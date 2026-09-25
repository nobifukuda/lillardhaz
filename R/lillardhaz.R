#' @importFrom stats pnorm dnorm qnorm sd coef logLik setNames
#' @importFrom pbivnorm pbivnorm
NULL

# ---------------------------------------------------------------------------
# Internal: piecewise-Gompertz survival/density, general K segments.
# theta: per-observation location index (vector, length n)
# slopes: per-segment slopes (vector, length K)
# nodes: interior nodes, ascending (vector, length K-1)
# t: per-observation duration (vector, length n)
# Direct port of _lillardhaz_pgomp.ado's recursion (Stata companion package).
# ---------------------------------------------------------------------------
.pgompertz_sf <- function(theta, slopes, nodes, t) {
  K <- length(slopes)
  n <- length(theta)

  levels_k <- numeric(K)
  startnodes_k <- numeric(K)
  cumlevel <- 0
  prevnode <- 0
  for (k in seq_len(K)) {
    levels_k[k] <- cumlevel
    startnodes_k[k] <- prevnode
    if (k < K) {
      thisnode <- nodes[k]
      seglen <- thisnode - prevnode
      cumlevel <- cumlevel + slopes[k] * seglen
      prevnode <- thisnode
    }
  }

  m <- rep(K, n)
  if (K > 1) {
    for (k in seq_len(K - 1)) {
      thisnode <- nodes[k]
      m[t < thisnode & m == K] <- k
    }
  }

  H <- numeric(n)
  h <- numeric(n)
  for (k in seq_len(K)) {
    lev_k <- levels_k[k]
    sl_k  <- slopes[k]
    sn_k  <- startnodes_k[k]

    if (k < K) {
      thisnode <- nodes[k]
      seglen <- thisnode - sn_k
      idx_full <- which(m > k)
      if (length(idx_full) > 0) {
        contrib <- if (abs(sl_k) < 1e-8) exp(lev_k) * seglen else
          exp(lev_k) / sl_k * (exp(sl_k * seglen) - 1)
        H[idx_full] <- H[idx_full] + exp(theta[idx_full]) * contrib
      }
    }

    idx_this <- which(m == k)
    if (length(idx_this) > 0) {
      tt <- t[idx_this] - sn_k
      contrib2 <- if (abs(sl_k) < 1e-8) exp(lev_k) * tt else
        exp(lev_k) / sl_k * (exp(sl_k * tt) - 1)
      H[idx_this] <- H[idx_this] + exp(theta[idx_this]) * contrib2
      h[idx_this] <- exp(theta[idx_this] + lev_k + sl_k * tt)
    }
  }

  S <- exp(-H)
  f <- h * S
  list(S = S, f = f, h = h)
}

# ---------------------------------------------------------------------------
# Internal: marginal S(t), f(t), and copula residual z for one hazard
# equation, dispatching on "lognormal" vs "pgompertz".
# ---------------------------------------------------------------------------
.hazard_marginal <- function(type, theta, t, lnsigma = NULL, slopes = NULL, nodes = NULL) {
  if (type == "lognormal") {
    sigma <- exp(lnsigma)
    z <- (log(t) - theta) / sigma
    S <- 1 - pnorm(z)
    f <- dnorm(z) / (sigma * t)
  } else {
    sf <- .pgompertz_sf(theta, slopes, nodes, t)
    S <- sf$S
    f <- sf$f
    z <- qnorm(pmin(pmax(1 - S, 1e-12), 1 - 1e-12))
  }
  list(S = S, f = f, z = z)
}

# ---------------------------------------------------------------------------
# Internal: unpack the parameter vector into named pieces, given the model
# layout (mirrors the equation ordering used by the Stata companion's d0
# evaluator: eq1 location+ancillary, eq2 location+ancillary, atanh_rho).
# ---------------------------------------------------------------------------
.param_layout <- function(eq1type, eq2type, cn1, cn2, K1 = NULL, K2 = NULL, corr = TRUE) {
  nm <- character(0)
  nm <- c(nm, paste0("eq1:", cn1))
  if (eq1type == "lognormal") nm <- c(nm, "ln_sigma1")
  else if (eq1type == "pgompertz") nm <- c(nm, paste0("s1_", seq_len(K1)))

  nm <- c(nm, paste0("eq2:", cn2))
  if (eq2type == "lognormal") nm <- c(nm, "ln_sigma2")
  else nm <- c(nm, paste0("s2_", seq_len(K2)))

  if (corr) nm <- c(nm, "atanh_rho")
  nm
}

.BAD_PENALTY <- -1e10

# Wraps .loglik_lillardhaz_core so that wild optimizer trial steps (e.g. an
# exploding pgompertz slope during a BFGS line search, which can overflow
# exp() into Inf/NaN and make pbivnorm() error on a non-finite input) return
# a large finite penalty instead of crashing maxLik -- the optimizer then
# simply steers away from that region, same fix used in the companion
# Python package for an analogous overflow.
.loglik_lillardhaz <- function(par, ...) {
  val <- tryCatch(
    suppressWarnings(.loglik_lillardhaz_core(par, ...)),
    error = function(e) NA_real_
  )
  if (!is.finite(val)) .BAD_PENALTY else val
}

.loglik_lillardhaz_core <- function(par, eq1type, eq2type, X1, X2,
                                     y1, d1, y2, d2, nodes1, nodes2, corr) {
  p1 <- ncol(X1); p2 <- ncol(X2)
  idx <- 1
  b1 <- par[idx:(idx + p1 - 1)]; idx <- idx + p1
  theta1 <- as.vector(X1 %*% b1)

  lnsigma1 <- NULL; slopes1 <- NULL
  if (eq1type == "lognormal") {
    lnsigma1 <- par[idx]; idx <- idx + 1
  } else if (eq1type == "pgompertz") {
    K1 <- length(nodes1) + 1
    slopes1 <- par[idx:(idx + K1 - 1)]; idx <- idx + K1
  }

  b2 <- par[idx:(idx + p2 - 1)]; idx <- idx + p2
  theta2 <- as.vector(X2 %*% b2)

  lnsigma2 <- NULL; slopes2 <- NULL
  if (eq2type == "lognormal") {
    lnsigma2 <- par[idx]; idx <- idx + 1
  } else {
    K2 <- length(nodes2) + 1
    slopes2 <- par[idx:(idx + K2 - 1)]; idx <- idx + K2
  }

  rho <- if (corr) tanh(par[idx]) else 0
  a <- sqrt(1 - rho^2)

  n <- length(y2)
  cont <- numeric(n)

  if (eq1type == "probit") {
    m2 <- .hazard_marginal(eq2type, theta2, y2, lnsigma2, slopes2, nodes2)
    f2 <- m2$f; z2 <- m2$z

    i <- d1 == 1 & d2 == 1
    cont[i] <- f2[i] * pnorm((theta1[i] + rho * z2[i]) / a)
    i <- d1 == 1 & d2 == 0
    cont[i] <- pnorm(theta1[i]) - pbivnorm(theta1[i], z2[i], rep(-rho, sum(i)))
    i <- d1 == 0 & d2 == 1
    cont[i] <- f2[i] * pnorm(-(theta1[i] + rho * z2[i]) / a)
    i <- d1 == 0 & d2 == 0
    cont[i] <- pbivnorm(-theta1[i], -z2[i], rep(-rho, sum(i)))
  } else {
    m1 <- .hazard_marginal(eq1type, theta1, y1, lnsigma1, slopes1, nodes1)
    m2 <- .hazard_marginal(eq2type, theta2, y2, lnsigma2, slopes2, nodes2)
    f1 <- m1$f; z1 <- m1$z
    f2 <- m2$f; z2 <- m2$z

    i <- d1 == 1 & d2 == 1
    cont[i] <- exp(-(z1[i]^2 - 2 * rho * z1[i] * z2[i] + z2[i]^2) / (2 * a^2)) / (2 * pi * a) *
      f1[i] / dnorm(z1[i]) * f2[i] / dnorm(z2[i])
    i <- d1 == 1 & d2 == 0
    cont[i] <- f1[i] * (1 - pnorm((z2[i] - rho * z1[i]) / a))
    i <- d1 == 0 & d2 == 1
    cont[i] <- f2[i] * (1 - pnorm((z1[i] - rho * z2[i]) / a))
    i <- d1 == 0 & d2 == 0
    cont[i] <- pbivnorm(-z1[i], -z2[i], rep(rho, sum(i)))
  }

  sum(log(pmax(cont, 1e-300)))
}

#' Fit a Lillard-style simultaneous-equations hazard/probit model
#'
#' Fits a two-equation simultaneous model in the style of Lillard (1993): a
#' pair of processes, each either a binary probit outcome or a continuous-
#' time hazard duration, linked through a single correlation parameter
#' between their underlying error terms and estimated jointly by maximum
#' likelihood via a Gaussian copula.
#'
#' @param eq1 one of \code{"probit"}, \code{"lognormal"}, \code{"pgompertz"}.
#' @param eq2 one of \code{"lognormal"}, \code{"pgompertz"}.
#' @param data a data frame.
#' @param y1 name of equation 1's outcome: the binary 0/1 variable if
#'   \code{eq1="probit"}, otherwise the duration variable.
#' @param d1 name of equation 1's failure indicator (ignored if
#'   \code{eq1="probit"}).
#' @param y2 name of equation 2's duration variable.
#' @param d2 name of equation 2's failure indicator.
#' @param x1 character vector of covariate names for equation 1's location
#'   index (may be empty for an intercept-only equation).
#' @param x2 character vector of covariate names for equation 2's location
#'   index.
#' @param nodes1,nodes2 interior nodes (ascending numeric vectors) for a
#'   piecewise-Gompertz equation 1/equation 2; required exactly when that
#'   equation is \code{"pgompertz"}.
#' @param corr logical; if \code{FALSE}, fixes the copula correlation at 0
#'   and fits the two equations independently (see the manual for why this
#'   reduces exactly to two separate univariate fits).
#' @param start optional named numeric vector of starting values.
#' @param method optimization method passed to \code{maxLik::maxLik}.
#'
#' @return An object of class \code{"lillardhaz"}.
#' @examples
#' set.seed(1)
#' n <- 2000
#' x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
#' t1 <- exp(1 - 0.3 * x1 + 0.5 * rnorm(n))
#' c1 <- -log(runif(n)) / 0.05
#' t2 <- exp(0.6 + 0.4 * x2 + 0.6 * rnorm(n))
#' c2 <- -log(runif(n)) / 0.05
#' d <- data.frame(
#'   time1 = pmin(t1, c1), event1 = as.integer(t1 <= c1),
#'   time2 = pmin(t2, c2), event2 = as.integer(t2 <= c2),
#'   x1 = x1, x2 = x2
#' )
#' fit <- fit_lillardhaz("lognormal", "lognormal", d,
#'                        y1 = "time1", d1 = "event1", y2 = "time2", d2 = "event2",
#'                        x1 = "x1", x2 = "x2")
#' summary(fit)
#' @export
fit_lillardhaz <- function(eq1, eq2, data, y1, d1 = NULL, y2, d2,
                            x1 = character(0), x2 = character(0),
                            nodes1 = NULL, nodes2 = NULL, corr = TRUE,
                            start = NULL, method = "BFGS") {
  eq1 <- match.arg(eq1, c("probit", "lognormal", "pgompertz"))
  eq2 <- match.arg(eq2, c("lognormal", "pgompertz"))
  if (eq1 == "pgompertz" && is.null(nodes1)) stop("nodes1 is required when eq1='pgompertz'")
  if (eq2 == "pgompertz" && is.null(nodes2)) stop("nodes2 is required when eq2='pgompertz'")

  data <- as.data.frame(data)
  X1 <- if (length(x1) > 0) cbind(1, as.matrix(data[, x1, drop = FALSE])) else matrix(1, nrow(data), 1)
  colnames(X1) <- c("_cons", x1)
  X2 <- if (length(x2) > 0) cbind(1, as.matrix(data[, x2, drop = FALSE])) else matrix(1, nrow(data), 1)
  colnames(X2) <- c("_cons", x2)

  y1v <- data[[y1]]
  d1v <- if (eq1 == "probit") NULL else data[[d1]]
  y2v <- data[[y2]]
  d2v <- data[[d2]]

  if (eq1 == "probit") d1v <- y1v  # loglik dispatch uses d1 as the y1 binary outcome when eq1=="probit"

  K1 <- if (eq1 == "pgompertz") length(nodes1) + 1 else NA
  K2 <- if (eq2 == "pgompertz") length(nodes2) + 1 else NA

  nm <- .param_layout(eq1, eq2, colnames(X1), colnames(X2), K1, K2, corr)

  if (is.null(start)) {
    start <- setNames(numeric(length(nm)), nm)
    if (eq1 != "probit") {
      start[grepl("^eq1:", nm)] <- 0
      start[1] <- mean(log(y1v))
      if (eq1 == "lognormal") start["ln_sigma1"] <- log(sd(log(y1v)))
      if (eq1 == "pgompertz") start[grepl("^s1_", nm)] <- 0.05
    } else {
      start[1] <- 0
    }
    start[grepl("^eq2:", nm)][1] <- mean(log(y2v))
    if (eq2 == "lognormal") start["ln_sigma2"] <- log(sd(log(y2v)))
    if (eq2 == "pgompertz") start[grepl("^s2_", nm)] <- 0.05
    if (corr) start["atanh_rho"] <- 0
  }

  fit <- maxLik::maxLik(
    logLik = .loglik_lillardhaz,
    start = start,
    method = method,
    eq1type = eq1, eq2type = eq2, X1 = X1, X2 = X2,
    y1 = y1v, d1 = d1v, y2 = y2v, d2 = d2v,
    nodes1 = nodes1, nodes2 = nodes2, corr = corr
  )

  b <- stats::coef(fit)
  names(b) <- nm
  se <- tryCatch(sqrt(diag(solve(-fit$hessian))), error = function(e) rep(NA_real_, length(b)))
  names(se) <- nm

  rho <- if (corr) tanh(b["atanh_rho"]) else 0
  rho_se <- if (corr) (1 - rho^2) * se["atanh_rho"] else NA_real_

  structure(list(
    coefficients = b, se = se, loglik = as.numeric(stats::logLik(fit)),
    fit = fit, eq1type = eq1, eq2type = eq2, nodes1 = nodes1, nodes2 = nodes2,
    corr = corr, rho = unname(rho), rho_se = unname(rho_se),
    x1 = x1, x2 = x2, y1name = y1, y2name = y2, n = nrow(data), call = match.call()
  ), class = "lillardhaz")
}

#' @export
print.lillardhaz <- function(x, ...) {
  cat("Lillard-style simultaneous-equations model\n")
  cat(sprintf("eq1: %s   eq2: %s   corr: %s\n", x$eq1type, x$eq2type, x$corr))
  cat(sprintf("Observations: %d   Log-likelihood: %.3f\n\n", x$n, x$loglik))
  print(summary(x))
  invisible(x)
}

#' @export
summary.lillardhaz <- function(object, ...) {
  z <- object$coefficients / object$se
  p <- 2 * pnorm(-abs(z))
  tab <- data.frame(
    Estimate = object$coefficients, `Std. Error` = object$se,
    `z value` = z, `Pr(>|z|)` = p, check.names = FALSE
  )
  if (object$corr) {
    tab <- rbind(tab, rho = c(object$rho, object$rho_se, NA, NA))
  }
  structure(list(table = tab, object = object), class = "summary.lillardhaz")
}

#' @export
print.summary.lillardhaz <- function(x, ...) {
  print(round(x$table, 4))
  invisible(x)
}

#' Predicted survival/density/index from a fitted lillardhaz model
#'
#' @param object a fitted \code{"lillardhaz"} object.
#' @param newdata optional data frame; defaults to the estimation sample via
#'   \code{object$call}.
#' @param type one of \code{"surv2"} (default), \code{"surv1"}, \code{"dens1"},
#'   \code{"dens2"}, \code{"pr1"}, \code{"xb1"}, \code{"xb2"}.
#' @param ... unused.
#' @export
predict.lillardhaz <- function(object, newdata, type = "surv2", ...) {
  type <- match.arg(type, c("surv2", "surv1", "dens1", "dens2", "pr1", "xb1", "xb2"))
  b <- object$coefficients
  x1 <- object$x1; x2 <- object$x2

  X1 <- if (length(x1) > 0) cbind(1, as.matrix(newdata[, x1, drop = FALSE])) else matrix(1, nrow(newdata), 1)
  X2 <- if (length(x2) > 0) cbind(1, as.matrix(newdata[, x2, drop = FALSE])) else matrix(1, nrow(newdata), 1)
  theta1 <- as.vector(X1 %*% b[grepl("^eq1:", names(b))])
  theta2 <- as.vector(X2 %*% b[grepl("^eq2:", names(b))])

  if (type == "xb1") return(theta1)
  if (type == "xb2") return(theta2)
  if (type == "pr1") {
    if (object$eq1type != "probit") stop("pr1 is only available when eq1='probit'")
    return(pnorm(theta1))
  }

  if (type %in% c("surv1", "dens1")) {
    if (object$eq1type == "probit") stop("surv1/dens1 are not available when eq1='probit'")
    m1 <- if (object$eq1type == "lognormal")
      .hazard_marginal("lognormal", theta1, newdata[[object$y1name]], lnsigma = b["ln_sigma1"])
    else
      .hazard_marginal("pgompertz", theta1, newdata[[object$y1name]],
                        slopes = b[grepl("^s1_", names(b))], nodes = object$nodes1)
    return(if (type == "surv1") m1$S else m1$f)
  }

  y2name <- object$y2name
  m2 <- if (object$eq2type == "lognormal")
    .hazard_marginal("lognormal", theta2, newdata[[y2name]], lnsigma = b["ln_sigma2"])
  else
    .hazard_marginal("pgompertz", theta2, newdata[[y2name]],
                      slopes = b[grepl("^s2_", names(b))], nodes = object$nodes2)
  if (type == "dens2") m2$f else m2$S
}
