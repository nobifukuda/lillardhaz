test_that("lognormal x lognormal, correlated, recovers true parameters", {
  set.seed(42)
  n <- 20000
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  b10 <- 1.0; b11 <- -0.3; lnsig1 <- -0.2
  b20 <- 0.6; b21 <- 0.4; lnsig2 <- 0.1
  rho_true <- 0.5

  u1 <- rnorm(n)
  u2 <- rho_true * u1 + sqrt(1 - rho_true^2) * rnorm(n)
  t1 <- exp(b10 + b11 * x1 + exp(lnsig1) * u1)
  t2 <- exp(b20 + b21 * x2 + exp(lnsig2) * u2)
  c1 <- -log(runif(n)) / 0.03
  c2 <- -log(runif(n)) / 0.03
  d <- data.frame(
    time1 = pmin(t1, c1), event1 = as.integer(t1 <= c1),
    time2 = pmin(t2, c2), event2 = as.integer(t2 <= c2),
    x1 = x1, x2 = x2
  )

  fit <- fit_lillardhaz("lognormal", "lognormal", d,
                         y1 = "time1", d1 = "event1", y2 = "time2", d2 = "event2",
                         x1 = "x1", x2 = "x2")

  b <- fit$coefficients
  expect_lt(abs(b["eq1:x1"] - b11), 0.05)
  expect_lt(abs(b["eq1:_cons"] - b10), 0.05)
  expect_lt(abs(b["eq2:x2"] - b21), 0.05)
  expect_lt(abs(b["ln_sigma1"] - lnsig1), 0.05)
  expect_lt(abs(b["ln_sigma2"] - lnsig2), 0.05)
  expect_lt(abs(fit$rho - rho_true), 0.05)
})

test_that("probit x lognormal, correlated, recovers true parameters", {
  set.seed(7)
  n <- 20000
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  g0 <- 0.3; g1 <- -0.5
  b20 <- 0.7; b21 <- 0.35; lnsig2 <- -0.15
  rho_true <- -0.4

  e1 <- rnorm(n)
  u2 <- rho_true * e1 + sqrt(1 - rho_true^2) * rnorm(n)
  y1 <- as.integer(g0 + g1 * x1 + e1 > 0)
  t2 <- exp(b20 + b21 * x2 + exp(lnsig2) * u2)
  c2 <- -log(runif(n)) / 0.03
  d <- data.frame(
    y1 = y1, time2 = pmin(t2, c2), event2 = as.integer(t2 <= c2), x1 = x1, x2 = x2
  )

  fit <- fit_lillardhaz("probit", "lognormal", d,
                         y1 = "y1", y2 = "time2", d2 = "event2",
                         x1 = "x1", x2 = "x2")

  b <- fit$coefficients
  expect_lt(abs(b["eq1:x1"] - g1), 0.05)
  expect_lt(abs(b["eq1:_cons"] - g0), 0.05)
  expect_lt(abs(b["eq2:x2"] - b21), 0.05)
  expect_lt(abs(fit$rho - rho_true), 0.05)
})

test_that("pgompertz x pgompertz, correlated, recovers true parameters", {
  set.seed(55)
  n <- 25000
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)

  b0 <- -1.3; b1 <- 0.3; alpha0 <- -0.4; alpha1 <- 0.35; alpha2 <- -0.1
  lev2_0 <- -0.9; lev2_1 <- -0.25; gamma0 <- -0.1; gamma1 <- 0.25; gamma2 <- -0.2
  rho_true <- 0.35

  L1a <- b0 + b1 * x1 + alpha0
  L1b_off <- alpha1 * 5
  L2a <- lev2_0 + lev2_1 * x2 + gamma0
  L2b_off <- gamma1 * 4

  e1 <- rnorm(n)
  e2 <- rho_true * e1 + sqrt(1 - rho_true^2) * rnorm(n)
  U1 <- pnorm(e1); U2 <- pnorm(e2)
  Htarget1 <- -log(1 - U1); Htarget2 <- -log(1 - U2)

  # With a negative terminal slope (alpha2, gamma2 < 0 here), the segment's
  # hazard decays toward 0 and its cumulative hazard is bounded as t -> Inf,
  # so a few extreme Htarget draws exceed what's reachable in finite time --
  # log() of a negative argument, i.e. NaN. This is a genuine defective-
  # distribution edge case (Gompertz with a negative slope), not a bug: those
  # observations simply never fail, so they're always censored (t = Inf).
  seg1a_true <- exp(L1a) / alpha1 * (exp(alpha1 * 5) - 1)
  t1 <- suppressWarnings(ifelse(Htarget1 < seg1a_true,
               log(1 + Htarget1 * alpha1 / exp(L1a)) / alpha1,
               5 + log(1 + (Htarget1 - seg1a_true) * alpha2 / exp(L1a + L1b_off)) / alpha2))
  t1[is.na(t1)] <- Inf

  seg2a_true <- exp(L2a) / gamma1 * (exp(gamma1 * 4) - 1)
  t2 <- suppressWarnings(ifelse(Htarget2 < seg2a_true,
               log(1 + Htarget2 * gamma1 / exp(L2a)) / gamma1,
               4 + log(1 + (Htarget2 - seg2a_true) * gamma2 / exp(L2a + L2b_off)) / gamma2))
  t2[is.na(t2)] <- Inf

  c1 <- -log(runif(n)) / 0.03
  c2 <- -log(runif(n)) / 0.03
  d <- data.frame(
    time1 = pmax(pmin(t1, c1), 1e-4), event1 = as.integer(t1 <= c1),
    time2 = pmax(pmin(t2, c2), 1e-4), event2 = as.integer(t2 <= c2),
    x1 = x1, x2 = x2
  )

  fit <- fit_lillardhaz("pgompertz", "pgompertz", d,
                         y1 = "time1", d1 = "event1", y2 = "time2", d2 = "event2",
                         x1 = "x1", x2 = "x2", nodes1 = 5, nodes2 = 4)

  b <- fit$coefficients
  expect_lt(abs(b["eq1:_cons"] - (b0 + alpha0)), 0.06)
  expect_lt(abs(b["eq1:x1"] - b1), 0.03)
  expect_lt(abs(b["s1_1"] - alpha1), 0.03)
  expect_lt(abs(b["s1_2"] - alpha2), 0.03)
  expect_lt(abs(b["eq2:_cons"] - (lev2_0 + gamma0)), 0.06)
  expect_lt(abs(b["eq2:x2"] - lev2_1), 0.03)
  expect_lt(abs(b["s2_1"] - gamma1), 0.03)
  expect_lt(abs(b["s2_2"] - gamma2), 0.03)
  expect_lt(abs(fit$rho - rho_true), 0.03)
})

test_that("nocorr reduces to two independent univariate fits", {
  set.seed(88)
  n <- 20000
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  b10 <- 1.0; b11 <- -0.3; lnsig1 <- -0.2
  b20 <- 0.6; b21 <- 0.4; lnsig2 <- 0.1

  t1 <- exp(b10 + b11 * x1 + exp(lnsig1) * rnorm(n))
  t2 <- exp(b20 + b21 * x2 + exp(lnsig2) * rnorm(n))
  c1 <- -log(runif(n)) / 0.03
  c2 <- -log(runif(n)) / 0.03
  d <- data.frame(
    time1 = pmin(t1, c1), event1 = as.integer(t1 <= c1),
    time2 = pmin(t2, c2), event2 = as.integer(t2 <= c2),
    x1 = x1, x2 = x2
  )

  fit <- fit_lillardhaz("lognormal", "lognormal", d,
                         y1 = "time1", d1 = "event1", y2 = "time2", d2 = "event2",
                         x1 = "x1", x2 = "x2", corr = FALSE)

  expect_identical(fit$corr, FALSE)
  b <- fit$coefficients
  expect_lt(abs(b["eq1:x1"] - b11), 0.05)
  expect_lt(abs(b["eq2:x2"] - b21), 0.05)
})
