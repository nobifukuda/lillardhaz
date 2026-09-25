# lillardhaz (R)

An R package for **simultaneous-equations hazard/probit models with a
Gaussian-copula correlation**, in the style of Lillard (1993): a pair of
processes — each a binary probit outcome or a continuous-time hazard
duration — linked through a single correlation parameter between their
underlying error terms, estimated jointly by maximum likelihood.

Four model families are supported via `eq1`/`eq2` in `fit_lillardhaz()`:

| eq1 | eq2 | typical use |
|---|---|---|
| `probit` | `lognormal` | a binary decision correlated with a log-normal AFT duration |
| `lognormal` | `lognormal` | two correlated log-normal AFT durations |
| `probit` | `pgompertz` | a binary decision correlated with a flexible piecewise-Gompertz duration |
| `pgompertz` | `pgompertz` | two correlated piecewise-Gompertz durations |

`pgompertz` is a piecewise-linear-in-time log-hazard (piecewise Gompertz)
with an arbitrary, user-specified number of segments (`nodes1`/`nodes2`).
Every combination can also be fit independently with `corr = FALSE` (rho
fixed at 0), which — as derived in the manual — reduces algebraically to
two separate univariate fits, giving a natural nested baseline for a
likelihood-ratio test of correlation.

A companion [Stata package](https://github.com/nobifukuda/lillardhaz-stata)
implements the identical models. See
[`docs/manual.html`](https://htmlpreview.github.io/?https://github.com/nobifukuda/lillardhaz/blob/main/docs/manual.html)
for the full model derivation (the joint likelihood for every eq1×eq2
combination, the piecewise-Gompertz closed-form hazard, the independence
reduction, and simulation-based verification results), and `?fit_lillardhaz`
for the function reference.

## Author

**Nobutaka Fukuda**, Tohoku University — <nobutaka.fukuda@tohoku.ac.jp>

## Installation

```r
# install.packages("devtools")
devtools::install_github("nobifukuda/lillardhaz")
```

(Once accepted to CRAN — see below — installation will simply be
`install.packages("lillardhaz")`.)

## Usage

```r
library(lillardhaz)

# Probit & log-normal hazard, correlated
fit <- fit_lillardhaz("probit", "lognormal", data,
                       y1 = "y1", y2 = "time2", d2 = "event2",
                       x1 = "x1", x2 = "x2")
summary(fit)
fit$rho; fit$rho_se

# Piecewise Gompertz & piecewise Gompertz, correlated, 2 interior nodes each side
fit <- fit_lillardhaz("pgompertz", "pgompertz", data,
                       y1 = "time1", d1 = "event1", y2 = "time2", d2 = "event2",
                       x1 = "x1", x2 = "x2", nodes1 = 5, nodes2 = 4)

# Same, but independent (nested test of rho = 0)
fit0 <- fit_lillardhaz("pgompertz", "pgompertz", data,
                        y1 = "time1", d1 = "event1", y2 = "time2", d2 = "event2",
                        x1 = "x1", x2 = "x2", nodes1 = 5, nodes2 = 4, corr = FALSE)
2 * (fit$loglik - fit0$loglik)   # LR test of rho = 0, ~chisq(1)

# Predictions
predict(fit, data, type = "surv2")   # fitted eq2 survival (default)
predict(fit, data, type = "dens2")   # fitted eq2 density
predict(fit, data, type = "xb1")     # eq1 linear index
```

## Verification

Every eq1×eq2 combination, correlated and independent, was verified by
simulating 20,000–25,000 observations under known parameters and confirming
`fit_lillardhaz()` recovers them; the scripts run as part of `R CMD check`
via `tests/testthat/test-fit.R`. See the manual's verification section for
a summary table. Optimizer robustness against overflow during BFGS trial
steps (mirroring a fix already used in the companion
[splitpopsurv](https://github.com/nobifukuda/splitpopsurv) Python package)
is documented there as well.

## License

MIT — see [LICENSE.md](LICENSE.md).

## References

Lillard, L. A. 1993. Simultaneous equations for hazards: Marriage duration
and fertility timing. *Journal of Econometrics* 56(1–2): 189–217.

Waite, L. J., & Lillard, L. A. 1991. Children and marital disruption.
*American Journal of Sociology* 96(4): 930–953.
