# lillardhaz 1.0.0

* Initial release. Fits two-equation simultaneous hazard/probit models
  (probit, log-normal, or piecewise-Gompertz equations) linked by a
  Gaussian-copula correlation, via `fit_lillardhaz()`.
* `print`, `summary`, and `predict` S3 methods for `"lillardhaz"` objects.
* Verified by simulation against known parameters for every eq1×eq2
  combination, correlated and independent (`corr = FALSE`).
