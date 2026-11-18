## Shared data-generating process for the missing-outcome tests.

## Small DGP: Y depends on A, X1, X2; the outcome is missing at random, with a
## probability that depends on X1 and A (so a complete-case analysis is not enough).
sim_missing_dgp <- function(n = 500, missing = TRUE) {
  X1 <- stats::rnorm(n)
  X2 <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * X1))
  Y <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.8 * X1 - 0.4 * X2))
  if (missing) {
    obs <- stats::rbinom(n, 1, stats::plogis(0.8 + 1.0 * X1 - 0.6 * A))
    Y[obs == 0] <- NA
  }
  list(Y = Y, A = A, W = data.frame(X1 = X1, X2 = X2))
}

## Truth for the tilted ATE with h = eta(g(X)), g(x) = plogis(0.5 * X1) (the true treatment
## propensity in sim_missing_dgp). eta = function(g) ...; eta = function(g) rep(1, length(g)) is the ATE.
true_tilted_missing <- function(eta, n_large = 4e5) {
  set.seed(777)
  X1 <- stats::rnorm(n_large); X2 <- stats::rnorm(n_large)
  Q1 <- stats::plogis(-0.3 + 0.5 + 0.8 * X1 - 0.4 * X2)
  Q0 <- stats::plogis(-0.3 + 0.8 * X1 - 0.4 * X2)
  h <- eta(stats::plogis(0.5 * X1))
  sum(h * (Q1 - Q0)) / sum(h)
}

true_ate_missing <- function(n_large = 4e5) true_tilted_missing(function(g) rep(1, length(g)), n_large)
