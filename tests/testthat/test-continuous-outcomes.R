## Tests for continuous outcomes (added in 0.2.0).
## A continuous Y is rescaled to [0, 1] with Ybounds, the TMLE is computed on that
## scale, and estimates / SEs / CIs are mapped back to the original scale.

## ---- Data-generating processes ---------------------------------------------

## Internal: Y = 1 + X1 - 0.5 X2 + A * (1 + 0.5 X1) + N(0, 1).
## Conditional effect tau(X) = 1 + 0.5 X1, so the ATE is 1.
sim_continuous_internal <- function(n = 1000) {
  X1 <- stats::rnorm(n)
  X2 <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.4 * X1 - 0.3 * X2))
  Y <- 1 + X1 - 0.5 * X2 + A * (1 + 0.5 * X1) + stats::rnorm(n)
  list(Y = Y, A = A, W = data.frame(X1 = X1, X2 = X2))
}

## True internally tilted effect: E[eta(g) tau(X)] / E[eta(g)], by Monte Carlo.
true_continuous_internal <- function(eta, n_large = 2e5) {
  set.seed(99)
  X1 <- stats::rnorm(n_large)
  X2 <- stats::rnorm(n_large)
  h <- eta(stats::plogis(0.4 * X1 - 0.3 * X2))
  sum(h * (1 + 0.5 * X1)) / sum(h)
}

## External: internal (D = 1) and external (D = 0) samples from one super-population,
## selection s(X) = plogis(0.8 X1 - 0.5 X2); same outcome model as above.
s_true <- function(X1, X2) stats::plogis(0.8 * X1 - 0.5 * X2)

sim_continuous_external <- function(n = 600, m = 600) {
  N_super <- 6 * (n + m)
  X1 <- stats::rnorm(N_super)
  X2 <- stats::rnorm(N_super)
  D_pool <- stats::rbinom(N_super, 1, s_true(X1, X2))
  i1 <- sample(which(D_pool == 1), n)
  i0 <- sample(which(D_pool == 0), m)
  W1 <- data.frame(X1 = X1[i1], X2 = X2[i1])
  A <- stats::rbinom(n, 1, stats::plogis(0.4 * W1$X1 - 0.3 * W1$X2))
  Y <- 1 + W1$X1 - 0.5 * W1$X2 + A * (1 + 0.5 * W1$X1) + stats::rnorm(n)
  list(Y = Y, A = A, D = c(rep(1L, n), rep(0L, m)),
       Ws = data.frame(X1 = c(X1[i1], X1[i0]), X2 = c(X2[i1], X2[i0])))
}

## True externally tilted effect: E[eta(s) tau(X) | D = 1] / E[eta(s) | D = 1].
true_continuous_external <- function(eta, n_large = 4e5) {
  set.seed(98)
  X1 <- stats::rnorm(n_large)
  X2 <- stats::rnorm(n_large)
  s <- s_true(X1, X2)
  D <- stats::rbinom(n_large, 1, s)
  h <- eta(s[D == 1])
  sum(h * (1 + 0.5 * X1[D == 1])) / sum(h)
}

## ---- Binary outcomes are left alone ----------------------------------------

test_that("a binary outcome is not rescaled", {
  set.seed(1)
  d <- sim_internal_dgp(300)
  spec <- tilt_spec("ATT")
  out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)
  expect_equal(out$Ybounds, c(0, 1))
  out01 <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                          Ybounds = c(0, 1))
  expect_equal(out01$TMLE, out$TMLE)
})

## ---- Results follow a change of units exactly ------------------------------

test_that("estimates and SEs scale with a linear change of units (all entry points)", {
  set.seed(2)
  d <- sim_continuous_internal(400)
  Y2 <- 10 + 5 * d$Y            # same outcome in different units
  spec <- tilt_spec("ATO")

  a <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)
  b <- estimate_tau_h_internal_tilted(Y2, d$A, d$W, spec$eta, spec$eta_dot)
  expect_equal(b$TMLE, 5 * a$TMLE, tolerance = 1e-6)
  expect_equal(b$TMLE_se, 5 * a$TMLE_se, tolerance = 1e-6)
  expect_equal(b$TMLE_CI, 5 * a$TMLE_CI, tolerance = 1e-6)
  expect_equal(b$AIPW, 5 * a$AIPW, tolerance = 1e-6)

  set.seed(3); a <- estimate_tau_h_fixed(d$Y, d$A, d$W, h = rep(1, 400), n_boot = 20)
  set.seed(3); b <- estimate_tau_h_fixed(Y2, d$A, d$W, h = rep(1, 400), n_boot = 20)
  expect_equal(b$TMLE, 5 * a$TMLE, tolerance = 1e-6)
  expect_equal(b$Bootstrap_SE, 5 * a$Bootstrap_SE, tolerance = 1e-6)

  set.seed(4)
  e <- sim_continuous_external(300, 300)
  spec_o <- tilt_spec_external("ATE-O")
  a <- estimate_tau_h_external_tilted(e$Y, e$A, e$D, e$Ws, spec_o$eta, spec_o$eta_dot)
  b <- estimate_tau_h_external_tilted(10 + 5 * e$Y, e$A, e$D, e$Ws, spec_o$eta, spec_o$eta_dot)
  expect_equal(b$TMLE, 5 * a$TMLE, tolerance = 1e-6)
  expect_equal(b$TMLE_se, 5 * a$TMLE_se, tolerance = 1e-6)
})

## ---- Approximately unbiased --------------------------------------------------

test_that("internally tilted TMLE is close to the truth for a continuous outcome", {
  set.seed(5)
  d <- sim_continuous_internal(2000)
  for (est in c("ATE", "ATT", "ATO")) {
    spec <- tilt_spec(est)
    out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)
    truth <- true_continuous_internal(spec$eta)
    expect_true(out$TMLE_converged)
    expect_lt(abs(out$TMLE - truth), 4 * out$TMLE_se)
    ## on the original scale, the SE of a difference in means of N(., 1)-ish data is ~0.05
    expect_gt(out$TMLE_se, 0.02)
    expect_lt(out$TMLE_se, 0.2)
  }
})

test_that("externally tilted TMLE is close to the truth for a continuous outcome", {
  set.seed(6)
  e <- sim_continuous_external(1500, 1500)
  for (est in c("ATE-I", "ATE-O")) {
    spec <- tilt_spec_external(est)
    out <- estimate_tau_h_external_tilted(e$Y, e$A, e$D, e$Ws, spec$eta, spec$eta_dot)
    truth <- true_continuous_external(spec$eta)
    expect_true(out$TMLE_converged)
    expect_lt(abs(out$TMLE - truth), 4 * out$TMLE_se)
  }
})

## ---- Works with the other features --------------------------------------------

test_that("continuous outcomes work with missing outcomes and the bootstrap", {
  set.seed(7)
  d <- sim_continuous_internal(400)
  spec <- tilt_spec("ATT")
  Y_miss <- d$Y
  Y_miss[stats::rbinom(400, 1, stats::plogis(1 + 0.5 * d$W$X1)) == 0] <- NA
  out <- estimate_tau_h_internal_tilted(Y_miss, d$A, d$W, spec$eta, spec$eta_dot)
  expect_true(is.finite(out$TMLE))
  expect_equal(out$Ybounds, range(Y_miss, na.rm = TRUE))

  b <- bootstrap_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, n_boot = 10)
  expect_true(is.finite(b$se) && b$se > 0)
  expect_equal(b$estimate,
               estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)$TMLE)

  e <- sim_continuous_external(200, 200)
  spec_o <- tilt_spec_external("ATE-O")
  b2 <- bootstrap_tau_h_external_tilted(e$Y, e$A, e$D, e$Ws, spec_o$eta, spec_o$eta_dot,
                                        n_boot = 10)
  expect_true(is.finite(b2$se) && b2$se > 0)
})

test_that("continuous outcomes work with SuperLearner", {
  skip_if_not_installed("SuperLearner")
  set.seed(8)
  d <- sim_continuous_internal(300)
  spec <- tilt_spec("ATE")
  out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                        Q_method = "SL", Q_SL.library = c("SL.glm", "SL.mean"))
  expect_true(is.finite(out$TMLE))
  expect_lt(abs(out$TMLE - 1), 0.5)
})

## ---- Input checks --------------------------------------------------------------

test_that("Ybounds and Y are validated", {
  set.seed(9)
  d <- sim_continuous_internal(100)
  spec <- tilt_spec("ATE")
  ## user bounds must contain every observed Y
  expect_error(estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                              Ybounds = c(0, 1)), "within `Ybounds`")
  ## wider known bounds are fine
  out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                        Ybounds = c(-20, 20))
  expect_equal(out$Ybounds, c(-20, 20))
  ## malformed bounds and a constant outcome are rejected
  expect_error(estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                              Ybounds = c(5, -5)), "`Ybounds` must be")
  expect_error(estimate_tau_h_internal_tilted(rep(3, 100), d$A, d$W, spec$eta, spec$eta_dot),
               "constant outcome")
  ## non-numeric Y is rejected
  expect_error(estimate_tau_h_internal_tilted(as.character(d$Y), d$A, d$W, spec$eta, spec$eta_dot),
               "`Y` must be")
  ## fit_Q() itself needs Y on the [0, 1] scale
  expect_error(fit_Q(d$Y, d$A, d$W), "needs `Y` in \\[0, 1\\]")
})
