## Tests for the full-refit bootstrap of the tilted estimators.

test_that("internal bootstrap returns the documented structure", {
  set.seed(1)
  d <- sim_internal_dgp(200)
  spec <- tilt_spec("ATT")
  b <- bootstrap_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, n_boot = 15)

  expect_named(b, c("estimate", "estimates", "se", "CI_normal", "CI_percentile",
                    "n_failed", "n_nonconverged"), ignore.order = TRUE)
  expect_length(b$estimates, 15)
  expect_true(is.finite(b$se) && b$se > 0)
  expect_length(b$CI_normal, 2)
  expect_lt(b$CI_normal[1], b$estimate)
  expect_gt(b$CI_normal[2], b$estimate)
  ## point estimate equals the plain entry-point result
  plain <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)
  expect_equal(b$estimate, plain$TMLE)
})

test_that("internal bootstrap SE is in the same ballpark as the EIF SE", {
  skip_on_cran()
  set.seed(2)
  d <- sim_internal_dgp(400)
  spec <- tilt_spec("ATO")
  b <- bootstrap_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, n_boot = 100)
  eif_se <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)$TMLE_se
  expect_gt(b$se / eif_se, 0.4)
  expect_lt(b$se / eif_se, 2.5)
})

test_that("bootstrap passes options through `...` and respects conf_level", {
  set.seed(3)
  d <- sim_internal_dgp(200)
  spec <- tilt_spec("ATC")
  b90 <- bootstrap_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, n_boot = 20,
                                         conf_level = 0.90, gbounds = c(0.05, 0.95))
  b99 <- bootstrap_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, n_boot = 20,
                                         conf_level = 0.99, gbounds = c(0.05, 0.95))
  expect_gt(diff(b99$CI_normal), diff(b90$CI_normal))
})

test_that("external bootstrap returns the documented structure", {
  set.seed(4)
  d <- sim_external_dgp(100, 200)
  spec <- tilt_spec_external("ATE-O")
  b <- bootstrap_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot, n_boot = 15)

  expect_length(b$estimates, 15)
  expect_true(is.finite(b$se) && b$se > 0)
  expect_lt(b$CI_normal[1], b$estimate)
  expect_gt(b$CI_normal[2], b$estimate)
  plain <- estimate_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot)
  expect_equal(b$estimate, plain$TMLE)
})

test_that("external bootstrap works with a separate Wg", {
  set.seed(5)
  d <- sim_external_dgp(100, 200)
  spec <- tilt_spec_external("ATE-I")
  Wg <- d$Ws[which(d$D == 1), c("X1", "X2"), drop = FALSE]
  b <- bootstrap_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot,
                                       Wg = Wg, n_boot = 10)
  expect_true(is.finite(b$se))
})

test_that("a healthy run reports n_failed = 0", {
  set.seed(6)
  d <- sim_internal_dgp(100)
  spec <- tilt_spec("ATT")
  b <- bootstrap_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, n_boot = 10)
  expect_equal(b$n_failed, 0L)
})
