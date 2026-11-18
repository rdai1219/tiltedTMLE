## Tests for review fixes 1-3 (Q bounds, g re-clamp, converged flag).
## Meant to sit next to the existing tests (uses helper-dgp.R) and to be run
## against the R_v2/ versions of nuisance-fitting.R, internal-tilted.R,
## external-tilted.R.

## ---- Fix 1: Q bounds -------------------------------------------------------

test_that("fit_Q() keeps Q inside Qbounds", {
  set.seed(1)
  d <- sim_internal_dgp(300)

  ## default bounds: strictly inside (0, 1)
  Qfit <- fit_Q(d$Y, d$A, d$W)
  for (q in Qfit[c("QAW", "Q1W", "Q0W")]) {
    expect_true(all(q >= 1e-4 & q <= 1 - 1e-4))
  }

  ## tight custom bounds are respected
  Qfit2 <- fit_Q(d$Y, d$A, d$W, Qbounds = c(0.3, 0.6))
  for (q in Qfit2[c("QAW", "Q1W", "Q0W")]) {
    expect_true(all(q >= 0.3 & q <= 0.6))
  }
})

test_that("predict_Qhat() applies Qbounds", {
  set.seed(2)
  d <- sim_internal_dgp(300)
  Qfit <- fit_Q(d$Y, d$A, d$W)
  p <- predict_Qhat(Qfit$fit, 1, d$W, Qbounds = c(0.3, 0.6))
  expect_true(all(p >= 0.3 & p <= 0.6))
})

test_that("TMLE returns finite results even when initial Q is pinned at the bounds", {
  set.seed(3)
  d <- sim_internal_dgp(300)
  Qfit <- fit_Q(d$Y, d$A, d$W, Qbounds = c(0.01, 0.99))
  gfit <- fit_g(d$A, d$W)
  spec <- tilt_spec("ATO")
  out <- tmle_tilted_h(d$Y, d$A, Qfit$Q1W, Qfit$Q0W, Qfit$QAW,
                       gfit$g1W, gfit$g0W, spec$eta, spec$eta_dot)
  expect_true(is.finite(out$estimate))
  expect_true(is.finite(out$se))
})

## ---- Fix 2: g re-clamp in the internal loop --------------------------------

test_that("tmle_tilted_h() keeps targeted g inside gbounds", {
  set.seed(4)
  d <- sim_internal_dgp(400)
  Qfit <- fit_Q(d$Y, d$A, d$W)
  gfit <- fit_g(d$A, d$W, gbounds = c(0.05, 0.95))
  for (est in c("ATT", "ATC", "ATO")) {
    spec <- tilt_spec(est)
    out <- tmle_tilted_h(d$Y, d$A, Qfit$Q1W, Qfit$Q0W, Qfit$QAW,
                         gfit$g1W, gfit$g0W, spec$eta, spec$eta_dot,
                         gbounds = c(0.05, 0.95))
    expect_true(all(out$g1W >= 0.05 & out$g1W <= 0.95), label = est)
    expect_equal(out$g0W, 1 - out$g1W)
  }
})

test_that("estimate_tau_h_internal_tilted() passes gbounds down to the loop", {
  set.seed(5)
  d <- sim_internal_dgp(400)
  spec <- tilt_spec("ATC")
  out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                        gbounds = c(0.1, 0.9))
  g <- out$tmle_fit$g1W
  expect_true(all(g >= 0.1 & g <= 0.9))
})

## ---- Fix 3: shared stopping rule and converged flag ------------------------

test_that("internal ATE converges in one iteration and reports converged = TRUE", {
  set.seed(6)
  d <- sim_internal_dgp(400)
  spec <- tilt_spec("ATE")
  out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                        truncate = FALSE)
  expect_equal(out$TMLE_iter, 1L)
  expect_true(out$TMLE_converged)
})

test_that("internal TMLE warns and sets converged = FALSE when max_iter is too small", {
  set.seed(7)
  d <- sim_internal_dgp(400)
  spec <- tilt_spec("ATT")
  expect_warning(
    out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                          max_iter = 1, tol = 1e-12),
    "did not converge"
  )
  expect_false(out$TMLE_converged)
  expect_equal(out$TMLE_iter, 1L)
})

test_that("internal ATT/ATC/ATO converge within the default max_iter", {
  set.seed(8)
  d <- sim_internal_dgp(400)
  for (est in c("ATT", "ATC", "ATO")) {
    spec <- tilt_spec(est)
    out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)
    expect_true(out$TMLE_converged, label = est)
  }
})

test_that("external TMLE reports converged and warns when max_iter is too small", {
  set.seed(9)
  d <- sim_external_dgp(300, 600)
  spec <- tilt_spec_external("ATE-O")

  out <- estimate_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot)
  expect_true(out$TMLE_converged)

  expect_warning(
    out2 <- estimate_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot,
                                           max_iter = 1, tol = 1e-12),
    "did not converge"
  )
  expect_false(out2$TMLE_converged)
})

test_that("external ATE-E still runs without error and keeps Q inside bounds", {
  set.seed(10)
  d <- sim_external_dgp(300, 600)
  spec <- tilt_spec_external("ATE-E")
  out <- estimate_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot)
  q <- c(out$tmle_fit$Q1W, out$tmle_fit$Q0W)
  expect_true(all(q >= 1e-4 & q <= 1 - 1e-4))
  expect_true(is.finite(out$TMLE))
})
