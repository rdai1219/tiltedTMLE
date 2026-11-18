## Tests for missing outcomes in estimate_tau_h_internal_tilted() (NA in Y, MAR given A and W).
## Uses sim_missing_dgp() / true_tilted_missing() from helper-missing-dgp.R.

test_that("with no missing outcomes, internal results are unchanged", {
  set.seed(1)
  d <- sim_missing_dgp(300, missing = FALSE)
  spec <- tilt_spec("ATT")
  new <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)
  expect_null(new$pifit)

  ## manual pipeline = the original (complete-data) code path
  Qfit <- fit_Q(d$Y, d$A, d$W)
  gfit <- fit_g(d$A, d$W)
  old <- estimate_tau_h_internal_tilted_given_fits(d$Y, d$A, Qfit, gfit, spec$eta, spec$eta_dot,
                                                   gbounds = c(0.025, 0.975))
  expect_equal(new$TMLE, old$TMLE)
  expect_equal(new$TMLE_se, old$TMLE_se)
  expect_equal(new$AIPW, old$AIPW)
})

test_that("eif_tilted_h() with Delta = 1 and pi = 1 equals the default", {
  set.seed(2)
  d <- sim_missing_dgp(200, missing = FALSE)
  Qfit <- fit_Q(d$Y, d$A, d$W); gfit <- fit_g(d$A, d$W)
  spec <- tilt_spec("ATO")
  eta_g <- spec$eta(gfit$g1W)
  tau <- sum(eta_g * (Qfit$Q1W - Qfit$Q0W)) / sum(eta_g)
  a <- eif_tilted_h(d$Y, d$A, eta_g, spec$eta_dot(gfit$g1W), Qfit$Q1W, Qfit$Q0W,
                    gfit$g1W, gfit$g0W, tau)
  b <- eif_tilted_h(d$Y, d$A, eta_g, spec$eta_dot(gfit$g1W), Qfit$Q1W, Qfit$Q0W,
                    gfit$g1W, gfit$g0W, tau,
                    Delta = rep(1L, 200), pi1W = rep(1, 200), pi0W = rep(1, 200))
  expect_equal(a, b)
})

test_that("internal tilted estimators run with NA in Y and converge", {
  set.seed(3)
  d <- sim_missing_dgp(500)
  expect_true(anyNA(d$Y))
  for (est in c("ATE", "ATT", "ATC", "ATO")) {
    spec <- tilt_spec(est)
    out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)
    expect_false(is.null(out$pifit))
    expect_true(out$TMLE_converged, label = est)
    for (nm in c("MLE", "IPW", "AIPW", "TMLE", "TMLE_se")) {
      expect_true(is.finite(out[[nm]]), label = paste(est, nm))
    }
    expect_lt(abs(mean(out$tmle_fit$EIF)), 1e-4)   # score equations (nearly) solved
  }
})

test_that("input checks for internal entry point with missing outcomes", {
  set.seed(4)
  d <- sim_missing_dgp(100)
  spec <- tilt_spec("ATT")
  Y_arm <- d$Y; Y_arm[d$A == 0] <- NA
  expect_error(estimate_tau_h_internal_tilted(Y_arm, d$A, d$W, spec$eta, spec$eta_dot),
               "No observed outcomes")
})

test_that("internal TMLE with missing outcomes (MAR) is approximately unbiased", {
  skip_on_cran()
  set.seed(5)
  B <- 40
  for (est in c("ATT", "ATO")) {
    spec <- tilt_spec(est)
    truth <- true_tilted_missing(spec$eta)
    ests <- numeric(B)
    for (b in seq_len(B)) {
      d <- sim_missing_dgp(600)
      ests[b] <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)$TMLE
    }
    expect_lt(abs(mean(ests) - truth), 0.04, label = est)
  }
})

test_that("internal bootstrap works with NA in Y", {
  set.seed(6)
  d <- sim_missing_dgp(300)
  spec <- tilt_spec("ATT")
  b <- bootstrap_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, n_boot = 10)
  expect_true(is.finite(b$se))
  expect_equal(b$n_failed, 0L)
})
