## Tests for missing outcomes in estimate_tau_h_external_tilted() (NA in Y, MAR given A and X).
## Uses sim_external_dgp() / approximate_true_tau_external() from helper-dgp.R.

## External data with outcomes missing at random among the INTERNAL rows:
## P(observed) depends on X1 and A (Y, A are internal-only, so Y gets NAs).
sim_external_missing <- function(n = 300, m = 600) {
  d <- sim_external_dgp(n, m)
  X1_int <- d$Ws$X1[which(d$D == 1)]
  obs <- stats::rbinom(n, 1, stats::plogis(0.8 + 0.8 * X1_int - 0.5 * d$A))
  d$Y[obs == 0] <- NA
  d
}

test_that("with no missing outcomes, external results match the manual (complete-data) pipeline", {
  set.seed(1)
  d <- sim_external_dgp(150, 300)
  idx1 <- which(d$D == 1)
  spec <- tilt_spec_external("ATE-O")
  new <- estimate_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot)
  expect_null(new$pifit)

  Ws <- d$Ws
  Qfit_int <- fit_Q(d$Y, d$A, Ws[idx1, , drop = FALSE])
  Qfit <- list(QAW = Qfit_int$QAW, Q1W = predict_Qhat(Qfit_int$fit, 1, Ws),
               Q0W = predict_Qhat(Qfit_int$fit, 0, Ws), fit = Qfit_int$fit, method = "glm")
  gfit <- fit_g(d$A, Ws[idx1, , drop = FALSE])
  sfit <- fit_g(d$D, Ws)
  g1W_all <- pmin(pmax(predict_ghat(gfit$fit, Ws), 0.025), 0.975)
  old <- estimate_tau_h_external_tilted_given_fits(d$Y, d$A, d$D, idx1, Qfit, gfit, sfit,
                                                   spec$eta, spec$eta_dot, g1W_all = g1W_all)
  expect_equal(new$TMLE, old$TMLE)
  expect_equal(new$TMLE_se, old$TMLE_se)
  expect_equal(new$IPW, old$IPW)
})

test_that("eif_external_tilted_h() with Delta = 1 and pi = 1 equals the default", {
  set.seed(2)
  d <- sim_external_dgp(100, 200)
  idx1 <- which(d$D == 1); n <- 100
  Ws <- d$Ws
  Qfit_int <- fit_Q(d$Y, d$A, Ws[idx1, , drop = FALSE])
  Q1W <- predict_Qhat(Qfit_int$fit, 1, Ws); Q0W <- predict_Qhat(Qfit_int$fit, 0, Ws)
  gfit <- fit_g(d$A, Ws[idx1, , drop = FALSE]); sfit <- fit_g(d$D, Ws)
  spec <- tilt_spec_external("ATE-O")
  eta_s <- spec$eta(sfit$g1W)
  tau <- sum(eta_s[idx1] * (Q1W[idx1] - Q0W[idx1])) / sum(eta_s[idx1])
  a <- eif_external_tilted_h(d$D, d$A, d$Y, idx1, eta_s, spec$eta_dot(sfit$g1W),
                             Q1W, Q0W, gfit$g1W, gfit$g0W, sfit$g1W, tau)
  b <- eif_external_tilted_h(d$D, d$A, d$Y, idx1, eta_s, spec$eta_dot(sfit$g1W),
                             Q1W, Q0W, gfit$g1W, gfit$g0W, sfit$g1W, tau,
                             Delta_int = rep(1L, n), pi1W_int = rep(1, n), pi0W_int = rep(1, n))
  expect_equal(a, b)
})

test_that("external tilted estimators run with NA in Y and converge", {
  set.seed(3)
  d <- sim_external_missing(300, 600)
  expect_true(anyNA(d$Y))
  for (est in c("ATE-I", "ATE-E", "ATE-O")) {
    spec <- tilt_spec_external(est)
    out <- estimate_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot)
    expect_false(is.null(out$pifit))
    expect_length(out$pifit$pi1W_all, 900)          # extrapolated to all N rows
    expect_true(out$TMLE_converged, label = est)
    for (nm in c("MLE", "IPW", "TMLE", "TMLE_se")) {
      expect_true(is.finite(out[[nm]]), label = paste(est, nm))
    }
    ## mean(EIF) is exactly 0 only for ATE-I; for ATE-E/ATE-O it is not forced to 0 by the
    ## algorithm (this is also true with complete data), so only check ATE-I.
    if (est == "ATE-I") expect_lt(abs(mean(out$tmle_fit$EIF)), 1e-6)
  }
})

test_that("input checks for external entry point with missing outcomes", {
  set.seed(4)
  d <- sim_external_missing(100, 200)
  spec <- tilt_spec_external("ATE-O")
  Y_arm <- d$Y; Y_arm[d$A == 0] <- NA
  expect_error(estimate_tau_h_external_tilted(Y_arm, d$A, d$D, d$Ws, spec$eta, spec$eta_dot),
               "No observed outcomes")
  Y_bad <- d$Y; Y_bad[which(!is.na(Y_bad))[1]] <- Inf
  expect_error(estimate_tau_h_external_tilted(Y_bad, d$A, d$D, d$Ws, spec$eta, spec$eta_dot),
               "`Y` must be")
})

test_that("external TMLE with missing outcomes (MAR) is approximately unbiased", {
  skip_on_cran()
  set.seed(5)
  B <- 30
  for (est in c("ATE-I", "ATE-O")) {
    spec <- tilt_spec_external(est)
    truth <- approximate_true_tau_external(spec$eta)
    ests <- numeric(B)
    for (b in seq_len(B)) {
      d <- sim_external_missing(300, 600)
      ests[b] <- estimate_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot)$TMLE
    }
    expect_lt(abs(mean(ests) - truth), 0.05, label = est)
  }
})

test_that("external bootstrap works with NA in Y", {
  set.seed(6)
  d <- sim_external_missing(150, 300)
  spec <- tilt_spec_external("ATE-O")
  b <- bootstrap_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot, n_boot = 8)
  expect_true(is.finite(b$se))
  expect_equal(b$n_failed, 0L)
})
