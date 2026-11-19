## Tests for cross-fitting (cv_folds, added in 0.3.0).
## Uses the DGPs in helper-dgp.R and helper-missing-dgp.R.

clamp <- function(x, b = c(0.025, 0.975)) pmin(pmax(x, b[1]), b[2])

## ---- Default: no cross-fitting ------------------------------------------------

test_that("cv_folds = 1 (the default) means no cross-fitting", {
  set.seed(1)
  d <- sim_internal_dgp(300)
  spec <- tilt_spec("ATT")
  out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)
  out1 <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, cv_folds = 1)
  expect_null(out$folds)
  expect_true("folds" %in% names(out))
  expect_equal(out1$TMLE, out$TMLE)
  expect_s3_class(out$Qfit$fit, "glm")
})

## ---- The nuisance predictions really are out-of-fold --------------------------

test_that("internal: every row is predicted by the model that did not see it", {
  set.seed(2)
  d <- sim_internal_dgp(300)
  spec <- tilt_spec("ATO")
  out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, cv_folds = 5)

  expect_length(out$folds, 300)
  expect_setequal(unique(out$folds), 1:5)
  expect_length(out$Qfit$fit, 5)
  expect_length(out$gfit$fit, 5)
  ## folds are balanced within each treatment arm
  counts <- table(out$folds, d$A)
  expect_true(all(apply(counts, 2, function(x) max(x) - min(x)) <= 1))

  for (v in 1:5) {
    te <- which(out$folds == v)
    ## the fold-v model was fit without the rows of fold v ...
    n_train <- sum(!is.na(d$Y[out$folds != v]))
    expect_equal(stats::nobs(out$Qfit$fit[[v]]), n_train)
    ## ... and supplies the predictions for exactly those rows
    expect_equal(out$Qfit$Q1W[te], predict_Qhat(out$Qfit$fit[[v]], 1, d$W[te, ]))
    expect_equal(out$Qfit$Q0W[te], predict_Qhat(out$Qfit$fit[[v]], 0, d$W[te, ]))
    expect_equal(out$gfit$g1W[te], clamp(predict_ghat(out$gfit$fit[[v]], d$W[te, ])))
  }
})

test_that("external: Q, g and s are out-of-fold over the pooled rows", {
  set.seed(3)
  e <- sim_external_dgp(300, 600)
  spec <- tilt_spec_external("ATE-O")
  out <- estimate_tau_h_external_tilted(e$Y, e$A, e$D, e$Ws, spec$eta, spec$eta_dot, cv_folds = 4)
  idx1 <- which(e$D == 1)

  expect_length(out$folds, 900)
  ## both internal and external rows are spread over all folds
  expect_setequal(unique(out$folds[e$D == 1]), 1:4)
  expect_setequal(unique(out$folds[e$D == 0]), 1:4)

  for (v in 1:4) {
    te <- which(out$folds == v)
    expect_equal(out$Qfit$Q1W[te], predict_Qhat(out$Qfit$fit[[v]], 1, e$Ws[te, ]))
    expect_equal(out$sfit$g1W[te], clamp(predict_ghat(out$sfit$fit[[v]], e$Ws[te, ])))
    te_int <- which(out$folds[idx1] == v)      # internal rows of fold v
    expect_equal(out$gfit$g1W[te_int],
                 clamp(predict_ghat(out$gfit$fit[[v]], e$Ws[idx1[te_int], ])))
  }
})

## ---- Estimates -----------------------------------------------------------------

test_that("cross-fitted TMLE is close to the truth and to the full-data TMLE", {
  set.seed(4)
  d <- sim_internal_dgp(1500)
  spec <- tilt_spec("ATT")
  full <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)
  cv <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, cv_folds = 5)
  truth <- approximate_true_tau_internal(spec$eta)
  expect_true(cv$TMLE_converged)
  expect_lt(abs(cv$TMLE - truth), 4 * cv$TMLE_se)
  ## with correctly specified glm models, cross-fitting changes little
  expect_lt(abs(cv$TMLE - full$TMLE), full$TMLE_se)

  set.seed(5)
  e <- sim_external_dgp(600, 1200)
  for (est in c("ATE-I", "ATE-O")) {
    sp <- tilt_spec_external(est)
    cv <- estimate_tau_h_external_tilted(e$Y, e$A, e$D, e$Ws, sp$eta, sp$eta_dot, cv_folds = 5)
    truth <- approximate_true_tau_external(sp$eta)
    expect_true(is.finite(cv$TMLE))
    expect_lt(abs(cv$TMLE - truth), 4 * cv$TMLE_se)
  }
})

test_that("results are reproducible with set.seed()", {
  d <- sim_internal_dgp(200)
  spec <- tilt_spec("ATO")
  set.seed(6); a <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, cv_folds = 3)
  set.seed(6); b <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot, cv_folds = 3)
  expect_identical(a$folds, b$folds)
  expect_equal(a$TMLE, b$TMLE)
})

## ---- Works with the other features ---------------------------------------------

test_that("cross-fitting works with fixed h, missing and continuous outcomes, and the bootstrap", {
  set.seed(7)
  d <- sim_internal_dgp(300)
  out <- estimate_tau_h_fixed(d$Y, d$A, d$W, h = rep(1, 300), n_boot = 10, cv_folds = 5)
  expect_true(is.finite(out$TMLE) && is.finite(out$Bootstrap_SE))
  expect_length(out$folds, 300)

  ## missing outcomes: pi is cross-fitted too
  m <- sim_missing_dgp(400)
  spec <- tilt_spec("ATT")
  out <- estimate_tau_h_internal_tilted(m$Y, m$A, m$W, spec$eta, spec$eta_dot, cv_folds = 5)
  expect_true(is.finite(out$TMLE))
  expect_length(out$pifit$fit, 5)

  ## continuous outcome
  Yc <- 1 + d$W$X1 + d$A * (1 + 0.5 * d$W$X1) + stats::rnorm(300)
  out <- estimate_tau_h_internal_tilted(Yc, d$A, d$W, spec$eta, spec$eta_dot, cv_folds = 5)
  expect_true(is.finite(out$TMLE))
  expect_equal(out$Ybounds, range(Yc))

  ## external with missing outcomes
  e <- sim_external_dgp(300, 300)
  Ye <- e$Y
  Ye[sample(300, 30)] <- NA
  sp <- tilt_spec_external("ATE-O")
  out <- estimate_tau_h_external_tilted(Ye, e$A, e$D, e$Ws, sp$eta, sp$eta_dot, cv_folds = 3)
  expect_true(is.finite(out$TMLE))
  expect_length(out$pifit$pi1W_all, 600)

  ## the full-refit bootstrap passes cv_folds on (new folds on every resample)
  b <- bootstrap_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                       n_boot = 5, cv_folds = 3)
  expect_true(is.finite(b$se))
})

test_that("cross-fitting works with SuperLearner", {
  skip_if_not_installed("SuperLearner")
  set.seed(8)
  d <- sim_internal_dgp(300)
  spec <- tilt_spec("ATE")
  out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                        Q_method = "SL", Q_SL.library = c("SL.glm", "SL.mean"),
                                        g_method = "SL", g_SL.library = c("SL.glm", "SL.mean"),
                                        cv_folds = 3)
  expect_true(is.finite(out$TMLE))
  expect_s3_class(out$Qfit$fit[[1]], "SuperLearner")
})

## ---- Input checks ------------------------------------------------------------------

test_that("invalid cv_folds is rejected", {
  set.seed(9)
  d <- sim_internal_dgp(100)
  spec <- tilt_spec("ATE")
  for (bad in list(0, 2.5, c(2, 3), 101, NA)) {
    expect_error(estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                                cv_folds = bad), "`cv_folds` must be")
  }
  expect_error(estimate_tau_h_fixed(d$Y, d$A, d$W, h = rep(1, 100), cv_folds = 0),
               "`cv_folds` must be")
})
