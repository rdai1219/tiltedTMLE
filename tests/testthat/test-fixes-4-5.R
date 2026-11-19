## Tests for review items 4-5 (fixed-h warning, input validation).
## Uses helper-dgp.R; run against the R_v2/ copies (incl. validate.R).

## ---- Item 4: warning when h is a function of the fitted g ------------------

test_that("estimate_tau_h_fixed() warns when h is a function of fitted g", {
  set.seed(1)
  d <- sim_internal_dgp(300)
  expect_warning(
    estimate_tau_h_fixed(d$Y, d$A, d$W, h = function(g1W, g0W) g1W, n_boot = 5),
    "ignores the eta_dot term"
  )
})

test_that("warn_fn_h = FALSE silences the warning, and a numeric h never warns", {
  set.seed(2)
  d <- sim_internal_dgp(300)
  expect_no_warning(
    estimate_tau_h_fixed(d$Y, d$A, d$W, h = function(g1W, g0W) g1W, n_boot = 5,
                         warn_fn_h = FALSE)
  )
  expect_no_warning(
    estimate_tau_h_fixed(d$Y, d$A, d$W, h = rep(1, length(d$Y)), n_boot = 5)
  )
})

## ---- Item 5: input validation ----------------------------------------------

test_that("non-numeric/infinite Y and missing A are rejected by all entry points", {
  set.seed(3)
  d <- sim_internal_dgp(100)
  spec <- tilt_spec("ATT")

  ## Y may be continuous since 0.2.0, so a value of 2 is now valid; Inf is not
  Y_bad <- d$Y; Y_bad[1] <- Inf
  A_na  <- d$A; A_na[1] <- NA

  expect_error(estimate_tau_h_fixed(Y_bad, d$A, d$W, h = rep(1, 100)), "`Y` must be")
  expect_error(estimate_tau_h_fixed(d$Y, A_na, d$W, h = rep(1, 100)), "`A` must be")
  expect_error(estimate_tau_h_internal_tilted(Y_bad, d$A, d$W, spec$eta, spec$eta_dot),
               "`Y` must be")
})

test_that("length mismatches and NA covariates are rejected", {
  set.seed(4)
  d <- sim_internal_dgp(100)
  spec <- tilt_spec("ATT")

  expect_error(estimate_tau_h_internal_tilted(d$Y[-1], d$A, d$W, spec$eta, spec$eta_dot),
               "same number of rows")

  W_na <- d$W; W_na$X1[3] <- NA
  expect_error(estimate_tau_h_internal_tilted(d$Y, d$A, W_na, spec$eta, spec$eta_dot),
               "missing values")
})

test_that("invalid gbounds are rejected", {
  set.seed(5)
  d <- sim_internal_dgp(100)
  spec <- tilt_spec("ATT")
  expect_error(estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                              gbounds = c(0.9, 0.1)),
               "gbounds")
  expect_error(estimate_tau_h_fixed(d$Y, d$A, d$W, h = rep(1, 100), gbounds = c(-0.1, 1)),
               "gbounds")
})

test_that("external entry point validates D, lengths, Wg and bounds", {
  set.seed(6)
  d <- sim_external_dgp(100, 200)
  spec <- tilt_spec_external("ATE-O")
  run <- function(...) estimate_tau_h_external_tilted(..., eta = spec$eta, eta_dot = spec$eta_dot)

  ## D not binary
  D_bad <- d$D; D_bad[1] <- 3L
  expect_error(run(d$Y, d$A, D_bad, d$Ws), "`D` must be")

  ## Y not internal-only (length != sum(D == 1))
  expect_error(run(c(d$Y, 0), c(d$A, 0), d$D, d$Ws), "internal-only")

  ## Ws rows != length(D)
  expect_error(run(d$Y, d$A, d$D, d$Ws[-1, ]), "one row per entry")

  ## Wg with wrong number of rows, or columns not in Ws
  expect_error(run(d$Y, d$A, d$D, d$Ws, Wg = d$Ws), "internal-only")
  Wg_bad <- d$Ws[which(d$D == 1), , drop = FALSE]; names(Wg_bad)[1] <- "Z"
  expect_error(run(d$Y, d$A, d$D, d$Ws, Wg = Wg_bad), "columns of `Ws`")

  ## bad sbounds
  expect_error(run(d$Y, d$A, d$D, d$Ws, sbounds = c(0.5, 0.5)), "sbounds")
})

test_that("valid inputs still run through the validated entry points", {
  set.seed(7)
  d <- sim_external_dgp(100, 200)
  spec <- tilt_spec_external("ATE-O")
  out <- estimate_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot)
  expect_true(is.finite(out$TMLE))

  di <- sim_internal_dgp(150)
  out2 <- estimate_tau_h_internal_tilted(di$Y, di$A, di$W,
                                         tilt_spec("ATO")$eta, tilt_spec("ATO")$eta_dot)
  expect_true(is.finite(out2$TMLE))
})
