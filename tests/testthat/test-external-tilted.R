test_that("tilt_spec_external() presets have the expected eta/eta_dot values", {
  s <- c(0.2, 0.5, 0.8)

  spec_i <- tilt_spec_external("ATE-I")
  expect_equal(spec_i$eta(s), rep(1, 3))
  expect_equal(spec_i$eta_dot(s), rep(0, 3))

  spec_e <- tilt_spec_external("ATE-E")
  expect_equal(spec_e$eta(s), (1 - s) / s)
  expect_equal(spec_e$eta_dot(s), -1 / s^2)

  spec_o <- tilt_spec_external("ATE-O")
  expect_equal(spec_o$eta(s), 1 - s)
  expect_equal(spec_o$eta_dot(s), rep(-1, 3))
})

test_that("ATE-I exactly matches the plain internal-only TMLE, ignoring external data", {
  ## Algebraic identity: with eta(s) = 1, theta_2 doesn't depend on s at
  ## all, so ATE-I should be bit-for-bit identical to running the plain
  ## internal-data-only ATE TMLE on the internal rows alone.
  set.seed(20)
  d <- sim_external_dgp(400, 800)
  idx1 <- which(d$D == 1)
  W_int <- d$Ws[idx1, , drop = FALSE]

  spec_ext <- tilt_spec_external("ATE-I")
  out_ext <- estimate_tau_h_external_tilted(
    d$Y, d$A, d$D, d$Ws, eta = spec_ext$eta, eta_dot = spec_ext$eta_dot,
    Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
    g_method = "glm", gform = A ~ X1 + X2 + X3,
    s_method = "glm", sform = A ~ X1 + X2 + X3,
    truncate_g = TRUE, gbounds = c(0.025, 0.975),
    truncate_s = TRUE, sbounds = c(0.025, 0.975)
  )

  out_int <- estimate_tau_h_fixed(
    d$Y, d$A, W_int, h = rep(1, length(d$Y)),
    Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
    g_method = "glm", gform = A ~ X1 + X2 + X3,
    truncate = TRUE, gbounds = c(0.025, 0.975), n_boot = 1
  )

  expect_equal(out_ext$TMLE, out_int$TMLE)
  expect_equal(out_ext$TMLE_se, out_int$TMLE_se)
})

test_that("sbounds re-clamping prevents the known ATE-E overflow failure mode", {
  ## Regression test for the specific failure this package's TMLE
  ## encountered during development: without re-clamping s after every
  ## fluctuation update, ATE-E's 1/s term can push Q to exactly 0/1 and
  ## crash fit_eps_h()'s glm() with "NA/NaN/Inf in 'y'". Loop across a
  ## handful of seeds since the failure was seed-dependent.
  for (seed in 100:104) {
    set.seed(seed)
    d <- sim_external_dgp(400, 800)
    spec <- tilt_spec_external("ATE-E")
    expect_no_error({
      out <- estimate_tau_h_external_tilted(
        d$Y, d$A, d$D, d$Ws, eta = spec$eta, eta_dot = spec$eta_dot,
        Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
        g_method = "glm", gform = A ~ X1 + X3,
        s_method = "glm", sform = A ~ X1 + X3,
        truncate_g = TRUE, gbounds = c(0.025, 0.975),
        truncate_s = TRUE, sbounds = c(0.025, 0.975)
      )
    })
  }
})

test_that("estimate_tau_h_external_tilted() point estimates are unbiased on average under correct specification", {
  skip_on_cran()
  set.seed(21)
  B <- 40
  for (est in c("ATE-I", "ATE-E", "ATE-O")) {
    spec <- tilt_spec_external(est)
    true_tau <- approximate_true_tau_external(spec$eta)

    ests <- numeric(B)
    for (b in seq_len(B)) {
      d <- sim_external_dgp(400, 800)
      out <- estimate_tau_h_external_tilted(
        d$Y, d$A, d$D, d$Ws, eta = spec$eta, eta_dot = spec$eta_dot,
        Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
        g_method = "glm", gform = A ~ X1 + X2 + X3,
        s_method = "glm", sform = A ~ X1 + X2 + X3,
        truncate_g = TRUE, gbounds = c(0.025, 0.975),
        truncate_s = TRUE, sbounds = c(0.025, 0.975)
      )
      ests[b] <- out$TMLE
    }
    bias <- mean(ests) - true_tau
    ## ATE-E is known to have much higher variance (1/s-term
    ## instability), so give it a looser tolerance than ATE-I/ATE-O.
    tol <- if (est == "ATE-E") 0.08 else 0.04
    expect_lt(abs(bias), tol, label = paste(est, "bias"))
  }
})
