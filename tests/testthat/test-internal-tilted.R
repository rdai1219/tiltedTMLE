test_that("tilt_spec() presets have the expected eta/eta_dot values", {
  g <- c(0.2, 0.5, 0.8)

  s_ate <- tilt_spec("ATE")
  expect_equal(s_ate$eta(g), rep(1, 3))
  expect_equal(s_ate$eta_dot(g), rep(0, 3))

  s_att <- tilt_spec("ATT")
  expect_equal(s_att$eta(g), g)
  expect_equal(s_att$eta_dot(g), rep(1, 3))

  s_atc <- tilt_spec("ATC")
  expect_equal(s_atc$eta(g), 1 - g)
  expect_equal(s_atc$eta_dot(g), rep(-1, 3))

  s_ato <- tilt_spec("ATO")
  expect_equal(s_ato$eta(g), g * (1 - g))
  expect_equal(s_ato$eta_dot(g), 1 - 2 * g)
})

test_that("ATE converges in a single iteration (eta_dot = 0 degenerate-fit guard)", {
  set.seed(10)
  d <- sim_internal_dgp(400)
  spec <- tilt_spec("ATE")
  out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, eta = spec$eta, eta_dot = spec$eta_dot,
                                Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
                                g_method = "glm", gform = A ~ X1 + X2 + X3,
                                truncate = FALSE)
  expect_equal(out$TMLE_iter, 1L)
})

test_that("ATT/ATC/ATO genuinely iterate (eta_dot != 0)", {
  set.seed(11)
  d <- sim_internal_dgp(400)
  for (est in c("ATT", "ATC", "ATO")) {
    spec <- tilt_spec(est)
    out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, eta = spec$eta, eta_dot = spec$eta_dot,
                                  Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
                                  g_method = "glm", gform = A ~ X1 + X2 + X3,
                                  truncate = FALSE)
    expect_gt(out$TMLE_iter, 1L)
  }
})

test_that("estimate_tau_h_internal_tilted() point estimates are unbiased on average under correct specification", {
  skip_on_cran()
  set.seed(12)
  B <- 60
  n <- 400
  for (est in c("ATE", "ATT", "ATC", "ATO")) {
    spec <- tilt_spec(est)
    true_tau <- approximate_true_tau_internal(spec$eta)

    ests <- numeric(B)
    for (b in seq_len(B)) {
      d <- sim_internal_dgp(n)
      out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, eta = spec$eta, eta_dot = spec$eta_dot,
                                    Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
                                    g_method = "glm", gform = A ~ X1 + X2 + X3,
                                    truncate = FALSE)
      ests[b] <- out$TMLE
    }
    bias <- mean(ests) - true_tau
    ## Loose bound: with B = 60 replications this is a sanity check,
    ## not a precise bias estimate -- just confirms no gross error
    ## (e.g. a sign flip or wrong nuisance plugged in) survived porting
    ## into the package.
    expect_lt(abs(bias), 0.05, label = paste(est, "bias"))
  }
})
