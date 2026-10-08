test_that("mle_h, ipw_h, aipw_h, tmle_h agree on a simple fixed-h case", {
  set.seed(1)
  d <- sim_internal_dgp(300)
  Qfit <- fit_Q(d$Y, d$A, d$W, method = "glm", Qform = Y ~ A * (X1 + X2 + X3))
  gfit <- fit_g(d$A, d$W, method = "glm", gform = A ~ X1 + X2 + X3, gbounds = c(0, 1))

  h <- rep(1, length(d$Y))
  mle_est <- mle_h(h, Qfit$Q1W, Qfit$Q0W)
  ipw_est <- ipw_h(d$Y, d$A, h, gfit$g1W)
  aipw_out <- aipw_h(d$Y, d$A, h, Qfit$Q1W, Qfit$Q0W, gfit$g1W)
  tmle_out <- tmle_h(d$Y, d$A, h, Qfit$QAW, Qfit$Q1W, Qfit$Q0W, gfit$g1W)

  expect_true(is.finite(mle_est))
  expect_true(is.finite(ipw_est))
  expect_true(is.finite(aipw_out$estimate) && aipw_out$se > 0)
  expect_true(is.finite(tmle_out$estimate) && tmle_out$se > 0)
  ## All four estimators should be in the same rough ballpark for a
  ## correctly-specified, well-behaved DGP.
  ests <- c(mle_est, ipw_est, aipw_out$estimate, tmle_out$estimate)
  expect_true(diff(range(ests)) < 0.2)
})

test_that("tmle_bootstrap_h returns a sane bootstrap distribution", {
  set.seed(2)
  d <- sim_internal_dgp(300)
  Qfit <- fit_Q(d$Y, d$A, d$W, method = "glm", Qform = Y ~ A * (X1 + X2 + X3))
  gfit <- fit_g(d$A, d$W, method = "glm", gform = A ~ X1 + X2 + X3, gbounds = c(0, 1))
  h <- rep(1, length(d$Y))
  tmle_out <- tmle_h(d$Y, d$A, h, Qfit$QAW, Qfit$Q1W, Qfit$Q0W, gfit$g1W)

  boot_out <- tmle_bootstrap_h(d$Y, d$A, h, Qfit$Q1W, Qfit$Q0W, gfit$g1W,
                                tmle_fit = tmle_out, n_boot = 50)

  expect_length(boot_out$estimates, 50)
  expect_true(all(is.finite(boot_out$estimates)))
  expect_true(boot_out$se > 0)
  expect_length(boot_out$CI_normal, 2)
  expect_lt(boot_out$CI_normal[1], boot_out$CI_normal[2])
})

test_that("estimate_tau_h_fixed() reproduces the same TMLE as calling tmle_h() directly", {
  set.seed(3)
  d <- sim_internal_dgp(300)
  out <- estimate_tau_h_fixed(d$Y, d$A, d$W, h = rep(1, length(d$Y)),
                               Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
                               g_method = "glm", gform = A ~ X1 + X2 + X3,
                               truncate = FALSE, n_boot = 20)

  Qfit <- fit_Q(d$Y, d$A, d$W, method = "glm", Qform = Y ~ A * (X1 + X2 + X3))
  gfit <- fit_g(d$A, d$W, method = "glm", gform = A ~ X1 + X2 + X3, gbounds = c(0, 1))
  direct <- tmle_h(d$Y, d$A, rep(1, length(d$Y)), Qfit$QAW, Qfit$Q1W, Qfit$Q0W, gfit$g1W)

  expect_equal(out$TMLE, direct$estimate)
  expect_equal(out$TMLE_se, direct$se)
})

test_that("estimate_tau_h_fixed() accepts h as a function of (g1W, g0W)", {
  set.seed(4)
  d <- sim_internal_dgp(300)
  out <- estimate_tau_h_fixed(d$Y, d$A, d$W, h = function(g1W, g0W) g1W,
                               Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
                               g_method = "glm", gform = A ~ X1 + X2 + X3,
                               truncate = FALSE, n_boot = 20)
  expect_true(is.finite(out$TMLE))
})
