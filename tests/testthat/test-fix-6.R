## Tests for review item 6: Q_args / g_args / s_args replace the shared `...`.

test_that("the three entry points no longer accept a stray `...`", {
  set.seed(1)
  d <- sim_internal_dgp(100)
  spec <- tilt_spec("ATT")
  ## an unknown argument is now an error instead of being silently forwarded
  expect_error(estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                              cvControl = list(V = 5)),
               "unused argument")
  expect_error(estimate_tau_h_fixed(d$Y, d$A, d$W, h = rep(1, 100), cvControl = list(V = 5)),
               "unused argument")
})

test_that("Q_args / g_args reach only their own fit (glm weights example)", {
  set.seed(2)
  d <- sim_internal_dgp(200)
  spec <- tilt_spec("ATT")

  ## glm() accepts `weights`; sending them only to Q must change Q but not g
  w <- sample(1:3, 200, replace = TRUE)
  base <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)
  wq   <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                         Q_args = list(weights = w))
  expect_equal(wq$gfit$g1W, base$gfit$g1W)          # g untouched
  expect_false(isTRUE(all.equal(wq$Qfit$Q1W, base$Qfit$Q1W)))  # Q changed

  wg <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                       g_args = list(weights = w))
  expect_equal(wg$Qfit$Q1W, base$Qfit$Q1W)          # Q untouched
  expect_false(isTRUE(all.equal(wg$gfit$g1W, base$gfit$g1W)))  # g changed
})

test_that("external entry point routes Q_args / g_args / s_args separately", {
  set.seed(3)
  d <- sim_external_dgp(100, 200)
  spec <- tilt_spec_external("ATE-O")
  n <- length(d$Y); N <- length(d$D)
  wn <- sample(1:3, n, replace = TRUE); wN <- sample(1:3, N, replace = TRUE)

  base <- estimate_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot)
  ws <- estimate_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot,
                                       s_args = list(weights = wN))
  ## only s changes
  expect_equal(ws$Qfit$Q1W, base$Qfit$Q1W)
  expect_equal(ws$gfit$g1W, base$gfit$g1W)
  expect_false(isTRUE(all.equal(ws$sfit$g1W, base$sfit$g1W)))

  wq <- estimate_tau_h_external_tilted(d$Y, d$A, d$D, d$Ws, spec$eta, spec$eta_dot,
                                       Q_args = list(weights = wn), g_args = list(weights = wn))
  expect_equal(wq$sfit$g1W, base$sfit$g1W)
})

test_that("default empty arg lists reproduce the previous results", {
  set.seed(4)
  d <- sim_internal_dgp(150)
  spec <- tilt_spec("ATO")
  a <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot)
  b <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                      Q_args = list(), g_args = list())
  expect_equal(a$TMLE, b$TMLE)
})
