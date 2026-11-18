## Tests for missing outcomes in estimate_tau_h_fixed() (NA in Y, MAR given A and W).

test_that("fit_pi() returns bounded probabilities for A = 1 and A = 0", {
  set.seed(1)
  d <- sim_missing_dgp(300)
  Delta <- as.integer(!is.na(d$Y))
  p <- fit_pi(Delta, d$A, d$W, pibounds = c(0.1, 0.95))
  expect_length(p$pi1W, 300)
  expect_true(all(p$pi1W >= 0.1 & p$pi1W <= 0.95))
  expect_true(all(p$pi0W >= 0.1 & p$pi0W <= 0.95))
  expect_equal(p$piAW, ifelse(d$A == 1, p$pi1W, p$pi0W))
})

test_that("with no missing outcomes, results are unchanged by the pi machinery", {
  set.seed(2)
  d <- sim_missing_dgp(300, missing = FALSE)
  h <- rep(1, 300)
  new <- estimate_tau_h_fixed(d$Y, d$A, d$W, h = h, n_boot = 5)
  expect_null(new$pifit)

  ## manual pipeline = the original (complete-data) code path
  Qfit <- fit_Q(d$Y, d$A, d$W)
  gfit <- fit_g(d$A, d$W)
  old <- estimate_tau_h_fixed_given_fits(d$Y, d$A, h, Qfit, gfit, n_boot = 5)
  expect_equal(new$TMLE, old$TMLE)
  expect_equal(new$TMLE_se, old$TMLE_se)
  expect_equal(new$AIPW, old$AIPW)
  expect_equal(new$IPW, old$IPW)
})

test_that("explicit pi = 1 and Delta = 1 reproduce tmle_h() defaults", {
  set.seed(3)
  d <- sim_missing_dgp(200, missing = FALSE)
  Qfit <- fit_Q(d$Y, d$A, d$W); gfit <- fit_g(d$A, d$W); h <- rep(1, 200)
  a <- tmle_h(d$Y, d$A, h, Qfit$QAW, Qfit$Q1W, Qfit$Q0W, gfit$g1W)
  b <- tmle_h(d$Y, d$A, h, Qfit$QAW, Qfit$Q1W, Qfit$Q0W, gfit$g1W,
              Delta = rep(1L, 200), pi1W = rep(1, 200), pi0W = rep(1, 200))
  expect_equal(a$estimate, b$estimate)
  expect_equal(a$EIF, b$EIF)
})

test_that("estimate_tau_h_fixed() runs with NA in Y and returns finite results", {
  set.seed(4)
  d <- sim_missing_dgp(400)
  expect_true(anyNA(d$Y))
  out <- estimate_tau_h_fixed(d$Y, d$A, d$W, h = rep(1, 400), n_boot = 20)
  expect_false(is.null(out$pifit))
  for (nm in c("MLE", "IPW", "AIPW", "TMLE", "TMLE_se", "Bootstrap_SE")) {
    expect_true(is.finite(out[[nm]]), label = nm)
  }
  ## the efficient influence function has mean zero at the TMLE (score equation solved)
  expect_lt(abs(mean(out$tmle_fit$EIF)), 1e-6)
})

test_that("input checks for missing outcomes", {
  set.seed(5)
  d <- sim_missing_dgp(100)

  ## observed Y values must still be 0/1
  Y_bad <- d$Y; Y_bad[which(!is.na(Y_bad))[1]] <- 2
  expect_error(estimate_tau_h_fixed(Y_bad, d$A, d$W, h = rep(1, 100)), "must be a 0/1")

  ## an arm with no observed outcome is an error
  Y_arm <- d$Y; Y_arm[d$A == 1] <- NA
  expect_error(estimate_tau_h_fixed(Y_arm, d$A, d$W, h = rep(1, 100)), "No observed outcomes")

  ## A must still be complete
  A_na <- d$A; A_na[1] <- NA
  expect_error(estimate_tau_h_fixed(d$Y, A_na, d$W, h = rep(1, 100)), "`A` must be")
})

test_that("TMLE with missing outcomes (MAR) is approximately unbiased and beats the naive estimate", {
  skip_on_cran()
  truth <- true_ate_missing()
  set.seed(6)
  B <- 40
  tmle_est <- naive_est <- numeric(B)
  for (b in seq_len(B)) {
    d <- sim_missing_dgp(600)
    out <- estimate_tau_h_fixed(d$Y, d$A, d$W, h = rep(1, 600), n_boot = 1)
    tmle_est[b] <- out$TMLE
    obs <- !is.na(d$Y)  # naive complete-case difference in means
    naive_est[b] <- mean(d$Y[obs & d$A == 1]) - mean(d$Y[obs & d$A == 0])
  }
  expect_lt(abs(mean(tmle_est) - truth), 0.03)
  expect_lt(abs(mean(tmle_est) - truth), abs(mean(naive_est) - truth))
})
