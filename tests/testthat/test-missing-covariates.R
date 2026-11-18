## Tests for add_missing_indicators() (indicator method for missing covariates).

test_that("continuous column: mean imputation + indicator", {
  W <- data.frame(age = c(50, NA, 70, 60))
  out <- add_missing_indicators(W)
  expect_equal(out$age, c(50, 60, 70, 60))        # mean of 50, 70, 60 = 60
  expect_equal(out$age_miss, c(0L, 1L, 0L, 0L))
  expect_named(out, c("age", "age_miss"))
})

test_that("binary column: majority-category imputation + indicator", {
  W <- data.frame(sex = c(1, 1, 0, NA, 1))         # majority = 1
  out <- add_missing_indicators(W)
  expect_equal(out$sex, c(1, 1, 0, 1, 1))
  expect_equal(out$sex_miss, c(0L, 0L, 0L, 1L, 0L))

  ## logical and 2-level factor columns are treated as binary too
  W2 <- data.frame(flag = c(TRUE, FALSE, FALSE, NA),
                   grp = factor(c("a", "b", "b", NA)))
  out2 <- add_missing_indicators(W2)
  expect_equal(out2$flag, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(as.character(out2$grp), c("a", "b", "b", "b"))
})

test_that("complete columns are untouched and get no indicator", {
  W <- data.frame(x1 = c(1, 2, 3), x2 = c(4, NA, 6))
  out <- add_missing_indicators(W)
  expect_equal(out$x1, W$x1)
  expect_named(out, c("x1", "x2", "x2_miss"))

  ## no missing data at all: returned as is
  W_full <- data.frame(a = 1:3, b = c(0.5, 0.6, 0.7))
  expect_equal(add_missing_indicators(W_full), W_full)
})

test_that("the input is not modified, and matrices are accepted", {
  W <- data.frame(x = c(1, NA, 3))
  W_copy <- W
  invisible(add_missing_indicators(W))
  expect_identical(W, W_copy)

  M <- matrix(c(1, NA, 3, 4, 5, 6), ncol = 2, dimnames = list(NULL, c("a", "b")))
  out <- add_missing_indicators(M)
  expect_named(out, c("a", "b", "a_miss"))
})

test_that("custom suffix works; all-NA column and name clashes are errors", {
  W <- data.frame(x = c(1, NA, 3))
  expect_named(add_missing_indicators(W, suffix = ".na"), c("x", "x.na"))

  expect_error(add_missing_indicators(data.frame(x = c(NA_real_, NA_real_))),
               "entirely missing")

  W_clash <- data.frame(x = c(1, NA, 3), x_miss = c(0, 0, 0))
  expect_error(add_missing_indicators(W_clash), "already in `W`")
})

test_that("estimators reject NA covariates, but run after add_missing_indicators()", {
  set.seed(1)
  d <- sim_internal_dgp(200)
  W_na <- d$W
  W_na$X1[sample(200, 20)] <- NA
  W_na$X3[sample(200, 10)] <- NA

  expect_error(estimate_tau_h_fixed(d$Y, d$A, W_na, h = rep(1, 200)), "add_missing_indicators")

  W_aug <- add_missing_indicators(W_na)
  expect_named(W_aug, c("X1", "X2", "X3", "X1_miss", "X3_miss"))
  out <- estimate_tau_h_fixed(d$Y, d$A, W_aug, h = rep(1, 200), n_boot = 10)
  expect_true(is.finite(out$TMLE))

  spec <- tilt_spec("ATT")
  out2 <- estimate_tau_h_internal_tilted(d$Y, d$A, W_aug, spec$eta, spec$eta_dot)
  expect_true(is.finite(out2$TMLE))
})
