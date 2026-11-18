## SuperLearner must work when the user has NOT attached it with library(SuperLearner):
## regression test for "object 'All' not found" (SuperLearner looked up its learners and the
## default screening function `All` from the package's environment, where they are not visible).
## Synthetic data only.

skip_on_cran()
skip_if_not_installed("SuperLearner")

## Make sure SuperLearner is NOT attached while these tests run (put it back afterwards)
if ("package:SuperLearner" %in% search()) {
  detach("package:SuperLearner", character.only = TRUE)
  withr::defer(library(SuperLearner), testthat::teardown_env())
}

sl_data <- function(n = 150) {
  set.seed(81)
  W <- data.frame(x1 = rnorm(n), x2 = rnorm(n), sex = factor(sample(c("F", "M"), n, TRUE)))
  A <- rbinom(n, 1, plogis(0.5 * W$x1))
  Y <- rbinom(n, 1, plogis(-0.3 + 0.5 * A + 0.8 * W$x1))
  list(Y = Y, A = A, W = W)
}

test_that("every learner the app offers works for the outcome and treatment models", {
  d <- sl_data(120)
  expect_false("package:SuperLearner" %in% search())
  for (lib in names(SL_PRESETS)) {
    pkg <- SL_PRESETS[[lib]]
    if (!is.na(pkg) && !requireNamespace(pkg, quietly = TRUE)) next
    q <- suppressWarnings(fit_Q(d$Y, d$A, d$W, method = "SL", SL.library = lib))
    g <- suppressWarnings(fit_g(d$A, d$W, method = "SL", SL.library = lib))
    expect_true(all(is.finite(q$Q1W)) && all(q$Q1W > 0 & q$Q1W < 1), label = paste(lib, "Q"))
    expect_true(all(is.finite(g$g1W)), label = paste(lib, "g"))
  }
})

test_that("all three entry points run with SuperLearner for Q and g", {
  d <- sl_data(150)
  lib <- c("SL.glm", "SL.mean")
  out_f <- estimate_tau_h_fixed(d$Y, d$A, d$W, h = rep(1, 150), n_boot = 3,
                                Q_method = "SL", Q_SL.library = lib, g_method = "SL", g_SL.library = lib)
  expect_true(is.finite(out_f$TMLE))

  spec <- tilt_spec("ATT")
  out_i <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                          Q_method = "SL", Q_SL.library = lib, g_method = "SL", g_SL.library = lib)
  expect_true(is.finite(out_i$TMLE))

  ## external: pooled covariates, the first 100 rows are internal
  set.seed(82)
  Ws <- rbind(d$W[1:100, ], transform(d$W[101:150, ], x1 = x1 + 0.3))
  D <- c(rep(1L, 100), rep(0L, 50))
  spec_e <- tilt_spec_external("ATE-O")
  out_e <- estimate_tau_h_external_tilted(d$Y[1:100], d$A[1:100], D, Ws, spec_e$eta, spec_e$eta_dot,
                                          Q_method = "SL", Q_SL.library = lib, g_method = "SL", g_SL.library = lib,
                                          s_method = "SL", s_SL.library = lib)
  expect_true(is.finite(out_e$TMLE))
})

test_that("SuperLearner also works for the missing-outcome model", {
  d <- sl_data(150)
  d$Y[1:20] <- NA
  spec <- tilt_spec("ATE")
  lib <- c("SL.glm", "SL.mean")
  out <- estimate_tau_h_internal_tilted(d$Y, d$A, d$W, spec$eta, spec$eta_dot,
                                        Q_method = "SL", Q_SL.library = lib,
                                        g_method = "SL", g_SL.library = lib,
                                        pi_method = "SL", pi_SL.library = lib)
  expect_true(is.finite(out$TMLE))
  expect_false(is.null(out$pifit))
})

test_that("the app's run step works with SuperLearner choices (the failing case from the app)", {
  set.seed(83)
  n <- 150
  df <- data.frame(y = sample(c("No", "Yes"), n, TRUE), a = sample(c("C", "T"), n, TRUE),
                   x1 = rnorm(n), x2 = rnorm(n), w = runif(n, 0.5, 2))
  p <- prepare_analysis_data(df, "fixed", "y", "a", c("x1", "x2"), "Yes", "T", h = "w")
  s <- build_analysis_spec(p, list(Q_method = "SL", Q_library = c("SL.glm", "SL.glmnet"),
                                   g_method = "SL", g_library = c("SL.glm", "SL.glmnet")))
  r <- run_analysis_main(p, s)
  expect_false(inherits(r, "run_error"), info = if (inherits(r, "run_error")) r$message else "")
  expect_true(is.finite(r$table$Estimate[r$table$Estimator == "TMLE"]))
})
