## Tests for the app's settings step: build_analysis_spec() and the fixed-h weight column.
## Synthetic data only.

prep_for <- function(setting, missing_Y = FALSE, with_h = FALSE) {
  set.seed(31)
  n <- 120
  df <- data.frame(y = sample(c("No", "Yes"), n, TRUE), a = sample(c("C", "T"), n, TRUE),
                   x1 = rnorm(n), x2 = rnorm(n), w = runif(n, 0.5, 2))
  if (missing_Y) df$y[1:8] <- NA
  if (setting == "external") {
    df$src <- rep(c("trial", "ext"), each = n / 2)
    df$y[df$src == "ext"] <- NA; df$a[df$src == "ext"] <- NA
  }
  prepare_analysis_data(df, setting, "y", "a", c("x1", "x2"), "Yes", "T",
                        d = if (setting == "external") "src", d_one = if (setting == "external") "trial",
                        h = if (with_h || setting == "fixed") "w")
}

test_that("fixed-h weight column: optional, validated, only for the fixed setting", {
  ## h(x) is required for the fixed setting, with a helpful message
  df0 <- data.frame(y = c("No", "Yes", "No", "Yes"), a = c("C", "T", "C", "T"), x = 1:4)
  expect_error(prepare_analysis_data(df0, "fixed", "y", "a", "x", "Yes", "T"), "choose the column that holds the tilting weights")
  p1 <- prep_for("fixed"); expect_equal(length(p1$h), 120); expect_true(all(p1$h > 0))
  expect_named(p1$summary$h, c("min", "mean", "max", "ess"))
  expect_lte(p1$summary$h$ess, 120)
  expect_null(prep_for("internal", with_h = TRUE)$h)          # ignored outside fixed
  df <- data.frame(y = c("No", "Yes", "No", "Yes"), a = c("C", "T", "C", "T"), x = 1:4, w = c(1, -1, 1, 1))
  expect_error(prepare_analysis_data(df, "fixed", "y", "a", "x", "Yes", "T", h = "w"), "no negative")
  df$w <- c("a", "b", "c", "d")
  expect_error(prepare_analysis_data(df, "fixed", "y", "a", "x", "Yes", "T", h = "w"), "must be numeric")
})

test_that("estimand_choices() lists the right estimands per method", {
  expect_equal(unname(estimand_choices("internal")), c("ATE", "ATT", "ATC", "ATO", "custom"))
  expect_equal(unname(estimand_choices("external")), c("ATE-I", "ATE-E", "ATE-O", "custom"))
  expect_equal(unname(estimand_choices("fixed")), "fixed")
})

test_that("internal plan: arguments match the package and the plan really runs", {
  p <- prep_for("internal")
  s <- build_analysis_spec(p, list(estimand = "ATT", g_lower = 0.05))
  expect_equal(s$fun, "estimate_tau_h_internal_tilted")
  expect_equal(s$args$gbounds, c(0.05, 0.95))
  expect_true(is.function(s$eta) && is.function(s$eta_dot))
  expect_true(any(grepl("ATT", s$description)))
  out <- do.call(estimate_tau_h_internal_tilted, c(list(p$Y, p$A, p$W, s$eta, s$eta_dot), s$args))
  expect_true(is.finite(out$TMLE))
})

test_that("external plan includes the selection model and runs", {
  p <- prep_for("external")
  s <- build_analysis_spec(p, list(estimand = "ATE-O", s_lower = 0.05))
  expect_equal(s$fun, "estimate_tau_h_external_tilted")
  expect_equal(s$args$sbounds, c(0.05, 0.95))
  expect_true(any(grepl("Selection model", s$description)))
  out <- do.call(estimate_tau_h_external_tilted, c(list(p$Y, p$A, p$D, p$Ws, s$eta, s$eta_dot), s$args))
  expect_true(is.finite(out$TMLE))
})

test_that("fixed plan runs, with the cheap built-in bootstrap unless the user asks for one", {
  p <- prep_for("fixed")
  s <- build_analysis_spec(p, list())
  expect_equal(s$fun, "estimate_tau_h_fixed")
  expect_equal(s$args$n_boot, 2L)
  out <- do.call(estimate_tau_h_fixed, c(list(p$Y, p$A, p$W, p$h), s$args))
  expect_true(is.finite(out$TMLE))
  s2 <- build_analysis_spec(p, list(bootstrap = TRUE, n_boot = 50))
  expect_equal(s2$args$n_boot, 50)
  expect_true(any(grepl("target step", s2$warnings)))
})

test_that("missing outcomes add the missing-outcome model to the plan", {
  p <- prep_for("internal", missing_Y = TRUE)
  s <- build_analysis_spec(p, list(estimand = "ATE", pi_lower = 0.1))
  expect_equal(s$args$pibounds, c(0.1, 1))
  expect_true(any(grepl("Missing-outcome model", s$description)))
  out <- do.call(estimate_tau_h_internal_tilted, c(list(p$Y, p$A, p$W, s$eta, s$eta_dot), s$args))
  expect_true(is.finite(out$TMLE))
  ## no missing outcomes: no pi arguments, no pi line
  s0 <- build_analysis_spec(prep_for("internal"), list())
  expect_null(s0$args$pi_method)
})

test_that("SuperLearner choices map to arguments and are validated", {
  skip_if_not_installed("SuperLearner")
  p <- prep_for("internal")
  s <- build_analysis_spec(p, list(Q_method = "SL", Q_library = c("SL.glm", "SL.mean")))
  expect_equal(s$args$Q_method, "SL"); expect_equal(s$args$Q_SL.library, c("SL.glm", "SL.mean"))
  expect_true(any(grepl("SuperLearner \\(glm, mean\\)", s$description)))
  expect_error(build_analysis_spec(p, list(Q_method = "SL", Q_library = character())), "at least one")
  expect_error(build_analysis_spec(p, list(Q_method = "SL", Q_library = "SL.nope")), "Unknown learner")
})

test_that("invalid settings give readable errors", {
  p <- prep_for("internal")
  expect_error(build_analysis_spec(p, list(estimand = "ATE-O")), "not available")
  expect_error(build_analysis_spec(p, list(g_lower = 0.6)), "lower bound")
  expect_error(build_analysis_spec(p, list(conf_level = 1.2)), "Confidence level")
  expect_error(build_analysis_spec(p, list(bootstrap = TRUE, n_boot = 5)), "between 20 and 5000")
  expect_error(build_analysis_spec(p, list(max_iter = 0)), "at least 1")
  expect_error(build_analysis_spec(p, list(tol = 0)), "positive")
})

test_that("hints appear for ATT and ATE-E without the bootstrap, and disappear with it", {
  p <- prep_for("internal")
  expect_true(any(grepl("consider turning on the bootstrap", build_analysis_spec(p, list(estimand = "ATT"))$warnings)))
  expect_length(build_analysis_spec(p, list(estimand = "ATT", bootstrap = TRUE))$warnings, 0)
  expect_length(build_analysis_spec(p, list(estimand = "ATE"))$warnings, 0)
  pe <- prep_for("external")
  expect_true(any(grepl("ATE-E", build_analysis_spec(pe, list(estimand = "ATE-E"))$warnings)))
})
