## Tests for the app's run step: run_analysis_main(), run_analysis_bootstrap() and the
## background wrapper. Synthetic data only.

prep_run <- function(setting, missing_Y = FALSE) {
  set.seed(51)
  n <- 160
  df <- data.frame(y = sample(c("No", "Yes"), n, TRUE), a = sample(c("C", "T"), n, TRUE),
                   x1 = rnorm(n), x2 = rnorm(n), w = runif(n, 0.5, 2))
  if (missing_Y) df$y[1:12] <- NA
  if (setting == "external") {
    df$src <- rep(c("trial", "ext"), each = n / 2)
    df$y[df$src == "ext"] <- NA; df$a[df$src == "ext"] <- NA
  }
  prepare_analysis_data(df, setting, "y", "a", c("x1", "x2"), "Yes", "T",
                        d = if (setting == "external") "src", d_one = if (setting == "external") "trial",
                        h = if (setting == "fixed") "w")
}

test_that("main fit: fixed, internal and external return tables and diagnostics", {
  for (setting in c("fixed", "internal", "external")) {
    p <- prep_run(setting)
    est <- switch(setting, fixed = "fixed", internal = "ATT", external = "ATE-O")
    s <- build_analysis_spec(p, list(estimand = est))
    r <- run_analysis_main(p, s)
    expect_false(inherits(r, "run_error"), label = setting)
    expect_true("TMLE" %in% r$table$Estimator, label = setting)
    expect_true(all(is.finite(r$table$Estimate)), label = setting)
    tm <- r$table[r$table$Estimator == "TMLE", ]
    expect_lt(tm$`CI lower`, tm$Estimate); expect_gt(tm$`CI upper`, tm$Estimate)
    expect_true(is.numeric(r$diagnostics$EIF), label = setting)
    expect_true(r$seconds >= 0)
    expect_match(r$package_version, "^[0-9]")
  }
})

test_that("main fit matches calling the package function directly", {
  p <- prep_run("internal")
  s <- build_analysis_spec(p, list(estimand = "ATO"))
  r <- run_analysis_main(p, s)
  spec_fns <- tilt_spec("ATO")
  direct <- do.call(estimate_tau_h_internal_tilted,
                    c(list(p$Y, p$A, p$W, spec_fns$eta, spec_fns$eta_dot), s$args))
  expect_equal(r$table$Estimate[r$table$Estimator == "TMLE"], direct$TMLE)
  expect_equal(r$table$SE[r$table$Estimator == "TMLE"], direct$TMLE_se)
})

test_that("confidence level is applied to the EIF interval", {
  p <- prep_run("internal")
  w90 <- run_analysis_main(p, build_analysis_spec(p, list(estimand = "ATT", conf_level = 0.90)))$table
  w99 <- run_analysis_main(p, build_analysis_spec(p, list(estimand = "ATT", conf_level = 0.99)))$table
  width <- function(t) { r <- t[t$Estimator == "TMLE", ]; r$`CI upper` - r$`CI lower` }
  expect_gt(width(w99), width(w90))
})

test_that("missing outcomes flow through the run, with diagnostics for pi", {
  p <- prep_run("internal", missing_Y = TRUE)
  r <- run_analysis_main(p, build_analysis_spec(p, list(estimand = "ATE")))
  expect_false(inherits(r, "run_error"))
  expect_length(r$diagnostics$pi1W, 160)
})

test_that("fixed method: bootstrap row appears only when requested", {
  p <- prep_run("fixed")
  r0 <- run_analysis_main(p, build_analysis_spec(p, list()))
  expect_false(any(grepl("bootstrap", r0$table$Estimator)))
  r1 <- run_analysis_main(p, build_analysis_spec(p, list(bootstrap = TRUE, n_boot = 30)))
  expect_true(any(grepl("target-step bootstrap", r1$table$Estimator)))
})

test_that("warnings are captured, and failures become a run_error (not a crash)", {
  p <- prep_run("internal")
  s <- build_analysis_spec(p, list(estimand = "ATT"))
  s_few <- s; s_few$args$max_iter <- 1; s_few$args$tol <- 1e-12     # forces non-convergence
  r <- run_analysis_main(p, s_few)
  expect_false(r$converged)
  expect_true(any(grepl("did not converge", r$warnings)))
  s_bad <- s; s_bad$args$Q_method <- "nope"
  e <- run_analysis_main(p, s_bad)
  expect_s3_class(e, "run_error"); expect_true(nzchar(e$message))
})

test_that("bootstrap stage: internal and external, seed makes it reproducible", {
  p <- prep_run("internal")
  s <- build_analysis_spec(p, list(estimand = "ATT", bootstrap = TRUE, n_boot = 20, seed = 7))
  b1 <- run_analysis_bootstrap(p, s); b2 <- run_analysis_bootstrap(p, s)
  expect_false(inherits(b1, "run_error"))
  expect_equal(b1$estimates, b2$estimates)                   # same seed, same result
  expect_equal(length(b1$estimates), 20)
  expect_true(is.finite(b1$se))

  pe <- prep_run("external")
  se <- build_analysis_spec(pe, list(estimand = "ATE-O", bootstrap = TRUE, n_boot = 20, seed = 3))
  be <- run_analysis_bootstrap(pe, se)
  expect_false(inherits(be, "run_error")); expect_true(is.finite(be$se))
})

test_that("background run: the worker process loads this package version and returns the bundle", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("pkgload")
  old_plan <- future::plan(future::multisession, workers = 1)
  on.exit(future::plan(old_plan), add = TRUE)
  p <- prep_run("internal")
  s <- build_analysis_spec(p, list(estimand = "ATT"))
  fut <- run_in_background("main", p, s)
  r <- future::value(fut)
  expect_false(inherits(r, "run_error"))
  local <- run_analysis_main(p, s)
  expect_equal(r$table$Estimate, local$table$Estimate)       # same code, same numbers
  expect_equal(r$package_version, as.character(utils::packageVersion("tiltedTMLE")))
})

test_that("a failed run reports where it failed (function names and message, never data)", {
  p <- prep_run("internal")
  s <- build_analysis_spec(p, list(estimand = "ATT"))
  s$args$Q_method <- "nope"
  e <- run_analysis_main(p, s)
  expect_s3_class(e, "run_error")
  expect_true(nzchar(e$message))
  expect_true(is.character(e$where) && length(e$where) > 0)
  expect_false(any(grepl("Yes|Treated", e$where)))          # call names only
  expect_true(any(grepl("estimate_tau_h_internal_tilted|match.arg", e$where)))
})
