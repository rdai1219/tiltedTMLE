## Tests for the app's results step: tables, diagnostics, plots, summary text and the
## reproduction script. Synthetic data only.

skip_if_not_installed("ggplot2")

## Prepared data, a plan and a main-fit bundle for one method
run_for <- function(setting, estimand, missing_Y = FALSE, settings = list()) {
  set.seed(61)
  n <- 160
  df <- data.frame(response = sample(c("No", "Yes"), n, TRUE), arm = sample(c("Control", "Treated"), n, TRUE),
                   x1 = rnorm(n), x2 = rnorm(n), site = "A", w = runif(n, 0.5, 2), stringsAsFactors = FALSE)
  df$x1[sample(n, 8)] <- NA
  if (missing_Y) df$response[1:12] <- NA
  if (setting == "external") {
    df$src <- rep(c("trial", "ext"), each = n / 2)
    df$response[df$src == "ext"] <- NA; df$arm[df$src == "ext"] <- NA
  }
  p <- prepare_analysis_data(df, setting, "response", "arm", c("x1", "x2", "site"), "Yes", "Treated",
                             d = if (setting == "external") "src", d_one = if (setting == "external") "trial",
                             h = if (setting == "fixed") "w")
  s <- build_analysis_spec(p, utils::modifyList(list(estimand = estimand), settings))
  list(df = df, prepared = p, spec = s, main = run_analysis_main(p, s))
}

test_that("effective sample size and bound share behave sensibly", {
  expect_equal(effective_sample_size(rep(1, 50)), 50)
  expect_lt(effective_sample_size(c(10, rep(1, 49))), 50)
  expect_equal(share_at_bounds(c(0.025, 0.5, 0.975, 0.4), c(0.025, 0.975)), 0.5)
  expect_equal(share_at_bounds(c(0.1, 0.9), c(0, 1)), 0)
})

test_that("results table adds bootstrap rows only when a bootstrap result is given", {
  r <- run_for("internal", "ATT", settings = list(bootstrap = TRUE, n_boot = 20, seed = 1))
  expect_equal(nrow(results_table(r$main)), 4)
  b <- run_analysis_bootstrap(r$prepared, r$spec)
  expect_equal(nrow(results_table(r$main, b)), 6)
  expect_equal(nrow(results_table(r$main, structure(list(message = "x"), class = "run_error"))), 4)
})

test_that("diagnostics table has the expected rows, per method", {
  for (cfg in list(c("fixed", "fixed"), c("internal", "ATO"), c("external", "ATE-O"))) {
    r <- run_for(cfg[1], cfg[2], missing_Y = cfg[1] != "external")
    d <- diagnostics_table(r$main)
    expect_named(d, c("Diagnostic", "Value"))
    expect_true(any(grepl("Effective sample size", d$Diagnostic)), label = cfg[1])
    expect_equal(any(grepl("Selection propensity", d$Diagnostic)), cfg[1] == "external", label = cfg[1])
    expect_equal(any(grepl("converged", d$Diagnostic)), cfg[1] != "fixed", label = cfg[1])
  }
})

test_that("all three plots build for every method", {
  for (cfg in list(c("fixed", "fixed"), c("internal", "ATT"), c("external", "ATE-O"))) {
    r <- run_for(cfg[1], cfg[2])
    for (f in list(plot_overlap, plot_weights, plot_estimators)) {
      p <- f(r$main)
      expect_s3_class(p, "ggplot")
      tmp <- tempfile(fileext = ".png")
      expect_no_error(suppressWarnings(ggplot2::ggsave(tmp, p, width = 5, height = 4)))
      expect_true(file.exists(tmp))
    }
  }
  r <- run_for("internal", "ATT", settings = list(bootstrap = TRUE, n_boot = 20, seed = 1))
  b <- run_analysis_bootstrap(r$prepared, r$spec)
  expect_s3_class(plot_estimators(r$main, b), "ggplot")
})

test_that("summary text carries the estimate, the level and the caveats", {
  r <- run_for("internal", "ATT", missing_Y = TRUE, settings = list(conf_level = 0.9))
  s <- summary_text(r$main)
  tm <- r$main$table[r$main$table$Estimator == "TMLE", ]
  expect_match(s$sentence, "ATT")
  expect_match(s$sentence, sprintf("%.3f", tm$Estimate), fixed = TRUE)
  expect_match(s$sentence, "90% CI")
  expect_true(any(grepl("no unmeasured confounding", s$caveats)))
  expect_true(any(grepl("12 outcome\\(s\\) were missing", s$caveats)))
  r0 <- run_for("internal", "ATT")
  expect_false(any(grepl("were missing", summary_text(r0$main)$caveats)))
})

test_that("the reproduction script reruns and gives the same TMLE (all three methods)", {
  for (cfg in list(c("fixed", "fixed"), c("internal", "ATO"), c("external", "ATE-O"))) {
    r <- run_for(cfg[1], cfg[2], missing_Y = cfg[1] == "internal")
    csv <- tempfile(fileext = ".csv")
    utils::write.csv(r$df, csv, row.names = FALSE, na = "")
    code <- reproduction_script(r$prepared, r$spec, file_name = csv)
    env <- new.env(parent = globalenv())
    suppressWarnings(eval(parse(text = code), envir = env))
    from_app <- r$main$table$Estimate[r$main$table$Estimator == "TMLE"]
    expect_equal(env$fit$TMLE, from_app, label = cfg[1])
  }
})

test_that("the reproduction script includes the bootstrap and the seed when requested", {
  r <- run_for("internal", "ATT", settings = list(bootstrap = TRUE, n_boot = 20, seed = 5))
  code <- reproduction_script(r$prepared, r$spec)
  expect_match(code, "set.seed(5)", fixed = TRUE)
  expect_match(code, "bootstrap_tau_h_internal_tilted")
  expect_match(code, "n_boot = 20")
  r2 <- run_for("internal", "ATT")
  expect_false(grepl("bootstrap_tau_h", reproduction_script(r2$prepared, r2$spec)))
})

test_that("the HTML report renders, contains the key pieces and no individual-level data", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  skip_on_cran()
  r <- run_for("internal", "ATT", settings = list(bootstrap = TRUE, n_boot = 20, seed = 2))
  b <- run_analysis_bootstrap(r$prepared, r$spec)
  out <- tempfile(fileext = ".html")
  render_report(r$main, b, reproduction_script(r$prepared, r$spec), out)
  expect_true(file.exists(out))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "tiltedTMLE analysis report")
  expect_match(html, "ATT")
  expect_match(html, "full-refit bootstrap")
  expect_match(html, "Effective sample size")
  expect_match(html, "estimate_tau_h_internal_tilted")
  ## summary-level only: no column from the raw data frame appears as a data column header
  expect_false(grepl("<th[^>]*>x1</th>", html))
})
