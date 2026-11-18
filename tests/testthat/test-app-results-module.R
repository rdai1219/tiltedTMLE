## Tests for the Shiny results module (inst/app/R/mod_results.R) via shiny::testServer(),
## with fake run reactives built from a real synthetic run.

skip_if_not_installed("shiny")
skip_if_not_installed("ggplot2")

module_env3 <- new.env(parent = environment(prepare_analysis_data))
sys.source(system.file("app", "R", "mod_results.R", package = "tiltedTMLE"), envir = module_env3)

make_run <- function(setting = "internal", estimand = "ATT") {
  set.seed(71)
  n <- 140
  df <- data.frame(y = sample(c("No", "Yes"), n, TRUE), a = sample(c("C", "T"), n, TRUE),
                   x1 = rnorm(n), x2 = rnorm(n))
  if (setting == "external") {
    df$src <- rep(c("trial", "ext"), each = n / 2); df$y[df$src == "ext"] <- NA; df$a[df$src == "ext"] <- NA
  }
  p <- prepare_analysis_data(df, setting, "y", "a", c("x1", "x2"), "Yes", "T",
                             d = if (setting == "external") "src", d_one = if (setting == "external") "trial")
  s <- build_analysis_spec(p, list(estimand = estimand))
  list(prepared = p, spec = s, main = run_analysis_main(p, s))
}

test_that("results module renders summary, table, diagnostics and script for a finished run", {
  r <- make_run()
  run <- list(inputs = shiny::reactive(list(prepared = r$prepared, spec = r$spec)),
              main = shiny::reactive(r$main), boot = shiny::reactive(NULL))
  shiny::testServer(module_env3$mod_results_server, args = list(run = run), {
    expect_match(as.character(output$tabs$html), "Reproduce &amp; report|Reproduce & report")
    expect_match(as.character(output$summary$html), "ATT")
    expect_match(as.character(output$table), "TMLE")
    expect_match(output$script, "estimate_tau_h_internal_tilted")
    expect_match(as.character(output$diag_table), "Effective sample size")
  })
})

test_that("results module shows a placeholder before any run and an error after a failed run", {
  none <- list(inputs = shiny::reactive(NULL), main = shiny::reactive(NULL), boot = shiny::reactive(NULL))
  shiny::testServer(module_env3$mod_results_server, args = list(run = none), {
    expect_error(output$tabs)                     # no tabs before a run
    expect_match(as.character(output$message$html), "Run an analysis")
  })
  r <- make_run()
  bad <- structure(list(message = "model failed"), class = "run_error")
  failed <- list(inputs = shiny::reactive(list(prepared = r$prepared, spec = r$spec)),
                 main = shiny::reactive(bad), boot = shiny::reactive(NULL))
  shiny::testServer(module_env3$mod_results_server, args = list(run = failed), {
    expect_error(output$tabs)                     # no tabs after a failed run
    expect_match(as.character(output$message$html), "model failed")
  })
})
