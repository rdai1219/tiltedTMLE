## Tests for the Shiny settings module (inst/app/R/mod_settings.R) via shiny::testServer(),
## with a synthetic `prepared` reactive standing in for the data module.

skip_if_not_installed("shiny")

module_env2 <- new.env(parent = environment(prepare_analysis_data))
sys.source(system.file("app", "R", "mod_settings.R", package = "tiltedTMLE"), envir = module_env2)

make_prepared <- function(setting = "internal", missing_Y = FALSE) {
  set.seed(41)
  n <- 100
  df <- data.frame(y = sample(c("No", "Yes"), n, TRUE), a = sample(c("C", "T"), n, TRUE),
                   x1 = rnorm(n), x2 = rnorm(n))
  if (missing_Y) df$y[1:6] <- NA
  if (setting == "external") {
    df$src <- rep(c("trial", "ext"), each = n / 2)
    df$y[df$src == "ext"] <- NA; df$a[df$src == "ext"] <- NA
  }
  prepare_analysis_data(df, setting, "y", "a", c("x1", "x2"), "Yes", "T",
                        d = if (setting == "external") "src", d_one = if (setting == "external") "trial")
}

test_that("settings module returns a plan that follows the inputs", {
  prepared <- shiny::reactive(make_prepared("internal"))
  shiny::testServer(module_env2$mod_settings_server, args = list(prepared = prepared), {
    session$setInputs(estimand = "ATO", Q_method = "glm", g_method = "glm", bootstrap = FALSE,
                      conf_level = 0.9, g_lower = 0.05, max_iter = 30, tol = 1e-5)
    s <- session$returned()
    expect_false(inherits(s, "spec_error"))
    expect_equal(s$fun, "estimate_tau_h_internal_tilted")
    expect_equal(s$estimand, "ATO")
    expect_equal(s$args$gbounds, c(0.05, 0.95))
    expect_equal(s$args$max_iter, 30)
    expect_equal(s$conf_level, 0.9)
  })
})

test_that("settings module turns invalid input into a readable error object", {
  prepared <- shiny::reactive(make_prepared("internal"))
  shiny::testServer(module_env2$mod_settings_server, args = list(prepared = prepared), {
    session$setInputs(estimand = "ATT", Q_method = "glm", g_method = "glm", bootstrap = TRUE,
                      n_boot = 5, conf_level = 0.95, g_lower = 0.025, max_iter = 50, tol = 1e-6)
    s <- session$returned()
    expect_s3_class(s, "spec_error")
    expect_match(s$message, "between 20 and 5000")
  })
})

test_that("external data: plan has a selection model and bounds", {
  prepared <- shiny::reactive(make_prepared("external"))
  shiny::testServer(module_env2$mod_settings_server, args = list(prepared = prepared), {
    session$setInputs(estimand = "ATE-O", Q_method = "glm", g_method = "glm", s_method = "glm",
                      bootstrap = FALSE, conf_level = 0.95, g_lower = 0.025, s_lower = 0.1,
                      max_iter = 50, tol = 1e-6)
    s <- session$returned()
    expect_equal(s$fun, "estimate_tau_h_external_tilted")
    expect_equal(s$args$sbounds, c(0.1, 0.9))
  })
})

test_that("missing outcomes: the plan includes the missing-outcome model", {
  prepared <- shiny::reactive(make_prepared("internal", missing_Y = TRUE))
  shiny::testServer(module_env2$mod_settings_server, args = list(prepared = prepared), {
    session$setInputs(estimand = "ATE", Q_method = "glm", g_method = "glm", pi_method = "glm",
                      bootstrap = FALSE, conf_level = 0.95, g_lower = 0.025, pi_lower = 0.05,
                      max_iter = 50, tol = 1e-6)
    s <- session$returned()
    expect_equal(s$args$pibounds, c(0.05, 1))
  })
})

test_that("settings module builds a custom-tilt plan and reports a bad formula readably", {
  prepared <- shiny::reactive(make_prepared("internal"))
  shiny::testServer(module_env2$mod_settings_server, args = list(prepared = prepared), {
    session$setInputs(estimand = "custom", custom_eta = "g*(1-g)", custom_deta = "", Q_method = "glm", g_method = "glm",
                      bootstrap = FALSE, conf_level = 0.95, g_lower = 0.025, max_iter = 50, tol = 1e-6)
    s <- session$returned()
    expect_false(inherits(s, "spec_error"))
    expect_equal(s$custom$var, "g")
    session$setInputs(custom_eta = "system('x')")
    expect_s3_class(session$returned(), "spec_error")
    expect_match(session$returned()$message, "not allowed")
  })
})
