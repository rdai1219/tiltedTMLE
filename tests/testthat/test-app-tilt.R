## Tests for custom tilting functions: make_tilt() and how they flow through the app's plan,
## run step and reproduction script. Synthetic data only.

test_that("make_tilt() reproduces the preset tilting functions and derives derivatives automatically", {
  g <- c(0.1, 0.4, 0.8)
  for (case in list(list("1", "ATE"), list("g", "ATT"), list("1-g", "ATC"), list("g*(1-g)", "ATO"))) {
    t <- make_tilt(case[[1]], "", "g")
    ref <- tilt_spec(case[[2]])
    expect_equal(t$eta(g), ref$eta(g), label = case[[2]])
    expect_equal(t$eta_dot(g), ref$eta_dot(g), label = case[[2]])
    expect_true(t$deta_auto)
  }
  ## external presets, variable s
  ts <- make_tilt("(1-s)/s", "", "s")
  expect_equal(ts$eta(c(0.2, 0.5)), tilt_spec_external("ATE-E")$eta(c(0.2, 0.5)))
  expect_equal(ts$eta_dot(c(0.2, 0.5)), tilt_spec_external("ATE-E")$eta_dot(c(0.2, 0.5)))
})

test_that("a typed derivative is used as given; non-differentiable formulas ask for one", {
  t <- make_tilt("g*(1-g)", "1-2*g", "g")
  expect_false(t$deta_auto)
  expect_equal(t$eta_dot(0.25), 0.5)
  expect_error(make_tilt("pmin(g, 1-g)", "", "g"), "type it in the derivative box")
  t2 <- make_tilt("pmin(g, 1-g)", "ifelse(g < 0.5, 1, -1)", "g")      # matching-weights style
  expect_equal(t2$eta(c(0.2, 0.9)), c(0.2, 0.1))
  expect_equal(t2$eta_dot(c(0.2, 0.9)), c(1, -1))
})

test_that("a constant eta returns one value per input", {
  t <- make_tilt("2", "", "g")
  expect_equal(t$eta(c(0.1, 0.5, 0.9)), rep(2, 3))
  expect_equal(t$eta_dot(c(0.1, 0.5, 0.9)), rep(0, 3))
})

test_that("invalid and unsafe formulas are rejected with readable messages", {
  expect_error(make_tilt("", "", "g"), "Enter a formula")
  expect_error(make_tilt("(g", "", "g"), "Cannot read")
  expect_error(make_tilt("h*(1-h)", "", "g"), "may only use the variable 'g'")
  expect_error(make_tilt("s", "", "g"), "may only use the variable 'g'")
  expect_error(make_tilt("g - 0.5", "", "g"), "must not be negative")
  expect_error(make_tilt("log(g - 0.5)", "", "g"), "not finite")
  expect_error(make_tilt("0*g", "", "g"), "zero everywhere")
  ## nothing but arithmetic and the allowed functions can run
  expect_error(make_tilt("system('echo hi')", "", "g"), "not allowed")
  expect_error(make_tilt("base::sqrt(g)", "", "g"), "not allowed")
  expect_error(make_tilt("{g}", "", "g"), "not allowed")
  expect_error(make_tilt("g <- 2", "", "g"), "not allowed")
  expect_error(make_tilt("x <- g", "", "g"))                      # also rejected (unknown variable)
  expect_error(make_tilt("function(g) g", "", "g"), "not allowed")
  expect_error(make_tilt("g", "unlink('a')", "g"), "not allowed")
})

## prepared data for a method
prep_t <- function(setting) {
  set.seed(91)
  n <- 160
  df <- data.frame(y = sample(c("No", "Yes"), n, TRUE), a = sample(c("C", "T"), n, TRUE),
                   x1 = rnorm(n), x2 = rnorm(n))
  if (setting == "external") {
    df$src <- rep(c("trial", "ext"), each = n / 2); df$y[df$src == "ext"] <- NA; df$a[df$src == "ext"] <- NA
  }
  prepare_analysis_data(df, setting, "y", "a", c("x1", "x2"), "Yes", "T",
                        d = if (setting == "external") "src", d_one = if (setting == "external") "trial")
}

test_that("a custom plan equals the preset it spells out (internal ATO and external ATE-O)", {
  p <- prep_t("internal")
  s_c <- build_analysis_spec(p, list(estimand = "custom", custom_eta = "g*(1-g)"))
  s_p <- build_analysis_spec(p, list(estimand = "ATO"))
  r_c <- run_analysis_main(p, s_c); r_p <- run_analysis_main(p, s_p)
  expect_false(inherits(r_c, "run_error"))
  expect_equal(r_c$table$Estimate, r_p$table$Estimate, tolerance = 1e-8)
  expect_equal(r_c$table$SE, r_p$table$SE, tolerance = 1e-8)
  expect_true(any(grepl("Tilting function: eta\\(g\\) = g\\*\\(1-g\\)", s_c$description)))
  expect_true(any(grepl("computed automatically", s_c$description)))

  pe <- prep_t("external")
  se_c <- build_analysis_spec(pe, list(estimand = "custom", custom_eta = "1-s"))
  se_p <- build_analysis_spec(pe, list(estimand = "ATE-O"))
  expect_equal(se_c$custom$var, "s")
  expect_equal(run_analysis_main(pe, se_c)$table$Estimate, run_analysis_main(pe, se_p)$table$Estimate, tolerance = 1e-8)
})

test_that("custom estimand is offered for internal and external methods, with readable plan errors", {
  expect_true("custom" %in% estimand_choices("internal"))
  expect_true("custom" %in% estimand_choices("external"))
  expect_false("custom" %in% estimand_choices("fixed"))
  p <- prep_t("internal")
  expect_error(build_analysis_spec(p, list(estimand = "custom", custom_eta = "q*g")), "may only use the variable 'g'")
  expect_error(build_analysis_spec(p, list(estimand = "custom", custom_eta = "")), "Enter a formula")
})

test_that("the bootstrap and the worker rebuild a custom tilting function", {
  p <- prep_t("internal")
  s <- build_analysis_spec(p, list(estimand = "custom", custom_eta = "g", bootstrap = TRUE, n_boot = 20, seed = 3))
  s_light <- s; s_light$eta <- NULL; s_light$eta_dot <- NULL            # what the background worker receives
  b <- run_analysis_bootstrap(p, s_light)
  expect_false(inherits(b, "run_error"))
  expect_true(is.finite(b$se))
})

test_that("the reproduction script writes a custom tilting function out and reruns to the same TMLE", {
  p <- prep_t("internal")
  df <- data.frame(y = ifelse(is.na(p$Y), NA, ifelse(p$Y == 1, "Yes", "No")), a = ifelse(p$A == 1, "T", "C"),
                   x1 = p$W$x1, x2 = p$W$x2)
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(df, csv, row.names = FALSE, na = "")
  s <- build_analysis_spec(p, list(estimand = "custom", custom_eta = "g*(1-g)"))
  code <- reproduction_script(p, s, file_name = csv)
  expect_match(code, "tilt <- list(eta = function(g) g*(1-g)", fixed = TRUE)
  expect_match(code, "eta_dot = function(g)", fixed = TRUE)
  env <- new.env(parent = globalenv())
  suppressWarnings(eval(parse(text = code), envir = env))
  from_app <- run_analysis_main(p, s)$table$Estimate[1]
  expect_equal(env$fit$TMLE, from_app)
})

test_that("summary text names the custom tilting function generically", {
  p <- prep_t("internal")
  s <- build_analysis_spec(p, list(estimand = "custom", custom_eta = "g"))
  expect_match(summary_text(run_analysis_main(p, s))$sentence, "your own tilting function")
})

test_that("plot_tilt() builds", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot_tilt(make_tilt("g*(1-g)", "", "g")), "ggplot")
})
