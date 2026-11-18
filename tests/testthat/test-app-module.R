## Tests for the Shiny data module (inst/app/R/mod_data.R) using shiny::testServer(),
## with a small synthetic CSV written to a temp file.

skip_if_not_installed("shiny")
skip_if_not_installed("readr")

## Load the module code into an environment whose parent is the package namespace
module_env <- new.env(parent = environment(prepare_analysis_data))
sys.source(system.file("app", "R", "mod_data.R", package = "tiltedTMLE"), envir = module_env)

make_csv <- function(external = FALSE) {
  set.seed(21)
  n <- 100
  df <- data.frame(
    outcome = sample(c("No", "Yes"), n, replace = TRUE),
    arm = sample(c("Control", "Treated"), n, replace = TRUE),
    age = round(rnorm(n, 60, 8)),
    sex = sample(c("F", "M"), n, replace = TRUE))
  df$age[1:4] <- NA                                      # some missing covariates
  if (external) {
    df$source <- rep(c("trial", "external"), each = n / 2)
    df$outcome[df$source == "external"] <- NA
    df$arm[df$source == "external"] <- NA
  }
  path <- tempfile(fileext = ".csv")
  utils::write.csv(df, path, row.names = FALSE)
  list(path = path, df = df)
}

test_that("module: upload, map columns and validate (single study)", {
  f <- make_csv()
  shiny::testServer(module_env$mod_data_server, {
    session$setInputs(file = data.frame(name = "d.csv", datapath = f$path, stringsAsFactors = FALSE),
                      setting = "internal", add_indicators = TRUE, show_preview = FALSE,
                      y = "outcome", y_one = "Yes", a = "arm", a_one = "Treated",
                      covariates = c("age", "sex"), validate = 1)
    p <- session$returned()
    expect_false(inherits(p, "prep_error"))
    expect_equal(length(p$Y), 100)
    expect_true("age_miss" %in% names(p$W))
    expect_equal(p$summary$n_rows, 100)
  })
})

test_that("module: changing the method after a check marks the data as needing a new check", {
  f <- make_csv()
  shiny::testServer(module_env$mod_data_server, {
    session$setInputs(file = data.frame(name = "d.csv", datapath = f$path, stringsAsFactors = FALSE),
                      setting = "internal", add_indicators = TRUE, show_preview = FALSE,
                      y = "outcome", y_one = "Yes", a = "arm", a_one = "Treated",
                      covariates = c("age", "sex"), validate = 1)
    expect_equal(session$returned()$setting, "internal")
    session$setInputs(setting = "external")                    # switch method, no new check
    p <- session$returned()
    expect_s3_class(p, "prep_error")
    expect_match(p$message, "Check data")
    session$setInputs(setting = "internal")                    # back to the checked choices
    expect_false(inherits(session$returned(), "prep_error"))
  })
})

test_that("module: external setting splits the pooled table", {
  f <- make_csv(external = TRUE)
  shiny::testServer(module_env$mod_data_server, {
    session$setInputs(file = data.frame(name = "d.csv", datapath = f$path, stringsAsFactors = FALSE),
                      setting = "external", add_indicators = TRUE, show_preview = FALSE,
                      y = "outcome", y_one = "Yes", a = "arm", a_one = "Treated",
                      d = "source", d_one = "trial", covariates = c("age", "sex"), validate = 1)
    p <- session$returned()
    expect_false(inherits(p, "prep_error"))
    expect_equal(length(p$Y), 50)
    expect_equal(nrow(p$Ws), 100)
  })
})

test_that("module: a bad choice returns a readable error, not a crash", {
  f <- make_csv()
  shiny::testServer(module_env$mod_data_server, {
    session$setInputs(file = data.frame(name = "d.csv", datapath = f$path, stringsAsFactors = FALSE),
                      setting = "internal", add_indicators = FALSE, show_preview = FALSE,
                      y = "outcome", y_one = "Yes", a = "arm", a_one = "Treated",
                      covariates = c("age", "sex"), validate = 1)
    p <- session$returned()
    expect_s3_class(p, "prep_error")
    expect_match(p$message, "missing values")
  })
})
