## Tests for the app's data step: prepare_analysis_data() and its helpers.
## All data here are synthetic (generated in the test).

## Small synthetic data frame with text-coded binary columns, as an uploaded file might have
make_app_df <- function(n = 120, external = FALSE) {
  set.seed(11)
  df <- data.frame(
    outcome = sample(c("No", "Yes"), n, replace = TRUE),
    arm     = sample(c("Control", "Treated"), n, replace = TRUE),
    age     = round(rnorm(n, 60, 8)),
    sex     = sample(c("F", "M"), n, replace = TRUE),
    `bmi value` = rnorm(n, 27, 4),          # name with a space on purpose
    w = runif(n, 0.5, 2),                   # tilting weights h(x) for the fixed method
    check.names = FALSE, stringsAsFactors = FALSE)
  if (external) {
    df$source <- rep(c("trial", "external"), c(n / 2, n / 2))
    df$outcome[df$source == "external"] <- NA   # external rows have no outcome/treatment
    df$arm[df$source == "external"] <- NA
  }
  df
}

## the fixed-h setting needs the weight column; this helper supplies it
prep_fx <- function(df, ...) prepare_analysis_data(df, "fixed", ..., h = "w")

test_that("binary_levels() and recode_binary() work and refuse bad input", {
  x <- c("No", "Yes", NA, "Yes")
  expect_equal(binary_levels(x), c("No", "Yes"))
  expect_equal(recode_binary(x, "Yes", "Y"), c(0L, 1L, NA, 1L))
  expect_error(recode_binary(x, "Maybe", "Y"), "choose which value counts as 1")
  expect_error(recode_binary(c("a", "b", "c"), "a", "Y"), "exactly 2 distinct")
  expect_error(recode_binary(c(1, 1, 1), "1", "Y"), "exactly 2 distinct")
  ## numeric 0/1 columns work too (compared as text)
  expect_equal(recode_binary(c(0, 1, 1, NA), "1", "A"), c(0L, 1L, 1L, NA))
})

test_that("fixed/internal setting: recoding, covariates and summary", {
  df <- make_app_df()
  p <- prepare_analysis_data(df, "internal", y = "outcome", a = "arm",
                             covariates = c("age", "sex", "bmi value"),
                             y_one = "Yes", a_one = "Treated")
  expect_equal(p$Y, as.integer(df$outcome == "Yes"))
  expect_equal(p$A, as.integer(df$arm == "Treated"))
  expect_named(p$W, c("age", "sex", "bmi.value"))           # safe name
  expect_true(is.factor(p$W$sex))                            # text -> factor
  expect_equal(p$summary$n_treated + p$summary$n_control, 120)
  expect_true(any(grepl("valid R names", p$notes)))
  ## the result is directly usable by an estimator
  out <- estimate_tau_h_fixed(p$Y, p$A, p$W, h = rep(1, length(p$Y)), n_boot = 5)
  expect_true(is.finite(out$TMLE))
})

test_that("missing outcomes are kept as NA and noted", {
  df <- make_app_df()
  df$outcome[1:10] <- NA
  p <- prepare_analysis_data(df, "internal", "outcome", "arm", c("age", "sex"), "Yes", "Treated")
  expect_equal(sum(is.na(p$Y)), 10)
  expect_equal(p$summary$n_outcome_missing, 10)
  expect_true(any(grepl("missing at random", p$notes)))
})

test_that("treatment must be complete; roles must be unique; columns must exist", {
  df <- make_app_df()
  df$arm[3] <- NA
  expect_error(prep_fx(df, "outcome", "arm", "age", "Yes", "Treated"),
               "`A` \\(treatment\\) must not have missing")
  df <- make_app_df()
  expect_error(prep_fx(df, "outcome", "arm", c("age", "arm"), "Yes", "Treated"),
               "only have one role")
  expect_error(prep_fx(df, "outcome", "arm", "nope", "Yes", "Treated"),
               "not found")
  expect_error(prep_fx(df, "outcome", "arm", character(), "Yes", "Treated"),
               "at least one covariate")
})

test_that("missing covariates: indicators added by default, error when switched off", {
  df <- make_app_df()
  df$age[1:6] <- NA
  df$sex[2] <- NA
  p <- prep_fx(df, "outcome", "arm", c("age", "sex"), "Yes", "Treated")
  expect_false(anyNA(p$W))
  expect_true(all(c("age_miss", "sex_miss") %in% names(p$W)))
  expect_equal(sum(p$W$age_miss), 6)
  expect_equal(unname(p$summary$covariate_missing[c("age", "sex")]), c(6L, 1L))
  expect_error(prep_fx(df, "outcome", "arm", c("age", "sex"), "Yes", "Treated",
                                     add_indicators = FALSE),
               "missing values")
})

test_that("constant covariates are dropped with a note", {
  df <- make_app_df(); df$site <- "A"
  p <- prep_fx(df, "outcome", "arm", c("age", "site"), "Yes", "Treated")
  expect_named(p$W, "age")
  expect_true(any(grepl("Dropped constant", p$notes)))
  expect_error(prep_fx(df, "outcome", "arm", "site", "Yes", "Treated"),
               "all were constant")
})

test_that("external setting: one table is split into internal-only Y, A and pooled Ws, D", {
  df <- make_app_df(120, external = TRUE)
  p <- prepare_analysis_data(df, "external", "outcome", "arm", c("age", "sex"), "Yes", "Treated",
                             d = "source", d_one = "trial")
  idx1 <- which(p$D == 1)
  expect_length(p$Y, 60); expect_length(p$A, 60)             # internal rows only
  expect_equal(nrow(p$Ws), 120); expect_length(p$D, 120)     # pooled
  expect_equal(p$Y, as.integer(df$outcome[idx1] == "Yes"))
  expect_equal(p$summary$n_external, 60)
  ## and it runs through the external estimator
  spec <- tilt_spec_external("ATE-O")
  out <- estimate_tau_h_external_tilted(p$Y, p$A, p$D, p$Ws, spec$eta, spec$eta_dot)
  expect_true(is.finite(out$TMLE))
})

test_that("external setting: blank Y on external rows is fine, blank Y on internal rows is a missing outcome", {
  df <- make_app_df(120, external = TRUE)
  df$outcome[5] <- NA                                   # an INTERNAL row with missing outcome
  p <- prepare_analysis_data(df, "external", "outcome", "arm", c("age", "sex"), "Yes", "Treated",
                             d = "source", d_one = "trial")
  expect_equal(p$summary$n_outcome_missing, 1)
})

test_that("external setting: checks on D and on factor levels", {
  df <- make_app_df(120, external = TRUE)
  expect_error(prepare_analysis_data(df, "external", "outcome", "arm", "age", "Yes", "Treated"),
               "Choose the column that marks")
  df_one <- df; df_one$source <- "trial"
  expect_error(prepare_analysis_data(df_one, "external", "outcome", "arm", "age", "Yes", "Treated",
                                     d = "source", d_one = "trial"),
               "exactly 2 distinct")
  ## a covariate level that exists only among external rows cannot be predicted
  df$sex[df$source == "external"][1] <- "X"
  expect_error(prepare_analysis_data(df, "external", "outcome", "arm", c("age", "sex"), "Yes", "Treated",
                                     d = "source", d_one = "trial"),
               "only in external rows")
})

test_that("covariates named like the models' internal variables (A, Y, Delta) are renamed", {
  df <- make_app_df()
  df$A <- rnorm(nrow(df)); df$Y <- rnorm(nrow(df))
  p <- prep_fx(df, "outcome", "arm", c("age", "A", "Y"), "Yes", "Treated")
  expect_true(all(c("A_cov", "Y_cov") %in% names(p$W)))
  expect_false(any(c("A", "Y") %in% names(p$W)))
  expect_true(any(grepl("renamed", p$notes)))
  out <- estimate_tau_h_fixed(p$Y, p$A, p$W, h = rep(1, length(p$Y)), n_boot = 3)
  expect_true(is.finite(out$TMLE))
})

test_that("fixed setting: the weight column is required, checked, and summarised", {
  df <- make_app_df()
  expect_error(prepare_analysis_data(df, "fixed", "outcome", "arm", "age", "Yes", "Treated"),
               "choose the column that holds the tilting weights")
  p <- prep_fx(df, "outcome", "arm", "age", "Yes", "Treated")
  expect_equal(p$h, df$w)
  expect_equal(p$summary$h$min, min(df$w))
  expect_lte(p$summary$h$ess, nrow(df))
  df$w[1] <- -1
  expect_error(prep_fx(df, "outcome", "arm", "age", "Yes", "Treated"), "no negative")
  df$w <- 0
  expect_error(prep_fx(df, "outcome", "arm", "age", "Yes", "Treated"), "positive sum")
  ## the weights column cannot also be a covariate
  df$w <- runif(nrow(df), 0.5, 2)
  expect_error(prep_fx(df, "outcome", "arm", c("age", "w"), "Yes", "Treated"), "only have one role")
})
