# Settings step of the Shiny app, as plain R functions (no shiny code, easy to test).
#
# build_analysis_spec() takes the prepared data (from prepare_analysis_data()) plus
# the user's choices, checks them, and returns an "analysis plan": exactly which
# package function to call with which arguments (Step 4 will run it), a short
# plain-language description, and any warnings.

# Learners the app offers for SuperLearner, and the R package each one needs.
SL_PRESETS <- c(SL.glm = NA, SL.mean = NA, SL.glmnet = "glmnet",
                SL.ranger = "ranger", SL.earth = "earth", SL.xgboost = "xgboost")

# Estimands each method offers (names = what the user sees).
estimand_choices <- function(setting) {
  switch(setting,
    fixed    = c("Tilted ATE with your weights h(x)" = "fixed"),
    internal = c("ATE" = "ATE", "ATT (treated)" = "ATT", "ATC (controls)" = "ATC", "ATO (overlap)" = "ATO",
                 "Custom: your own eta(g)" = "custom"),
    external = c("ATE-I (internal population)" = "ATE-I", "ATE-E (external population)" = "ATE-E",
                 "ATE-O (overlap population)" = "ATE-O", "Custom: your own eta(s)" = "custom"))
}

# Check that a nuisance model choice is usable: "glm", or "SL" with >= 1 learner whose
# package is installed. Returns the library to pass on.
check_learner <- function(method, library, label) {
  if (method == "glm") return(c("SL.glm"))        # library is ignored by glm
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("SuperLearner is not installed, so ", label, " must use glm.", call. = FALSE)
  }
  if (length(library) == 0) stop("Choose at least one SuperLearner learner for ", label, ".", call. = FALSE)
  bad <- setdiff(library, names(SL_PRESETS))
  if (length(bad) > 0) stop("Unknown learner(s) for ", label, ": ", paste(bad, collapse = ", "), call. = FALSE)
  need <- stats::na.omit(unname(SL_PRESETS[library]))
  not_installed <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(not_installed) > 0) {
    stop("Learner package(s) not installed for ", label, ": ", paste(not_installed, collapse = ", "), call. = FALSE)
  }
  library
}

# One lower truncation bound l gives bounds c(l, 1 - l); l = 0 means no truncation.
check_lower_bound <- function(l, label) {
  if (length(l) != 1 || is.na(l) || l < 0 || l >= 0.5) {
    stop(label, " lower bound must be a number from 0 (no truncation) up to, but not including, 0.5.", call. = FALSE)
  }
  c(l, 1 - l)
}

# prepared   output of prepare_analysis_data()
# settings   named list of the user's choices (see the defaults below)
# Returns an analysis plan; stops with a readable message on invalid settings.
build_analysis_spec <- function(prepared, settings = list()) {
  defaults <- list(
    estimand = NULL, custom_eta = "g*(1-g)", custom_deta = "",
    Q_method = "glm", Q_library = c("SL.glm", "SL.glmnet"),
    g_method = "glm", g_library = c("SL.glm", "SL.glmnet"),
    s_method = "glm", s_library = c("SL.glm", "SL.glmnet"),
    pi_method = "glm", pi_library = c("SL.glm", "SL.glmnet"),
    g_lower = 0.025, s_lower = 0.025, pi_lower = 0.025,
    bootstrap = FALSE, n_boot = 200, seed = NA, conf_level = 0.95,
    max_iter = 50, tol = 1e-6)
  st <- utils::modifyList(defaults, settings)
  method <- prepared$setting
  warnings <- character()

  # ---- estimand ----
  valid <- unname(estimand_choices(method))
  if (is.null(st$estimand)) st$estimand <- valid[1]
  if (!(st$estimand %in% valid)) {
    stop("Estimand '", st$estimand, "' is not available for the ", method, " method.", call. = FALSE)
  }

  # ---- nuisance models ----
  Q_lib <- check_learner(st$Q_method, st$Q_library, "the outcome model")
  g_lib <- check_learner(st$g_method, st$g_library, "the treatment model")
  s_lib <- if (method == "external") check_learner(st$s_method, st$s_library, "the selection model") else "SL.glm"
  has_missing_Y <- prepared$summary$n_outcome_missing > 0
  pi_lib <- if (has_missing_Y) check_learner(st$pi_method, st$pi_library, "the missing-outcome model") else "SL.glm"

  # ---- truncation ----
  gbounds <- check_lower_bound(st$g_lower, "Treatment propensity")
  sbounds <- if (method == "external") check_lower_bound(st$s_lower, "Selection propensity") else NULL
  pibounds <- if (has_missing_Y) { b <- check_lower_bound(st$pi_lower, "Missing-outcome propensity"); c(b[1], 1) } else c(0.025, 1)

  # ---- inference and iteration controls ----
  if (!is.numeric(st$conf_level) || st$conf_level <= 0.5 || st$conf_level >= 0.999) {
    stop("Confidence level must be between 0.5 and 0.999.", call. = FALSE)
  }
  if (isTRUE(st$bootstrap) && (is.na(st$n_boot) || st$n_boot < 20 || st$n_boot > 5000)) {
    stop("Number of bootstrap resamples must be between 20 and 5000.", call. = FALSE)
  }
  if (is.na(st$max_iter) || st$max_iter < 1) stop("Maximum iterations must be at least 1.", call. = FALSE)
  if (is.na(st$tol) || st$tol <= 0) stop("Tolerance must be positive.", call. = FALSE)

  # ---- arguments for the package function (names match the package's arguments) ----
  common <- list(
    Q_method = st$Q_method, Q_SL.library = Q_lib,
    g_method = st$g_method, g_SL.library = g_lib,
    gbounds = gbounds)          # truncation is on by default in every entry point
  if (has_missing_Y) {
    common <- c(common, list(pi_method = st$pi_method, pi_SL.library = pi_lib, pibounds = pibounds))
  }

  custom <- NULL
  if (method == "fixed") {
    # the fixed-h bootstrap only refits the target step and always runs; use a cheap
    # run when the user did not ask for it
    args <- c(common, list(n_boot = if (isTRUE(st$bootstrap)) st$n_boot else 2L))
    fun <- "estimate_tau_h_fixed"
    tilt <- NULL
  } else {
    tilt_var <- if (method == "internal") "g" else "s"      # eta is a function of g (internal) or s (external)
    tilt <- if (st$estimand == "custom") make_tilt(st$custom_eta, st$custom_deta, var = tilt_var)
            else if (method == "internal") tilt_spec(st$estimand) else tilt_spec_external(st$estimand)
    if (st$estimand == "custom") {
      custom <- list(eta = tilt$eta_text, deta = st$custom_deta, var = tilt_var,
                     deta_used = tilt$deta_text, deta_auto = tilt$deta_auto)
    }
    args <- c(common, list(max_iter = st$max_iter, tol = st$tol))
    if (method == "external") {
      args <- c(args, list(s_method = st$s_method, s_SL.library = s_lib, sbounds = sbounds))
    }
    fun <- if (method == "internal") "estimate_tau_h_internal_tilted" else "estimate_tau_h_external_tilted"
  }

  # ---- warnings the user should see before running ----
  uses_SL <- any(c(st$Q_method, st$g_method, if (method == "external") st$s_method,
                   if (has_missing_Y) st$pi_method) == "SL")
  if (uses_SL && isTRUE(st$bootstrap) && method != "fixed" && st$n_boot > 100) {
    warnings <- c(warnings, "SuperLearner with a large bootstrap can take a long time.")
  }
  if (st$estimand %in% c("ATT", "ATE-E") && !isTRUE(st$bootstrap)) {
    warnings <- c(warnings, paste0("For the ", st$estimand, ", EIF-based intervals can be too narrow in ",
                                   "moderate samples; consider turning on the bootstrap."))
  }
  if (method == "fixed" && isTRUE(st$bootstrap)) {
    warnings <- c(warnings, "For the fixed-h method the bootstrap refits only the target step (nuisance fits are held fixed).")
  }

  # ---- plain-language description ----
  learner_txt <- function(m, lib) if (m == "glm") "glm" else paste0("SuperLearner (", paste(sub("^SL\\.", "", lib), collapse = ", "), ")")
  lines <- c(
    paste0("Method: ", method, "; estimand: ", if (method == "fixed") "tilted ATE with your weights h(x)" else if (st$estimand == "custom") "custom tilting function" else st$estimand),
    if (!is.null(custom)) paste0("Tilting function: eta(", custom$var, ") = ", custom$eta, "; derivative = ", custom$deta_used,
                                 if (custom$deta_auto) " (computed automatically)" else ""),
    paste0("Outcome model: ", learner_txt(st$Q_method, Q_lib)),
    paste0("Treatment model: ", learner_txt(st$g_method, g_lib), " (truncated at ", gbounds[1], " / ", gbounds[2], ")"),
    if (method == "external") paste0("Selection model: ", learner_txt(st$s_method, s_lib), " (truncated at ", sbounds[1], " / ", sbounds[2], ")"),
    if (has_missing_Y) paste0("Missing-outcome model: ", learner_txt(st$pi_method, pi_lib), " (lower bound ", pibounds[1], ")"),
    paste0("Confidence level: ", format(100 * st$conf_level), "%; bootstrap: ",
           if (isTRUE(st$bootstrap)) paste0(st$n_boot, " resamples") else "off (EIF-based intervals only)"))

  list(method = method, estimand = st$estimand, fun = fun, args = args,
       eta = tilt$eta, eta_dot = tilt$eta_dot,
       custom = custom,
       bootstrap = list(on = isTRUE(st$bootstrap), n_boot = st$n_boot, seed = st$seed),
       conf_level = st$conf_level, description = lines, warnings = warnings)
}
