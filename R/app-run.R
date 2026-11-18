# Run step of the Shiny app, as plain R functions (no shiny code, easy to test).
#
# The app runs an analysis in two stages, each in a background R process:
#   stage 1  run_analysis_main():      the main fit (seconds): estimates, SEs, intervals,
#                                      diagnostics. Shown to the user immediately.
#   stage 2  run_analysis_bootstrap(): the full-refit bootstrap (internal / external
#                                      methods, only if requested); added when finished.
# Both take the prepared data (prepare_analysis_data()) and the plan (build_analysis_spec()),
# keep everything in memory, and return either a result bundle or a "run_error" object.

# ---- helpers ----

# Wald interval from an estimate and a standard error at the chosen confidence level.
wald_ci <- function(est, se, conf_level) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  c(est - z * se, est + z * se)
}

# One row of the estimates table.
est_row <- function(name, est, se = NA_real_, ci = c(NA_real_, NA_real_), note = "") {
  data.frame(Estimator = name, Estimate = est, SE = se, `CI lower` = ci[1], `CI upper` = ci[2],
             Note = note, check.names = FALSE, stringsAsFactors = FALSE)
}

# Tilting functions for the plan (rebuilt here, so no function needs to cross to the worker).
tilt_for <- function(spec) {
  if (spec$estimand == "custom") return(make_tilt(spec$custom$eta, spec$custom$deta, spec$custom$var))
  if (spec$method == "internal") tilt_spec(spec$estimand) else tilt_spec_external(spec$estimand)
}

# Names of the functions on the call stack (no arguments, so no data), used to show where an
# error happened. Wrapper functions are left out.
call_chain <- function() {
  nm <- vapply(sys.calls(), function(cl) paste(deparse(cl[[1]], nlines = 1), collapse = ""), "")
  nm <- nm[!nm %in% c("withCallingHandlers", "tryCatch", "tryCatchList", "tryCatchOne", "doTryCatch", "do.call", "eval", "eval.parent")]
  utils::tail(unique(nm), 8)
}

# Turn a caught error into a run_error object: message, the failing call and the call chain.
run_error_from <- function(e, chain) {
  structure(list(message = conditionMessage(e),
                 call = if (is.null(conditionCall(e))) NA_character_ else paste(deparse(conditionCall(e), nlines = 1), collapse = ""),
                 where = chain),
            class = "run_error")
}

# Run expr; collect its warnings (as text) instead of printing them.
with_warnings <- function(expr) {
  w <- character()
  value <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(value = value, warnings = unique(w))
}

# ---- stage 1: the main fit ----

# Returns a result bundle (a list) or a "run_error" object with a message.
run_analysis_main <- function(prepared, spec) {
  t0 <- Sys.time()
  chain <- NULL
  out <- tryCatch(withCallingHandlers({
    fit <- with_warnings(fit_main(prepared, spec))
    bundle <- fit$value
    bundle$warnings <- fit$warnings
    bundle$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    bundle$spec_description <- spec$description
    bundle$package_version <- as.character(utils::packageVersion("tiltedTMLE"))
    bundle$r_version <- R.version.string
    bundle
  }, error = function(e) chain <<- call_chain()),
  error = function(e) run_error_from(e, chain))
  out
}

# The estimator call for each method; builds the estimates table and the diagnostics
# the results step needs (propensity scores, influence function). Raw data are not kept.
fit_main <- function(prepared, spec) {
  cl <- spec$conf_level
  method <- spec$method

  if (method == "fixed") {
    out <- do.call(estimate_tau_h_fixed,
                   c(list(prepared$Y, prepared$A, prepared$W, prepared$h), spec$args))
    tab <- rbind(
      est_row("TMLE", out$TMLE, out$TMLE_se, wald_ci(out$TMLE, out$TMLE_se, cl), "EIF-based interval"),
      est_row("AIPW", out$AIPW, out$AIPW_se, wald_ci(out$AIPW, out$AIPW_se, cl), "benchmark"),
      est_row("IPW", out$IPW, note = "benchmark"),
      est_row("MLE (plug-in)", out$MLE, note = "benchmark"))
    boot <- NULL
    if (spec$bootstrap$on) {                       # the fixed-h bootstrap is cheap: done here
      tab <- rbind(tab,
        est_row("TMLE (target-step bootstrap)", out$TMLE, out$Bootstrap_SE,
                wald_ci(out$TMLE, out$Bootstrap_SE, cl), "normal interval; target step only"))
    }
    diagnostics <- list(A = prepared$A, g1W = out$gfit$g1W, h = out$h, EIF = out$tmle_fit$EIF,
                 pi1W = out$pifit$pi1W, pi0W = out$pifit$pi0W)
    converged <- NA; iter <- NA

  } else if (method == "internal") {
    tilt <- tilt_for(spec)
    out <- do.call(estimate_tau_h_internal_tilted,
                   c(list(prepared$Y, prepared$A, prepared$W, tilt$eta, tilt$eta_dot), spec$args))
    tab <- rbind(
      est_row("TMLE", out$TMLE, out$TMLE_se, wald_ci(out$TMLE, out$TMLE_se, cl), "EIF-based interval"),
      est_row("AIPW", out$AIPW, out$AIPW_se, wald_ci(out$AIPW, out$AIPW_se, cl), "benchmark (not targeted)"),
      est_row("IPW", out$IPW, note = "benchmark (not targeted)"),
      est_row("MLE (plug-in)", out$MLE, note = "benchmark (not targeted)"))
    diagnostics <- list(A = prepared$A, g1W = out$tmle_fit$g1W, g1W_initial = out$gfit$g1W,
                 h = out$h, EIF = out$tmle_fit$EIF, pi1W = out$pifit$pi1W, pi0W = out$pifit$pi0W)
    converged <- out$TMLE_converged; iter <- out$TMLE_iter

  } else {  # external
    tilt <- tilt_for(spec)
    out <- do.call(estimate_tau_h_external_tilted,
                   c(list(prepared$Y, prepared$A, prepared$D, prepared$Ws, tilt$eta, tilt$eta_dot), spec$args))
    tab <- rbind(
      est_row("TMLE", out$TMLE, out$TMLE_se, wald_ci(out$TMLE, out$TMLE_se, cl), "EIF-based interval"),
      est_row("IPW", out$IPW, note = "benchmark (not targeted)"),
      est_row("MLE (plug-in)", out$MLE, note = "benchmark (not targeted)"))
    diagnostics <- list(A = prepared$A, D = prepared$D, s1W = out$tmle_fit$s1W, s1W_initial = out$sfit$g1W,
                 g1W = out$gfit$g1W, h = out$h, EIF = out$tmle_fit$EIF,
                 pi1W = out$pifit$pi1W_all, pi0W = out$pifit$pi0W_all)
    converged <- out$TMLE_converged; iter <- out$TMLE_iter
  }

  list(method = method, estimand = spec$estimand, conf_level = cl, table = tab,
       converged = converged, iterations = iter, diagnostics = diagnostics,
       bounds = list(g = spec$args$gbounds, s = spec$args$sbounds),
       n_analysis = prepared$summary$n_analysis,
       n_missing_outcomes = prepared$summary$n_outcome_missing)
}

# ---- stage 2: the full-refit bootstrap (internal and external methods) ----

# Returns list(se, CI_normal, CI_percentile, estimates, n_failed, n_nonconverged, ...)
# or a "run_error".
run_analysis_bootstrap <- function(prepared, spec) {
  t0 <- Sys.time()
  chain <- NULL
  tryCatch(withCallingHandlers({
    if (!is.na(spec$bootstrap$seed)) set.seed(spec$bootstrap$seed)   # reproducible when a seed is given
    tilt <- tilt_for(spec)
    common <- c(list(n_boot = spec$bootstrap$n_boot, conf_level = spec$conf_level), spec$args)
    b <- with_warnings(
      if (spec$method == "internal") {
        do.call(bootstrap_tau_h_internal_tilted,
                c(list(prepared$Y, prepared$A, prepared$W, tilt$eta, tilt$eta_dot), common))
      } else {
        do.call(bootstrap_tau_h_external_tilted,
                c(list(prepared$Y, prepared$A, prepared$D, prepared$Ws, tilt$eta, tilt$eta_dot), common))
      })
    out <- b$value
    out$warnings <- b$warnings
    out$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    out$n_boot <- spec$bootstrap$n_boot
    out$seed <- spec$bootstrap$seed
    out
  }, error = function(e) chain <<- call_chain()),
  error = function(e) run_error_from(e, chain))
}

# ---- background execution ----

# Start `kind` ("main" or "bootstrap") in a separate R process and return a future.
# The worker loads the SAME tiltedTMLE as the calling session: from the source folder
# when the session uses pkgload::load_all() (development), otherwise the installed copy.
# Functions are not sent across (the plan's tilting functions are rebuilt in the worker).
run_in_background <- function(kind = c("main", "bootstrap"), prepared, spec,
                              pkg_path = getNamespaceInfo(asNamespace("tiltedTMLE"), "path")) {
  kind <- match.arg(kind)
  spec$eta <- NULL; spec$eta_dot <- NULL
  future::future({
    if (file.exists(file.path(pkg_path, "DESCRIPTION")) && !dir.exists(file.path(pkg_path, "Meta"))) {
      pkgload::load_all(pkg_path, quiet = TRUE)
    } else {
      loadNamespace("tiltedTMLE")
    }
    # look the function up in the loaded namespace (avoids `:::` inside the package)
    fun <- get(if (kind == "main") "run_analysis_main" else "run_analysis_bootstrap",
               envir = asNamespace("tiltedTMLE"))
    fun(prepared, spec)
  }, seed = TRUE)
}
