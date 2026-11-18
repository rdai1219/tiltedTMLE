# Results step of the Shiny app, as plain R functions (no shiny code, easy to test).
# They turn a run's result bundle (run_analysis_main() / run_analysis_bootstrap()) into
# the pieces the results view shows: an estimates table, diagnostics, three plots, a
# plain-language summary, and an R script that reproduces the analysis.
# Only summaries are used; the raw data are never part of the bundle.

utils::globalVariables(c("what", "value", "group", "panel", "w", "Estimator", "Estimate", "CI lower", "CI upper", "x"))

# What the estimand is called in sentences.
estimand_label <- function(estimand) {
  switch(estimand,
    fixed = "tilted average treatment effect with your weights h(x)",
    custom = "tilted average treatment effect with your own tilting function",
    ATE = "average treatment effect (ATE)", ATT = "average treatment effect in the treated (ATT)",
    ATC = "average treatment effect in the controls (ATC)", ATO = "average treatment effect in the overlap population (ATO)",
    `ATE-I` = "average treatment effect in the internal population (ATE-I)",
    `ATE-E` = "average treatment effect in the external population (ATE-E)",
    `ATE-O` = "average treatment effect in the overlap population (ATE-O)",
    estimand)
}

# Estimates table: the main table plus the full-refit bootstrap rows when available.
results_table <- function(main, boot = NULL) {
  tab <- main$table
  if (!is.null(boot) && !inherits(boot, "run_error")) {
    est <- tab$Estimate[tab$Estimator == "TMLE"]
    tab <- rbind(tab,
      est_row("TMLE (full-refit bootstrap, normal)", est, boot$se, boot$CI_normal, "SE from resampling"),
      est_row("TMLE (full-refit bootstrap, percentile)", est, boot$se, boot$CI_percentile, "percentile interval"))
  }
  tab
}

# Effective sample size of a set of weights: (sum h)^2 / sum h^2.
effective_sample_size <- function(h) sum(h)^2 / sum(h^2)

# Share of fitted probabilities sitting on a truncation bound.
share_at_bounds <- function(p, bounds) {
  if (is.null(bounds) || all(bounds == c(0, 1))) return(0)
  mean(p <= bounds[1] + 1e-8 | p >= bounds[2] - 1e-8)
}

# Short table of data and fit diagnostics.
diagnostics_table <- function(main) {
  d <- main$diagnostics
  h <- d$h
  rows <- list(
    c("Rows analysed", format(main$n_analysis)),
    c("Outcomes missing", format(main$n_missing_outcomes)),
    c("Effective sample size of the tilting weights",
      sprintf("%.1f of %d (%.0f%%)", effective_sample_size(h), length(h), 100 * effective_sample_size(h) / length(h))),
    c("Treatment propensity at truncation bounds", sprintf("%.1f%%", 100 * share_at_bounds(d$g1W, main$bounds$g))))
  if (!is.null(d$s1W)) {
    rows <- c(rows, list(c("Selection propensity at truncation bounds",
                           sprintf("%.1f%%", 100 * share_at_bounds(d$s1W, main$bounds$s)))))
  }
  if (!is.na(main$converged)) {
    rows <- c(rows, list(c("Targeting converged",
                           paste0(if (main$converged) "yes" else "NO", " (", main$iterations, " iterations)"))))
  }
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("Diagnostic", "Value")
  out
}

# ---- plots (ggplot2) ----

plot_theme <- function() ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())

# Propensity overlap: fitted propensities by group, truncation bounds dashed.
plot_overlap <- function(main) {
  d <- main$diagnostics
  if (main$method == "external") {
    pts <- rbind(
      data.frame(value = d$s1W, group = ifelse(d$D == 1, "Internal", "External"),
                 panel = "Selection propensity s(x)", stringsAsFactors = FALSE),
      data.frame(value = d$g1W, group = ifelse(d$A == 1, "Treated", "Control"),
                 panel = "Treatment propensity g(x), internal rows", stringsAsFactors = FALSE))
    bnd <- rbind(data.frame(panel = "Selection propensity s(x)", x = main$bounds$s),
                 data.frame(panel = "Treatment propensity g(x), internal rows", x = main$bounds$g))
  } else {
    pts <- data.frame(value = d$g1W, group = ifelse(d$A == 1, "Treated", "Control"),
                      panel = "Treatment propensity g(x)", stringsAsFactors = FALSE)
    bnd <- data.frame(panel = "Treatment propensity g(x)", x = main$bounds$g)
  }
  bnd <- bnd[bnd$x > 0 & bnd$x < 1, , drop = FALSE]          # (0, 1) = no truncation, nothing to draw
  p <- ggplot2::ggplot(pts, ggplot2::aes(x = value, fill = group)) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.55, bins = 30, colour = NA) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "Fitted probability", y = "Count", title = "Propensity overlap",
                  subtitle = "Little overlap between groups means the estimate relies on extrapolation") +
    plot_theme()
  if (nrow(bnd) > 0) {
    p <- p + ggplot2::geom_vline(data = bnd, ggplot2::aes(xintercept = x), linetype = "dashed", inherit.aes = FALSE)
  }
  p
}

# Tilting weights relative to their mean, with the effective sample size.
plot_weights <- function(main) {
  h <- main$diagnostics$h
  df <- data.frame(w = h / mean(h))
  ess <- effective_sample_size(h)
  ggplot2::ggplot(df, ggplot2::aes(x = w)) +
    ggplot2::geom_histogram(bins = 30, fill = "#4C78A8", colour = NA) +
    ggplot2::labs(x = "Tilting weight h(x) relative to its mean", y = "Count", title = "Tilting weights",
                  subtitle = sprintf("Effective sample size %.1f of %d (%.0f%%)", ess, length(h), 100 * ess / length(h))) +
    plot_theme()
}

# Forest-style comparison of the estimators (points only where no interval exists).
plot_estimators <- function(main, boot = NULL) {
  tab <- results_table(main, boot)
  tab$Estimator <- factor(tab$Estimator, levels = rev(tab$Estimator))
  ggplot2::ggplot(tab, ggplot2::aes(x = Estimate, y = Estimator)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = `CI lower`, xmax = `CI upper`), orientation = "y", width = 0.2, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.5, colour = "#4C78A8") +
    ggplot2::labs(x = "Estimated effect (difference in outcome probability)", y = NULL,
                  title = "Estimators compared",
                  subtitle = "Benchmarks without intervals are plotted as points") +
    plot_theme()
}

# ---- plain-language summary ----

# Returns list(sentence, caveats) as text.
summary_text <- function(main, boot = NULL) {
  tm <- main$table[main$table$Estimator == "TMLE", ]
  lvl <- format(100 * main$conf_level)
  sentence <- sprintf(
    "Using TMLE, the estimated %s is %.3f (%s%% CI %.3f to %.3f): the difference in outcome probability, treated minus control.",
    estimand_label(main$estimand), tm$Estimate, lvl, tm$`CI lower`, tm$`CI upper`)
  caveats <- c(
    paste("The estimate has a causal interpretation only if there is no unmeasured confounding,",
          "treatment (and, for the external methods, selection) has positive probability for every covariate pattern,",
          "and the models used are adequate."),
    if (main$n_missing_outcomes > 0)
      paste0(main$n_missing_outcomes, " outcome(s) were missing; the analysis assumes they are missing at random given treatment and covariates."),
    if (!is.na(main$converged) && !main$converged)
      "The targeting iteration did not converge; treat the TMLE result with caution.")
  list(sentence = sentence, caveats = caveats)
}

# ---- reproduction script ----

# Text of an R script that reruns the analysis with the package on the user's own file.
# Uses only exported functions, so it works for anyone with tiltedTMLE installed.
reproduction_script <- function(prepared, spec, file_name = "your_data.csv") {
  m <- prepared$mapping
  q <- function(x) encodeString(x, quote = "\"")   # escapes backslashes (Windows paths) and quotes
  vec <- function(x) paste0("c(", paste(q(x), collapse = ", "), ")")
  args_txt <- paste(deparse(spec$args, width.cutoff = 70L), collapse = "\n")

  # the tilting functions: a preset, or the user's own formula (written out as functions)
  tilt_lines <- function(preset_fun) {
    if (spec$estimand != "custom") return(paste0("tilt <- ", preset_fun, "(", q(spec$estimand), ")"))
    v <- spec$custom$var
    c(paste0("tilt <- list(eta = function(", v, ") ", spec$custom$eta, ",   # your tilting function"),
      paste0("             eta_dot = function(", v, ") ", spec$custom$deta_used, ")   # and its derivative"))
  }

  # Read the file the way the app does: blank cells count as missing (NA)
  excel <- tolower(tools::file_ext(file_name)) %in% c("xlsx", "xls")
  read_line <- if (excel) {
    paste0("dat <- as.data.frame(readxl::read_excel(", q(file_name), "))   # <- edit: path to your data file")
  } else {
    paste0("dat <- read.csv(", q(file_name), ", na.strings = c(\"\", \"NA\"))   # <- edit: path to your data file")
  }

  lines <- c(
    "# Reproduces the analysis from the tiltedTMLE app.",
    paste0("# tiltedTMLE ", as.character(utils::packageVersion("tiltedTMLE")), "; ", R.version.string),
    "library(tiltedTMLE)",
    "",
    read_line,
    "")

  if (m$setting == "external") {
    lines <- c(lines,
      paste0("D <- as.integer(as.character(dat[[", q(m$d), "]]) == ", q(m$d_one), ")   # 1 = internal (trial) rows"),
      "idx1 <- which(D == 1)",
      paste0("Y <- ifelse(is.na(dat[[", q(m$y), "]][idx1]), NA, as.integer(as.character(dat[[", q(m$y), "]][idx1]) == ", q(m$y_one), "))"),
      paste0("A <- as.integer(as.character(dat[[", q(m$a), "]][idx1]) == ", q(m$a_one), ")"))
  } else {
    lines <- c(lines,
      paste0("Y <- ifelse(is.na(dat[[", q(m$y), "]]), NA, as.integer(as.character(dat[[", q(m$y), "]]) == ", q(m$y_one), "))"),
      paste0("A <- as.integer(as.character(dat[[", q(m$a), "]]) == ", q(m$a_one), ")"))
  }
  lines <- c(lines,
    paste0("W <- dat[, ", vec(m$covariates), ", drop = FALSE]"),
    "names(W) <- make.names(names(W), unique = TRUE)",
    "W[] <- lapply(W, function(x) if (is.character(x)) factor(x) else x)")
  if (length(m$dropped_constant) > 0) {
    lines <- c(lines, paste0("W <- W[, setdiff(names(W), ", vec(m$dropped_constant), "), drop = FALSE]   # constant columns were dropped"))
  }
  if (isTRUE(m$add_indicators)) lines <- c(lines, "W <- add_missing_indicators(W)   # fill missing covariates, add indicators")
  lines <- c(lines, "", paste0("args <- ", args_txt), "")

  if (spec$method == "fixed") {
    h_line <- paste0("h <- dat[[", q(m$h), "]]   # tilting weights h(x), one per row")
    lines <- c(lines, h_line,
      "fit <- do.call(estimate_tau_h_fixed, c(list(Y, A, W, h), args))")
  } else if (spec$method == "internal") {
    lines <- c(lines, tilt_lines("tilt_spec"),
      "fit <- do.call(estimate_tau_h_internal_tilted, c(list(Y, A, W, tilt$eta, tilt$eta_dot), args))")
  } else {
    lines <- c(lines, tilt_lines("tilt_spec_external"),
      "fit <- do.call(estimate_tau_h_external_tilted, c(list(Y, A, D, W, tilt$eta, tilt$eta_dot), args))")
  }
  lines <- c(lines, "",
    "c(TMLE = fit$TMLE, SE = fit$TMLE_se)",
    paste0("fit$TMLE + c(-1, 1) * qnorm(", 1 - (1 - spec$conf_level) / 2, ") * fit$TMLE_se   # ", 100 * spec$conf_level, "% EIF-based interval"))

  if (spec$bootstrap$on && spec$method != "fixed") {
    seed_line <- if (is.na(spec$bootstrap$seed)) NULL else paste0("set.seed(", spec$bootstrap$seed, ")")
    first <- if (spec$method == "internal") "Y, A, W" else "Y, A, D, W"
    boot_fun <- if (spec$method == "internal") "bootstrap_tau_h_internal_tilted" else "bootstrap_tau_h_external_tilted"
    lines <- c(lines, "", "# full-refit bootstrap", seed_line,
      paste0("boot <- do.call(", boot_fun, ", c(list(", first, ", tilt$eta, tilt$eta_dot, n_boot = ",
             spec$bootstrap$n_boot, ", conf_level = ", spec$conf_level, "), args))"),
      "c(SE = boot$se); boot$CI_normal; boot$CI_percentile")
  }
  paste(lines, collapse = "\n")
}

# ---- HTML report ----

# Render the report template (inst/app/report.Rmd) to `output_file`. Summaries only.
# The template runs in an environment whose parent is the package namespace, so it can
# call the results functions above directly.
render_report <- function(main, boot, script, output_file) {
  template <- system.file("app", "report.Rmd", package = "tiltedTMLE")
  if (template == "") stop("Report template not found.", call. = FALSE)
  work <- tempfile("report_"); dir.create(work)
  file.copy(template, file.path(work, "report.Rmd"))
  output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
  rmarkdown::render(file.path(work, "report.Rmd"),
                    output_file = basename(output_file), output_dir = dirname(output_file),
                    params = list(main = main, boot = boot, script = script),
                    envir = new.env(parent = asNamespace("tiltedTMLE")), quiet = TRUE)
  invisible(output_file)
}
