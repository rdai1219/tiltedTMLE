# Custom tilting functions for the Shiny app, as plain R functions (no shiny code, easy to test).
#
# The estimators take a tilting function eta (of the treatment propensity g for the internal
# method, or of the selection propensity s for the external method) and its derivative eta_dot.
# The four named estimands (ATE, ATT, ATC, ATO) are just presets; this file lets the user type
# their own eta, for example  g * (1 - g)  for the overlap weights.
#
# Because the text is typed by a user, it is NOT evaluated as ordinary R code. It must use only
# the one variable (g or s), numbers and a short list of arithmetic / math functions; anything
# else is rejected before evaluation.

# Functions a typed formula may use.
TILT_ALLOWED <- c("+", "-", "*", "/", "^", "(", "exp", "log", "sqrt", "abs", "pmin", "pmax",
                  "ifelse", "<", ">", "<=", ">=", "==", "!=", "&", "|", "!", "plogis")

# An environment holding only the allowed functions (and nothing else) for evaluating a formula.
tilt_eval_env <- function() {
  env <- new.env(parent = emptyenv())
  for (nm in setdiff(TILT_ALLOWED, "plogis")) assign(nm, get(nm, envir = baseenv()), envir = env)
  assign("plogis", stats::plogis, envir = env)
  env
}

# Read a formula typed as text; stop with a readable message unless it only uses `var`,
# numbers and allowed functions.
parse_tilt_expr <- function(text, var, label) {
  text <- trimws(text)
  if (length(text) != 1 || is.na(text) || !nzchar(text)) {
    stop("Enter a formula for ", label, " (for example ", var, " * (1 - ", var, ")).", call. = FALSE)
  }
  expr <- tryCatch(str2lang(text), error = function(e) {
    stop("Cannot read the formula for ", label, ": ", conditionMessage(e), call. = FALSE)
  })
  bad_vars <- setdiff(all.vars(expr), var)
  if (length(bad_vars) > 0) {
    stop("The formula for ", label, " may only use the variable '", var, "' (found: ",
         paste(bad_vars, collapse = ", "), ").", call. = FALSE)
  }
  funs <- setdiff(all.names(expr, functions = TRUE, unique = TRUE), all.vars(expr))
  bad_funs <- setdiff(funs, TILT_ALLOWED)
  if (length(bad_funs) > 0) {
    stop("The formula for ", label, " uses something that is not allowed: ", paste(bad_funs, collapse = ", "),
         ". Allowed: ", paste(setdiff(TILT_ALLOWED, c("(", "plogis")), collapse = " "), " ( ) plogis.", call. = FALSE)
  }
  expr
}

# Turn a checked expression into a function of `var` that always returns one number per input.
tilt_function <- function(expr, var) {
  env <- tilt_eval_env()
  force(expr); force(var)
  function(x) {
    vals <- new.env(parent = env)
    assign(var, x, envir = vals)
    rep_len(as.numeric(eval(expr, vals)), length(x))
  }
}

# Build eta and eta_dot from typed text.
#   eta_text   formula for eta, e.g. "g * (1 - g)"
#   deta_text  formula for its derivative; blank = differentiate eta automatically (stats::D)
#   var        "g" (internal method) or "s" (external method)
# Returns list(eta, eta_dot, eta_text, deta_text, deta_auto, var). Stops with a readable message if
# a formula is invalid, not finite on (0, 1), or eta is negative.
make_tilt <- function(eta_text, deta_text = "", var = "g") {
  eta_expr <- parse_tilt_expr(eta_text, var, paste0("eta(", var, ")"))
  deta_auto <- is.null(deta_text) || is.na(deta_text) || !nzchar(trimws(deta_text))
  if (deta_auto) {
    deta_expr <- tryCatch(stats::D(eta_expr, var), error = function(e) {
      stop("The derivative of eta(", var, ") could not be computed automatically (", conditionMessage(e),
           "). Please type it in the derivative box.", call. = FALSE)
    })
    deta_expr <- parse_tilt_expr(paste(deparse(deta_expr, width.cutoff = 500L), collapse = ""), var, "the automatic derivative")
  } else {
    deta_expr <- parse_tilt_expr(deta_text, var, paste0("the derivative of eta(", var, ")"))
  }
  eta <- tilt_function(eta_expr, var)
  eta_dot <- tilt_function(deta_expr, var)

  # sanity checks on a grid of probabilities
  grid <- seq(0.01, 0.99, by = 0.01)
  e <- suppressWarnings(eta(grid)); d <- suppressWarnings(eta_dot(grid))   # NaN warnings become the messages below
  if (anyNA(e) || any(!is.finite(e))) stop("eta(", var, ") is not finite for some ", var, " in (0, 1).", call. = FALSE)
  if (any(e < -1e-12)) stop("eta(", var, ") must not be negative: it is a weight for each subject.", call. = FALSE)
  if (all(e == 0)) stop("eta(", var, ") is zero everywhere, so no subject would get any weight.", call. = FALSE)
  if (anyNA(d) || any(!is.finite(d))) stop("The derivative of eta(", var, ") is not finite for some ", var, " in (0, 1).", call. = FALSE)

  list(eta = eta, eta_dot = eta_dot, eta_text = trimws(eta_text),
       deta_text = if (deta_auto) paste(deparse(deta_expr, width.cutoff = 500L), collapse = "") else trimws(deta_text),
       deta_auto = deta_auto, var = var)
}

# Plot eta and its derivative over (0, 1), so the user can see what they typed.
plot_tilt <- function(tilt) {
  grid <- seq(0.01, 0.99, by = 0.01)
  df <- rbind(data.frame(x = grid, value = tilt$eta(grid), what = paste0("eta(", tilt$var, ")")),
              data.frame(x = grid, value = tilt$eta_dot(grid), what = paste0("derivative of eta(", tilt$var, ")")))
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = value)) +
    ggplot2::geom_line(colour = "#4C78A8") +
    ggplot2::facet_wrap(~what, scales = "free_y") +
    ggplot2::labs(x = if (tilt$var == "g") "treatment propensity g" else "selection propensity s", y = NULL) +
    ggplot2::theme_minimal(base_size = 9)
}
