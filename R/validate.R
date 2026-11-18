# Internal input checks, called at the top of the three raw-data entry points.
# Each check stops with a short, specific message instead of letting a bad
# input silently produce invalid fits.

# x must be a 0/1 vector with no NA.
check_binary <- function(x, name) {
  if (anyNA(x) || !all(x %in% c(0, 1))) {
    stop("`", name, "` must be a 0/1 vector with no missing values.", call. = FALSE)
  }
  invisible(TRUE)
}

# Bounds must satisfy 0 <= lo < hi <= 1.
check_bounds <- function(b, name) {
  if (length(b) != 2 || anyNA(b) || b[1] < 0 || b[2] > 1 || b[1] >= b[2]) {
    stop("`", name, "` must be a length-2 vector with 0 <= lower < upper <= 1.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Checks for the fixed-h and internally tilted entry points: (Y, A, W).
# allow_missing_Y = TRUE lets Y contain NA (missing outcomes, Delta = 0); the
# observed Y values must still be 0/1, and each arm needs some observed outcomes.
check_inputs <- function(Y, A, W, gbounds, allow_missing_Y = FALSE) {
  W <- as.data.frame(W)
  if (length(Y) != length(A) || length(Y) != nrow(W)) {
    stop("`Y`, `A` and `W` must have the same number of rows (got ",
         length(Y), ", ", length(A), ", ", nrow(W), ").", call. = FALSE)
  }
  check_binary(A, "A")
  if (allow_missing_Y) {
    obs <- !is.na(Y)
    check_binary(Y[obs], "Y")  # only the observed values need to be 0/1
    for (a in 0:1) {
      if (!any(obs & A == a)) {
        stop("No observed outcomes in treatment arm A = ", a, ".", call. = FALSE)
      }
    }
  } else {
    check_binary(Y, "Y")
  }
  if (anyNA(W)) stop("`W` must not contain missing values. Fill them in first with add_missing_indicators().", call. = FALSE)
  check_bounds(gbounds, "gbounds")
  invisible(TRUE)
}

# Checks for the externally tilted entry point.
# Y, A are internal-only (length n = sum(D == 1)); D and Ws are pooled (length N).
# allow_missing_Y = TRUE lets Y contain NA (see check_inputs()).
check_inputs_external <- function(Y, A, D, Ws, Wg, gbounds, sbounds, allow_missing_Y = FALSE) {
  check_binary(A, "A")
  check_binary(D, "D")
  Ws <- as.data.frame(Ws)
  n <- sum(D == 1)
  if (length(Y) != n || length(A) != n) {
    stop("`Y` and `A` must be internal-only: length sum(D == 1) = ", n,
         " (got ", length(Y), " and ", length(A), ").", call. = FALSE)
  }
  if (allow_missing_Y) {
    obs <- !is.na(Y)
    check_binary(Y[obs], "Y")  # only the observed values need to be 0/1
    for (a in 0:1) {
      if (!any(obs & A == a)) {
        stop("No observed outcomes in treatment arm A = ", a, ".", call. = FALSE)
      }
    }
  } else {
    check_binary(Y, "Y")
  }
  if (nrow(Ws) != length(D)) {
    stop("`Ws` must have one row per entry of `D` (got ", nrow(Ws),
         " rows, length(D) = ", length(D), ").", call. = FALSE)
  }
  if (anyNA(Ws)) stop("`Ws` must not contain missing values. Fill them in first with add_missing_indicators().", call. = FALSE)
  if (!is.null(Wg)) {
    Wg <- as.data.frame(Wg)
    if (nrow(Wg) != n) {
      stop("`Wg` must be internal-only: ", n, " rows (got ", nrow(Wg), ").", call. = FALSE)
    }
    if (!all(names(Wg) %in% names(Ws))) {
      stop("All column names of `Wg` must also be columns of `Ws`.", call. = FALSE)
    }
    if (anyNA(Wg)) stop("`Wg` must not contain missing values. Fill them in first with add_missing_indicators().", call. = FALSE)
  }
  check_bounds(gbounds, "gbounds")
  check_bounds(sbounds, "sbounds")
  invisible(TRUE)
}
