# Missing baseline covariates: the "indicator method".
#
# For every covariate with missing values we
#   (1) fill the gaps with a simple value, and
#   (2) add a 0/1 column saying where the gaps were.
# The package's estimators then run on the augmented covariate set W*,
# exactly as in Sec 5.1 of Fang, He, and Dai (2026) / Fang (2024).

# Internal helper: the most frequent non-missing value of x (first one on ties).
# Keeps the type of x (numeric stays numeric, factor stays factor, ...).
most_common <- function(x) {
  vals <- x[!is.na(x)]
  u <- unique(vals)
  u[which.max(tabulate(match(vals, u)))]
}

# Internal helper: TRUE if x is a "binary" covariate, i.e. has exactly 2
# distinct non-missing values (0/1 numeric, logical, 2-level factor, ...).
is_binary_cov <- function(x) {
  length(unique(x[!is.na(x)])) == 2
}

#' Fill in missing covariates and add missingness indicators
#'
#' Implements the indicator method for missing baseline covariates used in
#' the reference paper (Sec 5.1). For each column of `W` that has at least one
#' `NA`:
#' \itemize{
#'   \item the missing entries are imputed: a continuous (numeric,
#'     non-binary) column by its sample mean, and a binary column (or any
#'     other non-numeric column, e.g. a factor) by its most common value;
#'   \item a new 0/1 column named `<column>_miss` is appended, equal to `1`
#'     where the original value was missing.
#' }
#' Columns with no missing values are returned unchanged and get no
#' indicator (a constant indicator would break the default `glm()` formulas).
#'
#' Call this on the covariates BEFORE passing them to
#' [estimate_tau_h_fixed()], [estimate_tau_h_internal_tilted()] or
#' [estimate_tau_h_external_tilted()], which stop with an error if given `NA`s.
#' Only covariates are handled; missing outcomes or treatments are not.
#'
#' For the external setting, call this once on the full pooled `Ws` (so
#' internal and external rows share the same imputed values and indicators).
#'
#' @param W Covariate data frame (or matrix, coerced via [as.data.frame()]).
#' @param suffix Suffix for the indicator column names, default `"_miss"`.
#'
#' @return A data frame with the same rows as `W`: the original columns
#'   (missing values filled in) followed by one indicator column per
#'   covariate that had missing values.
#'
#' @references Fang (2024), "A sensitivity analysis to handle missing
#'   covariates", Journal of Causal Inference 10(1): 372-414.
#' @seealso [estimate_tau_h_fixed()], [estimate_tau_h_internal_tilted()],
#'   [estimate_tau_h_external_tilted()].
#' @examples
#' W <- data.frame(age = c(50, NA, 61, 45), sex = c(1, 0, NA, 1), bmi = c(22, 25, 31, 28))
#' add_missing_indicators(W)
#' @export
add_missing_indicators <- function(W, suffix = "_miss") {
  W <- as.data.frame(W)
  has_na <- vapply(W, anyNA, logical(1))   # which columns have any NA
  miss_cols <- names(W)[has_na]

  # Checks: nothing to impute from, or a name clash with the new indicators
  all_na <- vapply(W[miss_cols], function(x) all(is.na(x)), logical(1))
  if (any(all_na)) {
    stop("Column(s) entirely missing, cannot impute: ",
         paste(miss_cols[all_na], collapse = ", "), call. = FALSE)
  }
  new_names <- paste0(miss_cols, suffix)
  if (any(new_names %in% names(W))) {
    stop("Indicator column name(s) already in `W`: ",
         paste(intersect(new_names, names(W)), collapse = ", "),
         ". Use a different `suffix`.", call. = FALSE)
  }

  indicators <- list()
  for (col in miss_cols) {
    x <- W[[col]]
    is_miss <- is.na(x)

    # (1) impute: mean for continuous, most common value otherwise
    fill <- if (is.numeric(x) && !is_binary_cov(x)) mean(x, na.rm = TRUE) else most_common(x)
    x[is_miss] <- fill
    W[[col]] <- x

    # (2) indicator: 1 where the value was originally missing
    indicators[[paste0(col, suffix)]] <- as.integer(is_miss)
  }

  if (length(indicators) > 0) W <- cbind(W, as.data.frame(indicators))
  W
}
