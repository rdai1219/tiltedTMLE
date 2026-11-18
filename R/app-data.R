# Data step of the Shiny app, as plain R functions (no shiny code here, so they are
# easy to test). The app's data module calls prepare_analysis_data() and passes the
# result on to the estimators.
#
# What it does, in order:
#   1. recode binary columns to 0/1 (the user says which value counts as 1),
#   2. for the external setting, split one pooled table into the objects the package
#      expects (Y, A for internal rows only; D and Ws for all rows),
#   3. clean covariates (text -> factor, safe names, drop constant columns),
#   4. handle missing covariates (add indicators, or stop),
#   5. build a short summary the app can display (counts only, never raw rows).
#
# Missing OUTCOMES need no special step: they stay NA in Y and the estimators handle them.

# Distinct non-missing values of x as text, sorted: what the app offers as
# "the value that counts as 1".
binary_levels <- function(x) sort(unique(as.character(x[!is.na(x)])))

# Recode x to 0/1: the value equal to `one` (compared as text) becomes 1, the other
# value becomes 0, NA stays NA. Stops unless x has exactly 2 distinct values.
recode_binary <- function(x, one, name) {
  lev <- binary_levels(x)
  if (length(lev) != 2) {
    stop("`", name, "` must have exactly 2 distinct non-missing values (found ",
         length(lev), ").", call. = FALSE)
  }
  if (is.null(one) || !(one %in% lev)) {
    stop("For `", name, "`, choose which value counts as 1 (one of: ",
         paste(lev, collapse = ", "), ").", call. = FALSE)
  }
  ifelse(is.na(x), NA_integer_, as.integer(as.character(x) == one))
}

# Turn an uploaded table into analysis-ready objects.
#
# data        data.frame, one row per subject.
# setting     "fixed", "internal" or "external" (which entry point will be used).
# y, a        column names of the outcome and the treatment; y_one, a_one = value that is 1.
# covariates  column names of the covariates.
# d, d_one    external setting only: column marking internal (1) vs external (0) rows.
#             Y and A are read from the internal rows only; their values on external
#             rows are ignored.
# add_indicators  TRUE: fill missing covariates and add <name>_miss columns
#                 (see add_missing_indicators()); FALSE: stop if any are missing.
# h           fixed setting only, REQUIRED: the column holding the pre-specified tilting weights
#             h(x), one non-negative number per row (e.g. all 1s gives the plain ATE).
#
# Returns a list: setting, Y, A, (W for fixed/internal | Ws and D for external),
# h (fixed setting only: the weight vector), summary (counts only) and notes
# (messages to show the user).
prepare_analysis_data <- function(data, setting = c("fixed", "internal", "external"),
                                  y, a, covariates, y_one, a_one,
                                  d = NULL, d_one = NULL, add_indicators = TRUE, h = NULL) {
  setting <- match.arg(setting)
  data <- as.data.frame(data)
  notes <- character()
  if (setting != "fixed") h <- NULL            # a weight column only applies to the fixed setting

  # ---- 0. basic column checks ----
  roles <- c(y, a, covariates, if (setting == "external") d, h)
  absent <- setdiff(roles, names(data))
  if (length(absent) > 0) stop("Column(s) not found: ", paste(absent, collapse = ", "), call. = FALSE)
  if (length(covariates) == 0) stop("Choose at least one covariate.", call. = FALSE)
  if (anyDuplicated(roles)) stop("A column can only have one role (outcome, treatment, covariate, D or weights).", call. = FALSE)
  if (setting == "fixed" && is.null(h)) {
    stop("For the fixed-h method, choose the column that holds the tilting weights h(x): one non-negative ",
         "number per row (for example, 1 for everyone gives the plain ATE).", call. = FALSE)
  }

  # ---- 1. which rows carry Y and A? ----
  if (setting == "external") {
    if (is.null(d)) stop("Choose the column that marks internal vs external rows.", call. = FALSE)
    D <- recode_binary(data[[d]], d_one, "D")
    if (anyNA(D)) stop("`D` must not have missing values.", call. = FALSE)
    if (!all(c(0, 1) %in% D)) stop("`D` must contain both internal and external rows.", call. = FALSE)
    idx1 <- which(D == 1)                    # internal rows
  } else {
    D <- NULL
    idx1 <- seq_len(nrow(data))              # all rows
  }
  Y <- recode_binary(data[[y]][idx1], y_one, "Y")   # NA allowed = missing outcome
  A <- recode_binary(data[[a]][idx1], a_one, "A")
  if (anyNA(A)) stop("`A` (treatment) must not have missing values.", call. = FALSE)

  # ---- 1b. fixed setting: optional pre-specified weights h(x) ----
  h_vec <- NULL
  if (!is.null(h)) {
    h_vec <- data[[h]]
    if (!is.numeric(h_vec)) stop("The weight column must be numeric.", call. = FALSE)
    if (anyNA(h_vec) || any(h_vec < 0) || sum(h_vec) <= 0) {
      stop("The weight column must have no missing values, no negative values, and a positive sum.",
           call. = FALSE)
    }
  }

  # ---- 2. covariates: text -> factor, safe names, drop constant columns ----
  W <- data[covariates]
  names(W) <- make.names(names(W), unique = TRUE)   # spaces etc. would break model formulas
  renamed <- !identical(names(W), covariates)
  if (renamed) {
    notes <- c(notes, "Some covariate names were changed to valid R names.")
  }
  # the models use variables called A, Y and Delta internally; a covariate with one of these
  # names would clash, so rename it
  clash <- names(W) %in% c("A", "Y", "Delta")
  if (any(clash)) {
    notes <- c(notes, paste0("Covariate(s) ", paste(names(W)[clash], collapse = ", "),
                             " renamed (added _cov) because the name is used internally."))
    names(W)[clash] <- paste0(names(W)[clash], "_cov")
    renamed <- TRUE
  }
  W[] <- lapply(W, function(x) if (is.character(x)) factor(x) else x)
  constant <- vapply(W, function(x) length(unique(x[!is.na(x)])) <= 1, logical(1))
  dropped_constant <- names(W)[constant]
  if (any(constant)) {
    notes <- c(notes, paste0("Dropped constant covariate(s): ", paste(names(W)[constant], collapse = ", "), "."))
    W <- W[!constant]
  }
  if (ncol(W) == 0) stop("No usable covariates left (all were constant).", call. = FALSE)

  # External setting: Q is fit on internal rows and predicted for all rows, so a
  # factor level seen only in external rows cannot be predicted.
  if (setting == "external") {
    for (nm in names(W)[vapply(W, is.factor, logical(1))]) {
      unseen <- setdiff(unique(as.character(W[[nm]][-idx1])), unique(as.character(W[[nm]][idx1])))
      unseen <- unseen[!is.na(unseen)]
      if (length(unseen) > 0) {
        stop("Covariate `", nm, "` has value(s) only in external rows (",
             paste(unseen, collapse = ", "), "); merge or remove them.", call. = FALSE)
      }
    }
  }

  # ---- 3. missing covariates ----
  miss <- colSums(is.na(W))
  miss <- miss[miss > 0]
  if (length(miss) > 0) {
    if (add_indicators) {
      W <- add_missing_indicators(W)
      notes <- c(notes, paste0("Filled missing covariates and added indicator columns for: ",
                               paste(names(miss), collapse = ", "), "."))
    } else {
      stop("Covariates have missing values (", paste(names(miss), collapse = ", "),
           "). Turn on missingness indicators or fix the file.", call. = FALSE)
    }
  }

  # ---- 4. summary: counts only, never raw rows ----
  n_missing_Y <- sum(is.na(Y))
  if (n_missing_Y > 0) {
    notes <- c(notes, paste0(n_missing_Y, " missing outcome(s): assumed missing at random given treatment ",
                             "and covariates; the estimators adjust for this automatically."))
  }
  summary <- list(
    setting = setting,
    n_rows = nrow(data),
    n_analysis = length(Y),                              # internal rows (or all rows)
    n_external = if (setting == "external") sum(D == 0) else NA_integer_,
    n_treated = sum(A == 1), n_control = sum(A == 0),
    n_outcome_missing = n_missing_Y,
    outcome_missing_by_arm = c(control = sum(is.na(Y[A == 0])), treated = sum(is.na(Y[A == 1]))),
    outcome_mean_observed = mean(Y, na.rm = TRUE),
    h = if (is.null(h_vec)) NULL else list(min = min(h_vec), mean = mean(h_vec), max = max(h_vec),
                                           ess = effective_sample_size(h_vec)),
    covariates = names(W),
    covariate_missing = miss)

  # how the file was mapped: used to write a script that reproduces the analysis
  mapping <- list(setting = setting, y = y, y_one = y_one, a = a, a_one = a_one, covariates = covariates,
                  d = d, d_one = d_one, h = h, add_indicators = add_indicators,
                  dropped_constant = dropped_constant)
  out <- list(setting = setting, Y = Y, A = A, summary = summary, notes = notes, mapping = mapping)
  if (setting == "external") {
    out$D <- D
    out$Ws <- W                    # pooled covariates, all rows (internal rows are which(D == 1))
  } else {
    out$W <- W
  }
  if (setting == "fixed") out$h <- h_vec
  out
}
