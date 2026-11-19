# Continuous outcomes.
#
# All the TMLE machinery in this package works with an outcome in [0, 1]: the
# outcome regression Q is modelled on the logit scale and every targeting step is
# a logistic fluctuation. A continuous outcome is handled the standard way
# (Fang, He and Dai, 2026, Sec. 3.3.2; Gruber and van der Laan, 2012):
#
#   1. Choose bounds [a, b] containing every observed Y (default: the observed range).
#   2. Rescale: Y* = (Y - a) / (b - a), which lies in [0, 1].
#   3. Run the usual estimators on Y*. Q is fit with a quasi-binomial glm (or a
#      gaussian SuperLearner, clamped into [0, 1]), and the logistic fluctuation
#      uses a quasi-binomial working model, which is valid for any Y* in [0, 1].
#   4. Map back: every effect estimate, standard error and confidence limit is
#      multiplied by (b - a). The shift a cancels because tau is a difference of means.
#
# A binary 0/1 outcome gets bounds c(0, 1), so steps 2 and 4 leave it unchanged.

# Internal: TRUE if every non-missing value of Y is 0 or 1.
is_binary_outcome <- function(Y) all(Y[!is.na(Y)] %in% c(0, 1))

# Internal: resolve the outcome bounds c(a, b).
# Ybounds = NULL gives c(0, 1) for a binary outcome and the observed range otherwise.
outcome_bounds <- function(Y, Ybounds = NULL) {
  y <- Y[!is.na(Y)]
  if (is.null(Ybounds)) {
    Ybounds <- if (is_binary_outcome(y)) c(0, 1) else range(y)
  }
  if (length(Ybounds) != 2 || !all(is.finite(Ybounds)) || Ybounds[1] >= Ybounds[2]) {
    stop("`Ybounds` must be a length-2 vector of finite values with lower < upper ",
         "(a constant outcome cannot be analysed).", call. = FALSE)
  }
  if (any(y < Ybounds[1] | y > Ybounds[2])) {
    stop("All observed `Y` values must lie within `Ybounds` = [",
         Ybounds[1], ", ", Ybounds[2], "].", call. = FALSE)
  }
  Ybounds
}

# Internal: map Y from [a, b] onto [0, 1] (NA stays NA).
rescale_outcome <- function(Y, Ybounds) (Y - Ybounds[1]) / (Ybounds[2] - Ybounds[1])

# Internal: map the result of an estimate_tau_h_*_given_fits() call from the
# [0, 1] scale back to the original outcome scale. Effect estimates, standard
# errors and confidence limits are multiplied by (b - a). The fitted nuisance
# objects (Qfit, tmle_fit, ...) stay on the [0, 1] scale. Adds `Ybounds`.
unscale_results <- function(out, Ybounds) {
  k <- Ybounds[2] - Ybounds[1]
  fields <- c("MLE", "IPW", "AIPW", "AIPW_se", "TMLE", "TMLE_se", "TMLE_CI",
              "Bootstrap_SE", "Bootstrap_CI_normal", "Bootstrap_CI_percentile")
  for (f in intersect(fields, names(out))) out[[f]] <- k * out[[f]]
  out$Ybounds <- Ybounds
  out
}
