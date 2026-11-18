# Missing outcomes (assumed missing at random given A and W).
#
# Notation (same as the reference papers):
#   Delta = 1 if the outcome Y is observed, 0 if it is missing.
#         The package infers it from the data: Delta = 1 - is.na(Y).
#   pi_a(W) = P(Delta = 1 | A = a, W)   -- the "missingness propensity score".
#
# How it enters the TMLE:
#   * Q is fit on the rows with Delta = 1 only, then predicted for ALL rows.
#   * In every clever covariate / EIF term, the treatment propensity g_a(W) is
#     replaced by g_a(W) * pi_a(W)  (probability of getting arm a AND being observed),
#     and the residual (Y - Q) is multiplied by Delta.
# With no missing outcomes pi = 1 and everything reduces to the original estimators.

#' Fit the missing-outcome propensity pi_a(W) = P(Delta = 1 | A = a, W)
#'
#' Fits the indicator `Delta` (1 = outcome observed) on treatment and covariates,
#' by logistic regression or [SuperLearner::SuperLearner()], and returns the
#' predicted probability of being observed under `A = 1` and `A = 0`.
#'
#' @param Delta Binary (`0`/`1`) vector, `1` if the outcome is observed.
#' @param A Binary (`0`/`1`) treatment vector (observed for every row).
#' @param W Covariate data frame (or matrix, coerced via [as.data.frame()]).
#' @param method `"glm"` (default) or `"SL"`; see [fit_Q()].
#' @param piform Formula, used only when `method = "glm"`. Defaults to main
#'   effects, `Delta ~ A + W1 + W2 + ...`.
#' @param SL.library SuperLearner library, used only when `method = "SL"`.
#' @param pibounds Length-2 truncation bounds for the fitted probabilities,
#'   default `c(0.025, 1)` (guards against huge `1/pi` weights).
#' @param ... Passed through to `glm()` or `SuperLearner::SuperLearner()`.
#'
#' @return A list with components:
#' \item{pi1W}{Fitted `P(Delta = 1 | A = 1, W)`, truncated to `pibounds`.}
#' \item{pi0W}{Fitted `P(Delta = 1 | A = 0, W)`, truncated to `pibounds`.}
#' \item{piAW}{Fitted probability at the observed treatment.}
#' \item{fit}{The underlying fitted model object.}
#' \item{method}{`"glm"` or `"SL"`.}
#'
#' @seealso [estimate_tau_h_fixed()], which calls this when `Y` has `NA`s.
#' @examples
#' set.seed(1)
#' n <- 300
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, 0.5)
#' Delta <- rbinom(n, 1, stats::plogis(1 + 0.5 * W$X1 - 0.3 * A))
#' pifit <- fit_pi(Delta, A, W)
#' summary(pifit$pi1W)
#' @export
fit_pi <- function(Delta, A, W, method = c("glm", "SL"), piform = NULL,
                   SL.library = c("SL.glm"), pibounds = c(0.025, 1), ...) {
  method <- match.arg(method)
  W <- as.data.frame(W)

  if (method == "glm") {
    if (is.null(piform)) piform <- stats::as.formula(paste("Delta ~ A +", paste(names(W), collapse = "+")))
    fit <- stats::glm(piform, data = data.frame(Delta = Delta, A = A, W),
                      family = stats::binomial(), ...)
  } else {
    if (!requireNamespace("SuperLearner", quietly = TRUE)) {
      stop("Package 'SuperLearner' is required for method = \"SL\". ",
           "Install it with install.packages(\"SuperLearner\").", call. = FALSE)
    }
    fit <- run_superlearner(Y = Delta, X = data.frame(A = A, W),
                                      SL.library = SL.library, ...)
  }

  # Predict at A = 1 and A = 0 for ALL rows (predict_Qhat handles glm and SL;
  # c(0, 1) bounds there = no clamping, we truncate with pibounds below).
  pi1W <- clamp_bounds(predict_Qhat(fit, 1, W, Qbounds = c(0, 1)), pibounds)
  pi0W <- clamp_bounds(predict_Qhat(fit, 0, W, Qbounds = c(0, 1)), pibounds)
  piAW <- A * pi1W + (1 - A) * pi0W

  list(pi1W = pi1W, pi0W = pi0W, piAW = piAW, fit = fit, method = method)
}

# Internal helper (used by the raw-data entry points): fit Q on the rows with an
# OBSERVED outcome (Delta = 1), then predict it for ALL rows, so that
# Q(1, W), Q(0, W) and Q(A, W) exist even where Y is missing.
# Returns a list shaped like fit_Q(): QAW, Q1W, Q0W, fit, method.
fit_Q_observed <- function(Y, A, W, method, Qform, SL.library, Q_args = list()) {
  obs <- which(!is.na(Y))
  Qfit_obs <- do.call(fit_Q, c(list(Y[obs], A[obs], W[obs, , drop = FALSE], method = method,
                                    Qform = Qform, SL.library = SL.library), Q_args))
  Q1W <- predict_Qhat(Qfit_obs$fit, 1, W)
  Q0W <- predict_Qhat(Qfit_obs$fit, 0, W)
  list(QAW = A * Q1W + (1 - A) * Q0W, Q1W = Q1W, Q0W = Q0W,
       fit = Qfit_obs$fit, method = Qfit_obs$method)
}

# Internal helper: fit pi_a(W) = P(Delta = 1 | A = a, W) only when some outcomes
# are missing; returns NULL otherwise (meaning pi = 1 everywhere).
fit_pi_if_missing <- function(Y, A, W, method, piform, SL.library, pibounds, pi_args = list()) {
  Delta <- as.integer(!is.na(Y))
  if (all(Delta == 1)) return(NULL)
  do.call(fit_pi, c(list(Delta, A, W, method = method, piform = piform,
                         SL.library = SL.library, pibounds = pibounds), pi_args))
}
