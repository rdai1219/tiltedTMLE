#' Estimate a tilted ATE for a fixed h given already-fitted Q and g
#'
#' Core estimation step given ALREADY-FITTED `Q` and `g` (as returned by
#' [fit_Q()]/[fit_g()]) plus an `h(x)` specification. [estimate_tau_h_fixed()]
#' calls this after fitting `Q`/`g` from raw data; exposed separately so
#' that anything re-using the same `Q`/`g` fit across multiple `h`'s
#' (e.g. a simulation study looping over h-choices under one
#' nuisance-model scenario) doesn't have to refit.
#'
#' `h` can be specified either way:
#' \itemize{
#'   \item a numeric vector of length `N`: a FIXED, precomputed tilting
#'     function (e.g. a prespecified `h(x)`, or `h = 1` for the ATE).
#'   \item a function of `(g1W, g0W)`: a function of the FITTED
#'     propensity score (e.g. `function(g1W, g0W) g1W` for the ATT) --
#'     evaluated AFTER `g` is fit, so it automatically tracks whatever
#'     `g` (correct, misspecified, GLM, or SL) was actually estimated.
#' }
#'
#' @param Y Outcome vector.
#' @param A Binary treatment vector.
#' @param h Tilting-function specification; see Details.
#' @param Qfit A fitted-Q object, as returned by [fit_Q()].
#' @param gfit A fitted-g object, as returned by [fit_g()].
#' @param n_boot Number of target-step bootstrap resamples (see [tmle_bootstrap_h()]).
#' @param boot_strata Optional stratification vector for the bootstrap
#'   (defaults to `A`).
#'
#' @return A list with components `h`, `MLE`, `IPW`, `AIPW`, `AIPW_se`,
#'   `TMLE`, `TMLE_se`, `TMLE_CI`, `Bootstrap_SE`,
#'   `Bootstrap_CI_normal`, `Bootstrap_CI_percentile`, plus the fitted
#'   objects `Qfit`, `gfit`, `tmle_fit`, `boot_fit`.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' Qfit <- fit_Q(Y, A, W)
#' gfit <- fit_g(A, W)
#' out <- estimate_tau_h_fixed_given_fits(Y, A, h = rep(1, n), Qfit, gfit, n_boot = 50)
#' out$TMLE
#' out$TMLE_CI
#' @export
estimate_tau_h_fixed_given_fits <- function(Y, A, h, Qfit, gfit, n_boot = 300, boot_strata = NULL) {
  h_vec <- if (is.function(h)) h(gfit$g1W, gfit$g0W) else h
  stopifnot(length(h_vec) == length(Y))
  strata <- if (is.null(boot_strata)) A else boot_strata

  mle_est <- mle_h(h_vec, Qfit$Q1W, Qfit$Q0W)
  ipw_est <- ipw_h(Y, A, h_vec, gfit$g1W)
  aipw_out <- aipw_h(Y, A, h_vec, Qfit$Q1W, Qfit$Q0W, gfit$g1W)
  tmle_out <- tmle_h(Y, A, h_vec, Qfit$QAW, Qfit$Q1W, Qfit$Q0W, gfit$g1W)
  boot_out <- tmle_bootstrap_h(Y, A, h_vec, Qfit$Q1W, Qfit$Q0W, gfit$g1W,
                                tmle_fit = tmle_out, n_boot = n_boot, strata = strata)

  list(h = h_vec,
       MLE = mle_est,
       IPW = ipw_est,
       AIPW = aipw_out$estimate, AIPW_se = aipw_out$se,
       TMLE = tmle_out$estimate, TMLE_se = tmle_out$se, TMLE_CI = tmle_out$CI,
       Bootstrap_SE = boot_out$se,
       Bootstrap_CI_normal = boot_out$CI_normal,
       Bootstrap_CI_percentile = boot_out$CI_percentile,
       Qfit = Qfit, gfit = gfit, tmle_fit = tmle_out, boot_fit = boot_out)
}

#' Estimate a tilted ATE for a fixed h from raw data
#'
#' The single entry point: given raw data `(Y, A, W)`, an `h(x)`
#' specification (fixed vector, or `function(g1W, g0W)`), and choices for
#' how `Q`/`g` are fit, this fits `Q` and `g`, builds `h`, and returns
#' MLE/IPW/AIPW/TMLE point estimates plus TMLE's EIF-Wald and
#' target-step bootstrap inference.
#'
#' @param Y Outcome vector.
#' @param A Binary treatment vector.
#' @param W Covariate data frame.
#' @param h Tilting-function specification; see [estimate_tau_h_fixed_given_fits()].
#' @param Q_method,Qform,Q_SL.library Passed to [fit_Q()].
#' @param g_method,gform,g_SL.library Passed to [fit_g()].
#' @param truncate If `TRUE` (default), `gbounds` is used as given
#'   (default `c(0.025, 0.975)`). If `FALSE`, `gbounds` is forced to
#'   `c(0, 1)`, i.e. no truncation.
#' @param gbounds Propensity-score truncation bounds; see `truncate`.
#' @param n_boot,boot_strata Passed to [estimate_tau_h_fixed_given_fits()].
#' @param ... Passed through to [fit_Q()]/[fit_g()] (e.g. `cvControl` for SL).
#'
#' @return See [estimate_tau_h_fixed_given_fits()].
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#'
#' ## ATE: h(x) = 1
#' out_ate <- estimate_tau_h_fixed(Y, A, W, h = rep(1, n), n_boot = 50)
#' out_ate$TMLE
#'
#' ## ATT: h(x) = fitted propensity score itself
#' out_att <- estimate_tau_h_fixed(Y, A, W, h = function(g1W, g0W) g1W, n_boot = 50)
#' out_att$TMLE
#' @export
estimate_tau_h_fixed <- function(Y, A, W, h,
                                  Q_method = c("glm", "SL"), Qform = NULL, Q_SL.library = c("SL.glm"),
                                  g_method = c("glm", "SL"), gform = NULL, g_SL.library = c("SL.glm"),
                                  truncate = TRUE, gbounds = c(0.025, 0.975),
                                  n_boot = 300, boot_strata = NULL, ...) {
  Q_method <- match.arg(Q_method)
  g_method <- match.arg(g_method)
  if (!truncate) gbounds <- c(0, 1)

  Qfit <- fit_Q(Y, A, W, method = Q_method, Qform = Qform, SL.library = Q_SL.library, ...)
  gfit <- fit_g(A, W, method = g_method, gform = gform, SL.library = g_SL.library,
                gbounds = gbounds, ...)

  estimate_tau_h_fixed_given_fits(Y, A, h, Qfit, gfit, n_boot = n_boot, boot_strata = boot_strata)
}
