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
#' @param Y Binary (`0`/`1`) outcome vector; see [fit_Q()] Details.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param h Tilting-function specification; see Details.
#' @param Qfit A fitted-Q object, as returned by [fit_Q()].
#' @param gfit A fitted-g object, as returned by [fit_g()].
#' @param n_boot Number of target-step bootstrap resamples (see [tmle_bootstrap_h()]).
#' @param boot_strata Optional stratification vector for the bootstrap
#'   (defaults to `A`).
#' @param warn_fn_h If `TRUE` (default), warn when `h` is a function of the
#'   fitted propensity score. The fixed-h EIF treats `h` as known, so its
#'   SE/CI ignore the extra `eta_dot` term and are not valid for such `h`
#'   (e.g. ATT/ATC/ATO); use [estimate_tau_h_internal_tilted()] for those.
#' @param pifit Optional fitted missing-outcome propensity, as returned by [fit_pi()].
#'   `NULL` (default) means no missing outcomes (`Y` has no `NA`).
#'
#' @return A list with components:
#' \item{h}{The resolved numeric tilting-function weight vector.}
#' \item{MLE}{Plug-in (g-computation) point estimate; see [mle_h()].}
#' \item{IPW}{Inverse-probability-weighted point estimate; see [ipw_h()].}
#' \item{AIPW,AIPW_se}{Augmented IPW point estimate and its EIF-based
#'   standard error; see [aipw_h()].}
#' \item{TMLE,TMLE_se,TMLE_CI}{TMLE point estimate, EIF-based standard
#'   error, and Wald 95% confidence interval; see [tmle_h()].}
#' \item{Bootstrap_SE}{Target-step bootstrap standard error; see [tmle_bootstrap_h()].}
#' \item{Bootstrap_CI_normal}{Wald 95% confidence interval using the bootstrap SE.}
#' \item{Bootstrap_CI_percentile}{Percentile-based 95% confidence interval
#'   from the bootstrap distribution.}
#' \item{Qfit,gfit}{The fitted-Q and fitted-g objects passed in.}
#' \item{tmle_fit}{The full [tmle_h()] result.}
#' \item{boot_fit}{The full [tmle_bootstrap_h()] result.}
#'
#' @seealso [estimate_tau_h_fixed()], the raw-data entry point that calls
#'   this after fitting `Q`/`g`.
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
estimate_tau_h_fixed_given_fits <- function(Y, A, h, Qfit, gfit, n_boot = 300, boot_strata = NULL,
                                             warn_fn_h = TRUE, pifit = NULL) {
  if (is.function(h) && warn_fn_h) {
    warning("`h` is a function of the fitted propensity score, but the fixed-h standard error ",
            "treats h as known and ignores the eta_dot term, so SE/CI may be too narrow. ",
            "For ATT/ATC/ATO use estimate_tau_h_internal_tilted(). ",
            "Set warn_fn_h = FALSE to silence this warning.", call. = FALSE)
  }
  h_vec <- if (is.function(h)) h(gfit$g1W, gfit$g0W) else h
  stopifnot(length(h_vec) == length(Y))
  strata <- if (is.null(boot_strata)) A else boot_strata

  # Missing outcomes: Delta = 1 where Y is observed; pi = P(Delta = 1 | A, W).
  # Without a pifit (no missing outcomes) pi = 1 and nothing changes.
  Delta <- as.integer(!is.na(Y))
  pi1W <- if (is.null(pifit)) rep(1, length(Y)) else pifit$pi1W
  pi0W <- if (is.null(pifit)) rep(1, length(Y)) else pifit$pi0W

  mle_est <- mle_h(h_vec, Qfit$Q1W, Qfit$Q0W)
  ipw_est <- ipw_h(Y, A, h_vec, gfit$g1W, Delta = Delta, pi1W = pi1W, pi0W = pi0W)
  aipw_out <- aipw_h(Y, A, h_vec, Qfit$Q1W, Qfit$Q0W, gfit$g1W,
                     Delta = Delta, pi1W = pi1W, pi0W = pi0W)
  tmle_out <- tmle_h(Y, A, h_vec, Qfit$QAW, Qfit$Q1W, Qfit$Q0W, gfit$g1W,
                     Delta = Delta, pi1W = pi1W, pi0W = pi0W)
  boot_out <- tmle_bootstrap_h(Y, A, h_vec, Qfit$Q1W, Qfit$Q0W, gfit$g1W,
                                tmle_fit = tmle_out, n_boot = n_boot, strata = strata,
                                Delta = Delta, pi1W = pi1W, pi0W = pi0W)

  list(h = h_vec,
       MLE = mle_est,
       IPW = ipw_est,
       AIPW = aipw_out$estimate, AIPW_se = aipw_out$se,
       TMLE = tmle_out$estimate, TMLE_se = tmle_out$se, TMLE_CI = tmle_out$CI,
       Bootstrap_SE = boot_out$se,
       Bootstrap_CI_normal = boot_out$CI_normal,
       Bootstrap_CI_percentile = boot_out$CI_percentile,
       Qfit = Qfit, gfit = gfit, pifit = pifit, tmle_fit = tmle_out, boot_fit = boot_out)
}

#' Estimate a tilted ATE for a fixed h from raw data
#'
#' The single entry point: given raw data `(Y, A, W)`, an `h(x)`
#' specification (fixed vector, or `function(g1W, g0W)`), and choices for
#' how `Q`/`g` are fit, this fits `Q` and `g`, builds `h`, and returns
#' MLE/IPW/AIPW/TMLE point estimates plus TMLE's EIF-Wald and
#' target-step bootstrap inference.
#'
#' @param Y Binary (`0`/`1`) outcome vector; see [fit_Q()] Details.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param W Covariate data frame.
#' @param h Tilting-function specification; see [estimate_tau_h_fixed_given_fits()].
#' @param Q_method,Qform,Q_SL.library Passed to [fit_Q()] as `method`,
#'   `Qform`, `SL.library` (`Q_method` is `"glm"` or `"SL"`).
#' @param g_method,gform,g_SL.library Passed to [fit_g()] as `method`,
#'   `gform`, `SL.library` (`g_method` is `"glm"` or `"SL"`).
#' @param truncate If `TRUE` (default), `gbounds` is used as given
#'   (default `c(0.025, 0.975)`). If `FALSE`, `gbounds` is forced to
#'   `c(0, 1)`, i.e. no truncation.
#' @param gbounds Propensity-score truncation bounds; see `truncate`.
#' @param n_boot,boot_strata,warn_fn_h Passed to [estimate_tau_h_fixed_given_fits()].
#' @param pi_method,piform,pi_SL.library,pibounds Passed to [fit_pi()] as `method`,
#'   `piform`, `SL.library`, `pibounds`. Used only when `Y` has missing values.
#' @param Q_args,g_args,pi_args Named lists of extra arguments passed ONLY to [fit_Q()] /
#'   [fit_g()] / [fit_pi()] (e.g. `Q_args = list(cvControl = list(V = 5))` for SL).
#'   Replaces the old shared `...`, so settings meant for one fit can no longer leak
#'   into another.
#'
#' @section Missing outcomes:
#' `Y` may contain `NA`. Outcomes are assumed missing at random given `A` and `W`.
#' Let `Delta = 1` if `Y` is observed. `Q` is fit on the observed rows and predicted
#' for all rows; the missingness propensity `pi_a(W) = P(Delta = 1 | A = a, W)` is fit
#' with [fit_pi()]; and `g` is replaced by `g * pi` in the clever covariate and the
#' EIF (see the Remark on missing data in Fang, He, and Dai, 2026). `A` and `W` must
#' be complete; see [add_missing_indicators()] for missing covariates.
#'
#' @return See [estimate_tau_h_fixed_given_fits()].
#'
#' @references Fang, He, and Dai, "Targeted Learning for Externally
#'   Controlled Trials and Beyond" (AbbVie, 2026).
#' @seealso [estimate_tau_h_internal_tilted()], [estimate_tau_h_external_tilted()]
#'   for the other two tilting-function settings this package implements,
#'   [fit_Q()], [fit_g()] for the underlying nuisance fits.
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
#' out_att <- estimate_tau_h_fixed(Y, A, W, h = function(g1W, g0W) g1W, n_boot = 50,
#'                               warn_fn_h = FALSE)
#' out_att$TMLE
#' @export
estimate_tau_h_fixed <- function(Y, A, W, h,
                                  Q_method = c("glm", "SL"), Qform = NULL, Q_SL.library = c("SL.glm"),
                                  g_method = c("glm", "SL"), gform = NULL, g_SL.library = c("SL.glm"),
                                  truncate = TRUE, gbounds = c(0.025, 0.975),
                                  n_boot = 300, boot_strata = NULL, warn_fn_h = TRUE,
                                  pi_method = c("glm", "SL"), piform = NULL,
                                  pi_SL.library = c("SL.glm"), pibounds = c(0.025, 1),
                                  Q_args = list(), g_args = list(), pi_args = list()) {
  Q_method <- match.arg(Q_method)
  g_method <- match.arg(g_method)
  pi_method <- match.arg(pi_method)
  if (!truncate) gbounds <- c(0, 1)
  check_inputs(Y, A, W, gbounds, allow_missing_Y = TRUE)
  check_bounds(pibounds, "pibounds")
  W <- as.data.frame(W)

  # Q fit on observed outcomes, predicted for all rows; pi fit only if Y has NAs
  Qfit <- fit_Q_observed(Y, A, W, method = Q_method, Qform = Qform,
                         SL.library = Q_SL.library, Q_args = Q_args)
  pifit <- fit_pi_if_missing(Y, A, W, method = pi_method, piform = piform,
                             SL.library = pi_SL.library, pibounds = pibounds, pi_args = pi_args)

  # Treatment propensity g (A is observed for every row)
  gfit <- do.call(fit_g, c(list(A, W, method = g_method, gform = gform,
                                SL.library = g_SL.library, gbounds = gbounds), g_args))

  estimate_tau_h_fixed_given_fits(Y, A, h, Qfit, gfit, n_boot = n_boot, boot_strata = boot_strata,
                                  warn_fn_h = warn_fn_h, pifit = pifit)
}
