#' Efficient influence function for an internally tilted ATE
#'
#' Evaluated at any `(Q, g, tau)`; used both for the Step-3 convergence
#' check and for the final standard error at convergence inside
#' [tmle_tilted_h()].
#'
#' @param Y Outcome vector on the `[0, 1]` scale (binary `0`/`1`, or a
#'   continuous outcome rescaled to `[0, 1]`); see [fit_Q()] Details.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param eta_g `eta(g1W)`.
#' @param eta_dot_g `eta_dot(g1W)`.
#' @param Q1W,Q0W Fitted outcome regression at `A = 1` and `A = 0`.
#' @param g1W,g0W Fitted treatment propensity at `A = 1` and `A = 0`.
#' @param tau Current estimate of the tilted ATE.
#' @param Delta Indicator that `Y` is observed; default `1 - is.na(Y)`.
#' @param pi1W,pi0W Missing-outcome propensity `P(Delta = 1 | A = a, W)`;
#'   default `1` (no missing outcomes). See [fit_pi()].
#'
#' @details With missing outcomes (missing at random given `A`, `W`) only the
#'   residual term changes: `(A/g1W - (1-A)/g0W) * (Y - Q)` becomes
#'   `Delta * (A/(g1W*pi1W) - (1-A)/(g0W*pi0W)) * (Y - Q)`. The other two
#'   terms involve only `A` and `W`, which are fully observed.
#'
#' @return Numeric vector, the efficient influence function evaluated at
#'   each observation.
#' @seealso [tmle_tilted_h()], which calls this at every iteration.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' Qfit <- fit_Q(Y, A, W)
#' gfit <- fit_g(A, W)
#' spec <- tilt_spec("ATT")
#' eta_g <- spec$eta(gfit$g1W)
#' tau <- sum(eta_g * (Qfit$Q1W - Qfit$Q0W)) / sum(eta_g)
#' EIF <- eif_tilted_h(Y, A, eta_g, spec$eta_dot(gfit$g1W),
#'                      Qfit$Q1W, Qfit$Q0W, gfit$g1W, gfit$g0W, tau)
#' mean(EIF)
#' @export
eif_tilted_h <- function(Y, A, eta_g, eta_dot_g, Q1W, Q0W, g1W, g0W, tau,
                         Delta = as.integer(!is.na(Y)),
                         pi1W = rep(1, length(Y)), pi0W = rep(1, length(Y))) {
  theta2 <- mean(eta_g)
  Y0 <- ifelse(Delta == 1, Y, 0)   # missing Y never used: multiplied by Delta = 0
  # residual term: only observed outcomes, weighted by 1 / P(A = a and observed | W)
  eta_g / theta2 * Delta * (A / (g1W * pi1W) * (Y0 - Q1W) -
                              (1 - A) / (g0W * pi0W) * (Y0 - Q0W)) +
    eta_g / theta2 * (Q1W - Q0W - tau) +
    eta_dot_g / theta2 * (Q1W - Q0W - tau) * (A - g1W)
}

#' Iterative TMLE for an internally tilted ATE
#'
#' Implements Sec 2.3 of the reference paper: an iterative TMLE that
#' alternately targets the outcome regression `Q` and the treatment
#' propensity `g`, since the tilting function `h(x) = eta(g(x))` depends
#' on `g` itself.
#'
#' Each iteration: Step 2(a) targets `Q` via a logistic fluctuation with
#' clever covariate `eta(g)/g1W` (`A = 1`) / `-eta(g)/g0W` (`A = 0`); Step
#' 2(b) targets `g` via a logistic fluctuation with clever covariate
#' `eta_dot(g) * (Q1W - Q0W - tau)`; Step 2(c) recomputes `tau`. The loop
#' (Step 3) stops once the empirical mean of the EIF is below `tol` or
#' `max_iter` is reached. `eta_dot = 0` (the plain ATE) collapses Step
#' 2(b) to a no-op via [fit_eps_h()]'s existing `eps <- 0` guard for a
#' degenerate (all-zero-covariate) fit, so the ATE converges in a single
#' iteration without a special-cased branch.
#'
#' @param Y Outcome vector on the `[0, 1]` scale (binary `0`/`1`, or a
#'   continuous outcome rescaled to `[0, 1]`); see [fit_Q()] Details.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param Q1W,Q0W,QAW Fitted outcome regression at `A = 1`, `A = 0`, and
#'   the observed treatment.
#' @param g1W,g0W Fitted treatment propensity at `A = 1` and `A = 0`.
#' @param eta,eta_dot Tilting function and its derivative, each a
#'   function of `g`; see [tilt_spec()].
#' @param gbounds Truncation bounds re-applied to `g1` after EVERY fluctuation
#'   update (Step 2(b)), mirroring `sbounds` in [tmle_external_tilted_h()].
#'   Default `c(0, 1)` is a no-op; [estimate_tau_h_internal_tilted()] passes
#'   its own `gbounds`.
#' @param Qbounds Bounds re-applied to `Q` after every update; see [fit_Q()].
#' @param max_iter Maximum number of targeting iterations.
#' @param tol Convergence tolerance on the empirical mean of the EIF.
#' @param Delta,pi1W,pi0W Missing-outcome indicator and propensity, as in
#'   [eif_tilted_h()]. With missing outcomes the `Q` clever covariate becomes
#'   `eta(g)/(g1W*pi1W)` / `-eta(g)/(g0W*pi0W)` and `Q` is targeted on the rows
#'   with `Delta = 1`. `pi` is not targeted (it does not enter the tilt
#'   `eta(g)`), and the `g` step is unchanged because `A` is always observed.
#'
#' @return A list with components:
#' \item{estimate}{TMLE point estimate of the internally tilted ATE.}
#' \item{se}{EIF-based standard error, at convergence.}
#' \item{CI}{Wald 95% confidence interval, `estimate +/- 1.96 * se`.}
#' \item{iter}{Number of targeting iterations used (`<= max_iter`).}
#' \item{converged}{`TRUE` if the stopping rule was met before `max_iter`.}
#' \item{EIF}{Efficient influence function at convergence, evaluated at each observation.}
#' \item{Q1W,Q0W,QAW}{Targeted outcome regression at `A = 1`, `A = 0`, and the observed treatment.}
#' \item{g1W,g0W}{Targeted treatment propensity at `A = 1` and `A = 0`.}
#'
#' @seealso [eif_tilted_h()] for the EIF this iterates to convergence,
#'   [estimate_tau_h_internal_tilted()] for the raw-data entry point.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' Qfit <- fit_Q(Y, A, W)
#' gfit <- fit_g(A, W)
#' spec <- tilt_spec("ATT")
#' out <- tmle_tilted_h(Y, A, Qfit$Q1W, Qfit$Q0W, Qfit$QAW, gfit$g1W, gfit$g0W,
#'                       spec$eta, spec$eta_dot)
#' out$estimate
#' out$iter
#' @export
tmle_tilted_h <- function(Y, A, Q1W, Q0W, QAW, g1W, g0W, eta, eta_dot,
                           gbounds = c(0, 1), Qbounds = c(1e-4, 1 - 1e-4),
                           max_iter = 50, tol = 1e-6,
                           Delta = as.integer(!is.na(Y)),
                           pi1W = rep(1, length(Y)), pi0W = rep(1, length(Y))) {
  N <- length(Y)
  obs <- which(Delta == 1)   # rows with an observed outcome
  g1 <- g1W
  g0 <- g0W
  Q1 <- Q1W
  Q0 <- Q0W
  QA <- QAW

  eta_g <- eta(g1)
  tau <- sum(eta_g * (Q1 - Q0)) / sum(eta_g)
  EIF <- eif_tilted_h(Y, A, eta_g, eta_dot(g1), Q1, Q0, g1, g0, tau, Delta, pi1W, pi0W)

  iter_used <- 0
  converged <- FALSE
  for (k in seq_len(max_iter)) {
    iter_used <- k

    eta_g <- eta(g1)
    H1_A1 <- eta_g / (g1 * pi1W)     # clever covariate: g replaced by g * pi
    H1_A0 <- -eta_g / (g0 * pi0W)
    H1_obs <- A * H1_A1 + (1 - A) * H1_A0
    eps1 <- fit_eps_h(Y, H1_obs, stats::qlogis(QA), indices = obs)  # observed rows only
    QA <- stats::plogis(stats::qlogis(QA) + eps1 * H1_obs)
    Q1 <- stats::plogis(stats::qlogis(Q1) + eps1 * H1_A1)
    Q0 <- stats::plogis(stats::qlogis(Q0) + eps1 * H1_A0)
    # Fix 1: keep Q finite on the logit scale
    QA <- clamp_bounds(QA, Qbounds)
    Q1 <- clamp_bounds(Q1, Qbounds)
    Q0 <- clamp_bounds(Q0, Qbounds)

    H2 <- eta_dot(g1) * (Q1 - Q0 - tau)
    eps2 <- fit_eps_h(A, H2, stats::qlogis(g1))
    g1 <- stats::plogis(stats::qlogis(g1) + eps2 * H2)
    g1 <- clamp_bounds(g1, gbounds)  # Fix 2: re-clamp g after fluctuation
    g0 <- 1 - g1

    eta_g <- eta(g1)
    tau <- sum(eta_g * (Q1 - Q0)) / sum(eta_g)

    EIF <- eif_tilted_h(Y, A, eta_g, eta_dot(g1), Q1, Q0, g1, g0, tau, Delta, pi1W, pi0W)
    # Fix 3: same stopping rule as the external TMLE (EIF mean ~ 0, or no more movement)
    if (abs(mean(EIF)) < tol || (abs(eps1) < tol && abs(eps2) < tol)) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    warning("tmle_tilted_h() did not converge in ", max_iter,
            " iterations; consider a larger max_iter or tol.", call. = FALSE)
  }

  se <- sqrt(mean(EIF^2) / N)
  ci <- tau + c(-1, 1) * stats::qnorm(0.975) * se

  list(estimate = tau, se = se, CI = ci, iter = iter_used, converged = converged, EIF = EIF,
       Q1W = Q1, Q0W = Q0, QAW = QA, g1W = g1, g0W = g0)
}

#' Estimate an internally tilted ATE given already-fitted Q and g
#'
#' MLE/IPW/AIPW here are the un-targeted benchmark estimators of Sec 2.1,
#' evaluated at `h = eta(initial fitted g)` -- i.e. NOT re-targeted the
#' way TMLE's `g` is. TMLE itself iteratively targets both `Q` and `g`
#' via [tmle_tilted_h()] (Sec 2.3).
#'
#' @param Y Outcome vector on the `[0, 1]` scale (binary `0`/`1`, or a
#'   continuous outcome rescaled to `[0, 1]`); see [fit_Q()] Details.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param Qfit A fitted-Q object, as returned by [fit_Q()].
#' @param gfit A fitted-g object, as returned by [fit_g()].
#' @param eta,eta_dot Tilting function and its derivative; see [tilt_spec()].
#' @param gbounds Re-clamp bounds for `g` inside the loop; see [tmle_tilted_h()].
#' @param max_iter,tol Passed to [tmle_tilted_h()].
#' @param pifit Optional fitted missing-outcome propensity, as returned by
#'   [fit_pi()]. `NULL` (default) means no missing outcomes.
#'
#' @return A list with components:
#' \item{h}{The resolved tilting weight vector, `eta(gfit$g1W)`.}
#' \item{MLE}{Plug-in (g-computation) benchmark estimate; see [mle_h()].}
#' \item{IPW}{Inverse-probability-weighted benchmark estimate; see [ipw_h()].}
#' \item{AIPW,AIPW_se}{Augmented IPW benchmark estimate and its EIF-based
#'   standard error; see [aipw_h()].}
#' \item{TMLE,TMLE_se,TMLE_CI}{Iterative TMLE point estimate, EIF-based
#'   standard error, and Wald 95% confidence interval; see [tmle_tilted_h()].}
#' \item{TMLE_iter}{Number of targeting iterations the TMLE used.}
#' \item{TMLE_converged}{Whether the TMLE met its stopping rule.}
#' \item{Qfit,gfit}{The fitted-Q and fitted-g objects passed in.}
#' \item{tmle_fit}{The full [tmle_tilted_h()] result.}
#'
#' @seealso [estimate_tau_h_internal_tilted()], the raw-data entry point
#'   that calls this after fitting `Q`/`g`.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' Qfit <- fit_Q(Y, A, W)
#' gfit <- fit_g(A, W)
#' spec <- tilt_spec("ATT")
#' out <- estimate_tau_h_internal_tilted_given_fits(Y, A, Qfit, gfit, spec$eta, spec$eta_dot)
#' out$TMLE
#' @export
estimate_tau_h_internal_tilted_given_fits <- function(Y, A, Qfit, gfit, eta, eta_dot,
                                                        gbounds = c(0, 1),
                                                        max_iter = 50, tol = 1e-6,
                                                        pifit = NULL) {
  h_vec <- eta(gfit$g1W)

  # Missing outcomes: Delta = 1 where Y is observed; pi = P(Delta = 1 | A, W).
  # Without a pifit (no missing outcomes) pi = 1 and nothing changes.
  Delta <- as.integer(!is.na(Y))
  pi1W <- if (is.null(pifit)) rep(1, length(Y)) else pifit$pi1W
  pi0W <- if (is.null(pifit)) rep(1, length(Y)) else pifit$pi0W

  mle_est <- mle_h(h_vec, Qfit$Q1W, Qfit$Q0W)
  ipw_est <- ipw_h(Y, A, h_vec, gfit$g1W, Delta = Delta, pi1W = pi1W, pi0W = pi0W)
  aipw_out <- aipw_h(Y, A, h_vec, Qfit$Q1W, Qfit$Q0W, gfit$g1W,
                     Delta = Delta, pi1W = pi1W, pi0W = pi0W)

  tmle_out <- tmle_tilted_h(Y, A, Qfit$Q1W, Qfit$Q0W, Qfit$QAW,
                             gfit$g1W, gfit$g0W, eta, eta_dot,
                             gbounds = gbounds, max_iter = max_iter, tol = tol,
                             Delta = Delta, pi1W = pi1W, pi0W = pi0W)

  list(h = h_vec,
       MLE = mle_est,
       IPW = ipw_est,
       AIPW = aipw_out$estimate, AIPW_se = aipw_out$se,
       TMLE = tmle_out$estimate, TMLE_se = tmle_out$se, TMLE_CI = tmle_out$CI,
       TMLE_iter = tmle_out$iter, TMLE_converged = tmle_out$converged,
       Qfit = Qfit, gfit = gfit, pifit = pifit, tmle_fit = tmle_out)
}

#' Estimate an internally tilted ATE from raw data
#'
#' The single entry point: given raw data `(Y, A, W)` and a tilting
#' function pair `(eta, eta_dot)` (e.g. from [tilt_spec()]), fits `Q`
#' and `g`, then returns MLE/IPW/AIPW benchmarks plus the iterative
#' TMLE and its EIF-based inference.
#'
#' @param Y Outcome vector: binary (`0`/`1`) or continuous, with `NA` for a
#'   missing outcome. A continuous outcome is rescaled internally; see `Ybounds`.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param W Covariate data frame.
#' @param eta,eta_dot Tilting function and its derivative; see [tilt_spec()].
#' @param Q_method,Qform,Q_SL.library Passed to [fit_Q()] as `method`,
#'   `Qform`, `SL.library` (`Q_method` is `"glm"` or `"SL"`).
#' @param g_method,gform,g_SL.library Passed to [fit_g()] as `method`,
#'   `gform`, `SL.library` (`g_method` is `"glm"` or `"SL"`).
#' @param truncate If `TRUE` (default), `gbounds` is used as given
#'   (default `c(0.025, 0.975)`). If `FALSE`, `gbounds` is forced to
#'   `c(0, 1)`, i.e. no truncation.
#' @param gbounds Propensity-score truncation bounds; see `truncate`.
#' @param max_iter,tol Passed to [tmle_tilted_h()].
#' @param pi_method,piform,pi_SL.library,pibounds Passed to [fit_pi()]; used only
#'   when `Y` has missing values.
#' @param Ybounds Optional bounds `c(a, b)` for a continuous outcome. `Y` is
#'   rescaled to `(Y - a) / (b - a)` before fitting and the results are mapped
#'   back. `NULL` (default) uses `c(0, 1)` for a binary outcome and the observed
#'   range of `Y` otherwise. Known limits (e.g. the range of a rating scale)
#'   can be supplied instead; they must contain every observed `Y`.
#' @param Q_args,g_args,pi_args Named lists of extra arguments passed ONLY to [fit_Q()] /
#'   [fit_g()] / [fit_pi()] (e.g. `Q_args = list(cvControl = list(V = 5))` for SL).
#'   Replaces the old shared `...`, so settings meant for one fit can no longer leak
#'   into another.
#'
#' @section Continuous outcomes:
#' `Y` may be continuous. Following Fang, He, and Dai (2026, Sec. 3.3.2), `Y` is
#' rescaled to `[0, 1]` using `Ybounds`, the TMLE is computed on that scale
#' (quasi-binomial outcome model and logistic fluctuation), and the point
#' estimates, standard errors and confidence limits are multiplied by `b - a`,
#' so they are reported on the original outcome scale. The returned list also
#' contains `Ybounds`; the fitted nuisance objects (`Qfit`, `tmle_fit`, ...)
#' stay on the rescaled `[0, 1]` scale. A binary outcome is not rescaled.
#'
#' @section Missing outcomes:
#' `Y` may contain `NA` (outcome missing at random given `A` and `W`). `Q` is fit on
#' the observed rows and predicted for all rows, `pi_a(W) = P(Delta = 1 | A = a, W)`
#' is fit with [fit_pi()], and `g` is replaced by `g * pi` in the `Q` clever covariate
#' and the EIF; see [tmle_tilted_h()]. `A` and `W` must be complete.
#'
#' @return See [estimate_tau_h_internal_tilted_given_fits()].
#'
#' @references Fang, "Targeted Learning for Tilted Average Treatment
#'   Effects" (AbbVie, 2026), Sec. 2.
#' @seealso [tilt_spec()] for the four named internal estimands,
#'   [estimate_tau_h_fixed()], [estimate_tau_h_external_tilted()] for
#'   the other two tilting-function settings this package implements.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#'
#' ## ATT: tilt the ATE by h(x) = eta(g(x)) = g(x)
#' spec <- tilt_spec("ATT")
#' out <- estimate_tau_h_internal_tilted(Y, A, W, spec$eta, spec$eta_dot)
#' out$TMLE
#' out$TMLE_iter
#' @export
estimate_tau_h_internal_tilted <- function(Y, A, W, eta, eta_dot,
                                            Q_method = c("glm", "SL"), Qform = NULL, Q_SL.library = c("SL.glm"),
                                            g_method = c("glm", "SL"), gform = NULL, g_SL.library = c("SL.glm"),
                                            truncate = TRUE, gbounds = c(0.025, 0.975),
                                            max_iter = 50, tol = 1e-6,
                                            pi_method = c("glm", "SL"), piform = NULL,
                                            pi_SL.library = c("SL.glm"), pibounds = c(0.025, 1),
                                            Ybounds = NULL,
                                            Q_args = list(), g_args = list(), pi_args = list()) {
  Q_method <- match.arg(Q_method)
  g_method <- match.arg(g_method)
  pi_method <- match.arg(pi_method)
  if (!truncate) gbounds <- c(0, 1)
  check_inputs(Y, A, W, gbounds, allow_missing_Y = TRUE)
  check_bounds(pibounds, "pibounds")
  W <- as.data.frame(W)

  # Continuous outcome: work with Y rescaled to [0, 1]; results are mapped back
  # at the end (a binary outcome has Ybounds = c(0, 1) and is unchanged).
  Ybounds <- outcome_bounds(Y, Ybounds)
  Y <- rescale_outcome(Y, Ybounds)

  # Q fit on observed outcomes, predicted for all rows; pi fit only if Y has NAs
  Qfit <- fit_Q_observed(Y, A, W, method = Q_method, Qform = Qform,
                         SL.library = Q_SL.library, Q_args = Q_args)
  pifit <- fit_pi_if_missing(Y, A, W, method = pi_method, piform = piform,
                             SL.library = pi_SL.library, pibounds = pibounds, pi_args = pi_args)
  gfit <- do.call(fit_g, c(list(A, W, method = g_method, gform = gform,
                                SL.library = g_SL.library, gbounds = gbounds), g_args))

  out <- estimate_tau_h_internal_tilted_given_fits(Y, A, Qfit, gfit, eta, eta_dot,
                                                   gbounds = gbounds,
                                                   max_iter = max_iter, tol = tol,
                                                   pifit = pifit)
  unscale_results(out, Ybounds)
}
