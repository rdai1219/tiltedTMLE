#' Efficient influence function for an internally tilted ATE
#'
#' Evaluated at any `(Q, g, tau)`; used both for the Step-3 convergence
#' check and for the final standard error at convergence inside
#' [tmle_tilted_h()].
#'
#' @param Y Binary (`0`/`1`) outcome vector; see [fit_Q()] Details.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param eta_g `eta(g1W)`.
#' @param eta_dot_g `eta_dot(g1W)`.
#' @param Q1W,Q0W Fitted outcome regression at `A = 1` and `A = 0`.
#' @param g1W,g0W Fitted treatment propensity at `A = 1` and `A = 0`.
#' @param tau Current estimate of the tilted ATE.
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
eif_tilted_h <- function(Y, A, eta_g, eta_dot_g, Q1W, Q0W, g1W, g0W, tau) {
  theta2 <- mean(eta_g)
  eta_g / theta2 * (A / g1W * (Y - Q1W) - (1 - A) / g0W * (Y - Q0W)) +
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
#' @param Y Binary (`0`/`1`) outcome vector; see [fit_Q()] Details.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param Q1W,Q0W,QAW Fitted outcome regression at `A = 1`, `A = 0`, and
#'   the observed treatment.
#' @param g1W,g0W Fitted treatment propensity at `A = 1` and `A = 0`.
#' @param eta,eta_dot Tilting function and its derivative, each a
#'   function of `g`; see [tilt_spec()].
#' @param max_iter Maximum number of targeting iterations.
#' @param tol Convergence tolerance on the empirical mean of the EIF.
#'
#' @return A list with components:
#' \item{estimate}{TMLE point estimate of the internally tilted ATE.}
#' \item{se}{EIF-based standard error, at convergence.}
#' \item{CI}{Wald 95% confidence interval, `estimate +/- 1.96 * se`.}
#' \item{iter}{Number of targeting iterations used (`<= max_iter`).}
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
                           max_iter = 50, tol = 1e-6) {
  N <- length(Y)
  g1 <- g1W
  g0 <- g0W
  Q1 <- Q1W
  Q0 <- Q0W
  QA <- QAW

  eta_g <- eta(g1)
  tau <- sum(eta_g * (Q1 - Q0)) / sum(eta_g)
  EIF <- eif_tilted_h(Y, A, eta_g, eta_dot(g1), Q1, Q0, g1, g0, tau)

  iter_used <- 0
  for (k in seq_len(max_iter)) {
    iter_used <- k

    eta_g <- eta(g1)
    H1_A1 <- eta_g / g1
    H1_A0 <- -eta_g / g0
    H1_obs <- A * H1_A1 + (1 - A) * H1_A0
    eps1 <- fit_eps_h(Y, H1_obs, stats::qlogis(QA))
    QA <- stats::plogis(stats::qlogis(QA) + eps1 * H1_obs)
    Q1 <- stats::plogis(stats::qlogis(Q1) + eps1 * H1_A1)
    Q0 <- stats::plogis(stats::qlogis(Q0) + eps1 * H1_A0)

    H2 <- eta_dot(g1) * (Q1 - Q0 - tau)
    eps2 <- fit_eps_h(A, H2, stats::qlogis(g1))
    g1 <- stats::plogis(stats::qlogis(g1) + eps2 * H2)
    g0 <- 1 - g1

    eta_g <- eta(g1)
    tau <- sum(eta_g * (Q1 - Q0)) / sum(eta_g)

    EIF <- eif_tilted_h(Y, A, eta_g, eta_dot(g1), Q1, Q0, g1, g0, tau)
    if (abs(mean(EIF)) < tol) break
  }

  se <- sqrt(mean(EIF^2) / N)
  ci <- tau + c(-1, 1) * stats::qnorm(0.975) * se

  list(estimate = tau, se = se, CI = ci, iter = iter_used, EIF = EIF,
       Q1W = Q1, Q0W = Q0, QAW = QA, g1W = g1, g0W = g0)
}

#' Estimate an internally tilted ATE given already-fitted Q and g
#'
#' MLE/IPW/AIPW here are the un-targeted benchmark estimators of Sec 2.1,
#' evaluated at `h = eta(initial fitted g)` -- i.e. NOT re-targeted the
#' way TMLE's `g` is. TMLE itself iteratively targets both `Q` and `g`
#' via [tmle_tilted_h()] (Sec 2.3).
#'
#' @param Y Binary (`0`/`1`) outcome vector; see [fit_Q()] Details.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param Qfit A fitted-Q object, as returned by [fit_Q()].
#' @param gfit A fitted-g object, as returned by [fit_g()].
#' @param eta,eta_dot Tilting function and its derivative; see [tilt_spec()].
#' @param max_iter,tol Passed to [tmle_tilted_h()].
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
                                                        max_iter = 50, tol = 1e-6) {
  h_vec <- eta(gfit$g1W)

  mle_est <- mle_h(h_vec, Qfit$Q1W, Qfit$Q0W)
  ipw_est <- ipw_h(Y, A, h_vec, gfit$g1W)
  aipw_out <- aipw_h(Y, A, h_vec, Qfit$Q1W, Qfit$Q0W, gfit$g1W)

  tmle_out <- tmle_tilted_h(Y, A, Qfit$Q1W, Qfit$Q0W, Qfit$QAW,
                             gfit$g1W, gfit$g0W, eta, eta_dot,
                             max_iter = max_iter, tol = tol)

  list(h = h_vec,
       MLE = mle_est,
       IPW = ipw_est,
       AIPW = aipw_out$estimate, AIPW_se = aipw_out$se,
       TMLE = tmle_out$estimate, TMLE_se = tmle_out$se, TMLE_CI = tmle_out$CI,
       TMLE_iter = tmle_out$iter,
       Qfit = Qfit, gfit = gfit, tmle_fit = tmle_out)
}

#' Estimate an internally tilted ATE from raw data
#'
#' The single entry point: given raw data `(Y, A, W)` and a tilting
#' function pair `(eta, eta_dot)` (e.g. from [tilt_spec()]), fits `Q`
#' and `g`, then returns MLE/IPW/AIPW benchmarks plus the iterative
#' TMLE and its EIF-based inference.
#'
#' @param Y Binary (`0`/`1`) outcome vector; see [fit_Q()] Details.
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
#' @param ... Passed through to [fit_Q()]/[fit_g()] (e.g. `cvControl` for SL).
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
                                            max_iter = 50, tol = 1e-6, ...) {
  Q_method <- match.arg(Q_method)
  g_method <- match.arg(g_method)
  if (!truncate) gbounds <- c(0, 1)

  Qfit <- fit_Q(Y, A, W, method = Q_method, Qform = Qform, SL.library = Q_SL.library, ...)
  gfit <- fit_g(A, W, method = g_method, gform = gform, SL.library = g_SL.library,
                gbounds = gbounds, ...)

  estimate_tau_h_internal_tilted_given_fits(Y, A, Qfit, gfit, eta, eta_dot,
                                             max_iter = max_iter, tol = tol)
}
