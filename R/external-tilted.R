#' Efficient influence function for an externally tilted ATE
#'
#' Evaluated at any `(Q, g, s, tau)`; used both for the Step-3
#' convergence check and for the final standard error at convergence
#' inside [tmle_external_tilted_h()].
#'
#' @param D Length-`N` pooled selection indicator (`1` = internal, `0` = external).
#' @param A_int,Y_int Length-`n` internal-only treatment and binary
#'   (`0`/`1`) outcome vectors; see [fit_Q()] Details.
#' @param idx1 Integer vector, `which(D == 1)`.
#' @param eta_s,eta_dot_s Length-`N` pooled vectors, `eta(s1W)` and `eta_dot(s1W)`.
#' @param Q1W,Q0W Length-`N` pooled fitted outcome regression at `A = 1`
#'   and `A = 0` (extrapolated onto external rows).
#' @param g1W_int,g0W_int Length-`n` internal-only fitted treatment
#'   propensity at `A = 1` and `A = 0`.
#' @param s1W Length-`N` pooled fitted selection propensity.
#' @param tau Current estimate of the tilted ATE.
#'
#' @return Numeric vector of length `N`, the efficient influence function
#'   evaluated at each pooled observation.
#' @seealso [tmle_external_tilted_h()], which calls this at every iteration.
#' @examples
#' set.seed(1)
#' n <- 150; m <- 150
#' Ws <- data.frame(X1 = rnorm(n + m))
#' D <- c(rep(1L, n), rep(0L, m))
#' A_int <- rbinom(n, 1, stats::plogis(0.3 * Ws$X1[1:n]))
#' Y_int <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A_int + 0.3 * Ws$X1[1:n]))
#' idx1 <- which(D == 1)
#'
#' Qfit_int <- fit_Q(Y_int, A_int, Ws[idx1, , drop = FALSE])
#' Q1W <- predict_Qhat(Qfit_int$fit, 1, Ws)
#' Q0W <- predict_Qhat(Qfit_int$fit, 0, Ws)
#' gfit <- fit_g(A_int, Ws[idx1, , drop = FALSE])
#' sfit <- fit_g(D, Ws)
#'
#' spec <- tilt_spec_external("ATE-O")
#' eta_s <- spec$eta(sfit$g1W)
#' tau <- sum(eta_s[idx1] * (Q1W[idx1] - Q0W[idx1])) / sum(eta_s[idx1])
#' EIF <- eif_external_tilted_h(D, A_int, Y_int, idx1, eta_s, spec$eta_dot(sfit$g1W),
#'                               Q1W, Q0W, gfit$g1W, gfit$g0W, sfit$g1W, tau)
#' mean(EIF)
#' @export
eif_external_tilted_h <- function(D, A_int, Y_int, idx1, eta_s, eta_dot_s,
                                   Q1W, Q0W, g1W_int, g0W_int, s1W, tau) {
  N <- length(D)
  theta2 <- mean(D * eta_s)

  term1 <- numeric(N)
  term1[idx1] <- eta_s[idx1] * (A_int / g1W_int * (Y_int - Q1W[idx1]) -
                                 (1 - A_int) / g0W_int * (Y_int - Q0W[idx1]))
  term1 <- term1 / theta2

  term2 <- (s1W * eta_s / theta2) * (Q1W - Q0W - tau)
  term3 <- ((eta_s + s1W * eta_dot_s) / theta2) * (Q1W - Q0W - tau) * (D - s1W)

  term1 + term2 + term3
}

#' Iterative TMLE for an externally tilted ATE
#'
#' Implements Sec 3.3 of the reference paper: an iterative TMLE that
#' alternately targets the outcome regression `Q` (fit on internal rows,
#' applied to all `N` pooled rows) and the selection propensity `s` (fit
#' on all `N` pooled rows), since the tilting function `h(x) = eta(s(x))`
#' depends on `s` itself. The treatment propensity `g` is fit once,
#' upstream, and is never re-targeted inside this loop.
#'
#' Each iteration: Step 2(a) targets `Q` via a logistic fluctuation fit on
#' internal rows only but applied to all `N` rows; Step 2(b) targets `s`
#' via a logistic fluctuation fit on all `N` rows; Step 2(c) recomputes
#' `tau`. The loop stops once the empirical mean of the EIF is below
#' `tol` or `max_iter` is reached. `eta_dot = 0` (ATE-I) collapses Step
#' 2(b) to a no-op via [fit_eps_h()]'s existing `eps <- 0` guard for a
#' degenerate (all-zero-covariate) fit.
#'
#' `g1W_int`/`g0W_int` (length `n`, internal rows only) are used for the
#' EIF's term 1 and are held fixed throughout -- `g` is fit once,
#' upstream, and never re-targeted. `g1W_all`/`g0W_all` (length `N`, `g`
#' extrapolated onto ALL pooled covariates, including external rows)
#' are needed to evaluate the `Q` clever covariate at every row, since
#' `Q(Xs,1|D=1)`/`Q(Xs,0|D=1)` must be extrapolated to external
#' covariates for Step 2(b)'s clever covariate and the EIF.
#'
#' `sbounds` re-clamps `s1` after EVERY fluctuation update (Step 2(b)),
#' not just at the initial fit, because the logistic fluctuation itself
#' is unbounded and, for `eta(s)` with a `1/s` term (e.g. ATE-E), can
#' walk `s1` back below the lower truncation bound within a few
#' iterations even when the initial fit respected it, eventually
#' overflowing `Q` to exactly `0`/`1` on the logit scale and crashing
#' the next iteration's logistic regression. `sbounds = c(0, 1)` is a
#' no-op clamp, matching [fit_g()]'s `gbounds` convention.
#'
#' @param D Length-`N` pooled selection indicator.
#' @param A_int,Y_int Length-`n` internal-only treatment and binary
#'   (`0`/`1`) outcome vectors; see [fit_Q()] Details.
#' @param idx1 Integer vector, `which(D == 1)`.
#' @param Q1W,Q0W Length-`N` pooled fitted outcome regression at `A = 1`
#'   and `A = 0`.
#' @param QAW_int Length-`n` internal-only fitted outcome regression at
#'   the observed treatment.
#' @param g1W_int,g0W_int Length-`n` internal-only fitted treatment
#'   propensity at `A = 1` and `A = 0`.
#' @param g1W_all,g0W_all Length-`N` pooled fitted treatment propensity
#'   (extrapolated onto external rows), at `A = 1` and `A = 0`.
#' @param s1W Length-`N` pooled fitted selection propensity.
#' @param eta,eta_dot Tilting function and its derivative, each a
#'   function of `s`; see [tilt_spec_external()].
#' @param sbounds Truncation bounds re-applied to `s1` after every
#'   fluctuation update; see Details.
#' @param max_iter Maximum number of targeting iterations.
#' @param tol Convergence tolerance.
#'
#' @return A list with components:
#' \item{estimate}{TMLE point estimate of the externally tilted ATE.}
#' \item{se}{EIF-based standard error, at convergence.}
#' \item{CI}{Wald 95% confidence interval, `estimate +/- 1.96 * se`.}
#' \item{iter}{Number of targeting iterations used (`<= max_iter`).}
#' \item{EIF}{Efficient influence function at convergence, evaluated at
#'   each pooled observation.}
#' \item{Q1W,Q0W}{Targeted pooled outcome regression at `A = 1` and `A = 0`.}
#' \item{QAW_int}{Targeted internal-only outcome regression at the observed treatment.}
#' \item{s1W}{Targeted pooled selection propensity.}
#'
#' @seealso [eif_external_tilted_h()] for the EIF this iterates to
#'   convergence, [estimate_tau_h_external_tilted()] for the raw-data
#'   entry point.
#' @examples
#' set.seed(1)
#' n <- 150; m <- 150
#' Ws <- data.frame(X1 = rnorm(n + m))
#' D <- c(rep(1L, n), rep(0L, m))
#' A_int <- rbinom(n, 1, stats::plogis(0.3 * Ws$X1[1:n]))
#' Y_int <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A_int + 0.3 * Ws$X1[1:n]))
#' idx1 <- which(D == 1)
#'
#' Qfit_int <- fit_Q(Y_int, A_int, Ws[idx1, , drop = FALSE])
#' Q1W <- predict_Qhat(Qfit_int$fit, 1, Ws)
#' Q0W <- predict_Qhat(Qfit_int$fit, 0, Ws)
#' gfit <- fit_g(A_int, Ws[idx1, , drop = FALSE])
#' g1W_all <- predict_ghat(gfit$fit, Ws)
#' sfit <- fit_g(D, Ws)
#'
#' spec <- tilt_spec_external("ATE-O")
#' out <- tmle_external_tilted_h(D, A_int, Y_int, idx1, Q1W, Q0W, Qfit_int$QAW,
#'                                gfit$g1W, gfit$g0W, g1W_all, 1 - g1W_all,
#'                                sfit$g1W, spec$eta, spec$eta_dot)
#' out$estimate
#' out$iter
#' @export
tmle_external_tilted_h <- function(D, A_int, Y_int, idx1,
                                    Q1W, Q0W, QAW_int,
                                    g1W_int, g0W_int = 1 - g1W_int,
                                    g1W_all, g0W_all = 1 - g1W_all,
                                    s1W, eta, eta_dot,
                                    sbounds = c(0.025, 0.975),
                                    max_iter = 50, tol = 1e-6) {
  N <- length(D)
  Q1 <- Q1W
  Q0 <- Q0W
  QA_int <- QAW_int
  s1 <- s1W

  eta_s <- eta(s1)
  tau <- sum(eta_s[idx1] * (Q1[idx1] - Q0[idx1])) / sum(eta_s[idx1])
  EIF <- eif_external_tilted_h(D, A_int, Y_int, idx1, eta_s, eta_dot(s1),
                                Q1, Q0, g1W_int, g0W_int, s1, tau)

  iter_used <- 0
  for (k in seq_len(max_iter)) {
    iter_used <- k

    eta_s <- eta(s1)
    H1_A1 <- eta_s / g1W_all
    H1_A0 <- -eta_s / g0W_all
    H1_obs_int <- A_int * H1_A1[idx1] + (1 - A_int) * H1_A0[idx1]
    eps1 <- fit_eps_h(Y_int, H1_obs_int, stats::qlogis(QA_int))
    QA_int <- stats::plogis(stats::qlogis(QA_int) + eps1 * H1_obs_int)
    Q1 <- stats::plogis(stats::qlogis(Q1) + eps1 * H1_A1)
    Q0 <- stats::plogis(stats::qlogis(Q0) + eps1 * H1_A0)

    H2 <- (eta_s + s1 * eta_dot(s1)) * (Q1 - Q0 - tau)
    eps2 <- fit_eps_h(D, H2, stats::qlogis(s1))
    s1 <- stats::plogis(stats::qlogis(s1) + eps2 * H2)
    s1 <- pmin(pmax(s1, sbounds[1]), sbounds[2])

    eta_s <- eta(s1)
    tau <- sum(eta_s[idx1] * (Q1[idx1] - Q0[idx1])) / sum(eta_s[idx1])

    EIF <- eif_external_tilted_h(D, A_int, Y_int, idx1, eta_s, eta_dot(s1),
                                  Q1, Q0, g1W_int, g0W_int, s1, tau)
    if (abs(mean(EIF)) < tol) break

    if (abs(eps1) < tol && abs(eps2) < tol) break
  }

  se <- sqrt(mean(EIF^2) / N)
  ci <- tau + c(-1, 1) * stats::qnorm(0.975) * se

  list(estimate = tau, se = se, CI = ci, iter = iter_used, EIF = EIF,
       Q1W = Q1, Q0W = Q0, QAW_int = QA_int, s1W = s1)
}

#' Estimate an externally tilted ATE given already-fitted Q, g, and s
#'
#' MLE/IPW here are the un-targeted benchmark estimators of Sec 3.1,
#' evaluated at `h = eta(initial fitted s)` -- i.e. NOT re-targeted the
#' way TMLE's `s` is -- and computed on the internal (`D = 1`) rows
#' only, since `Y`/`A` are only observed there. No AIPW benchmark is
#' defined in Sec 3, so none is included here.
#'
#' Length contract: `Y`, `A` are length `n` (internal rows only); `D`,
#' `sfit$g1W`, `g1W_all` are length `N` (pooled); `gfit$g1W` is length
#' `n` (fit on internal rows only); `Qfit$Q1W`/`Q0W` must be length `N`
#' -- i.e. the outcome regression fit on internal rows but with
#' predictions EXTRAPOLATED onto all `N` pooled covariates (see
#' [predict_Qhat()] and [estimate_tau_h_external_tilted()]).
#'
#' @param Y,A Length-`n` internal-only binary (`0`/`1`) outcome and
#'   treatment vectors; see [fit_Q()] Details.
#' @param D Length-`N` pooled selection indicator.
#' @param idx1 Integer vector, `which(D == 1)`.
#' @param Qfit A fitted-Q object (length-`N` `Q1W`/`Q0W`, extrapolated
#'   onto pooled covariates; length-`n` `QAW`).
#' @param gfit A fitted-g object (length-`n`, internal rows only).
#' @param sfit A fitted-s object (length-`N`, pooled), as returned by
#'   [fit_g()] called with `D` in place of the treatment indicator.
#' @param eta,eta_dot Tilting function and its derivative; see [tilt_spec_external()].
#' @param g1W_all,g0W_all Length-`N` pooled fitted treatment propensity.
#' @param sbounds,max_iter,tol Passed to [tmle_external_tilted_h()].
#'
#' @return A list with components:
#' \item{h}{The resolved tilting weight vector (internal rows only), `eta(sfit$g1W[idx1])`.}
#' \item{MLE}{Plug-in (g-computation) benchmark estimate; see [mle_h()].}
#' \item{IPW}{Inverse-probability-weighted benchmark estimate; see [ipw_h()].
#'   (No AIPW benchmark is defined in Sec 3 of the reference paper, so none is included.)}
#' \item{TMLE,TMLE_se,TMLE_CI}{Iterative TMLE point estimate, EIF-based
#'   standard error, and Wald 95% confidence interval; see [tmle_external_tilted_h()].}
#' \item{TMLE_iter}{Number of targeting iterations the TMLE used.}
#' \item{Qfit,gfit,sfit}{The fitted-Q, fitted-g, and fitted-s objects passed in.}
#' \item{tmle_fit}{The full [tmle_external_tilted_h()] result.}
#'
#' @seealso [estimate_tau_h_external_tilted()], the raw-data entry point
#'   that calls this after fitting `Q`/`g`/`s`.
#' @examples
#' set.seed(1)
#' n <- 150; m <- 150
#' Ws <- data.frame(X1 = rnorm(n + m))
#' D <- c(rep(1L, n), rep(0L, m))
#' A_int <- rbinom(n, 1, stats::plogis(0.3 * Ws$X1[1:n]))
#' Y_int <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A_int + 0.3 * Ws$X1[1:n]))
#' idx1 <- which(D == 1)
#'
#' Qfit_int <- fit_Q(Y_int, A_int, Ws[idx1, , drop = FALSE])
#' Qfit <- list(QAW = Qfit_int$QAW,
#'              Q1W = predict_Qhat(Qfit_int$fit, 1, Ws),
#'              Q0W = predict_Qhat(Qfit_int$fit, 0, Ws),
#'              fit = Qfit_int$fit, method = Qfit_int$method)
#' gfit <- fit_g(A_int, Ws[idx1, , drop = FALSE])
#' g1W_all <- predict_ghat(gfit$fit, Ws)
#' sfit <- fit_g(D, Ws)
#'
#' spec <- tilt_spec_external("ATE-O")
#' out <- estimate_tau_h_external_tilted_given_fits(Y_int, A_int, D, idx1, Qfit, gfit, sfit,
#'                                                    spec$eta, spec$eta_dot, g1W_all = g1W_all)
#' out$TMLE
#' @export
estimate_tau_h_external_tilted_given_fits <- function(Y, A, D, idx1,
                                                        Qfit, gfit, sfit,
                                                        eta, eta_dot,
                                                        g1W_all, g0W_all = 1 - g1W_all,
                                                        sbounds = c(0.025, 0.975),
                                                        max_iter = 50, tol = 1e-6) {
  h_vec_int <- eta(sfit$g1W[idx1])

  mle_est <- mle_h(h_vec_int, Qfit$Q1W[idx1], Qfit$Q0W[idx1])
  ipw_est <- ipw_h(Y, A, h_vec_int, gfit$g1W)

  tmle_out <- tmle_external_tilted_h(D, A, Y, idx1,
                                      Qfit$Q1W, Qfit$Q0W, Qfit$QAW,
                                      gfit$g1W, gfit$g0W,
                                      g1W_all, g0W_all,
                                      sfit$g1W, eta, eta_dot,
                                      sbounds = sbounds,
                                      max_iter = max_iter, tol = tol)

  list(h = h_vec_int,
       MLE = mle_est,
       IPW = ipw_est,
       TMLE = tmle_out$estimate, TMLE_se = tmle_out$se, TMLE_CI = tmle_out$CI,
       TMLE_iter = tmle_out$iter,
       Qfit = Qfit, gfit = gfit, sfit = sfit, tmle_fit = tmle_out)
}

#' Estimate an externally tilted ATE from raw data
#'
#' The single entry point: given raw pooled data (`Y`, `A` observed
#' only for `D = 1`; `D`; `Ws` covariates for all `N = n + m` rows) and
#' a tilting function pair `(eta, eta_dot)` (e.g. from
#' [tilt_spec_external()]), fits `Q` and `g` on the internal (`D = 1`)
#' rows and `s` on all `N` rows, then returns MLE/IPW benchmarks plus
#' the iterative TMLE and its EIF-based inference.
#'
#' `Y`, `A` are length-`n` vectors (internal rows only); `D` is
#' length-`N`; `Ws` is an `N`-row covariate data frame used for both
#' `s(x)` and, via `Wg`, `g(x)` (`Wg` defaults to the same covariates,
#' restricted to the internal rows, but may be supplied separately per
#' the reference paper's convention that the confounder set for `g`,
#' `X^g`, may be a strict subset of the covariate set for `s` and `Q`,
#' `X^s`).
#'
#' @param Y,A Length-`n` internal-only binary (`0`/`1`) outcome and
#'   treatment vectors; see [fit_Q()] Details.
#' @param D Length-`N` pooled selection indicator.
#' @param Ws Length-`N` pooled covariate data frame.
#' @param eta,eta_dot Tilting function and its derivative; see [tilt_spec_external()].
#' @param Wg Optional covariate data frame for fitting `g`, restricted
#'   to internal rows (defaults to `Ws` restricted to the internal rows).
#' @param Q_method,Qform,Q_SL.library Passed to [fit_Q()] as `method`,
#'   `Qform`, `SL.library` (`Q_method` is `"glm"` or `"SL"`).
#' @param g_method,gform,g_SL.library Passed to [fit_g()] for `g`, as
#'   `method`, `gform`, `SL.library` (`g_method` is `"glm"` or `"SL"`).
#' @param s_method,sform,s_SL.library Passed to [fit_g()] for `s`, as
#'   `method`, `gform`, `SL.library` (`s_method` is `"glm"` or `"SL"`).
#' @param truncate_g If `TRUE` (default), `gbounds` is used as given.
#'   If `FALSE`, `gbounds` is forced to `c(0, 1)`.
#' @param gbounds Treatment-propensity truncation bounds.
#' @param truncate_s If `TRUE` (default), `sbounds` is used as given.
#'   If `FALSE`, `sbounds` is forced to `c(0, 1)`.
#' @param sbounds Selection-propensity truncation bounds.
#' @param max_iter,tol Passed to [tmle_external_tilted_h()].
#' @param ... Passed through to [fit_Q()]/[fit_g()] (e.g. `cvControl` for SL).
#'
#' @return See [estimate_tau_h_external_tilted_given_fits()].
#'
#' @references Fang, "Targeted Learning for Tilted Average Treatment
#'   Effects" (AbbVie, 2026), Sec. 3.
#' @seealso [tilt_spec_external()] for the three named external estimands,
#'   [estimate_tau_h_fixed()], [estimate_tau_h_internal_tilted()] for
#'   the other two tilting-function settings this package implements.
#' @examples
#' set.seed(1)
#' n <- 150; m <- 150
#' Ws <- data.frame(X1 = rnorm(n + m))
#' D <- c(rep(1L, n), rep(0L, m))
#' A <- rbinom(n, 1, stats::plogis(0.3 * Ws$X1[1:n]))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * Ws$X1[1:n]))
#'
#' ## ATE-O: tilt by h(x) = eta(s(x)) = 1 - s(x), the selection-overlap population
#' spec <- tilt_spec_external("ATE-O")
#' out <- estimate_tau_h_external_tilted(Y, A, D, Ws, spec$eta, spec$eta_dot)
#' out$TMLE
#' out$TMLE_iter
#' @export
estimate_tau_h_external_tilted <- function(Y, A, D, Ws, eta, eta_dot,
                                            Wg = NULL,
                                            Q_method = c("glm", "SL"), Qform = NULL, Q_SL.library = c("SL.glm"),
                                            g_method = c("glm", "SL"), gform = NULL, g_SL.library = c("SL.glm"),
                                            s_method = c("glm", "SL"), sform = NULL, s_SL.library = c("SL.glm"),
                                            truncate_g = TRUE, gbounds = c(0.025, 0.975),
                                            truncate_s = TRUE, sbounds = c(0.025, 0.975),
                                            max_iter = 50, tol = 1e-6, ...) {
  Q_method <- match.arg(Q_method)
  g_method <- match.arg(g_method)
  s_method <- match.arg(s_method)
  if (!truncate_g) gbounds <- c(0, 1)
  if (!truncate_s) sbounds <- c(0, 1)

  Ws <- as.data.frame(Ws)
  idx1 <- which(D == 1)
  if (is.null(Wg)) Wg <- Ws[idx1, , drop = FALSE]

  Qfit_int <- fit_Q(Y, A, Ws[idx1, , drop = FALSE], method = Q_method, Qform = Qform,
                     SL.library = Q_SL.library, ...)
  gfit <- fit_g(A, Wg, method = g_method, gform = gform, SL.library = g_SL.library,
                gbounds = gbounds, ...)
  sfit <- fit_g(D, Ws, method = s_method, gform = sform, SL.library = s_SL.library,
                gbounds = sbounds, ...)

  Qfit <- list(QAW = Qfit_int$QAW,
               Q1W = predict_Qhat(Qfit_int$fit, 1, Ws),
               Q0W = predict_Qhat(Qfit_int$fit, 0, Ws),
               fit = Qfit_int$fit, method = Qfit_int$method)
  Wg_all <- Ws[, names(Wg), drop = FALSE]
  g1W_all <- pmin(pmax(predict_ghat(gfit$fit, Wg_all), gbounds[1]), gbounds[2])

  estimate_tau_h_external_tilted_given_fits(Y, A, D, idx1, Qfit, gfit, sfit,
                                             eta, eta_dot,
                                             g1W_all = g1W_all,
                                             sbounds = sbounds,
                                             max_iter = max_iter, tol = tol)
}
