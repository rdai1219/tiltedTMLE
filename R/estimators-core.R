#' Plug-in (g-computation) estimator of a tilted ATE
#'
#' \deqn{\hat\tau_h^{MLE} = \sum_i h(X_i)[\hat Q(X_i,1)-\hat Q(X_i,0)] / \sum_i h(X_i)}
#'
#' @param h Numeric tilting-function weight vector, one entry per observation.
#' @param Q1W,Q0W Fitted outcome regression at `A = 1` and `A = 0`.
#'
#' @return Numeric point estimate of \eqn{\tau_h}.
#' @seealso [ipw_h()], [aipw_h()] for the other Sec 2.1/3.1 benchmark
#'   estimators, [tmle_h()] for the targeted (TMLE) estimator.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' Qfit <- fit_Q(Y, A, W)
#' mle_h(h = rep(1, n), Qfit$Q1W, Qfit$Q0W)
#' @export
mle_h <- function(h, Q1W, Q0W) {
  sum(h * (Q1W - Q0W)) / sum(h)
}

#' Inverse-probability-weighted estimator of a tilted ATE
#'
#' @param Y Outcome vector on the `[0, 1]` scale (binary `0`/`1`, or a
#'   continuous outcome rescaled to `[0, 1]`); see [fit_Q()] Details.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param h Numeric tilting-function weight vector.
#' @param g1W,g0W Fitted treatment propensity at `A = 1` (and, by default,
#'   `1 - g1W` at `A = 0`).
#' @param Delta Indicator that `Y` is observed. Default `1 - is.na(Y)`; `Y`
#'   may be `NA` where `Delta = 0`.
#' @param pi1W,pi0W Missing-outcome propensity `P(Delta = 1 | A = a, W)` at
#'   `A = 1` and `A = 0` (see [fit_pi()]). Default `1`, i.e. no missing outcomes.
#'
#' @return Numeric point estimate of \eqn{\tau_h}.
#' @seealso [mle_h()], [aipw_h()] for the other benchmark estimators,
#'   [tmle_h()] for the targeted (TMLE) estimator.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' gfit <- fit_g(A, W)
#' ipw_h(Y, A, h = rep(1, n), gfit$g1W)
#' @export
ipw_h <- function(Y, A, h, g1W, g0W = 1 - g1W,
                  Delta = as.integer(!is.na(Y)),
                  pi1W = rep(1, length(Y)), pi0W = rep(1, length(Y))) {
  Y0 <- ifelse(Delta == 1, Y, 0)   # missing outcomes get weight 0 through Delta
  sum(h * (Delta * A * Y0 / (g1W * pi1W) - Delta * (1 - A) * Y0 / (g0W * pi0W))) / sum(h)
}

#' Augmented inverse-probability-weighted (AIPW) estimator of a tilted ATE
#'
#' Doubly-robust benchmark estimator, with an EIF-based standard error of
#' the same form used by [tmle_h()].
#'
#' @inheritParams ipw_h
#' @param Q1W,Q0W Fitted outcome regression at `A = 1` and `A = 0`.
#'
#' @return A list with components:
#' \item{estimate}{Point estimate of \eqn{\tau_h}.}
#' \item{se}{EIF-based standard error.}
#'
#' @seealso [mle_h()], [ipw_h()] for the other benchmark estimators,
#'   [tmle_h()] for the targeted (TMLE) estimator.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' Qfit <- fit_Q(Y, A, W)
#' gfit <- fit_g(A, W)
#' aipw_h(Y, A, h = rep(1, n), Qfit$Q1W, Qfit$Q0W, gfit$g1W)
#' @export
aipw_h <- function(Y, A, h, Q1W, Q0W, g1W, g0W = 1 - g1W,
                   Delta = as.integer(!is.na(Y)),
                   pi1W = rep(1, length(Y)), pi0W = rep(1, length(Y))) {
  N <- length(Y)
  h_star <- h / mean(h)
  Y0 <- ifelse(Delta == 1, Y, 0)
  # Residual terms only count where the outcome is observed (Delta = 1),
  # weighted by 1 / (P(A = a | W) * P(Delta = 1 | A = a, W)).
  terms <- Q1W - Q0W +
    Delta * (A / (g1W * pi1W) * (Y0 - Q1W) - (1 - A) / (g0W * pi0W) * (Y0 - Q0W))
  estimate <- mean(h_star * terms)
  EIF <- h_star * (terms - estimate)
  list(estimate = estimate, se = sqrt(mean(EIF^2) / N))
}

#' Fit the TMLE target-step fluctuation parameter
#'
#' The shared target-step logistic regression
#' \eqn{\mathrm{logit}[Q(X,A;\epsilon)] = \mathrm{logit}[\hat Q(X,A)] + \epsilon H(X,A)}
#' used throughout the package's TMLE implementations. Factored out so
#' [tmle_h()] and [tmle_bootstrap_h()] (and the internally/externally
#' tilted TMLE cores) can all call it on different subsets of rows
#' without duplicating logic.
#'
#' @param Y Response vector with values in `[0, 1]` -- the outcome `Y`
#'   (binary, or continuous rescaled to `[0, 1]`) when targeting `Q`, or a
#'   selection/treatment indicator (`D`/`A`) when targeting `g`/`s`. The
#'   fluctuation is fit by quasi-binomial regression, which gives the usual
#'   logistic-regression estimate for binary `Y` and remains valid for
#'   continuous `Y` in `[0, 1]`.
#' @param H The clever covariate, same length as `Y`.
#' @param logitQAW Offset: the current fit on the logit scale, same length as `Y`.
#' @param indices Row indices of `Y`/`H`/`logitQAW` to fit on (defaults to all).
#'
#' @return Numeric scalar: the fitted fluctuation coefficient `eps`, or
#'   `0` if the fit is degenerate or non-finite (guards against
#'   degenerate/separated fits, e.g. in a bootstrap resample).
#' @seealso [tmle_h()], [tmle_bootstrap_h()] for callers that target `Q`;
#'   [tmle_tilted_h()] and [tmle_external_tilted_h()] for callers that
#'   also use this to target `g`/`s`.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' Qfit <- fit_Q(Y, A, W)
#' gfit <- fit_g(A, W)
#' H <- A / gfit$g1W - (1 - A) / gfit$g0W
#' fit_eps_h(Y, H, stats::qlogis(Qfit$QAW))
#' @export
fit_eps_h <- function(Y, H, logitQAW, indices = seq_along(Y)) {
  dat <- data.frame(Y = Y[indices], H = H[indices], off = logitQAW[indices])
  # quasibinomial: same estimate as binomial for 0/1 Y, and also valid for Y in [0, 1]
  fit <- suppressWarnings(stats::glm(Y ~ -1 + H + offset(off), data = dat,
                                     family = stats::quasibinomial()))
  eps <- unname(stats::coef(fit)["H"])
  if (!is.finite(eps)) eps <- 0
  eps
}

#' Targeted maximum likelihood estimator (TMLE) of a tilted ATE for a fixed h
#'
#' Implements the initial step + target step (single logistic
#' fluctuation) and EIF-based standard error for a tilting function `h`
#' that is either a fixed, precomputed vector or has already been
#' evaluated at a fitted propensity score (see [estimate_tau_h_fixed()]
#' for the higher-level entry point that handles both cases).
#'
#' The target step runs one logistic regression of `Y` on the clever
#' covariate `H = h_star * (A/g1W - (1-A)/g0W)` (with `h_star = h /
#' mean(h)`), with the current fit on the logit scale as an offset (see
#' [fit_eps_h()]). `H`/`logitQAW` are returned so [tmle_bootstrap_h()] can
#' reuse them exactly on each bootstrap resample, instead of recomputing.
#'
#' @param Y Outcome vector on the `[0, 1]` scale (binary `0`/`1`, or a
#'   continuous outcome rescaled to `[0, 1]`); see [fit_Q()] Details.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param h Numeric tilting-function weight vector.
#' @param QAW Fitted outcome regression at the observed treatment.
#' @param Q1W,Q0W Fitted outcome regression at `A = 1` and `A = 0`.
#' @param g1W,g0W Fitted treatment propensity at `A = 1` (and, by default,
#'   `1 - g1W` at `A = 0`).
#' @param Delta,pi1W,pi0W Missing-outcome indicator and propensity
#'   `P(Delta = 1 | A = a, W)`; see [ipw_h()]. With missing outcomes the clever
#'   covariate uses `g * pi` in place of `g`, the target step is fit on the rows
#'   with `Delta = 1`, and the EIF residual is multiplied by `Delta`.
#'
#' @return A list with components:
#' \item{estimate}{TMLE point estimate of \eqn{\tau_h}.}
#' \item{se}{EIF-based standard error.}
#' \item{CI}{Wald 95% confidence interval, `estimate +/- 1.96 * se`.}
#' \item{eps}{Fitted target-step fluctuation coefficient.}
#' \item{Q1W_star,Q0W_star}{Targeted outcome regression at `A = 1` and `A = 0`.}
#' \item{EIF}{Efficient influence function, evaluated at each observation.}
#' \item{H}{The clever covariate used in the target step; reused by [tmle_bootstrap_h()].}
#' \item{logitQAW}{The initial fit on the logit scale (the fluctuation's offset); reused by [tmle_bootstrap_h()].}
#'
#' @seealso [fit_eps_h()] for the underlying fluctuation fitter,
#'   [tmle_bootstrap_h()] for bootstrap inference reusing this fit,
#'   [estimate_tau_h_fixed()] for the typical calling context.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' Qfit <- fit_Q(Y, A, W)
#' gfit <- fit_g(A, W)
#' out <- tmle_h(Y, A, h = rep(1, n), Qfit$QAW, Qfit$Q1W, Qfit$Q0W, gfit$g1W)
#' out$estimate
#' out$CI
#' @export
tmle_h <- function(Y, A, h, QAW, Q1W, Q0W, g1W, g0W = 1 - g1W,
                   Delta = as.integer(!is.na(Y)),
                   pi1W = rep(1, length(Y)), pi0W = rep(1, length(Y))) {
  N <- length(Y)
  h_star <- h / mean(h)
  Y0 <- ifelse(Delta == 1, Y, 0)         # Y with missing entries set to 0 (they are never used)
  gpi1 <- g1W * pi1W                     # P(A = 1 and observed | W)
  gpi0 <- g0W * pi0W                     # P(A = 0 and observed | W)
  HAW <- A / gpi1 - (1 - A) / gpi0
  H <- h_star * HAW
  logitQAW <- stats::qlogis(QAW)

  # Target step: fit eps on the rows with an observed outcome only
  eps <- fit_eps_h(Y, H, logitQAW, indices = which(Delta == 1))

  QAW_star <- stats::plogis(logitQAW + eps * H)
  Q1W_star <- stats::plogis(stats::qlogis(Q1W) + eps * (h_star / gpi1))
  Q0W_star <- stats::plogis(stats::qlogis(Q0W) - eps * (h_star / gpi0))

  # Plug-in over ALL rows (Q is predicted for everyone, including missing Y)
  estimate <- sum(h * (Q1W_star - Q0W_star)) / sum(h)
  EIF <- h_star * (Q1W_star - Q0W_star - estimate) + Delta * h_star * HAW * (Y0 - QAW_star)
  se <- sqrt(mean(EIF^2) / N)
  ci <- estimate + c(-1, 1) * stats::qnorm(0.975) * se

  list(estimate = estimate, se = se, CI = ci, eps = eps,
       Q1W_star = Q1W_star, Q0W_star = Q0W_star, EIF = EIF,
       H = H, logitQAW = logitQAW)
}

#' Target-step-only stratified bootstrap for a tilted-ATE TMLE
#'
#' Initial nuisance fits (`Q1W`, `Q0W`, `g1W`, `h`) are held fixed; only
#' the target-step fluctuation `eps` is refit on each resample.
#'
#' @param Y Outcome vector on the `[0, 1]` scale (binary `0`/`1`, or a
#'   continuous outcome rescaled to `[0, 1]`); see [fit_Q()] Details.
#' @param A Binary (`0`/`1`) treatment vector.
#' @param h Numeric tilting-function weight vector.
#' @param Q1W,Q0W Fitted outcome regression at `A = 1` and `A = 0`.
#' @param g1W,g0W Fitted treatment propensity at `A = 1` (and, by default,
#'   `1 - g1W` at `A = 0`).
#' @param tmle_fit A fitted-TMLE object, as returned by [tmle_h()], whose
#'   `H`/`logitQAW` components are reused on each resample.
#' @param n_boot Number of bootstrap resamples.
#' @param strata Vector controlling what's held fixed in resample
#'   composition (defaults to arm-only stratification via `A`).
#' @param subset Row indices eligible for resampling at all (defaults to
#'   everything).
#' @param Delta,pi1W,pi0W Missing-outcome indicator and propensity, as in
#'   [tmle_h()]. Whole rows are resampled; `eps` is refit on the resampled rows
#'   with `Delta = 1`.
#'
#' @return A list with components:
#' \item{estimates}{Numeric vector of length `n_boot`, the bootstrap
#'   distribution of \eqn{\hat\tau_h}.}
#' \item{se}{Standard deviation of `estimates`.}
#' \item{CI_normal}{Wald 95% confidence interval centered at `tmle_fit$estimate`,
#'   using the bootstrap `se`.}
#' \item{CI_percentile}{Percentile-based 95% confidence interval, the
#'   2.5th/97.5th quantiles of `estimates`.}
#'
#' @seealso [tmle_h()], whose fit this function resamples from.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' Qfit <- fit_Q(Y, A, W)
#' gfit <- fit_g(A, W)
#' h <- rep(1, n)
#' tmle_out <- tmle_h(Y, A, h, Qfit$QAW, Qfit$Q1W, Qfit$Q0W, gfit$g1W)
#' boot_out <- tmle_bootstrap_h(Y, A, h, Qfit$Q1W, Qfit$Q0W, gfit$g1W,
#'                               tmle_fit = tmle_out, n_boot = 50)
#' boot_out$se
#' boot_out$CI_normal
#' @export
tmle_bootstrap_h <- function(Y, A, h, Q1W, Q0W, g1W, g0W = 1 - g1W,
                              tmle_fit, n_boot = 300,
                              strata = A, subset = seq_along(Y),
                              Delta = as.integer(!is.na(Y)),
                              pi1W = rep(1, length(Y)), pi0W = rep(1, length(Y))) {
  h_star <- h / mean(h)
  H <- tmle_fit$H
  logitQAW <- tmle_fit$logitQAW
  gpi1 <- g1W * pi1W
  gpi0 <- g0W * pi0W

  strata_id <- split(subset, strata[subset])

  boot_est <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    # sample(idx, ...) would misbehave for a length-1 stratum, so index explicitly
    boot_id <- unlist(lapply(strata_id, function(idx) idx[sample.int(length(idx), replace = TRUE)]))
    obs_id <- boot_id[Delta[boot_id] == 1]   # only observed outcomes inform eps
    eps_b <- fit_eps_h(Y, H, logitQAW, indices = obs_id)
    Q1_b <- stats::plogis(stats::qlogis(Q1W) + eps_b * (h_star / gpi1))
    Q0_b <- stats::plogis(stats::qlogis(Q0W) - eps_b * (h_star / gpi0))
    boot_est[b] <- sum(h * (Q1_b - Q0_b)) / sum(h)
  }

  se <- stats::sd(boot_est)
  list(estimates = boot_est, se = se,
       CI_normal = tmle_fit$estimate + c(-1, 1) * stats::qnorm(0.975) * se,
       CI_percentile = unname(stats::quantile(boot_est, c(0.025, 0.975))))
}
