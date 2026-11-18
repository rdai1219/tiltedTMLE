# Full-refit nonparametric bootstrap for the internally / externally tilted TMLEs.
# Each resample re-fits Q, g (and s) AND re-runs the iterative targeting, so the
# bootstrap captures the variability of the nuisance estimation (unlike the
# target-step-only bootstrap used for fixed h in [tmle_bootstrap_h()]).

# Internal: stratified resample of row positions 1..length(strata), with replacement.
boot_indices <- function(strata) {
  unlist(lapply(split(seq_along(strata), strata),
                function(i) i[sample.int(length(i), replace = TRUE)]),
         use.names = FALSE)
}

# Internal: turn a vector of bootstrap estimates into SE + two CIs.
boot_summary <- function(est, point, conf_level) {
  ok <- is.finite(est)
  e <- est[ok]
  if (length(e) < 2) stop("Fewer than 2 successful bootstrap resamples.", call. = FALSE)
  alpha <- 1 - conf_level
  se <- stats::sd(e)
  list(estimate = point,
       estimates = est,
       se = se,
       CI_normal = point + c(-1, 1) * stats::qnorm(1 - alpha / 2) * se,
       CI_percentile = unname(stats::quantile(e, c(alpha / 2, 1 - alpha / 2))),
       n_failed = sum(!ok))
}

#' Full-refit bootstrap for an internally tilted TMLE
#'
#' Resamples the rows of `(Y, A, W)` with replacement (stratified by `strata`,
#' default the treatment arm `A`), and on every resample re-fits `Q` and `g`
#' and re-runs the iterative TMLE via [estimate_tau_h_internal_tilted()]. The
#' bootstrap SE therefore reflects nuisance-estimation variability, which the
#' EIF-based SE and the target-step bootstrap of [tmle_bootstrap_h()] do not
#' fully capture when flexible learners are used.
#'
#' @param Y,A,W,eta,eta_dot As in [estimate_tau_h_internal_tilted()].
#' @param n_boot Number of bootstrap resamples.
#' @param strata Vector (length `n`) defining resampling strata; default `A`.
#' @param conf_level Confidence level, default 0.95.
#' @param ... Further arguments passed to [estimate_tau_h_internal_tilted()]
#'   (e.g. `Q_method`, `gbounds`, `Q_args`).
#'
#' @return A list with components:
#' \item{estimate}{TMLE point estimate on the original data.}
#' \item{estimates}{Length-`n_boot` bootstrap TMLE estimates (`NA` for failed resamples).}
#' \item{se}{Bootstrap standard error (SD of the successful `estimates`).}
#' \item{CI_normal}{Wald interval centered at `estimate` using the bootstrap `se`.}
#' \item{CI_percentile}{Percentile interval from the bootstrap distribution.}
#' \item{n_failed}{Number of resamples that errored or gave a non-finite estimate.}
#' \item{n_nonconverged}{Number of successful resamples whose TMLE hit `max_iter`.}
#'
#' @seealso [bootstrap_tau_h_external_tilted()] for the external analogue.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' spec <- tilt_spec("ATT")
#' b <- bootstrap_tau_h_internal_tilted(Y, A, W, spec$eta, spec$eta_dot, n_boot = 20)
#' b$se
#' b$CI_percentile
#' @export
bootstrap_tau_h_internal_tilted <- function(Y, A, W, eta, eta_dot, n_boot = 200,
                                             strata = A, conf_level = 0.95, ...) {
  W <- as.data.frame(W)
  point <- estimate_tau_h_internal_tilted(Y, A, W, eta, eta_dot, ...)$TMLE

  est <- rep(NA_real_, n_boot)
  nonconv <- 0L
  for (b in seq_len(n_boot)) {
    id <- boot_indices(strata)
    # A failed resample (e.g. separation in glm) is recorded as NA, not fatal
    out <- tryCatch(
      suppressWarnings(estimate_tau_h_internal_tilted(Y[id], A[id], W[id, , drop = FALSE],
                                                      eta, eta_dot, ...)),
      error = function(e) NULL)
    if (!is.null(out)) {
      est[b] <- out$TMLE
      if (!isTRUE(out$TMLE_converged)) nonconv <- nonconv + 1L
    }
  }

  res <- boot_summary(est, point, conf_level)
  res$n_nonconverged <- nonconv
  res
}

#' Full-refit bootstrap for an externally tilted TMLE
#'
#' Resamples the internal rows (stratified by `strata`, default the treatment
#' arm `A`) and the external rows separately, with replacement, keeping the
#' sample sizes `n` and `m` fixed. On every resample `Q`, `g`, `s` are re-fit
#' and the iterative TMLE is re-run via [estimate_tau_h_external_tilted()].
#'
#' @param Y,A Length-`n` internal-only outcome and treatment.
#' @param D,Ws,eta,eta_dot As in [estimate_tau_h_external_tilted()].
#' @param Wg Optional internal-only covariates for `g`, as in
#'   [estimate_tau_h_external_tilted()].
#' @param n_boot,conf_level See [bootstrap_tau_h_internal_tilted()].
#' @param strata Length-`n` vector for stratifying the internal resample; default `A`.
#' @param ... Further arguments passed to [estimate_tau_h_external_tilted()].
#'
#' @return Same components as [bootstrap_tau_h_internal_tilted()].
#' @seealso [bootstrap_tau_h_internal_tilted()].
#' @examples
#' set.seed(1)
#' n <- 150; m <- 150
#' Ws <- data.frame(X1 = rnorm(n + m))
#' D <- c(rep(1L, n), rep(0L, m))
#' A <- rbinom(n, 1, stats::plogis(0.3 * Ws$X1[1:n]))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * Ws$X1[1:n]))
#' spec <- tilt_spec_external("ATE-O")
#' b <- bootstrap_tau_h_external_tilted(Y, A, D, Ws, spec$eta, spec$eta_dot, n_boot = 20)
#' b$se
#' @export
bootstrap_tau_h_external_tilted <- function(Y, A, D, Ws, eta, eta_dot, Wg = NULL,
                                             n_boot = 200, strata = A,
                                             conf_level = 0.95, ...) {
  Ws <- as.data.frame(Ws)
  idx1 <- which(D == 1)
  idx0 <- which(D == 0)
  point <- estimate_tau_h_external_tilted(Y, A, D, Ws, eta, eta_dot, Wg = Wg, ...)$TMLE

  est <- rep(NA_real_, n_boot)
  nonconv <- 0L
  for (b in seq_len(n_boot)) {
    p <- boot_indices(strata)                          # positions among internal rows
    q <- idx0[sample.int(length(idx0), replace = TRUE)]  # external rows
    Ws_b <- Ws[c(idx1[p], q), , drop = FALSE]
    D_b <- c(rep(1L, length(p)), rep(0L, length(q)))
    Wg_b <- if (is.null(Wg)) NULL else as.data.frame(Wg)[p, , drop = FALSE]
    out <- tryCatch(
      suppressWarnings(estimate_tau_h_external_tilted(Y[p], A[p], D_b, Ws_b, eta, eta_dot,
                                                      Wg = Wg_b, ...)),
      error = function(e) NULL)
    if (!is.null(out)) {
      est[b] <- out$TMLE
      if (!isTRUE(out$TMLE_converged)) nonconv <- nonconv + 1L
    }
  }

  res <- boot_summary(est, point, conf_level)
  res$n_nonconverged <- nonconv
  res
}
