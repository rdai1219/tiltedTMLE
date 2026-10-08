#' Fit the outcome regression Q(x, a)
#'
#' Fits an outcome regression of a binary `Y` on treatment `A` and
#' covariates `W`, via logistic regression or [SuperLearner::SuperLearner()],
#' and returns predicted values at the observed treatment, at `A = 1`,
#' and at `A = 0`.
#'
#' @param Y Binary outcome vector.
#' @param A Binary treatment vector.
#' @param W Covariate data frame (or matrix, coerced via [as.data.frame()]).
#' @param method `"glm"` (default) or `"SL"` for [SuperLearner::SuperLearner()].
#' @param Qform Optional model formula for `method = "glm"`. Defaults to
#'   `Y ~ A * (W1 + W2 + ...)` over all columns of `W`.
#' @param SL.library Character vector of SuperLearner wrapper names, used
#'   when `method = "SL"`.
#' @param ... Passed through to `glm()` or `SuperLearner::SuperLearner()`
#'   (e.g. `cvControl` for SL).
#'
#' @return A list with components `QAW`, `Q1W`, `Q0W` (predicted
#'   probabilities at the observed, at `A = 1`, and at `A = 0` treatment,
#'   respectively), `fit` (the underlying fitted model object), and
#'   `method`.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' Qfit <- fit_Q(Y, A, W)
#' str(Qfit[c("QAW", "Q1W", "Q0W")])
#' @export
fit_Q <- function(Y, A, W, method = c("glm", "SL"),
                   Qform = NULL, SL.library = c("SL.glm"), ...) {
  method <- match.arg(method)
  W <- as.data.frame(W)

  if (method == "glm") {
    if (is.null(Qform)) Qform <- stats::as.formula(paste("Y ~ A*(", paste(names(W), collapse = "+"), ")"))
    fit <- stats::glm(Qform, data = data.frame(Y = Y, A = A, W), family = stats::binomial())
    QAW <- stats::predict(fit, newdata = data.frame(A = A, W), type = "response")
    Q1W <- stats::predict(fit, newdata = data.frame(A = 1, W), type = "response")
    Q0W <- stats::predict(fit, newdata = data.frame(A = 0, W), type = "response")
  } else {
    if (!requireNamespace("SuperLearner", quietly = TRUE)) {
      stop("Package 'SuperLearner' is required for method = \"SL\". ",
           "Install it with install.packages(\"SuperLearner\").", call. = FALSE)
    }
    WA <- data.frame(A = A, W)
    fit <- SuperLearner::SuperLearner(Y = Y, X = WA, family = stats::binomial(),
                                       SL.library = SL.library, ...)
    QAW <- as.numeric(stats::predict(fit, newdata = WA)$pred)
    Q1W <- as.numeric(stats::predict(fit, newdata = data.frame(A = 1, W))$pred)
    Q0W <- as.numeric(stats::predict(fit, newdata = data.frame(A = 0, W))$pred)
  }

  list(QAW = QAW, Q1W = Q1W, Q0W = Q0W, fit = fit, method = method)
}

#' Fit a propensity (or selection) score g(x) = P(A = 1 | X)
#'
#' Fits a binary treatment (or, when called with `D` in place of `A`, a
#' selection indicator) on covariates `W`, via logistic regression or
#' [SuperLearner::SuperLearner()], with optional truncation of the fitted
#' probabilities away from 0 and 1.
#'
#' @param A Binary treatment (or selection) indicator vector.
#' @param W Covariate data frame (or matrix, coerced via [as.data.frame()]).
#' @param method `"glm"` (default) or `"SL"` for [SuperLearner::SuperLearner()].
#' @param gform Optional model formula for `method = "glm"`. Defaults to
#'   `A ~ W1 + W2 + ...` over all columns of `W`.
#' @param SL.library Character vector of SuperLearner wrapper names, used
#'   when `method = "SL"`.
#' @param gbounds Length-2 vector of truncation bounds applied to the
#'   fitted probabilities. `c(0.025, 0.975)` (default) truncates;
#'   `c(0, 1)` is a no-op clamp (fitted probabilities from `glm`/`SL` are
#'   always strictly inside `(0, 1)`), which is how "no truncation" is
#'   realized by callers.
#' @param ... Passed through to `glm()` or `SuperLearner::SuperLearner()`.
#'
#' @return A list with components `g1W` (truncated fitted `P(A = 1 | W)`),
#'   `g0W` (`1 - g1W`), `fit` (the underlying fitted model object), and
#'   `method`.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' gfit <- fit_g(A, W)
#' summary(gfit$g1W)
#' @export
fit_g <- function(A, W, method = c("glm", "SL"),
                   gform = NULL, SL.library = c("SL.glm"), gbounds = c(0.025, 0.975), ...) {
  method <- match.arg(method)
  W <- as.data.frame(W)

  if (method == "glm") {
    if (is.null(gform)) gform <- stats::as.formula(paste("A ~", paste(names(W), collapse = "+")))
    fit <- stats::glm(gform, data = data.frame(A = A, W), family = stats::binomial())
    g1W <- stats::predict(fit, newdata = W, type = "response")
  } else {
    if (!requireNamespace("SuperLearner", quietly = TRUE)) {
      stop("Package 'SuperLearner' is required for method = \"SL\". ",
           "Install it with install.packages(\"SuperLearner\").", call. = FALSE)
    }
    fit <- SuperLearner::SuperLearner(Y = A, X = W, family = stats::binomial(),
                                       SL.library = SL.library, ...)
    g1W <- as.numeric(stats::predict(fit, newdata = W)$pred)
  }

  g1W <- pmin(pmax(g1W, gbounds[1]), gbounds[2])
  list(g1W = g1W, g0W = 1 - g1W, fit = fit, method = method)
}

#' Predict P(A = 1 | X) from a fitted propensity-score model
#'
#' Dispatches prediction on either a `glm` or a `SuperLearner` fit object
#' returned by [fit_g()], since `predict.SuperLearner()` doesn't accept
#' `type = "response"` and returns a list (unlike `predict.glm()`). Lets
#' any caller (e.g. extrapolating a fitted g onto new covariates, or a
#' simulation's own Monte Carlo truth calculation) predict from a
#' [fit_g()] result without caring which method produced it.
#'
#' @param g_fit The `fit` component of a [fit_g()] result (a `glm` or
#'   `SuperLearner` object).
#' @param newdata Covariate data frame to predict at.
#'
#' @return Numeric vector of predicted probabilities.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' gfit <- fit_g(A, W)
#' predict_ghat(gfit$fit, newdata = W[1:5, ])
#' @export
predict_ghat <- function(g_fit, newdata) {
  if (inherits(g_fit, "SuperLearner")) {
    as.numeric(stats::predict(g_fit, newdata = newdata)$pred)
  } else {
    as.numeric(stats::predict(g_fit, newdata = newdata, type = "response"))
  }
}

#' Predict Q(x, a) at a fixed treatment arm from a fitted outcome-regression model
#'
#' The outcome-regression analog of [predict_ghat()]: dispatches on
#' whether the fit object is a `glm` or `SuperLearner`, and predicts at a
#' fixed arm `A_new` (either a scalar, recycled, or a vector of the same
#' length as `nrow(W_new)`) over new covariates. Used to extrapolate an
#' outcome regression fit on one subset of rows (e.g. an internal trial
#' arm) onto a larger pooled covariate set (e.g. including external,
#' covariate-only rows).
#'
#' @param Q_fit The `fit` component of a [fit_Q()] result (a `glm` or
#'   `SuperLearner` object).
#' @param A_new Fixed treatment arm (`0` or `1`) to predict at.
#' @param W_new Covariate data frame to predict at.
#'
#' @return Numeric vector of predicted probabilities.
#' @examples
#' set.seed(1)
#' n <- 200
#' W <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
#' A <- rbinom(n, 1, stats::plogis(0.3 * W$X1))
#' Y <- rbinom(n, 1, stats::plogis(-0.3 + 0.5 * A + 0.3 * W$X1))
#' Qfit <- fit_Q(Y, A, W)
#' predict_Qhat(Qfit$fit, A_new = 1, W_new = W[1:5, ])
#' @export
predict_Qhat <- function(Q_fit, A_new, W_new) {
  W_new <- as.data.frame(W_new)
  if (inherits(Q_fit, "SuperLearner")) {
    as.numeric(stats::predict(Q_fit, newdata = data.frame(A = A_new, W_new))$pred)
  } else {
    as.numeric(stats::predict(Q_fit, newdata = data.frame(A = A_new, W_new), type = "response"))
  }
}
