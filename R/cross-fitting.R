# Cross-fitting (sample splitting) for the nuisance models.
#
# Without cross-fitting, Q, g (and s, pi) are fit on all the data and then used,
# on that same data, in the targeting step and the EIF-based standard error. With
# flexible learners (e.g. a SuperLearner library with random forests or boosting)
# the fits can partly memorise the data, which makes the residuals too small and
# the confidence intervals too narrow.
#
# With cv_folds = V > 1 the rows are split into V folds. For each fold v, every
# nuisance model is fit on the other V - 1 folds and used to predict the rows in
# fold v ("out-of-fold" predictions). The usual targeting step and EIF-based
# inference are then run on these out-of-fold predictions, pooled over all rows
# (as in CV-TMLE; Zheng and van der Laan, 2011). Folds are stratified so that
# every fold has a similar mix of treated / control (and internal / external) rows.
# Fold assignment is random: use set.seed() for reproducible results.

# Internal: check the cv_folds argument (1 = no cross-fitting).
check_cv_folds <- function(cv_folds, n) {
  if (length(cv_folds) != 1 || !is.numeric(cv_folds) || is.na(cv_folds) ||
      cv_folds != round(cv_folds) || cv_folds < 1 || cv_folds > n) {
    stop("`cv_folds` must be a whole number between 1 (no cross-fitting) and the ",
         "number of rows.", call. = FALSE)
  }
  invisible(TRUE)
}

# Internal: assign each row to one of V folds, balancing the folds within each
# level of `strata`. Returns an integer vector of fold ids 1..V.
make_folds <- function(strata, V) {
  folds <- integer(length(strata))
  for (s in unique(strata)) {
    rows <- which(strata == s)
    folds[rows] <- sample(rep_len(seq_len(V), length(rows)))
  }
  folds
}

# Internal: generic out-of-fold prediction.
# For each fold v, fit_fun(train_rows) fits a model on the rows NOT in fold v and
# pred_fun(fit, test_rows) predicts for the rows IN fold v (a vector, or a matrix
# with one column per prediction). Returns the predictions for all rows (as a
# matrix) and the list of V fitted models.
cross_predict <- function(folds, fit_fun, pred_fun) {
  V <- max(folds)
  fits <- vector("list", V)
  pred <- NULL
  for (v in seq_len(V)) {
    test <- which(folds == v)
    train <- which(folds != v)
    fits[[v]] <- fit_fun(train)
    p <- as.matrix(pred_fun(fits[[v]], test))
    if (is.null(pred)) pred <- matrix(NA_real_, length(folds), ncol(p))
    pred[test, ] <- p
  }
  list(pred = pred, fits = fits)
}

# Internal: cross-fitted Q, g and (if Y has NAs) pi for the fixed-h and internally
# tilted entry points. Returns list(Qfit, gfit, pifit), shaped like the
# full-data fits, except that `fit` holds the list of V per-fold models.
crossfit_nuisance_internal <- function(Y, A, W, folds,
                                       Q_method, Qform, Q_SL.library, Q_args,
                                       g_method, gform, g_SL.library, gbounds, g_args,
                                       pi_method, piform, pi_SL.library, pibounds, pi_args) {
  # Outcome regression: fit on the observed outcomes of the training folds,
  # predict Q(1, W) and Q(0, W) for the held-out fold
  Qcv <- cross_predict(folds,
    fit_fun = function(tr) {
      fit_Q_observed(Y[tr], A[tr], W[tr, , drop = FALSE], method = Q_method, Qform = Qform,
                     SL.library = Q_SL.library, Q_args = Q_args)$fit
    },
    pred_fun = function(fit, te) {
      cbind(predict_Qhat(fit, 1, W[te, , drop = FALSE]), predict_Qhat(fit, 0, W[te, , drop = FALSE]))
    })
  Q1W <- Qcv$pred[, 1]
  Q0W <- Qcv$pred[, 2]
  Qfit <- list(QAW = A * Q1W + (1 - A) * Q0W, Q1W = Q1W, Q0W = Q0W,
               fit = Qcv$fits, method = Q_method)

  # Treatment propensity score g
  gcv <- cross_predict(folds,
    fit_fun = function(tr) {
      do.call(fit_g, c(list(A[tr], W[tr, , drop = FALSE], method = g_method, gform = gform,
                            SL.library = g_SL.library, gbounds = gbounds), g_args))$fit
    },
    pred_fun = function(fit, te) predict_ghat(fit, W[te, , drop = FALSE]))
  g1W <- clamp_bounds(gcv$pred[, 1], gbounds)
  gfit <- list(g1W = g1W, g0W = 1 - g1W, fit = gcv$fits, method = g_method)

  # Missingness propensity pi_a(W) = P(Delta = 1 | A = a, W), only if some Y are missing
  pifit <- NULL
  if (anyNA(Y)) {
    Delta <- as.integer(!is.na(Y))
    picv <- cross_predict(folds,
      fit_fun = function(tr) {
        do.call(fit_pi, c(list(Delta[tr], A[tr], W[tr, , drop = FALSE], method = pi_method,
                               piform = piform, SL.library = pi_SL.library,
                               pibounds = pibounds), pi_args))$fit
      },
      pred_fun = function(fit, te) {
        cbind(predict_Qhat(fit, 1, W[te, , drop = FALSE], Qbounds = c(0, 1)),
              predict_Qhat(fit, 0, W[te, , drop = FALSE], Qbounds = c(0, 1)))
      })
    pi1W <- clamp_bounds(picv$pred[, 1], pibounds)
    pi0W <- clamp_bounds(picv$pred[, 2], pibounds)
    pifit <- list(pi1W = pi1W, pi0W = pi0W, piAW = A * pi1W + (1 - A) * pi0W,
                  fit = picv$fits, method = pi_method)
  }

  list(Qfit = Qfit, gfit = gfit, pifit = pifit)
}

# Internal: cross-fitted Q, g, s and (if Y has NAs) pi for the externally tilted
# entry point. Folds are over all N pooled rows. In each training set, Q, g and pi
# are fit on its internal rows and s on all its rows; predictions are made for
# every pooled row (internal and external) in the held-out fold.
# Returns list(Qfit, gfit, sfit, pifit, g1W_all), shaped like the full-data fits.
crossfit_nuisance_external <- function(Y, A, D, Ws, Wg, folds,
                                       Q_method, Qform, Q_SL.library, Q_args,
                                       g_method, gform, g_SL.library, gbounds, g_args,
                                       s_method, sform, s_SL.library, sbounds, s_args,
                                       pi_method, piform, pi_SL.library, pibounds, pi_args) {
  idx1 <- which(D == 1)
  Wg_all <- Ws[, names(Wg), drop = FALSE]
  # Positions (within the internal-only Y, A, Wg) of the internal rows among `rows`
  internal_pos <- function(rows) match(intersect(rows, idx1), idx1)

  # Outcome regression: internal training rows with an observed outcome
  Qcv <- cross_predict(folds,
    fit_fun = function(tr) {
      p <- internal_pos(tr)
      fit_Q_observed(Y[p], A[p], Ws[idx1[p], , drop = FALSE], method = Q_method, Qform = Qform,
                     SL.library = Q_SL.library, Q_args = Q_args)$fit
    },
    pred_fun = function(fit, te) {
      cbind(predict_Qhat(fit, 1, Ws[te, , drop = FALSE]), predict_Qhat(fit, 0, Ws[te, , drop = FALSE]))
    })
  Q1W <- Qcv$pred[, 1]
  Q0W <- Qcv$pred[, 2]
  Qfit <- list(QAW = A * Q1W[idx1] + (1 - A) * Q0W[idx1], Q1W = Q1W, Q0W = Q0W,
               fit = Qcv$fits, method = Q_method)

  # Treatment propensity g: internal training rows, predicted for all pooled rows
  gcv <- cross_predict(folds,
    fit_fun = function(tr) {
      p <- internal_pos(tr)
      do.call(fit_g, c(list(A[p], Wg[p, , drop = FALSE], method = g_method, gform = gform,
                            SL.library = g_SL.library, gbounds = gbounds), g_args))$fit
    },
    pred_fun = function(fit, te) predict_ghat(fit, Wg_all[te, , drop = FALSE]))
  g1W_all <- clamp_bounds(gcv$pred[, 1], gbounds)
  gfit <- list(g1W = g1W_all[idx1], g0W = 1 - g1W_all[idx1], fit = gcv$fits, method = g_method)

  # Selection propensity s: all training rows (internal and external)
  scv <- cross_predict(folds,
    fit_fun = function(tr) {
      do.call(fit_g, c(list(D[tr], Ws[tr, , drop = FALSE], method = s_method, gform = sform,
                            SL.library = s_SL.library, gbounds = sbounds), s_args))$fit
    },
    pred_fun = function(fit, te) predict_ghat(fit, Ws[te, , drop = FALSE]))
  s1W <- clamp_bounds(scv$pred[, 1], sbounds)
  sfit <- list(g1W = s1W, g0W = 1 - s1W, fit = scv$fits, method = s_method)

  # Missingness propensity pi: internal training rows, predicted for all pooled rows
  pifit <- NULL
  if (anyNA(Y)) {
    Delta <- as.integer(!is.na(Y))
    picv <- cross_predict(folds,
      fit_fun = function(tr) {
        p <- internal_pos(tr)
        do.call(fit_pi, c(list(Delta[p], A[p], Ws[idx1[p], , drop = FALSE], method = pi_method,
                               piform = piform, SL.library = pi_SL.library,
                               pibounds = pibounds), pi_args))$fit
      },
      pred_fun = function(fit, te) {
        cbind(predict_Qhat(fit, 1, Ws[te, , drop = FALSE], Qbounds = c(0, 1)),
              predict_Qhat(fit, 0, Ws[te, , drop = FALSE], Qbounds = c(0, 1)))
      })
    pi1W_all <- clamp_bounds(picv$pred[, 1], pibounds)
    pi0W_all <- clamp_bounds(picv$pred[, 2], pibounds)
    pifit <- list(pi1W = pi1W_all[idx1], pi0W = pi0W_all[idx1],
                  piAW = A * pi1W_all[idx1] + (1 - A) * pi0W_all[idx1],
                  pi1W_all = pi1W_all, pi0W_all = pi0W_all,
                  fit = picv$fits, method = pi_method)
  }

  list(Qfit = Qfit, gfit = gfit, sfit = sfit, pifit = pifit, g1W_all = g1W_all)
}
