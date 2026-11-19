# tiltedTMLE 0.3.0

## New features

* Cross-fitting. New argument `cv_folds` in `estimate_tau_h_fixed()`,
  `estimate_tau_h_internal_tilted()` and `estimate_tau_h_external_tilted()` (and, through
  `...`, the full-refit bootstrap functions). With `cv_folds = V > 1` the rows are split into
  `V` folds (stratified by treatment arm and, in the external setting, by internal/external
  membership); every nuisance model (`Q`, `g`, `s`, `pi`) is fit on `V - 1` folds and predicts
  the held-out fold, and the TMLE and its EIF-based inference use these out-of-fold
  predictions. This keeps the standard errors valid with flexible SuperLearner libraries.
  The default `cv_folds = 1` gives the same results as before. Results gain `folds`
  (the fold of each row, or `NULL`), and with cross-fitting the `fit` element of each
  fitted-nuisance object is the list of per-fold models.

# tiltedTMLE 0.2.0

## New features

* Continuous outcomes. `Y` may now be continuous in `estimate_tau_h_fixed()`,
  `estimate_tau_h_internal_tilted()`, `estimate_tau_h_external_tilted()` and the two
  full-refit bootstrap functions. Following Fang, He and Dai (2026, Sec. 3.3.2), `Y` is
  rescaled to `[0, 1]` with the new argument `Ybounds` (default: `c(0, 1)` for a binary
  outcome, the observed range otherwise), the TMLE is computed on that scale, and the
  estimates, standard errors and confidence intervals are mapped back to the original
  scale. The results gain a `Ybounds` element. Works together with missing outcomes,
  SuperLearner and the bootstrap (which keeps the same bounds on every resample).

## Changes

* `fit_Q()` fits a quasi-binomial glm (or a gaussian SuperLearner) when `Y` is not 0/1,
  and stops if `Y` is outside `[0, 1]`.
* The targeting step `fit_eps_h()` uses a quasi-binomial working model. For 0/1 data this
  gives exactly the same estimate as before.
* Input checks now accept any finite numeric `Y`; a non-numeric or infinite `Y` is an error.
* The Shiny app still handles binary outcomes only.

# tiltedTMLE 0.1.0

## New features

* `bootstrap_tau_h_internal_tilted()` and `bootstrap_tau_h_external_tilted()`:
  full-refit nonparametric bootstrap (nuisance models and iterative targeting
  are re-run on every resample), returning a bootstrap SE and normal /
  percentile confidence intervals.
* `add_missing_indicators()`: indicator method for missing baseline covariates
  (mean imputation for continuous, most common value for binary, plus a
  `<name>_miss` 0/1 column).

* Missing outcomes (missing at random given `A` and `W`): `Y` may contain `NA` in
  all three entry points (`estimate_tau_h_fixed()`, `estimate_tau_h_internal_tilted()`,
  `estimate_tau_h_external_tilted()`). `Q` is fit on the
  observed rows, the missingness propensity `P(Delta = 1 | A, W)` is fit with the new
  `fit_pi()`, and `g` is replaced by `g * pi` in the clever covariate and the EIF.
  In the external setting `pi` is fit on the internal rows and extrapolated to the external rows.

* `run_app()`: an interactive Shiny app (inst/app). Upload a CSV or Excel file, map the columns,
  choose the method (fixed h, internally tilted, externally tilted), estimand, nuisance learners
  (glm or SuperLearner), truncation and bootstrap options, then run the analysis in a background
  R process. Results show the estimates (TMLE and benchmarks), propensity-overlap and tilting-weight
  diagnostics, an estimator comparison plot and a plain-language summary, with downloads (CSV, PNG,
  an R script that reproduces the analysis, and an HTML report). Only summaries are shown or
  downloaded, and nothing is written to disk by the app. Needs the suggested packages shiny, bslib,
  readr, readxl, future, promises, ggplot2 (and rmarkdown for the report).

* App: the internally and externally tilted methods accept a **custom tilting function**. Type eta(g)
  (internal) or eta(s) (external), for example `g*(1-g)`; the derivative is computed automatically or
  can be typed. Formulas are checked and evaluated in a restricted environment (only the one variable,
  numbers and a short list of arithmetic and math functions), shown in a preview plot, and written into
  the downloadable reproduction script.
* App: for the fixed-h method the tilting weights h(x) are now a required column with help text, and the
  data check summarises them (min / mean / max and effective sample size).

## Improvements

* `fit_Q()` and `predict_Qhat()` gain `Qbounds` (default `c(1e-4, 1 - 1e-4)`);
  `Q` is also re-bounded after every targeting update, so `qlogis(Q)` stays
  finite.
* `tmle_tilted_h()` gains `gbounds` and re-clamps `g` after every targeting
  update, mirroring `sbounds` in the external version.
  `estimate_tau_h_internal_tilted()` passes its `gbounds` down.
* The internal and external TMLE loops now share one stopping rule and return
  `converged`; the `*_given_fits` functions return `TMLE_converged`. A warning
  is issued if `max_iter` is reached without convergence.
* The raw-data entry points validate their inputs (0/1 `Y`/`A`/`D`, matching
  lengths, no `NA`, sensible bounds) and give specific error messages.
* `estimate_tau_h_fixed()` warns when `h` is a function of the fitted
  propensity score, because its standard error then ignores the `eta_dot` term
  (new argument `warn_fn_h` to silence it).

## Breaking changes

* The shared `...` argument of `estimate_tau_h_fixed()`,
  `estimate_tau_h_internal_tilted()` and `estimate_tau_h_external_tilted()` is
  replaced by `Q_args`, `g_args` (and `s_args` for the external version), so
  extra arguments reach only the intended fit.
* `NA` in `Y`, `A`, `D`, `W`, `Ws` or `Wg` now raises an error instead of
  silently dropping rows; see `add_missing_indicators()` for covariates.

## Bug fixes

* App: after changing the method or column choices, the data must be checked again; the Settings tab
  now says so instead of silently keeping the previous method's options. Column choices (including the
  "counts as 1" values) are kept when the method changes.
* App: a covariate named A, Y or Delta clashed with the models' internal variable names; such columns are
  now renamed with a _cov suffix.
* SuperLearner (`method = "SL"`) failed with "object 'All' not found" unless the user had run
  `library(SuperLearner)`: SuperLearner looks up its learners and the default screening function by name
  from the calling environment, which cannot see them inside the package. All SuperLearner calls (outcome,
  treatment/selection and missing-outcome models) now go through one helper that points SuperLearner at
  its own namespace. A regression test runs every learner offered by the app with SuperLearner not attached.
* `fit_Q()` and `fit_g()` documented that `...` is passed to `glm()` but did
  not pass it; it now does (with `method = "glm"`).

# tiltedTMLE 0.0.1

* Initial version: fixed-h, internally tilted and externally tilted TMLE.
