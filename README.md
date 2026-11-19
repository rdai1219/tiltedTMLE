# tiltedTMLE

Targeted maximum likelihood estimation (TMLE) for tilted average treatment
effect (ATE) estimands, including the plain ATE/ATT/ATC/ATO, and ATE-I/ATE-E/ATE-O
estimands for externally controlled trials.

## Installation

```r
# install.packages("remotes")  # if not already installed
remotes::install_github("rdai1219/tiltedTMLE")
```

## Overview

`tiltedTMLE` estimates estimands of the form

```
tau_h = E[h(X) * (Q(X,1) - Q(X,0))] / E[h(X)]
```

for a tilting function `h(x)` that is either fixed and known, or itself a
function of an estimated nuisance that must be targeted alongside the outcome
regression. The outcome `Y` can be binary (0/1) or continuous. There are three
entry points, one per setting:

- `estimate_tau_h_fixed()` -- `h(x)` is a fixed, precomputed weight (e.g.
  `h(x) = 1` for the plain ATE), estimated with a single TMLE targeting step.
- `estimate_tau_h_internal_tilted()` -- `h(x) = eta(g(x))`, a known function
  of the treatment propensity score, within one study population (ATE / ATT
  / ATC / ATO). `g` is iteratively targeted.
- `estimate_tau_h_external_tilted()` -- `h(x) = eta(s(x))`, a known function
  of the selection propensity score between an internal trial and an
  external covariate-only population (ATE-I / ATE-E / ATE-O). The selection
  score `s` is iteratively targeted.

## Quick start

```r
library(tiltedTMLE)

set.seed(42)
n <- 300  # internal trial
m <- 300  # external, covariate-only

W1 <- rnorm(n + m)
W2 <- rnorm(n + m)
D  <- c(rep(1L, n), rep(0L, m))  # 1 = internal trial, 0 = external
Ws <- data.frame(X1 = W1, X2 = W2)

A <- rbinom(n, 1, plogis(0.3 * W1[1:n] - 0.2 * W2[1:n]))
Y <- rbinom(n, 1, plogis(-0.3 + 0.5 * A + 0.3 * W1[1:n] - 0.2 * W2[1:n]))
W <- Ws[1:n, , drop = FALSE]  # internal-only covariates

# Fixed h: the plain ATE
out_ate <- estimate_tau_h_fixed(Y, A, W, h = rep(1, n))
out_ate$TMLE
out_ate$TMLE_CI

# Internally tilted: the ATT (h(x) = eta(g(x)) = g(x))
spec_att <- tilt_spec("ATT")
out_att <- estimate_tau_h_internal_tilted(Y, A, W, spec_att$eta, spec_att$eta_dot)
out_att$TMLE

# Externally tilted: the overlap population (ATE-O, h(x) = eta(s(x)) = 1 - s(x))
spec_o <- tilt_spec_external("ATE-O")
out_ateo <- estimate_tau_h_external_tilted(Y, A, D, Ws, spec_o$eta, spec_o$eta_dot)
out_ateo$TMLE
```

Each entry point returns a list including `TMLE` (point estimate),
`TMLE_se`/`TMLE_CI` (EIF-based Wald inference), un-targeted benchmark
estimators (`MLE`, `IPW`, and `AIPW` where defined), and the underlying
fitted-nuisance objects for further inspection. See `?estimate_tau_h_fixed`,
`?estimate_tau_h_internal_tilted`, and `?estimate_tau_h_external_tilted` for
the full set of arguments (SuperLearner nuisance fitting, propensity
truncation bounds, bootstrap options), or run `vignette("tiltedTMLE")` for a
longer walkthrough.

## Continuous outcomes

A continuous `Y` is passed in the same way. It is rescaled to `[0, 1]` using `Ybounds`
(by default the observed range), the TMLE is computed on that scale, and the estimate,
standard error and confidence interval are reported on the original scale.

```r
Y_cont <- 2 + 1.5 * A + W$X1 + rnorm(n)
out_c <- estimate_tau_h_internal_tilted(Y_cont, A, W, spec_att$eta, spec_att$eta_dot)
out_c$TMLE
out_c$Ybounds

# Known limits (e.g. a 0-100 rating scale) can be given instead of the observed range
# estimate_tau_h_internal_tilted(score, A, W, spec_att$eta, spec_att$eta_dot, Ybounds = c(0, 100))
```

## Cross-fitting with flexible learners

With flexible SuperLearner libraries, set `cv_folds` (e.g. 5): each nuisance model is fit
on the other folds and predicts the held-out fold, which keeps the standard errors and
confidence intervals valid. Fold assignment is random, so set a seed.

```r
set.seed(2026)
out_cv <- estimate_tau_h_internal_tilted(
  Y, A, W, spec_att$eta, spec_att$eta_dot,
  Q_method = "SL", Q_SL.library = c("SL.glm", "SL.ranger", "SL.xgboost"),
  g_method = "SL", g_SL.library = c("SL.glm", "SL.ranger", "SL.xgboost"),
  cv_folds = 5)
```

## Bootstrap inference for the tilted estimators

The EIF-based Wald interval can be too narrow at moderate sample sizes. For the
internally and externally tilted estimators there is a full-refit bootstrap:
every resample re-fits the nuisance models and re-runs the iterative TMLE.

```r
b <- bootstrap_tau_h_internal_tilted(Y, A, W, spec_att$eta, spec_att$eta_dot, n_boot = 200)
b$se             # bootstrap standard error
b$CI_percentile  # percentile interval (also b$CI_normal)

b2 <- bootstrap_tau_h_external_tilted(Y, A, D, Ws, spec_o$eta, spec_o$eta_dot, n_boot = 200)
```

## Missing baseline covariates

The estimators require complete data. Missing *outcomes* (`NA` in `Y`, missing at random given `A` and `W`) are supported in all three
entry points: the package fits `P(outcome observed | A, W)` and adjusts the TMLE automatically.

Missing *covariates* can be handled with
the indicator method (mean/majority-value imputation plus a `<name>_miss`
indicator column) before calling any entry point. Missing treatments are not supported.

```r
W_aug <- add_missing_indicators(W)   # for the external setting, apply to the pooled Ws
out <- estimate_tau_h_internal_tilted(Y, A, W_aug, spec_att$eta, spec_att$eta_dot)
```

## Shiny app

```r
tiltedTMLE::run_app()
```

For the internally and externally tilted methods you can pick ATE / ATT / ATC / ATO (ATE-I / ATE-E / ATE-O) or
**Custom** and type your own tilting function eta(g) (or eta(s)); for the fixed-h method you supply a column of
weights h(x).

Upload a data file (CSV or Excel), map the outcome, treatment and covariates (and, for the external
setting, the internal/external indicator), choose the estimand and models, and run. The app runs the
analysis in a background R process, shows the estimates with diagnostics and plots, and lets you
download the estimates (CSV), the plots (PNG), an R script that reproduces the analysis, and an HTML
report. It keeps data in memory only. Needs the suggested packages `shiny`, `bslib`, `readr`,
`readxl`, `future`, `promises` and `ggplot2` (plus `rmarkdown` for the report).

## Notes

- `estimate_tau_h_fixed()` accepts `h` as a function of the fitted propensity
  score (e.g. the ATT), but then treats `h` as known, so its standard errors are
  not valid for such `h` and a warning is issued. Use
  `estimate_tau_h_internal_tilted()` for the ATT/ATC/ATO.
- Extra arguments for the nuisance fits go in `Q_args`, `g_args` and `s_args`
  (named lists), e.g. `Q_args = list(cvControl = list(V = 5))` for SuperLearner.
- Each tilted fit reports `TMLE_converged`; a warning is given if the iteration
  hits `max_iter`.

## Documentation

- `vignette("tiltedTMLE")` -- a short getting-started walkthrough of all
  three entry points.
- `help(package = "tiltedTMLE")` -- full function reference.
- `TMLE_tilted_algorithms.pdf` (in this repo) -- the underlying derivations,
  efficient influence functions, and algorithm pseudocode.

## License

MIT (c) Yixin Fang, Weili He, Ran Dai
