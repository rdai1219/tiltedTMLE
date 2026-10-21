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
regression. There are three entry points, one per setting:

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

## Documentation

- `vignette("tiltedTMLE")` -- a short getting-started walkthrough of all
  three entry points.
- `help(package = "tiltedTMLE")` -- full function reference.
- `TMLE_tilted_algorithms.pdf` (in this repo) -- the underlying derivations,
  efficient influence functions, and algorithm pseudocode.

## License

MIT (c) Yixin Fang, Weili He, Ran Dai
