#' tiltedTMLE: Targeted Maximum Likelihood Estimation for Tilted Average Treatment Effects
#'
#' Implements targeted maximum likelihood estimators (TMLE) for a broad
#' class of tilted average treatment effect (ATE) estimands, covering two
#' related but distinct tilting-function settings. Fang, He, and Dai,
#' "Targeted Learning for Externally Controlled Trials and Beyond" (AbbVie,
#' 2026), treats the tilting function `h(x)` as a known, fixed weight
#' prespecified by the user, estimated with a single targeting step. Fang,
#' "Targeted Learning for Tilted Average Treatment Effects" (AbbVie, 2026),
#' instead treats `h(x)` as a known function of an unknown nuisance
#' function that must itself be estimated and iteratively targeted
#' alongside the outcome regression.
#'
#' @section Fixed/precomputed tilting functions:
#' [estimate_tau_h_fixed()] / [tmle_h()]: `h(x)` is either a fixed, precomputed
#' vector (e.g. `h(x) = 1` for the ATE), or a function evaluated once at
#' the fitted propensity score and never re-targeted (e.g. `h(x) =` the
#' estimated propensity score itself for the ATT).
#'
#' @section Internally tilted estimands:
#' [estimate_tau_h_internal_tilted()] / [tmle_tilted_h()]: ATE, ATT, ATC, ATO,
#' obtained by tilting the ATE by `h(x) = eta(g(x))`, a known function of
#' the treatment propensity score `g(x)`, within a single study population
#' (see [tilt_spec()]). `g(x)` is iteratively targeted alongside the
#' outcome regression.
#'
#' @section Externally tilted estimands:
#' [estimate_tau_h_external_tilted()] / [tmle_external_tilted_h()]:
#' ATE-I, ATE-E, ATE-O, obtained by tilting by `h(x) = eta(s(x))`, a known
#' function of the selection propensity score `s(x)` between an internal
#' trial population and an external covariate-only population (see
#' [tilt_spec_external()]). `s(x)` is iteratively targeted alongside the
#' outcome regression.
#'
#' @keywords internal
"_PACKAGE"
