#' Internally tilted ATE presets
#'
#' Returns `eta(g)` and its derivative `eta_dot(g)` for the four
#' internally tilted estimands of Section 2 of the reference paper, where
#' `g` is the treatment propensity score.
#'
#' @param estimand One of `"ATE"`, `"ATT"`, `"ATC"`, `"ATO"`.
#'
#' @return A list with components:
#' \item{eta}{Function of `g`, the tilting weight `eta(g)`.}
#' \item{eta_dot}{Function of `g`, the derivative `d/dg eta(g)`.}
#'
#' @seealso [estimate_tau_h_internal_tilted()], the entry point this is
#'   normally passed to.
#' @examples
#' spec <- tilt_spec("ATT")
#' g <- c(0.2, 0.5, 0.8)
#' spec$eta(g)
#' spec$eta_dot(g)
#' @export
tilt_spec <- function(estimand = c("ATE", "ATT", "ATC", "ATO")) {
  estimand <- match.arg(estimand)
  switch(estimand,
    ATE = list(eta = function(g) rep(1, length(g)),
               eta_dot = function(g) rep(0, length(g))),
    ATT = list(eta = function(g) g,
               eta_dot = function(g) rep(1, length(g))),
    ATC = list(eta = function(g) 1 - g,
               eta_dot = function(g) rep(-1, length(g))),
    ATO = list(eta = function(g) g * (1 - g),
               eta_dot = function(g) 1 - 2 * g)
  )
}

#' Externally tilted ATE presets
#'
#' Returns `eta(s)` and its derivative `eta_dot(s)` for the three
#' externally tilted estimands of Section 3 of the reference paper, where
#' `s` is the selection propensity score between an internal trial
#' population and an external covariate-only population.
#'
#' @param estimand One of `"ATE-I"` (internal population), `"ATE-E"`
#'   (external population), `"ATE-O"` (overlap population).
#'
#' @return A list with components:
#' \item{eta}{Function of `s`, the tilting weight `eta(s)`.}
#' \item{eta_dot}{Function of `s`, the derivative `d/ds eta(s)`.}
#'
#' @seealso [estimate_tau_h_external_tilted()], the entry point this is
#'   normally passed to.
#' @examples
#' spec <- tilt_spec_external("ATE-O")
#' s <- c(0.2, 0.5, 0.8)
#' spec$eta(s)
#' spec$eta_dot(s)
#' @export
tilt_spec_external <- function(estimand = c("ATE-I", "ATE-E", "ATE-O")) {
  estimand <- match.arg(estimand)
  switch(estimand,
    `ATE-I` = list(eta = function(s) rep(1, length(s)),
                   eta_dot = function(s) rep(0, length(s))),
    `ATE-E` = list(eta = function(s) (1 - s) / s,
                   eta_dot = function(s) -1 / s^2),
    `ATE-O` = list(eta = function(s) 1 - s,
                   eta_dot = function(s) rep(-1, length(s)))
  )
}
