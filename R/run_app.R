#' Launch the tiltedTMLE Shiny app
#'
#' Starts an interactive app for running the package's estimators on an uploaded
#' data file: upload a CSV or Excel file, map the columns (outcome, treatment,
#' covariates, and for the external setting the internal/external indicator), and
#' check the data. The app keeps data in memory only; nothing is written to disk by
#' the app itself.
#'
#' Requires the suggested packages `shiny`, `bslib`, `readr` and `readxl`.
#'
#' @param launch.browser Open the app in a browser? Defaults to `TRUE` in an
#'   interactive session.
#' @param ... Passed to [shiny::runApp()] (for example `port = 4321`).
#'
#' @return Called for its side effect (starts the app); returns when it is stopped.
#' @examples
#' if (interactive()) run_app()
#' @export
run_app <- function(launch.browser = interactive(), ...) {
  for (pkg in c("shiny", "bslib", "readr", "readxl")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for the app. Install it with install.packages(\"",
           pkg, "\").", call. = FALSE)
    }
  }
  app_dir <- system.file("app", package = "tiltedTMLE")
  if (app_dir == "") stop("Could not find the app folder; try reinstalling tiltedTMLE.", call. = FALSE)
  shiny::runApp(app_dir, launch.browser = launch.browser, ...)
}
