# tiltedTMLE Shiny app. Launch with tiltedTMLE::run_app().
# Shiny automatically loads every file in the R/ folder next to this file.

if (!"tiltedTMLE" %in% loadedNamespaces()) library(tiltedTMLE)

# The module code calls the package's internal helpers, so make them visible here.
for (nm in c("prepare_analysis_data", "binary_levels", "build_analysis_spec", "estimand_choices",
             "SL_PRESETS", "run_in_background", "est_row",
             "results_table", "diagnostics_table", "plot_overlap", "plot_weights", "plot_estimators",
             "summary_text", "reproduction_script", "render_report", "make_tilt", "plot_tilt")) {
  assign(nm, get(nm, envir = asNamespace("tiltedTMLE")))
}

options(shiny.maxRequestSize = 50 * 1024^2)   # allow uploads up to 50 MB

# Analyses run in separate background R processes (at most 2 at a time on this machine),
# so the app stays responsive during long SuperLearner or bootstrap runs.
future::plan(future::multisession, workers = 2)

ui <- bslib::page_sidebar(
  title = "tiltedTMLE: tilted treatment effects with TMLE",
  sidebar = bslib::sidebar(width = 380, mod_data_ui("data")),
  bslib::navset_card_tab(
    id = "tabs",
    bslib::nav_panel("1. Data", mod_data_output_ui("data")),
    bslib::nav_panel("2. Settings",
                     bslib::layout_columns(col_widths = c(6, 6),
                                           mod_settings_ui("settings"),
                                           shiny::tagList(mod_settings_output_ui("settings"),
                                                          shiny::hr(),
                                                          mod_run_ui("run")))),
    bslib::nav_panel("3. Results", mod_results_ui("results"))
  )
)

server <- function(input, output, session) {
  prepared <- mod_data_server("data")
  spec <- mod_settings_server("settings", prepared)
  run <- mod_run_server("run", prepared, spec)
  mod_results_server("results", run)
}

shiny::shinyApp(ui, server)
