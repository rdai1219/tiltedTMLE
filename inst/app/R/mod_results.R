# Results module: shows a finished run in three sub-tabs.
#   Estimates            plain-language summary, estimates table, estimator comparison plot, CSV
#   Diagnostics          overlap and weight plots, diagnostics table, PNG downloads
#   Reproduce & report   R script that reruns the analysis, HTML report
# All numbers, plots, text and files come from plain functions in R/app-results.R; this file
# is only the user interface. Only summaries are shown or downloaded, never raw rows.
#
# mod_results_server(id, run) takes the list returned by mod_run_server().

mod_results_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(shiny::uiOutput(ns("message")), shiny::uiOutput(ns("tabs")))
}

# The three sub-tabs (built by the server once a run has finished; ids need the module namespace)
results_tabs <- function(ns) {
  bslib::navset_card_tab(
    bslib::nav_panel("Estimates",
      shiny::uiOutput(ns("summary")),
      shiny::tableOutput(ns("table")),
      shiny::plotOutput(ns("plot_estimators"), height = "300px", fill = FALSE),
      shiny::uiOutput(ns("notes")),
      shiny::downloadButton(ns("dl_csv"), "Download estimates (CSV)")),
    bslib::nav_panel("Diagnostics",
      shiny::tableOutput(ns("diag_table")),
      shiny::plotOutput(ns("plot_overlap"), height = "380px", fill = FALSE),
      shiny::downloadButton(ns("dl_overlap"), "Overlap plot (PNG)"),
      shiny::plotOutput(ns("plot_weights"), height = "300px", fill = FALSE),
      shiny::downloadButton(ns("dl_weights"), "Weights plot (PNG)"),
      shiny::downloadButton(ns("dl_estimators"), "Estimator plot (PNG)")),
    bslib::nav_panel("Reproduce & report",
      shiny::p("This R script reruns the analysis with the tiltedTMLE package on your own file."),
      shiny::verbatimTextOutput(ns("script")),
      shiny::downloadButton(ns("dl_script"), "Download R script"),
      shiny::downloadButton(ns("dl_report"), "Download HTML report")))
}

mod_results_server <- function(id, run) {
  shiny::moduleServer(id, function(input, output, session) {
    # A finished main fit (NULL until a run has produced one)
    main <- shiny::reactive({
      if (is.null(run$inputs())) return(NULL)
      m <- tryCatch(run$main(), error = function(e) NULL)
      if (is.null(m)) NULL else m
    })
    ok_main <- shiny::reactive({ m <- main(); shiny::req(m, !inherits(m, "run_error")); m })
    # The bootstrap result, only if it exists and succeeded
    boot_ok <- shiny::reactive({
      b <- tryCatch(run$boot(), error = function(e) NULL)
      if (is.null(b) || inherits(b, "run_error")) NULL else b
    })

    # The tabs exist only after a successful main fit
    output$tabs <- shiny::renderUI({ ok_main(); results_tabs(session$ns) })

    output$message <- shiny::renderUI({
      m <- main()
      if (is.null(m)) return(shiny::p("Run an analysis on the Settings tab to see results here."))
      if (inherits(m, "run_error")) {
        shiny::div(class = "alert alert-danger", shiny::strong("The analysis failed: "), m$message,
                   shiny::tags$details(shiny::tags$summary("Technical details (no data included)"),
                                       shiny::tags$p("Failing call: ", if (is.null(m$call) || is.na(m$call)) "unknown" else m$call),
                                       shiny::tags$p("Call chain: ", paste(m$where, collapse = " > "))))
      }
    })

    # ---- Estimates tab ----
    output$summary <- shiny::renderUI({
      s <- summary_text(ok_main(), boot_ok())
      shiny::tagList(shiny::p(shiny::strong(s$sentence)),
                     shiny::tags$small(shiny::tags$ul(lapply(s$caveats, shiny::tags$li))))
    })
    output$table <- shiny::renderTable(results_table(ok_main(), boot_ok()), digits = 4, na = "")
    output$plot_estimators <- shiny::renderPlot(plot_estimators(ok_main(), boot_ok()))
    output$notes <- shiny::renderUI({
      m <- ok_main(); b <- boot_ok(); be <- tryCatch(run$boot(), error = function(e) NULL)
      items <- c(
        paste0("Main fit took ", round(m$seconds, 1), " s."),
        m$warnings,
        if (inherits(be, "run_error")) paste0("Bootstrap failed: ", be$message),
        if (!is.null(b)) paste0("Bootstrap: ", b$n_failed, " failed and ", b$n_nonconverged,
                                " non-converged resample(s) of ", length(b$estimates), "; took ", round(b$seconds), " s."))
      shiny::tags$ul(lapply(items, shiny::tags$li))
    })

    # ---- Diagnostics tab ----
    output$diag_table <- shiny::renderTable(diagnostics_table(ok_main()))
    output$plot_overlap <- shiny::renderPlot(plot_overlap(ok_main()))
    output$plot_weights <- shiny::renderPlot(plot_weights(ok_main()))

    # ---- Reproduce & report tab ----
    script_text <- shiny::reactive({
      shiny::req(run$inputs(), ok_main())
      i <- run$inputs(); reproduction_script(i$prepared, i$spec)
    })
    output$script <- shiny::renderText(script_text())

    # ---- downloads (summaries only) ----
    png_download <- function(plot_fun, file_name) {
      shiny::downloadHandler(
        filename = function() file_name,
        content = function(file) ggplot2::ggsave(file, plot_fun(), width = 7, height = 5, dpi = 150))
    }
    output$dl_csv <- shiny::downloadHandler(
      filename = function() "tiltedTMLE_estimates.csv",
      content = function(file) utils::write.csv(results_table(ok_main(), boot_ok()), file, row.names = FALSE, na = ""))
    output$dl_overlap <- png_download(function() plot_overlap(ok_main()), "propensity_overlap.png")
    output$dl_weights <- png_download(function() plot_weights(ok_main()), "tilting_weights.png")
    output$dl_estimators <- png_download(function() plot_estimators(ok_main(), boot_ok()), "estimator_comparison.png")
    output$dl_script <- shiny::downloadHandler(
      filename = function() "reproduce_analysis.R",
      content = function(file) writeLines(script_text(), file))
    output$dl_report <- shiny::downloadHandler(
      filename = function() "tiltedTMLE_report.html",
      content = function(file) render_report(ok_main(), boot_ok(), script_text(), file))
  })
}
