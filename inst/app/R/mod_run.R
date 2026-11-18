# Run module: a "Run analysis" button that starts the plan in a background R process.
#   stage 1 (main fit, seconds)  -> results appear straight away
#   stage 2 (bootstrap, if on)   -> added when it finishes
# The calculations live in R/app-run.R; this file is only the user interface.
# Everything stays in memory. One run at a time per session.
#
# mod_run_server(id, prepared, spec) returns a list of reactives for Step 5:
#   main()  the stage-1 result bundle (or a "run_error"),
#   boot()  the stage-2 bootstrap result (or a "run_error"; NULL if not requested),
#   inputs() the frozen data + plan the run used.

# Button and status (shown next to the plan)
mod_run_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    bslib::input_task_button(ns("run"), "Run analysis"),
    shiny::uiOutput(ns("status"))
  )
}

mod_run_server <- function(id, prepared, spec) {
  shiny::moduleServer(id, function(input, output, session) {
    # Same package version in the worker as in this session (see run_in_background())
    pkg_path <- getNamespaceInfo(asNamespace("tiltedTMLE"), "path")

    main_task <- shiny::ExtendedTask$new(function(prepared, spec, pkg_path) {
      promises::as.promise(run_in_background("main", prepared, spec, pkg_path))
    })
    boot_task <- shiny::ExtendedTask$new(function(prepared, spec, pkg_path) {
      promises::as.promise(run_in_background("bootstrap", prepared, spec, pkg_path))
    })
    bslib::bind_task_button(main_task, "run")

    inputs <- shiny::reactiveVal(NULL)       # data and plan frozen at click time
    started <- shiny::reactiveVal(NULL)
    boot_wanted <- shiny::reactive({
      i <- inputs(); !is.null(i) && i$spec$bootstrap$on && i$spec$method != "fixed"
    })

    shiny::observeEvent(input$run, {
      p <- prepared(); s <- spec()
      if (inherits(p, "prep_error") || inherits(s, "spec_error")) {
        shiny::showNotification("Fix the data or settings first.", type = "error"); return()
      }
      if (identical(boot_task$status(), "running")) {
        shiny::showNotification("The previous bootstrap is still running.", type = "warning"); return()
      }
      inputs(list(prepared = p, spec = s))   # freeze: later edits do not affect this run
      started(Sys.time())
      main_task$invoke(p, s, pkg_path)
    })

    # Start the bootstrap once the main fit has succeeded
    shiny::observeEvent(main_task$result(), {
      r <- main_task$result()
      if (boot_wanted() && !inherits(r, "run_error")) {
        i <- inputs(); boot_task$invoke(i$prepared, i$spec, pkg_path)
      }
    })

    # Results as reactives (errors from the background process become run_error objects)
    as_error <- function(e) structure(list(message = conditionMessage(e)), class = "run_error")
    main <- shiny::reactive(tryCatch(main_task$result(), error = as_error))
    boot <- shiny::reactive({
      shiny::req(boot_wanted())
      tryCatch(boot_task$result(), error = as_error)
    })

    # Status line with elapsed time while running
    output$status <- shiny::renderUI({
      m <- main_task$status(); b <- boot_task$status()
      running <- identical(m, "running") || (boot_wanted() && identical(b, "running")) ||
        (boot_wanted() && identical(m, "success") && identical(b, "initial"))
      if (running) shiny::invalidateLater(1000)
      elapsed <- if (!is.null(started())) round(as.numeric(difftime(Sys.time(), started(), units = "secs"))) else 0
      txt <- if (identical(m, "running")) paste0("Running the main fit... ", elapsed, " s")
             else if (running) paste0("Main fit done; running the bootstrap... ", elapsed, " s")
             else if (identical(m, "error") || identical(b, "error")) "The run stopped with an error (see Results)."
             else if (identical(m, "success")) "Finished."
             else "Ready."
      shiny::p(shiny::em(txt))
    })

    list(main = main, boot = boot, inputs = inputs)
  })
}
