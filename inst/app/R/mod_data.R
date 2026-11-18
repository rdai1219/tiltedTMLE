# Data module: upload a file, map columns, validate. Everything stays in memory;
# nothing is written to disk by the app. The heavy lifting is done by the plain R
# function tiltedTMLE:::prepare_analysis_data() (R/app-data.R), so this file is
# only the user interface.
#
# Returns (from mod_data_server): a reactive holding the prepared data, or an
# "error" object with a message. Step 3 (settings) and Step 4 (run) will use it.

# ---- UI pieces ----

# Inputs shown in the sidebar
mod_data_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h5("1. Upload data"),
    shiny::fileInput(ns("file"), "CSV or Excel file", accept = c(".csv", ".xlsx", ".xls")),
    shiny::radioButtons(
      ns("setting"), "Method",
      choices = c("Fixed h (single study)" = "fixed",
                  "Internally tilted (single study)" = "internal",
                  "Externally tilted (trial + external covariates)" = "external")),
    shiny::uiOutput(ns("mapping")),
    shiny::checkboxInput(ns("add_indicators"),
                         "Fill missing covariates and add missingness indicators", TRUE),
    shiny::checkboxInput(ns("show_preview"), "Show a preview of the first rows (displays raw data)", FALSE),
    shiny::actionButton(ns("validate"), "Check data", class = "btn-primary")
  )
}

# Results shown in the main panel
mod_data_output_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("status")),
    shiny::tableOutput(ns("summary_table")),
    shiny::tableOutput(ns("covariate_table")),
    shiny::uiOutput(ns("preview_ui"))
  )
}

# ---- Server ----

mod_data_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {

    # Read the uploaded file (kept in memory only)
    raw <- shiny::reactive({
      shiny::req(input$file)
      path <- input$file$datapath
      ext <- tolower(tools::file_ext(input$file$name))
      if (ext == "csv") {
        as.data.frame(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
      } else if (ext %in% c("xlsx", "xls")) {
        as.data.frame(readxl::read_excel(path))
      } else {
        stop("Please upload a .csv, .xlsx or .xls file.")
      }
    })

    # Column mapping controls, built from the file's column names
    output$mapping <- shiny::renderUI({
      ns <- session$ns
      cols <- names(raw())
      # keep earlier choices when the controls are rebuilt (for example after changing the method)
      keep <- function(id) {
        v <- shiny::isolate(input[[id]])
        if (!is.null(v) && length(v) > 0 && all(v %in% cols)) v else NULL
      }
      shiny::tagList(
        shiny::h5("2. Map columns"),
        shiny::selectInput(ns("y"), "Outcome (binary)", cols, selected = keep("y")),
        shiny::uiOutput(ns("y_one_ui")),
        shiny::selectInput(ns("a"), "Treatment (binary)", cols, selected = keep("a")),
        shiny::uiOutput(ns("a_one_ui")),
        # method-specific controls stay in place and are only shown when they apply, so a change of
        # method never rebuilds (and never resets) the other choices
        shiny::conditionalPanel(
          sprintf("input['%s'] == 'external'", ns("setting")),
          shiny::selectInput(ns("d"), "Internal vs external indicator (binary)", cols, selected = keep("d")),
          shiny::uiOutput(ns("d_one_ui"))),
        shiny::selectizeInput(ns("covariates"), "Covariates", cols, selected = keep("covariates"), multiple = TRUE),
        shiny::conditionalPanel(
          sprintf("input['%s'] == 'fixed'", ns("setting")),
          shiny::selectInput(ns("h"), "Tilting weights h(x) (required)", c("Choose a column..." = "", cols), selected = keep("h")),
          shiny::helpText("One non-negative number per row: the weight h(x) that defines your target population.",
                          "For example, 1 for everyone gives the plain ATE; 1 for subjects aged 65+ and 0 for the rest",
                          "targets that subgroup. Rows with a larger weight count more."))
      )
    })

    # Remember the "counts as 1" choices: those inputs disappear while the controls are rebuilt
    # (for example after changing the method), so keep the last non-empty value.
    remembered <- shiny::reactiveValues()
    for (nm in c("y_one", "a_one", "d_one")) {
      local({
        id <- nm
        shiny::observe({ v <- input[[id]]; if (!is.null(v)) remembered[[id]] <- v })
      })
    }

    # "Which value counts as 1?" pickers, shown once a column is chosen
    one_picker <- function(col_input, label, id_name) {
      shiny::renderUI({
        shiny::req(input[[col_input]], raw())
        v <- utils::head(binary_levels(raw()[[input[[col_input]]]]), 20)
        old <- shiny::isolate(remembered[[id_name]])
        shiny::selectInput(session$ns(id_name), label, v, selected = if (!is.null(old) && old %in% v) old else NULL)
      })
    }
    output$y_one_ui <- one_picker("y", "Outcome value that counts as 1 (event)", "y_one")
    output$a_one_ui <- one_picker("a", "Treatment value that counts as 1 (treated)", "a_one")
    output$d_one_ui <- one_picker("d", "Value that marks INTERNAL (trial) rows", "d_one")

    # The choices that define the analysis data (a change in any of them needs a new check)
    current_choices <- function() {
      list(setting = input$setting, y = input$y, y_one = input$y_one, a = input$a, a_one = input$a_one,
           covariates = input$covariates,
           d = if (identical(input$setting, "external")) input$d else NULL,
           d_one = if (identical(input$setting, "external")) input$d_one else NULL,
           h = if (identical(input$setting, "fixed")) input$h else NULL,
           add_indicators = input$add_indicators)
    }
    checked_choices <- shiny::reactiveVal(NULL)

    # Run the validation when the user clicks "Check data"
    validated <- shiny::eventReactive(input$validate, {
      checked_choices(current_choices())
      tryCatch(
        prepare_analysis_data(
          raw(), setting = input$setting,
          y = input$y, a = input$a, covariates = input$covariates,
          y_one = input$y_one, a_one = input$a_one,
          d = if (input$setting == "external") input$d else NULL,
          d_one = if (input$setting == "external") input$d_one else NULL,
          add_indicators = input$add_indicators,
          h = if (input$setting == "fixed" && !is.null(input$h) && nzchar(input$h)) input$h else NULL),
        error = function(e) structure(list(message = conditionMessage(e)), class = "prep_error"))
    })

    # What the rest of the app sees: the checked data, unless the user changed the method or the
    # column choices since the last check, in which case it asks for a new check.
    prepared <- shiny::reactive({
      p <- validated()
      if (!identical(checked_choices(), current_choices())) {
        return(structure(list(message = "The method or column choices changed since the last check. Click \"Check data\" again."),
                         class = "prep_error"))
      }
      p
    })

    output$status <- shiny::renderUI({
      p <- prepared()
      if (inherits(p, "prep_error")) {
        shiny::div(class = "alert alert-danger", shiny::strong("Please fix: "), p$message)
      } else {
        shiny::div(class = "alert alert-success",
                   shiny::strong("Data look usable."),
                   if (length(p$notes) > 0) shiny::tags$ul(lapply(p$notes, shiny::tags$li)))
      }
    })

    output$summary_table <- shiny::renderTable({
      p <- prepared(); shiny::req(!inherits(p, "prep_error")); s <- p$summary
      tab <- data.frame(
        Item = c("Rows in file", "Rows analysed", if (s$setting == "external") "External rows",
                 "Treated", "Control", "Outcome missing (treated / control)", "Outcome mean (observed)"),
        Value = c(s$n_rows, s$n_analysis, if (s$setting == "external") s$n_external,
                  s$n_treated, s$n_control,
                  paste0(s$outcome_missing_by_arm[["treated"]], " / ", s$outcome_missing_by_arm[["control"]]),
                  format(round(s$outcome_mean_observed, 3), nsmall = 3)),
        stringsAsFactors = FALSE)
      if (!is.null(s$h)) {   # fixed method: show what the weights look like
        tab <- rbind(tab, data.frame(
          Item = c("h(x): min / mean / max", "h(x): effective sample size"),
          Value = c(paste(format(round(c(s$h$min, s$h$mean, s$h$max), 3), nsmall = 3), collapse = " / "),
                    sprintf("%.1f of %d", s$h$ess, s$n_analysis)),
          stringsAsFactors = FALSE))
      }
      tab
    })

    output$covariate_table <- shiny::renderTable({
      p <- prepared(); shiny::req(!inherits(p, "prep_error")); s <- p$summary
      miss <- as.integer(s$covariate_missing[s$covariates])
      data.frame(Covariate = s$covariates,
                 `Missing before filling` = ifelse(is.na(miss), 0L, miss),
                 check.names = FALSE)
    })

    # Optional raw preview, off by default (privacy)
    output$preview_ui <- shiny::renderUI({
      shiny::req(input$show_preview)
      shiny::tagList(shiny::h5("Preview (first 10 rows)"), shiny::tableOutput(session$ns("preview")))
    })
    output$preview <- shiny::renderTable(utils::head(raw(), 10))

    prepared
  })
}
