# Settings module: choose the estimand, nuisance learners, truncation and inference
# options. It does no estimation; it turns the choices into an "analysis plan" with
# tiltedTMLE:::build_analysis_spec() (R/app-spec.R), which Step 4 will run.
#
# mod_settings_server(id, prepared) returns a reactive holding the plan, or an error
# object (class "spec_error") with a readable message. `prepared` is the reactive
# returned by mod_data_server().

# Learner picker for one nuisance model: glm, or SuperLearner with checkboxes
learner_ui <- function(ns, prefix, label) {
  shiny::tagList(
    shiny::radioButtons(ns(paste0(prefix, "_method")), label, c("glm" = "glm", "SuperLearner" = "SL"),
                        inline = TRUE),
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'SL'", ns(paste0(prefix, "_method"))),
      shiny::checkboxGroupInput(ns(paste0(prefix, "_library")), NULL,
                                choices = names(SL_PRESETS), selected = c("SL.glm", "SL.glmnet")))
  )
}

mod_settings_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("estimand_ui")),
    shiny::uiOutput(ns("custom_ui")),
    shiny::h5("Nuisance models"),
    learner_ui(ns, "Q", "Outcome model"),
    learner_ui(ns, "g", "Treatment model"),
    shiny::uiOutput(ns("s_ui")),
    shiny::uiOutput(ns("pi_ui")),
    shiny::h5("Inference"),
    shiny::checkboxInput(ns("bootstrap"), "Bootstrap confidence intervals", FALSE),
    shiny::conditionalPanel(
      sprintf("input['%s']", ns("bootstrap")),
      shiny::numericInput(ns("n_boot"), "Number of resamples", 200, min = 20, max = 5000, step = 50),
      shiny::numericInput(ns("seed"), "Random seed (blank = random)", NA)),
    shiny::numericInput(ns("conf_level"), "Confidence level", 0.95, min = 0.5, max = 0.999, step = 0.01),
    shiny::tags$details(
      shiny::tags$summary("Advanced"),
      shiny::numericInput(ns("g_lower"), "Treatment propensity: lower truncation bound (0 = none)", 0.025, min = 0, max = 0.49, step = 0.005),
      shiny::uiOutput(ns("s_lower_ui")),
      shiny::uiOutput(ns("pi_lower_ui")),
      shiny::numericInput(ns("max_iter"), "Maximum targeting iterations (tilted methods)", 50, min = 1),
      shiny::numericInput(ns("tol"), "Convergence tolerance (tilted methods)", 1e-6, min = 1e-12))
  )
}

mod_settings_output_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("plan"))
}

mod_settings_server <- function(id, prepared) {
  shiny::moduleServer(id, function(input, output, session) {
    ok_data <- shiny::reactive({
      p <- prepared(); shiny::req(!inherits(p, "prep_error")); p
    })

    # Estimand choices depend on the method chosen in the data step
    output$estimand_ui <- shiny::renderUI({
      p <- ok_data()
      shiny::tagList(
        shiny::h5("Estimand"),
        shiny::selectInput(session$ns("estimand"), NULL, estimand_choices(p$setting)))
    })
    # Custom tilting function: shown only when "Custom" is the chosen estimand
    output$custom_ui <- shiny::renderUI({
      p <- ok_data(); shiny::req(p$setting != "fixed", identical(input$estimand, "custom"))
      v <- if (p$setting == "internal") "g" else "s"
      what <- if (p$setting == "internal") "treatment propensity g(x)" else "selection propensity s(x)"
      shiny::tagList(
        shiny::textInput(session$ns("custom_eta"), paste0("eta(", v, ") ="), paste0(v, "*(1-", v, ")")),
        shiny::textInput(session$ns("custom_deta"), paste0("derivative of eta(", v, ") (blank = compute automatically)"), ""),
        shiny::helpText(paste0("h(x) = eta(", v, "), where ", v, " is the ", what, ". Use only the variable ", v,
                               ", numbers and + - * / ^ ( ) exp log sqrt abs pmin pmax ifelse plogis and comparisons. ",
                               "Examples: 1 (ATE), ", v, " (ATT), 1-", v, " (ATC), ", v, "*(1-", v, ") (ATO).")),
        # the wrapper stops the surrounding flex layout from squeezing the plot
        shiny::div(style = "flex: none; min-height: 250px;",
                   shiny::plotOutput(session$ns("eta_plot"), height = "240px", fill = FALSE)))
    })
    output$eta_plot <- shiny::renderPlot({
      p <- ok_data(); shiny::req(identical(input$estimand, "custom"))
      v <- if (p$setting == "internal") "g" else "s"
      plot_tilt(make_tilt(input$custom_eta, input$custom_deta, v))
    })

    # Selection model and its bound only exist for the external method
    output$s_ui <- shiny::renderUI({
      shiny::req(ok_data()$setting == "external"); learner_ui(session$ns, "s", "Selection model")
    })
    output$s_lower_ui <- shiny::renderUI({
      shiny::req(ok_data()$setting == "external")
      shiny::numericInput(session$ns("s_lower"), "Selection propensity: lower truncation bound (0 = none)", 0.025, min = 0, max = 0.49, step = 0.005)
    })
    # Missing-outcome model only matters when some outcomes are missing
    output$pi_ui <- shiny::renderUI({
      shiny::req(ok_data()$summary$n_outcome_missing > 0)
      learner_ui(session$ns, "pi", "Missing-outcome model")
    })
    output$pi_lower_ui <- shiny::renderUI({
      shiny::req(ok_data()$summary$n_outcome_missing > 0)
      shiny::numericInput(session$ns("pi_lower"), "Missing-outcome propensity: lower bound", 0.025, min = 0, max = 0.49, step = 0.005)
    })

    # Build the plan whenever data or settings change
    spec <- shiny::reactive({
      p <- ok_data()
      st <- list(estimand = input$estimand, custom_eta = input$custom_eta, custom_deta = input$custom_deta,
                 Q_method = input$Q_method, Q_library = input$Q_library,
                 g_method = input$g_method, g_library = input$g_library,
                 s_method = input$s_method, s_library = input$s_library,
                 pi_method = input$pi_method, pi_library = input$pi_library,
                 g_lower = input$g_lower, s_lower = input$s_lower, pi_lower = input$pi_lower,
                 bootstrap = isTRUE(input$bootstrap), n_boot = input$n_boot, seed = input$seed,
                 conf_level = input$conf_level, max_iter = input$max_iter, tol = input$tol)
      st <- Filter(Negate(is.null), st)           # inputs not shown yet fall back to defaults
      tryCatch(build_analysis_spec(p, st),
               error = function(e) structure(list(message = conditionMessage(e)), class = "spec_error"))
    })

    output$plan <- shiny::renderUI({
      p <- tryCatch(prepared(), error = function(e) NULL)
      if (is.null(p)) return(shiny::div(class = "alert alert-info", "Upload your data and click \"Check data\" on the Data tab first."))
      if (inherits(p, "prep_error")) return(shiny::div(class = "alert alert-warning", shiny::strong("Data step: "), p$message))
      s <- spec()
      if (inherits(s, "spec_error")) {
        return(shiny::div(class = "alert alert-danger", shiny::strong("Please fix: "), s$message))
      }
      shiny::tagList(
        shiny::h5("Analysis plan"),
        shiny::tags$ul(lapply(s$description, shiny::tags$li)),
        if (length(s$warnings) > 0)
          shiny::div(class = "alert alert-warning", shiny::tags$ul(lapply(s$warnings, shiny::tags$li))))
    })

    spec
  })
}
