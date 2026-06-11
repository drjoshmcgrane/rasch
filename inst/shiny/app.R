# RaschR Shiny GUI
# ---------------------------------------------------------------------------
# Run from a folder containing RaschR.R, RaschR_plots.R and this file:
#   shiny::runApp("app.R")
# Requires only the 'shiny' package; all analysis and plotting is base R.
# ---------------------------------------------------------------------------
library(shiny)
# Use the installed package if available, otherwise source the loose script files.
if (requireNamespace("RaschR", quietly = TRUE)) {
  library(RaschR)
} else {
  source("RaschR.R"); source("RaschR_plots.R")
}

# Built-in demo data (10 polytomous items, one deliberately disordered, plus a group).
.demo_data <- function(seed = 11, Np = 1200) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  mvec <- rep(c(2, 3), length.out = 10)
  tau_true <- lapply(mvec, function(m) sort(rnorm(m, 0, 0.9)))
  tau_true[[2]] <- c(1.2, -1.3, 0.6)                       # disordered item
  th <- rnorm(Np, 0, 1.4); grp <- rep(c("ref", "foc"), each = Np / 2)
  shift <- numeric(Np); 
  X <- sapply(seq_along(mvec), function(i) {
    sft <- if (i == 5) ifelse(grp == "foc", 0.9, 0) else 0  # uniform DIF planted on item 5
    sapply(seq_len(Np), function(n) sample(0:mvec[i], 1, prob = simP(th[n] - sft[n], tau_true[[i]])))
  })
  colnames(X) <- sprintf("Q%02d", seq_along(mvec))
  data.frame(X, group = grp, check.names = FALSE)
}

ui <- fluidPage(
  tags$head(tags$style(HTML(
    ".btn-primary{background:#1f3b5c;border-color:#1f3b5c} h2{color:#1f3b5c}"))),
  titlePanel("RaschR — pairwise Rasch analysis"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("file", "Response data (CSV)", accept = ".csv"),
      checkboxInput("demo", "Use built-in demo data", TRUE),
      selectInput("model", "Model", c("Partial credit (PCM)" = "PCM", "Rating scale (RSM)" = "RSM")),
      selectInput("solver", "Solver (PCM only)", c("Least squares" = "LS", "Reciprocal averaging" = "RA")),
      sliderInput("ng", "Class intervals", min = 2, max = 16, value = 8),
      uiOutput("group_ui"),
      actionButton("run", "Run analysis", class = "btn-primary"),
      tags$hr(),
      tags$small("Estimation is pairwise conditional; person estimates are Warm's WLE.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Items & fit", br(),
                 verbatimTextOutput("overview"), tableOutput("items")),
        tabPanel("Thresholds", br(), verbatimTextOutput("thr")),
        tabPanel("ICC", br(), uiOutput("sel_icc"), plotOutput("icc", height = "460px")),
        tabPanel("Category curves", br(), uiOutput("sel_ccc"), plotOutput("ccc", height = "460px")),
        tabPanel("Person-item map", br(), plotOutput("pimap", height = "620px")),
        tabPanel("Dimensionality", br(), verbatimTextOutput("dim"), tableOutput("loadings")),
        tabPanel("Local dependence", br(), plotOutput("rcor", height = "460px"), tableOutput("rpairs")),
        tabPanel("DIF", br(), helpText("Requires a group column selected in the sidebar."),
                 tableOutput("dif"))
      )
    )
  )
)

server <- function(input, output, session) {

  raw_data <- reactive({
    if (isTRUE(input$demo) || is.null(input$file)) .demo_data()
    else read.csv(input$file$datapath, check.names = FALSE)
  })

  output$group_ui <- renderUI({
    df <- raw_data()
    selectInput("group", "DIF group column (optional)", c("(none)", names(df)),
                selected = if ("group" %in% names(df)) "group" else "(none)")
  })

  analysis <- eventReactive(input$run, {
    df <- raw_data(); grpcol <- input$group
    if (!is.null(grpcol) && grpcol != "(none)" && grpcol %in% names(df)) {
      g <- df[[grpcol]]; X <- as.matrix(df[, setdiff(names(df), grpcol), drop = FALSE])
    } else { g <- NULL; X <- as.matrix(df) }
    X <- apply(X, 2, function(col) suppressWarnings(as.integer(col)))
    keep <- apply(X, 2, function(col) length(unique(col[!is.na(col)])) > 1)
    X <- X[, keep, drop = FALSE]
    fit <- rasch(X, model = input$model, solver = input$solver, n_groups = input$ng)
    list(fit = fit, group = g)
  })

  output$overview <- renderPrint({
    a <- analysis(); f <- a$fit
    cat(sprintf("Model: %s   Items: %d   Persons: %d   Class intervals: %d\n",
                f$model, nrow(f$items), length(f$theta_person), input$ng))
    cat(sprintf("Person Separation Index (PSI): %.3f   Separation ratio: %.2f\n",
                f$psi$PSI, f$psi$separation))
    cat(sprintf("Total item-trait chi-square: %.1f on %d df (p = %.4f)\n",
                f$total_chisq, f$total_df,
                pchisq(f$total_chisq, f$total_df, lower.tail = FALSE)))
    cat(sprintf("Items flagged misfitting (Bonferroni): %d\n", sum(f$items$misfit)))
  })

  output$items <- renderTable({
    f <- analysis()$fit
    d <- f$items
    data.frame(Item = d$item, Max = d$max, Location = round(d$location, 3),
               InfitMS = round(d$infit_ms, 2), OutfitMS = round(d$outfit_ms, 2),
               FitResid = round(d$fit_resid, 2), ChiSq = round(d$chisq, 1),
               p = signif(d$chisq_p, 3), Misfit = ifelse(d$misfit, "*", ""))
  })

  output$thr <- renderPrint({
    for (d in analysis()$fit$thresholds_diag) {
      cat(sprintf("%-8s  ordered=%-5s  thresholds = [%s]\n", d$item, d$ordered,
                  paste(sprintf("%.2f", d$thresholds), collapse = ", ")))
      if (length(d$never_modal_categories))
        cat(sprintf("           categories never most-probable: %s\n",
                    paste(d$never_modal_categories, collapse = ", ")))
    }
  })

  output$sel_icc <- renderUI(selectInput("icc_item", "Item", analysis()$fit$items$item))
  output$sel_ccc <- renderUI(selectInput("ccc_item", "Item", analysis()$fit$items$item))
  output$icc <- renderPlot({ req(input$icc_item); plot_icc(analysis()$fit, input$icc_item) })
  output$ccc <- renderPlot({ req(input$ccc_item); plot_ccc(analysis()$fit, input$ccc_item) })
  output$pimap <- renderPlot(plot_pimap(analysis()$fit))

  output$dim <- renderPrint({
    dt <- dimensionality_test(analysis()$fit)
    if (!is.null(dt$note)) { cat(dt$note); return(invisible()) }
    cat(sprintf("Proportion of significant person t-tests: %.3f  (95%% CI %.3f to %.3f)\n",
                dt$prop_significant, dt$ci[1], dt$ci[2]))
    cat(sprintf("First residual eigenvalue: %.2f\n", dt$first_eigenvalue))
    cat(sprintf("Verdict: %s\n", if (dt$multidimensional)
      "lower CI exceeds 5% — unidimensionality is questionable" else
      "consistent with unidimensionality"))
  })
  output$loadings <- renderTable({
    pc <- residual_pca(analysis()$fit)$loadings
    data.frame(Item = pc$item, PC1_loading = round(pc$pc1_loading, 3))
  })

  output$rcor <- renderPlot(plot_resid_cor(analysis()$fit))
  output$rpairs <- renderTable({
    fl <- residual_correlations(analysis()$fit)$flagged
    if (!nrow(fl)) return(data.frame(Note = "no item pairs exceed the flag threshold"))
    data.frame(Item_A = fl$item_a, Item_B = fl$item_b,
               ResidCor = round(fl$resid_cor, 3), Excess = round(fl$excess, 3))
  })

  output$dif <- renderTable({
    a <- analysis(); req(!is.null(a$group))
    d <- dif_anova(a$fit, group = a$group, n_groups = input$ng)
    data.frame(Item = d$item,
               F_uniform = round(d$F_uniform, 2), p_uniform = signif(d$p_uniform, 3),
               Uniform_DIF = ifelse(d$uniform_DIF, "*", ""),
               F_nonuniform = round(d$F_nonuniform, 2), p_nonuniform = signif(d$p_nonuniform, 3),
               Nonuniform_DIF = ifelse(d$nonuniform_DIF, "*", ""))
  })
}

shinyApp(ui, server)
