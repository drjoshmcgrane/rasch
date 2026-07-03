# rmt Shiny GUI
# ---------------------------------------------------------------------------
# A modern bslib interface to the full rmt analysis: data upload with ID,
# person-factor, and item column nomination; pairwise conditional ML
# estimation (Andrich & Luo 2003); the complete test-of-fit suite;
# every diagnostic plot with per-plot PNG and PDF downloads; and one-click
# export of all tables and plots as a ZIP archive.
# Launch with rmt::run_app(), or shiny::runApp() from this folder.
# ---------------------------------------------------------------------------
library(shiny)
library(bslib)
library(DT)

if (requireNamespace("rmt", quietly = TRUE)) {
  library(rmt)
} else {
  rdir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
  if (dir.exists(rdir)) {
    for (f in list.files(rdir, "\\.R$", full.names = TRUE)) source(f)
  } else stop("Install rmt, or run the app from inst/shiny in the source tree")
}

# --- demo data: 10 polytomous items, one disordered, DIF on Q05 -------------
.demo_data <- function(seed = 11, Np = 1200) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  mvec <- rep(c(2, 3), length.out = 10)
  tau_true <- lapply(mvec, function(m) sort(rnorm(m, 0, 0.9)))
  tau_true[[2]] <- c(1.2, -1.3, 0.6)                       # disordered item
  th <- rnorm(Np, 0, 1.4)
  grp <- rep(c("reference", "focal"), each = Np / 2)
  sex <- sample(c("female", "male"), Np, replace = TRUE)
  X <- sapply(seq_along(mvec), function(i) {
    sft <- if (i == 5) ifelse(grp == "focal", 0.9, 0) else numeric(Np)  # uniform DIF
    sapply(seq_len(Np), function(n) sample(0:mvec[i], 1, prob = simP(th[n] - sft[n], tau_true[[i]])))
  })
  colnames(X) <- sprintf("Q%02d", seq_along(mvec))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
             group = grp, sex = sex, check.names = FALSE)
}

# dichotomous demo: 15 multiple-choice items (raw A-D responses), DIF planted
# on I05 by group, and I07 deliberately miskeyed (true correct C, key says A)
.demo_dich <- function(seed = 41, Np = 1000) {
  set.seed(seed)
  d <- seq(-2, 2, length.out = 15)
  grp <- rep(c("reference", "focal"), each = Np / 2)
  sex <- sample(c("female", "male"), Np, replace = TRUE)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(seq_along(d), function(i) {
    sft <- if (i == 5) ifelse(grp == "focal", 0.8, 0) else 0
    correct <- if (i == 7) "C" else "A"
    ok <- rbinom(Np, 1, plogis(th - d[i] - sft))
    ifelse(ok == 1, correct,
           sample(setdiff(c("A", "B", "C", "D"), correct), Np, replace = TRUE))
  })
  colnames(X) <- sprintf("I%02d", seq_along(d))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
             group = grp, sex = sex, check.names = FALSE)
}

# the demo key: all "A" (so I07 is the discoverable miskey)
.demo_dich_key <- function()
  setNames(rep("A", 15), sprintf("I%02d", 1:15))

# paired-comparison demo: 8 essays compared pairwise by 10 judges, with
# judge J09 answering at random (discoverable in the judge fit table)
.demo_btl <- function(seed = 47, reps = 22) {
  set.seed(seed)
  beta <- setNames(seq(-1.4, 1.4, length.out = 8), sprintf("E%02d", 1:8))
  pr <- t(utils::combn(names(beta), 2))
  d <- data.frame(object_a = rep(pr[, 1], each = reps),
                  object_b = rep(pr[, 2], each = reps))
  d$judge <- sprintf("J%02d", sample(1:10, nrow(d), replace = TRUE))
  p <- plogis(beta[d$object_a] - beta[d$object_b])
  p[d$judge == "J09"] <- 0.5
  d$winner <- ifelse(runif(nrow(d)) < p, d$object_a, d$object_b)
  d[sample(nrow(d)), ]
}

# rating scale demo: common step structure, item locations vary
.demo_rsm <- function(seed = 51, Np = 1000) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  loc <- seq(-1.2, 1.2, length.out = 8)
  step <- c(-0.9, 0.0, 0.9)
  grp <- rep(c("reference", "focal"), each = Np / 2)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(loc, function(b) sapply(th, function(t)
    sample(0:3, 1, prob = simP(t, b + step))))
  colnames(X) <- sprintf("R%02d", seq_along(loc))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
             group = grp, check.names = FALSE)
}

# long-format rated demo: 5 items, 6 raters (one erratic), incomplete design
.demo_long <- function(seed = 21, Np = 250) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  persons <- sprintf("P%04d", seq_len(Np)); raters <- paste0("Rater_", 1:6)
  th <- setNames(rnorm(Np, 0, 1.3), persons)
  rho <- setNames(c(-0.9, -0.4, -0.1, 0.1, 0.4, 0.9), raters)
  tau <- list(Essay = c(-1.2, 0.2, 1.1), Argument = c(-0.8, 0.5, 1.3),
              Evidence = c(-1.5, -0.2, 0.9), Style = c(-0.6, 0.4, 1.2),
              Mechanics = c(-1.0, 0.0, 1.0))
  d <- expand.grid(person = persons, item = names(tau), rater = raters,
                   stringsAsFactors = FALSE)
  seen <- unlist(lapply(persons, function(p) paste(p, sample(raters, 3))))
  d <- d[paste(d$person, d$rater) %in% seen, ]
  d$score <- mapply(function(p, i, r) {
    if (r == "Rater_6" && runif(1) < 0.2) return(sample(0:3, 1))  # erratic rater
    sample(0:3, 1, prob = simP(th[p], tau[[i]] + rho[r]))
  }, d$person, d$item, d$rater)
  rownames(d) <- NULL
  d
}

# frames demo: 2 person groups x 3 item sets with distinct units
.demo_efrm <- function(seed = 31, per_g = 350) {
  set.seed(seed)
  simP <- function(th, tau, r) { x <- 0:length(tau); p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p) }
  glev <- c("year5", "year7"); grp <- rep(glev, each = per_g); Np <- length(grp)
  phi <- c(year5 = 0.8, year7 = 1.25)
  sets <- rep(c("Number", "Algebra", "Space"), each = 6)
  alpha <- c(Number = 0.75, Algebra = 1.0, Space = 4 / 3)
  th <- rnorm(Np, 0, 1.3) + ifelse(grp == "year7", 0.5, 0)
  d <- as.numeric(sapply(c(-0.3, 0.1, 0.2), function(m) m + seq(-1.2, 1.2, length.out = 6)))
  X <- sapply(seq_along(sets), function(i) sapply(seq_len(Np), function(n)
    sample(0:2, 1, prob = simP(th[n], d[i] + c(-0.5, 0.5),
                               alpha[sets[i]] * phi[grp[n]]))))
  colnames(X) <- sprintf("%s_%02d", sets, seq_along(sets))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X, year_group = grp,
             check.names = FALSE)
}

NONE <- "(none)"
# defined here (not mid-server) because observers created early reference it
FACTORIAL <- "(all factors: factorial)"

# p-values as text: "%.3f" alone prints a misleading 0.000 for tiny p
fmt_p <- function(p)
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))

theme <- bs_theme(
  version = 5, bg = "#f8fafc", fg = "#0f172a",
  primary = "#2563eb", secondary = "#64748b",
  success = "#0f766e", danger = "#dc2626", warning = "#f59e0b",
  "navbar-bg" = "#0f172a", "card-border-color" = "#e2e8f0",
  "border-radius" = "0.65rem", "font-size-base" = "0.95rem"
)

css <- HTML("
  .card { box-shadow: 0 1px 3px rgba(15,23,42,.08); }
  .card-header { background: #fff; font-weight: 600; border-bottom: 1px solid #e2e8f0; }
  .navbar-brand { font-weight: 700; letter-spacing: .02em; }
  .value-box-title { font-size: .72rem; text-transform: uppercase; letter-spacing: .04em; white-space: nowrap; }
  .value-box-value { font-size: 1.45rem; }
  pre, .shiny-text-output { white-space: pre-wrap; font-size: .82rem; }
  .btn-xs { padding: .1rem .5rem; font-size: .75rem; }
  .form-label { font-weight: 600; font-size: .85rem; }
  table.dataTable { font-size: .85rem; }
")

# Card with a plot and PNG/PDF download buttons in the header. The body is
# non-fillable so flex sizing can never compress the fixed-height plot
# (the cause of the squashed plots), and a percentage height is avoided
# because it races the layout and renders a zero-height device.
plotCard <- function(id, title, height = "560px") {
  card(
    full_screen = TRUE,
    card_header(div(class = "d-flex justify-content-between align-items-center",
      span(title),
      div(class = "btn-group",
          downloadButton(paste0(id, "_png"), "PNG", class = "btn-outline-secondary btn-xs"),
          downloadButton(paste0(id, "_pdf"), "PDF", class = "btn-outline-secondary btn-xs")))),
    card_body(plotOutput(id, height = height), padding = 8, fillable = FALSE)
  )
}

tableCard <- function(id, title, note = NULL) {
  card(
    full_screen = TRUE,
    card_header(div(class = "d-flex justify-content-between align-items-center",
      span(title),
      downloadButton(paste0(id, "_csv"), "CSV", class = "btn-outline-secondary btn-xs"))),
    card_body(if (!is.null(note)) p(class = "text-muted small mb-2", note),
              DTOutput(id), padding = 12)
  )
}

ui <- page_navbar(
  title = span("rmt"),
  theme = theme,
  # normal scrolling pages: never compress content to fit the viewport
  fillable = FALSE,
  header = tags$head(tags$style(css)),

  # ----------------------------------------------------------------- DATA --
  nav_panel("Data",
    layout_sidebar(
      sidebar = sidebar(width = 330,
        h6("Data source"),
        fileInput("file", NULL, accept = c(".csv", ".txt", ".tsv"),
                  buttonLabel = "Browse…", placeholder = "CSV / TSV file"),
        selectInput("demo_choice", "Or pick an example dataset",
                    c("(none)" = "none",
                      "Multiple choice, dichotomous" = "dich",
                      "Polytomous (PCM)" = "pcm",
                      "Rating scale (RSM)" = "rsm",
                      "Ratings, long format (MFRM)" = "mfrm",
                      "Item sets x groups (EFRM)" = "efrm",
                      "Paired comparisons (BTL)" = "btl")),
        radioButtons("model_type", "Model",
                     c("Dichotomous" = "dich",
                       "Partial credit (PCM)" = "pcm",
                       "Rating scale (RSM)" = "rsm",
                       "Many-facet (MFRM)" = "mfrm",
                       "Extended frames (EFRM)" = "efrm",
                       "Paired comparisons (BTL)" = "btl")),
        hr(),
        conditionalPanel("['dich','pcm','rsm'].indexOf(input.model_type) > -1",
          h6("Column roles"),
          selectInput("id_col", "ID variable", NONE),
          selectizeInput("factor_cols", "Person factors (DIF groups)", NULL,
                         multiple = TRUE,
                         options = list(placeholder = "none selected")),
          selectizeInput("item_cols", "Item columns", NULL, multiple = TRUE,
                         options = list(placeholder = "all remaining columns")),
          hr(),
          fileInput("key_file", "Multiple-choice key (CSV: item,key — use \"A/C\" for double keys — or item,option,score for polytomous option scoring)",
                    accept = ".csv", placeholder = "optional"),
          fileInput("anchor_file", "Anchors for equating (CSV: item,k,tau)",
                    accept = ".csv", placeholder = "optional"),
          radioButtons("anchor_type", "Anchor as",
                       c("Individual thresholds" = "individual",
                         "Average item locations" = "average"),
                       inline = TRUE),
          p(class = "text-muted small mt-1",
            "Anchors match by item name; rows for items not present are ignored. Individual anchoring fixes each listed threshold; average anchoring fixes each item's mean location (thresholds stay free). Save an anchor file from the Items page of a previous analysis.")
        ),
        conditionalPanel("input.model_type == 'efrm'",
          h6("Column roles"),
          selectInput("ef_id", "ID variable", NONE),
          selectInput("ef_group", "Person group column ((none) = single group; units differ by item set only)",
                      NONE),
          selectizeInput("ef_items", "Item columns", NULL, multiple = TRUE,
                         options = list(placeholder = "all remaining columns")),
          fileInput("ef_sets", "Item-set map (CSV: item,set)",
                    accept = ".csv", placeholder = "optional"),
          checkboxInput("ef_prefix", "Infer sets from item-name prefix", TRUE),
          selectInput("ef_se", "Standard errors",
                      c("Hybrid (fast)" = "hybrid",
                        "Full person bootstrap (slow, exact)" = "bootstrap")),
          numericInput("ef_reps", "Bootstrap replicates", value = 200,
                       min = 50, step = 50),
          p(class = "text-muted small",
            "Each item-set by group cell is a frame with its own unit. Group units come from the person-free pairwise comparisons; set units from persons common to the sets.")
        ),
        conditionalPanel("input.model_type == 'mfrm'",
          h6("Column roles (one row per response)"),
          selectInput("lp_person", "Person column", NONE),
          selectInput("lp_item", "Item column", NONE),
          selectInput("lp_score", "Score column", NONE),
          selectizeInput("lp_facets", "Facet columns (e.g. rater)", NULL,
                         multiple = TRUE,
                         options = list(placeholder = "choose at least one")),
          selectInput("lp_interaction", "Item-by-facet interaction (optional)", NONE),
          p(class = "text-muted small",
            "Each item x facet combination is calibrated jointly; facet severities are reported with SEs and fit. An interaction lets one facet be more or less severe on particular items.")
        ),
        conditionalPanel("input.model_type == 'btl'",
          h6("Column roles (one comparison per row)"),
          selectInput("bt_a", "Object A column", NONE),
          selectInput("bt_b", "Object B column", NONE),
          selectInput("bt_win", "Winner column", NONE),
          selectInput("bt_judge", "Judge column (optional)", NONE),
          selectInput("bt_count", "Count column (optional)", NONE),
          radioButtons("bt_ties", "Ties",
                       c("Drop" = "drop", "Half a win each" = "half")),
          p(class = "text-muted small",
            "The Bradley-Terry-Luce model: the conditional (person-free) form of the dichotomous Rasch model, estimated by the same conventions. A judge column enables the judge fit table and clusters the standard errors by judge. Results appear on the BTL tab.")
        ),
        conditionalPanel("['dich','pcm'].indexOf(input.model_type) > -1",
          radioButtons("thr_mode", "Threshold estimation",
                       c("Free thresholds" = "free",
                         "Principal components (Andrich)" = "pc")),
          conditionalPanel("input.thr_mode == 'pc'",
            selectInput("pc_rank", "Components",
                        c("Location only" = "1",
                          "+ spread (equal spread)" = "2",
                          "+ skewness" = "3",
                          "+ kurtosis (full PC)" = "4"),
                        selected = "4"),
            p(class = "text-muted small",
              "Thresholds follow a polynomial trend across categories; useful with sparse categories. Anchors cannot be combined with this option."))),
        checkboxInput("ng_auto", "Automatic class intervals (at least 50 per interval)", TRUE),
        conditionalPanel("!input.ng_auto",
          sliderInput("ng", "Class intervals", min = 2, max = 16, value = 8)),
        numericInput("run_adjN",
                     "Adjust chi-square to sample size N (optional)",
                     value = NA, min = 50),
        h6("Estimation"),
        numericInput("maxit", "Maximum iterations", value = 60, min = 5, step = 5),
        numericInput("tol", "Convergence criterion", value = 1e-8,
                     min = 1e-12, step = 1e-8),
        actionButton("run", "Run analysis", class = "btn-primary w-100 btn-lg mt-2"),
        conditionalPanel("output.has_override",
          uiOutput("override_status"),
          actionButton("reset_override", "Reset overrides",
                       class = "btn-outline-warning w-100 mt-1")),
        p(class = "text-muted small mt-3",
          "Estimation: pairwise conditional maximum likelihood (Andrich & Luo 2003).",
          "Person measures: Warm weighted likelihood.")
      ),
      card(card_header("Data preview"),
           card_body(uiOutput("data_info"), DTOutput("preview"), padding = 12))
    )
  ),

  # -------------------------------------------------------------- SUMMARY --
  nav_panel("Summary",
    uiOutput("vboxes"),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      card(card_header("Test of fit"), card_body(verbatimTextOutput("fit_summary"))),
      card(card_header("Targeting & reliability"), card_body(verbatimTextOutput("targeting")))
    ),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      card(
        full_screen = TRUE,
        card_header(div(class = "d-flex justify-content-between align-items-center",
          span("Score-to-measure table"),
          downloadButton("score_tbl_csv", "CSV", class = "btn-outline-secondary btn-xs"))),
        card_body(
          p(class = "text-muted small mb-2",
            "Location and SE for every raw score (complete responders), with the frequency and cumulative percentage at each score."),
          div(class = "d-flex gap-3 flex-wrap",
            selectInput("st_method", "Estimator",
                        c("WLE (Warm)" = "wle", "MLE" = "mle"), width = "170px"),
            selectInput("st_extremes", "Extreme scores",
                        c("Model" = "model",
                          "Extrapolated (geometric)" = "extrapolated"),
                        width = "200px")),
          DTOutput("score_tbl"), padding = 12)),
      tableCard("thr_tbl", "Thresholds with standard errors")
    ),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      card(
        full_screen = TRUE,
        card_header(div(class = "d-flex justify-content-between align-items-center",
          span("Traditional statistics (CTT)"),
          downloadButton("ctt_tbl_csv", "CSV", class = "btn-outline-secondary btn-xs"))),
        card_body(uiOutput("ctt_head"), DTOutput("ctt_tbl"), padding = 12)),
      uiOutput("lr_ui")
    )
  ),

  # ---------------------------------------------------------------- ITEMS --
  nav_panel("Items",
    div(class = "mb-2 d-flex align-items-center gap-3 flex-wrap",
        checkboxInput("show_obs",
                      "Show observed points on the category and threshold curves",
                      TRUE, width = "420px"),
        downloadButton("dl_anchors", "Save anchors (CSV: item,k,tau)",
                       class = "btn-outline-secondary mb-3")),
    tableCard("items_tbl", "Item statistics",
              "Click a row to inspect that item's chi-square detail and curves below. Location and SE from the pairwise conditional likelihood; fit residual ~ N(0,1) under fit; item-trait chi-square over class intervals; ANOVA F over class-interval cells; misfit flag uses BH-adjusted probabilities."),
    uiOutput("pc_comp_ui"),
    card(full_screen = TRUE,
      card_header(div(class = "d-flex justify-content-between align-items-center",
        span("Chi-square detail (selected item)"),
        div(class = "btn-group",
            downloadButton("chisq_int_csv", "Intervals CSV",
                           class = "btn-outline-secondary btn-xs"),
            downloadButton("chisq_cat_csv", "Categories CSV",
                           class = "btn-outline-secondary btn-xs")))),
      card_body(
        uiOutput("chisq_caption"),
        DTOutput("chisq_int_tbl"),
        h6("Response categories by class interval", class = "mt-3"),
        DTOutput("chisq_cat_tbl"), padding = 12)),
    layout_columns(col_widths = 12,
      plotCard("icc", "Item characteristic curve"),
      plotCard("ccc", "Category probability curves")),
    layout_columns(col_widths = 12,
      plotCard("tpc", "Threshold probability curves"),
      plotCard("cfreq", "Category frequencies")),
    layout_columns(col_widths = 12,
      tableCard("distractor_tbl", "Distractor analysis",
                "Multiple-choice analyses only (provide a key). Locations use the rest measure; a distractor whose takers are abler than the keyed option's flags a possible miskey."),
      plotCard("distractor_plot", "Option curves"),
      card(card_header("Polytomous option scoring (Andrich & Styles 2011)"),
           card_body(
             p(class = "text-muted",
               "Propose partial credit for informative distractors from the rest-measure evidence. Review substantively, download, edit if needed, and upload as the key (item,option,score) to refit."),
             layout_columns(col_widths = c(3, 3, 6),
               numericInput("rescore_min_n", "Min takers", 20, min = 5, step = 5),
               numericInput("rescore_z", "Separation z", 1.96, min = 0.5, step = 0.1),
               actionButton("rescore_go", "Propose option scores",
                            class = "btn-primary mt-4")),
             DT::DTOutput("rescore_tbl"),
             downloadButton("dl_rescore", "Download proposed key CSV",
                            class = "btn-sm"))))
  ),

  # -------------------------------------------------------------- PERSONS --
  nav_panel("Persons",
    tableCard("person_tbl", "Person estimates",
              "Warm WLE location and SE per person, with raw score, fit statistics, and your ID and factor columns. Click a row to draw that person's characteristic curve below."),
    layout_columns(col_widths = 12,
      plotCard("pcc", "Person characteristic curve (selected person)"),
      card(
        full_screen = TRUE,
        card_header(div(class = "d-flex justify-content-between align-items-center",
          span("Fit residual distribution"),
          div(class = "btn-group",
              downloadButton("rdist_png", "PNG", class = "btn-outline-secondary btn-xs"),
              downloadButton("rdist_pdf", "PDF", class = "btn-outline-secondary btn-xs")))),
        card_body(
          div(class = "d-flex gap-4 flex-wrap",
            radioButtons("rd_what", NULL,
                         c("Items" = "items", "Persons" = "persons"),
                         selected = "persons", inline = TRUE),
            radioButtons("rd_stat", NULL,
                         c("Fit residual (log-transformed)" = "fit_resid",
                           "Natural" = "natural"), inline = TRUE)),
          plotOutput("rdist", height = "520px"), padding = 8, fillable = FALSE))),
    layout_columns(col_widths = 12,
      plotCard("pfit", "Person fit"),
      plotCard("pim_p", "Person-item threshold distribution"))
  ),

  # ----------------------------------------------------------------- TEST --
  nav_panel("Test plots",
    layout_columns(col_widths = 12,
      plotCard("thrmap", "Threshold map"),
      plotCard("imap", "Item map: location by fit residual")),
    layout_columns(col_widths = 12,
      plotCard("tcc", "Test characteristic curve"),
      plotCard("tif", "Test information & SEM")),
    plotCard("guttman", "Guttman scalogram", height = "640px")
  ),

  # ------------------------------------------------------------------ DIF --
  nav_panel("DIF",
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        selectInput("dif_factor", "Person factor", NONE),
        conditionalPanel("input.dif_factor == '(all factors: factorial)'",
          radioButtons("dif_effects", "Model",
                       c("Full factorial" = "factorial",
                         "Main effects only" = "main"))),
        selectInput("dif_item", "Item for ICC by group", NONE),
        numericInput("dif_alpha", "Significance level (alpha)", value = 0.05,
                     min = 0.001, max = 0.5, step = 0.01),
        selectInput("dif_padj", "Multiplicity adjustment",
                    c("Benjamini-Hochberg" = "BH",
                      "Bonferroni" = "bonferroni",
                      "None" = "none")),
        p(class = "text-muted small",
          "ANOVA of standardised residuals: factor effects = uniform DIF; factor x class-interval terms = non-uniform DIF. Probabilities are adjusted across items by the chosen method. With several factors, choose the factorial option to model them jointly: significant interactions supersede their main effects, and Tukey HSD compares the levels of each significant group term."),
        hr(),
        conditionalPanel("input.dif_factor != '(all factors: factorial)'",
          actionButton("make_split", "Resolve: split this item by this factor",
                       class = "btn-outline-primary w-100"),
          p(class = "text-muted small mt-2",
            "Replaces the selected item with one item per group level (each level keeps only its own responses) and re-analyses; the split locations quantify the DIF. Splitting works one factor at a time; choose a single factor above."))),
      tableCard("dif_tbl", "DIF analysis of variance"),
      tableCard("dif_tukey_tbl", "Tukey HSD comparisons",
                "Pairwise level comparisons for significant, non-superseded group terms (factorial mode)."),
      card(card_header("DIF size in logits (practical significance)"),
           card_body(
             p(class = "text-muted",
               "Resolves the selected item by the selected factor (in factorial mode: by every significant, non-superseded group term) and reports pairwise location differences in logits with Holm familywise adjustment. Differences of at least the criterion are flagged as practically significant."),
             layout_columns(col_widths = c(4, 4, 4),
               numericInput("dif_size_flag", "Practical criterion (logits)",
                            0.5, min = 0.1, step = 0.1),
               numericInput("dif_size_minn", "Min responders per level", 20,
                            min = 5, step = 5),
               actionButton("dif_size_go", "Compute DIF size",
                            class = "btn-primary mt-4")),
             DT::DTOutput("dif_size_tbl"),
             downloadButton("dl_dif_size", "Download CSV", class = "btn-sm"))),
      plotCard("dif_icc", "ICC by group (DIF plot)")
    )
  ),

  # --------------------------------------------------------------- FACETS --
  nav_panel("Facets",
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        selectInput("facet_sel", "Facet", NONE),
        p(class = "text-muted small",
          "Severities from the joint calibration (positive = more severe). Pooled fit residuals beyond +/-2.5 flag inconsistent levels. Long-format analyses only.")),
      tableCard("facet_tbl", "Facet severities and fit"),
      plotCard("facet_plot", "Severity caterpillar plot"),
      tableCard("facet_int_tbl", "Item-by-facet interactions",
                "Shown when the analysis was run in interactive facet mode; gamma is the extra severity of a level on a particular item.")
    )
  ),

  # -------------------------------------------------------------- EQUATING --
  nav_panel("Equating",
    layout_sidebar(
      sidebar = sidebar(width = 300, open = "always",
        radioButtons("eq_source", "Reference",
                     c("Uploaded calibration CSV" = "csv",
                       "A kept fit from Compare" = "kept")),
        conditionalPanel("input.eq_source == 'csv'",
          fileInput("eq_file", "Reference calibration (CSV: item,location,se)",
                    accept = ".csv")),
        conditionalPanel("input.eq_source == 'kept'",
          selectInput("eq_kept", "Kept fit", NONE)),
        radioButtons("eq_shift", "Scale alignment",
                     c("Allow a shift between origins" = "mean",
                       "Compare raw locations (anchored scales)" = "none")),
        downloadButton("dl_calib", "Save current calibration (CSV)",
                       class = "btn-outline-secondary w-100"),
        p(class = "text-muted small mt-2",
          "Common items (matched by name) are tested against the shifted identity line; flagged items show drift and weaken the equating link. Save a calibration now to equate a future analysis against it.")),
      tableCard("eq_tbl", "Common-item comparison"),
      plotCard("eq_plot", "Equating plot")
    )
  ),

  # --------------------------------------------------------------- FRAMES --
  nav_panel("Frames",
    layout_sidebar(
      sidebar = sidebar(width = 290, open = "always",
        selectInput("frame_item", "Item for ICC across frames", NONE),
        p(class = "text-muted small",
          "Units rho = alpha (set) x phi (group) on a common arbitrary scale. Within a frame all curves are parallel; across frames they fan with the unit. Extended frame of reference analyses only.")),
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(7, 5)),
        tableCard("frame_tbl", "Frames: units, origins, pooled fit"),
        div(tableCard("phi_tbl", "Person group units (phi)"),
            tableCard("alpha_tbl", "Item set units (alpha) and locations"))),
      layout_columns(col_widths = 12,
        plotCard("frame_plot", "Frame units"),
        plotCard("frame_icc", "ICC across frames")),
      card(card_header("Equal-unit comparison"),
           card_body(verbatimTextOutput("efrm_cmp")))
    )
  ),

  # ------------------------------------------------------- DIMENSIONALITY --
  nav_panel("Dimensionality",
    layout_sidebar(
      sidebar = sidebar(width = 300, open = "always",
        h6("t-test item subsets"),
        selectizeInput("dim_pos", "Subset A", NULL, multiple = TRUE,
                       options = list(placeholder = "default: positive PC1 loadings")),
        selectizeInput("dim_neg", "Subset B", NULL, multiple = TRUE,
                       options = list(placeholder = "default: negative PC1 loadings")),
        actionButton("dim_apply", "Run t-test with these subsets",
                     class = "btn-outline-primary w-100"),
        p(class = "text-muted small mt-2",
          "Leave both empty (and press the button) to return to the first-contrast split. Persons extreme on either subset are excluded; the proportion of significant tests carries an exact binomial confidence interval.")),
      layout_columns(col_widths = 12,
        card(card_header("Unidimensionality t-test (Smith)"),
             card_body(verbatimTextOutput("dim_txt"))),
        plotCard("scree", "Scree of the residual components")),
      layout_columns(col_widths = 12,
        plotCard("pca_plot", "Residual first contrast"),
        tableCard("loadings_tbl", "Component loadings (first 10)")),
      tableCard("eigen_tbl", "Residual eigenvalues (first 10)"),
      card(
        full_screen = TRUE,
        card_header(div(class = "d-flex justify-content-between align-items-center",
          span("Magnitude of multidimensionality (Andrich 2016)"),
          downloadButton("dm_tbl_csv", "CSV", class = "btn-outline-secondary btn-xs"))),
        card_body(
          p(class = "text-muted small",
            "Compares reliability with all items treated as independent (run1) against the subtest analysis in which each subset becomes one polytomous super-item. c is the unique-variance loading, rho the latent correlation between the subsets, and A the proportion of common variance. Uses the manual subsets above if set, otherwise the current PC1 split; every item must belong to a subset."),
          actionButton("dm_run", "Estimate from current subsets",
                       class = "btn-outline-primary"),
          DTOutput("dm_tbl"), padding = 12))
    )
  ),

  # ------------------------------------------------------ LOCAL DEPENDENCE --
  nav_panel("Local dependence",
    layout_columns(col_widths = 12,
      plotCard("rcor", "Residual correlations", height = "640px"),
      div(
        numericInput("ld_flag",
                     "Flag threshold (excess over the average residual correlation)",
                     value = 0.2, min = 0.05, max = 0.9, step = 0.05,
                     width = "420px"),
        tableCard("rpairs_tbl", "Flagged dependent pairs",
                  "Pairs more than the flag threshold above the average off-diagonal residual correlation."),
        card(card_header("Subtest (combine dependent items)"),
          card_body(
            p(class = "text-muted small",
              "Select two or more items to merge into one polytomous super-item and re-analyse; the dependence is absorbed into the subtest."),
            selectizeInput("subtest_items", NULL, NULL, multiple = TRUE,
                           options = list(placeholder = "items to combine")),
            actionButton("make_subtest", "Combine and re-analyse",
                         class = "btn-outline-primary w-100"),
            uiOutput("subtest_status"))))),
    layout_columns(col_widths = 12,
      card(
        full_screen = TRUE,
        card_header(div(class = "d-flex justify-content-between align-items-center",
          span("Response dependence magnitude (Andrich & Kreiner)"),
          downloadButton("dep_tbl_csv", "CSV", class = "btn-outline-secondary btn-xs"))),
        card_body(
          p(class = "text-muted small",
            "Resolves the dependent item by the categories of the independent item and re-analyses; d is the size of the dependence in logits, half the split of the resolved thresholds. Both items must share the same maximum score."),
          div(class = "d-flex gap-3 flex-wrap align-items-end",
            selectInput("dep_item", "Dependent item", NONE, width = "190px"),
            selectInput("ind_item", "Independent item", NONE, width = "190px"),
            div(class = "mb-3",
                actionButton("run_dep", "Estimate d", class = "btn-outline-primary"))),
          verbatimTextOutput("dep_txt"),
          DTOutput("dep_tbl"), padding = 12)),
      card(
        full_screen = TRUE,
        card_header(div(class = "d-flex justify-content-between align-items-center",
          span("Spread test (LUB)"),
          downloadButton("spread_tbl_csv", "CSV", class = "btn-outline-secondary btn-xs"))),
        card_body(
          actionButton("run_spread", "Run spread test",
                       class = "btn-outline-primary mb-2"),
          DTOutput("spread_tbl"),
          p(class = "text-muted small mt-2",
            "Spread below the least upper bound indicates dependence among subtest members (Andrich 1985). Polytomous items only; typically applied after combining items into a subtest."),
          padding = 12)))
  ),

  # ------------------------------------------------------------- GUESSING --
  nav_panel("Guessing",
    layout_sidebar(
      sidebar = sidebar(width = 300, open = "always",
        numericInput("guess_chance", "Chance success probability",
                     value = 0.25, min = 0.05, max = 0.95, step = 0.05),
        selectizeInput("guess_anchors", "Anchor items (common origin)", NULL,
                       multiple = TRUE,
                       options = list(placeholder = "automatic: least-affected third")),
        actionButton("run_guess", "Run tailored analysis",
                     class = "btn-primary w-100"),
        p(class = "text-muted small mt-2",
          "The tailored procedure of Andrich, Marais and Humphry (2012): every response whose modelled success probability falls below the chance level is set to missing and the test is re-calibrated on a common origin. Difficult items becoming harder in the tailored calibration signals guessing. Dichotomous analyses only.")),
      layout_columns(col_widths = 12,
        card(card_header("Tailored analysis"),
             card_body(verbatimTextOutput("guess_txt"))),
        tableCard("guess_tbl", "Initial vs tailored calibration",
                  "shift = tailored minus origin-equated location; z > 1.96 flags items significantly harder after tailoring (a guessing signature).")),
      plotCard("guess_plot", "Tailored vs origin-equated calibrations")
    )
  ),

  # ------------------------------------------------------------------ BTL --
  nav_panel("BTL",
    layout_columns(col_widths = 12,
      card(card_header("Paired comparisons (Bradley-Terry-Luce)"),
           card_body(uiOutput("btl_boxes"),
                     verbatimTextOutput("btl_summary"))),
      tableCard("btl_obj_tbl", "Object locations and fit",
                "Conditional (person-free) estimation with sum-zero identification and sandwich standard errors; the fit residual is the log-of-mean-square statistic over each object's comparisons (Andrich & Marais 2019)."),
      plotCard("btl_plot", "Object caterpillar"),
      tableCard("btl_pairs_tbl", "Pairwise goodness of fit",
                "Observed against expected win proportions for every pair; the total chi-square tests the BTL structure."),
      tableCard("btl_judges_tbl", "Judge fit",
                "Available when a judge column is nominated; an erratic judge carries a large positive fit residual, exactly as an erratic person does."))
  ),

  # -------------------------------------------------------------- COMPARE --
  nav_panel("Compare",
    layout_sidebar(
      sidebar = sidebar(width = 300, open = "always",
        actionButton("keep_fit", "Keep current fit for comparison",
                     class = "btn-primary w-100"),
        actionButton("clear_fits", "Clear kept fits",
                     class = "btn-outline-secondary w-100 mt-2"),
        selectInput("cmp_ref", "Reference fit (two_delta_ll)", NONE),
        p(class = "text-muted small mt-3",
          "Run an analysis, keep it, change the model or settings, run again, and keep that too. For fits of the same data the pairwise conditional log-likelihoods are compared directly (descriptive, composite likelihood; most meaningful for nested structures such as RSM inside PCM). Across different data preparations, compare the calibration-free columns: chi-square per df, fit residual SDs (ideal 1), PSI, and alpha.")),
      tableCard("cmp_tbl", "Model comparison",
                "Reference for two_delta_ll is the fit chosen in the sidebar.")
    )
  ),

  # --------------------------------------------------------------- EXPORT --
  nav_panel("Export",
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      card(card_header("Download everything"),
        card_body(
          p("One archive containing every table (CSV), every plot (in the formats chosen), and a plain-text analysis summary."),
          checkboxGroupInput("exp_formats", "Plot formats",
                             c("PNG" = "png", "PDF" = "pdf"), selected = c("png", "pdf"),
                             inline = TRUE),
          checkboxInput("exp_items", "Include the per-item plot set (ICC, categories, thresholds, frequencies for every item)", TRUE),
          downloadButton("dl_zip", "Download all results (ZIP)", class = "btn-primary btn-lg"))),
      card(card_header("What is included"),
        card_body(tags$ul(
          tags$li("Item statistics with the item ANOVA fit table, thresholds with SEs, person estimates (with ID and factors), score-to-measure table with score frequencies"),
          tags$li("Chi-square class-interval detail for every item, traditional (CTT) statistics, and the principal components table when estimated"),
          tags$li("Residual correlations, flagged dependent pairs, PCA loadings, category frequencies, DIF ANOVA for every factor"),
          tags$li("Person-item distribution, threshold map, TCC, TIF, item and person fit maps, item and person fit residual distributions, residual heatmap, PCA plot"),
          tags$li("Per-item ICC, category curves, threshold curves, and frequency charts"),
          tags$li("For many-facet analyses: facet severities with SEs and fit, structural item thresholds, and severity caterpillar plots"),
          tags$li("summary.txt with the full test-of-fit report"))))
    )
  )
)

server <- function(input, output, session) {

  # ------------------------------------------------------------- data in --
  # picking an example dataset also selects the matching model; uploading a
  # file clears the example selection
  observeEvent(input$demo_choice, {
    if (!identical(input$demo_choice, "none"))
      updateRadioButtons(session, "model_type", selected = input$demo_choice)
  }, ignoreInit = TRUE)
  observeEvent(input$file,
    updateSelectInput(session, "demo_choice", selected = "none"))

  raw_data <- reactive({
    if (!identical(input$demo_choice %||% "none", "none"))
      return(switch(input$demo_choice,
                    dich = .demo_dich(), rsm = .demo_rsm(),
                    mfrm = .demo_long(), efrm = .demo_efrm(),
                    btl = .demo_btl(), .demo_data()))
    req(input$file)
    ext <- tolower(tools::file_ext(input$file$name))
    sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
    read.csv(input$file$datapath, sep = sep, check.names = FALSE,
             stringsAsFactors = FALSE)
  })

  anchors_in <- reactive({
    if (is.null(input$anchor_file)) return(NULL)
    a <- tryCatch(read.csv(input$anchor_file$datapath, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    if (is.null(a) || !all(c("item", "k", "tau") %in% names(a))) {
      showNotification("Anchor file needs columns item, k, tau - ignored.",
                       type = "warning")
      return(NULL)
    }
    a
  })

  observeEvent(raw_data(), {
    df <- raw_data(); nm <- names(df)
    guess_id <- nm[grepl("^id$|_id$|^person", tolower(nm))][1]
    guess_fac <- intersect(nm, c("group", "sex", "gender", "site", "country", "age_group"))
    updateSelectInput(session, "id_col", choices = c(NONE, nm),
                      selected = if (!is.na(guess_id)) guess_id else NONE)
    updateSelectizeInput(session, "factor_cols", choices = nm, selected = guess_fac)
    updateSelectizeInput(session, "item_cols", choices = nm,
                         selected = setdiff(nm, c(guess_id, guess_fac)))
    # long-format guesses
    g_per <- nm[grepl("person|candidate|student|^id$|_id$", tolower(nm))][1]
    g_itm <- nm[grepl("item|task|criterion|question", tolower(nm))][1]
    g_sco <- nm[grepl("score|rating|grade|mark", tolower(nm))][1]
    g_fac <- setdiff(nm[grepl("rater|judge|marker|occasion|time", tolower(nm))],
                     c(g_per, g_itm, g_sco))
    updateSelectInput(session, "lp_person", choices = c(NONE, nm),
                      selected = if (!is.na(g_per)) g_per else NONE)
    updateSelectInput(session, "lp_item", choices = c(NONE, nm),
                      selected = if (!is.na(g_itm)) g_itm else NONE)
    updateSelectInput(session, "lp_score", choices = c(NONE, nm),
                      selected = if (!is.na(g_sco)) g_sco else NONE)
    updateSelectizeInput(session, "lp_facets", choices = nm, selected = g_fac)
    # frames layout guesses
    g_grp <- nm[grepl("group|year|grade|cohort|class$", tolower(nm))][1]
    updateSelectInput(session, "ef_id", choices = c(NONE, nm),
                      selected = if (!is.na(guess_id)) guess_id else NONE)
    updateSelectInput(session, "ef_group", choices = c(NONE, nm),
                      selected = if (!is.na(g_grp)) g_grp else NONE)
    updateSelectizeInput(session, "ef_items", choices = nm,
                         selected = setdiff(nm, c(guess_id, g_grp)))
    # paired-comparison guesses
    g_a <- nm[grepl("^a$|object_a|left|first|option_a", tolower(nm))][1]
    g_b <- nm[grepl("^b$|object_b|right|second|option_b", tolower(nm))][1]
    g_w <- nm[grepl("win|preferred|chosen|better", tolower(nm))][1]
    g_j <- nm[grepl("judge|rater|marker", tolower(nm))][1]
    g_c <- nm[grepl("^count$|^n$|freq", tolower(nm))][1]
    updateSelectInput(session, "bt_a", choices = c(NONE, nm),
                      selected = if (!is.na(g_a)) g_a else NONE)
    updateSelectInput(session, "bt_b", choices = c(NONE, nm),
                      selected = if (!is.na(g_b)) g_b else NONE)
    updateSelectInput(session, "bt_win", choices = c(NONE, nm),
                      selected = if (!is.na(g_w)) g_w else NONE)
    updateSelectInput(session, "bt_judge", choices = c(NONE, nm),
                      selected = if (!is.na(g_j)) g_j else NONE)
    updateSelectInput(session, "bt_count", choices = c(NONE, nm),
                      selected = if (!is.na(g_c)) g_c else NONE)
  })

  # item -> set map: uploaded CSV wins; otherwise infer from the item-name
  # prefix (the part before trailing digits/separators)
  ef_setmap <- reactive({
    its <- if (length(input$ef_items)) input$ef_items else
      setdiff(names(raw_data()),
              c(if (!is.null(input$ef_id) && input$ef_id != NONE) input$ef_id,
                if (!is.null(input$ef_group) && input$ef_group != NONE) input$ef_group))
    if (!is.null(input$ef_sets)) {
      mp <- tryCatch(read.csv(input$ef_sets$datapath, stringsAsFactors = FALSE),
                     error = function(e) NULL)
      if (!is.null(mp) && all(c("item", "set") %in% names(mp))) {
        out <- setNames(as.character(mp$set), mp$item)
        miss <- setdiff(its, names(out))
        if (length(miss)) out[miss] <- "(rest)"
        return(out[its])
      }
      showNotification("Item-set CSV needs columns item,set - using prefixes.",
                       type = "warning")
    }
    if (isTRUE(input$ef_prefix)) {
      pref <- sub("[_. -]*[0-9]+$", "", its)
      pref[pref == ""] <- "(rest)"
      return(setNames(pref, its))
    }
    setNames(rep("all", length(its)), its)
  })

  observeEvent(input$lp_facets, {
    sel <- if (!is.null(input$lp_interaction) &&
               input$lp_interaction %in% input$lp_facets)
      input$lp_interaction else NONE
    updateSelectInput(session, "lp_interaction",
                      choices = c(NONE, input$lp_facets), selected = sel)
  }, ignoreNULL = FALSE)

  # keep item choices free of the chosen ID / factor columns
  observeEvent(c(input$id_col, input$factor_cols), {
    df <- raw_data(); nm <- names(df)
    taken <- c(if (!is.null(input$id_col) && input$id_col != NONE) input$id_col,
               input$factor_cols)
    sel <- setdiff(if (length(input$item_cols)) input$item_cols else nm, taken)
    updateSelectizeInput(session, "item_cols", choices = setdiff(nm, taken),
                         selected = sel)
  }, ignoreInit = TRUE)

  output$data_info <- renderUI({
    if (identical(input$demo_choice %||% "none", "none") && is.null(input$file))
      return(p(class = "text-muted",
               "Upload a CSV/TSV file, or pick an example dataset in the sidebar, to begin."))
    df <- raw_data()
    p(class = "text-muted",
      sprintf("%d rows x %d columns. Nominate the column roles in the sidebar, then run the analysis. Missing responses may be left blank or coded as -1; any negative score is read as missing.",
              nrow(df), ncol(df)))
  })
  output$preview <- renderDT({
    datatable(head(raw_data(), 200), rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
  })

  # ----------------------------------------------------------------- fit --
  override_fit <- reactiveVal(NULL)
  override_desc <- reactiveVal(NULL)
  # clear any subtest/split override as soon as a fresh run is requested;
  # fit() short-circuits on the override, so analysis() cannot clear it itself
  observeEvent(input$run, { override_fit(NULL); override_desc(NULL) },
               priority = 10)
  output$has_override <- reactive(!is.null(override_fit()))
  outputOptions(output, "has_override", suspendWhenHidden = FALSE)
  output$override_status <- renderUI({
    if (is.null(override_desc())) return(NULL)
    p(class = "text-warning small mb-1 mt-2",
      paste("Active override -", override_desc()))
  })
  observeEvent(input$reset_override, {
    override_fit(NULL); override_desc(NULL)
    showNotification("Override cleared; showing the base analysis.",
                     type = "message", duration = 5)
  })

  analysis <- eventReactive(input$run, {
    df <- raw_data()
    # automatic class intervals pass NULL; rasch() resolves the rule and
    # reports the value used in fit$n_groups
    ng <- if (isTRUE(input$ng_auto %||% TRUE)) NULL else input$ng
    # chi-square sample-size adjustment is applied inside the fit, so every
    # tab and export sees the same adjusted statistics
    adjN <- if (!is.null(input$run_adjN) && !is.na(input$run_adjN) &&
                input$run_adjN > 0) input$run_adjN else NA
    withProgress(message = "Estimating (pairwise conditional ML)…", value = 0.3, {
      fit <- tryCatch({
        if (identical(input$model_type, "btl")) {
          if (any(c(input$bt_a, input$bt_b, input$bt_win) == NONE))
            stop("nominate the object A, object B, and winner columns")
          btl(df, object_a = input$bt_a, object_b = input$bt_b,
              winner = input$bt_win,
              judge = if (!is.null(input$bt_judge) && input$bt_judge != NONE)
                input$bt_judge else NULL,
              count = if (!is.null(input$bt_count) && input$bt_count != NONE)
                input$bt_count else NULL,
              ties = input$bt_ties %||% "drop",
              maxit = max(5, input$maxit %||% 60),
              tol = max(1e-12, input$tol %||% 1e-8))
        } else if (identical(input$model_type, "efrm")) {
          rasch_efrm(df,
                     item_sets = ef_setmap(),
                     groups = if (is.null(input$ef_group) ||
                                  input$ef_group == NONE)
                       rep("(all)", nrow(df)) else input$ef_group,
                     id = if (!is.null(input$ef_id) && input$ef_id != NONE)
                       input$ef_id else NULL,
                     items = names(ef_setmap()),
                     n_groups = ng, adjust_N = adjN,
                     maxit = max(5, input$maxit %||% 60),
                     tol = max(1e-12, input$tol %||% 1e-8),
                     se_method = input$ef_se %||% "hybrid",
                     boot_reps = if (!is.null(input$ef_reps) &&
                                     !is.na(input$ef_reps))
                       max(50, input$ef_reps) else NULL)
        } else if (identical(input$model_type, "mfrm")) {
          if (any(c(input$lp_person, input$lp_item, input$lp_score) == NONE) ||
              !length(input$lp_facets))
            stop("nominate the person, item, score, and at least one facet column")
          rasch_mfrm(df, person = input$lp_person, item = input$lp_item,
                     score = input$lp_score, facets = input$lp_facets,
                     n_groups = ng, adjust_N = adjN,
                     interaction = if (!is.null(input$lp_interaction) &&
                                       input$lp_interaction != NONE)
                       input$lp_interaction else NULL,
                     maxit = max(5, input$maxit %||% 60),
                     tol = max(1e-12, input$tol %||% 1e-8))
        } else {
          idc <- if (!is.null(input$id_col) && input$id_col != NONE) input$id_col else NULL
          fac <- if (length(input$factor_cols)) input$factor_cols else NULL
          its <- if (length(input$item_cols)) input$item_cols else NULL
          # multiple-choice key: uploaded CSV, or the demo key for the demo
          mc_key <- NULL
          if (!is.null(input$key_file)) {
            kf <- tryCatch(read.csv(input$key_file$datapath,
                                    stringsAsFactors = FALSE),
                           error = function(e) NULL)
            if (!is.null(kf) && (all(c("item", "key") %in% names(kf)) ||
                                 all(c("item", "option", "score") %in% names(kf))))
              mc_key <- kf
            else showNotification("Key CSV needs columns item,key or item,option,score - ignored.",
                                  type = "warning")
          } else if (identical(input$demo_choice, "dich")) {
            mc_key <- .demo_dich_key()
          }
          # anchors match by item name; rows for absent items are ignored
          anc <- anchors_in()
          if (!is.null(anc)) {
            cand <- if (!is.null(its)) its else setdiff(names(df), c(idc, fac))
            present <- as.character(anc$item) %in% cand
            if (!all(present))
              showNotification(sprintf("%d anchor row(s) ignored (items not in this dataset)",
                                       sum(!present)), type = "warning")
            anc <- anc[present, , drop = FALSE]
            if (!nrow(anc)) anc <- NULL
            # average anchoring: collapse to one mean-location anchor per item
            if (!is.null(anc) && identical(input$anchor_type, "average")) {
              mu <- tapply(anc$tau, as.character(anc$item), mean)
              anc <- data.frame(item = names(mu), k = NA, tau = as.numeric(mu))
            }
          }
          # principal-components (Andrich) threshold estimation; PCM only,
          # and not combinable with anchors
          pcc <- if (input$model_type %in% c("dich", "pcm") &&
                     identical(input$thr_mode, "pc"))
            as.integer(input$pc_rank %||% "4") else NULL
          if (!is.null(pcc) && !is.null(anc)) {
            showNotification("Anchors are ignored under principal-components threshold estimation.",
                             type = "warning", duration = 8)
            anc <- NULL
          }
          f0 <- rasch(df, model = if (identical(input$model_type, "rsm")) "RSM" else "PCM",
                      id = idc, factors = fac, items = its,
                      n_groups = ng, adjust_N = adjN, anchors = anc,
                      key = mc_key, pc_components = pcc,
                      maxit = max(5, input$maxit %||% 60),
                      tol = max(1e-12, input$tol %||% 1e-8))
          if (identical(input$model_type, "dich") && any(f0$m > 1L))
            showNotification("Some items have more than two categories; they were fitted with partial credit thresholds.",
                             type = "warning", duration = 10)
          f0
        }
      }, error = function(e) e)
    })
    if (inherits(fit, "error")) {
      showNotification(paste("Analysis failed:", conditionMessage(fit)),
                       type = "error", duration = NULL)
      return(NULL)
    }
    # routine handling notes are informational; only real problems warn
    if (length(fit$notes))
      showNotification(paste(fit$notes, collapse = "\n"), type = "message",
                       duration = 8)
    conv <- if (!is.null(fit$est)) fit$est$converged else fit$converged
    if (!isTRUE(conv))
      showNotification("Estimation did not converge; consider raising the maximum iterations or loosening the convergence criterion.",
                       type = "warning", duration = NULL)
    override_fit(NULL); override_desc(NULL)
    # paired-comparison results live on their own tab; the Rasch tabs
    # suspend while a BTL analysis is current
    if (inherits(fit, "rmt_btl")) {
      btl_fit(fit)
      try(nav_select("BTL", session = session), silent = TRUE)
      return(NULL)
    }
    btl_fit(NULL)
    try(nav_select("Summary", session = session), silent = TRUE)
    fit
  })
  btl_fit <- reactiveVal(NULL)
  fit <- reactive({
    f <- override_fit()
    if (is.null(f)) f <- analysis()
    req(f); f
  })

  observeEvent(fit(), {
    its <- fit()$items$item
    updateSelectInput(session, "dif_item", choices = its, selected = its[1])
    updateSelectizeInput(session, "subtest_items", choices = its, selected = character(0))
    fac <- names(fit()$factors)
    dif_choices <- if (length(fac)) c(fac, FACTORIAL) else NONE
    updateSelectInput(session, "dif_factor", choices = dif_choices,
                      selected = dif_choices[1])
    fs <- if (inherits(fit(), "rasch_mfrm")) fit()$facet_spec else NONE
    updateSelectInput(session, "facet_sel", choices = fs, selected = fs[1])
    fi <- if (inherits(fit(), "rasch_efrm"))
      unique(fit()$virtual_map$item) else NONE
    updateSelectInput(session, "frame_item", choices = fi, selected = fi[1])
    updateSelectizeInput(session, "dim_pos", choices = its, selected = character(0))
    updateSelectizeInput(session, "dim_neg", choices = its, selected = character(0))
    updateSelectInput(session, "dep_item", choices = its,
                      selected = its[min(2L, length(its))])
    updateSelectInput(session, "ind_item", choices = its, selected = its[1])
    updateSelectizeInput(session, "guess_anchors", choices = its,
                         selected = character(0))
    # results computed on request belong to the fit they came from
    lr_res(NULL); dep_res(NULL); spread_res(NULL); dm_res(NULL); guess_res(NULL)
  })

  observeEvent(input$make_subtest, {
    req(length(input$subtest_items) >= 2)
    res <- tryCatch(combine_items(fit(), list(input$subtest_items)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Subtest failed:", conditionMessage(res)), type = "error")
    } else {
      override_fit(res)
      override_desc(paste("subtest:", paste(input$subtest_items, collapse = " + ")))
      showNotification("Re-analysed with the subtest in place. Reset the override (Data page) or run again to return to the base fit.",
                       type = "message", duration = 8)
    }
  })
  output$subtest_status <- renderUI({
    if (is.null(override_desc())) return(NULL)
    p(class = "text-success small mt-2", paste("Active:", override_desc()))
  })

  observeEvent(input$make_split, {
    f <- fit()
    req(input$dif_item %in% f$items$item,
        !is.null(f$factors), input$dif_factor %in% names(f$factors))
    res <- tryCatch(split_items(f, input$dif_item, by = input$dif_factor),
                    error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Split failed:", conditionMessage(res)), type = "error")
    } else {
      override_fit(res)
      override_desc(sprintf("split: item %s by %s", input$dif_item, input$dif_factor))
      showNotification(
        sprintf("Re-analysed with %s split by %s. Reset the override (Data page) or run again to return to the base fit.",
                input$dif_item, input$dif_factor),
        type = "message", duration = 8)
    }
  })

  sel_item <- reactive({
    f <- fit()
    i <- input$items_tbl_rows_selected
    if (length(i)) f$items$item[i] else f$items$item[1]
  })

  # ------------------------------------------------------- plot plumbing --
  register_plot <- function(id, fun, w = 9, h = 6) {
    output[[id]] <- renderPlot(fun(), res = 96)
    for (fmt in c("png", "pdf")) local({
      fmt_ <- fmt
      output[[paste0(id, "_", fmt_)]] <- downloadHandler(
        filename = function() paste0("rmt_", id, ".", fmt_),
        content = function(file) {
          # 300 dpi PNG (and vector PDF) for publication
          if (fmt_ == "png") png(file, width = w, height = h, units = "in", res = 300)
          else pdf(file, width = w, height = h)
          fun(); dev.off()
        })
    })
  }
  register_table <- function(id, fun, dt_fun) {
    output[[id]] <- renderDT(dt_fun())
    output[[paste0(id, "_csv")]] <- downloadHandler(
      filename = function() paste0("rmt_", id, ".csv"),
      content = function(file) write.csv(fun(), file, row.names = FALSE))
  }
  num_dt <- function(d, digits = 3, ...) {
    num <- vapply(d, is.numeric, TRUE)
    datatable(d, rownames = FALSE, ...,
              options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
      formatRound(names(d)[num], digits)
  }

  # -------------------------------------------------------------- summary --
  output$vboxes <- renderUI({
    f <- fit()
    layout_column_wrap(width = "165px", fill = FALSE, class = "mb-3",
      value_box("Persons", nrow(f$X), theme = "primary"),
      value_box("Items", ncol(f$X), theme = "primary"),
      value_box("PSI", sprintf("%.3f", f$psi$PSI), theme = "success",
                p(class = "small mb-0", sprintf("%.3f no extremes", f$psi_noext$PSI))),
      value_box("Alpha", sprintf("%.3f", f$alpha$alpha), theme = "success",
                p(class = "small mb-0",
                  if (isFALSE(f$alpha$applicable))
                    sprintf("complete cases only (n = %d)", f$alpha$n)
                  else sprintf("n = %d complete", f$alpha$n))),
      value_box("Item-trait p", fmt_p(f$total_chisq_p),
                theme = if (f$total_chisq_p < 0.05) "danger" else "secondary"),
      value_box("Power of fit", f$power_of_fit, theme = "secondary")
    )
  })

  output$fit_summary <- renderPrint({
    f <- fit(); ss <- f$summary_stats
    cat(sprintf("Model: %s  |  Estimation: pairwise conditional ML (%s, %d iterations)\n",
                f$model, if (f$est$converged) "converged" else "NOT CONVERGED",
                f$est$iterations))
    cat(sprintf("Total item-trait chi-square: %.3f on %d df, p = %s  (%d class intervals)\n",
                f$total_chisq, f$total_df, fmt_p(f$total_chisq_p), f$n_groups))
    cat(sprintf("Item fit residual:   mean %6.2f  SD %5.2f  skew %5.2f  kurt %5.2f\n",
                f$item_fit_summary$mean, f$item_fit_summary$sd,
                f$item_fit_summary$skewness, f$item_fit_summary$kurtosis))
    cat(sprintf("Person fit residual: mean %6.2f  SD %5.2f  skew %5.2f  kurt %5.2f  (ideal 0, 1)\n",
                f$person_fit_summary$mean, f$person_fit_summary$sd,
                f$person_fit_summary$skewness, f$person_fit_summary$kurtosis))
    cat(sprintf("Fit-location correlation: items %.3f, persons %.3f; cell df factor %.3f\n",
                ss$cor_item_fit_location, ss$cor_person_fit_location,
                ss$df_factor))
    cat(sprintf("Items flagged misfitting (BH-adjusted): %d of %d\n",
                sum(f$items$misfit, na.rm = TRUE), nrow(f$items)))
    dis <- names(which(vapply(f$thresholds_diag, function(d)
      !d$ordered && length(d$thresholds) > 1, TRUE)))
    cat("Disordered thresholds:", if (length(dis)) paste(dis, collapse = ", ") else "none", "\n")
    if (length(f$notes)) cat("Notes:", paste(f$notes, collapse = "; "), "\n")
  })

  output$targeting <- renderPrint({
    f <- fit(); t <- f$targeting; ss <- f$summary_stats
    cat(sprintf("Person location: mean %6.3f  SD %.3f  skew %5.2f  kurt %5.2f\n",
                ss$person_location$mean, ss$person_location$sd,
                ss$person_location$skewness, ss$person_location$kurtosis))
    cat(sprintf("  without extremes: mean %.3f  SD %.3f\n",
                ss$person_location_noext$mean, ss$person_location_noext$sd))
    cat(sprintf("Item location:   mean %6.3f (constrained)  SD %.3f  skew %5.2f  kurt %5.2f\n",
                ss$item_location$mean, ss$item_location$sd,
                ss$item_location$skewness, ss$item_location$kurtosis))
    cat(sprintf("Threshold range:      %6.3f to %.3f\n",
                t$threshold_range[1], t$threshold_range[2]))
    cat(sprintf("Persons beyond thresholds: %.1f%% below, %.1f%% above\n",
                100 * t$prop_below, 100 * t$prop_above))
    cat(sprintf("\nPSI: %.3f (separation %.3f)\n", f$psi$PSI, f$psi$separation))
    cat(sprintf("PSI without extremes: %.3f (n = %d)\n", f$psi_noext$PSI, f$psi_noext$n))
    cat(sprintf("Item separation index: %.3f\n", f$isi$PSI))
    cat(sprintf("Cronbach alpha: %.3f%s\n", f$alpha$alpha,
                if (isFALSE(f$alpha$applicable))
                  sprintf(" - complete cases only (n = %d)", f$alpha$n)
                else sprintf(" (n = %d complete cases)", f$alpha$n)))
  })

  # score table with the chosen estimator and extreme-score treatment; the
  # geometric extrapolation needs >= 4 score points, so failures fall back
  score_dat <- reactive({
    f <- fit()
    if (is.null(f$score_table)) return(NULL)
    method <- input$st_method %||% "wle"
    ext <- input$st_extremes %||% "model"
    tryCatch(score_table(f, method = method, extremes = ext),
             error = function(e) {
               showNotification(paste0("Score table (", method, ", ", ext,
                                       "): ", conditionMessage(e),
                                       " - showing model extremes."),
                                type = "warning", duration = 8)
               score_table(f, method = method, extremes = "model")
             })
  })
  register_table("score_tbl", function() {
    f <- fit()
    if (!is.null(f$score_table)) score_dat() else f$score_curves
  }, function() {
    f <- fit()
    if (!is.null(f$score_table)) {
      d <- score_dat()
      d$extrapolated <- ifelse(d$extrapolated, "*", "")
      datatable(d, rownames = FALSE,
                options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
        formatRound(c("theta", "se"), 3) |>
        formatRound("cum_pct", 1)
    } else {
      validate(need(!is.null(f$score_curves),
                    "No score conversion available for this fit."))
      datatable(f$score_curves, rownames = FALSE,
                caption = "Raw scores are not sufficient under unequal frame units; per-group expected-score curves replace the score table.",
                options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
        formatRound(c("theta", "expected_score", "sem"), 3)
    }
  })

  # traditional (CTT) statistics
  ctt_res <- reactive(tryCatch(ctt_table(fit()), error = function(e) e))
  output$ctt_head <- renderUI({
    ct <- ctt_res()
    if (inherits(ct, "error"))
      return(p(class = "text-muted",
               paste("Traditional statistics unavailable:", conditionMessage(ct))))
    p(class = "small mb-2", HTML(sprintf(
      "Raw score mean <b>%.2f</b>, SD <b>%.2f</b>; alpha <b>%.3f</b>; SEM <b>%.2f</b> (one value for all persons) &mdash; complete cases n = %d.",
      ct$mean, ct$sd, ct$alpha, ct$sem, ct$n)))
  })
  output$ctt_tbl <- renderDT({
    ct <- ctt_res()
    validate(need(!inherits(ct, "error"),
                  "No traditional statistics for this fit."))
    num_dt(ct$table)
  })
  output$ctt_tbl_csv <- downloadHandler(
    filename = function() "rmt_ctt_tbl.csv",
    content = function(file) {
      ct <- ctt_res(); req(!inherits(ct, "error"))
      write.csv(ct$table, file, row.names = FALSE)
    })

  # likelihood-ratio test of PCM against the rating parameterisation; only
  # meaningful for a PCM fit whose items share a common maximum score > 1
  lr_res <- reactiveVal(NULL)
  output$lr_ui <- renderUI({
    f <- fit()
    if (!identical(f$model, "PCM") || length(unique(f$m)) != 1L ||
        max(f$m) < 2L) return(NULL)
    card(card_header("Likelihood-ratio test (PCM vs rating)"),
      card_body(
        p(class = "text-muted small",
          "Refits the same data with the rating (common threshold structure) parameterisation and compares the pairwise conditional log-likelihoods. A non-significant outcome supports adopting the simpler rating model; use the adjusted statistic for inference."),
        actionButton("run_lr", "Run likelihood-ratio test",
                     class = "btn-outline-primary"),
        verbatimTextOutput("lr_txt")))
  })
  observeEvent(input$run_lr, {
    f <- fit()
    r <- withProgress(message = "Refitting with the rating parameterisation…",
                      value = 0.4,
                      tryCatch(lr_test(f), error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("LR test failed:", conditionMessage(r)),
                       type = "error", duration = NULL)
    else lr_res(r)
  })
  output$lr_txt <- renderPrint({
    r <- lr_res()
    validate(need(!is.null(r), "Press the button to run the test."))
    cat(sprintf("Raw composite ChiSq %.3f on %d df, p = %s (conventional display; anticonservative)\n",
                r$chisq, r$df, fmt_p(r$p)))
    if (is.finite(r$chisq_adj))
      cat(sprintf("Adjusted ChiSq %.3f, p = %s (Kent 1982 calibration)\n",
                  r$chisq_adj, fmt_p(r$p_adj)))
    else cat("Adjusted (Kent 1982) statistic unavailable for this fit.\n")
    cat(sprintf("log-likelihood (pairwise composite): PCM %.3f, rating %.3f\n",
                r$loglik_pcm, r$loglik_rsm))
  })
  register_table("thr_tbl", function() {
    f <- fit(); d <- f$thresholds
    data.frame(item = f$items$item[d$item], threshold = d$k,
               tau = d$tau, se = d$se)
  }, function() {
    f <- fit(); d <- f$thresholds
    num_dt(data.frame(item = f$items$item[d$item], threshold = d$k,
                      tau = d$tau, se = d$se))
  })

  # ---------------------------------------------------------------- items --
  # any chi-square sample-size adjustment is applied inside the fit (the
  # adjust_N run control), so the table shows the fit's own statistics
  register_table("items_tbl", function() fit()$items, function() {
    d <- fit()$items
    d$misfit <- ifelse(d$misfit, "*", "")
    num_dt(d, selection = "single")
  })

  # per-class-interval breakdown of the selected item's chi-square
  chisq_res <- reactive(chisq_detail(fit(), sel_item()))
  output$chisq_caption <- renderUI({
    cd <- chisq_res()
    p(class = "small mb-2", HTML(sprintf(
      "<b>%s</b> (location %.3f): total ChiSq <b>%.3f</b> on %d df, p = %s; whole-sample mean = %.3f. Intervals with fewer than 2 responders carry no chi-square contribution.",
      cd$item, cd$location, cd$chisq, cd$df, fmt_p(cd$p), cd$ave)))
  })
  output$chisq_int_tbl <- renderDT({
    d <- chisq_res()$intervals
    d$used <- ifelse(d$used, "", "(excluded)")
    num_dt(d)
  })
  output$chisq_cat_tbl <- renderDT(num_dt(chisq_res()$categories))
  output$chisq_int_csv <- downloadHandler(
    filename = function() paste0("rmt_chisq_intervals_", sel_item(), ".csv"),
    content = function(file)
      write.csv(chisq_res()$intervals, file, row.names = FALSE))
  output$chisq_cat_csv <- downloadHandler(
    filename = function() paste0("rmt_chisq_categories_", sel_item(), ".csv"),
    content = function(file)
      write.csv(chisq_res()$categories, file, row.names = FALSE))

  # principal-components estimates (only for pc_components fits)
  output$pc_comp_ui <- renderUI({
    if (is.null(fit()$est$components)) return(NULL)
    tableCard("pc_tbl", "Principal components (location/spread/skewness/kurtosis)",
              "Andrich principal-components threshold estimates with standard errors; NA where an item's number of thresholds does not support the component.")
  })
  register_table("pc_tbl", function() fit()$est$components, function() {
    validate(need(!is.null(fit()$est$components),
                  "Run with principal-components threshold estimation to see the components."))
    num_dt(fit()$est$components)
  })

  register_plot("icc",  function() plot_icc(fit(), sel_item()))
  register_plot("ccc",  function()
    plot_ccc(fit(), sel_item(), observed = isTRUE(input$show_obs)))
  register_plot("tpc",  function()
    plot_threshold_prob(fit(), sel_item(), observed = isTRUE(input$show_obs)))
  register_plot("cfreq", function() plot_catfreq(fit(), sel_item()))
  mc_dat <- reactive({
    f <- fit()
    validate(need(!is.null(f$mc),
                  "Provide a multiple-choice key (CSV: item,key) to see distractor analysis."))
    distractor_analysis(f)
  })
  register_table("distractor_tbl", function() mc_dat(), function() {
    d <- mc_dat()
    d$keyed <- ifelse(d$keyed, "*", "")
    d$flag <- ifelse(d$flag, "MISKEY?", "")
    num_dt(d)
  })
  register_plot("distractor_plot", function() {
    f <- fit()
    validate(need(!is.null(f$mc),
                  "Provide a multiple-choice key (CSV: item,key) to see option curves."))
    it <- if (sel_item() %in% colnames(f$mc$raw)) sel_item() else
      colnames(f$mc$raw)[1]
    plot_distractors(f, it)
  })

  # polytomous option-scoring proposal (Andrich & Styles 2011)
  rescore_res <- reactiveVal(NULL)
  observeEvent(input$rescore_go, {
    f <- fit()
    if (is.null(f$mc)) {
      showNotification("Provide a multiple-choice key first.", type = "warning")
      return()
    }
    res <- tryCatch(distractor_rescore(f,
                                       min_n = max(2, input$rescore_min_n %||% 20),
                                       z = max(0.1, input$rescore_z %||% 1.96)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Rescore proposal failed:", conditionMessage(res)),
                       type = "error")
      return()
    }
    rescore_res(res)
    n_cred <- sum(res$option_scores$score > 0 &
                  !res$evidence$keyed[match(paste(res$option_scores$item,
                                                  res$option_scores$option),
                                            paste(res$evidence$item,
                                                  res$evidence$option))])
    showNotification(sprintf("%d distractor(s) proposed for partial credit.",
                             n_cred), type = "message")
  })
  output$rescore_tbl <- DT::renderDT({
    res <- rescore_res()
    validate(need(!is.null(res), "Run the proposal to see the evidence table."))
    d <- res$evidence
    d$keyed <- ifelse(d$keyed, "*", "")
    num_dt(d)
  })
  output$dl_rescore <- downloadHandler(
    filename = function() "option_scores.csv",
    content = function(file) {
      res <- rescore_res()
      if (is.null(res)) stop("run the proposal first")
      write.csv(res$option_scores, file, row.names = FALSE)
    })

  # -------------------------------------------------------------- persons --
  register_table("person_tbl", function() fit()$person, function() {
    d <- fit()$person
    datatable(d, rownames = FALSE, filter = "top", selection = "single",
              options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
      formatRound(names(d)[vapply(d, is.numeric, TRUE) &
                           !names(d) %in% c("raw", "max_raw", "n_items", "class_interval")], 3)
  })
  register_plot("pcc", function() {
    f <- fit()
    n <- input$person_tbl_rows_selected
    n <- if (length(n)) n[1] else 1L
    plot_pcc(f, n)
  })
  register_plot("rdist", function()
    plot_resid_dist(fit(), what = input$rd_what %||% "persons",
                    statistic = input$rd_stat %||% "fit_resid"))
  register_plot("pfit",  function() plot_person_fit(fit()))
  register_plot("pim_p", function() plot_pimap(fit()))

  # ------------------------------------------------------------ test plots --
  register_plot("thrmap", function() plot_threshold_map(fit()), h = 7)
  register_plot("imap",   function() plot_item_map(fit()))
  register_plot("tcc",    function() plot_tcc(fit()))
  register_plot("tif",    function() plot_tif(fit()))
  register_plot("guttman", function() plot_guttman(fit()), h = 7)

  # ------------------------------------------------------------------ DIF --
  # (the FACTORIAL constant is defined at the top of the file: observers
  # created before this section reference it)
  dif_alpha <- reactive({
    a <- input$dif_alpha
    if (is.null(a) || is.na(a) || a <= 0 || a >= 1) 0.05 else a
  })
  dif_res <- reactive({
    f <- fit(); req(!is.null(f$factors), length(names(f$factors)) > 0)
    dif_anova(f, p_adjust = input$dif_padj %||% "BH", alpha = dif_alpha())
  })
  dif_fact <- reactive({
    f <- fit(); req(!is.null(f$factors), length(names(f$factors)) >= 1)
    dif_anova_factorial(f, p_adjust = input$dif_padj %||% "BH",
                        alpha = dif_alpha(),
                        effects = input$dif_effects %||% "factorial")
  })
  register_table("dif_tbl", function() {
    if (identical(input$dif_factor, FACTORIAL)) dif_fact()$terms else dif_res()
  }, function() {
    if (identical(input$dif_factor, FACTORIAL)) {
      d <- dif_fact()$terms
      d$significant <- ifelse(d$significant, "*", "")
      d$superseded <- ifelse(d$superseded, "(superseded)", "")
      num_dt(d) |> formatRound(c("p", "p_adj"), 3)
    } else {
      d <- dif_res()
      if (!is.null(input$dif_factor) && input$dif_factor %in% d$factor)
        d <- d[d$factor == input$dif_factor, ]
      d$uniform_DIF <- ifelse(d$uniform_DIF, "*", "")
      d$nonuniform_DIF <- ifelse(d$nonuniform_DIF, "*", "")
      num_dt(d) |> formatRound(c("p_uniform", "p_nonuniform",
                                "p_uniform_adj", "p_nonuniform_adj"), 3)
    }
  })
  register_table("dif_tukey_tbl", function() dif_fact()$tukey, function() {
    validate(need(identical(input$dif_factor, FACTORIAL),
                  "Choose the factorial option in the sidebar to see Tukey HSD comparisons."))
    tk <- dif_fact()$tukey
    if (!nrow(tk))
      return(datatable(data.frame(note = "no significant group terms to compare"),
                       rownames = FALSE, options = list(dom = "t")))
    num_dt(tk) |> formatRound("p_tukey", 3)
  })
  register_plot("dif_icc", function() {
    f <- fit()
    req(input$dif_item %in% f$items$item,
        !is.null(f$factors), input$dif_factor %in% names(f$factors))
    plot_icc(f, input$dif_item, group = input$dif_factor)
  })

  # DIF magnitude in logits: single factor -> the selected item and factor;
  # factorial -> sizes for every significant, non-superseded group term
  dif_size_res <- reactiveVal(NULL)
  observeEvent(input$dif_size_go, {
    f <- fit()
    flg <- max(0.05, input$dif_size_flag %||% 0.5)
    mn <- max(2, input$dif_size_minn %||% 20)
    res <- tryCatch({
      if (identical(input$dif_factor, FACTORIAL)) {
        fa <- dif_anova_factorial(f, sizes = TRUE,
                                  p_adjust = input$dif_padj %||% "BH",
                                  alpha = dif_alpha(),
                                  effects = input$dif_effects %||% "factorial")
        sz <- fa$sizes
        if (is.null(sz) || !nrow(sz))
          stop("no significant, non-superseded group terms to size")
        sz$practical <- abs(sz$difference) >= flg
        sz
      } else {
        req(input$dif_item %in% f$items$item,
            input$dif_factor %in% names(f$factors))
        ds <- dif_size(f, input$dif_item, by = input$dif_factor,
                       flag_logits = flg, min_n = mn)
        cbind(item = ds$item, term = ds$by, ds$pairs)
      }
    }, error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("DIF size:", conditionMessage(res)),
                       type = "warning")
      dif_size_res(NULL)
    } else dif_size_res(res)
  })
  output$dif_size_tbl <- DT::renderDT({
    d <- dif_size_res()
    validate(need(!is.null(d), "Compute DIF sizes to see the logit-scale comparisons."))
    d$significant <- ifelse(d$significant, "*", "")
    d$practical <- ifelse(d$practical, "PRACTICAL", "")
    num_dt(d)
  })
  output$dl_dif_size <- downloadHandler(
    filename = function() "dif_sizes.csv",
    content = function(file) {
      d <- dif_size_res()
      if (is.null(d)) stop("compute DIF sizes first")
      write.csv(d, file, row.names = FALSE)
    })

  # ------------------------------------------------------------------- BTL --
  bfit <- reactive({
    validate(need(!is.null(btl_fit()),
                  "Choose the paired-comparisons (BTL) model on the Data tab and run."))
    btl_fit()
  })
  output$btl_boxes <- renderUI({
    f <- bfit()
    layout_column_wrap(width = "165px", fill = FALSE, class = "mb-3",
      value_box("Objects", nrow(f$objects), theme = "primary"),
      value_box("Comparisons", sprintf("%.0f", f$n_comparisons), theme = "primary"),
      if (!is.null(f$judges))
        value_box("Judges", nrow(f$judges), theme = "primary"),
      value_box("Object separation", sprintf("%.3f", f$osi$PSI),
                theme = "success"),
      value_box("Pairwise fit p", fmt_p(f$total_p),
                theme = if (f$total_p < 0.05) "danger" else "success"))
  })
  output$btl_summary <- renderText({
    f <- bfit()
    paste0(sprintf("Conditional ML: %s in %d iterations; sandwich SEs%s.\n",
                   if (f$converged) "converged" else "NOT converged",
                   f$iterations,
                   if (f$clustered) " clustered by judge" else ""),
           sprintf("Pairwise chi-square %.2f on %d df, p = %s.\n",
                   f$total_chisq, f$total_df, fmt_p(f$total_p)),
           sprintf("Log-likelihood %.3f.", f$loglik),
           if (length(f$notes)) paste0("\nNotes: ",
                                       paste(f$notes, collapse = "; ")) else "")
  })
  register_table("btl_obj_tbl", function() bfit()$objects,
                 function() num_dt(bfit()$objects))
  register_table("btl_pairs_tbl", function() bfit()$pairs,
                 function() num_dt(bfit()$pairs))
  register_table("btl_judges_tbl", function() {
    validate(need(!is.null(bfit()$judges), "No judge column was nominated."))
    bfit()$judges
  }, function() {
    validate(need(!is.null(bfit()$judges), "No judge column was nominated."))
    d <- bfit()$judges
    d$misfit <- ifelse(!is.na(d$fit_resid) & abs(d$fit_resid) > 2.5, "*", "")
    num_dt(d)
  })
  register_plot("btl_plot", function() plot_btl(bfit()))

  # ---------------------------------------------------------------- facets --
  facet_dat <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm"),
                  "Run a long-format (many-facet) analysis to see facet results."))
    req(input$facet_sel %in% f$facet_spec)
    f$facet_effects[[input$facet_sel]]
  })
  register_table("facet_tbl", function() facet_dat(), function() {
    d <- facet_dat()
    datatable(d, rownames = FALSE,
              options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
      formatRound(setdiff(names(d)[vapply(d, is.numeric, TRUE)], "n"), 3)
  })
  register_plot("facet_plot", function() {
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm"),
                  "Run a long-format (many-facet) analysis to see facet results."))
    req(input$facet_sel %in% f$facet_spec)
    plot_facets(f, input$facet_sel)
  })
  facet_int <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm") && !is.null(f$interaction),
                  "Run a long-format analysis with an item-by-facet interaction selected."))
    f$interaction_effects
  })
  register_table("facet_int_tbl", function() facet_int(), function() {
    d <- facet_int()
    d$significant <- ifelse(abs(d$gamma) > 1.96 * d$se, "*", "")
    num_dt(d)
  })

  # ----------------------------------------------------------------- equating --
  eq_ref <- reactive({
    req(input$eq_file)
    a <- tryCatch(read.csv(input$eq_file$datapath, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    validate(need(!is.null(a) && all(c("item", "location") %in% names(a)),
                  "The reference CSV needs columns item, location (and ideally se)."))
    a
  })
  # reference: an uploaded calibration CSV, or a fit kept on the Compare page
  eq_reference <- reactive({
    if (identical(input$eq_source, "kept")) {
      k <- kept_fits()
      validate(need(length(k) >= 1,
                    "Keep a fit on the Compare page to use it as the equating reference."))
      validate(need(!is.null(input$eq_kept) && input$eq_kept %in% names(k),
                    "Choose a kept fit in the sidebar."))
      k[[input$eq_kept]]
    } else {
      validate(need(!is.null(input$eq_file),
                    "Upload a reference calibration (item, location, se) to equate against."))
      eq_ref()
    }
  })
  eq_res <- reactive(equate_tests(fit(), eq_reference(), shift = input$eq_shift))
  register_table("eq_tbl", function() eq_res()$table, function() {
    eq <- eq_res()
    d <- eq$table
    d$drift <- ifelse(d$drift, "*", "")
    num_dt(d) |> formatRound(c("p", "p_adj"), 3)
  })
  register_plot("eq_plot", function()
    plot_equate(fit(), eq_reference(), shift = input$eq_shift))
  output$dl_anchors <- downloadHandler(
    filename = function() format(Sys.time(), "rmt_anchors_%Y%m%d_%H%M.csv"),
    content = function(file) {
      f <- fit()
      thr <- f$thresholds
      write.csv(data.frame(item = f$items$item[thr$item], k = thr$k,
                           tau = thr$tau), file, row.names = FALSE)
    })

  output$dl_calib <- downloadHandler(
    filename = function() format(Sys.time(), "rmt_calibration_%Y%m%d_%H%M.csv"),
    content = function(file) {
      f <- fit()
      write.csv(data.frame(item = f$items$item, location = f$items$location,
                           se = f$items$se), file, row.names = FALSE)
    })

  # ---------------------------------------------------------------- frames --
  efrm_fit <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_efrm"),
                  "Run a frames (extended frame of reference) analysis to see results here."))
    f
  })
  register_table("frame_tbl", function() efrm_fit()$frames,
                 function() num_dt(efrm_fit()$frames))
  register_table("phi_tbl", function() efrm_fit()$phi_table,
                 function() num_dt(efrm_fit()$phi_table))
  register_table("alpha_tbl", function()
    merge(efrm_fit()$alpha_table, efrm_fit()$set_table[, c("set", "mu", "n_items")],
          by = "set"),
    function() num_dt(merge(efrm_fit()$alpha_table,
                            efrm_fit()$set_table[, c("set", "mu", "n_items")],
                            by = "set")))
  register_plot("frame_plot", function() plot_frames(efrm_fit()))
  register_plot("frame_icc", function() {
    f <- efrm_fit()
    req(input$frame_item %in% f$virtual_map$item)
    plot_icc_frames(f, input$frame_item)
  })
  output$efrm_cmp <- renderPrint({
    f <- efrm_fit(); cmp <- f$efrm_vs_rasch
    cat(sprintf("Pairwise conditional log-likelihood: frames model %.3f, equal units %.3f\n",
                cmp$ll_efrm, cmp$ll_equal))
    cat(sprintf("2 x improvement: %.3f with %d extra unit parameter(s)\n",
                cmp$two_delta_ll, cmp$extra_parameters))
    cat("(composite likelihood: descriptive; informative for ",
        cmp$informative_for, ")\n", sep = "")
    if (!is.null(cmp$unit_tests)) {
      cat("\nWald tests of the units (H0: unit = 1):\n")
      print(cmp$unit_tests, digits = 3, row.names = FALSE)
    }
    cat(sprintf("\nItem fit residual SD under the frames model: %.3f\n",
                f$item_fit_summary$sd))
  })

  # --------------------------------------------------------- dimensionality --
  dim_subsets <- reactiveVal(NULL)
  observeEvent(input$run, dim_subsets(NULL), priority = 9)
  observeEvent(input$dim_apply, {
    if (length(input$dim_pos) >= 2 && length(input$dim_neg) >= 2) {
      if (length(intersect(input$dim_pos, input$dim_neg))) {
        showNotification("The two subsets must be disjoint.", type = "error")
      } else dim_subsets(list(pos = input$dim_pos, neg = input$dim_neg))
    } else if (!length(input$dim_pos) && !length(input$dim_neg)) {
      dim_subsets(NULL)
      showNotification("Reset to the first-contrast split.", type = "message")
    } else {
      showNotification("Nominate at least two items in each subset (or leave both empty).",
                       type = "warning")
    }
  })
  dim_res <- reactive({
    s <- dim_subsets()
    if (is.null(s)) dimensionality_test(fit())
    else dimensionality_test(fit(), items_positive = s$pos,
                             items_negative = s$neg)
  })
  output$dim_txt <- renderPrint({
    dt <- dim_res()
    if (!is.null(dt$note)) { cat(dt$note); return(invisible()) }
    cat(sprintf("Item split: %s\n", dt$split))
    cat(sprintf("First residual eigenvalue: %.3f\n", dt$first_eigenvalue))
    cat(sprintf("Significant person t-tests: %.1f%%  (exact 95%% CI %.1f%% to %.1f%%, n = %d)\n",
                100 * dt$prop_significant, 100 * dt$ci[1], 100 * dt$ci[2], dt$n))
    cat(sprintf("Persons excluded (extreme on a subset): %d\n", dt$n_excluded_extreme))
    cat(sprintf("Verdict: %s\n", if (dt$multidimensional)
      "lower CI exceeds 5% - unidimensionality is questionable"
      else "consistent with unidimensionality"))
    if (!is.null(dt$paired_t))
      cat(sprintf("Paired t-test of subset means: mean difference %.3f, t = %.2f (df %.0f), p = %s\n",
                  dt$paired_t$mean_difference, dt$paired_t$t,
                  dt$paired_t$df, fmt_p(dt$paired_t$p)))
    cat("\nSubset A items:\n ", paste(dt$items_positive, collapse = ", "), "\n")
    cat("Subset B items:\n ", paste(dt$items_negative, collapse = ", "), "\n")
  })

  # magnitude of multidimensionality (Andrich 2016): needs every item in a
  # subset; the PC1 split satisfies this by construction, manual subsets may not
  dm_res <- reactiveVal(NULL)
  observeEvent(input$dm_run, {
    f <- fit()
    s <- dim_subsets()
    if (is.null(s)) {
      dr <- dim_res()
      if (!is.null(dr$note)) {
        showNotification(paste("No usable subsets:", dr$note), type = "warning")
        return()
      }
      s <- list(pos = dr$items_positive, neg = dr$items_negative)
    }
    allit <- c(s$pos, s$neg)
    if (!setequal(allit, f$items$item) || anyDuplicated(allit) > 0) {
      left <- setdiff(f$items$item, allit)
      showNotification(paste0(
        "The magnitude estimate needs every item assigned to exactly one subset. ",
        if (length(left)) paste0("Unassigned: ", paste(left, collapse = ", "), ". ") else "",
        "Adjust subsets A and B (or leave both empty for the PC1 split) and re-run the t-test first."),
        type = "warning", duration = 10)
      return()
    }
    r <- withProgress(message = "Subtest re-analysis…", value = 0.4,
                      tryCatch(dimensionality_magnitude(f, list(s$pos, s$neg)),
                               error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Magnitude estimate failed:", conditionMessage(r)),
                       type = "error", duration = NULL)
    else dm_res(r)
  })
  output$dm_tbl <- renderDT({
    r <- dm_res()
    validate(need(!is.null(r),
                  "Press the button; the current subsets (manual, or the PC1 split) are combined into super-items and the two reliability calculations compared."))
    num_dt(r$table)
  })
  output$dm_tbl_csv <- downloadHandler(
    filename = function() "rmt_dimensionality_magnitude.csv",
    content = function(file) {
      r <- dm_res(); req(!is.null(r))
      write.csv(r$table, file, row.names = FALSE)
    })
  register_plot("scree", function() plot_scree(fit()))
  register_table("loadings_tbl", function() residual_pca(fit())$loadings_matrix,
                 function() num_dt(residual_pca(fit())$loadings_matrix))
  register_table("eigen_tbl", function() residual_pca(fit())$eigen_table,
                 function() num_dt(residual_pca(fit())$eigen_table))
  register_plot("pca_plot", function() plot_pca(fit()))

  # -------------------------------------------------------- local dependence --
  register_plot("rcor", function() plot_resid_cor(fit()), w = 8, h = 8)
  ld_res <- reactive({
    fl <- input$ld_flag
    if (is.null(fl) || is.na(fl) || fl <= 0) fl <- 0.2
    residual_correlations(fit(), flag = fl)
  })
  register_table("rpairs_tbl", function() ld_res()$flagged,
                 function() {
    fl <- ld_res()$flagged
    if (!nrow(fl)) fl <- data.frame(note = "no item pairs exceed the flag threshold")
    num_dt(fl)
  })

  # response dependence magnitude (Andrich & Kreiner resolved-item refit)
  dep_res <- reactiveVal(NULL)
  observeEvent(input$run_dep, {
    f <- fit()
    req(input$dep_item %in% f$items$item, input$ind_item %in% f$items$item)
    r <- withProgress(message = "Resolving and re-analysing…", value = 0.4,
                      tryCatch(dependence_magnitude(f,
                                                    dependent = input$dep_item,
                                                    independent = input$ind_item),
                               error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Dependence estimate failed:", conditionMessage(r)),
                       type = "error", duration = NULL)
    else dep_res(r)
  })
  output$dep_txt <- renderPrint({
    r <- dep_res()
    validate(need(!is.null(r),
                  "Choose the dependent and independent items and press the button."))
    cat(sprintf("Dependence of %s on %s: d = %.3f logits (se %.3f), z = %.2f, p = %s\n",
                r$dependent, r$independent, r$d, r$se, r$z, fmt_p(r$p)))
  })
  output$dep_tbl <- renderDT({
    r <- dep_res()
    validate(need(!is.null(r), ""))
    num_dt(r$thresholds)
  })
  output$dep_tbl_csv <- downloadHandler(
    filename = function() "rmt_dependence_thresholds.csv",
    content = function(file) {
      r <- dep_res(); req(!is.null(r))
      write.csv(r$thresholds, file, row.names = FALSE)
    })

  # spread test against Andrich's least upper bounds
  spread_res <- reactiveVal(NULL)
  observeEvent(input$run_spread, {
    r <- withProgress(message = "Principal-components refit…", value = 0.4,
                      tryCatch(spread_test(fit()), error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Spread test failed:", conditionMessage(r)),
                       type = "error", duration = NULL)
    else spread_res(r)
  })
  output$spread_tbl <- renderDT({
    r <- spread_res()
    validate(need(!is.null(r),
                  "Press the button to run the spread test (needs at least one polytomous item, e.g. after combining a subtest)."))
    d <- r
    d$dependent <- ifelse(d$dependent, "*", "")
    num_dt(d)
  })
  output$spread_tbl_csv <- downloadHandler(
    filename = function() "rmt_spread_test.csv",
    content = function(file) {
      r <- spread_res(); req(!is.null(r))
      write.csv(r, file, row.names = FALSE)
    })

  # ---------------------------------------------------------------- guessing --
  guess_res <- reactiveVal(NULL)
  observeEvent(input$run_guess, {
    f <- fit()
    if (max(f$m) != 1L) {
      showNotification("Tailored analysis applies to dichotomous analyses only.",
                       type = "warning")
      return()
    }
    ch <- input$guess_chance
    if (is.null(ch) || is.na(ch) || ch <= 0 || ch >= 1) ch <- 0.25
    anc <- if (length(input$guess_anchors)) input$guess_anchors else NULL
    r <- withProgress(message = "Tailored analysis (three re-analyses)…",
                      value = 0.3,
                      tryCatch(tailored_analysis(f, chance = ch,
                                                 anchor_items = anc),
                               error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Tailored analysis failed:", conditionMessage(r)),
                       type = "error", duration = NULL)
    else guess_res(r)
  })
  output$guess_txt <- renderPrint({
    validate(need(max(fit()$m) == 1L,
                  "Tailored analysis applies to dichotomous (multiple-choice) data only."))
    r <- guess_res()
    validate(need(!is.null(r),
                  "Set the chance level in the sidebar and press the button."))
    cat(sprintf("Responses set to missing (P below chance %.2f): %d\n",
                r$chance, r$n_removed))
    cat("Anchor items for the common origin:",
        paste(r$anchor_items, collapse = ", "), "\n")
    up <- sum(r$table$z > 1.96, na.rm = TRUE)
    cat(sprintf("Items significantly harder in the tailored analysis (z > 1.96): %d of %d\n",
                up, nrow(r$table)))
    cat("Verdict:", if (up > 0) "guessing indicated"
        else "no guessing signature", "\n")
  })
  register_table("guess_tbl", function() {
    r <- guess_res(); req(!is.null(r)); r$table
  }, function() {
    validate(need(max(fit()$m) == 1L,
                  "Tailored analysis applies to dichotomous (multiple-choice) data only."))
    r <- guess_res()
    validate(need(!is.null(r), "Run the tailored analysis to see the comparison."))
    num_dt(r$table)
  })
  register_plot("guess_plot", function() {
    validate(need(max(fit()$m) == 1L,
                  "Tailored analysis applies to dichotomous (multiple-choice) data only."))
    r <- guess_res()
    validate(need(!is.null(r), "Run the tailored analysis to see the equating plot."))
    plot_equate(r$tailored, r$origin_equated, shift = "none")
  })

  # ---------------------------------------------------------------- compare --
  kept_fits <- reactiveVal(list())
  observeEvent(input$keep_fit, {
    f <- fit()
    k <- kept_fits()
    lab <- sprintf("%d_%s", length(k) + 1L, f$model)
    k[[lab]] <- f
    kept_fits(k)
    updateSelectInput(session, "cmp_ref", choices = names(k),
                      selected = if (!is.null(input$cmp_ref) &&
                                     input$cmp_ref %in% names(k))
                        input$cmp_ref else names(k)[1])
    updateSelectInput(session, "eq_kept", choices = names(k),
                      selected = if (!is.null(input$eq_kept) &&
                                     input$eq_kept %in% names(k))
                        input$eq_kept else names(k)[length(k)])
    showNotification(sprintf("Kept '%s' (%d fit(s) held).", lab, length(k)),
                     type = "message", duration = 5)
  })
  observeEvent(input$clear_fits, {
    kept_fits(list())
    updateSelectInput(session, "cmp_ref", choices = NONE, selected = NONE)
    updateSelectInput(session, "eq_kept", choices = NONE, selected = NONE)
    showNotification("Cleared kept fits.", type = "message", duration = 4)
  })
  cmp_res <- reactive({
    k <- kept_fits()
    validate(need(length(k) >= 2,
                  "Keep at least two fits (run, keep, change settings, run, keep) to compare."))
    ref <- if (!is.null(input$cmp_ref) && input$cmp_ref %in% names(k))
      input$cmp_ref else 1
    as.data.frame(do.call(compare_fits, c(k, list(reference = ref))))
  })
  register_table("cmp_tbl", function() cmp_res(), function() {
    d <- cmp_res()
    d$same_data <- ifelse(d$same_data, "yes", "no")
    num_dt(d)
  })

  # ----------------------------------------------------------------- export --
  output$dl_zip <- downloadHandler(
    filename = function() format(Sys.time(), "rmt_results_%Y%m%d_%H%M.zip"),
    content = function(file) {
      f <- fit()
      tmp <- file.path(tempdir(), paste0("rmt_", as.integer(Sys.time())))
      withProgress(message = "Writing all tables and plots…", value = 0.4, {
        save_outputs(f, tmp,
                     formats = if (length(input$exp_formats)) input$exp_formats else "png",
                     item_plots = isTRUE(input$exp_items))
      })
      owd <- setwd(tmp); on.exit(setwd(owd), add = TRUE)
      utils::zip(zipfile = file, files = list.files(".", recursive = TRUE),
                 flags = "-r9Xq")
    })
}

shinyApp(ui, server)
