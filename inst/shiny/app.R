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
library(bsicons)

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

# value-box guard: NULL, NA, NaN, and Inf must never reach a conditional or
# a sprintf; such values display as an em dash on a neutral theme
finite1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

# helpers for the "R code for this analysis" disclosure
qstr <- function(x) paste0('"', x, '"')
qvec <- function(x)
  if (length(x) == 1L) qstr(x) else
    paste0("c(", paste(qstr(x), collapse = ", "), ")")

theme <- bs_theme(
  version = 5, preset = "shiny",
  bg = "#f8fafc", fg = "#0f172a",
  primary = "#2563eb", secondary = "#64748b",
  success = "#0f766e", danger = "#dc2626", warning = "#d97706",
  base_font = font_google("Inter"),
  heading_font = font_google("Inter"),
  code_font = font_google("JetBrains Mono"),
  "navbar-bg" = "#0f172a",
  "headings-font-weight" = "600",
  "font-size-base" = "0.925rem"
)

css <- HTML("
  .navbar-brand { font-weight: 700; letter-spacing: .02em; }
  .card-header { font-weight: 600; }
  .value-box-title { font-size: .72rem; text-transform: uppercase; letter-spacing: .04em; white-space: nowrap; }
  .value-box-value { font-size: 1.45rem; }
  pre, .shiny-text-output { white-space: pre-wrap; font-size: .82rem; }
  .btn-xs { padding: .1rem .5rem; font-size: .75rem; }
  .form-label { font-weight: 600; font-size: .85rem; }
  /* APA-style tables: tabular numerals, no vertical rules, no zebra,
     a strong rule under the header row */
  table.dataTable { font-size: .85rem; }
  table.dataTable td { font-variant-numeric: tabular-nums; }
  table.dataTable thead th {
    font-weight: 600;
    border-bottom: 2px solid var(--bs-emphasis-color) !important;
  }
  table.dataTable th, table.dataTable td {
    border-left: none !important; border-right: none !important;
  }
  table.dataTable.table-striped > tbody > tr:nth-of-type(odd) > * {
    --bs-table-accent-bg: transparent;
  }
  .table-note { color: var(--bs-secondary-color); font-size: .8rem; }
  .empty-state { text-align: center; padding: 3rem 1rem; color: var(--bs-secondary-color); }
  .shiny-output-error-validation {
    text-align: center; padding: 2rem 1rem; color: var(--bs-secondary-color);
  }
  .nav-status .badge { font-weight: 500; }
")

# Card with a plot and PNG/PDF download buttons in the header. The body is
# non-fillable so flex sizing can never compress the fixed-height plot
# (the cause of the squashed plots), and a percentage height is avoided
# because it races the layout and renders a zero-height device.
# data-bs-theme is pinned to light because base plots draw on white.
plotCard <- function(id, title, height = "560px") {
  card(
    full_screen = TRUE,
    `data-bs-theme` = "light",
    card_header(div(class = "d-flex justify-content-between align-items-center",
      span(title),
      div(class = "btn-group",
          downloadButton(paste0(id, "_png"), "PNG", class = "btn-outline-secondary btn-xs"),
          downloadButton(paste0(id, "_pdf"), "PDF", class = "btn-outline-secondary btn-xs")))),
    card_body(plotOutput(id, height = height), padding = 8, fillable = FALSE)
  )
}

# `info` adds a header tooltip defining the key statistic; `footer` takes a
# small UI slot rendered under the table (dynamic interpretation notes)
tableCard <- function(id, title, note = NULL, info = NULL, footer = NULL,
                      controls = NULL) {
  card(
    full_screen = TRUE,
    card_header(div(class = "d-flex justify-content-between align-items-center",
      span(title,
           if (!is.null(info))
             tooltip(bs_icon("info-circle", class = "ms-1 text-secondary"), info)),
      div(class = "d-flex align-items-center gap-3",
          controls,
          downloadButton(paste0(id, "_csv"), "CSV",
                         class = "btn-outline-secondary btn-xs")))),
    card_body(if (!is.null(note)) p(class = "text-muted small mb-2", note),
              DTOutput(id),
              if (!is.null(footer)) div(class = "table-note mt-2", footer),
              padding = 12)
  )
}


# compact header switch revealing every column of a curated table
cols_switch <- function(id)
  div(class = "small text-secondary",
      input_switch(id, "All columns", value = FALSE))

# card header with an info-circle tooltip (for non-table cards)
info_header <- function(title, info)
  card_header(span(title,
    tooltip(bs_icon("info-circle", class = "ms-1 text-secondary"), info)))

# ---------------------------------------------------------------------------
# Panels are built as objects and assembled into the workflow-ordered navbar
# (with Structure / Invariance / More menus) at the end of the UI section.
# ----------------------------------------------------------------- DATA --
panel_data <- nav_panel("Data", icon = bs_icon("database"),
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
        accordion(
          id = "run_settings", multiple = TRUE,
          open = c("Data roles", "Model"),
          accordion_panel("Model", icon = bs_icon("diagram-2"),
            radioButtons("model_type", NULL,
                         c("Dichotomous" = "dich",
                           "Partial credit (PCM)" = "pcm",
                           "Rating scale (RSM)" = "rsm",
                           "Many-facet (MFRM)" = "mfrm",
                           "Extended frames (EFRM)" = "efrm",
                           "Paired comparisons (BTL)" = "btl"))),
          accordion_panel("Data roles", icon = bs_icon("table"),
            conditionalPanel("['dich','pcm','rsm'].indexOf(input.model_type) > -1",
              selectInput("id_col", "ID variable", NONE),
              selectizeInput("factor_cols", "Person factors (DIF groups)", NULL,
                             multiple = TRUE,
                             options = list(placeholder = "none selected")),
              selectizeInput("item_cols", "Item columns", NULL, multiple = TRUE,
                             options = list(placeholder = "all remaining columns"))
            ),
            conditionalPanel("input.model_type == 'efrm'",
              selectInput("ef_id", "ID variable", NONE),
              selectInput("ef_group", "Person group column ((none) = single group; units differ by item set only)",
                          NONE),
              selectizeInput("ef_items", "Item columns", NULL, multiple = TRUE,
                             options = list(placeholder = "all remaining columns")),
              fileInput("ef_sets", "Item-set map (CSV: item,set)",
                        accept = ".csv", placeholder = "optional"),
              checkboxInput("ef_prefix", "Infer sets from item-name prefix", TRUE),
              p(class = "text-muted small",
                "Each item-set by group cell is a frame with its own unit. Group units come from the person-free pairwise comparisons; set units from persons common to the sets.")
            ),
            conditionalPanel("input.model_type == 'mfrm'",
              h6("One row per response"),
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
              h6("One comparison per row"),
              selectInput("bt_a", "Object A column", NONE),
              selectInput("bt_b", "Object B column", NONE),
              selectInput("bt_win", "Winner column", NONE),
              selectInput("bt_judge", "Judge column (optional)", NONE),
              selectInput("bt_count", "Count column (optional)", NONE),
              radioButtons("bt_ties", "Ties",
                           c("Drop" = "drop", "Half a win each" = "half")),
              p(class = "text-muted small",
                "The Bradley-Terry-Luce model: the conditional (person-free) form of the dichotomous Rasch model, estimated by the same conventions. A judge column enables the judge fit table and clusters the standard errors by judge. Results appear on the BTL tab.")
            )),
          accordion_panel("Estimation options", icon = bs_icon("gear"),
            checkboxInput("ng_auto", "Automatic class intervals (at least 50 per interval)", TRUE),
            conditionalPanel("!input.ng_auto",
              sliderInput("ng", "Class intervals", min = 2, max = 16, value = 8)),
            numericInput("run_adjN",
                         "Adjust chi-square to sample size N (optional)",
                         value = NA, min = 50),
            numericInput("maxit", "Maximum iterations", value = 60, min = 5, step = 5),
            numericInput("tol", "Convergence criterion", value = 1e-8,
                         min = 1e-12, step = 1e-8),
            conditionalPanel("input.model_type == 'efrm'",
              selectInput("ef_se", "Standard errors",
                          c("Hybrid (fast)" = "hybrid",
                            "Full person bootstrap (slow, exact)" = "bootstrap")),
              numericInput("ef_reps", "Bootstrap replicates", value = 200,
                           min = 50, step = 50))),
          accordion_panel("Scoring & anchors", icon = bs_icon("key"),
            conditionalPanel("['dich','pcm','rsm'].indexOf(input.model_type) > -1",
              fileInput("key_file", "Multiple-choice key (CSV: item,key — use \"A/C\" for double keys — or item,option,score for polytomous option scoring)",
                        accept = ".csv", placeholder = "optional"),
              fileInput("anchor_file", "Anchors for equating (CSV: item,k,tau)",
                        accept = ".csv", placeholder = "optional"),
              radioButtons("anchor_type", "Anchor as",
                           c("Individual thresholds" = "individual",
                             "Average item locations" = "average"),
                           inline = TRUE),
              p(class = "text-muted small mt-1",
                "Anchors match by item name; rows for items not present are ignored. Individual anchoring fixes each listed threshold; average anchoring fixes each item's mean location (thresholds stay free). Save an anchor file from the Items page of a previous analysis.")),
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
                  "Thresholds follow a polynomial trend across categories; useful with sparse categories. Anchors cannot be combined with this option."))))
        ),
        input_task_button("run", "Estimate", icon = bs_icon("play-fill"),
                          type = "primary", class = "w-100 btn-lg mt-2"),
        conditionalPanel("output.has_override",
          uiOutput("override_status"),
          actionButton("reset_override", "Reset overrides",
                       class = "btn-outline-warning w-100 mt-1")),
        p(class = "text-muted small mt-3",
          "Estimation: pairwise conditional maximum likelihood (Andrich & Luo 2003).",
          "Person measures: Warm weighted likelihood.")
      ),
      uiOutput("data_main")
    )
  )

# -------------------------------------------------------------- SUMMARY --
panel_summary <- nav_panel("Summary", icon = bs_icon("clipboard-data"),
    uiOutput("vboxes"),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      card(info_header("Test of fit",
             "The total item-trait chi-square tests the invariance of item ordering across the trait; a significant result means at least one item's difficulty is not invariant across class intervals (Andrich & Marais 2019). The cell df factor scales each item's chi-square degrees of freedom by the proportion of class-interval cells with enough responders to contribute."),
           card_body(verbatimTextOutput("fit_summary"))),
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
  )

# ---------------------------------------------------------------- ITEMS --
panel_items <- nav_panel("Items", icon = bs_icon("list-check"),
    uiOutput("items_vboxes"),
    div(class = "mb-2 d-flex align-items-center gap-3 flex-wrap",
        checkboxInput("show_obs",
                      "Show observed points on the category and threshold curves",
                      TRUE, width = "420px"),
        downloadButton("dl_anchors", "Save anchors (CSV: item,k,tau)",
                       class = "btn-outline-secondary mb-3")),
    layout_columns(col_widths = c(7, 5),
      tableCard("items_tbl", "Item statistics",
        controls = cols_switch("items_full"),
                "Click a row to explore that item on the right. Fit residual ~ N(0,1) under fit; misfit flag uses BH-adjusted chi-square probabilities.",
                info = "Fit residual: log-of-mean-square statistic, approximately N(0,1) under fit; |values| > 2.5 are conventionally flagged (Andrich & Marais 2019).",
                footer = uiOutput("items_note")),
      navset_card_underline(
        title = uiOutput("sel_item_title", inline = TRUE),
        full_screen = TRUE,
        nav_panel("Curve",
                  plotOutput("icc", height = "440px"),
                  div(class = "text-end",
                      downloadButton("icc_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("icc_pdf", "PDF", class = "btn-outline-secondary btn-xs"))),
        nav_panel("Categories",
                  plotOutput("ccc", height = "440px"),
                  div(class = "text-end",
                      downloadButton("ccc_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("ccc_pdf", "PDF", class = "btn-outline-secondary btn-xs"))),
        nav_panel("Thresholds",
                  plotOutput("tpc", height = "440px"),
                  div(class = "text-end",
                      downloadButton("tpc_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("tpc_pdf", "PDF", class = "btn-outline-secondary btn-xs"))),
        nav_panel("Frequencies",
                  plotOutput("cfreq", height = "440px"),
                  div(class = "text-end",
                      downloadButton("cfreq_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("cfreq_pdf", "PDF", class = "btn-outline-secondary btn-xs"))),
        nav_panel("Chi-square",
                  uiOutput("chisq_caption"),
                  DTOutput("chisq_int_tbl"),
                  h6("Response categories by class interval", class = "mt-3"),
                  DTOutput("chisq_cat_tbl"),
                  div(class = "text-end mt-2",
                      downloadButton("chisq_int_csv", "Intervals CSV",
                                     class = "btn-outline-secondary btn-xs"),
                      downloadButton("chisq_cat_csv", "Categories CSV",
                                     class = "btn-outline-secondary btn-xs"))))),
    uiOutput("pc_comp_ui"),
    conditionalPanel("output.has_mc == true",
    layout_columns(col_widths = 12,
      tableCard("distractor_tbl", "Distractor analysis",
                "Locations use the rest measure; a distractor whose takers are abler than the keyed option's flags a possible miskey."),
      plotCard("distractor_plot", "Option curves"),
      card(card_header("Polytomous option scoring (Andrich & Styles 2011)"),
           card_body(
             p(class = "text-muted",
               "Propose partial credit for informative distractors from the rest-measure evidence. Review substantively, download, edit if needed, and upload as the key (item,option,score) to refit."),
             layout_columns(col_widths = c(3, 3, 3, 3),
               numericInput("rescore_min_n", "Min takers", 20, min = 5, step = 5),
               numericInput("rescore_z", "Separation z", 1.96, min = 0.5, step = 0.1),
               actionButton("rescore_go", "Propose option scores",
                            class = "btn-primary mt-4"),
               div(class = "mt-4", cols_switch("rescore_full"))),
             DT::DTOutput("rescore_tbl"),
             downloadButton("dl_rescore", "Download proposed key CSV",
                            class = "btn-sm")))))
  )

# -------------------------------------------------------------- PERSONS --
panel_persons <- nav_panel("Persons", icon = bs_icon("people"),
    uiOutput("persons_vboxes"),
    tableCard("person_tbl", "Person estimates",
        controls = cols_switch("persons_full"),
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
  )

# ----------------------------------------------------------------- TEST --
panel_test <- nav_panel("Test", icon = bs_icon("graph-up"),
    layout_columns(col_widths = 12,
      plotCard("thrmap", "Threshold map"),
      plotCard("imap", "Item map: location by fit residual")),
    layout_columns(col_widths = 12,
      plotCard("tcc", "Test characteristic curve"),
      plotCard("tif", "Test information & SEM")),
    plotCard("guttman", "Guttman scalogram", height = "640px")
  )

# ------------------------------------------------------------------ DIF --
panel_dif <- nav_panel("DIF", icon = bs_icon("sliders"),
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
      tableCard("dif_tbl", "DIF analysis of variance",
        controls = cols_switch("dif_full"),
                info = "ANOVA of standardised residuals: a significant factor effect indicates uniform DIF, a significant factor-by-class-interval interaction indicates non-uniform DIF (Andrich & Marais 2019).",
                footer = uiOutput("dif_note")),
      tableCard("dif_tukey_tbl", "Tukey HSD comparisons",
                "Pairwise level comparisons for significant, non-superseded group terms (factorial mode)."),
      card(card_header("DIF size in logits (practical significance)"),
           card_body(
             p(class = "text-muted",
               "Resolves the selected item by the selected factor (in factorial mode: by every significant, non-superseded group term) and reports pairwise location differences in logits with Holm familywise adjustment. Differences of at least the criterion are flagged as practically significant."),
             layout_columns(col_widths = c(3, 3, 3, 3),
               numericInput("dif_size_flag", "Practical criterion (logits)",
                            0.5, min = 0.1, step = 0.1),
               numericInput("dif_size_minn", "Min responders per level", 20,
                            min = 5, step = 5),
               actionButton("dif_size_go", "Compute DIF size",
                            class = "btn-primary mt-4"),
               div(class = "mt-4", cols_switch("difsize_full"))),
             DT::DTOutput("dif_size_tbl"),
             downloadButton("dl_dif_size", "Download CSV", class = "btn-sm"))),
      plotCard("dif_icc", "ICC by group (DIF plot)")
    )
  )

# --------------------------------------------------------------- FACETS --
panel_facets <- nav_panel("Facets", icon = bs_icon("person-badge"),
    layout_sidebar(
      sidebar = sidebar(width = 280, open = "always",
        selectInput("facet_sel", "Facet", NONE),
        p(class = "text-muted small",
          "Severities from the joint calibration (positive = more severe). Pooled fit residuals beyond +/-2.5 flag inconsistent levels. Long-format analyses only.")),
      tableCard("facet_tbl", "Facet severities and fit",
        controls = cols_switch("facets_full")),
      plotCard("facet_plot", "Severity caterpillar plot"),
      tableCard("facet_int_tbl", "Item-by-facet interactions",
                "Shown when the analysis was run in interactive facet mode; gamma is the extra severity of a level on a particular item.")
    )
  )

# -------------------------------------------------------------- EQUATING --
panel_equating <- nav_panel("Equating", icon = bs_icon("arrow-left-right"),
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
      tableCard("eq_tbl", "Common-item comparison",
                controls = cols_switch("eq_full")),
      plotCard("eq_plot", "Equating plot")
    )
  )

# --------------------------------------------------------------- FRAMES --
panel_frames <- nav_panel("Frames", icon = bs_icon("grid-3x3"),
    layout_sidebar(
      sidebar = sidebar(width = 290, open = "always",
        selectInput("frame_item", "Item for ICC across frames", NONE),
        p(class = "text-muted small",
          "Units rho = alpha (set) x phi (group) on a common arbitrary scale. Within a frame all curves are parallel; across frames they fan with the unit. Extended frame of reference analyses only.")),
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(7, 5)),
        tableCard("frame_tbl", "Frames: units, origins, pooled fit",
                  controls = cols_switch("frames_full")),
        div(tableCard("phi_tbl", "Person group units (phi)"),
            tableCard("alpha_tbl", "Item set units (alpha) and locations"))),
      layout_columns(col_widths = 12,
        plotCard("frame_plot", "Frame units"),
        plotCard("frame_icc", "ICC across frames")),
      card(card_header("Equal-unit comparison"),
           card_body(verbatimTextOutput("efrm_cmp")))
    )
  )

# ------------------------------------------------------- DIMENSIONALITY --
panel_dim <- nav_panel("Dimensionality", icon = bs_icon("diagram-3"),
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
        card(info_header("Unidimensionality t-test (Smith)",
               "Each person is measured separately on the two item subsets and the estimates compared by t-test; unidimensionality is questioned when clearly more than 5% of tests are significant."),
             card_body(verbatimTextOutput("dim_txt"))),
        plotCard("scree", "Scree of the residual components")),
      layout_columns(col_widths = 12,
        plotCard("pca_plot", "Residual first contrast"),
        tableCard("loadings_tbl", "Component loadings (first 10)",
                  controls = cols_switch("load_full"))),
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
  )

# ------------------------------------------------------ LOCAL DEPENDENCE --
panel_ld <- nav_panel("Local dependence", icon = bs_icon("link-45deg"),
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
  )

# ------------------------------------------------------------- GUESSING --
panel_guess <- nav_panel("Guessing", icon = bs_icon("question-diamond"),
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
  )

# ------------------------------------------------------------------ BTL --
panel_btl <- nav_panel("BTL", icon = bs_icon("trophy"),
    layout_columns(col_widths = 12,
      card(card_header("Paired comparisons (Bradley-Terry-Luce)"),
           card_body(uiOutput("btl_boxes"),
                     verbatimTextOutput("btl_summary"))),
      tableCard("btl_obj_tbl", "Object locations and fit",
        controls = cols_switch("btl_full"),
                "Conditional (person-free) estimation with sum-zero identification and sandwich standard errors; the fit residual is the log-of-mean-square statistic over each object's comparisons (Andrich & Marais 2019)."),
      plotCard("btl_plot", "Object caterpillar"),
      tableCard("btl_pairs_tbl", "Pairwise goodness of fit",
                "Observed against expected win proportions for every pair; the total chi-square tests the BTL structure."),
      tableCard("btl_judges_tbl", "Judge fit",
                "Available when a judge column is nominated; an erratic judge carries a large positive fit residual, exactly as an erratic person does."))
  )

# -------------------------------------------------------------- COMPARE --
panel_compare <- nav_panel("Compare", icon = bs_icon("columns-gap"),
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
                "Reference for the log-likelihood comparison is the fit chosen in the sidebar.",
                controls = cols_switch("cmp_full"))
    )
  )

# --------------------------------------------------------------- EXPORT --
panel_export <- nav_panel("Export", icon = bs_icon("download"),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      card(card_header(span(bs_icon("file-earmark-text"),
                            " Analysis report (single HTML file)")),
        card_body(
          p("A self-contained HTML report of the current analysis: the summary statistics, every diagnostic table, and the test-level plots embedded as images. One portable file, ready to e-mail or archive."),
          p(class = "text-muted small",
            "Available for Rasch analyses (dichotomous, PCM, RSM, MFRM, EFRM); paired-comparison (BTL) fits are not covered."),
          downloadButton("dl_report", "Download report (HTML)",
                         class = "btn-primary btn-lg", icon = icon("file")))),
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

# ------------------------------------------------------------ ASSEMBLY --
# Workflow order: data -> summary -> items -> persons -> test, then the
# structure, invariance, and utility menus; status chips and the dark-mode
# toggle sit at the right of the navbar.
ui <- page_navbar(
  id = "nav",
  title = span("rmt"),
  theme = theme,
  # normal scrolling pages: never compress content to fit the viewport
  fillable = FALSE,
  header = tagList(
    tags$head(tags$style(css)),
    busyIndicatorOptions(spinner_type = "ring2")),
  panel_data,
  panel_summary,
  panel_items,
  panel_persons,
  panel_test,
  nav_menu("Structure",
    panel_dim,
    panel_ld,
    panel_guess),
  nav_menu("Invariance",
    panel_dif,
    panel_equating,
    panel_facets,
    panel_frames),
  nav_menu("More",
    panel_btl,
    panel_compare,
    panel_export),
  nav_spacer(),
  nav_item(uiOutput("nav_status")),
  nav_item(downloadLink("dl_report_nav", label = bs_icon("file-earmark-text"),
                        class = "nav-link px-2",
                        title = "Analysis report (HTML)")),
  nav_item(input_dark_mode(id = "app_mode"))
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

  # Data page main area: an empty-state hero before any data is loaded, the
  # summary strip + preview + R-code disclosure once data is in
  .demo_labels <- c(dich = "Multiple choice, dichotomous",
                    pcm = "Polytomous (PCM)",
                    rsm = "Rating scale (RSM)",
                    mfrm = "Ratings, long format (MFRM)",
                    efrm = "Item sets x groups (EFRM)",
                    btl = "Paired comparisons (BTL)")
  .demo_chip_labels <- c(dich = "Multiple choice", pcm = "Polytomous (PCM)",
                         rsm = "Rating scale", mfrm = "Ratings (MFRM)",
                         efrm = "Frames (EFRM)", btl = "Paired comparisons")
  output$data_main <- renderUI({
    if (identical(input$demo_choice %||% "none", "none") && is.null(input$file)) {
      div(class = "mx-auto", style = "max-width: 760px; margin-top: 8vh;",
        card(
          card_body(class = "empty-state",
            bs_icon("clipboard-data", size = "3rem", class = "text-primary mb-2"),
            h2("Welcome to rmt"),
            p(class = "lead mb-4",
              "Rasch measurement: pairwise conditional estimation, the complete test-of-fit suite, and every diagnostic table and plot — with one-click export."),
            div(class = "d-flex justify-content-center gap-2 flex-wrap mb-4",
              actionButton("hero_demo", "Try an example dataset",
                           icon = icon("table"), class = "btn-primary btn-lg"),
              tags$button(class = "btn btn-outline-primary btn-lg", type = "button",
                          onclick = "document.getElementById('file').click();",
                          bs_icon("upload"), " Upload data")),
            p(class = "text-muted small mb-2", "Or start from a specific example:"),
            div(class = "d-flex flex-wrap gap-1 justify-content-center",
              lapply(names(.demo_chip_labels), function(k)
                actionButton(paste0("demo_chip_", k), .demo_chip_labels[[k]],
                             class = "btn-outline-secondary btn-sm"))))))
    } else {
      tagList(
        uiOutput("data_strip"),
        card(card_header("Data preview"),
             card_body(uiOutput("data_info"), DTOutput("preview"), padding = 12)),
        accordion(id = "rcode_acc", open = FALSE, class = "mt-3",
          accordion_panel("R code for this analysis", icon = bs_icon("code-slash"),
            p(class = "text-muted small mb-2",
              "The exact rmt call reproducing the current run; updates on every estimation."),
            verbatimTextOutput("rcode_fit"))))
    }
  })
  observeEvent(input$hero_demo,
    updateSelectInput(session, "demo_choice", selected = "pcm"))
  lapply(c("dich", "pcm", "rsm", "mfrm", "efrm", "btl"), function(k)
    observeEvent(input[[paste0("demo_chip_", k)]],
      updateSelectInput(session, "demo_choice", selected = k)))
  output$data_strip <- renderUI({
    df <- raw_data()
    vals <- as.matrix(df)
    miss <- 100 * mean(is.na(vals) | trimws(vals) == "", na.rm = FALSE)
    layout_column_wrap(width = "160px", fill = FALSE, class = "mb-3",
      value_box("Rows", format(nrow(df), big.mark = ","),
                showcase = bs_icon("list-ol"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Columns", ncol(df), showcase = bs_icon("layout-three-columns"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Missing", sprintf("%.1f%%", miss),
                showcase = bs_icon("droplet-half"),
                showcase_layout = "left center",
                theme = if (miss > 20) "warning" else "secondary"))
  })
  output$data_info <- renderUI({
    df <- raw_data()
    p(class = "text-muted",
      sprintf("%d rows x %d columns. Nominate the column roles in the sidebar, then press Estimate. Missing responses may be left blank or coded as -1; any negative score is read as missing.",
              nrow(df), ncol(df)))
  })
  output$preview <- renderDT({
    datatable(head(raw_data(), 200), rownames = FALSE, style = "bootstrap5",
              class = "table-sm compact hover order-column",
              options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
  })
  output$rcode_fit <- renderText({
    validate(need(!is.null(rcode_str()),
                  "Run an analysis to see the reproducible R code."))
    rcode_str()
  })

  # ----------------------------------------------------------------- fit --
  override_fit <- reactiveVal(NULL)
  override_desc <- reactiveVal(NULL)
  # the exact rmt call reproducing the current run (built alongside the fit)
  rcode_str <- reactiveVal(NULL)
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
    # reproducible-code pieces (spliced into the branch-specific call below)
    src_line <- if (!identical(input$demo_choice %||% "none", "none"))
      paste0("# dat: the \"", .demo_labels[[input$demo_choice]],
             "\" example dataset generated by the app")
    else
      sprintf('dat <- read.csv("%s"%s, check.names = FALSE)',
              input$file$name,
              if (tolower(tools::file_ext(input$file$name)) %in%
                  c("tsv", "txt")) ', sep = "\\t"' else "")
    code_args_common <- c(
      if (!is.null(ng)) paste0("n_groups = ", ng),
      if (!is.na(adjN)) paste0("adjust_N = ", adjN))
    code_est <- c(paste0("maxit = ", max(5, input$maxit %||% 60)),
                  paste0("tol = ", format(max(1e-12, input$tol %||% 1e-8))))
    code_call <- NULL
    code_notes <- character(0)
    withProgress(message = "Estimating (pairwise conditional ML)…", value = 0.3, {
      fit <- tryCatch({
        if (identical(input$model_type, "btl")) {
          if (any(c(input$bt_a, input$bt_b, input$bt_win) == NONE))
            stop("nominate the object A, object B, and winner columns")
          code_call <- paste0("fit <- btl(dat,\n  ", paste(c(
            paste0("object_a = ", qstr(input$bt_a)),
            paste0("object_b = ", qstr(input$bt_b)),
            paste0("winner = ", qstr(input$bt_win)),
            if (!is.null(input$bt_judge) && input$bt_judge != NONE)
              paste0("judge = ", qstr(input$bt_judge)),
            if (!is.null(input$bt_count) && input$bt_count != NONE)
              paste0("count = ", qstr(input$bt_count)),
            paste0("ties = ", qstr(input$bt_ties %||% "drop")),
            code_est), collapse = ",\n  "), ")")
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
          sm <- ef_setmap()
          code_call <- paste0("fit <- rasch_efrm(dat,\n  ", paste(c(
            paste0("item_sets = ",
                   paste(deparse(sm), collapse = "\n    ")),
            if (is.null(input$ef_group) || input$ef_group == NONE)
              'groups = rep("(all)", nrow(dat))'
            else paste0("groups = ", qstr(input$ef_group)),
            if (!is.null(input$ef_id) && input$ef_id != NONE)
              paste0("id = ", qstr(input$ef_id)),
            paste0("items = ", qvec(names(sm))),
            code_args_common,
            code_est,
            paste0("se_method = ", qstr(input$ef_se %||% "hybrid")),
            if (!is.null(input$ef_reps) && !is.na(input$ef_reps))
              paste0("boot_reps = ", max(50, input$ef_reps))),
            collapse = ",\n  "), ")")
          rasch_efrm(df,
                     item_sets = sm,
                     groups = if (is.null(input$ef_group) ||
                                  input$ef_group == NONE)
                       rep("(all)", nrow(df)) else input$ef_group,
                     id = if (!is.null(input$ef_id) && input$ef_id != NONE)
                       input$ef_id else NULL,
                     items = names(sm),
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
          code_call <- paste0("fit <- rasch_mfrm(dat,\n  ", paste(c(
            paste0("person = ", qstr(input$lp_person)),
            paste0("item = ", qstr(input$lp_item)),
            paste0("score = ", qstr(input$lp_score)),
            paste0("facets = ", qvec(input$lp_facets)),
            code_args_common,
            if (!is.null(input$lp_interaction) && input$lp_interaction != NONE)
              paste0("interaction = ", qstr(input$lp_interaction)),
            code_est), collapse = ",\n  "), ")")
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
          code_notes <- c(
            if (!is.null(anc))
              "# anchor rows for items not in the data are dropped before fitting",
            if (!is.null(anc) && identical(input$anchor_type, "average"))
              "# average anchoring: the anchor file is collapsed to one mean location per item")
          code_call <- paste0("fit <- rasch(dat,\n  ", paste(c(
            paste0("model = ",
                   qstr(if (identical(input$model_type, "rsm")) "RSM" else "PCM")),
            if (!is.null(idc)) paste0("id = ", qstr(idc)),
            if (!is.null(fac)) paste0("factors = ", qvec(fac)),
            if (!is.null(its)) paste0("items = ", qvec(its)),
            code_args_common,
            if (!is.null(anc) && !is.null(input$anchor_file))
              paste0("anchors = read.csv(", qstr(input$anchor_file$name), ")"),
            if (!is.null(mc_key) && !is.null(input$key_file))
              paste0("key = read.csv(", qstr(input$key_file$name), ")")
            else if (!is.null(mc_key))
              'key = setNames(rep("A", 15), sprintf("I%02d", 1:15))',
            if (!is.null(pcc)) paste0("pc_components = ", pcc),
            code_est), collapse = ",\n  "), ")")
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
    if (!is.null(code_call))
      rcode_str(paste(c("library(rmt)", "", src_line,
                        if (length(code_notes)) c("", code_notes), "",
                        code_call), collapse = "\n"))
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
      try(nav_select("nav", "BTL", session = session), silent = TRUE)
      return(NULL)
    }
    btl_fit(NULL)
    try(nav_select("nav", "Summary", session = session), silent = TRUE)
    fit
  })
  btl_fit <- reactiveVal(NULL)
  fit <- reactive({
    f <- override_fit()
    if (is.null(f)) f <- analysis()
    req(f); f
  })

  output$has_mc <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    !is.null(f) && !is.null(f$mc)
  })
  outputOptions(output, "has_mc", suspendWhenHidden = FALSE)

  output$sel_item_title <- renderUI(span(class = "fw-semibold",
    tryCatch(sel_item(), error = function(e) "Selected item")))

  # only offer the pages that apply to the current analysis: Facets needs a
  # many-facet fit, Frames an extended-frames fit, BTL a paired-comparison
  # analysis, and Guessing a dichotomous one. Everything else stays.
  observe({
    f <- tryCatch(fit(), error = function(e) NULL)
    bf <- btl_fit()
    show <- function(target, on) {
      fun <- if (isTRUE(on)) nav_show else nav_hide
      try(fun("nav", target, session = session), silent = TRUE)
    }
    show("Facets", inherits(f, "rasch_mfrm"))
    show("Frames", inherits(f, "rasch_efrm"))
    show("BTL", !is.null(bf))
    show("Guessing", !is.null(f) && !inherits(f, "rasch_mfrm") &&
           !inherits(f, "rasch_efrm") && max(f$m) == 1L)
    rasch_on <- !is.null(f)
    for (tgt in c("Summary", "Items", "Persons", "Test", "Dimensionality",
                  "Local dependence", "DIF", "Equating"))
      show(tgt, rasch_on)
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
  # APA-leaning DT wrapper: Bootstrap 5 skin, right-aligned numerics, paging
  # controls only when the table needs them. `fit_col` colours fit residuals
  # beyond |2.5| with the theme danger colour; `p_bold` bolds p-values < .05.
  # curated display columns: fit objects carry every statistic, but the
  # tables show a readable core; the per-table "detailed columns" switch
  # reveals the rest (CSV downloads always contain everything)
  CORE <- list(
    items = c("item", "location", "se", "fit_resid", "infit_ms", "outfit_ms",
              "chisq", "df", "p_adj", "misfit"),
    person = c("id", "raw", "max_raw", "theta", "se", "extreme", "fit_resid"),
    dif = c("factor", "item", "F_uniform", "p_uniform_adj", "eta2_uniform",
            "F_nonuniform", "p_nonuniform_adj", "eta2_nonuniform",
            "uniform_DIF", "nonuniform_DIF"),
    facet = c("level", "severity", "se", "n", "fit_resid"),
    btl_obj = c("object", "location", "se", "comparisons", "wins", "fit_resid"),
    btl_judge = c("judge", "n", "fit_resid", "misfit"),
    equate = c("item", "location_1", "location_2", "adj_difference", "t",
               "p_adj", "drift"),
    dif_size = c("item", "term", "level_a", "level_b", "difference", "lower",
                 "upper", "p_adj", "significant", "practical"),
    frames = c("set", "group", "rho", "se_log_rho", "origin", "fit_resid",
               "n_responses"),
    compare = c("label", "model", "persons", "items", "two_delta_ll",
                "chisq_per_df", "item_fit_sd", "person_fit_sd", "PSI", "alpha"),
    rescore = c("item", "option", "keyed", "n", "prop", "mean_location",
                "z_sep", "proposed"))
  curate <- function(d, which, full = FALSE, extra = NULL) {
    if (isTRUE(full)) return(d)
    keep <- c(CORE[[which]], extra)
    d[, intersect(keep, names(d)), drop = FALSE]
  }
  # display headers for the tables; downloads always keep the raw names
  DISPLAY_NAMES <- c(
    fit_resid = "Fit resid", fit_resid_pooled = "Pooled fit resid",
    natural_resid = "Natural resid", infit_ms = "Infit MS",
    outfit_ms = "Outfit MS", infit_z = "Infit z", outfit_z = "Outfit z",
    se = "SE", theta = "Location", max_raw = "Max score", raw = "Raw score",
    n_items = "Items", chisq = "Chi-sq", df_fit = "Fit df", p = "p",
    p_adj = "Adj. p", p_bonf = "Bonf. p", p_anova = "ANOVA p",
    F_anova = "ANOVA F", F_uniform = "F (uniform)",
    F_nonuniform = "F (non-unif.)", p_uniform_adj = "Adj. p (uniform)",
    p_nonuniform_adj = "Adj. p (non-unif.)",
    eta2_uniform = "η² (uniform)",
    eta2_nonuniform = "η² (non-unif.)",
    eta2_partial = "Partial η²",
    uniform_DIF = "Uniform", nonuniform_DIF = "Non-uniform",
    mean_location = "Mean location", point_biserial = "Point-biserial",
    se_location = "SE", z_sep = "Separation z",
    alpha_drop = "α if deleted", item_total = "Item-total r",
    item_rest = "Item-rest r", di = "Discrimination", cum_pct = "Cum. %",
    exp_prop = "Expected", obs_prop = "Observed", obs_mean = "Observed mean",
    exp_value = "Expected value", theta_mean = "Mean location",
    theta_max = "Max location", chisq_per_df = "Chi-sq/df",
    two_delta_ll = "2Δ log-lik", se_log_phi = "SE (log φ)",
    se_log_alpha = "SE (log α)", se_log_rho = "SE (log ρ)",
    mu = "Origin", comparisons = "Comparisons",
    obs_p = "Observed", est_p = "Expected", obs_t = "Threshold prop.")
  # p-value columns render as "<0.001" / 3 dp on the client, so sorting
  # still uses the raw value; detection runs on the ORIGINAL column names
  P_COL_RE <- "^p$|^p_|_p$|^prob$|p_tukey|p_anova|p_adj|p_bonf|p_uniform|p_nonuniform"
  P_RENDER <- DT::JS("function(data,type,row){ if(type==='display'){ if(data===null||data==='') return ''; var x=Number(data); return x<0.001 ? '&lt;0.001' : x.toFixed(3);} return data; }")
  num_dt <- function(d, digits = 3, fit_col = NULL, p_bold = NULL, ...) {
    orig <- names(d)
    num <- vapply(d, is.numeric, TRUE)
    # integer-valued columns (counts, whole-number df) show no decimals;
    # fractional df columns fail the test and keep the 3-dp rounding
    intcol <- vapply(d, function(v)
      is.numeric(v) && all(is.na(v) | v == round(v)), TRUE)
    pcol <- num & grepl(P_COL_RE, orig)
    opts <- list(pageLength = 15, scrollX = TRUE,
                 dom = if (nrow(d) > 12) "tip" else "t")
    cdefs <- list()
    if (any(num))
      cdefs[[length(cdefs) + 1L]] <- list(className = "dt-right",
                                          targets = which(num) - 1L)
    for (j in which(pcol))
      cdefs[[length(cdefs) + 1L]] <- list(targets = j - 1L, render = P_RENDER)
    if (length(cdefs)) opts$columnDefs <- cdefs
    # display-only renaming; formatting targets are column positions, so
    # they stay tied to the original names computed above
    hit <- orig %in% names(DISPLAY_NAMES)
    names(d)[hit] <- DISPLAY_NAMES[orig[hit]]
    dt <- datatable(d, rownames = FALSE, style = "bootstrap5",
                    class = "table-sm compact hover order-column", ...,
                    options = opts)
    rnd <- which(num & !intcol & !pcol)
    if (length(rnd)) dt <- formatRound(dt, rnd, digits)
    for (fc in which(orig %in% fit_col))
      dt <- formatStyle(dt, fc, color = styleInterval(
        c(-2.5, 2.5), c("var(--bs-danger)", "inherit", "var(--bs-danger)")))
    for (pc in which(orig %in% p_bold))
      dt <- formatStyle(dt, pc,
                        fontWeight = styleInterval(0.05, c("bold", "normal")))
    dt
  }

  # ------------------------------------------------- navbar status chips --
  # compact badges once a fit exists: model, N persons, N items, PSI
  # (objects/comparisons/OSI for a paired-comparison fit); nothing before
  # the first run
  output$nav_status <- renderUI({
    chip <- function(txt, kind = "secondary")
      span(class = paste0("badge text-bg-", kind), txt)
    b <- btl_fit()
    if (!is.null(b)) {
      osi <- b$osi$PSI
      return(div(class = "nav-status d-flex align-items-center gap-1 px-2",
        chip("BTL", "primary"),
        chip(paste(nrow(b$objects), "objects")),
        chip(sprintf("%.0f comparisons", b$n_comparisons)),
        chip(if (finite1(osi)) sprintf("OSI %.2f", osi) else "OSI —",
             if (!finite1(osi)) "secondary"
             else if (osi >= 0.7) "success" else "danger")))
    }
    f <- override_fit()
    if (is.null(f)) f <- tryCatch(analysis(), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    psi <- f$psi$PSI
    div(class = "nav-status d-flex align-items-center gap-1 px-2",
      chip(f$model, "primary"),
      chip(paste(nrow(f$X), "persons")),
      chip(paste(ncol(f$X), "items")),
      chip(if (finite1(psi)) sprintf("PSI %.2f", psi) else "PSI —",
           if (!finite1(psi)) "secondary"
           else if (psi >= 0.7) "success" else "danger"))
  })

  # -------------------------------------------------------------- summary --
  output$vboxes <- renderUI({
    f <- fit()
    layout_column_wrap(width = "185px", fill = FALSE, class = "mb-3",
      value_box("Persons", nrow(f$X), showcase = bs_icon("people"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Items", ncol(f$X), showcase = bs_icon("list-check"),
                showcase_layout = "left center", theme = "primary"),
      value_box(span("PSI",
                     tooltip(bs_icon("info-circle", class = "ms-1"),
                             "Person Separation Index: the proportion of variance in person estimates not attributable to measurement error; 0.7 is a conventional minimum for distinguishing groups of persons (Andrich & Marais 2019).")),
                if (finite1(f$psi$PSI)) sprintf("%.3f", f$psi$PSI) else "—",
                showcase = bs_icon("speedometer2"),
                showcase_layout = "left center",
                theme = if (!finite1(f$psi$PSI)) "secondary"
                        else if (f$psi$PSI >= 0.7) "success" else "warning",
                p(class = "small mb-0",
                  if (finite1(f$psi_noext$PSI))
                    sprintf("%.3f no extremes", f$psi_noext$PSI)
                  else "— no extremes")),
      value_box("Alpha",
                if (finite1(f$alpha$alpha)) sprintf("%.3f", f$alpha$alpha)
                else "—",
                showcase = bs_icon("calculator"),
                showcase_layout = "left center",
                theme = if (!finite1(f$alpha$alpha)) "secondary"
                        else if (f$alpha$alpha >= 0.7) "success" else "warning",
                p(class = "small mb-0",
                  if (isFALSE(f$alpha$applicable))
                    sprintf("complete cases only (n = %d)", f$alpha$n)
                  else sprintf("n = %d complete", f$alpha$n))),
      value_box("Item-trait p",
                if (finite1(f$total_chisq_p)) fmt_p(f$total_chisq_p) else "—",
                showcase = bs_icon("clipboard-check"),
                showcase_layout = "left center",
                theme = if (!finite1(f$total_chisq_p)) "secondary"
                        else if (f$total_chisq_p < 0.05) "danger" else "success"),
      value_box("Power of fit", f$power_of_fit,
                showcase = bs_icon("lightning-charge"),
                showcase_layout = "left center", theme = "secondary")
    )
  })

  output$fit_summary <- renderPrint({
    f <- fit(); ss <- f$summary_stats
    cat(sprintf("Model: %s  |  Estimation: pairwise conditional ML (%s, %d iterations)\n",
                f$model, if (f$est$converged) "converged" else "NOT CONVERGED",
                f$est$iterations))
    cat(sprintf("Total item-trait chi-square: %.3f on %d df, p = %s  (%d class intervals)\n",
                f$total_chisq, f$total_df, fmt_p(f$total_chisq_p), f$n_groups))
    cat(sprintf("Item fit residual:   mean %6.2f  SD %5.2f  skew %5.2f  kurt %5.2f  (ideal 0, 1)\n",
                f$item_fit_summary$mean, f$item_fit_summary$sd,
                f$item_fit_summary$skewness, f$item_fit_summary$kurtosis))
    cat(sprintf("Person fit residual: mean %6.2f  SD %5.2f  skew %5.2f  kurt %5.2f  (ideal 0, 1)\n",
                f$person_fit_summary$mean, f$person_fit_summary$sd,
                f$person_fit_summary$skewness, f$person_fit_summary$kurtosis))
    cat(sprintf("Fit-location correlation: items %.3f, persons %.3f\n",
                ss$cor_item_fit_location, ss$cor_person_fit_location))
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
      if ("extrapolated" %in% names(d)) {
        if (!isTRUE(any(d$extrapolated))) d$extrapolated <- NULL
        else d$extrapolated <- ifelse(d$extrapolated, "*", "")
      }
      datatable(d, rownames = FALSE, style = "bootstrap5",
                class = "table-sm compact hover order-column",
                options = list(pageLength = 15, scrollX = TRUE,
                               dom = if (nrow(d) > 12) "tip" else "t")) |>
        formatRound(c("theta", "se"), 3) |>
        formatRound("cum_pct", 1)
    } else {
      validate(need(!is.null(f$score_curves),
                    "No score conversion available for this fit."))
      datatable(f$score_curves, rownames = FALSE, style = "bootstrap5",
                class = "table-sm compact hover order-column",
                caption = "Raw scores are not sufficient under unequal frame units; per-group expected-score curves replace the score table.",
                options = list(pageLength = 15, scrollX = TRUE, dom = "tip")) |>
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
    card(info_header("Likelihood-ratio test (PCM vs rating)",
           "Compares the partial credit model against the more parsimonious rating parameterisation with common thresholds; a non-significant result supports the rating model."),
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
    cat(sprintf("Raw composite chi-square %.3f on %d df, p = %s (conventional display; anticonservative)\n",
                r$chisq, r$df, fmt_p(r$p)))
    if (is.finite(r$chisq_adj))
      cat(sprintf("Adjusted chi-square %.3f, p = %s (Kent 1982 calibration)\n",
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
  output$items_vboxes <- renderUI({
    f <- fit()
    mis <- sum(f$items$misfit, na.rm = TRUE)
    dis <- sum(vapply(f$thresholds_diag, function(d)
      !d$ordered && length(d$thresholds) > 1, TRUE))
    layout_column_wrap(width = "200px", fill = FALSE, class = "mb-3",
      value_box("Items", nrow(f$items), showcase = bs_icon("list-check"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Misfitting items", mis,
                showcase = bs_icon("exclamation-triangle"),
                showcase_layout = "left center",
                theme = if (mis > 0) "danger" else "success"),
      value_box("Disordered thresholds", dis,
                showcase = bs_icon("arrow-down-up"),
                showcase_layout = "left center",
                theme = if (dis > 0) "warning" else "success"))
  })
  output$items_note <- renderUI({
    f <- fit(); d <- f$items
    dis <- names(which(vapply(f$thresholds_diag, function(x)
      !x$ordered && length(x$thresholds) > 1, TRUE)))
    sprintf("Note. %d of %d items beyond |fit residual| 2.5; %d flagged by adjusted chi-square; disordered thresholds: %s.",
            sum(abs(d$fit_resid) > 2.5, na.rm = TRUE), nrow(d),
            sum(d$misfit, na.rm = TRUE),
            if (length(dis)) paste(dis, collapse = ", ") else "none")
  })
  # any chi-square sample-size adjustment is applied inside the fit (the
  # adjust_N run control), so the table shows the fit's own statistics
  register_table("items_tbl", function() fit()$items, function() {
    d <- curate(fit()$items, "items", full = isTRUE(input$items_full),
                extra = if (length(unique(fit()$m)) > 1) "max")
    d$misfit <- ifelse(d$misfit, "*", "")
    num_dt(d, selection = "single", fit_col = "fit_resid",
           p_bold = c("p_adj", "p_anova"))
  })

  # per-class-interval breakdown of the selected item's chi-square
  chisq_res <- reactive(chisq_detail(fit(), sel_item()))
  output$chisq_caption <- renderUI({
    cd <- chisq_res()
    p(class = "small mb-2", HTML(sprintf(
      "<b>%s</b> (location %.3f): total chi-square <b>%.3f</b> on %d df, p = %s; whole-sample mean = %.3f. Intervals with fewer than 2 responders carry no chi-square contribution.",
      cd$item, cd$location, cd$chisq, cd$df, fmt_p(cd$p), cd$ave)))
  })
  output$chisq_int_tbl <- renderDT({
    d <- chisq_res()$intervals
    d$Excluded <- ifelse(d$used, "", "*")
    d$used <- NULL
    d$theta_max <- NULL   # in the CSV download, not the default display
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
    d <- curate(res$evidence, "rescore", full = isTRUE(input$rescore_full))
    if ("keyed" %in% names(d)) d$keyed <- ifelse(d$keyed, "*", "")
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
  output$persons_vboxes <- renderUI({
    f <- fit(); d <- f$person
    mis <- sum(abs(d$fit_resid) > 2.5, na.rm = TRUE)
    layout_column_wrap(width = "200px", fill = FALSE, class = "mb-3",
      value_box("Persons", nrow(d), showcase = bs_icon("people"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Extreme scores", sum(d$extreme, na.rm = TRUE),
                showcase = bs_icon("arrows-expand"),
                showcase_layout = "left center", theme = "secondary"),
      value_box("Misfitting persons", mis,
                showcase = bs_icon("exclamation-triangle"),
                showcase_layout = "left center",
                theme = if (mis > 0) "danger" else "success"))
  })
  register_table("person_tbl", function() fit()$person, function() {
    fac <- names(fit()$factors)
    d <- curate(fit()$person, "person", full = isTRUE(input$persons_full),
                extra = fac)
    dt <- datatable(d, rownames = FALSE, filter = "top", selection = "single",
                    style = "bootstrap5",
                    class = "table-sm compact hover order-column",
                    options = list(pageLength = 15, scrollX = TRUE, dom = "tip")) |>
      formatRound(names(d)[vapply(d, is.numeric, TRUE) &
                           !names(d) %in% c("raw", "max_raw", "n_items", "class_interval")], 3)
    if ("fit_resid" %in% names(d))
      dt <- formatStyle(dt, "fit_resid", color = styleInterval(
        c(-2.5, 2.5), c("var(--bs-danger)", "inherit", "var(--bs-danger)")))
    dt
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
      num_dt(d)
    } else {
      d <- dif_res()
      if (!is.null(input$dif_factor) && input$dif_factor %in% d$factor)
        d <- d[d$factor == input$dif_factor, ]
      d <- curate(d, "dif", full = isTRUE(input$dif_full))
      d$uniform_DIF <- ifelse(d$uniform_DIF, "*", "")
      d$nonuniform_DIF <- ifelse(d$nonuniform_DIF, "*", "")
      num_dt(d)
    }
  })
  output$dif_note <- renderUI({
    if (identical(input$dif_factor, FACTORIAL)) {
      d <- dif_fact()$terms
      sup <- sum(d$superseded, na.rm = TRUE)
      sprintf("Note. %d of %d terms significant after adjustment%s.",
              sum(d$significant, na.rm = TRUE), nrow(d),
              if (sup) sprintf(" (%d superseded by an interaction)", sup) else "")
    } else {
      d <- dif_res()
      parts <- vapply(split(d, d$factor), function(g)
        sprintf("%s: %d uniform, %d non-uniform", g$factor[1],
                sum(g$uniform_DIF, na.rm = TRUE),
                sum(g$nonuniform_DIF, na.rm = TRUE)), "")
      paste0("Note. Items flagged per factor - ",
             paste(parts, collapse = "; "), ".")
    }
  })
  register_table("dif_tukey_tbl", function() dif_fact()$tukey, function() {
    validate(need(identical(input$dif_factor, FACTORIAL),
                  "Choose the factorial option in the sidebar to see Tukey HSD comparisons."))
    tk <- dif_fact()$tukey
    if (!nrow(tk))
      return(datatable(data.frame(note = "no significant group terms to compare"),
                       rownames = FALSE, style = "bootstrap5",
                       class = "table-sm compact",
                       options = list(dom = "t")))
    num_dt(tk)
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
    d <- curate(d, "dif_size", full = isTRUE(input$difsize_full))
    if ("significant" %in% names(d))
      d$significant <- ifelse(d$significant, "*", "")
    if ("practical" %in% names(d))
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
                  "Run a paired-comparisons (BTL) analysis from the Data page to see results here."))
    btl_fit()
  })
  output$btl_boxes <- renderUI({
    f <- bfit()
    layout_column_wrap(width = "165px", fill = FALSE, class = "mb-3",
      value_box("Objects", nrow(f$objects), theme = "primary"),
      value_box("Comparisons", sprintf("%.0f", f$n_comparisons), theme = "primary"),
      if (!is.null(f$judges))
        value_box("Judges", nrow(f$judges), theme = "primary"),
      value_box("Object separation",
                if (finite1(f$osi$PSI)) sprintf("%.3f", f$osi$PSI) else "—",
                theme = if (!finite1(f$osi$PSI)) "secondary"
                        else if (f$osi$PSI >= 0.7) "success" else "danger"),
      value_box("Pairwise fit p",
                if (finite1(f$total_p)) fmt_p(f$total_p) else "—",
                theme = if (!finite1(f$total_p)) "secondary"
                        else if (f$total_p < 0.05) "danger" else "success"))
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
                 function() num_dt(curate(bfit()$objects, "btl_obj",
                                          full = isTRUE(input$btl_full)),
                                   fit_col = "fit_resid"))
  register_table("btl_pairs_tbl", function() bfit()$pairs,
                 function() {
    d <- bfit()$pairs
    d$chisq <- NULL   # residual^2; redundant on screen, kept in the CSV
    num_dt(d)
  })
  register_table("btl_judges_tbl", function() {
    validate(need(!is.null(bfit()$judges), "No judge column was nominated."))
    bfit()$judges
  }, function() {
    validate(need(!is.null(bfit()$judges), "No judge column was nominated."))
    d <- bfit()$judges
    d$misfit <- ifelse(!is.na(d$fit_resid) & abs(d$fit_resid) > 2.5, "*", "")
    num_dt(curate(d, "btl_judge", full = isTRUE(input$btl_full)),
           fit_col = "fit_resid")
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
    d <- curate(facet_dat(), "facet", full = isTRUE(input$facets_full))
    datatable(d, rownames = FALSE, style = "bootstrap5",
              class = "table-sm compact hover order-column",
              options = list(pageLength = 15, scrollX = TRUE,
                             dom = if (nrow(d) > 12) "tip" else "t")) |>
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
    d <- curate(eq_res()$table, "equate", full = isTRUE(input$eq_full))
    if ("drift" %in% names(d)) d$drift <- ifelse(d$drift, "*", "")
    num_dt(d)
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
                 function() num_dt(curate(efrm_fit()$frames, "frames",
                                          full = isTRUE(input$frames_full))))
  register_table("phi_tbl", function() efrm_fit()$phi_table,
                 function() num_dt(efrm_fit()$phi_table))
  # keep the fit's original set order: merge() sorts by the key, so it is
  # re-matched to fit$set_table$set
  efrm_alpha_tbl <- reactive({
    f <- efrm_fit()
    d <- merge(f$alpha_table, f$set_table[, c("set", "mu", "n_items")],
               by = "set", sort = FALSE)
    d <- d[stats::na.omit(match(f$set_table$set, d$set)), , drop = FALSE]
    rownames(d) <- NULL
    d
  })
  register_table("alpha_tbl", function() efrm_alpha_tbl(),
                 function() num_dt(efrm_alpha_tbl()))
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
                 function() {
    d <- residual_pca(fit())$loadings_matrix
    if (!isTRUE(input$load_full))
      d <- d[, intersect(c("item", "PC1", "PC2", "PC3"), names(d)),
             drop = FALSE]
    num_dt(d)
  })
  register_table("eigen_tbl", function() residual_pca(fit())$eigen_table,
                 function() {
    d <- residual_pca(fit())$eigen_table
    d$proportion <- 100 * d$proportion
    d$cumulative <- 100 * d$cumulative
    names(d)[match(c("proportion", "cumulative"), names(d))] <-
      c("Proportion %", "Cumulative %")
    num_dt(d) |> formatRound(c("Proportion %", "Cumulative %"), 1)
  })
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
                  "Run a dichotomous (multiple-choice) analysis to use the tailored guessing procedure."))
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
                  "Run a dichotomous (multiple-choice) analysis to use the tailored guessing procedure."))
    r <- guess_res()
    validate(need(!is.null(r), "Run the tailored analysis to see the comparison."))
    num_dt(r$table)
  })
  register_plot("guess_plot", function() {
    validate(need(max(fit()$m) == 1L,
                  "Run a dichotomous (multiple-choice) analysis to use the tailored guessing procedure."))
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
    if ("same_data" %in% names(d))
      d$same_data <- ifelse(d$same_data, "yes", "no")
    num_dt(curate(d, "compare", full = isTRUE(input$cmp_full)))
  })

  # ----------------------------------------------------------------- export --
  # single-file HTML report; one content function feeds both the Export-tab
  # button and the navbar icon link (Rasch fits only; BTL is notified)
  report_content <- function(file) {
    if (!is.null(btl_fit())) {
      showNotification("The HTML report covers Rasch analyses; paired-comparison (BTL) fits are not yet supported.",
                       type = "warning", duration = 8)
      stop("report unavailable for a BTL fit")
    }
    f <- override_fit()
    if (is.null(f)) f <- tryCatch(analysis(), error = function(e) NULL)
    if (is.null(f)) {
      showNotification("Run an analysis first, then download the report.",
                       type = "warning", duration = 8)
      stop("no fit to report")
    }
    withProgress(message = "Building the HTML report…", value = 0.4,
                 report_html(f, file))
  }
  output$dl_report <- downloadHandler(
    filename = function() "rmt_report.html", content = report_content)
  output$dl_report_nav <- downloadHandler(
    filename = function() "rmt_report.html", content = report_content)

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
