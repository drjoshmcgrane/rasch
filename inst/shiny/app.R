# rasch Shiny GUI
# ---------------------------------------------------------------------------
# A modern bslib interface to the full rasch analysis: data upload with ID,
# person-factor, and item column nomination; pairwise conditional ML
# estimation (Zwinderman 1995); the complete test-of-fit suite;
# diagnostic plots with per-plot PNG and PDF downloads; and a ZIP archive of
# the fitted-model outputs supported by save_outputs().
# Launch with rasch::run_app(), or shiny::runApp() from this folder.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(bsicons)
})

# Prefer the source tree only when this is genuinely the in-tree app. An
# installed app has a package/R directory too, but not the original rasch.R
# source file. This distinction prevents shiny::runApp("inst/shiny") from
# silently loading an older installed release while a newer tree is being
# developed.
.source_candidates <- c(file.path("..", "..", "R"), "R")
.source_hits <- .source_candidates[
  file.exists(file.path(.source_candidates, "rasch.R"))]
.rasch_source_dir <- if (length(.source_hits))
  normalizePath(.source_hits[1], mustWork = TRUE) else NULL

if (!is.null(.rasch_source_dir)) {
  .rasch_source_root <- normalizePath(dirname(.rasch_source_dir),
                                      mustWork = TRUE)
  .rasch_loaded_root <- if ("rasch" %in% loadedNamespaces()) {
    normalizePath(getNamespaceInfo(asNamespace("rasch"), "path"),
                  mustWork = FALSE)
  } else {
    ""
  }
  # Re-sourcing app.R during tests or development must not reload an already
  # active namespace: doing so invalidates mocked bindings and any objects
  # created from that namespace. A direct source-tree launch still loads the
  # tree, including its compiled code, through pkgload.
  if (!identical(.rasch_loaded_root, .rasch_source_root)) {
    if (!requireNamespace("pkgload", quietly = TRUE))
      stop("Running the app from a source tree requires pkgload; install rasch ",
           "first, or install pkgload for development use")
    suppressWarnings(pkgload::load_all(.rasch_source_root, quiet = TRUE))
  } else if (!"package:rasch" %in% search()) {
    suppressPackageStartupMessages(library(rasch))
  }
} else if (requireNamespace("rasch", quietly = TRUE)) {
  library(rasch)
} else {
  stop("Install rasch, or run the app from inst/shiny in the source tree")
}

# Canonical, succinct explainers shared by every result-card helper.  Locate
# the file both when Shiny has made inst/shiny the working directory and when
# the source-tree app is launched from the package root.
.help_candidates <- c("help.R", file.path("inst", "shiny", "help.R"),
                      system.file("shiny", "help.R", package = "rasch"))
.help_file <- .help_candidates[nzchar(.help_candidates) &
                                 file.exists(.help_candidates)][1]
if (is.na(.help_file)) stop("The Shiny explainer registry (help.R) is missing")
source(.help_file, local = TRUE)

# Load the bundled examples from their own source file. The same file is used
# by .app_example_data(), which keeps the code shown by the app executable.
.example_candidates <- c("examples.R", file.path("inst", "shiny", "examples.R"),
                         system.file("shiny", "examples.R", package = "rasch"))
.example_file <- .example_candidates[nzchar(.example_candidates) &
                                       file.exists(.example_candidates)][1]
if (is.na(.example_file)) stop("The bundled app examples are missing")
source(.example_file, local = TRUE)

NONE <- "(none)"
# the sentinel VALUE stays "(none)" (the server compares against it), but it
# is always displayed as "None"; selects with no meaningful pre-fit choice
# use empty choices plus a selectize placeholder instead of a sentinel row
NONE_CH <- c(None = "(none)")

# null-coalescing helper: defined here so the app does not depend on the
# base R version that introduced it (R >= 4.4)
`%||%` <- function(a, b) if (is.null(a)) b else a

# A source column is a symbol, not text to be interpolated into an R formula
# or code expression. deparse(as.name()) handles reserved words, leading
# dot-digits, spaces, colons and embedded backticks without guessing which
# names happen to be syntactic.
.app_quote_name <- function(x) {
  if (!is.character(x) || anyNA(x) || any(!nzchar(trimws(x))))
    stop("column names must be non-empty character values", call. = FALSE)
  vapply(x, function(z)
    paste(deparse(as.name(z), backtick = TRUE), collapse = ""), "")
}

# Interaction controls use opaque values. A colon-joined value is ambiguous:
# the literal predictor "a:b" paired with "c" otherwise collides with "a"
# paired with "b:c". The displayed labels remain readable, while the map
# retains the two exact source names.
.app_explanatory_interactions <- function(main) {
  if (!is.character(main) || anyNA(main) || any(!nzchar(trimws(main))) ||
      anyDuplicated(main))
    stop("explanatory main effects must have unique non-empty names",
         call. = FALSE)
  pairs <- if (length(main) >= 2L)
    utils::combn(main, 2L, simplify = FALSE) else list()
  ids <- sprintf("interaction_%04d", seq_along(pairs))
  labels <- vapply(pairs, paste, "", collapse = " × ")
  ambiguous <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[ambiguous] <- paste0(labels[ambiguous], " [", ids[ambiguous], "]")
  list(map = stats::setNames(pairs, ids),
       choices = stats::setNames(ids, labels))
}

.app_explanatory_formula <- function(main, interactions = list()) {
  if (!is.character(main) || !length(main) || anyNA(main) ||
      any(!nzchar(trimws(main))) || anyDuplicated(main))
    stop("choose at least one explanatory main effect", call. = FALSE)
  if (!is.list(interactions) || any(!vapply(interactions, function(z)
    is.character(z) && length(z) == 2L && !anyNA(z) &&
      all(nzchar(trimws(z))) && all(z %in% main) && z[1L] != z[2L], logical(1))))
    stop("the explanatory interaction selection is invalid", call. = FALSE)
  rhs <- c(lapply(main, as.name), lapply(interactions, function(z)
    call(":", as.name(z[1L]), as.name(z[2L]))))
  rhs <- Reduce(function(a, b) call("+", a, b), rhs)
  stats::as.formula(call("~", rhs), env = parent.frame())
}

# Planned MFRM DIF is item-pooled by default. Keep its selector on the same
# unit while the item table continues to display observed item-by-facet cells.
.app_dif_item_choices <- function(fit) {
  out <- if (inherits(fit, "rasch_mfrm") && !is.null(fit$virtual_map))
    unique(as.character(fit$virtual_map$item)) else
      as.character(fit$items$item)
  out[!is.na(out) & nzchar(out)]
}

# A DataTable selection can briefly retain the preceding table's row number
# while a new fit is rendering. Resolve only one valid row; otherwise use the
# first row rather than passing NA or an out-of-range index to a plot.
.app_selected_row <- function(index, n, default = 1L) {
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) ||
      n != floor(n) || n < 1L)
    stop("`n` must be one positive whole row count", call. = FALSE)
  if (length(default) != 1L || !is.numeric(default) || !is.finite(default) ||
      default != floor(default) || default < 1L || default > n)
    stop("`default` must be one available row number", call. = FALSE)
  if (length(index) != 1L || !is.numeric(index) || !is.finite(index) ||
      index != floor(index) || index < 1L || index > n)
    return(as.integer(default))
  as.integer(index)
}

# Wald results now carry a t statistic and its reference degrees of freedom.
# Projects saved by earlier releases carry the same estimate as `z` without a
# df field. Prefer the current schema whenever it is present, but keep the
# legacy statistic readable after a project is restored.
.app_wald_reference <- function(x) {
  if (!is.list(x) && !is.data.frame(x))
    stop("a Wald result must be a list or data frame", call. = FALSE)
  nms <- names(x)
  statistic <- if ("t" %in% nms) "t" else if ("z" %in% nms) "z" else NA_character_
  value <- if (!is.na(statistic)) x[[statistic]] else NA_real_
  if (!is.numeric(value) || length(value) != 1L) value <- NA_real_
  df <- if (identical(statistic, "t") && "df" %in% nms) x[["df"]]
        else if (identical(statistic, "z")) Inf else NA_real_
  if (!is.numeric(df) || length(df) != 1L) df <- NA_real_
  label <- statistic
  if (identical(statistic, "t") && !is.na(df))
    label <- sprintf("t(%s)", if (is.infinite(df)) "Inf" else
      format(df, digits = 4L, trim = TRUE))
  list(statistic = statistic, label = label, value = value, df = df)
}

.app_wald_columns <- function(x) {
  nms <- names(x)
  if ("t" %in% nms) intersect(c("t", "df"), nms)
  else intersect("z", nms)
}

# Common-item equating works with unconstrained ordinary Rasch calibrations,
# not with explanatory, response-cell, or paired-comparison calibrations.
.app_equating_fit <- function(fit) {
  inherits(fit, "rasch") &&
    !inherits(fit, c("rasch_explanatory", "rasch_btl", "rasch_mfrm",
                    "rasch_efrm"))
}

# The Local page is for sequence-dependent exposure and carry-over. A static
# first-position term also has row-level dependence_data, but no judgment
# history to display on that page.
.app_has_btl_history <- function(fit) {
  inherits(fit, "rasch_btl") && !is.null(fit$dependence_data) &&
    any(c("exposure", "carry_over") %in% names(fit$dependence_data))
}

# Launchers differ in what the app's environment can see: a development
# load resolves the package internals through the search path, while
# shiny::runApp() on an installed copy resolves exports only. The two
# project helpers are internal, so each is taken from wherever it can be
# found -- the inherited scope first, the loaded namespace otherwise -- and
# the app behaves the same everywhere, including a source tree where the
# package is not installed.
.rasch_internal <- function(name) {
  if (exists(name, inherits = TRUE)) get(name, inherits = TRUE)
  else utils::getFromNamespace(name, "rasch")
}
.read_app_project <- .rasch_internal(".read_app_project")
.save_app_project <- .rasch_internal(".save_app_project")
.seal_app_project <- .rasch_internal(".seal_app_project")
.classical_design_applicable <-
  .rasch_internal(".classical_design_applicable")
.validate_dif_bootstrap <- .rasch_internal(".validate_dif_bootstrap")
.scree_analysis <- .rasch_internal(".scree_analysis")
.validate_scree_result <- .rasch_internal(".validate_scree_result")
.estimation_label <- .rasch_internal(".estimation_label")
.validate_dimensionality_test <-
  .rasch_internal(".validate_dimensionality_test")
.validate_btl_dimensionality <-
  .rasch_internal(".validate_btl_dimensionality")
.validate_frame_invariance <-
  .rasch_internal(".validate_frame_invariance")
.sim_explanatory_departure <-
  .rasch_internal(".sim_explanatory_departure")
.has_repeated_residual_units <-
  .rasch_internal(".has_repeated_residual_units")
.fit_boot_signature <- .rasch_internal(".fit_boot_signature")
.fit_boot_signature_matches <-
  .rasch_internal(".fit_boot_signature_matches")
.fit_boot_md5 <- .rasch_internal(".fit_boot_md5")
.sim_planted_count <- .rasch_internal(".sim_planted_count")

# Controls that determine the fitted analysis. Project files retain these
# separately from display-only choices so reopening an analysis also restores
# the inputs needed to estimate it again. File-backed inputs are embedded as
# parsed data in the project resources below rather than as expired temporary
# paths.
.project_radio_inputs <- c(
  "model_type", "lp_layout", "lp_structure", "rasch_calibration",
  "thr_structure", "thr_mode", "exp_level", "bt_thr", "bt_ties",
  "anchor_type", "btlef_se", "dif_effects", "bdif_effects", "inv_se",
  "eq_source", "eq_shift")
.project_select_inputs <- c(
  "id_col", "ef_id", "ef_group", "lp_person", "lp_item", "lp_score",
  "lp_interaction", "bt_a", "bt_b", "bt_win", "bt_judge", "bt_count",
  "pc_rank", "ng", "ef_se", "ef_workers", "btlef_panel",
  "btlef_workers", "pca_component", "dim_workers", "wright_renderer",
  "wright_type", "wright_person_panels", "wright_item_panels",
  "wright_person_style")
.project_selectize_inputs <- c(
  "factor_cols", "item_cols", "ef_items", "lp_facets", "lp_items_wide",
  "bt_margin", "bt_response", "bt_order", "bt_jfactors", "exp_main",
  "exp_interactions", "bdif_factors", "dim_pos", "dim_neg", "eq_kept")
.project_checkbox_inputs <- c(
  "ef_prefix", "bt_position", "ng_auto", "eq_csv_independent",
  "eq_kept_independent", "bt_eq_independent")
.project_numeric_inputs <- c(
  "maxit", "tol", "ef_reps", "ef_seed", "btlef_boot",
  "btlef_seed", "dif_alpha", "dif_boot_B", "dif_boot_seed",
  "bdif_alpha", "bdif_boot_B", "bdif_boot_seed", "inv_boot", "inv_seed",
  "dim_boot_B", "dim_boot_seed")

# Inputs that define the base calibration. A saved analysis may also retain
# later display and diagnostic choices, but these values must continue to
# describe the data and fitted object that were actually estimated.
.project_fit_inputs <- c(
  "model_type", "rasch_calibration", "thr_structure", "thr_mode",
  "exp_level", "bt_thr", "bt_ties", "anchor_type", "lp_layout",
  "lp_structure", "id_col", "ef_id", "ef_group", "lp_person",
  "lp_item", "lp_score", "lp_interaction", "bt_a", "bt_b", "bt_win",
  "bt_judge", "bt_count", "pc_rank", "ng", "ef_se", "ef_workers",
  "factor_cols", "item_cols", "ef_items", "lp_facets", "lp_items_wide",
  "bt_margin", "bt_response", "bt_order", "exp_main",
  "exp_interactions", "ef_prefix", "bt_position", "ng_auto", "maxit",
  "tol", "ef_reps", "ef_seed")

.collect_app_settings <- function(input) {
  ids <- unique(c(.project_radio_inputs, .project_select_inputs,
                  .project_selectize_inputs, .project_checkbox_inputs,
                  .project_numeric_inputs))
  current <- shiny::reactiveValuesToList(input, all.names = TRUE)
  # Predictor-type controls are generated from the uploaded metadata and do
  # not have fixed IDs. Retain their type, reference and ordinal-order values.
  dynamic <- grep("^exp_(type|ref|order)_[0-9]+$", names(current), value = TRUE)
  ids <- unique(c(ids, dynamic))
  out <- current[intersect(ids, names(current))]
  out[!vapply(out, is.null, logical(1))]
}

.collect_fit_settings <- function(input) {
  settings <- .collect_app_settings(input)
  dynamic <- grep("^exp_(type|ref|order)_[0-9]+$", names(settings),
                  value = TRUE)
  settings[intersect(names(settings), c(.project_fit_inputs, dynamic))]
}

.app_fit_source <- function(fit) {
  source <- attr(fit, "rasch_app_source", exact = TRUE)
  if (is.list(source) && is.data.frame(source$data) &&
      is.list(source$settings) && is.list(source$resources) &&
      is.list(source$simulation)) source else NULL
}

# Frame estimation is a downstream extension of the fitted CJ calibration.
# Resolve its data roles from that calibration's frozen run specification,
# not from sidebar controls that may have changed since it was fitted.
.app_btl_frame_source <- function(fit, current_data) {
  if (!inherits(fit, "rasch_btl"))
    stop("`fit` must be a Comparative Judgement calibration", call. = FALSE)
  source <- .app_fit_source(fit)
  if (is.null(source))
    stop(paste("the Comparative Judgement calibration predates saved run",
               "metadata; refit it before adding frames"), call. = FALSE)
  current <- as.data.frame(current_data, check.names = FALSE)
  if (!identical(source$data, current))
    stop(paste("the source data have changed since the Comparative Judgement",
               "calibration was fitted; refit it before adding frames"),
         call. = FALSE)

  settings <- source$settings
  column <- function(id, optional = FALSE) {
    value <- settings[[id]]
    absent <- is.null(value) || identical(value, NONE) ||
      (is.character(value) && length(value) == 1L && !is.na(value) &&
         !nzchar(value))
    if (optional && absent) return(NULL)
    if (length(value) != 1L || !is.character(value) || is.na(value) ||
        !nzchar(value) || !value %in% names(source$data))
      stop(sprintf("the fitted calibration has no valid `%s` data role", id),
           call. = FALSE)
    value
  }
  ties <- settings$bt_ties %||% "drop"
  if (length(ties) != 1L || !is.character(ties) || is.na(ties) ||
      !ties %in% c("drop", "half"))
    stop("the fitted calibration has no valid tie rule", call. = FALSE)

  list(data = source$data,
       object_a = column("bt_a"), object_b = column("bt_b"),
       winner = column("bt_win", optional = TRUE),
       judge = column("bt_judge", optional = TRUE),
       count = column("bt_count", optional = TRUE),
       response = column("bt_response", optional = TRUE),
       margin = column("bt_margin", optional = TRUE), ties = ties)
}

# Judge-group DIF is a downstream calculation on a fitted comparison design.
# A changed sidebar column is safe only when it carries exactly the same judge
# IDs, in the same row order, as the role used by the calibration.  This keeps
# equivalent metadata aliases usable while refusing a permuted or replacement
# judge role that would silently regroup the old comparisons.
.app_btl_dif_judge_role <- function(source, current_data, current_col) {
  if (is.null(source))
    stop(paste("the Comparative Judgement calibration predates saved run",
               "metadata; refit it before running judge-group DIF"),
         call. = FALSE)
  settings <- source$settings
  fitted_col <- settings[["bt_judge"]]
  absent <- is.null(fitted_col) || identical(fitted_col, NONE) ||
    (is.character(fitted_col) && length(fitted_col) == 1L &&
       !is.na(fitted_col) && !nzchar(fitted_col))
  if (absent || length(fitted_col) != 1L || !is.character(fitted_col) ||
      is.na(fitted_col) || !nzchar(fitted_col) ||
      !fitted_col %in% names(source$data))
    stop("the fitted Comparative Judgement calibration has no valid judge role",
         call. = FALSE)
  if (length(current_col) != 1L || !is.character(current_col) ||
      is.na(current_col) || !nzchar(current_col) ||
      identical(current_col, NONE) || !current_col %in% names(current_data))
    stop("judge-group DIF needs the judge column nominated on the Data page",
         call. = FALSE)
  fitted_ids <- as.character(source$data[[fitted_col]])
  current_ids <- as.character(current_data[[current_col]])
  if (length(fitted_ids) != length(current_ids) ||
      !identical(fitted_ids, current_ids))
    stop(paste("the selected judge column does not reproduce the fitted judge",
               "IDs; choose an equivalent metadata column or refit the",
               "Comparative Judgement calibration before running DIF"),
         call. = FALSE)
  list(fitted_col = fitted_col, current_col = current_col,
       ids = current_ids)
}

.restore_app_settings <- function(session, settings) {
  if (!is.list(settings) || !length(settings)) return(invisible(NULL))
  apply_ids <- function(ids, fun, argument = "selected") {
    for (id in intersect(ids, names(settings))) {
      args <- list(session = session, inputId = id)
      args[[argument]] <- settings[[id]]
      do.call(fun, args)
    }
  }
  apply_ids(.project_radio_inputs, shiny::updateRadioButtons)
  apply_ids(.project_select_inputs, shiny::updateSelectInput)
  apply_ids(.project_selectize_inputs, shiny::updateSelectizeInput)
  apply_ids(.project_checkbox_inputs, shiny::updateCheckboxInput,
            argument = "value")
  apply_ids(.project_numeric_inputs, shiny::updateNumericInput,
            argument = "value")
  dynamic <- grep("^exp_(type|ref)_[0-9]+$", names(settings), value = TRUE)
  apply_ids(dynamic, shiny::updateSelectInput)
  orders <- grep("^exp_order_[0-9]+$", names(settings), value = TRUE)
  apply_ids(orders, shiny::updateTextInput, argument = "value")
  invisible(NULL)
}


.efrm_detected_cores <- if (requireNamespace("rasch", quietly = TRUE))
  getFromNamespace(".efrm_available_workers", "rasch")() else if (
    exists(".efrm_available_workers", mode = "function"))
  .efrm_available_workers() else 1L
.efrm_worker_values <- seq_len(max(1L, min(4L, .efrm_detected_cores)))
.efrm_worker_choices <- stats::setNames(
  .efrm_worker_values,
  paste(.efrm_worker_values,
        ifelse(.efrm_worker_values == 1L, "worker", "workers")))

.wrightmap_available <- requireNamespace("WrightMap", quietly = TRUE)
.wrightmap_item_panels <- .wrightmap_available &&
  "item.groups" %in% names(formals(WrightMap::wrightMap))

# Match the package's collision-safe crossed-factor cells. This matters when
# a level itself contains the display separator (for example a colon).
app_factor_cells <- function(x, sep = ":") {
  fun <- if (exists(".factor_cells", mode = "function", inherits = TRUE))
    get(".factor_cells", mode = "function", inherits = TRUE)
  else getFromNamespace(".factor_cells", "rasch")
  fun(x, sep = sep)
}

app_factor_keys <- function(x) {
  fun <- if (exists(".factor_keys", mode = "function", inherits = TRUE))
    get(".factor_keys", mode = "function", inherits = TRUE)
  else getFromNamespace(".factor_keys", "rasch")
  fun(x)
}

# sanitise a proportion-type input (alpha level, chance probability): fall
# back to the default outside the open interval (0, 1)
clamp01 <- function(x, default)
  if (is.null(x) || is.na(x) || x <= 0 || x >= 1) default else x

# p-values as text: "%.3f" alone prints a misleading 0.000 for tiny p
# write.csv encodes a small double in scientific notation, so a downloaded
# table would disagree with the screen beside it. Raise scipen for the call:
# the file keeps full precision and loses the exponent.
# A raw preview is the one table that does not go through num_dt, so its
# doubles reach DataTables as serialised JSON and a small one renders in
# scientific notation. Round the non-integer numerics the same way.
round_preview <- function(dt, d) {
  num <- names(d)[vapply(d, function(v)
    is.numeric(v) && !all(is.na(v) | v == round(v)), TRUE)]
  if (length(num)) DT::formatRound(dt, num, 3) else dt
}

# The package refuses some analyses on some designs deliberately -- residual
# components on structurally disjoint columns, differential functioning on the
# factor that defines a frame -- and signals those with a rasch_refusal
# condition. A refusal is a result and belongs in the app's own validation
# voice; anything else is a fault and should still arrive in red.
soft <- function(expr) {
  tryCatch(expr,
           rasch_refusal = function(e) validate(need(FALSE, conditionMessage(e))))
}

write_csv_plain <- function(d, file) {
  op <- options(scipen = 999)
  on.exit(options(op), add = TRUE)
  write.csv(d, file, row.names = FALSE)
}

fmt_p <- function(p)
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))

# value-box guard: NULL, NA, NaN, and Inf must never reach a conditional or
# a sprintf; such values display as an em dash on a neutral theme
finite1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

inference_count <- function(p, alpha = 0.05) {
  tested <- sum(is.finite(p)); total <- length(p)
  flagged <- if (tested) sum(p[is.finite(p)] < alpha) else NA_integer_
  unavailable <- total - tested
  text <- if (!tested) "Unavailable"
    else if (!unavailable) sprintf("%d of %d", flagged, tested)
    else sprintf("%d of %d tested (%d unavailable)",
                 flagged, tested, unavailable)
  list(flagged = flagged, tested = tested, unavailable = unavailable,
       total = total, text = text)
}

# p-values phrased for prose ("p = 0.018" / "p < 0.001"), matching fmt_p's
# thresholds; used by the curated stat rows
p_lab <- function(p) {
  if (!finite1(p)) "p = NA"
  else if (p < 0.001) "p < 0.001"
  else sprintf("p = %.3f", p)
}

# curated stat-box rows (Summary page): one label-value line per statistic
stat_row <- function(label, value)
  div(class = "stat-row",
      span(class = "stat-label", label),
      span(class = "stat-value", value))
stat_rows <- function(...) div(class = "stat-rows", ...)

# Compact headline metrics.  These replace the large saturated value boxes:
# context remains visible, but results no longer consume the first several
# phone screens or give every statistic the visual weight of a warning.
metric_tile <- function(id, label, value, detail = NULL, icon = NULL,
                        status = c("neutral", "good", "bad", "warning",
                                   "accent", "person", "item")) {
  status <- match.arg(status)
  div(class = paste("metric-tile", paste0("metric-", status)),
    div(class = "metric-heading",
      if (!is.null(icon)) div(class = "metric-icon", glyph(icon)),
      span(class = "metric-label", label),
      info_icon(app_help(id), paste("About", label))),
    div(class = "metric-value", value),
    if (!is.null(detail)) div(class = "metric-detail", detail))
}
metric_grid <- function(...) div(class = "metric-grid mb-3", ...)

# measurement-themed value-box glyphs: inline SVG stroked with currentColor,
# so each glyph inherits its box's text colour in light and dark themes
.glyph_body <- list(
  # bell curve: the person distribution
  distribution = '<path d="M2 20 C7 20 8.5 5.5 12 5.5 S17 20 22 20"/>',
  # logit scale with alternating major/minor ticks
  ruler = '<path d="M2 15 H22"/><path d="M5 15 V9"/><path d="M9.5 15 V11.5"/><path d="M14 15 V9"/><path d="M18.5 15 V11.5"/>',
  # two distinct person distributions: separation / reliability
  separation = '<path d="M1 20 C4.5 20 5 10 7.5 10 S10.5 20 14 20"/><path d="M10 20 C13.5 20 14 5.5 16.5 5.5 S19.5 20 23 20"/>',
  alpha = '<text x="12" y="17.5" text-anchor="middle" font-size="17" font-style="italic" font-family="Georgia, serif" fill="currentColor" stroke="none">&#945;</text>',
  chisq = '<text x="11" y="17" text-anchor="middle" font-size="14" font-style="italic" font-family="Georgia, serif" fill="currentColor" stroke="none">&#967;&#178;</text>',
  # magnifier over a curve: person separation quality
  power = '<circle cx="10" cy="10" r="6.2"/><path d="M14.6 14.6 L20.5 20.5"/><path d="M6.8 11.5 Q10 6 13.2 11.5"/>',
  # data matrix rows / columns
  grid = '<rect x="3" y="4" width="18" height="16" rx="1.5"/><path d="M3 9.3 H21 M3 14.6 H21"/>',
  columns = '<rect x="3" y="4" width="18" height="16" rx="1.5"/><path d="M9 4 V20 M15 4 V20"/>',
  # 2x2 grid with a crossed-out cell: missing responses
  missing = '<rect x="3" y="4" width="18" height="16" rx="1.5"/><path d="M3 12 H21 M12 4 V20"/><path d="M14.5 14.5 L18.5 18.5 M18.5 14.5 L14.5 18.5" stroke-width="1.3"/>',
  # location span between two ends of the scale
  range = '<path d="M4 6 V18 M20 6 V18 M4 12 H20"/>',
  # a point off the trend: misfit
  outlier = '<path d="M3 17 C9 15.5 15 14.5 21 12.5"/><circle cx="16.5" cy="5.5" r="1.7" fill="currentColor" stroke="none"/>',
  # crossing solid/dashed step lines: reversed thresholds
  disorder = '<path d="M4 7 H10 L16 17 H20"/><path d="M4 17 H10 L16 7 H20" stroke-dasharray="2.6 2.2"/>',
  # two objects with a double-headed arrow: a paired comparison
  pair = '<circle cx="5.2" cy="12" r="2.7"/><circle cx="18.8" cy="12" r="2.7"/><path d="M9.3 12 H14.7 M11.2 10 L9.2 12 L11.2 14 M12.8 10 L14.8 12 L12.8 14"/>',
  # judge's balance
  balance = '<path d="M12 5 V19 M6 5 H18"/><path d="M6 5 L3.5 11 H8.5 Z"/><path d="M18 5 L15.5 11 H20.5 Z"/><path d="M9 19 H15"/>',
  # object locations as a podium
  podium = '<rect x="3.5" y="11" width="5" height="9" rx="1"/><rect x="9.5" y="6" width="5" height="14" rx="1"/><rect x="15.5" y="14" width="5" height="6" rx="1"/>')
glyph <- function(name)
  HTML(sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" style="width:100%%;height:100%%">%s</svg>',
               .glyph_body[[name]]))

# helpers for the "R code for this analysis" disclosure
qstr <- function(x) encodeString(as.character(x), quote = '"')
qvec <- function(x)
  if (length(x) == 1L) qstr(x) else
    paste0("c(", paste(qstr(x), collapse = ", "), ")")

theme <- bs_theme(
  version = 5, preset = "shiny",
  bg = "#f8fafc", fg = "#0f172a",
  primary = "#2563eb", secondary = "#64748b",
  success = "#0f766e", danger = "#dc2626", warning = "#d97706",
  base_font = font_collection("Inter", "Segoe UI", "Roboto",
                              "Helvetica Neue", "Arial", "sans-serif"),
  heading_font = font_collection("Inter", "Segoe UI", "Roboto",
                                 "Helvetica Neue", "Arial", "sans-serif"),
  code_font = font_collection("JetBrains Mono", "SFMono-Regular", "Consolas",
                              "Liberation Mono", "monospace"),
  "navbar-bg" = "#0f172a",
  "headings-font-weight" = "600",
  "font-size-base" = "0.925rem"
)

css <- HTML("
  .navbar-brand { font-weight: 700; letter-spacing: .02em; }
  /* two-line wordmark: the package name over the spelled-out tagline */
  .app-brand { display: inline-flex; flex-direction: column;
               justify-content: center; align-items: flex-start;
               text-align: left; line-height: 1.06; }
  .app-brand-name { font-weight: 700; font-size: 1.05rem; letter-spacing: .02em; }
  .app-brand-sub { font-weight: 400; font-size: .70rem; opacity: .6;
                   letter-spacing: .03em; }
  .card-header { font-weight: 600; }
  .value-box-title { font-size: .72rem; text-transform: uppercase; letter-spacing: .04em; white-space: nowrap; }
  .value-box-value { font-size: 1.45rem; }
  pre, .shiny-text-output { white-space: pre-wrap; font-size: .82rem; }
  .btn-xs { padding: .1rem .5rem; font-size: .75rem; }
  .form-label { font-weight: 600; font-size: .85rem; }
  /* APA-style tables: tabular numerals, no vertical rules, no zebra,
     a strong rule under the header row */
  table.dataTable { font-size: .85rem; width: auto !important; max-width: 100%; }
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
  /* curated stat boxes (Summary page): label-value rows with a hairline
     separator; tabular numerals keep the values aligned in both themes */
  .stat-rows { display: flex; flex-direction: column; }
  .stat-row { display: flex; justify-content: space-between; align-items: baseline; gap: 1rem; padding: .45rem 0; border-bottom: 1px solid var(--bs-border-color-translucent); }
  .stat-row:last-child { border-bottom: 0; }
  .stat-label { color: var(--bs-secondary-color); font-size: .85rem; }
  .stat-value { font-weight: 600; font-variant-numeric: tabular-nums; text-align: right; }
  .stat-head { color: var(--bs-secondary-color); font-size: .85rem; margin-bottom: .35rem; }
  .metric-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(145px, 1fr));
    gap: .65rem; }
  .metric-tile { min-width: 0; min-height: 88px; padding: .7rem .8rem;
    border: 1px solid var(--bs-border-color); border-top-width: 3px;
    border-radius: var(--bs-border-radius-lg); background: var(--bs-body-bg);
    box-shadow: 0 1px 2px rgba(15, 23, 42, .04); }
  .metric-neutral { border-top-color: var(--bs-border-color); }
  .metric-accent { border-top-color: var(--bs-primary); }
  /* the diagnostic pair: green when a check passes, red when it wants
     attention. Amber is not used for it, because amber now marks the item
     side of the model and the two readings would be indistinguishable. */
  .metric-good { border-top-color: var(--bs-success); }
  .metric-bad { border-top-color: var(--bs-danger); }
  .metric-warning { border-top-color: var(--bs-warning); }
  /* the two sides of the model, in the colours the plots already use for
     them: persons blue, items amber (see the Wright and person-item maps) */
  .metric-person { border-top-color: #2563eb; }
  .metric-item { border-top-color: #f59e0b; }
  .metric-heading { display: flex; align-items: center; min-width: 0;
    color: var(--bs-secondary-color); }
  .metric-icon { width: 1.15rem; height: 1.15rem; flex: 0 0 auto; margin-right: .35rem; }
  .metric-label { min-width: 0; overflow: hidden; text-overflow: ellipsis;
    white-space: nowrap; font-size: .73rem; font-weight: 650;
    text-transform: uppercase; letter-spacing: .035em; }
  .metric-value { margin-top: .28rem; font-size: 1.35rem; line-height: 1.1;
    font-weight: 650; font-variant-numeric: tabular-nums; }
  .metric-detail { margin-top: .2rem; color: var(--bs-secondary-color);
    font-size: .72rem; line-height: 1.25; white-space: nowrap;
    overflow: hidden; text-overflow: ellipsis; }
  @media (max-width: 575.98px) {
    .metric-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); gap: .5rem; }
    .metric-tile { min-height: 88px; padding: .6rem .65rem; }
    .metric-heading { align-items: flex-start; }
    .metric-label, .metric-detail { white-space: normal; overflow: visible;
      text-overflow: clip; line-height: 1.2; }
    .metric-value { font-size: 1.2rem; }
  }
  /* card headers: title left, action chips right, with a small gap
     between chips (the chips row right-aligns even when there is no
     title, via margin-left:auto) */
  .rasch-card-header { display: flex; align-items: center; justify-content: space-between; gap: .5rem; }
  .rasch-chips { display: flex; align-items: center; flex-wrap: wrap; gap: .35rem; margin-left: auto; }
  .rasch-info-button { display: inline-flex; align-items: center; justify-content: center;
    width: 1.3rem; height: 1.3rem; border: 1px solid rgba(var(--bs-primary-rgb), .16);
    background: var(--bs-primary-bg-subtle); color: var(--bs-primary);
    padding: 0; margin-left: .25rem; border-radius: 50%; line-height: 1; }
  .rasch-info-button:hover, .rasch-info-button:focus-visible {
    color: #fff; background: var(--bs-primary); border-color: var(--bs-primary);
    outline: none; box-shadow: 0 0 0 .18rem rgba(var(--bs-primary-rgb), .18); }
  .rasch-explainer { max-width: 22rem; font-size: .84rem; line-height: 1.42; }
  .rasch-accordion-info { display: flex; justify-content: flex-end;
    margin-bottom: .5rem; }
  .rasch-control-button { display: inline-flex; align-items: center; gap: .3rem; }
  .rasch-axis-popover { width: 250px; }
  .rasch-axis-popover .form-group, .rasch-axis-popover .shiny-input-container {
    width: 100%; margin-bottom: .65rem;
  }
  .rasch-axis-custom { display: grid; grid-template-columns: 1fr 1fr; gap: .5rem; }
  .rasch-axis-custom .form-group, .rasch-axis-custom .shiny-input-container {
    margin-bottom: 0;
  }
  .rasch-section-toolbar { display: flex; justify-content: flex-end;
    align-items: center; gap: .35rem; margin: .5rem 0; }
  .rasch-control-column, .rasch-result-column { min-width: 0; }
  .rasch-control-card { position: sticky; top: .75rem; }
  .rasch-control-card > .card-body { padding: .9rem; }
  @media (max-width: 991.98px) {
    .rasch-control-card { position: static; }
  }
  .rasch-plot-toolbar { display: flex; align-items: flex-end; flex-wrap: wrap;
    gap: .65rem; width: 100%; }
  .rasch-plot-toolbar .form-group,
  .rasch-plot-toolbar .shiny-input-container { margin-bottom: 0; }
  .rasch-plot-toolbar .rasch-icc-compare { flex: 1 1 240px; max-width: 320px; }
  .rasch-plot-toolbar .rasch-class-intervals { flex: 0 0 116px; }
  /* inline form controls that sit on a flex row with buttons: strip the
     bottom margin Shiny's containers carry */
  .rasch-inline-check .form-group, .rasch-inline-check .shiny-input-container,
  .rasch-inline-check .checkbox,
  .rasch-inline-select .form-group, .rasch-inline-select .shiny-input-container {
    margin-bottom: 0; width: auto;
  }
  .rasch-inline-select select.form-select { padding: .15rem 1.6rem .15rem .5rem; font-size: .78rem; }
  .rasch-predictor-types { display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: .5rem; }
  .rasch-predictor-type { min-width: 0; padding: .55rem .6rem;
    border: 1px solid var(--bs-border-color); border-radius: .5rem;
    background: var(--bs-tertiary-bg); }
  .rasch-predictor-type .form-group,
  .rasch-predictor-type .shiny-input-container { width: 100%; margin-bottom: .35rem; }
  .rasch-predictor-type .shiny-input-container:last-child { margin-bottom: 0; }
  /* collapsed advanced-settings disclosure inside the sidebar accordion */
  .rasch-advanced { margin-top: .5rem; }
  .rasch-advanced summary { cursor: pointer; font-size: .8rem; font-weight: 600; color: var(--bs-secondary-color); margin-bottom: .35rem; }
  .empty-state { text-align: center; padding: 3rem 1rem; color: var(--bs-secondary-color); }
  .shiny-output-error-validation {
    text-align: center; padding: 2rem 1rem; color: var(--bs-secondary-color);
  }
  .rasch-nav-summary { display: inline-flex; align-items: center; gap: .45rem;
    padding: .28rem .6rem; border: 1px solid rgba(255,255,255,.18);
    border-radius: 999px; color: rgba(255,255,255,.78); font-size: .76rem;
    line-height: 1; white-space: nowrap; }
  .rasch-nav-model { color: #fff; font-weight: 650; }
  .rasch-nav-sep { opacity: .35; }
  .rasch-nav-good { color: #99f6e4; }
  .rasch-nav-warn { color: #fde68a; }
  @media (max-width: 991.98px) {
    .rasch-nav-summary .rasch-nav-secondary { display: none; }
  }
  .analysis-pipeline { color: var(--bs-body-color); }
  .analysis-pipeline-steps { display: flex; align-items: center; flex-wrap: wrap;
    gap: .25rem; }
  .analysis-pipeline-arrow { width: .7rem; height: .7rem;
    color: var(--bs-secondary-color); }
  /* the DT bottom elements (info + pager) float; without clearance they
     collide with the collapsed R-code footer and can render outside the
     card. Clear the footer, give the wrapper self-clearing bottom room. */
  div.dataTables_wrapper { padding-bottom: .25rem; }
  div.dataTables_wrapper::after { content: ''; display: block; clear: both; }
  .rcode { clear: both; margin-top: .5rem; }
  .rcode summary { cursor: pointer; font-size: .78rem; color: var(--bs-secondary-color); }
  .rcode pre { font-size: .78rem; background: var(--bs-tertiary-bg); border: 1px solid var(--bs-border-color); border-radius: 6px; padding: .5rem .75rem; margin: .35rem 0 0; white-space: pre-wrap; }
  .rcode-copy { float: right; margin-top: .35rem; }
  /* floating hover-identification tooltip on the person-fit and item-fit-map
     scatter plots: a small dark pill that follows the cursor, positioned in
     CSS pixels from nearPoints()/coords_css. Colours are fixed (not themed)
     so it reads the same over the plot's white background in both app
     themes; pointer-events: none keeps it from stealing the hover it reports. */
  .rasch-hover-tip {
    position: absolute; pointer-events: none; z-index: 20;
    background: rgba(33, 37, 41, .85); color: #fff;
    font-size: .78rem; padding: .25rem .5rem; border-radius: 4px;
    white-space: nowrap; transform: translate(-50%, -100%);
  }
  /* Every fill item is allowed to shrink inside its grid cell. This prevents
     one wide table or plot from widening the entire page. */
  html, body { max-width: 100%; overflow-x: hidden; }
  .card, .card-body, .bslib-grid, .bslib-sidebar-layout,
  .main, .tab-content, .tab-pane { min-width: 0; max-width: 100%; }
  .shiny-plot-output, .shiny-image-output { width: 100% !important;
    max-width: 100%; }
  div.dataTables_scroll, div.dataTables_scrollBody { max-width: 100%; }
  @media (max-width: 575.98px) {
    .card-body { padding: .65rem !important; }
    .rasch-card-header { align-items: flex-start; flex-wrap: wrap; }
    .rasch-chips { width: 100%; justify-content: flex-end; }
    .form-label { font-size: .8rem; }
    table.dataTable { font-size: .78rem; }
    .shiny-plot-output { height: min(68vh, 520px) !important;
      min-height: 300px; }
    .empty-state { padding: 1.75rem .65rem; }
    .rasch-hover-tip { max-width: 82vw; overflow: hidden;
      text-overflow: ellipsis; }
  }
")

# collapsed per-output "R code" footer (jamovi-style syntax mode): shows the
# exact rasch call reproducing the output, updating with the current selections;
# the server registers a matching `<id>_code` renderText for every output
rcode_details <- function(id)
  tags$details(class = "rcode",
    tags$summary(bs_icon("code-slash"), " R code"),
    tags$button(class = "btn btn-outline-secondary btn-xs rcode-copy",
                type = "button",
                onclick = sprintf(
                  "navigator.clipboard.writeText(document.getElementById('%s_code').innerText)",
                  id),
                "Copy"),
    verbatimTextOutput(paste0(id, "_code"), placeholder = FALSE))

# Accessible information popover used throughout the app.  A real button gives
# keyboard and touch users the same access as mouse users; the concise copy is
# hidden until requested so the initial interface remains quiet.
info_icon <- function(info, label = "About this result") {
  if (is.null(info) || !length(info) || !nzchar(info)) return(NULL)
  popover(
    tags$button(type = "button", class = "rasch-info-button",
                `aria-label` = label,
                bs_icon("info-circle")),
    div(class = "rasch-explainer", info),
    placement = "auto"
  )
}

info_label <- function(label, info) {
  span(label, info_icon(info, paste("About", tolower(label))))
}

# Accordion headings are buttons. Their explainer therefore belongs inside
# the opened panel rather than inside the heading, where it would create an
# invalid button-within-button and break the surrounding navigation markup.
accordion_info <- function(info, label = "About this section")
  div(class = "rasch-accordion-info", info_icon(info, label))

# Compact axis settings. Presets cover ordinary use; exact limits remain
# available without occupying the analysis page. Automatic ranges are derived
# from the fitted person and threshold distributions by the server.
axis_control <- function(id, standard = c(-5, 5), wide = c(-8, 8),
                         label = "Axes") {
  mode_id <- paste0(id, "_mode")
  popover(
    tags$button(type = "button",
      class = "btn btn-outline-secondary btn-sm rasch-control-button",
      bs_icon("sliders"), label,
      `aria-label` = paste("Change", tolower(label))),
    div(class = "rasch-axis-popover",
      selectInput(mode_id, "Horizontal range",
        choices = stats::setNames(
          c("auto", "standard", "wide", "custom"),
          c("Automatic",
            sprintf("Standard (%g to %g)", standard[1], standard[2]),
            sprintf("Wide (%g to %g)", wide[1], wide[2]), "Custom")),
        selected = "auto", width = "100%"),
      conditionalPanel(
        sprintf("input['%s'] === 'custom'", mode_id),
        div(class = "rasch-axis-custom",
          numericInput(paste0(id, "_min"), "Minimum", standard[1], step = .5),
          numericInput(paste0(id, "_max"), "Maximum", standard[2], step = .5)))
    ),
    title = "Plot axes", placement = "auto"
  )
}

# A responsive controls-and-results layout that does not depend on bslib's
# collapsible sidebar JavaScript. Controls stack above results on small
# windows and remain visible beside them on larger screens.
app_layout_sidebar <- function(sidebar, ...) {
  controls <- do.call(tagList, sidebar$children)
  layout_columns(
    div(class = "rasch-control-column",
        card(class = "rasch-control-card", card_body(controls))),
    div(class = "rasch-result-column", ...),
    col_widths = breakpoints(sm = 12, lg = c(4, 8), xl = c(3, 9)),
    gap = "1rem")
}

# card header as a full-width flex bar: title (when given) on the left,
# action chips right-aligned with a small gap between them. Cards that sit
# inside an accordion panel that already names them pass title = NULL: the
# header then renders buttons-only (the info icon stays in the chips row).
card_header_bar <- function(title = NULL, buttons = NULL, info = NULL)
  card_header(class = "rasch-card-header",
    if (!is.null(title)) span(title, if (!is.null(info)) info_icon(info)),
    div(class = "rasch-chips",
        if (is.null(title) && !is.null(info)) info_icon(info),
        buttons))

# Card with a plot and PNG/PDF download buttons in the header. The body is
# non-fillable so flex sizing can never compress the fixed-height plot
# (the cause of the squashed plots), and a percentage height is avoided
# because it races the layout and renders a zero-height device.
# data-bs-theme is pinned to light because base plots draw on white.
# `info` adds a header tooltip; `controls` takes plot-display inputs rendered
# in a footer below the plot; `extra` takes further header buttons (e.g.
# the batch all-persons downloads on the kidmap card). title = NULL
# renders a buttons-only header (for cards inside named accordion panels).
# `hover = TRUE` opts a single-panel base-graphics plot into point
# identification: the plotOutput records hover coordinates (throttled) into
# input$<id>_hover, and the plot is wrapped in a position:relative div that
# floats a uiOutput(<id>_tip) tooltip (rendered server-side via nearPoints())
# over it. Left FALSE (the default) for every other plot card.
plotCard <- function(id, title = NULL, height = "560px", info = NULL,
                     controls = NULL, extra = NULL, hover = FALSE,
                     controls_top = FALSE) {
  info <- app_help(id, info)
  plot_out <- if (isTRUE(hover))
    div(style = "position: relative",
        plotOutput(id, height = height,
                   hover = hoverOpts(paste0(id, "_hover"), delay = 60,
                                     delayType = "throttle")),
        uiOutput(paste0(id, "_tip")))
  else plotOutput(id, height = height)
  # controls_top puts the toolbar between the header and the plot. The
  # footer is right for display tweaks read after the figure; a control that
  # CHANGES what is drawn -- and reveals further controls when switched --
  # belongs above a tall plot, or its consequences appear out of sight.
  card(
    id = paste0("card_", id),
    full_screen = TRUE,
    `data-bs-theme` = "light",
    card_header_bar(title, info = info, buttons = tagList(
      downloadButton(paste0(id, "_png"), "PNG", class = "btn-outline-secondary btn-xs"),
      downloadButton(paste0(id, "_pdf"), "PDF", class = "btn-outline-secondary btn-xs"),
      extra)),
    if (!is.null(controls) && isTRUE(controls_top))
      card_body(class = "rasch-plot-toolbar py-2", controls,
                padding = 8, fillable = FALSE),
    card_body(plot_out, rcode_details(id),
              padding = 8, fillable = FALSE),
    if (!is.null(controls) && !isTRUE(controls_top))
      card_footer(class = "rasch-plot-toolbar", controls)
  )
}

# `info` adds a header tooltip defining the key statistic; `footer` takes a
# small UI slot rendered under the table (dynamic interpretation notes);
# title = NULL renders a buttons-only header
tableCard <- function(id, title = NULL, note = NULL, info = NULL,
                      footer = NULL, controls = NULL) {
  info <- app_help(id, info)
  card(
    id = paste0("card_", id),
    full_screen = TRUE,
    card_header_bar(title, info = info, buttons = tagList(
      controls,
      downloadButton(paste0(id, "_csv"), "CSV",
                     class = "btn-outline-secondary btn-xs"))),
    # non-fillable body: DT outputs are fill items, and flex sizing inside a
    # natural-height card crops the last rows and the info/pager strip
    card_body(if (!is.null(note)) p(class = "text-muted small mb-2", note),
              DTOutput(id),
              if (!is.null(footer)) div(class = "table-note mt-2", footer),
              rcode_details(id),
              padding = 12, fillable = FALSE)
  )
}


# curated stat-box card: the body is a uiOutput of label-value rows built by
# the server; the CSV chip downloads the COMPLETE summary table (never the
# curated display), and the code footer names the call that builds it
statCard <- function(id, title = NULL, info = NULL, footer = NULL) {
  info <- app_help(id, info)
  card(
    id = paste0("card_", id),
    full_screen = TRUE,
    card_header_bar(title, info = info,
      buttons = downloadButton(paste0(id, "_csv"), "CSV",
                               class = "btn-outline-secondary btn-xs")),
    card_body(uiOutput(id),
              if (!is.null(footer)) div(class = "table-note mt-2", footer),
              rcode_details(id),
              padding = 12, fillable = FALSE)
  )
}

# bslib otherwise gives a full-screen nav card a random four-digit DOM id.
# A large app can then contain two cards with the same id, so retain the nav
# input id and give the surrounding card its own deterministic id.
stable_navset_card_underline <- function(..., id, header = NULL, footer = NULL,
                                         full_screen = FALSE) {
  out <- navset_card_underline(..., id = id, header = header, footer = footer,
                               full_screen = full_screen)
  if (isTRUE(full_screen)) out$attribs$id <- paste0("card_", id)
  out
}

# bslib assigns accordions without an explicit id, and each accordion panel,
# random four-digit DOM ids. With this many panels, collisions are plausible
# and can link one heading to another panel's body. Keep bslib's markup while
# replacing those ids deterministically.
.accordion_counter <- 0L
accordion <- function(..., id = NULL, open = NULL, multiple = TRUE,
                      class = NULL, width = NULL, height = NULL) {
  if (is.null(id)) {
    .accordion_counter <<- .accordion_counter + 1L
    id <- sprintf("rasch-accordion-%03d", .accordion_counter)
  }
  bslib::accordion(..., id = id, open = open, multiple = multiple,
                   class = class, width = width, height = height)
}

.accordion_panel_counter <- 0L
accordion_panel <- function(title, ..., value = title, icon = NULL) {
  .accordion_panel_counter <<- .accordion_panel_counter + 1L
  id <- sprintf("rasch-accordion-panel-%03d", .accordion_panel_counter)
  out <- bslib::accordion_panel(title, ..., value = value, icon = icon)
  button <- out$children[[1L]]$children[[1L]]
  button$attribs[["data-bs-target"]] <- paste0("#", id)
  button$attribs[["aria-controls"]] <- id
  out$children[[1L]]$children[[1L]] <- button
  out$children[[2L]]$attribs$id <- id
  out
}

# compact header switch revealing every column of a curated table
cols_switch <- function(id)
  div(class = "small text-secondary",
      input_switch(id, "All columns", value = FALSE))

# card header with an info-circle tooltip (for non-table cards)
info_header <- function(title, info)
  card_header(span(title, info_icon(info, paste("About", title))))

# ---------------------------------------------------------------------------
# Panels are built as objects and assembled into the workflow-ordered navbar
# (with Independence / Invariance / More menus) at the end of the UI section.
# ----------------------------------------------------------------- DATA --
panel_data <- nav_panel("Data", value = "p_data", icon = bs_icon("database"),
    app_layout_sidebar(
      sidebar = sidebar(width = 330,
        h6("Data source"),
        fileInput("file", NULL, accept = c(".csv", ".txt", ".tsv"),
                  buttonLabel = "Browse…", placeholder = "CSV / TSV file"),
        selectInput("demo_choice", "Example dataset",
                    c("None" = "none",
                      "Multiple choice, dichotomous" = "dich",
                      "Polytomous (PCM)" = "pcm",
                      "Rating scale (RSM)" = "rsm",
                      "Multiple Ratings (MFRM)" = "mfrm",
                      "Extended Frames (EFRM)" = "efrm",
                      "Comparative Judgement" = "btl")),
        accordion(
          id = "run_settings", multiple = TRUE,
          open = c("Data roles", "Model"),
          accordion_panel("Model", icon = bs_icon("diagram-2"),
            radioButtons("model_type", NULL,
                         c("Rasch" = "rasch",
                           "Multiple Ratings (MFRM)" = "mfrm",
                           "Extended Frames (EFRM)" = "efrm",
                           "Comparative Judgement" = "btl"))),
          accordion_panel("Data roles", icon = bs_icon("table"),
            conditionalPanel("input.model_type == 'rasch'",
              selectInput("id_col", "ID variable", NONE_CH),
              selectizeInput("factor_cols", "Person factors (DIF groups)", NULL,
                             multiple = TRUE,
                             options = list(placeholder = "none")),
              selectizeInput("item_cols", "Item columns", NULL, multiple = TRUE,
                             options = list(placeholder = "all remaining"))
            ),
            conditionalPanel("input.model_type == 'efrm'",
              selectInput("ef_id", "ID variable", NONE_CH),
              selectInput("ef_group",
                          span("Person group column",
                               info_icon("None treats all persons as one group, so the units differ by item set only.")),
                          NONE_CH),
              selectizeInput("ef_items", "Item columns", NULL, multiple = TRUE,
                             options = list(placeholder = "all remaining")),
              fileInput("ef_sets", info_label("Item-set map (CSV: item,set)",
                        paste("Each item-set by group cell is a frame with its own unit.",
                              "Group units use person-free pairwise comparisons;",
                              "set units use persons common to the sets.")),
                        accept = ".csv", placeholder = "optional"),
              checkboxInput("ef_prefix", "Infer sets from item-name prefix", TRUE)
            ),
            conditionalPanel("input.model_type == 'mfrm'",
              radioButtons("lp_layout", info_label("Data layout",
                           paste("Wide data have one row per person-by-facet",
                                 "combination and one column per item. Long data",
                                 "have person, item and score columns.")),
                           c("Items in columns (wide)" = "wide",
                             "One score per row (long)" = "long")),
              selectInput("lp_person", "Person column", NONE_CH),
              conditionalPanel("input.lp_layout == 'long'",
                selectInput("lp_item", "Item column", NONE_CH),
                selectInput("lp_score", "Score column", NONE_CH)),
              selectizeInput("lp_facets", info_label("Facet columns (e.g. rater)",
                             "Items and facet levels are calibrated jointly; facet severities carry standard errors and fit statistics."), NULL,
                             multiple = TRUE,
                             options = list(placeholder = "choose at least one")),
              conditionalPanel("input.lp_layout == 'wide'",
                selectizeInput("lp_items_wide", "Item columns", NULL,
                               multiple = TRUE,
                               options = list(placeholder = "all remaining"))),
              radioButtons("lp_structure", info_label("Facet structure",
                           paste("Additive estimates one severity per facet level.",
                                 "Interactive also estimates item-by-facet effects.")),
                           c("Additive" = "additive",
                             "Interactive (item-by-facet)" = "interactive")),
              conditionalPanel("input.lp_structure == 'interactive'",
                selectInput("lp_interaction", "Interacting facet", NULL))
            ),
            conditionalPanel("input.model_type == 'btl'",
              h6("One comparison per row"),
              selectInput("bt_a", "Object A column", NONE_CH),
              selectInput("bt_b", "Object B column", NONE_CH),
              conditionalPanel("!input.bt_response",
                selectInput("bt_win", "Winner column", NONE_CH),
                selectizeInput("bt_margin",
                               span("Margin of win (optional)",
                                    info_icon("Extent of the win (a little / much) as an ordered factor or increasing values. \"tie\"/\"draw\" in the winner column marks a tie; other unmatched values drop the row.")),
                               NULL,
                               options = list(placeholder = "none — dichotomous"))),
              selectizeInput("bt_response",
                             span("Polytomous response (optional)",
                                  info_icon(paste("Ordered preference for object A",
                                    "(worst to best, or scores 0..m); overrides the winner column.",
                                    "Ties belong in a middle category. A winner column fits",
                                    "the Bradley-Terry-Luce model; this fits its adjacent-categories extension."))),
                             NULL,
                             options = list(placeholder = "none — use winner")),
              selectInput("bt_judge",
                          span("Judge column (optional)",
                               info_icon("Enables the judge fit table and clusters standard errors by judge.")),
                          NONE_CH),
              conditionalPanel("input.bt_judge && input.bt_judge != '(none)'",
                selectizeInput("bt_order",
                               span("Judgment order (optional)",
                                    info_icon("Each judge's judgment sequence (timestamps or ranks); enables the exposure and carry-over dependence analysis.")),
                               NULL,
                               options = list(placeholder = "none")),
                selectizeInput("bt_jfactors",
                               span("Judge factors (optional)",
                                    info_icon("Judge groupings (constant within judge) for DIF by judge group.")),
                               NULL, multiple = TRUE,
                               options = list(placeholder = "none"))),
              checkboxInput("bt_position",
                            span("First-position advantage",
                                 info_icon("Object A is the first-presented of each pair; estimates the positional advantage (Davidson and Beaver 1977).")),
                            FALSE),
              conditionalPanel("input.rasch_calibration != 'explanatory'",
                fileInput("bt_anchor_file",
                          span("Anchor objects (CSV: object, location)",
                               info_icon(paste("Holds the named objects at their given locations",
                                 "and estimates the rest around them. Values are treated as",
                                 "fixed, so uncertainty from the earlier calibration is not",
                                 "included."))),
                          accept = ".csv", placeholder = "optional")),
              conditionalPanel("!input.bt_response && !input.bt_margin",
                radioButtons("bt_ties", "Ties",
                             c("Drop" = "drop", "Half a win each" = "half")))
            )),
          accordion_panel("Estimation options", icon = bs_icon("gear"),
            conditionalPanel("input.model_type == 'rasch' || input.model_type == 'btl'",
              radioButtons("rasch_calibration", info_label("Calibration",
                           paste("Free calibration estimates item thresholds or object locations directly.",
                                 "Explanatory calibration expresses them as functions",
                                 "of observed characteristics.")),
                           c("Free" = "free", "Explanatory" = "explanatory")),
              conditionalPanel("input.model_type == 'rasch' && input.rasch_calibration != 'explanatory'",
                radioButtons("thr_structure", info_label("Threshold structure",
                             paste("Dichotomous items need no setting. The rating scale",
                                   "model requires equal maximum scores; lr_test() compares",
                                   "it with the partial credit model.")),
                             c("Partial credit (item-specific)" = "pcm",
                               "Rating scale (common across items)" = "rsm")),
                conditionalPanel("input.thr_structure == 'pcm'",
                  radioButtons("thr_mode", "Threshold estimation",
                               c("Free thresholds" = "free",
                                 "Principal components (Andrich)" = "pc")),
                  conditionalPanel("input.thr_mode == 'pc'",
                    selectInput("pc_rank", info_label("Components",
                                paste("Constrains thresholds to a polynomial trend",
                                      "across categories, which can stabilise sparse",
                                      "categories. This option cannot be combined with anchors.")),
                                c("Location only" = "1",
                                  "+ spread (equal spread)" = "2",
                                  "+ skewness" = "3",
                                  "+ kurtosis (full PC)" = "4"),
                                selected = "4")))),
              conditionalPanel("input.rasch_calibration == 'explanatory'",
                fileInput("exp_predictors", info_label("Predictor metadata (CSV)",
                          paste("Use one row per item or object.",
                                "Threshold-level Rasch files also identify each threshold.")),
                          accept = ".csv", placeholder = "names and predictor columns"),
                conditionalPanel("input.model_type == 'rasch'",
                  radioButtons("exp_level", info_label("Predictor level",
                               paste("Item predictors are repeated over that item's thresholds.",
                                     "Threshold predictors may differ within an item.")),
                               c("Item locations" = "item",
                                 "Individual thresholds" = "threshold"))),
                uiOutput("exp_predictor_types"),
                uiOutput("exp_formula_controls"))),
            conditionalPanel(
              "input.model_type == 'btl' && (input.bt_response || input.bt_margin)",
              radioButtons("bt_thr", info_label("Threshold structure",
                           paste("Principal components pool the symmetric thresholds",
                                 "to a spread component so sparse categories borrow",
                                 "strength; free estimation fits each threshold pair.")),
                           c("Free symmetric" = "free",
                             "Principal components (spread)" = "pc"))),
            checkboxInput("ng_auto", info_label("Automatic class intervals",
                          "Chooses the largest practical interval count while retaining about 50 non-extreme observations per interval."), TRUE),
            conditionalPanel("!input.ng_auto",
              selectInput("ng", "Class intervals", choices = 2:16,
                          selected = 8, selectize = FALSE)),
            conditionalPanel("input.model_type == 'efrm'",
              selectInput("ef_se", info_label("Standard errors",
                          paste("Hybrid combines the stage-wise covariance estimates.",
                                "The full person bootstrap repeats the complete fit and is slower.")),
                          c("Hybrid" = "hybrid",
                            "Full person bootstrap (slow)" = "bootstrap")),
              numericInput("ef_reps", "Bootstrap replicates", value = 300,
                           min = 50, step = 50),
              numericInput("ef_seed", info_label("Bootstrap seed",
                           "Reproduces the same resamples and results for any selected worker count."),
                           value = 1, min = 0, step = 1),
              selectInput("ef_workers", info_label("Parallel workers",
                          paste("Runs independent bootstrap replicates concurrently.",
                                "Defaults to four, or fewer when the system limit is lower.")),
                          choices = .efrm_worker_choices,
                          selected = max(.efrm_worker_values))),
            tags$details(class = "rasch-advanced",
              tags$summary("Advanced"),
              numericInput("maxit", "Maximum iterations", value = 60, min = 5, step = 5),
              numericInput("tol", "Convergence criterion", value = 1e-8,
                           min = 1e-12, step = 1e-8),
              conditionalPanel("input.model_type == 'btl'",
                selectInput("bt_count", info_label("Count column (optional)",
                            paste("A row may represent several identical comparisons.",
                                  "Counts greater than one cannot be used with judgment order,",
                                  "which needs one row per comparison.")), NONE_CH)))),
          conditionalPanel("input.model_type == 'rasch'",
            accordion_panel("Scoring & anchors", icon = bs_icon("key"),
              fileInput("key_file",
                        span("Scoring key (CSV)",
                             info_icon("Columns item,key for a multiple-choice key — use \"A/C\" for a double key — or item,option,score for polytomous option scoring.")),
                        accept = ".csv", placeholder = "optional"),
              conditionalPanel("input.rasch_calibration != 'explanatory'",
                fileInput("anchor_file", info_label("Anchors for equating (CSV: item,k,tau)",
                          paste("Every named item must be present. Individual anchoring fixes",
                                "each threshold; location anchoring fixes each item's mean",
                                "while its thresholds remain free; average anchoring (RUMM's",
                                "average item anchoring) estimates every item free and shifts",
                                "the calibration so the anchor items' mean location equals the",
                                "mean anchor value. Supplied values are treated as fixed; their",
                                "earlier calibration uncertainty is not included.")),
                          accept = ".csv", placeholder = "optional"),
                radioButtons("anchor_type", "Anchor as",
                             c("Individual thresholds" = "individual",
                               "Item locations" = "average",
                               "Average of the anchor set" = "set_average"),
                             inline = TRUE))))
        ),
        input_task_button("run", "Estimate", icon = bs_icon("play-fill"),
                          type = "primary", class = "w-100 btn-lg mt-2"),
        uiOutput("efrm_job_controls"),
        conditionalPanel("output.has_override",
          uiOutput("override_status"),
          layout_columns(col_widths = c(6, 6), gap = "0.35rem",
            actionButton("undo_override", "Undo last change",
                         class = "btn-outline-secondary w-100"),
            actionButton("reset_override", "Reset all changes",
                         class = "btn-outline-warning w-100"))),
        div(class = "text-center mt-2",
            info_icon(paste("Item parameters use pairwise conditional maximum",
                            "likelihood (Zwinderman 1995); the",
                            "principal-components option follows Andrich &",
                            "Luo (2003). Person locations use Warm weighted",
                            "likelihood."),
                      "About estimation"))
      ),
      uiOutput("data_main")
    )
  )

# -------------------------------------------------------------- SUMMARY --
# bottom row built by the server: the likelihood-ratio card only applies to
# a PCM fit whose items share a common maximum score, so it hides otherwise
.lr_card <- function()
  card(info_header("PCM and rating scale comparison",
         "Compares the partial credit model with the rating scale model. Use the adjusted composite-likelihood statistic for inference."),
    card_body(
      input_task_button("run_lr", "Run likelihood-ratio test",
                        type = "primary"),
      verbatimTextOutput("lr_txt"),
      rcode_details("lr")))

panel_summary <- nav_panel("Summary", value = "p_summary", icon = bs_icon("clipboard-data"),
    # Rasch fits (hidden while a paired-comparison fit is active)
    conditionalPanel("output.is_btl != true",
    uiOutput("vboxes"),
    rcode_details("vboxes"),
    # stat-box cards sit inside plain divs: the grid row would otherwise
    # stretch them to equal height and pad the shorter card mid-row
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      div(statCard("fitsum_tbl", "Test of fit",
        info = "Compares observed and expected scores over class intervals. Structural models report response-cell fit. The CSV also contains fit-residual moments and fit-location correlations.",
        footer = uiOutput("fitsum_notes"))),
      div(statCard("targeting_tbl", "Targeting & reliability",
        info = "Compares the calibration-threshold and person distributions and reports the applicable reliability indices."))
    ),
    # server-rendered: the likelihood-ratio card only when it applies
    uiOutput("summary_bottom"),
    accordion(id = "test_acc", open = FALSE, class = "mb-3",
      accordion_panel("Test characteristic curve", value = "test_tcc",
        plotCard("tcc")),
      accordion_panel("Test information", value = "test_tif",
        plotCard("tif")),
      accordion_panel("Guttman scalogram", value = "test_guttman",
        plotCard("guttman", height = "640px", hover = TRUE))),
    div(class = "rasch-plot-toolbar mb-3",
        axis_control("ts_axis", standard = c(-6, 6), wide = c(-8, 8)))),
    # paired-comparison (BTL) fits: the headline value boxes and the
    # test-of-fit summary table. This block and the override block below are
    # siblings of the Rasch-only panel above, not children: nested inside it,
    # `is_btl == true` contradicts the ancestor's `is_btl != true` and the
    # BTL summary can never render.
    conditionalPanel("output.is_btl == true",
      uiOutput("btl_boxes"),
      rcode_details("btl_boxes"),
      layout_columns(col_widths = breakpoints(sm = 12, xl = 6),
        div(statCard("btl_fitsum_tbl", "Test of fit",
          info = "The pairwise chi-square compares observed and expected responses for each object pair. The object separation index estimates the proportion of observed location variance not attributable to error.",
          footer = uiOutput("btl_fitsum_notes")))),
      conditionalPanel("output.can_btl_boot == true",
        div(class = "d-flex align-items-end gap-2 flex-wrap mb-3",
          input_task_button("btl_boot_run", "Bootstrap the fit statistics",
                            icon = bs_icon("arrow-repeat"),
                            class = "btn-outline-primary"),
          div(style = "width: 110px;",
              numericInput("btl_boot_B", "Replicates", value = 999,
                           min = 99, max = 9999, step = 100)),
          div(style = "width: 90px;",
              numericInput("btl_boot_seed", "Seed", value = 1, step = 1)),
          span(class = "pb-2",
            info_icon(paste("Generates outcomes on the fitted comparison design",
                            "and refits the model. It calibrates whole-model,",
                            "pair, object and judge fit. Use at least 999",
                            "replicates for a final analysis. Inference is",
                            "withheld if fewer than 90% of refits are usable."))),
          uiOutput("btl_boot_state", inline = TRUE),
          uiOutput("btl_boot_job_controls", inline = TRUE)))),
    conditionalPanel("output.has_override == true",
      accordion(id = "change_acc", open = FALSE,
        accordion_panel("Changes from the original fit",
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)),
            tableCard("change_est_tbl", "Item or object estimates",
              info = "Original and active common-scale locations. The active fit is the fit used by every downstream table and plot; added or removed parameters are retained so the transformation is explicit."),
            conditionalPanel("output.is_btl != true",
              tableCard("change_person_tbl", "Person estimates",
                info = "Original and active person locations for the same response rows. This shows how a DIF split or dependence superitem changes person measurement as well as the item calibration.")))))
  ))

# ---------------------------------------------------------------- ITEMS --
panel_items <- nav_panel("Items", value = "p_items", icon = bs_icon("list-check"),
    # Rasch fits (hidden while a paired-comparison fit is active)
    conditionalPanel("output.is_btl != true",
    uiOutput("items_vboxes"),
    rcode_details("items_vboxes"),
    # The asymptotic fit statistics grow anticonservative with the sample;
    # the bootstrap replaces their reference distributions with ones formed
    # by refitting the model to data it generated. Post-estimation and
    # costly (B refits), so it runs on request and clears on a new fit.
    conditionalPanel("output.can_boot == true",
      div(class = "d-flex align-items-end gap-2 flex-wrap mb-3",
        input_task_button("boot_run", "Bootstrap the fit statistics",
                          icon = bs_icon("arrow-repeat"),
                          class = "btn-outline-primary"),
        div(style = "width: 110px;",
            numericInput("boot_B", "Replicates", value = 999,
                         min = 99, max = 9999, step = 100)),
        div(style = "width: 90px;",
            numericInput("boot_seed", "Seed", value = 1, step = 1)),
        span(class = "pb-2",
          info_icon(paste("Refits data generated under the fitted model. It",
                          "calibrates item fit and score-conditional person fit;",
                          "adjusted probabilities appear in both tables. They",
                          "control each statistic family under the fitted global",
                          "null, not invariant items or persons when other parts",
                          "of the model misfit. Use at",
                          "least 999 replicates for a final analysis and more",
                          "when the test has many items. Inference is withheld",
                          "if fewer than 90% of refits are usable; sparse",
                          "polytomous categories can cause this."))),
        uiOutput("boot_state", inline = TRUE),
        uiOutput("boot_job_controls", inline = TRUE))),
    conditionalPanel("output.anchor_download_available == true",
      div(class = "mb-2 d-flex justify-content-end",
          downloadButton("dl_anchors", "Save anchors (CSV: item,k,tau)",
                         class = "btn-outline-secondary btn-sm"))),
    conditionalPanel("output.has_structural_items == true",
      div(class = "mb-3",
        tableCard("structural_items_tbl", "Common-scale item estimates",
          controls = cols_switch("structural_items_full"),
          info = paste("Item locations on the model's common scale.",
                       "Item-by-frame or item-by-facet response cells are",
                       "shown separately below because they describe fit",
                       "rather than additional items.")))),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      tableCard("items_tbl", uiOutput("items_table_title", inline = TRUE),
        controls = cols_switch("items_full"),
                "Click a row to explore that item on the right. Fit residual ~ N(0,1) under fit.",
                info = paste("For Extended Frames and Multiple Ratings, each",
                             "row is an observed item-by-frame or item-by-facet",
                             "cell; otherwise each row is an item. Marking uses",
                             "the displayed fit criteria as working guides."),
                footer = uiOutput("items_note")),
      stable_navset_card_underline(
        id = "items_nav",
        # Keep the tab strip quiet. The selected item and plot explanation sit
        # above the result; display settings and batch downloads sit below it.
        header = div(class = "d-flex align-items-center gap-3 flex-wrap",
          uiOutput("sel_item_title", inline = TRUE),
          uiOutput("drop_item_ui", inline = TRUE),
          conditionalPanel(
            "input.items_nav == 'ICC'",
            info_icon(app_help("icc"), "About this plot")),
          conditionalPanel("input.items_nav == 'Categories'",
                           info_icon(app_help("ccc"), "About this plot")),
          conditionalPanel("input.items_nav == 'Thresholds'",
                           info_icon(app_help("tpc"), "About this plot")),
          conditionalPanel("input.items_nav == 'Frequencies'",
                           info_icon(app_help("cfreq"), "About this plot"))),
        footer = div(class = "rasch-plot-toolbar",
          conditionalPanel(
            "input.items_nav == 'ICC'",
            div(class = "rasch-icc-compare",
              selectizeInput(
                "icc_compare_items",
                info_label("Compare items",
                  paste("Adds up to seven items to the same plot. The item",
                        "selected in the table remains the primary item.")),
                choices = NULL, multiple = TRUE,
                options = list(maxItems = 7,
                               placeholder = "Add up to seven items"),
                width = "100%"))),
          conditionalPanel(
            "input.items_nav != 'Frequencies' && input.items_nav != 'Chi-square'",
            div(class = "d-flex align-items-end gap-2 flex-wrap",
              div(class = "rasch-inline-check pb-1",
                  checkboxInput(
                    "show_obs",
                    info_label("Observed points",
                      paste("Adds class-interval observations to the model",
                            "curves. Clear it to show the model alone.")),
                    TRUE, width = "auto")),
              div(class = "rasch-class-intervals",
                  selectInput(
                    "ex_ng",
                    info_label("Class intervals",
                      paste("Sets the number of person-location groups used",
                            "for the observed points. It does not change the",
                            "model curve.")),
                    choices = 2:16, selected = 8, selectize = FALSE,
                    width = "100%")),
              axis_control("ex_axis", standard = c(-5, 5), wide = c(-8, 8)))),
          conditionalPanel("input.items_nav != 'Chi-square'",
            class = "ms-auto",
            div(class = "rasch-chips",
                downloadButton("items_all_pdf", "All (PDF)",
                               class = "btn-outline-secondary btn-xs"),
                downloadButton("items_all_zip", "All (PNG)",
                               class = "btn-outline-secondary btn-xs")))),
        full_screen = TRUE,
        nav_panel("ICC",
                  plotOutput("icc", height = "440px"),
                  div(class = "text-end",
                      downloadButton("icc_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("icc_pdf", "PDF", class = "btn-outline-secondary btn-xs")),
                  rcode_details("icc")),
        nav_panel("Categories",
                  plotOutput("ccc", height = "440px"),
                  div(class = "text-end",
                      downloadButton("ccc_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("ccc_pdf", "PDF", class = "btn-outline-secondary btn-xs")),
                  rcode_details("ccc")),
        nav_panel("Thresholds",
                  plotOutput("tpc", height = "440px"),
                  div(class = "text-end",
                      downloadButton("tpc_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("tpc_pdf", "PDF", class = "btn-outline-secondary btn-xs")),
                  rcode_details("tpc")),
        nav_panel("Frequencies",
                  plotOutput("cfreq", height = "440px"),
                  div(class = "text-end",
                      downloadButton("cfreq_png", "PNG", class = "btn-outline-secondary btn-xs"),
                      downloadButton("cfreq_pdf", "PDF", class = "btn-outline-secondary btn-xs")),
                  rcode_details("cfreq")),
        nav_panel("Chi-square",
                  uiOutput("chisq_caption"),
                  h6(span("Class intervals",
                          info_icon(app_help("chisq_int_tbl"),
                                    "About this table"))),
                  DTOutput("chisq_int_tbl"),
                  h6(span("Response categories by class interval",
                          info_icon(app_help("chisq_cat_tbl"),
                                    "About this table")), class = "mt-3"),
                  DTOutput("chisq_cat_tbl"),
                  div(class = "text-end mt-2",
                      downloadButton("chisq_int_csv", "Intervals CSV",
                                     class = "btn-outline-secondary btn-xs"),
                      downloadButton("chisq_cat_csv", "Categories CSV",
                                     class = "btn-outline-secondary btn-xs")),
                  rcode_details("chisq")))),
    accordion(id = "items_acc", open = FALSE, class = "mt-3 mb-3",
      accordion_panel("Threshold map", value = "items_thrmap",
        plotCard("thrmap")),
      accordion_panel("Item fit map", value = "items_imap",
        plotCard("imap", hover = TRUE)),
      accordion_panel("Fit residual distribution", value = "items_rdist",
        plotCard("rdist_i")),
      accordion_panel("Traditional statistics", value = "items_ctt",
        card(
          id = "card_ctt_tbl",
          full_screen = TRUE,
          card_header_bar(
            info = app_help("ctt_tbl"),
            buttons = downloadButton("ctt_tbl_csv", "CSV",
                                     class = "btn-outline-secondary btn-xs")),
          card_body(uiOutput("ctt_head"), DTOutput("ctt_tbl"),
                    rcode_details("ctt_tbl"),
                    padding = 12, fillable = FALSE)))),
    uiOutput("pc_comp_ui"),
    conditionalPanel("output.has_mc == true",
    layout_columns(col_widths = 12,
      tableCard("distractor_tbl", "Distractor analysis",
                "Locations use the rest measure; a distractor whose takers are abler than the keyed option's flags a possible miskey."),
      plotCard("distractor_plot", "Option curves"),
      card(card_header_bar("Polytomous option scoring (Andrich & Styles 2011)",
             info = app_help("rescore_tbl"),
             buttons = conditionalPanel("output.has_rescore == true",
               downloadButton("dl_rescore", "Key CSV",
                              class = "btn-outline-secondary btn-xs"))),
           card_body(fillable = FALSE,
             layout_columns(col_widths = c(3, 3, 3, 3),
               numericInput("rescore_min_n", "Min takers", 20, min = 5, step = 5),
               numericInput("rescore_z", "Separation z", 1.96, min = 0.5, step = 0.1),
               div(class = "mt-4",
                   input_task_button("rescore_go", "Propose option scores",
                                     type = "primary")),
               conditionalPanel("output.has_rescore == true", class = "mt-4",
                                cols_switch("rescore_full"))),
             conditionalPanel("output.has_rescore != true",
               p(class = "text-muted small mb-0",
                 "Select Propose option scores to see the rest-measure evidence and proposed key.")),
             conditionalPanel("output.has_rescore == true",
               DT::DTOutput("rescore_tbl"),
               rcode_details("rescore_tbl"))))))),
    # paired-comparison (BTL) fits: the object side of the analysis
    conditionalPanel("output.is_btl == true",
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
        tableCard("btl_obj_tbl", "Object locations and fit",
          controls = cols_switch("btl_full"),
          info = paste("Object locations and standard errors, with fit pooled",
                       "over each object's comparisons. After a fit bootstrap,",
                       "the adjusted probabilities refer each statistic to its",
                       "fitted global null across objects. Click a row to plot",
                       "the object.")),
        plotCard("btl_occ", "Object characteristic curve",
          info = "Expected response for the selected object over opponent location, with observed means for sufficiently observed opponents.",
          extra = downloadButton("btl_occ_all_pdf", "All (PDF)",
                                 class = "btn-outline-secondary btn-xs"))),
      accordion(id = "btl_items_acc", open = FALSE,
                class = "mt-3 mb-3",
        accordion_panel("Object caterpillar", value = "btl_caterpillar",
          plotCard("btl_plot")),
        # polytomous (ordinal) fits only: hidden entirely for dichotomous fits
        conditionalPanel("output.btl_graded == true",
          accordion_panel("Symmetric thresholds", value = "btl_thresholds",
            tableCard("btl_thr_tbl",
                      note = "Adjacent-categories thresholds of the polytomous structure, constrained symmetric (tau_k = -tau_(m+1-k)) so the model is invariant to presentation order.")),
          accordion_panel("Threshold components", value = "btl_components",
            tableCard("btl_comp_tbl",
                      note = "Spread is the linear component; the skewness component is structurally zero under presentation-order symmetry. Under the PC structure kurtosis is constrained to zero.")),
          accordion_panel("Category probability curves", value = "btl_catcurves",
            plotCard("btl_cats",
              info = "The probability of each polytomous response category as a function of the location difference between the two objects; the paired-comparison counterpart of a polytomous item's category curves."))),
        accordion_panel("Pairwise fit", value = "btl_pairs",
          tableCard("btl_pairs_tbl",
                    info = paste("Observed and expected win proportions, or mean",
                                 "ordered responses, for every object pair.",
                                 "The fit bootstrap adds fitted-null adjusted",
                                 "probabilities for the pair chi-squares.")))))
  )

# --------------------------------------------------------- EXPLANATORY --
panel_explanatory <- nav_panel("Explanatory", value = "p_explanatory",
    icon = bs_icon("bezier2"),
    uiOutput("expl_boxes"),
    rcode_details("expl_boxes"),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(5, 7)),
      tableCard("expl_test_tbl", "Model comparison",
        info = paste("Compares the active explanatory restrictions with the",
                     "free calibration of the same responses. Inference uses",
                     "the Kent adjustment for pairwise composite likelihood.",
                     "Calibration R-squared is the proportion of well-determined",
                     "free threshold or object variation reproduced by the",
                     "explanatory model.")),
      tableCard("expl_coef_tbl", "Predictor effects",
        info = paste("Estimated effects of the nominated item, threshold or object",
                     "characteristics in logits. Holm adjustment covers the",
                     "displayed coefficient family."))),
    plotCard("expl_calibration", "Explanatory and free calibration",
      info = paste("Each point compares an active explanatory threshold or object",
                   "location with its freely estimated counterpart. The diagonal denotes",
                   "agreement; labels identify parameters with the largest",
                   "departures.")),
    card(class = "mb-3",
      card_header_bar("Fixed-departure diagnostics",
        info = paste("Each candidate is added separately to the active model.",
                     "Item location moves all thresholds together; threshold",
                     "structure changes their relative positions; CJ uses fixed",
                     "object-location departures. Probabilities",
                     "use Kent calibration and Holm adjustment."),
        buttons = downloadButton("expl_diag_tbl_csv", "CSV",
                                 class = "btn-outline-secondary btn-xs")),
      card_body(fillable = FALSE,
        p(class = "text-muted small mb-2",
          "Select a row to approve that fixed departure and repeat the complete calibration."),
        DTOutput("expl_diag_tbl"),
        uiOutput("expl_selected_note"),
        div(class = "d-flex gap-2 justify-content-end mt-2",
          input_task_button("expl_relax", "Add selected departure and refit",
                            type = "primary")),
        rcode_details("expl_diag_tbl"))),
    conditionalPanel("output.has_expl_relaxations == true",
      tableCard("expl_relax_tbl", "Fixed departures",
        info = paste("The fixed departures included in the active calibration.",
                     "Every downstream table and plot uses this recalibrated",
                     "model; use Undo or Reset on the Data page to reverse it."))))

# -------------------------------------------------------------- PERSONS --
panel_persons <- nav_panel("Persons", value = "p_persons", icon = bs_icon("people"),
    # Rasch fits (hidden while a paired-comparison fit is active)
    conditionalPanel("output.is_btl != true",
    uiOutput("persons_vboxes"),
    rcode_details("persons_vboxes"),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      tableCard("person_tbl", "Person estimates",
          controls = cols_switch("persons_full"),
          info = paste("Warm WLE location and standard error, raw score and fit",
                       "statistics for each person. After a fit bootstrap, the",
                       "table also gives score-conditional probabilities with",
                       "fitted-null adjustment across persons. Click a row to draw",
                       "that person's kidmap.")),
      plotCard("kidmap", "Kidmap",
        info = "The person diagnostic map (Wright, Mead & Ludlow 1980): thresholds the person achieved print to the right of the logit axis, thresholds not achieved to the left; the dashed line inside its confidence band is the person location. Achieved thresholds above the band and unachieved thresholds below it are unexpected responses.",
        controls = div(class = "d-flex align-items-center gap-1 me-1",
          span(class = "small text-secondary", "Confidence"),
          div(class = "rasch-inline-select",
              selectInput("kid_level", NULL,
                          c("90%" = "0.9", "95%" = "0.95", "99%" = "0.99"),
                          selected = "0.95", width = "85px"))),
        extra = tagList(
          downloadButton("kidmap_all_pdf", "All (PDF)",
                         class = "btn-outline-secondary btn-xs"),
          downloadButton("kidmap_all_zip", "All (PNG)",
                         class = "btn-outline-secondary btn-xs")))),
    accordion(id = "persons_acc", open = FALSE, class = "mt-3",
      accordion_panel("Person fit", value = "persons_pfit",
        plotCard("pfit", hover = TRUE)),
      accordion_panel("Fit residual distribution", value = "persons_rdist",
        plotCard("rdist_p")),
      accordion_panel("Externally weighted estimates", value = "persons_weights",
        accordion_info(
          "A secondary person measure based on externally imposed relative weights. The fitted calibration and the ordinary estimates used for fit, reliability, targeting and DIF do not change."),
        layout_columns(col_widths = breakpoints(sm = 12, lg = c(4, 8)),
          card(card_body(
            radioButtons("person_weight_level",
              info_label("Weights apply to",
                "Item weights apply to every response cell for that item. Set weights apply one relative weight to each nominated set."),
              c("Items" = "item", "Item sets" = "set"),
              selected = "item", inline = TRUE),
            fileInput("person_weights_file", "Weights CSV", accept = ".csv",
              buttonLabel = "Choose CSV…", placeholder = "No file selected"),
            p(class = "text-muted small",
              "Items: item, weight. Sets: item, set, weight. For an Extended Frames fit, set, weight is sufficient."),
            input_task_button("person_weights_go", "Calculate weighted estimates",
                              type = "primary", class = "w-100"))),
          tableCard("person_weight_tbl", "Weighted person estimates",
            info = "Person locations and sandwich standard errors from the externally weighted score. Zero-weight items are omitted. These are supplementary estimates, not replacements for the fitted Rasch measures."))))),
    # paired-comparison (BTL) fits: the judges are the persons here, so the
    # page carries their fit and their transitivity consistency -- the two
    # judge-level lenses (offered only when a judge column was nominated)
    conditionalPanel("output.is_btl == true && output.has_judges == true",
      accordion(id = "btl_judge_acc", open = FALSE,
        accordion_panel(
          title = "Judge fit",
          value = "btl_judge_fit",
          accordion_info(
            "The fit residual, infit, and outfit are calculated over each judge's comparisons."),
          layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
            tableCard("btl_judges_tbl",
              controls = cols_switch("btl_judges_full"),
              info = paste("Fit statistics pooled over each judge's comparisons.",
                           "A fit bootstrap adds fitted-null adjusted probabilities.",
                           "Click a row to map unexpected judgements.")),
            plotCard("btl_judge_map", title = "Unexpected judgements",
                     info = paste("Matchup residuals for the selected judge.",
                                  "Red matchups favour the weaker object and pass the Holm-adjusted familywise rule.",
                                  "Tied locations have no directional flag."),
                     height = "460px", hover = TRUE))),
        accordion_panel(
          title = "Judge consistency",
          value = "btl_judge_consistency",
          accordion_info(
            "Consistency is one minus the judge's circular-triad rate divided by the chance rate. A value of one is a transitive order; zero is the random-tournament benchmark."),
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(5, 7)),
            tableCard("btl_trans_judges_tbl", title = "Consistency by judge",
                      note = "Judges sorted least consistent first."),
            plotCard("btl_judge_consist", title = "Consistency dotplot",
                     info = "Each judge against the chance line at 0 and the clean-order line at 1. Judges near or below chance contribute little signal.",
                     height = "460px")))))
  )

# ------------------------------------------------------------ TARGETING --
panel_targeting <- nav_panel("Targeting", value = "p_targeting", icon = bs_icon("bullseye"),
    # paired-comparison (BTL) fits: the person-item targeting displays need a
    # person distribution that paired comparisons do not produce, so they hide
    # and the design-information analogues take their place
    conditionalPanel("output.is_btl == true",
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
        tableCard("btl_info_tbl", "Design information"),
        plotCard("btl_targeting_plot", "Design information and targeting")),
      conditionalPanel("output.active_btlef != true && output.btl_history_model != true",
      accordion(class = "mt-3", open = FALSE,
        accordion_panel("Next most informative pairs", value = "btl_next_panel",
          layout_columns(col_widths = breakpoints(sm = 12, md = c(4, 8)),
            div(
              numericInput("btl_next_n", "Pairs to recommend", value = 10,
                           min = 3, max = 50),
              checkboxInput("btl_next_wse",
                info_label("Prioritise poorly measured objects",
                  paste("Ranks candidate pairs by the information in one further",
                        "comparison. Priority weighting favours pairs expected to",
                        "reduce total location uncertainty most (Pollitt 2012).")),
                TRUE)),
            tableCard("btl_next_tbl", "Recommended comparisons",
              info = "Ranks pairs by a variance-reduction approximation, or by information without weighting. The stronger object is presented first; fitted position effects are retained. Reliability from an adaptive sample can be optimistic (Bramley 2015)."))))),
      conditionalPanel("output.active_btlef != true && output.btl_history_model == true",
        div(class = "alert alert-light border mt-3",
            "Pair recommendations require a specified judge and comparison history for this model. The observed-design information remains available.")),
      conditionalPanel("output.active_btlef == true",
        div(class = "alert alert-light border mt-3 d-flex align-items-center gap-2",
            bs_icon("info-circle"),
            span("Pair recommendations need a target panel and object set. Undo the frame adjustment to rank equal-unit pairs.")))),
    conditionalPanel("output.is_btl != true",
      layout_columns(col_widths = 12,
        plotCard("pim_p", "Person-item threshold distribution"),
        plotCard("wright", "Wright map", height = "640px",
          controls_top = TRUE,
          controls = tagList(
            div(class = "rasch-inline-select",
              selectInput("wright_renderer",
                info_label("Renderer",
                  paste("Native uses the package plot. WrightMap uses the",
                        "optional WrightMap package and enables panel layouts.")),
                choices = c("Native" = "native", "WrightMap" = "wrightmap"),
                selected = "native", width = "auto")),
            conditionalPanel("input.wright_renderer === 'wrightmap'",
              div(class = "rasch-inline-select",
                selectInput("wright_type",
                  info_label("Item estimates",
                    paste("Thresholds show every category transition; locations",
                          "show one mean threshold location per item.")),
                  choices = c("Thresholds" = "thresholds",
                              "Locations" = "locations"),
                  selected = "thresholds", width = "auto")),
              div(class = "rasch-inline-select",
                selectInput("wright_person_panels",
                  info_label("Person panels",
                    paste("Optionally separates the person distribution by a",
                          "retained person factor. The default is one panel.")),
                  choices = c("One panel" = ""), selected = "", width = "auto")),
              div(class = "rasch-inline-select",
                selectInput("wright_item_panels",
                  info_label("Item panels",
                    paste("Optionally separates items using an uploaded map,",
                          "or structural response columns by their fitted set,",
                          "group, frame, or facet.")),
                  choices = c("One panel" = ""), selected = "", width = "auto")),
              conditionalPanel("input.wright_item_panels === 'uploaded'",
                div(class = "rasch-inline-select",
                  fileInput("wright_item_map",
                    info_label("Item panel map (CSV)",
                      paste("Two columns, item and panel. Every displayed item",
                            "must occur once.")),
                    accept = ".csv", width = "220px"))),
              div(class = "rasch-inline-select",
                selectInput("wright_person_style",
                  info_label("Person display",
                    "Shows each person panel as a histogram or density curve."),
                  choices = c("Histogram" = "histogram", "Density" = "density"),
                  selected = "histogram", width = "auto")))
          ))),
      div(class = "rasch-plot-toolbar mb-3",
        div(class = "rasch-class-intervals",
          selectInput("tg_bins",
            info_label("Histogram bins",
              paste("Sets the grouping of the person and threshold",
                    "distributions. It does not change the estimates.")),
            choices = seq(10, 60, 5), selected = 35, selectize = FALSE,
            width = "100%")),
        div(class = "rasch-class-intervals",
          selectInput("tg_group",
            info_label("Person group",
              paste("Restricts the person distribution to one level of a",
                    "fitted person factor. The whole sample is shown by",
                    "default; the selection is named on the plot.")),
            choices = c("All persons" = ""), selected = "",
            selectize = FALSE, width = "100%")),
        div(class = "rasch-class-intervals",
          selectInput("tg_items",
            info_label("Item set",
              paste("Restricts the threshold distribution to one item set.",
                    "The whole instrument is shown by default; the",
                    "selection is named on the plot.")),
            choices = c("All items" = ""), selected = "",
            selectize = FALSE, width = "100%")),
        div(class = "rasch-inline-check pb-1",
          checkboxInput("tg_information",
            info_label("Test information",
              paste("Adds conditional test information on a separate",
                    "right-hand scale. Distinct administrable designs are",
                    "shown separately.")),
            value = FALSE, width = "auto")),
        axis_control("tg_axis", standard = c(-5, 5), wide = c(-8, 8),
                     label = "Plot axes")))
  )

# ------------------------------------------------------------------ DIF --
panel_dif <- nav_panel("DIF", value = "p_dif", icon = bs_icon("sliders"),
    # Rasch fits: person-factor DIF (hidden while a BTL fit is active)
    conditionalPanel("output.is_btl != true",
    app_layout_sidebar(
      sidebar = sidebar(width = 280, open = "desktop",
        conditionalPanel("output.dif_multifactor == true",
          radioButtons("dif_effects", info_label("Model",
                       paste("Main effects tests factors jointly. With interactions",
                             "also tests whether one factor's DIF changes over levels",
                             "of another; a supported interaction supersedes its",
                             "constituent main effects.")),
                       c("Main effects" = "main",
                         "With interactions" = "factorial"))),
        numericInput("dif_alpha",
                     info_label("Significance level (alpha)",
                       "Applied to Holm-adjusted probabilities across the DIF family."),
                     value = 0.05,
                     min = 0.001, max = 0.5, step = 0.01),
        conditionalPanel("output.dif_refit_available == true",
          hr(),
          h6(span("Resolve DIF",
                  info_icon(paste("Resolve selected splits a uniform-DIF item by",
                                  "the selected factor term. Automatic resolution",
                                  "uses the selected factor model and proceeds",
                                  "one item at a time while preserving a",
                                  "viable anchor set. Non-uniform DIF requires item",
                                  "review because a location split does not fit a",
                                  "change in discrimination."),
                            "About DIF resolution"))),
          input_task_button("make_split", "Resolve the selected item",
                            type = "primary", class = "w-100"),
          input_task_button("resolve_all", "Resolve all DIF automatically",
                            type = "primary", class = "w-100 mt-2"),
          conditionalPanel("output.has_override_dif",
            actionButton("reset_split", "Undo this change",
                         class = "btn-outline-warning w-100 mt-2"))),
        conditionalPanel("output.dif_refit_available != true",
          uiOutput("dif_refit_note"))),
      accordion(id = "dif_acc", open = "dif_anova",
        accordion_panel("DIF analysis of variance", value = "dif_anova",
          layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
            tableCard("dif_tbl",
              controls = cols_switch("dif_full"),
                      info = "ANOVA of standardised residuals: a significant factor effect indicates uniform DIF, a significant factor-by-class-interval interaction indicates non-uniform DIF (Andrich & Marais 2019). Click a row to see that item's characteristic curves by group (right) and, below, the pairwise comparisons for the selected item and term.",
                      footer = uiOutput("dif_note")),
            plotCard("dif_icc", "Characteristic curves by group",
              info = paste("The item characteristic curve with observed",
                           "class-interval means for each level of the selected",
                           "factor term. For Extended Frames, model curves are",
                           "shown by frame and the observed points are separated",
                           "by the non-frame DIF factor."))),
        # collapsed by default: the DT output suspends while hidden, so the
        # per-item terms table is only computed when the panel is first opened
        accordion_panel("Full ANOVA table", value = "dif_full_panel",
          tableCard("dif_full_tbl",
                    note = "The complete per-item ANOVA: every model term with its df, sums of squares, mean squares, F, and adjusted probability.")),
        accordion_panel("Bootstrap sensitivity", value = "dif_boot_panel",
          conditionalPanel("output.can_dif_boot != true",
            accordion_info(paste(
              "This analysis needs an estimable DIF model and a fitted",
              "model whose null response structure can be reproduced."))),
          conditionalPanel("output.can_dif_boot == true",
            card(card_body(fillable = FALSE,
              uiOutput("dif_boot_explainer"),
              layout_columns(col_widths = c(4, 4, 4),
                numericInput("dif_boot_B", "Replicates", 999,
                             min = 99, max = 9999, step = 100),
                numericInput("dif_boot_seed", "Seed", 2026,
                             min = 0, step = 1),
                div(class = "mt-4",
                  input_task_button("dif_boot_run", "Run sensitivity analysis",
                                    type = "primary", class = "w-100"))),
              uiOutput("dif_boot_job_controls", inline = TRUE),
              uiOutput("dif_boot_state"))),
            conditionalPanel("output.has_dif_boot == true",
              tableCard("dif_boot_tbl", "Bootstrap DIF probabilities",
                controls = cols_switch("dif_boot_full"),
                note = paste("Adjusted probabilities refer to the complete",
                             "item-by-term family. Read them beside, not in",
                             "place of, the primary DIF table."))))),
        accordion_panel("Post-hoc comparisons", value = "dif_pairwise",
          card(card_header_bar(info = app_help("dif_posthoc_tbl")),
               card_body(fillable = FALSE,
                 layout_columns(col_widths = c(4, 4, 4),
                   numericInput("dif_size_flag", "Practical criterion (logits)",
                                0.5, min = 0.1, step = 0.1),
                   numericInput("dif_size_minn", "Min responders", 20,
                                min = 5, step = 5),
                   div(class = "mt-4 d-flex justify-content-end align-items-start",
                       downloadButton("dl_dif_posthoc", "CSV",
                                      class = "btn-outline-secondary btn-xs"))),
                 uiOutput("dif_posthoc_heading"),
                 uiOutput("dif_posthoc_note"),
                 DT::DTOutput("dif_posthoc_tbl"),
                 rcode_details("dif_posthoc_tbl"),
                 conditionalPanel("output.dif_selected_interaction == true",
                   hr(),
                   h6(span("Resolved cell comparisons",
                           info_icon(app_help("dif_size_tbl"),
                                     "About this table"))),
                 uiOutput("dif_levels_note"),
                   DT::DTOutput("dif_size_tbl"),
                   rcode_details("dif_size_tbl"))))),
        # shown only after an automatic run: the trace of the splits that
        # resolved the DIF (the resolved fit is the active override)
        accordion_panel("Automatic resolution", value = "dif_resolve",
          conditionalPanel("output.has_resolve == true",
            card(
              card_header_bar(
                info = app_help("resolve_tbl"),
                buttons = downloadButton("resolve_tbl_csv", "CSV",
                                         class = "btn-outline-secondary btn-xs")),
              card_body(fillable = FALSE,
                uiOutput("resolve_summary"),
                DT::DTOutput("resolve_tbl"),
                rcode_details("resolve_tbl"))))),
        accordion_panel(
          title = "Planned contrasts",
          value = "dif_contrasts",
          card(card_body(fillable = FALSE,
                 accordion_info(app_help("contr_tbl"), "About this table"),
                 conditionalPanel("output.dif_followup_available != true",
                   p(class = "text-muted small mb-0",
                     "Resolved logit contrasts are unavailable for Extended Frames because an ordinary item split would discard the fitted frame units.")),
                 conditionalPanel("output.dif_followup_available == true",
                   layout_columns(col_widths = c(4, 4, 4),
                     selectizeInput("pc_items", "Items", NULL, multiple = TRUE,
                                    options = list(placeholder = "all items")),
                     div(class = "form-text mt-4",
                         info_label("Repeated measures",
                           "The person identifier stored with the fitted model is used automatically.")),
                     div(class = "mt-4",
                         input_task_button("pc_run", "Derive and test contrasts",
                                           type = "primary"))),
                   conditionalPanel("output.has_contr != true",
                     p(class = "text-muted small mb-0",
                       "Select Derive and test contrasts to see the contrast family and its tests.")),
                   conditionalPanel("output.has_contr == true",
                     uiOutput("contr_family"),
                     div(class = "d-flex justify-content-end align-items-center gap-3",
                         cols_switch("contr_full"),
                         downloadButton("contr_tbl_csv", "CSV",
                                        class = "btn-outline-secondary btn-xs")),
                     DT::DTOutput("contr_tbl"),
                     rcode_details("contr_tbl")))))))
    ))),
    # paired-comparison (BTL) fits: differential object functioning by
    # judge group (Bradley-Terry counterpart of the person-factor analysis)
    conditionalPanel("output.is_btl == true",
    app_layout_sidebar(
      sidebar = sidebar(width = 280, open = "desktop",
        selectizeInput("bdif_factors", "Judge factors", NULL, multiple = TRUE,
                       options = list(placeholder = "nominate judge factors on the Data page")),
        conditionalPanel("output.bdif_multifactor == true",
          radioButtons("bdif_effects", info_label("Model",
                       paste("Main effects tests judge factors jointly. With",
                             "interactions also tests whether object DIF changes",
                             "across combinations of judge factors.")),
                       c("Main effects" = "main",
                         "With interactions" = "factorial"))),
        numericInput("bdif_alpha", "Significance level (alpha)", value = 0.05,
                     min = 0.001, max = 0.5, step = 0.01),
        input_task_button("bdif_run", "Run DIF analysis",
                          type = "primary", class = "w-100")),
      accordion(id = "bdif_acc", open = "bdif_anova",
        accordion_panel("DIF analysis of variance", value = "bdif_anova",
          layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
            tableCard("bdif_anova_tbl",
              info = "ANOVA of each object's standardised residuals. A judge-group effect indicates uniform DIF; a group-by-opponent-band interaction indicates non-uniform DIF. Probabilities are Holm-adjusted. Click a row to show the characteristic curves by judge group.",
              footer = uiOutput("bdif_notes")),
            plotCard("bdif_occ", "Characteristic curves by group",
              info = "The object characteristic curve with the observed mean response per opponent overlaid separately for each judge group: the graphical display of DIF for the object of the selected table row.",
              hover = TRUE))),
        accordion_panel("DIF magnitude in logits", value = "bdif_size_panel",
          tableCard("bdif_sizes_tbl",
            note = "Pairwise differences between resolved group locations, in logits, with Holm-adjusted probabilities. Differences of at least 0.5 logits are flagged as practically significant.")),
        accordion_panel("Bootstrap sensitivity", value = "bdif_boot_panel",
          card(card_body(fillable = FALSE,
            accordion_info(paste(
              "Draws outcomes from the fitted comparison model while",
              "retaining judges and the comparison design, then refits the",
              "model and repeats the complete DIF analysis. The adjusted",
              "residual ANOVA remains the primary result. Bootstrap familywise",
              "probabilities refer to the fitted global invariant null; they",
              "do not provide strong control when another object has DIF.")),
            layout_columns(col_widths = c(4, 4, 4),
              numericInput("bdif_boot_B", "Replicates", 999,
                           min = 99, max = 9999, step = 100),
              numericInput("bdif_boot_seed", "Seed", 2026,
                           min = 0, step = 1),
              div(class = "mt-4",
                input_task_button("bdif_boot_run", "Run sensitivity analysis",
                                  type = "primary", class = "w-100"))),
            uiOutput("bdif_boot_job_controls", inline = TRUE),
            uiOutput("bdif_boot_state"))),
          conditionalPanel("output.has_bdif_boot == true",
            tableCard("bdif_boot_tbl", "Bootstrap DIF probabilities",
              controls = cols_switch("bdif_boot_full"),
              note = paste("Adjusted probabilities refer to the complete",
                           "object-by-term family. Read them beside, not in",
                           "place of, the primary DIF table."))))
    ))
  ))

# --------------------------------------------------------------- FACETS --
panel_facets <- nav_panel("Facets", value = "p_facets", icon = bs_icon("person-badge"),
    app_layout_sidebar(
      sidebar = sidebar(width = 280, open = "desktop",
        selectizeInput("facet_sel", info_label("Facet",
                       paste("Shows jointly calibrated severities for the selected",
                             "facet. Positive values denote greater severity;",
                             "large fit residuals indicate inconsistency.")), NULL,
                       options = list(placeholder = "run a Multiple Ratings analysis"))),
      tableCard("facet_tbl", "Facet severities and fit",
        controls = cols_switch("facets_full"),
        footer = uiOutput("facet_structure_note")),
      conditionalPanel("output.has_interaction == true",
        tableCard("facet_int_omnibus", "Item-by-facet omnibus test",
                  "Joint Wald test of the complete interaction family. Read this before the Holm-adjusted exploratory cell contrasts."),
        tableCard("facet_int_tbl", "Item-by-facet interactions",
                  "Estimated because the interactive facet structure was chosen in the data roles; gamma is the extra severity of a level on a particular item. Cell p-values are Holm-adjusted and are exploratory follow-ups to the omnibus test.")),
      plotCard("facet_plot", "Severity caterpillar plot")
    )
  )

# -------------------------------------------------------------- EQUATING --
panel_equating <- nav_panel("Equating", value = "p_equating", icon = bs_icon("arrow-left-right"),
    # paired-comparison (BTL) fits: common-object equating against a banked
    # calibration (standards maintenance across panels or years); the Rasch
    # common-item machinery hides while a BTL fit is active
    conditionalPanel("output.is_btl == true",
      app_layout_sidebar(
        sidebar = sidebar(width = 300, open = "desktop",
          fileInput("bt_eq_file", info_label("Reference calibration",
                    paste("Common objects link by name. A CSV with marginal",
                          "standard errors is descriptive unless the bank is",
                          "fixed or carries its joint covariance.")),
                    accept = ".csv"),
          checkboxInput("bt_eq_independent",
                        "Independent judges/comparisons", TRUE),
          uiOutput("btl_eq_summary")),
        layout_columns(col_widths = 12,
          tableCard("btl_eq_tbl", "Common-object comparison",
            info = "Each common object is compared with the shifted identity line. The shift is precision-weighted when at least two common objects have usable standard errors; otherwise it is an unweighted descriptive mean. Inferential columns are withheld unless the reference carries its joint covariance or is fixed. Where available, p_adj is the multiplicity-adjusted drift p-value, shown red below 0.05."),
          plotCard("btl_eq_plot", "Equating plot",
            info = "The two calibrations' common-object locations against the shifted identity line, with per-object error bars and a dotted guide band at the average pooled precision; objects that drift after the multiplicity adjustment are highlighted and labelled.",
            hover = TRUE))
      )
    ),
    conditionalPanel("output.is_btl != true",
    app_layout_sidebar(
      sidebar = sidebar(width = 300, open = "desktop",
        radioButtons("eq_source", info_label("Reference",
                     paste("Common items link by name. Drift inference requires",
                           "independent sampling units and joint location covariance;",
                           "a CSV with marginal standard errors is descriptive unless fixed.")),
                     c("Uploaded calibration CSV" = "csv",
                       "A kept fit from Compare" = "kept")),
        conditionalPanel("input.eq_source == 'csv'",
          fileInput("eq_file", "Reference calibration (CSV: item,location; optional se,max)",
                    accept = ".csv"),
          checkboxInput("eq_csv_independent",
                        "Independent sampling units", TRUE)),
        conditionalPanel("input.eq_source == 'kept'",
          selectizeInput("eq_kept", "Kept fit", NULL,
                         options = list(placeholder = "keep a fit on the Compare page")),
          checkboxInput("eq_kept_independent",
                        "Independent sampling units", FALSE)),
        radioButtons("eq_shift", "Scale alignment",
                     c("Allow a shift between origins" = "mean",
                       "Compare raw locations (anchored scales)" = "none")),
        downloadButton("dl_calib", "Save current calibration (CSV)",
                       class = "btn-outline-secondary w-100")),
      card(
        id = "card_eq_tbl",
        full_screen = TRUE,
        card_header_bar("Common-item comparison",
          info = app_help("eq_tbl"),
          buttons = conditionalPanel("output.has_eq == true",
            div(class = "rasch-chips",
                cols_switch("eq_full"),
                downloadButton("eq_tbl_csv", "CSV",
                               class = "btn-outline-secondary btn-xs")))),
        card_body(
          conditionalPanel("output.has_eq != true",
            p(class = "text-muted small mb-0",
              "Upload a reference calibration (or choose a kept fit) to see the common-item comparison.")),
          conditionalPanel("output.has_eq == true",
            DTOutput("eq_tbl"), rcode_details("eq_tbl")),
          padding = 12, fillable = FALSE)),
      card(
        id = "card_eq_plot",
        full_screen = TRUE,
        `data-bs-theme` = "light",
        card_header_bar("Equating plot",
          info = app_help("eq_plot"),
          buttons = conditionalPanel("output.has_eq == true",
            div(class = "rasch-chips",
                downloadButton("eq_plot_png", "PNG", class = "btn-outline-secondary btn-xs"),
                downloadButton("eq_plot_pdf", "PDF", class = "btn-outline-secondary btn-xs")))),
        card_body(
          conditionalPanel("output.has_eq != true",
            p(class = "text-muted small mb-0",
              "Upload a reference calibration (or choose a kept fit) to see the equating plot.")),
          conditionalPanel("output.has_eq == true",
            # hand-built card (not plotCard()): the same hover = TRUE wiring
            # plotCard() adds -- a hoverOpts plotOutput in a position:relative
            # div, with the tip floated over it by output$eq_plot_tip.
            div(style = "position: relative",
                plotOutput("eq_plot", height = "560px",
                           hover = hoverOpts("eq_plot_hover", delay = 60,
                                             delayType = "throttle")),
                uiOutput("eq_plot_tip")),
            rcode_details("eq_plot")),
          padding = 8, fillable = FALSE))
    ))
  )

# --------------------------------------------------------------- FRAMES --
panel_frames <- nav_panel("Extended Frames", value = "p_frames", icon = bs_icon("grid-3x3"),
    # paired-comparison (BTL) fits: the paired-comparison extended frame of
    # reference (judge-panel by object-set cells, this package's extension of
    # Humphry's model to the Bradley-Terry-Luce family); the Rasch EFRM
    # machinery below hides while a BTL fit is active
    conditionalPanel("output.is_btl == true",
      app_layout_sidebar(
        sidebar = sidebar(width = 340, open = "desktop",
          selectInput("btlef_panel", "Judge-panel column", NULL),
          fileInput("btlef_sets_file", info_label("Object sets (CSV: object, set)",
                    "If omitted, sets are inferred from the part of each object name before its trailing digits."),
                    accept = ".csv", placeholder = "optional"),
          numericInput("btlef_boot", info_label("Bootstrap replicates",
                       paste("More replicates give more stable standard errors",
                             "and probabilities but take longer.")), value = 200,
                       min = 50, max = 999),
          numericInput("btlef_seed", info_label("Bootstrap seed",
                       "Reproduces the same judge resamples for every selected worker count."),
                       value = 1, min = 0, step = 1),
          selectInput("btlef_workers", info_label("Parallel workers",
                      paste("Runs judge-bootstrap refits concurrently.",
                            "Defaults to four, or fewer when the system limit is lower.")),
                      choices = .efrm_worker_choices,
                      selected = max(.efrm_worker_values)),
          radioButtons("btlef_se", info_label("Standard errors",
                       paste("The bootstrap carries uncertainty through the frame",
                             "linking stages. Conditional standard errors are faster",
                             "but omit calibration uncertainty.")),
                       c("Judge bootstrap (recommended)" = "judge_bootstrap",
                         "Parametric bootstrap" = "bootstrap",
                         "Conditional (fast, understates)" = "conditional")),
          input_task_button("btlef_run", "Estimate frame units",
                            type = "primary", class = "w-100"),
          uiOutput("btlef_job_controls")),
        conditionalPanel("output.has_btlef != true",
          card(card_body(p(class = "text-muted small mb-0",
            "Choose the judge-panel column in the sidebar and press Estimate frame units to see results.")))),
        conditionalPanel("output.has_btlef == true",
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)),
            tableCard("btlef_phi_tbl", "Panel units (phi)",
              info = "The panel's discriminating power on the common scale; Wald test against log phi = 0."),
            conditionalPanel("output.btlef_multiset == true",
              tableCard("btlef_units_tbl", "Set units (alpha) and origins (kappa)",
                info = "Identified from the cross-set comparisons alone (no distributional assumption about the objects); the reference set is fixed at alpha = 1, kappa = 0."))),
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(7, 5)),
            plotCard("btlef_units_plot", "Frame units", height = "460px",
              info = "Caterpillar plot of the log units with reference-appropriate 95% intervals; the reference (one) is marked."),
            tableCard("btlef_frames_tbl", "Frame fit",
              info = "Each judge-panel by object-set cell holding within-set comparisons: rho = phi / alpha is the cell's common-scale discrimination, n_comparisons the comparison count, and fit_resid the pooled fit residual.")),
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)),
            tableCard("btlef_cmp_tbl", "Frame model comparison",
              info = app_help("btlef_cmp_tbl")),
            tableCard("btlef_omnibus_tbl", "Omnibus tests of equal units",
              info = app_help("btlef_omnibus_tbl"))),
          uiOutput("btlef_note"))
      )
    ),
    conditionalPanel("output.is_btl != true",
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(7, 5)),
        tableCard("frame_tbl", "Extended frames: units, origins, pooled fit",
                  controls = cols_switch("frames_full")),
        div(tableCard("phi_tbl", "Person group units (phi)"),
            tableCard("alpha_tbl", "Item set units (alpha) and locations"))),
      layout_columns(col_widths = 12,
        div(class = "rasch-section-toolbar",
            div(class = "rasch-inline-select",
              selectInput("inv_se",
                info_label("Uncertainty",
                  paste("Conditional provides fast location tests and descriptive",
                        "discrimination comparisons. Bootstrap resamples whole",
                        "person rows within person group, preserving their item-set",
                        "patterns, and provides inference for both.")),
                c("Conditional" = "conditional", "Bootstrap" = "bootstrap"),
                selected = "conditional", width = "150px")),
            conditionalPanel("input.inv_se == 'bootstrap'",
              numericInput("inv_boot",
                info_label("Replicates",
                  paste("Number of complete whole-person refits.",
                        "Rows are resampled within person group and keep",
                        "their observed item-set pattern.")),
                200, min = 30, max = 1000, step = 50, width = "110px"),
              numericInput("inv_seed",
                info_label("Seed",
                  "Fixes the bootstrap resamples so the result can be reproduced."),
                1, min = 1, step = 1, width = "100px"))),
        tableCard("frame_inv_summary_tbl", "Item invariance summary",
                  info = app_help("frame_inv_summary_tbl"))),
      layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
        tableCard("frame_inv_loc_tbl", "Item invariance: locations",
                  info = app_help("frame_inv_loc_tbl")),
        tableCard("frame_inv_disc_tbl", "Item invariance: discrimination",
                  info = app_help("frame_inv_disc_tbl"))),
      layout_columns(col_widths = 12,
        plotCard("frame_plot", "Frame units"),
        plotCard("frame_icc", "ICC across frames",
          controls = div(class = "rasch-icc-compare",
            selectizeInput("frame_item",
              info_label("Item",
                paste("Selects the item whose characteristic curve is drawn",
                      "across the fitted frames.")),
              NULL, options = list(placeholder = "run an Extended Frames analysis"),
              width = "100%")))),
      layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)),
        tableCard("efrm_cmp_tbl", "Frame model comparison",
          info = app_help("efrm_cmp_tbl")),
        tableCard("efrm_omnibus_tbl", "Omnibus tests of equal units",
          info = app_help("efrm_omnibus_tbl")))
    )
  )

# -------------------------------------------------- INDEPENDENCE: TRAIT --
panel_dim <- nav_panel("Trait", value = "p_dim", icon = bs_icon("diagram-3"),
    # paired-comparison (BTL) fits: the Rasch residual-PCA suite needs a
    # persons x items residual matrix that paired comparisons do not produce,
    # so it hides and the pair-structure analogues take its place
    conditionalPanel("output.is_btl == true",
      accordion(id = "btl_dim_acc", open = "btl_dim_swirl",
        accordion_panel(
          title = "Residual dimensions",
          value = "btl_dim_swirl",
          accordion_info(
            "Observed and fitted expected points are compared within each object pair, retaining fitted thresholds and position or history effects. Their residual log-odds matrix is decomposed into rotational planes, or bimensions (Gower 1977)."),
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)),
            plotCard("btl_scree", title = "Bimension strengths",
                     info = "Each bimension's strength and, when available, the simulated mean and finite-simulation 5% upper reference. A leading bar above the reference suggests residual structure under the conditional independence assumption.",
                     height = "460px"),
            plotCard("btl_dim_map", title = "Leading residual map",
                     info = "Objects in the leading residual plane. Coordinates show the pattern, not its magnitude or significance; use the scree plot for the conditional reference. Point size grows with location on the main scale.",
                     height = "460px")),
          tableCard("btl_bimensions_tbl", title = "Bimensions",
                    note = "Strength, share of the total residual, and, when available, the noise reference for the leading bimension.")),
        accordion_panel(
          title = "Preference loops",
          value = "btl_dim_loops",
          accordion_info(
            "A circular triad occurs when A is preferred to B, B to C, and C to A. One quarter is the random-tournament benchmark (Kendall and Babington Smith 1940)."),
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(7, 5)),
            tableCard("btl_trans_tbl", title = "Transitivity summary",
                      note = "Circular triads out of the complete triples, the loop rate against the 25% chance rate, and Kendall's coefficient of consistency when every pair was compared."),
            plotCard("btl_involve_plot", title = "Objects in loops",
                     info = "How many circular triads each object sits in - the objects whose order is least stable, and the likeliest seat of a second attribute.",
                     height = "460px"))))),
    conditionalPanel("output.is_btl != true",
    accordion(id = "dim_acc", open = "dim_components",
      accordion_panel("Residual components", value = "dim_components",
        layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)),
          tableCard("loadings_tbl", title = "Loadings",
                    note = "First 10 components shown."),
          plotCard("pca_biplot", title = "Biplot (PC1 vs PC2)",
                   height = "auto"))),
      accordion_panel("Scree plot", value = "dim_scree",
        plotCard("scree")),
      accordion_panel(
        title = "Unidimensionality t-test",
        value = "dim_ttest",
        accordion_info(
          paste("Smith's test compares each person's estimates from two item",
                "subsets. A split fixed in advance uses the binomial interval;",
                "a residual-derived split needs bootstrap calibration for an",
                "inferential verdict.")),
        layout_columns(col_widths = breakpoints(sm = 12, xl = c(4, 8)),
          div(
            h6(span("t-test item subsets",
                    info_icon(paste("Leave both subsets empty to use opposing",
                                    "loadings on the selected component. Persons",
                                    "extreme on either subset are excluded."),
                              "About item subsets"))),
            div(class = "mb-2 d-flex align-items-center gap-2",
              span(class = "small text-secondary", "Automatic split component"),
              div(class = "rasch-inline-select",
                  selectInput("pca_component", NULL, choices = 1, selected = 1,
                              width = "80px"))),
            selectizeInput("dim_pos", "Subset A", NULL, multiple = TRUE,
                           options = list(placeholder = "positive loadings on the selected component")),
            selectizeInput("dim_neg", "Subset B", NULL, multiple = TRUE,
                           options = list(placeholder = "negative loadings on the selected component")),
            numericInput("dim_boot_B", info_label(
              "Bootstrap replicates",
              paste("Leave at zero for a descriptive automatic-split result.",
                    "A positive value calibrates the data-driven split;",
                    "99 or more is suitable for inference.")),
              value = 0, min = 0, step = 99),
            div(class = "d-flex gap-2",
              div(class = "flex-fill",
                selectInput("dim_workers", info_label(
                  "Parallel workers",
                  "Bootstrap results reproduce across worker counts."),
                  choices = .efrm_worker_choices,
                  selected = max(.efrm_worker_values))),
              div(class = "flex-fill",
                numericInput("dim_boot_seed", info_label(
                  "Random seed", "Used only for bootstrap calibration."),
                  value = 1, min = 0, step = 1))),
            input_task_button("dim_apply", "Run t-test",
                              type = "primary", class = "w-100")),
          card(card_body(verbatimTextOutput("dim_txt"), rcode_details("dim"))))),
      accordion_panel("Magnitude of multidimensionality", value = "dim_magnitude",
        card(
          id = "card_dm_tbl",
          full_screen = TRUE,
          card_header_bar(info = app_help("dm_tbl")),
          card_body(
            div(input_task_button("dm_run", "Estimate from current subsets",
                                  type = "primary")),
            conditionalPanel("output.has_dm != true",
              p(class = "text-muted small mb-0 mt-2",
                "Select Estimate from current subsets to see the reliability comparison.")),
            conditionalPanel("output.has_dm == true",
              div(class = "d-flex justify-content-end",
                  downloadButton("dm_tbl_csv", "CSV",
                                 class = "btn-outline-secondary btn-xs")),
              DTOutput("dm_tbl"), rcode_details("dm_tbl")),
            padding = 12, fillable = FALSE))),
      accordion_panel("Eigenvalues", value = "dim_eigen",
        tableCard("eigen_tbl", note = "First 10 eigenvalues shown."))))
  )

# -------------------------------------------------- INDEPENDENCE: LOCAL --
panel_ld <- nav_panel("Local", value = "p_ld", icon = bs_icon("link-45deg"),
    # paired-comparison (BTL) fits: within-judge dependence estimated from
    # the judgment order; the Rasch Q3 suite hides while a BTL fit is active
    conditionalPanel("output.is_btl == true",
      conditionalPanel("output.has_btl_dep != true",
        card(card_body(p(class = "text-muted small mb-0",
          "Nominate a judgment-order column in the Data roles to estimate within-judge dependence.")))),
      conditionalPanel("output.has_btl_dep == true",
        layout_columns(col_widths = breakpoints(sm = 12, lg = c(5, 7)),
          tableCard("btl_dep_tbl", "Order effects",
            info = "Exposure is the effect of having seen an object before. Carry-over is the effect of earlier verdicts involving that object. If fitted, first-position advantage is a static presentation effect rather than within-judge dependence. All effects are estimated jointly with the object locations.",
            note = "An effect estimated on few informative comparisons carries a wide standard error; the plot and the comparison table below show which comparisons drive it."),
          plotCard("btl_dep_plot", "Dependence effect",
            info = "Observed residual departures are binned by the history covariate and shown with the fitted effect. Bin counts show how much information supports the estimate.",
            controls = div(class = "d-flex align-items-center gap-1 me-1",
              span(class = "small text-secondary", "Effect"),
              div(class = "rasch-inline-select",
                  selectInput("btl_dep_effect", NULL,
                              c("Exposure" = "exposure",
                                "Carry-over" = "carry_over"),
                              width = "130px"))))),
        accordion(class = "mt-3", open = FALSE,
          accordion_panel("Comparison covariates", value = "btl_dep_comps_panel",
            tableCard("btl_dep_comps",
              note = "Every comparison in judgment order with its exposure and carry-over covariates; a comparison is informative for an effect when its covariate is non-zero."))))),
    conditionalPanel("output.is_btl != true",
    accordion(id = "ld_acc", open = "ld_cormat",
      accordion_panel("Residual Correlations (Q3 statistics)", value = "ld_cormat",
        numericInput("ld_flag",
                     info_label("Flag threshold",
                       paste("Flags |Q3| or Q3* at or above this value. No threshold",
                             "is a universal critical value; 0.2 is a common",
                             "screening value rather than a test.")),
                     value = 0.2, min = 0.05, max = 0.9, step = 0.05,
                     width = "420px"),
        layout_columns(col_widths = 12,
          tableCard("cormat_q3_tbl", title = "Q3 correlations",
                    info = "The residual correlation of every item or response-cell pair, with 1.00 on the diagonal (Yen 1984). Pairs with |Q3| at or above the threshold are red (Yen 1993)."),
          plotCard("rcor_q3", height = "auto", hover = TRUE)),
        layout_columns(col_widths = 12,
          tableCard("cormat_q3s_tbl", title = "Adjusted Q3 (Q3*)",
                    info = "Each item or response-cell Q3 less the average off-diagonal Q3: 0 marks local independence. Pairs with Q3* at or above the threshold are red (Christensen, Makransky & Horton 2017)."),
          plotCard("rcor_q3s", height = "auto", hover = TRUE))),
      accordion_panel("Response dependence magnitude", value = "ld_dep",
        card(
          id = "card_dep_tbl",
          full_screen = TRUE,
          card_header_bar(info = app_help("dep_tbl")),
          card_body(
            conditionalPanel("output.dep_magnitude_available == true",
              div(class = "d-flex gap-3 flex-wrap align-items-end",
                selectizeInput("dep_item", "Dependent item", NULL, width = "190px",
                               options = list(placeholder = "run an analysis first")),
                selectizeInput("ind_item", "Independent item", NULL, width = "190px",
                               options = list(placeholder = "run an analysis first")),
                div(class = "mb-3",
                    input_task_button("run_dep", "Estimate d", type = "primary")))),
            conditionalPanel("output.dep_magnitude_available != true",
              p(class = "text-muted small mb-0",
                "Response-dependence magnitude is unavailable for mutually exclusive Extended Frames.")),
            conditionalPanel("output.has_dep != true",
              p(class = "text-muted small mb-0",
                "Select Estimate d to see the resolved magnitudes.")),
            conditionalPanel("output.has_dep == true",
              div(class = "d-flex justify-content-end",
                  downloadButton("dep_tbl_csv", "CSV",
                                 class = "btn-outline-secondary btn-xs")),
              verbatimTextOutput("dep_txt"),
              DTOutput("dep_tbl"), rcode_details("dep_tbl")),
            padding = 12, fillable = FALSE))),
      accordion_panel("Subtest (combine dependent items)", value = "ld_subtest",
        card(
          card_body(
            conditionalPanel("output.subtest_available == true",
              selectizeInput("subtest_items", info_label("Items to combine",
                             paste("Merges two or more items into a polytomous",
                                   "super-item and refits, absorbing their local",
                                   "dependence into the subtest.")), NULL, multiple = TRUE,
                             options = list(placeholder = "items to combine")),
              div(input_task_button("make_subtest", "Combine and re-analyse",
                                    type = "primary")),
              conditionalPanel("output.has_override_subtest",
                div(class = "mt-2",
                  actionButton("reset_subtest", "Undo this change",
                               class = "btn-outline-warning w-100")))),
            conditionalPanel("output.subtest_available != true",
              p(class = "text-muted small mb-0",
                "Form subtests in the source data before fitting Multiple Ratings or Extended Frames.")),
            uiOutput("subtest_status")))),
      accordion_panel("Spread test (LUB)", value = "ld_spread",
        card(
          id = "card_spread_tbl",
          full_screen = TRUE,
          card_header_bar(info = app_help("spread_tbl")),
          card_body(
            conditionalPanel("output.spread_available == true",
              div(class = "d-flex flex-wrap gap-2 align-items-end mb-2",
                  numericInput("spread_alpha",
                    info_label("Significance level",
                      "Applied after adjustment to the one-sided probabilities."),
                    value = 0.05, min = 0.001, max = 0.25, step = 0.01,
                    width = "160px"),
                  div(class = "mb-3",
                    input_task_button("run_spread", "Run spread test",
                                      type = "primary")))),
            conditionalPanel("output.has_spread != true",
              p(class = "text-muted small mb-0",
                "Run after forming a superitem to compare its spread with the binomial bound.")),
            conditionalPanel("output.has_spread == true",
              div(class = "d-flex justify-content-end",
                  downloadButton("spread_tbl_csv", "CSV",
                                 class = "btn-outline-secondary btn-xs")),
              DTOutput("spread_tbl"),
              rcode_details("spread_tbl")),
            padding = 12, fillable = FALSE)))))
  )

# ------------------------------------------------------------- GUESSING --
panel_guess <- nav_panel("Guessing", value = "p_guess", icon = bs_icon("question-diamond"),
    app_layout_sidebar(
      sidebar = sidebar(width = 300, open = "desktop",
        numericInput("guess_chance", info_label("Chance success probability",
                     paste("Responses with modelled success probability below",
                           "this value are removed before recalibration. Bootstrap",
                           "inference repeats the complete procedure.")),
                     value = 0.25, min = 0.05, max = 0.95, step = 0.05),
        selectizeInput("guess_anchors", "Anchor items (common origin)", NULL,
                       multiple = TRUE,
                       options = list(placeholder = "automatic: least-affected third")),
        checkboxInput("guess_bootstrap", "Person-bootstrap inference", FALSE),
        conditionalPanel("input.guess_bootstrap == true",
          numericInput("guess_boot_reps", info_label("Bootstrap replicates",
                       "Holm-adjusted tests need at least 40 usable draws per item. Allow more for failed refits; fewer draws give descriptive uncertainty only."), 999,
                       min = 50, step = 50),
          numericInput("guess_seed", "Random seed", 1, min = 0, step = 1)),
        input_task_button("run_guess", "Run tailored analysis",
                          type = "primary", class = "w-100")),
      layout_columns(col_widths = 12,
        card(card_header("Tailored analysis"),
             card_body(verbatimTextOutput("guess_txt"))),
        tableCard("guess_tbl", "Initial vs tailored calibration",
                  "shift = tailored minus origin-equated location. Shifts are descriptive by default; bootstrap percentile intervals and Holm-adjusted p-values appear only when person-bootstrap inference is selected.")),
      plotCard("guess_plot", "Tailored vs origin-equated calibrations", hover = TRUE)
    )
  )

# -------------------------------------------------------------- COMPARE --
panel_compare <- nav_panel("Compare", value = "p_compare", icon = bs_icon("columns-gap"),
    app_layout_sidebar(
      sidebar = sidebar(width = 300, open = "desktop",
        input_task_button("keep_fit", "Keep current fit for comparison",
                          type = "primary", class = "w-100"),
        actionButton("clear_fits", "Clear kept fits",
                     class = "btn-outline-secondary w-100 mt-2"),
        selectizeInput("cmp_ref", info_label("Reference fit",
                       paste("Information criteria are comparable only when fits",
                             "use the same data. For different preparations,",
                             "compare fit, targeting and reliability descriptively.")), NULL,
                       options = list(placeholder = "keep at least two fits"))),
      card(
        id = "card_cmp_tbl",
        full_screen = TRUE,
        card_header_bar("Model comparison",
          info = app_help("cmp_tbl"),
          buttons = conditionalPanel("output.has_cmp == true",
            div(class = "rasch-chips",
                cols_switch("cmp_full"),
                downloadButton("cmp_tbl_csv", "CSV",
                               class = "btn-outline-secondary btn-xs")))),
        card_body(
          conditionalPanel("output.has_cmp != true",
            p(class = "text-muted small mb-0",
              "Keep at least two fits (run, keep, change the settings, run and keep again) to see the comparison.")),
          conditionalPanel("output.has_cmp == true",
            p(class = "text-muted small mb-2",
              "Reference for the log-likelihood comparison is the fit chosen in the sidebar."),
            uiOutput("cmp_family_note"),
            DTOutput("cmp_tbl"), rcode_details("cmp_tbl")),
          padding = 12, fillable = FALSE))
    )
  )

# --------------------------------------------------------------- EXPORT --
panel_export <- nav_panel("Export", value = "p_export", icon = bs_icon("download"),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
      card(info_header("Save or reopen the analysis",
        "The .rasch file stores the source data, fitted model, downstream transformations, comparison fits and reproducible R call. Reopen it in this app to continue the analysis."),
        card_body(
          downloadButton("dl_project", "Save analysis (.rasch)",
                         class = "btn-primary", icon = bs_icon("floppy")),
          fileInput("project_file", NULL, accept = ".rasch",
                    buttonLabel = "Open analysis…",
                    placeholder = "Choose a saved .rasch file"))),
      card(info_header("Analysis report",
        "A formatted report of the active fit. HTML is self-contained; Word is editable; PDF requires a LaTeX installation such as TinyTeX."),
        card_body(
          radioButtons("report_format", NULL,
                       c("HTML" = "html", "Word" = "docx", "PDF" = "pdf"),
                       selected = "html", inline = TRUE),
          downloadButton("dl_report", "Download report",
                         class = "btn-primary", icon = bs_icon("file-earmark-text")))),
      card(info_header("Model results archive",
        paste("A ZIP containing the model tables supported by save_outputs(),",
              "plots in the selected formats, and a plain-text summary.",
              "On-demand app analyses without an export format remain in the",
              "saved .rasch analysis and in their table-specific downloads.")),
        card_body(
          checkboxGroupInput("exp_formats", "Plot formats",
                             c("PNG" = "png", "PDF" = "pdf"), selected = c("png", "pdf"),
                             inline = TRUE),
          checkboxInput("exp_items", "Include all item or object plots", TRUE),
          downloadButton("dl_zip", "Download model results (ZIP)",
                         class = "btn-primary")))
    )
  )

# --------------------------------------------------------------- SIMULATE --
# Generate data from any model family with dial-in departures, load it as the
# current dataset, and run the analysis to watch the matching diagnostic fire.
panel_simulate <- nav_panel("Simulate", value = "p_simulate", icon = bs_icon("dice-5"),
  app_layout_sidebar(
    sidebar = sidebar(width = 400, open = "desktop",
      radioButtons("sim_layout", "Data type", c(
        "Rasch" = "rasch",
        "Explanatory Rasch" = "rasch_exp",
        "Comparative Judgement" = "btl",
        "Explanatory Comparative Judgement" = "btl_exp",
        "Multiple Ratings (MFRM)" = "mfrm",
        "Extended Frames (EFRM)" = "efrm",
        "Extended Frames (Comparative Judgement)" = "btl_efrm")),
      # ---------------------------------------------------------- Rasch ----
      conditionalPanel("input.sim_layout == 'rasch' || input.sim_layout == 'rasch_exp'",
        sliderInput("sr_persons", "Persons", 100, 2000, 600, 50),
        sliderInput("sr_items", "Items", 5, 40, 15, 1),
        radioButtons("sr_model", "Model", inline = TRUE,
          c("Dichotomous" = "dichotomous", "Partial credit" = "PCM",
            "Rating scale" = "RSM")),
        conditionalPanel("input.sr_model != 'dichotomous'",
          sliderInput("sr_cats", "Categories", 3, 6, 4, 1)),
        accordion(open = FALSE, accordion_panel(
          title = "Population", value = "sr_pop",
          sliderInput("sr_mean", "Person mean (logits)", -2, 2, 0, 0.25),
          sliderInput("sr_sd", "Person SD (logits)", 0.4, 2.5, 1, 0.1),
          selectInput("sr_dist", "Person distribution", c(
            "Normal" = "normal", "Uniform" = "uniform",
            "Skewed" = "skew", "Bimodal" = "bimodal")),
          sliderInput("sr_diff", "Item difficulty range", -4, 4,
                      c(-2.5, 2.5), 0.25))),
        conditionalPanel("input.sim_layout == 'rasch_exp'",
          accordion(open = FALSE, accordion_panel(
            title = "Explanatory structure", value = "sr_explanatory",
            sliderInput("sx_cont", "Continuous-predictor effect", -2, 2, 1, 0.1),
            sliderInput("sx_cat", "Categorical-predictor effect", -2, 2, 0.6, 0.1),
            checkboxInput("sx_interaction", "Include their interaction", FALSE),
            conditionalPanel("input.sx_interaction == true",
              sliderInput("sx_int", "Interaction effect", -2, 2, 0.6, 0.1)),
            sliderInput("sx_depart", "Fixed departure for one item", 0, 2, 0, 0.1)))),
        accordion(open = FALSE, accordion_panel(
          title = span(bs_icon("bug"), " Misfit to plant"), value = "misfit",
          sliderInput("sr_over", "Over-discriminating items", 0, 5, 0, 1),
          sliderInput("sr_under", "Under-discriminating items", 0, 5, 0, 1),
          conditionalPanel("input.sr_model == 'dichotomous'",
            checkboxInput("sr_guess", "Guessing (lower asymptote)", FALSE)),
          checkboxInput("sr_2d", "Second dimension", FALSE),
          conditionalPanel("input.sr_2d",
            sliderInput("sr_rho", "correlation between dimensions", 0, 0.9, 0.3, 0.1)),
          checkboxInput("sr_dep", "Local dependence (item pair)", FALSE),
          checkboxInput("sr_dif", "DIF between two groups", FALSE),
          conditionalPanel("input.sr_dif",
            sliderInput("sr_difmag", "DIF magnitude (logits)", 0.3, 2, 1, 0.1)),
          conditionalPanel("input.sr_model != 'dichotomous'",
            checkboxInput("sr_style", "Response style", FALSE),
            conditionalPanel("input.sr_style",
              selectInput("sr_styletype", "style",
                c("Extreme categories" = "extreme",
                  "Middle categories" = "middle")))),
          sliderInput("sr_speeded", "Speededness (not-reached)", 0, 0.6, 0, 0.05),
          sliderInput("sr_careless", "Careless responders", 0, 0.3, 0, 0.02),
          sliderInput("sr_missing", "Missing data", 0, 0.3, 0, 0.02)))),
      # ---------------------------------------------------------- BTL ----
      conditionalPanel("input.sim_layout == 'btl' || input.sim_layout == 'btl_exp'",
        sliderInput("sb_objects", "Objects", 4, 20, 8, 1),
        sliderInput("sb_judges", "Judges", 3, 30, 12, 1),
        sliderInput("sb_reps", "Comparisons per pair", 5, 60, 25, 5),
        radioButtons("sb_model", "Verdict", inline = TRUE,
          c("Winner" = "dichotomous", "Polytomous margin" = "polytomous")),
        conditionalPanel("input.sb_model == 'polytomous'",
          sliderInput("sb_cats", "Categories", 3, 6, 4, 1)),
        sliderInput("sb_objsd", "Object-location SD", 0.5, 2.5, 1, 0.1),
        conditionalPanel("input.sim_layout == 'btl_exp'",
          accordion(open = FALSE, accordion_panel(
            title = "Explanatory structure", value = "sb_explanatory",
            sliderInput("sbe_cont", "Continuous-predictor effect", -2, 2, 1, 0.1),
            sliderInput("sbe_cat", "Categorical-predictor effect", -2, 2, 0.6, 0.1),
            checkboxInput("sbe_interaction", "Include their interaction", FALSE),
            conditionalPanel("input.sbe_interaction == true",
              sliderInput("sbe_int", "Interaction effect", -2, 2, 0.6, 0.1)),
            sliderInput("sbe_depart", "Fixed departure for one object", 0, 2, 0, 0.1)))),
        accordion(open = FALSE, accordion_panel(
          title = span(bs_icon("bug"), " Misfit to plant"), value = "misfit",
          sliderInput("sb_erratic", "Erratic judges", 0, 0.4, 0, 0.05),
          checkboxInput("sb_2a", "Second object attribute (two judge camps)", FALSE),
          conditionalPanel("input.sb_2a",
            sliderInput("sb_rho", "attribute correlation", 0, 0.9, 0.2, 0.1)),
          checkboxInput("sb_dep", "Within-judge dependence (order effects)", FALSE)))),
      # ---------------------------------------------------------- MFRM ----
      conditionalPanel("input.sim_layout == 'mfrm'",
        sliderInput("sm_persons", "Persons", 30, 300, 80, 10),
        sliderInput("sm_items", "Items", 3, 12, 5, 1),
        sliderInput("sm_raters", "Raters", 3, 15, 6, 1),
        sliderInput("sm_cats", "Categories", 3, 6, 4, 1),
        sliderInput("sm_thsd", "Person SD (logits)", 0.5, 2.5, 1.2, 0.1),
        sliderInput("sm_itemsd", "Item-difficulty SD", 0.3, 2, 1, 0.1),
        accordion(open = FALSE, accordion_panel(
          title = span(bs_icon("bug"), " Misfit to plant"), value = "misfit",
          sliderInput("sm_sev", "Rater-severity SD", 0, 1.5, 0.6, 0.1),
          sliderInput("sm_erratic", "Erratic raters", 0, 0.4, 0, 0.05),
          sliderInput("sm_halo", "Halo raters", 0, 0.4, 0, 0.05),
          checkboxInput("sm_int", "Rater-by-item interaction", FALSE)))),
      # ---------------------------------------------------------- EFRM ----
      conditionalPanel("input.sim_layout == 'efrm'",
        sliderInput("se_pergroup", "Persons per group", 100, 800, 300, 50),
        sliderInput("se_items", "Items per set", 5, 15, 8, 1),
        sliderInput("se_sets", "Item sets", 2, 4, 2, 1),
        sliderInput("se_groups", "Person groups", 2, 4, 2, 1),
        sliderInput("se_cats", "Response categories", 2, 6, 2, 1),
        sliderInput("se_thsd", "Person SD (logits)", 0.8, 2.5, 1.3, 0.1),
        sliderInput("se_setratio", "Set-unit ratio", 1, 2, 1.3, 0.05),
        sliderInput("se_grpratio", "Group-unit ratio", 1, 2, 1, 0.05),
        accordion(open = FALSE, accordion_panel(
          title = span(bs_icon("bug"), " Misfit to plant"), value = "se_misfit",
          checkboxInput("se_drift", "Item drift in one group", FALSE),
          conditionalPanel("input.se_drift == true",
            sliderInput("se_driftmag", "Item drift (logits)", 0.3, 2, 1, 0.1)),
          sliderInput("se_careless", "Careless responders", 0, 0.3, 0, 0.02),
          sliderInput("se_missing", "Missing data", 0, 0.3, 0, 0.02)))),
      # ----------------------------------------------------- BTL-EFRM ----
      conditionalPanel("input.sim_layout == 'btl_efrm'",
        sliderInput("sbf_objects", "Objects per set", 3, 12, 6, 1),
        sliderInput("sbf_sets", "Object sets", 2, 4, 2, 1),
        sliderInput("sbf_judges", "Judges per panel", 3, 15, 6, 1),
        sliderInput("sbf_panels", "Judge panels", 2, 4, 2, 1),
        sliderInput("sbf_within", "Comparisons per within-set pair", 5, 40, 20, 5),
        sliderInput("sbf_cross", "Comparisons per cross-set pair", 5, 40, 20, 5),
        sliderInput("sbf_objsd", "Object-location SD", 0.5, 2.5, 1, 0.1),
        sliderInput("sbf_setratio", "Set-unit ratio", 1, 2, 1.3, 0.05),
        sliderInput("sbf_panelratio", "Panel-unit ratio", 1, 2, 1, 0.05),
        sliderInput("sbf_origin", "Set-origin span", 0, 2, 0.5, 0.1),
        accordion(open = FALSE, accordion_panel(
          title = span(bs_icon("bug"), " Misfit to plant"), value = "sbf_misfit",
          sliderInput("sbf_erratic", "Erratic judges", 0, 0.4, 0, 0.05)))),
      hr(),
      numericInput("sim_seed", "Random seed", 1, min = 1),
      input_task_button("sim_go", "Simulate & load", type = "primary",
                        class = "w-100")),
    # ----------------------------------------------------------- main ----
    card(card_body(
      p(class = "lead mb-1", "Generate data with known parameters and optional departures."),
      p(class = "text-muted small mb-0",
        "Pick a data type, set the size and any departures, then select Simulate & load. The data loads with its roles assigned; go to Data and select Estimate to fit it. The planted values are listed below and attached to the data.")),
    uiOutput("sim_truth"),
    conditionalPanel("output.sim_ready == 'yes'", rcode_details("sim_truth")),
    conditionalPanel("output.sim_has_recovery == 'yes'",
      card(class = "mt-3", full_screen = TRUE,
        id = "card_sim_recovery",
        card_header(span(bs_icon("bullseye"), " Generating and fitted parameters",
                         info_icon(app_help("sim_recovery_tbl"),
                                   "About this table"))),
        card_body(
          p(class = "text-muted small mb-2",
            "Planted and fitted values after estimation. Locations are mean-centred because the model fixes only the origin."),
          uiOutput("sim_recovery_note"),
          DTOutput("sim_recovery_tbl"),
          plotOutput("sim_recovery_plot", height = "300px"),
          rcode_details("sim_recovery"), padding = 10))),
    conditionalPanel("output.sim_ready == 'yes'",
      card(class = "mt-3",
        card_header(span(bs_icon("code-slash"), " Reproducible code")),
        card_body(verbatimTextOutput("sim_code"), padding = 10)),
      card(class = "mt-3", full_screen = TRUE,
        id = "card_sim_preview",
        card_header_bar("Preview of the loaded data (first rows)",
          info = app_help("sim_preview"),
          buttons = div(class = "d-flex gap-1",
            downloadButton("sim_data_csv", "Data (CSV)",
                           class = "btn-outline-secondary btn-xs"),
            downloadButton("sim_bundle_zip", "Data + setup (ZIP)",
                           class = "btn-outline-secondary btn-xs"))),
        card_body(DTOutput("sim_preview"), rcode_details("sim_preview"),
                  padding = 8))))
  ))

# ------------------------------------------------------------ ASSEMBLY --

# Workflow order: data -> summary -> items -> persons -> test, then the
# independence, invariance, and utility menus (the two requirements of
# measurement); status chips and the dark-mode toggle sit at the right of
# the navbar.
ui <- page_navbar(
  id = "nav",
  title = span(class = "app-brand",
               span("rasch", class = "app-brand-name"),
               span("measurement theory", class = "app-brand-sub")),
  theme = theme,
  # normal scrolling pages: never compress content to fit the viewport
  fillable = FALSE,
  header = tagList(
    tags$head(tags$style(css),
      # nav visibility by data-value: shiny::hideTab (behind bslib's
      # nav_hide) does not reach nav_panels nested inside a nav_menu
      # dropdown, so the server toggles entries itself through this handler;
      # it covers top-level links, dropdown items, and menu toggles alike
      tags$script(HTML("
        Shiny.addCustomMessageHandler('rasch-nav-vis', function(msg) {
          document.querySelectorAll('.navbar a[data-value]').forEach(function(a) {
            if (a.getAttribute('data-value') !== msg.value) return;
            var li = a.closest('li');
            (li || a).style.display = msg.show ? '' : 'none';
          });
        });
      "))),
    busyIndicatorOptions(spinner_type = "ring2")),
  panel_data,
  panel_summary,
  panel_items,
  panel_explanatory,
  panel_persons,
  panel_targeting,
  nav_menu("Independence", value = "menu_independence",
    panel_ld,
    panel_dim),
  nav_menu("Invariance", value = "menu_invariance",
    panel_dif,
    panel_equating,
    panel_guess,
    panel_facets,
    panel_frames),
  nav_menu("More", value = "menu_more",
    panel_simulate,
    panel_compare,
    panel_export),
  nav_spacer(),
  nav_item(uiOutput("nav_status")),
  nav_item(downloadLink("dl_report_nav", label = bs_icon("file-earmark-text"),
                        class = "nav-link px-2",
                        title = "Analysis report (HTML)")),
  nav_item(input_dark_mode())
)

server <- function(input, output, session) {

  # ------------------------------------------------------------- data in --
  # A monotonically increasing token identifies the active analysis context.
  # Background workers retain the token and their input data at launch; their
  # results are accepted only while both still match. This prevents a worker
  # started for one analysis from updating a later fit or reopened project.
  analysis_context <- reactiveVal(0L)
  advance_analysis_context <- function() {
    z <- isolate(analysis_context())
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) ||
        z >= .Machine$integer.max) z <- 0L
    z <- as.integer(z) + 1L
    analysis_context(z)
    invisible(z)
  }

  # simulated data (from the Simulate page) takes precedence over a demo or
  # an upload until one of those replaces it
  sim_data <- reactiveVal(NULL)
  sim_truth_val <- reactiveVal(NULL)
  sim_code_val <- reactiveVal(NULL)
  sim_predictors_val <- reactiveVal(NULL)
  sim_interactions_val <- reactiveVal(character(0))
  restored_project_name <- reactiveVal(NULL)
  restored_project_settings <- reactiveVal(list())
  restored_project_resources <- reactiveVal(list())
  person_weight_state <- reactiveVal(NULL)
  restored_dimensionality <- reactiveVal(NULL)
  restored_subtest <- reactiveVal(NULL)
  restored_invariance <- reactiveVal(NULL)
  # generation stamps: which simulation is loaded, and which one the current
  # fit was estimated on (recovery only renders when they agree)
  sim_gen <- reactiveVal(0L)
  fitted_sim_gen <- reactiveVal(NULL)
  # picking an example dataset also selects the matching model; uploading a
  # file clears the example selection
  observeEvent(input$demo_choice, {
    dc <- input$demo_choice
    invalidate_source_results("The example dataset selection changed")
    if (!identical(dc, "none")) {
      sim_data(NULL); sim_truth_val(NULL); sim_code_val(NULL)
      sim_predictors_val(NULL); sim_interactions_val(character(0))
      restored_project_name(NULL)
      restored_project_settings(list()); restored_project_resources(list())
      updateRadioButtons(session, "model_type",
                         selected = if (dc %in% c("dich", "pcm", "rsm"))
                           "rasch" else dc)
      if (dc %in% c("dich", "pcm", "rsm"))
        updateRadioButtons(session, "thr_structure",
                           selected = if (identical(dc, "rsm")) "rsm" else "pcm")
    }
  }, ignoreInit = TRUE)
  observeEvent(input$file, {
    invalidate_source_results("A new data file was uploaded")
    sim_data(NULL); sim_truth_val(NULL); sim_code_val(NULL)
    sim_predictors_val(NULL); sim_interactions_val(character(0))
    restored_project_name(NULL)
    restored_project_settings(list()); restored_project_resources(list())
    updateSelectInput(session, "demo_choice", selected = "none")
  })
  observeEvent(input$exp_predictors, {
    if (!is.null(input$exp_predictors)) {
      sim_predictors_val(NULL)
      sim_interactions_val(character(0))
      # Dynamic predictor controls from a reopened project belong to its
      # embedded metadata, not to a newly uploaded predictor file.
      s <- restored_project_settings()
      old_dynamic <- grep("^exp_(type|ref|order)_[0-9]+$", names(s),
                          value = TRUE)
      if (length(old_dynamic)) {
        s[old_dynamic] <- NULL
        restored_project_settings(s)
      }
    }
  }, ignoreInit = TRUE)

  raw_data <- reactive({
    if (!is.null(sim_data())) return(sim_data())
    if (!identical(input$demo_choice %||% "none", "none"))
      return(switch(input$demo_choice,
                    dich = .demo_dich(), rsm = .demo_rsm(),
                    mfrm = .demo_mfrm(), efrm = .demo_efrm(),
                    btl = .demo_btl(), .demo_data()))
    req(input$file)
    ext <- tolower(tools::file_ext(input$file$name))
    sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
    read.csv(input$file$datapath, sep = sep, check.names = FALSE,
             stringsAsFactors = FALSE)
  })

  # ---- Simulate page: build the call, generate, and load as current data --
  observeEvent(input$sim_go, {
    lay <- input$sim_layout
    seed_raw <- suppressWarnings(as.numeric(input$sim_seed %||% NA_real_))
    if (length(seed_raw) != 1L || !is.finite(seed_raw) || seed_raw < 0 ||
        seed_raw != floor(seed_raw) || seed_raw > .Machine$integer.max) {
      showNotification(
        "Seed must be one non-negative whole number within the integer range",
        type = "warning")
      return()
    }
    seed <- as.integer(seed_raw)
    sim_call <- NULL
    predictors <- NULL
    wanted_interactions <- character(0)
    nums <- function(x) paste(format(as.numeric(x), digits = 15, trim = TRUE,
                                     scientific = FALSE), collapse = ", ")
    d <- tryCatch(withProgress(message = "Simulating…", value = 0.5, {
      if (lay %in% c("rasch", "rasch_exp")) {
        I <- input$sr_items; disc <- rep(1, I)
        no <- min(input$sr_over, I); nu <- min(input$sr_under, I - no)
        if (no + nu >= I && no + nu > 0)
          stop("Over- and under-discriminating items must leave at least one ",
               "item with the reference discrimination.")
        if (no > 0) disc[seq_len(no)] <- 2.8
        if (nu > 0) disc[no + seq_len(nu)] <- 0.4
        guess <- if (isTRUE(input$sr_guess) && input$sr_model == "dichotomous") {
          g <- rep(0, I); g[seq_len(max(1L, round(0.4 * I)))] <- 0.25; g } else 0
        k2 <- max(2L, round(0.4 * I))
        s2 <- if (isTRUE(input$sr_2d))
          list(items = sprintf("I%02d", (I - k2 + 1L):I), rho = input$sr_rho) else NULL
        dep <- if (isTRUE(input$sr_dep) && I >= 2)
          list(pairs = list(c("I01", "I02")), strength = 2.2) else NULL
        dif <- if (isTRUE(input$sr_dif))
          list(items = sprintf("I%02d", max(1L, round(I / 2))),
               uniform = input$sr_difmag) else NULL
        rstyle <- if (isTRUE(input$sr_style) && input$sr_model != "dichotomous")
          list(type = input$sr_styletype %||% "extreme", prop = 0.25) else NULL
        difficulty <- input$sr_diff %||% c(-2.5, 2.5)
        if (lay == "rasch_exp") {
          item <- sprintf("I%02d", seq_len(I))
          exposure <- seq(-1, 1, length.out = I)
          type <- rep(c("A", "B"), length.out = I)
          difficulty <- (input$sx_cont %||% 1) * exposure +
            (input$sx_cat %||% 0.6) * (type == "B")
          if (isTRUE(input$sx_interaction)) {
            difficulty <- difficulty + (input$sx_int %||% 0.6) *
              exposure * (type == "B")
            wanted_interactions <- "exposure:type"
          }
          design_formula <- if (isTRUE(input$sx_interaction))
            ~ exposure * type else ~ exposure + type
          depart_item <- item[max(1L, round(I / 2))]
          if ((input$sx_depart %||% 0) > 0) {
            departure <- .sim_explanatory_departure(
              model.matrix(design_formula), input$sx_depart)
            difficulty <- difficulty + departure$values
            depart_item <- item[departure$index]
          }
          difficulty <- difficulty - mean(difficulty)
          predictors <- data.frame(item = item, exposure = exposure,
                                   type = type, stringsAsFactors = FALSE)
        }
        out <- simulate_rasch(input$sr_persons, I, model = input$sr_model,
          n_categories = input$sr_cats %||% 4L,
          theta_mean = input$sr_mean %||% 0, theta_sd = input$sr_sd %||% 1,
          theta_dist = input$sr_dist %||% "normal",
          difficulty = difficulty,
          discrimination = disc,
          guessing = guess, second_dim = s2, dependence = dep, dif = dif,
          n_groups = if (!is.null(dif)) 2L else 1L,
          careless = input$sr_careless, response_style = rstyle,
          speeded = input$sr_speeded %||% 0, missing = input$sr_missing,
          seed = seed)
        a <- c(sprintf("n_persons = %d, n_items = %d", input$sr_persons, I),
          sprintf('model = "%s"', input$sr_model),
          sprintf("theta_mean = %s, theta_sd = %s", input$sr_mean, input$sr_sd),
          if (input$sr_dist != "normal") sprintf('theta_dist = "%s"', input$sr_dist),
          sprintf("difficulty = c(%s)", nums(difficulty)),
          if (input$sr_model != "dichotomous")
            sprintf("n_categories = %d", input$sr_cats %||% 4L),
          if (input$sr_over > 0 || input$sr_under > 0)
            sprintf("discrimination = c(%s)", nums(disc)),
          if (any(guess != 0)) sprintf("guessing = c(%s)", nums(guess)),
          if (!is.null(s2)) sprintf('second_dim = list(items = c(%s), rho = %s)',
            paste0('"', s2$items, '"', collapse = ", "), s2$rho),
          if (!is.null(dep))
            'dependence = list(pairs = list(c("I01", "I02")), strength = 2.2)',
          if (!is.null(dif)) sprintf(
            'dif = list(items = "%s", uniform = %s), n_groups = 2',
            dif$items, dif$uniform),
          if (!is.null(rstyle)) sprintf(
            'response_style = list(type = "%s", prop = 0.25)', rstyle$type),
          if (input$sr_speeded > 0) sprintf("speeded = %s", input$sr_speeded),
          if (input$sr_careless > 0) sprintf("careless = %s", input$sr_careless),
          if (input$sr_missing > 0) sprintf("missing = %s", input$sr_missing))
        base_call <- sprintf("simulate_rasch(\n  %s,\n  seed = %d)",
                             paste(Filter(nzchar, a), collapse = ",\n  "), seed)
        if (lay == "rasch_exp") {
          sim_call <- paste0("local({\n  predictors <- data.frame(\n",
            "    item = ", qvec(predictors$item), ",\n",
            "    exposure = c(", nums(predictors$exposure), "),\n",
            "    type = ", qvec(predictors$type), ",\n",
            "    stringsAsFactors = FALSE)\n  d <- ",
            gsub("\n", "\n  ", base_call),
            "\n  attr(d, \"predictors\") <- predictors\n",
            "  truth <- attr(d, \"truth\")\n",
            "  truth$predictors <- predictors\n",
            "  truth$explanatory_formula <- \"",
            if (length(wanted_interactions))
              "~ exposure + type + exposure:type" else "~ exposure + type",
            "\"\n",
            if ((input$sx_depart %||% 0) > 0) sprintf(
              paste0(
                "  truth$planted <- c(truth$planted, \"fixed explanatory departure for %s (%.2f logits)\")\n",
                "  truth$departure_types <- c(truth$departure_types, \"a fixed explanatory item departure\")\n",
                "  truth$explanatory_departure <- list(target = \"%s\", level = \"item\")\n"),
              depart_item, input$sx_depart, depart_item) else "",
            "  attr(d, \"truth\") <- truth\n  d\n})")
          tr <- attr(out, "truth")
          tr$predictors <- predictors
          tr$explanatory_formula <- if (length(wanted_interactions))
            "~ exposure + type + exposure:type" else "~ exposure + type"
          if ((input$sx_depart %||% 0) > 0) {
            tr$planted <- c(tr$planted,
              sprintf("fixed explanatory departure for %s (%.2f logits)",
                      depart_item, input$sx_depart))
            tr$departure_types <- c(tr$departure_types,
                                    "a fixed explanatory item departure")
            tr$explanatory_departure <- list(target = depart_item,
                                             level = "item")
          }
          attr(out, "truth") <- tr
          attr(out, "predictors") <- predictors
        } else sim_call <- base_call
        out
      } else if (lay %in% c("btl", "btl_exp")) {
        s2 <- if (isTRUE(input$sb_2a)) list(rho = input$sb_rho) else NULL
        dep <- if (isTRUE(input$sb_dep)) list(exposure = 0.6, carry_over = 1) else NULL
        object_locations <- NULL
        if (lay == "btl_exp") {
          K <- input$sb_objects; object <- sprintf("O%d", seq_len(K))
          exposure <- seq(-1, 1, length.out = K)
          type <- rep(c("A", "B"), length.out = K)
          object_locations <- (input$sbe_cont %||% 1) * exposure +
            (input$sbe_cat %||% 0.6) * (type == "B")
          if (isTRUE(input$sbe_interaction)) {
            object_locations <- object_locations + (input$sbe_int %||% 0.6) *
              exposure * (type == "B")
            wanted_interactions <- "exposure:type"
          }
          design_formula <- if (isTRUE(input$sbe_interaction))
            ~ exposure * type else ~ exposure + type
          depart_object <- object[max(1L, round(K / 2))]
          if ((input$sbe_depart %||% 0) > 0) {
            departure <- .sim_explanatory_departure(
              model.matrix(design_formula), input$sbe_depart)
            object_locations <- object_locations + departure$values
            depart_object <- object[departure$index]
          }
          predictors <- data.frame(object = object, exposure = exposure,
                                   type = type, stringsAsFactors = FALSE)
        }
        out <- simulate_btl(input$sb_objects, input$sb_judges, input$sb_reps,
          model = input$sb_model, n_categories = input$sb_cats %||% 4L,
          object_sd = input$sb_objsd %||% 1,
          object_locations = object_locations,
          second_attribute = s2, erratic_judges = input$sb_erratic,
          dependence = dep, seed = seed)
        base_call <- sprintf(
          'simulate_btl(%d, %d, %d, model = "%s"%s, object_sd = %s%s%s%s%s,\n  seed = %d)',
          input$sb_objects, input$sb_judges, input$sb_reps, input$sb_model,
          if (identical(input$sb_model, "polytomous"))
            sprintf(", n_categories = %d", input$sb_cats %||% 4L) else "",
          input$sb_objsd,
          if (!is.null(object_locations))
            sprintf(", object_locations = c(%s)", nums(object_locations)) else "",
          if (input$sb_erratic > 0)
            sprintf(", erratic_judges = %s", input$sb_erratic) else "",
          if (isTRUE(input$sb_2a))
            sprintf(", second_attribute = list(rho = %s)", input$sb_rho) else "",
          if (isTRUE(input$sb_dep))
            ", dependence = list(exposure = 0.6, carry_over = 1)" else "", seed)
        if (lay == "btl_exp") {
          sim_call <- paste0("local({\n  predictors <- data.frame(\n",
            "    object = ", qvec(predictors$object), ",\n",
            "    exposure = c(", nums(predictors$exposure), "),\n",
            "    type = ", qvec(predictors$type), ",\n",
            "    stringsAsFactors = FALSE)\n  d <- ",
            gsub("\n", "\n  ", base_call),
            "\n  attr(d, \"predictors\") <- predictors\n",
            "  truth <- attr(d, \"truth\")\n",
            "  truth$predictors <- predictors\n",
            "  truth$explanatory_formula <- \"",
            if (length(wanted_interactions))
              "~ exposure + type + exposure:type" else "~ exposure + type",
            "\"\n",
            if ((input$sbe_depart %||% 0) > 0) sprintf(
              paste0(
                "  truth$planted <- c(truth$planted, \"fixed explanatory departure for %s (%.2f logits)\")\n",
                "  truth$departure_types <- c(truth$departure_types, \"a fixed explanatory object departure\")\n",
                "  truth$explanatory_departure <- list(target = \"%s\", level = \"object\")\n"),
              depart_object, input$sbe_depart, depart_object) else "",
            "  attr(d, \"truth\") <- truth\n  d\n})")
          tr <- attr(out, "truth")
          tr$predictors <- predictors
          tr$explanatory_formula <- if (length(wanted_interactions))
            "~ exposure + type + exposure:type" else "~ exposure + type"
          if ((input$sbe_depart %||% 0) > 0) {
            tr$planted <- c(tr$planted,
              sprintf("fixed explanatory departure for %s (%.2f logits)",
                      depart_object, input$sbe_depart))
            tr$departure_types <- c(tr$departure_types,
                                    "a fixed explanatory object departure")
            tr$explanatory_departure <- list(target = depart_object,
                                             level = "object")
          }
          attr(out, "truth") <- tr
          attr(out, "predictors") <- predictors
        } else sim_call <- base_call
        out
      } else if (lay == "mfrm") {
        # The simulator assigns erratic raters from R1 upwards. Plant the
        # interaction on the first remaining rater so selecting both controls
        # does not put a model-based bias on a rater whose responses are then
        # replaced by random ratings.
        n_erratic <- .sim_planted_count(input$sm_erratic %||% 0,
                                        input$sm_raters)
        non_erratic <- setdiff(sprintf("R%d", seq_len(input$sm_raters)),
                               sprintf("R%d", seq_len(n_erratic)))
        if (isTRUE(input$sm_int) && !length(non_erratic))
          stop("a rater-by-item interaction needs at least one non-erratic rater")
        interaction_rater <- if (length(non_erratic)) non_erratic[1L] else NA_character_
        intr <- if (isTRUE(input$sm_int))
          list(rater = interaction_rater, item = "I2", bias = 1.8) else NULL
        out <- simulate_mfrm(input$sm_persons, input$sm_items, input$sm_raters,
          n_categories = input$sm_cats, theta_sd = input$sm_thsd %||% 1.2,
          item_sd = input$sm_itemsd %||% 1, rater_severity_sd = input$sm_sev,
          erratic_raters = input$sm_erratic, halo = input$sm_halo %||% 0,
          interaction = intr, seed = seed)
        sim_call <- sprintf(
          'simulate_mfrm(%d, %d, %d, n_categories = %d, theta_sd = %s, item_sd = %s, rater_severity_sd = %s%s%s%s,\n  seed = %d)',
          input$sm_persons, input$sm_items, input$sm_raters, input$sm_cats,
          input$sm_thsd %||% 1.2, input$sm_itemsd %||% 1, input$sm_sev,
          if (input$sm_erratic > 0)
            sprintf(", erratic_raters = %s", input$sm_erratic) else "",
          if (input$sm_halo > 0) sprintf(", halo = %s", input$sm_halo) else "",
          if (isTRUE(input$sm_int))
            sprintf(', interaction = list(rater = "%s", item = "I2", bias = 1.8)',
                    interaction_rater) else "",
          seed)
        out
      } else if (lay == "efrm") {
        drift <- if (isTRUE(input$se_drift))
          list(items = sprintf("S1I%02d", max(1L, round(input$se_items / 2))),
               group = "g2", shift = input$se_driftmag) else NULL
        out <- simulate_efrm(input$se_pergroup, input$se_items, input$se_sets,
          input$se_groups, set_unit_ratio = input$se_setratio,
          group_unit_ratio = input$se_grpratio,
          n_categories = input$se_cats %||% 2L,
          theta_sd = input$se_thsd %||% 1.3, item_drift = drift,
          careless = input$se_careless %||% 0,
          missing = input$se_missing %||% 0, seed = seed)
        sim_call <- sprintf(
          'simulate_efrm(%d, %d, %d, %d, set_unit_ratio = %s, group_unit_ratio = %s, n_categories = %d, theta_sd = %s%s%s%s,\n  seed = %d)',
          input$se_pergroup, input$se_items, input$se_sets, input$se_groups,
          input$se_setratio, input$se_grpratio, input$se_cats %||% 2L,
          input$se_thsd %||% 1.3,
          if (!is.null(drift)) sprintf(
            ', item_drift = list(items = "%s", group = "g2", shift = %s)',
            drift$items, drift$shift) else "",
          if ((input$se_careless %||% 0) > 0)
            sprintf(", careless = %s", input$se_careless) else "",
          if ((input$se_missing %||% 0) > 0)
            sprintf(", missing = %s", input$se_missing) else "", seed)
        out
      } else {
        span_units <- function(ratio, n) exp(seq(0, log(ratio), length.out = n))
        panel_units <- span_units(input$sbf_panelratio, input$sbf_panels)
        set_units <- span_units(input$sbf_setratio, input$sbf_sets)
        set_origins <- seq(0, input$sbf_origin, length.out = input$sbf_sets)
        out <- simulate_btl_efrm(
          input$sbf_objects, input$sbf_sets, input$sbf_judges,
          input$sbf_panels, input$sbf_within, input$sbf_cross,
          panel_units = panel_units, set_units = set_units,
          set_origins = set_origins, object_sd = input$sbf_objsd,
          erratic_judges = input$sbf_erratic, seed = seed)
        sim_call <- sprintf(
          paste0('simulate_btl_efrm(%d, %d, %d, %d, %d, %d,\n',
            '  panel_units = c(%s), set_units = c(%s),\n',
            '  set_origins = c(%s), object_sd = %s%s, seed = %d)'),
          input$sbf_objects, input$sbf_sets, input$sbf_judges,
          input$sbf_panels, input$sbf_within, input$sbf_cross,
          nums(panel_units), nums(set_units), nums(set_origins),
          input$sbf_objsd,
          if (input$sbf_erratic > 0)
            sprintf(", erratic_judges = %s", input$sbf_erratic) else "", seed)
        out
      }
    }), error = function(e) e)
    if (inherits(d, "error")) {
      showNotification(paste("Simulation failed:", conditionMessage(d)),
                       type = "error", duration = 10); return()
    }
    sim_truth_val(attr(d, "truth"))
    sim_predictors_val(predictors)
    sim_interactions_val(wanted_interactions)
    sim_code_val(sim_call)
    restored_project_name(NULL)
    restored_project_settings(list())
    restored_project_resources(list())
    updateSelectInput(session, "demo_choice", selected = "none")
    model_selected <- if (lay %in% c("rasch", "rasch_exp")) "rasch" else
      if (lay %in% c("btl", "btl_exp", "btl_efrm")) "btl" else lay
    updateRadioButtons(session, "model_type", selected = model_selected)
    if (lay %in% c("rasch", "rasch_exp"))
      updateRadioButtons(session, "thr_structure",
        selected = if (input$sr_model == "RSM") "rsm" else "pcm")
    updateRadioButtons(session, "rasch_calibration",
      selected = if (lay %in% c("rasch_exp", "btl_exp")) "explanatory" else "free")
    if (lay == "rasch_exp")
      updateRadioButtons(session, "exp_level", selected = "item")
    # simulated many-facet data is long (person, item, rater, score)
    if (lay == "mfrm")
      updateRadioButtons(session, "lp_layout", selected = "long")
    invalidate_source_results("A new simulation was loaded")
    sim_gen(sim_gen() + 1L)
    sim_data(as.data.frame(d))   # plain frame -> raw_data() -> role guessing
    showNotification("Simulated data loaded. Go to Data and select Estimate.",
                     type = "message", duration = 8)
  })
  output$sim_ready <- reactive(if (!is.null(sim_truth_val())) "yes" else "no")
  outputOptions(output, "sim_ready", suspendWhenHidden = FALSE)
  output$sim_truth <- renderUI({
    tr <- sim_truth_val(); if (is.null(tr)) return(NULL)
    card(class = "mt-3",
      card_header(span(bs_icon("check-circle-fill", class = "text-success"),
                       " Loaded")),
      card_body(
        p(class = "mb-1", strong(tr$description)),
        if (length(tr$planted)) tagList(
          p(class = "mb-1 small text-muted", "Generating values and departures:"),
          tags$ul(class = "small mb-1", lapply(tr$planted, tags$li)))
        else p(class = "small text-muted", "Model-conforming (no departures)."),
        p(class = "small mb-0", bs_icon("arrow-right-circle"),
          " Go to Data and select Estimate to examine the matching diagnostic.")))
  })
  output$sim_preview <- renderDT({
    req(!is.null(sim_data()))
    d <- head(sim_data(), 12)
    round_preview(
      datatable(d, rownames = FALSE, style = "bootstrap5", class = "table-sm compact",
                options = list(dom = "t", scrollX = TRUE)), d)
  })
  sim_file_stem <- function() paste0(
    "rasch_sim_", sim_truth_val()$layout %||% "data", "_",
    format(Sys.Date(), "%Y%m%d"))
  output$sim_data_csv <- downloadHandler(
    filename = function() paste0(sim_file_stem(), ".csv"),
    content = function(file) {
      req(!is.null(sim_data()))
      write_csv_plain(sim_data(), file)
    })
  write_sim_bundle <- function(file) {
    req(!is.null(sim_data()), nzchar(sim_code_val() %||% ""))
    bundle <- tempfile("rasch-simulation-")
    dir.create(bundle)
    on.exit(unlink(bundle, recursive = TRUE, force = TRUE), add = TRUE)
    write_csv_plain(sim_data(), file.path(bundle, "data.csv"))
    writeLines(c("# Recreate the simulated data and its truth attribute",
                 "library(rasch)",
                 paste0("data <- ", sim_code_val())),
               file.path(bundle, "simulation.R"))
    saveRDS(sim_truth_val(), file.path(bundle, "truth.rds"))
    if (!is.null(sim_predictors_val()))
      write_csv_plain(sim_predictors_val(),
                      file.path(bundle, "predictors.csv"))
    writeLines(c(
      "data.csv contains the simulated analysis data.",
      "simulation.R recreates the classed dataset and its truth attribute.",
      "truth.rds contains the exact generating values and planted departures.",
      if (!is.null(sim_predictors_val()))
        "predictors.csv contains the explanatory item or object metadata."
      else NULL), file.path(bundle, "README.txt"))
    old <- setwd(bundle); on.exit(setwd(old), add = TRUE)
    utils::zip(file, list.files(".", all.files = FALSE), flags = "-9Xq")
    invisible(file)
  }
  output$sim_bundle_zip <- downloadHandler(
    filename = function() paste0(sim_file_stem(), ".zip"),
    content = write_sim_bundle)
  output$sim_code <- renderText(sim_code_val() %||% "")
  # parameter recovery: the fit matching the simulated layout, compared to the
  # attached truth (re-attached, since sim_data() is stored as a plain frame)
  sim_recovery_val <- reactive({
    tr <- sim_truth_val(); req(!is.null(tr))
    # only compare a fit estimated on THIS simulation; a fit left over from
    # an earlier one would render misleading recovery numbers
    req(identical(fitted_sim_gen(), sim_gen()))
    f <- if (identical(tr$layout, "btl")) btl_fit() else
      if (identical(tr$layout, "btl_efrm")) {
        z <- bfit()
        req(inherits(z, "rasch_btl_efrm"))
        z
      } else fit_or_null()
    req(!is.null(f))
    obj <- sim_data(); attr(obj, "truth") <- tr
    tryCatch(sim_recovery(f, obj), error = function(e) NULL)
  })
  output$sim_has_recovery <- reactive({
    ok <- tryCatch(!is.null(sim_recovery_val()), error = function(e) FALSE)
    if (ok) "yes" else "no"
  })
  outputOptions(output, "sim_has_recovery", suspendWhenHidden = FALSE)
  output$sim_recovery_tbl <- renderDT({
    r <- sim_recovery_val(); req(!is.null(r)); num_dt(r$summary)
  })
  output$sim_recovery_note <- renderUI({
    r <- sim_recovery_val(); req(!is.null(r))
    if (is.null(r$note)) return(NULL)
    p(class = "text-warning-emphasis small mb-2", r$note)
  })
  output$sim_recovery_plot <- renderPlot({
    r <- sim_recovery_val(); req(!is.null(r)); plot_recovery(r)
  }, res = 96)

  anchors_in <- reactive({
    a <- if (!is.null(input$anchor_file))
      tryCatch(read.csv(input$anchor_file$datapath,
                        check.names = FALSE, stringsAsFactors = FALSE),
               error = function(e)
                 stop("could not read the anchor CSV: ",
                      conditionMessage(e), call. = FALSE))
    else restored_project_resources()$anchors
    if (is.null(a)) return(NULL)
    if (anyDuplicated(names(a)))
      stop("the anchor CSV has duplicate column names", call. = FALSE)
    if (!all(c("item", "k", "tau") %in% names(a)))
      stop("the anchor CSV needs columns item, k, tau", call. = FALSE)
    if (!nrow(a)) stop("the anchor CSV contains no anchor rows", call. = FALSE)
    missing_item <- is.na(a$item)
    a$item <- trimws(as.character(a$item))
    a$item[missing_item] <- NA_character_
    if (anyNA(a$item) || any(!nzchar(a$item)))
      stop("every anchor row needs a non-blank item name", call. = FALSE)
    a
  })

  key_in <- reactive({
    kf <- if (!is.null(input$key_file))
      tryCatch(read.csv(input$key_file$datapath,
                        check.names = FALSE, stringsAsFactors = FALSE),
               error = function(e)
                 stop("could not read the scoring-key CSV: ",
                      conditionMessage(e), call. = FALSE))
    else restored_project_resources()$key
    if (is.null(kf)) {
      if (identical(input$demo_choice, "dich")) return(.demo_dich_key())
      return(NULL)
    }
    if (anyDuplicated(names(kf)))
      stop("the scoring-key CSV has duplicate column names")
    if (!(all(c("item", "key") %in% names(kf)) ||
          all(c("item", "option", "score") %in% names(kf))))
      stop("the scoring-key CSV needs columns item,key or item,option,score")
    if (!nrow(kf)) stop("the scoring-key CSV contains no key rows")
    kf
  })

  exp_predictors_raw <- reactive({
    if (!is.null(sim_predictors_val())) {
      p <- sim_predictors_val()
    } else if (!is.null(input$exp_predictors)) {
      p <- tryCatch(read.csv(input$exp_predictors$datapath,
                             check.names = FALSE, stringsAsFactors = FALSE),
                    error = function(e) e)
      if (inherits(p, "error"))
        stop("could not read the predictor metadata: ", conditionMessage(p))
    } else {
      p <- restored_project_resources()$predictors
      if (is.null(p)) req(FALSE)
    }
    key <- if (identical(input$model_type, "btl")) "object" else "item"
    if (!key %in% names(p))
      stop("predictor metadata needs a ", key, " column")
    p
  })
  exp_predictor_vars <- reactive(setdiff(names(exp_predictors_raw()),
    c("item", "object", "threshold", "threshold_number")))
  exp_predictor_type <- function(p, nm) {
    j <- match(nm, exp_predictor_vars())
    id <- paste0("exp_type_", j)
    input[[id]] %||% restored_project_settings()[[id]] %||%
      if (is.numeric(p[[nm]])) "continuous" else "categorical"
  }
  exp_level_order <- function(p, nm) {
    j <- match(nm, exp_predictor_vars())
    id <- paste0("exp_order_", j)
    raw <- input[[id]] %||% restored_project_settings()[[id]] %||% ""
    lev <- trimws(strsplit(raw, "[>,;]", perl = TRUE)[[1]])
    lev <- lev[nzchar(lev)]
    observed <- unique(as.character(p[[nm]][!is.na(p[[nm]])]))
    if (!length(lev)) lev <- if (is.numeric(p[[nm]]))
      as.character(sort(unique(p[[nm]][!is.na(p[[nm]])]))) else observed
    if (anyDuplicated(lev) || !setequal(lev, observed))
      stop("ordinal order for ", nm,
           " must list every observed level exactly once")
    lev
  }
  exp_category_levels <- function(p, nm) {
    observed <- unique(as.character(p[[nm]][!is.na(p[[nm]])]))
    lev <- if (is.numeric(p[[nm]]))
      as.character(sort(unique(p[[nm]][!is.na(p[[nm]])]))) else
        sort(observed)
    if (!length(lev)) stop("predictor ", nm, " has no observed values")
    j <- match(nm, exp_predictor_vars())
    id <- paste0("exp_ref_", j)
    ref <- input[[id]] %||% restored_project_settings()[[id]] %||% lev[1]
    if (!ref %in% lev) stop("categorical reference for ", nm,
                             " is not an observed level")
    c(ref, setdiff(lev, ref))
  }
  exp_predictors_in <- reactive({
    p <- exp_predictors_raw()
    for (nm in exp_predictor_vars()) {
      type <- exp_predictor_type(p, nm)
      if (type == "ordinal") {
        p[[nm]] <- ordered(as.character(p[[nm]]),
                           levels = exp_level_order(p, nm))
      } else if (type == "categorical") {
        p[[nm]] <- factor(p[[nm]], levels = exp_category_levels(p, nm))
      } else {
        original_na <- is.na(p[[nm]])
        value <- suppressWarnings(as.numeric(as.character(p[[nm]])))
        if (anyNA(value) && any(!original_na & is.na(value)))
          stop("continuous predictor ", nm, " contains non-numeric values")
        p[[nm]] <- value
      }
    }
    p
  })
  output$exp_predictor_types <- renderUI({
    if (is.null(input$exp_predictors) && is.null(sim_predictors_val()) &&
        is.null(restored_project_resources()$predictors)) return(NULL)
    p <- tryCatch(exp_predictors_raw(), error = function(e) NULL)
    if (is.null(p)) return(NULL)
    vars <- exp_predictor_vars()
    tagList(
      div(class = "small text-muted mb-2",
          info_label("Predictor types",
            paste("Categorical predictors compare levels with a selected reference.",
                  "Ordinal predictors estimate adjacent changes along the",
                  "declared order. Continuous predictors estimate a linear",
                  "effect per unit."))),
      div(class = "rasch-predictor-types", lapply(vars, function(nm) {
        j <- match(nm, vars)
        type <- exp_predictor_type(p, nm)
        default <- if (is.numeric(p[[nm]]))
          sort(unique(p[[nm]][!is.na(p[[nm]])])) else
            unique(as.character(p[[nm]][!is.na(p[[nm]])]))
        div(class = "rasch-predictor-type",
          selectInput(paste0("exp_type_", j), nm,
            choices = c("Categorical" = "categorical",
                        "Ordinal" = "ordinal",
                        "Continuous" = "continuous"),
            selected = type),
          if (type == "categorical")
            selectInput(paste0("exp_ref_", j), "Reference level",
              choices = exp_category_levels(p, nm),
              selected = (input[[paste0("exp_ref_", j)]] %||%
                            exp_category_levels(p, nm)[1])),
          if (type == "ordinal")
            textInput(paste0("exp_order_", j), paste("Order for", nm),
              value = paste(default, collapse = " > "),
              placeholder = "lowest > middle > highest"))
      })))
  })
  exp_predictor_code <- reactive({
    p <- exp_predictors_raw()
    vapply(exp_predictor_vars(), function(nm) {
      type <- exp_predictor_type(p, nm)
      if (type == "ordinal")
        return(sprintf(
          "predictors[[%s]] <- ordered(predictors[[%s]], levels = %s)",
          qstr(nm), qstr(nm), qvec(exp_level_order(p, nm))))
      if (type == "categorical")
        return(sprintf("predictors[[%s]] <- factor(predictors[[%s]], levels = %s)",
                       qstr(nm), qstr(nm), qvec(exp_category_levels(p, nm))))
      sprintf("predictors[[%s]] <- as.numeric(predictors[[%s]])",
              qstr(nm), qstr(nm))
    }, character(1))
  })
  output$exp_formula_controls <- renderUI({
    if (is.null(input$exp_predictors) && is.null(sim_predictors_val()) &&
        is.null(restored_project_resources()$predictors))
      return(p(class = "text-muted small",
               "Upload predictor metadata to choose effects and interactions."))
    p <- tryCatch(exp_predictors_raw(), error = function(e) NULL)
    if (is.null(p)) return(p(class = "text-danger small",
                             "The predictor metadata could not be read."))
    is_cj <- identical(input$model_type, "btl")
    vars <- setdiff(names(p), c("item", "object", "threshold",
                                "threshold_number"))
    # `threshold` is supplied by the estimator even for item-level metadata.
    # It is selected by default only when the nominated responses appear
    # polytomous.
    its <- input$item_cols
    poly <- FALSE
    if (length(its)) {
      d <- raw_data()[, its, drop = FALSE]
      poly <- any(vapply(d, function(v)
        suppressWarnings(max(as.numeric(as.character(v)), na.rm = TRUE)) > 1,
        logical(1)))
    }
    choices <- c(vars, if (!is_cj) c("threshold", "threshold_number"))
    selected <- c(vars, if (!is_cj && poly) "threshold")
    tagList(
      selectizeInput("exp_main", info_label("Main effects",
        if (is_cj) "Predictors of object location." else
          paste("Predictors of item or threshold location. Threshold denotes",
                "a categorical within-item threshold effect; threshold_number",
                "fits one linear effect per threshold number.")),
        choices = choices, selected = selected, multiple = TRUE,
        options = list(placeholder = "choose at least one")),
      selectizeInput("exp_interactions", info_label("Two-way interactions",
        paste("Allows one selected predictor's effect to differ over another.",
              "Both main effects remain in the model.")),
        choices = character(0), multiple = TRUE,
        options = list(placeholder = "none")))
  })
  exp_interaction_map <- reactiveVal(list())
  observeEvent(input$exp_main, {
    z <- input$exp_main %||% character(0)
    spec <- .app_explanatory_interactions(z)
    old_map <- exp_interaction_map()
    requested <- unique(c(
      input$exp_interactions %||% character(0),
      restored_project_settings()$exp_interactions %||% character(0),
      sim_interactions_val()))
    wanted <- list()
    for (value in requested) {
      pair <- old_map[[value]] %||% spec$map[[value]]
      # Schema-2 projects and bundled simulations created before opaque
      # interaction ids stored the readable colon-joined label. Retain that
      # unambiguous legacy case; new projects never rely on it.
      if (is.null(pair)) {
        hit <- which(vapply(spec$map, function(p)
          identical(paste(p, collapse = ":"), value), logical(1)))
        if (length(hit) == 1L) pair <- spec$map[[hit]]
      }
      if (!is.null(pair)) wanted[[length(wanted) + 1L]] <- pair
    }
    selected <- names(spec$map)[vapply(spec$map, function(pair)
      any(vapply(wanted, identical, logical(1), y = pair)), logical(1))]
    exp_interaction_map(spec$map)
    updateSelectizeInput(session, "exp_interactions",
                         choices = spec$choices, selected = selected,
                         server = TRUE)
  }, ignoreNULL = FALSE)

  exp_formula <- reactive({
    main <- input$exp_main
    ids <- input$exp_interactions %||% character(0)
    map <- exp_interaction_map()
    if (length(ids) && any(!ids %in% names(map)))
      stop("the explanatory interaction selection is stale; choose it again")
    .app_explanatory_formula(main, unname(map[ids]))
  })

  # paired-comparison anchors: a two-column CSV (object, location) that places
  # this calibration on an existing scale, parsed defensively like anchors_in()
  bt_anchors_in <- reactive({
    a <- if (!is.null(input$bt_anchor_file))
      tryCatch(read.csv(input$bt_anchor_file$datapath,
                        check.names = FALSE, stringsAsFactors = FALSE),
               error = function(e)
                 stop("could not read the paired-comparison anchor CSV: ",
                      conditionMessage(e), call. = FALSE))
    else restored_project_resources()$bt_anchors
    if (is.null(a)) return(NULL)
    if (anyDuplicated(names(a)))
      stop("the paired-comparison anchor CSV has duplicate column names",
           call. = FALSE)
    if (!all(c("object", "location") %in% names(a)))
      stop("the paired-comparison anchor CSV needs columns object, location",
           call. = FALSE)
    if (!nrow(a))
      stop("the paired-comparison anchor CSV contains no anchor rows",
           call. = FALSE)
    missing_object <- is.na(a$object)
    a$object <- trimws(as.character(a$object))
    a$object[missing_object] <- NA_character_
    if (anyNA(a$object) || any(!nzchar(a$object)))
      stop("every paired-comparison anchor row needs a non-blank object name",
           call. = FALSE)
    a
  })

  observeEvent(raw_data(), {
    df <- raw_data(); nm <- names(df)
    guess_id <- nm[grepl("^id$|_id$|^person", tolower(nm))][1]
    guess_fac <- intersect(nm, c("group", "sex", "gender", "site", "country", "age_group"))
    updateSelectInput(session, "id_col", choices = c(NONE_CH, nm),
                      selected = if (!is.na(guess_id)) guess_id else NONE)
    updateSelectizeInput(session, "factor_cols", choices = nm, selected = guess_fac)
    updateSelectizeInput(session, "item_cols", choices = nm,
                         selected = setdiff(nm, c(guess_id, guess_fac)))
    # MFRM role guesses
    g_per <- nm[grepl("person|candidate|student|^id$|_id$", tolower(nm))][1]
    g_itm <- nm[grepl("item|task|criterion|question", tolower(nm))][1]
    g_sco <- nm[grepl("score|rating|grade|mark", tolower(nm))][1]
    g_fac <- setdiff(nm[grepl("rater|judge|marker|occasion|time", tolower(nm))],
                     c(g_per, g_itm, g_sco))
    updateSelectInput(session, "lp_person", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_per)) g_per else NONE)
    updateSelectInput(session, "lp_item", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_itm)) g_itm else NONE)
    updateSelectInput(session, "lp_score", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_sco)) g_sco else NONE)
    updateSelectizeInput(session, "lp_facets", choices = nm, selected = g_fac)
    # wide layout: item columns = the remaining columns (like the Rasch
    # item_cols guess), excluding the guessed person and facet columns
    updateSelectizeInput(session, "lp_items_wide", choices = nm,
                         selected = setdiff(nm, c(g_per, g_fac)))
    # frames layout guesses
    g_grp <- nm[grepl("group|year|grade|cohort|class$", tolower(nm))][1]
    updateSelectInput(session, "ef_id", choices = c(NONE_CH, nm),
                      selected = if (!is.na(guess_id)) guess_id else NONE)
    updateSelectInput(session, "ef_group", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_grp)) g_grp else NONE)
    updateSelectizeInput(session, "ef_items", choices = nm,
                         selected = setdiff(nm, c(guess_id, g_grp)))
    # paired-comparison guesses
    g_a <- nm[grepl("^a$|object_a|left|first|option_a", tolower(nm))][1]
    g_b <- nm[grepl("^b$|object_b|right|second|option_b", tolower(nm))][1]
    g_w <- nm[grepl("win|preferred|chosen|better", tolower(nm))][1]
    g_j <- nm[grepl("judge|rater|marker", tolower(nm))][1]
    g_c <- nm[grepl("^count$|^n$|freq", tolower(nm))][1]
    g_o <- nm[grepl("^t$|order|seq|trial|round|^time$", tolower(nm))][1]
    # judge factors are columns constant within each judge (and not another
    # role) -- exactly the shape of a panel or rater-group variable, so they
    # are offered for judge-group DIF by default
    g_f <- character(0)
    if (!is.na(g_j)) {
      cand <- setdiff(nm, stats::na.omit(c(g_a, g_b, g_w, g_j, g_o, g_c)))
      if (length(cand))
        g_f <- cand[vapply(cand, function(cn) all(tapply(as.character(df[[cn]]),
                     df[[g_j]], function(v) length(unique(v)) == 1L)), TRUE)]
    }
    updateSelectInput(session, "bt_a", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_a)) g_a else NONE)
    updateSelectInput(session, "bt_b", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_b)) g_b else NONE)
    updateSelectInput(session, "bt_win", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_w)) g_w else NONE)
    # polytomous margin column, if present (a winner and a response are mutually
    # exclusive; auto-detecting the response keeps polytomous data one-click)
    g_resp <- nm[grepl("^response$|^grade$|^rating$", tolower(nm))][1]
    updateSelectizeInput(session, "bt_response", choices = c("", nm),
                         selected = if (!is.na(g_resp)) g_resp else "")
    updateSelectizeInput(session, "bt_margin", choices = c("", nm),
                         selected = "")
    updateSelectInput(session, "bt_judge", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_j)) g_j else NONE)
    # empty choice = the "none" placeholder (clearable)
    updateSelectizeInput(session, "bt_order", choices = c("", nm),
                         selected = if (!is.na(g_o)) g_o else "")
    saved_jf <- names(restored_project_resources()$btl_dif_factors$factors)
    updateSelectizeInput(session, "bt_jfactors", choices = unique(c(nm, saved_jf)),
                         selected = g_f)
    # judge-panel column for the paired-comparison frames extension: the same
    # judge-constant candidates as the judge factors above, one of them chosen
    updateSelectInput(session, "btlef_panel",
                      choices = if (length(g_f)) g_f else character(0),
                      selected = if (length(g_f)) g_f[1] else character(0))
    updateSelectInput(session, "bt_count", choices = c(NONE_CH, nm),
                      selected = if (!is.na(g_c)) g_c else NONE)
  })

  read_frame_map <- function(upload, id_col, units, label) {
    mp <- tryCatch(
      read.csv(upload$datapath, stringsAsFactors = FALSE,
               check.names = FALSE),
      error = function(e) e)
    if (inherits(mp, "error"))
      validate(need(FALSE, paste0("Could not read the ", label,
                                 "-set CSV: ", conditionMessage(mp))))
    validate(need(!anyDuplicated(names(mp)),
                  paste0("The ", label, "-set CSV has duplicate column names.")))
    validate(need(all(c(id_col, "set") %in% names(mp)),
                  paste0("The ", label, "-set CSV needs columns ",
                         id_col, " and set.")))
    ids <- as.character(mp[[id_col]])
    sets <- as.character(mp$set)
    validate(need(!anyNA(ids) && all(nzchar(trimws(ids))),
                  paste0("Every row of the ", label,
                         "-set CSV needs a non-blank ", id_col, ".")))
    validate(need(!anyNA(sets) && all(nzchar(trimws(sets))),
                  paste0("Every row of the ", label,
                         "-set CSV needs a non-blank set.")))
    ids <- trimws(ids)
    sets <- trimws(sets)
    validate(need(!anyDuplicated(ids),
                  paste0("The ", label, "-set CSV assigns ", id_col,
                         " values more than once: ",
                         paste(unique(ids[duplicated(ids)]), collapse = ", "), ".")))
    extra <- setdiff(ids, units)
    validate(need(!length(extra),
                  paste0("The ", label, "-set CSV includes unknown ", id_col,
                         " values: ", paste(utils::head(extra, 5L),
                                            collapse = ", "), ".")))
    stats::setNames(sets, ids)
  }

  # item -> set map: uploaded CSV wins; otherwise infer from the item-name
  # prefix (the part before trailing digits/separators)
  ef_setmap <- reactive({
    its <- if (length(input$ef_items)) input$ef_items else
      setdiff(names(raw_data()),
              c(if (!is.null(input$ef_id) && input$ef_id != NONE) input$ef_id,
                if (!is.null(input$ef_group) && input$ef_group != NONE) input$ef_group))
    if (!is.null(input$ef_sets)) {
      out <- read_frame_map(input$ef_sets, "item", its, "item")
      miss <- setdiff(its, names(out))
      if (length(miss)) out[miss] <- "(rest)"
      return(out[its])
    }
    saved <- restored_project_resources()$ef_setmap
    if (!is.null(saved)) {
      miss <- setdiff(its, names(saved))
      validate(need(!length(miss),
                    paste0("The saved item-set map does not cover: ",
                           paste(utils::head(miss, 5L), collapse = ", "), ".")))
      return(saved[its])
    }
    if (isTRUE(input$ef_prefix)) {
      pref <- sub("[_. -]*[0-9]+$", "", its)
      pref[pref == ""] <- "(rest)"
      return(setNames(pref, its))
    }
    setNames(rep("all", length(its)), its)
  })

  # the interacting facet (shown in interactive facet mode) is chosen from
  # the nominated facet columns; a single facet is preselected automatically
  observeEvent(input$lp_facets, {
    fs <- input$lp_facets
    sel <- if (!is.null(input$lp_interaction) && input$lp_interaction %in% fs)
      input$lp_interaction else if (length(fs)) fs[1] else character(0)
    updateSelectInput(session, "lp_interaction",
                      choices = if (length(fs)) fs else character(0),
                      selected = sel)
  }, ignoreNULL = FALSE)

  # keep the wide-mode item choices free of the chosen person / facet columns
  observeEvent(c(input$lp_person, input$lp_facets), {
    df <- raw_data(); nm <- names(df)
    taken <- c(if (!is.null(input$lp_person) && input$lp_person != NONE)
                 input$lp_person,
               input$lp_facets)
    sel <- setdiff(if (length(input$lp_items_wide)) input$lp_items_wide else nm,
                   taken)
    updateSelectizeInput(session, "lp_items_wide",
                         choices = setdiff(nm, taken), selected = sel)
  }, ignoreInit = TRUE)

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
                    mfrm = "Multiple Ratings (MFRM)",
                    efrm = "Extended Frames (EFRM)",
                    btl = "Comparative Judgement")
  .demo_chip_labels <- c(dich = "Multiple choice", pcm = "Polytomous (PCM)",
                         rsm = "Rating scale", mfrm = "Multiple Ratings",
                         efrm = "Extended Frames", btl = "Comparative Judgement")
  output$data_main <- renderUI({
    if (is.null(sim_data()) &&
        identical(input$demo_choice %||% "none", "none") &&
        is.null(input$file)) {
      div(class = "mx-auto", style = "max-width: 760px; margin-top: 8vh;",
        card(
          card_body(class = "empty-state",
            bs_icon("clipboard-data", size = "3rem", class = "text-primary mb-2"),
            h2("Welcome to rasch"),
            p(class = "lead mb-4",
              "Fit and evaluate models within Rasch Measurement Theory."),
            # Reopening a saved analysis is a way to START one, so it belongs
            # beside uploading rather than only under More: the file input it
            # triggers stays on the Export panel, next to the download that
            # writes the file.
            div(class = "d-flex justify-content-center gap-2 flex-wrap mb-4",
              tags$button(class = "btn btn-outline-primary btn-lg", type = "button",
                          onclick = "document.getElementById('file').click();",
                          bs_icon("upload"), " Upload data"),
              tags$button(class = "btn btn-outline-primary btn-lg", type = "button",
                          onclick = "document.getElementById('project_file').click();",
                          bs_icon("folder2-open"), " Open saved analysis")),
            p(class = "text-muted small mb-2", "Or start from a specific example:"),
            div(class = "d-flex flex-wrap gap-1 justify-content-center",
              lapply(names(.demo_chip_labels), function(k)
                actionButton(paste0("demo_chip_", k), .demo_chip_labels[[k]],
                             class = "btn-outline-secondary btn-sm"))))))
    } else {
      tagList(
        uiOutput("data_strip"),
        rcode_details("data_strip"),
        card(card_header(span("Data preview",
                              info_icon(app_help("preview"),
                                        "About this table"))),
             card_body(uiOutput("data_info"), DTOutput("preview"),
                       rcode_details("preview"),
                       padding = 12, fillable = FALSE)),
        accordion(id = "rcode_acc", open = FALSE, class = "mt-3",
          accordion_panel("R code for this analysis", icon = bs_icon("code-slash"),
            p(class = "text-muted small mb-2",
              "The exact rasch call reproducing the current run; updates on every estimation."),
            verbatimTextOutput("rcode_fit"))))
    }
  })
  lapply(c("dich", "pcm", "rsm", "mfrm", "efrm", "btl"), function(k)
    observeEvent(input[[paste0("demo_chip_", k)]],
      updateSelectInput(session, "demo_choice", selected = k)))
  output$data_strip <- renderUI({
    df <- raw_data()
    vals <- as.matrix(df)
    miss <- 100 * mean(is.na(vals) | trimws(vals) == "", na.rm = FALSE)
    layout_column_wrap(width = "160px", fill = FALSE, class = "mb-3",
      value_box("Rows", format(nrow(df), big.mark = ","),
                showcase = glyph("grid"),
                showcase_layout = "left center", theme = "primary"),
      value_box("Columns", ncol(df), showcase = glyph("columns"),
                showcase_layout = "left center", theme = "warning"),
      value_box("Missing", sprintf("%.1f%%", miss),
                showcase = glyph("missing"),
                showcase_layout = "left center",
                theme = if (miss > 20) "warning" else "secondary"))
  })
  output$data_info <- renderUI({
    df <- raw_data()
    p(class = "text-muted",
      sprintf("%d rows x %d columns.%s Nominate the column roles in the sidebar, then press Estimate. Missing responses may be left blank or coded as -1; any negative score is read as missing.",
              nrow(df), ncol(df),
              if (nrow(df) > 200) " First 200 rows shown in the preview." else ""))
  })
  output$preview <- renderDT({
    d <- head(raw_data(), 200)
    round_preview(
      datatable(d, rownames = FALSE, style = "bootstrap5",
                class = "table-sm compact hover order-column",
                options = list(pageLength = 25, scrollX = TRUE, dom = "tip")), d)
  })
  output$rcode_fit <- renderText({
    validate(need(!is.null(current_rcode()),
                  "Run an analysis to see the reproducible R code."))
    current_rcode()
  })

  # ----------------------------------------------------------------- fit --
  # Ordered structural changes applied to the base fit.  Each entry carries
  # the resulting fitted object, so every downstream output reads one active
  # analysis while the app can still show, undo and eventually serialise the
  # complete transformation history.
  analysis_steps <- reactiveVal(list())
  btl_analysis_steps <- reactiveVal(list())
  restoring_project <- reactiveVal(FALSE)
  active_step <- function() {
    h <- analysis_steps()
    if (length(h)) h[[length(h)]] else NULL
  }
  active_step_type <- function() {
    s <- active_step()
    if (is.null(s)) NULL else s$type
  }
  push_analysis_step <- function(type, label, fitted, details = list(),
                                 code = NULL) {
    h <- analysis_steps()
    h[[length(h) + 1L]] <- list(
      type = type, label = label, fit = fitted, details = details,
      code = code,
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )
    analysis_steps(h)
    advance_analysis_context()
    invisible(fitted)
  }
  clear_analysis_steps <- function() analysis_steps(list())
  active_btl_step <- function() {
    h <- btl_analysis_steps()
    if (length(h)) h[[length(h)]] else NULL
  }
  push_btl_analysis_step <- function(type, label, fitted, details = list(),
                                     code = NULL) {
    h <- btl_analysis_steps()
    h[[length(h) + 1L]] <- list(
      type = type, label = label, fit = fitted, details = details, code = code,
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
    btl_analysis_steps(h)
    advance_analysis_context()
    invisible(fitted)
  }
  clear_btl_analysis_steps <- function() btl_analysis_steps(list())
  undo_analysis_step <- function(notify = TRUE) {
    if (!is.null(tryCatch(btl_fit(), error = function(e) NULL)) &&
        length(btl_analysis_steps())) {
      h <- btl_analysis_steps(); removed <- h[[length(h)]]
      btl_analysis_steps(h[-length(h)])
      if (identical(removed$type, "btl_frames"))
        try(btlef_res(NULL), silent = TRUE)
      advance_analysis_context()
      if (isTRUE(notify))
        showNotification(paste("Undid", removed$label), type = "message",
                         duration = 5)
      return(invisible(TRUE))
    }
    h <- analysis_steps()
    if (!length(h)) return(invisible(FALSE))
    removed <- h[[length(h)]]
    h <- h[-length(h)]
    analysis_steps(h)
    advance_analysis_context()
    if (identical(removed$type, "dif_auto")) resolve_res(NULL)
    if (isTRUE(notify))
      showNotification(paste("Undid", removed$label), type = "message",
                       duration = 5)
    invisible(TRUE)
  }
  # the exact rasch call reproducing the current run (built alongside the fit)
  rcode_str <- reactiveVal(NULL)
  # The main disclosure is a complete analysis, not merely the first fit.
  # Structural changes are appended in the order applied, and the untouched
  # fit is retained for the before/after tables.
  current_rcode <- reactive({
    base <- rcode_str()
    if (is.null(base)) return(NULL)
    is_bt <- !is.null(tryCatch(btl_fit(), error = function(e) NULL))
    h <- if (is_bt) btl_analysis_steps() else analysis_steps()
    step_code <- vapply(h, function(s) s$code %||% "", character(1))
    step_code <- step_code[nzchar(step_code)]
    if (!length(step_code)) return(base)
    original <- if (is_bt) "original_bt <- bt" else "original_fit <- fit"
    paste(c(base, "", original, "", step_code), collapse = "\n")
  })
  output$has_override <- reactive(length(analysis_steps()) > 0L ||
                                    length(btl_analysis_steps()) > 0L)
  outputOptions(output, "has_override", suspendWhenHidden = FALSE)
  # Local undo buttons appear only when that page created the most recent step.
  output$has_override_subtest <- reactive(
    identical(active_step_type(), "superitem"))
  output$has_override_dif <- reactive(
    isTRUE(active_step_type() %in% c("dif_split", "dif_auto")))
  outputOptions(output, "has_override_subtest", suspendWhenHidden = FALSE)
  outputOptions(output, "has_override_dif", suspendWhenHidden = FALSE)
  output$override_status <- renderUI({
    is_bt <- !is.null(tryCatch(btl_fit(), error = function(e) NULL))
    h <- if (is_bt) btl_analysis_steps() else analysis_steps()
    if (!length(h)) return(NULL)
    div(class = "analysis-pipeline mt-2 mb-1",
      div(class = "small fw-semibold mb-1", "Active analysis"),
      div(class = "analysis-pipeline-steps",
        span(class = "badge text-bg-light border",
             if (is_bt) "Equal-unit fit" else "Original fit"),
        lapply(h, function(s) tagList(
          bs_icon("chevron-right", class = "analysis-pipeline-arrow"),
          span(class = "badge text-bg-warning", s$label)))))
  })
  # Reset returns every page to the base analysis. Local buttons undo only the
  # most recent page-specific change, preserving any earlier transformations.
  reset_to_original <- function() {
    clear_analysis_steps(); clear_btl_analysis_steps(); resolve_res(NULL)
    try(btlef_res(NULL), silent = TRUE)
    advance_analysis_context()
    showNotification("Reset to the original data; showing the base analysis.",
                     type = "message", duration = 5)
  }
  observeEvent(input$undo_override, undo_analysis_step())
  observeEvent(input$reset_override, reset_to_original())
  observeEvent(input$reset_split,   undo_analysis_step())
  observeEvent(input$reset_subtest, undo_analysis_step())

  # estimation runs as a side-effecting observer that stores the completed
  # fit in fit_val; analysis() is a pure accessor, so every reader keeps
  # working unchanged
  fit_val <- reactiveVal(NULL)
  analysis <- reactive(fit_val())
  btl_fit <- reactiveVal(NULL)
  clear_btl_fit_results <- function() {
    bdif_res(NULL)
    bdif_meta(NULL)
    btlef_res(NULL)
    invisible(NULL)
  }
  # shared estimation-control resolution (used by every model branch and
  # the reproducible-code footers)
  est_opts <- reactive(list(maxit = max(5, input$maxit %||% 60),
                            tol = max(1e-12, input$tol %||% 1e-8)))
  data_source_code <- reactive({
    if (!is.null(sim_data()) && nzchar(sim_code_val() %||% ""))
      return(paste0("dat <- ", sim_code_val()))
    if (!is.null(sim_data()) && nzchar(restored_project_name() %||% ""))
      return(paste0("project <- readRDS(", qstr(restored_project_name()),
                    ")\ndat <- project$data"))
    if (identical(input$demo_choice %||% "none", "pcm"))
      return(paste(
        "dat <- simulate_rasch(",
        "  n_persons = 600,",
        "  n_items = 12,",
        "  model = \"PCM\",",
        "  n_categories = 4,",
        "  difficulty = c(-1.5, 1.5),",
        "  disordered = \"I04\",",
        "  dependence = list(pairs = list(c(\"I10\", \"I11\")), strength = 1.3),",
        "  dif = list(items = \"I08\", uniform = 0.8),",
        "  n_groups = 3,",
        "  seed = 17",
        ")", sep = "\n"))
    if (identical(input$demo_choice %||% "none", "btl"))
      return(paste(
        "dat <- simulate_btl(",
        "  n_objects = 8,",
        "  n_judges = 48,",
        "  reps_per_pair = 84,",
        "  erratic_judges = 2 / 48,",
        "  dependence = list(exposure = 0.7, carry_over = 0),",
        "  seed = 2",
        ")",
        "judge_number <- as.integer(sub(\"^J\", \"\", dat$judge))",
        "dat$panel <- factor(ifelse(judge_number %% 2L,",
        "                           \"panel A\", \"panel B\"))",
        "dat$experience <- factor(ifelse(judge_number <= 24,",
        "                                \"experienced\", \"novice\"))",
        sep = "\n"))
    if (!identical(input$demo_choice %||% "none", "none"))
      return(sprintf("dat <- rasch:::.app_example_data(%s)",
                     qstr(input$demo_choice)))
    req(input$file)
    sprintf('dat <- read.csv(%s%s, check.names = FALSE, stringsAsFactors = FALSE)',
            qstr(input$file$name),
            if (tolower(tools::file_ext(input$file$name)) %in%
                c("tsv", "txt")) ', sep = "\\t"' else "")
  })

  # Store a completed fit from either the ordinary in-process route or the
  # cancellable EFRM worker. Keeping this in one function prevents the two
  # execution routes from diverging in navigation, notes, reproducible code,
  # convergence handling, and downstream active-state propagation.
  complete_fit <- function(fit, code_call, code_notes, src_line,
                           simulation_stamp = NULL, source_state = NULL) {
    if (inherits(fit, "error")) {
      showNotification(paste("Analysis failed:", conditionMessage(fit)),
                       type = "error", duration = 10)
      return(invisible(NULL))
    }
    if (!is.null(source_state))
      attr(fit, "rasch_app_source") <- source_state
    if (!is.null(code_call))
      rcode_str(paste(c("library(rasch)", "", src_line,
                        if (length(code_notes)) c("", code_notes), "",
                        code_call), collapse = "\n"))
    if (length(fit$notes))
      # each note is an independent statement, and the notification renders
      # as HTML, where a newline is only whitespace: joining on one runs
      # separate notes into a single unreadable sentence
      showNotification(
        if (length(fit$notes) == 1L) fit$notes
        else tags$ul(class = "mb-0 ps-3",
                     lapply(fit$notes, function(n) tags$li(n))),
        type = "message", duration = 8)
    conv <- if (!is.null(fit$est)) fit$est$converged else fit$converged
    if (!isTRUE(conv)) {
      msg <- if (inherits(fit, "rasch_efrm") &&
                 isTRUE(fit$est$stage1_converged))
        paste("The conditional calibration converged, but the item-set link",
              "did not. Check the common-person design, targeting and score",
              "range within each linked set.")
      else paste("Estimation did not converge; consider raising the maximum",
                 "iterations or loosening the convergence criterion.")
      showNotification(msg, type = "warning", duration = 10)
    }
    # Replace the preceding analysis only after the new fit has succeeded.
    # A failed or cancelled background fit therefore leaves its fitted object
    # and complete structural history available.
    advance_analysis_context()
    clear_analysis_steps()
    clear_btl_analysis_steps()
    clear_btl_fit_results()
    # A weighted person table belongs to one completed calibration.  Clear it
    # here rather than relying only on the ordinary-fit observer: when a Rasch
    # analysis is replaced by Comparative Judgement, fit() becomes unavailable
    # before that observer can run and the old table would otherwise survive.
    person_weight_state(NULL)
    fitted_sim_gen(simulation_stamp)
    if (inherits(fit, "rasch_btl")) {
      btl_fit(fit); fit_val(NULL)
      try(nav_select("nav", "p_summary", session = session), silent = TRUE)
      return(invisible(NULL))
    }
    btl_fit(NULL)
    try(nav_select("nav", "p_summary", session = session), silent = TRUE)
    fit_val(NULL)
    fit_val(fit)
    invisible(fit)
  }

  # EFRM uncertainty runs in a separate R process. The main Shiny process can
  # therefore continue to receive input, update a real progress bar, and stop
  # the worker immediately when the analyst presses Cancel.
  efrm_job <- reactiveVal(NULL)
  output$efrm_job_controls <- renderUI({
    if (is.null(efrm_job())) return(NULL)
    actionButton("cancel_efrm", "Cancel EFRM estimation",
                 icon = bs_icon("x-circle"),
                 class = "btn-outline-danger w-100 mt-2")
  })
  outputOptions(output, "efrm_job_controls", suspendWhenHidden = FALSE)

  close_efrm_job <- function(st, remove_files = TRUE) {
    if (!is.null(st$progress)) try(st$progress$close(), silent = TRUE)
    if (isTRUE(remove_files))
      unlink(c(st$progress_file, st$log_file), force = TRUE)
    invisible(NULL)
  }

  stop_efrm_process <- function(process) {
    stopped <- tryCatch({
      process$kill_tree()
      TRUE
    }, error = function(e) FALSE)
    if (!stopped && process$is_alive()) try(process$kill(), silent = TRUE)
    invisible(NULL)
  }

  background_job_is_current <- function(st) {
    if (!is.list(st) ||
        !identical(st$context, isolate(analysis_context()))) return(FALSE)
    current_data <- tryCatch(isolate(raw_data()), error = function(e) NULL)
    if (!identical(st$data, current_data)) return(FALSE)
    if (isTRUE(st$check_btl_fit) &&
        !identical(st$base_fit, isolate(btl_fit()))) return(FALSE)
    TRUE
  }

  cancel_efrm_job <- function() {
    st <- isolate(efrm_job())
    if (is.null(st)) return(invisible(FALSE))
    if (st$process$is_alive()) stop_efrm_process(st$process)
    close_efrm_job(st)
    efrm_job(NULL)
    invisible(TRUE)
  }

  observeEvent(input$cancel_efrm, {
    if (!cancel_efrm_job()) return(invisible(NULL))
    showNotification("EFRM estimation cancelled.", type = "message",
                     duration = 5)
  })

  session$onSessionEnded(function() {
    st <- isolate(efrm_job())
    if (!is.null(st) && st$process$is_alive())
      stop_efrm_process(st$process)
    if (!is.null(st)) close_efrm_job(st)
  })

  observe({
    st <- efrm_job()
    if (is.null(st)) return()
    invalidateLater(250, session)
    if (file.exists(st$progress_file)) {
      z <- tryCatch(strsplit(readLines(st$progress_file, warn = FALSE)[1],
                             "\t", fixed = TRUE)[[1]],
                    error = function(e) character(0))
      if (length(z) == 3L) {
        current <- suppressWarnings(as.numeric(z[2])); total <-
          suppressWarnings(as.numeric(z[3])); stage <- z[1]
        fraction <- if (is.finite(current) && is.finite(total) && total > 0)
          pmin(pmax(current / total, 0), 1) else 0
        value <- switch(stage,
          "conditional calibration" = 0.02 + 0.08 * fraction,
          "linking bootstrap" = if (identical(st$se_method, "bootstrap"))
            0.10 + 0.40 * fraction else 0.10 + 0.84 * fraction,
          "full person bootstrap" = 0.50 + 0.44 * fraction,
          "finalising" = 0.99, 0.02)
        detail <- if (stage %in% c("linking bootstrap",
                                   "full person bootstrap") && total > 0)
          sprintf("%s: %d of %d", stage, round(current), round(total)) else stage
        st$progress$set(value = value, detail = detail)
      }
    }
    if (st$process$is_alive()) return()

    result <- tryCatch(st$process$get_result(), error = function(e) e)
    close_efrm_job(st)
    efrm_job(NULL)
    if (!background_job_is_current(st)) {
      showNotification(paste(
        "The completed EFRM result was not used because the active data or",
        "analysis changed while it was running."),
        type = "warning", duration = 8)
      return()
    }
    if (inherits(result, "error")) {
      complete_fit(result, st$code_call, character(0), st$src_line,
                   st$simulation_stamp)
      return()
    }
    if (length(result$warnings))
      showNotification(paste(result$warnings, collapse = "\n"),
                       type = "warning", duration = 10)
    complete_fit(result$value, st$code_call, character(0), st$src_line,
                 st$simulation_stamp, st$source_state)
  })

  observeEvent(input$run, {
    if (!is.null(boot_job())) {
      showNotification("A fit bootstrap is already running. Cancel it before starting another analysis.",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    if (!is.null(efrm_job())) {
      showNotification("An EFRM estimation is already running. Cancel it before starting another analysis.",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    if (!is.null(btlef_job())) {
      showNotification("A frame estimation is already running. Cancel it before starting another analysis.",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    df <- raw_data()
    # Freeze the source of the calibration at launch. The data and controls
    # in the sidebar may subsequently change while the fitted object remains
    # available; a project file and any downstream frame model must still use
    # the source that produced that object.
    run_source <- list(
      data = as.data.frame(df, check.names = FALSE),
      settings = .collect_fit_settings(input),
      resources = list(anchors = NULL, bt_anchors = NULL, key = NULL,
                       predictors = NULL, ef_setmap = NULL),
      simulation = list(
        data = sim_data(), truth = sim_truth_val(), code = sim_code_val(),
        predictors = sim_predictors_val(),
        interactions = sim_interactions_val(), generation = sim_gen(),
        fitted_generation = if (!is.null(sim_data())) sim_gen() else NULL))
    # automatic class intervals pass NULL; rasch() resolves the rule and
    # reports the value used in fit$n_groups
    ng <- if (isTRUE(input$ng_auto %||% TRUE)) NULL else input$ng
    # reproducible-code pieces (spliced into the branch-specific call below)
    src_line <- data_source_code()
    code_args_common <- c(
      if (!is.null(ng)) paste0("n_groups = ", ng))
    eo <- est_opts()
    code_est <- c(paste0("maxit = ", eo$maxit),
                  paste0("tol = ", format(eo$tol)))
    code_call <- NULL
    code_notes <- character(0)

    if (identical(input$model_type, "efrm")) {
      sm <- ef_setmap()
      run_source$resources["ef_setmap"] <- list(sm)
      reps_raw <- suppressWarnings(as.numeric(input$ef_reps))
      if (length(reps_raw) != 1L || !is.finite(reps_raw) ||
          reps_raw != floor(reps_raw) || reps_raw < 50L ||
          reps_raw > .Machine$integer.max) {
        showNotification(paste("EFRM bootstrap replicates must be one whole",
                               "number of at least 50"),
                         type = "warning", duration = 7)
        return(invisible(NULL))
      }
      reps <- as.integer(reps_raw)
      workers_raw <- suppressWarnings(as.numeric(input$ef_workers))
      if (length(workers_raw) != 1L || !is.finite(workers_raw) ||
          workers_raw != floor(workers_raw) ||
          !workers_raw %in% .efrm_worker_values) {
        showNotification("EFRM workers must be one of the available whole-number choices",
                         type = "warning", duration = 7)
        return(invisible(NULL))
      }
      workers <- as.integer(workers_raw)
      seed_raw <- suppressWarnings(as.numeric(input$ef_seed))
      if (length(seed_raw) != 1L || !is.finite(seed_raw) || seed_raw < 0 ||
          seed_raw != floor(seed_raw) || seed_raw > .Machine$integer.max) {
        showNotification(paste("EFRM seed must be one non-negative whole",
                               "number within the integer range"),
                         type = "warning", duration = 7)
        return(invisible(NULL))
      }
      seed <- as.integer(seed_raw)
      group_arg <- if (is.null(input$ef_group) || input$ef_group == NONE)
        rep("(all)", nrow(df)) else input$ef_group
      id_arg <- if (!is.null(input$ef_id) && input$ef_id != NONE)
        input$ef_id else NULL
      code_call <- paste0("fit <- rasch_efrm(dat,\n  ", paste(c(
        paste0("item_sets = ", paste(deparse(sm), collapse = "\n    ")),
        if (is.null(input$ef_group) || input$ef_group == NONE)
          'groups = rep("(all)", nrow(dat))'
        else paste0("groups = ", qstr(input$ef_group)),
        if (!is.null(id_arg)) paste0("id = ", qstr(input$ef_id)),
        paste0("items = ", qvec(names(sm))), code_args_common, code_est,
        paste0("se_method = ", qstr(input$ef_se %||% "hybrid")),
        if (!is.null(reps)) paste0("boot_reps = ", reps),
        paste0("workers = ", workers), paste0("seed = ", seed)),
        collapse = ",\n  "), ")")
      fit_args <- list(
        data = df, item_sets = sm, groups = group_arg, id = id_arg,
        items = names(sm), n_groups = ng,
        maxit = eo$maxit, tol = eo$tol,
        se_method = input$ef_se %||% "hybrid", boot_reps = reps,
        workers = workers, seed = seed)
      progress_file <- tempfile("rasch-efrm-", fileext = ".progress")
      log_file <- tempfile("rasch-efrm-", fileext = ".log")
      progress_bar <- shiny::Progress$new(session, min = 0, max = 1)
      progress_bar$set(message = "Estimating Extended Frames",
                       detail = "conditional calibration", value = 0.02)
      worker <- tryCatch(callr::r_bg(
        function(args, progress_path, source_dir) {
          if (!is.null(source_dir)) {
            if (!requireNamespace("pkgload", quietly = TRUE))
              stop("pkgload is needed for a source-tree background analysis")
            pkgload::load_all(dirname(source_dir), quiet = TRUE)
          }
          progress_fun <- function(stage, current, total) {
            writeLines(paste(stage, current, total, sep = "\t"),
                       progress_path)
          }
          args$progress <- progress_fun
          warnings <- character(0)
          value <- withCallingHandlers(
            do.call(if (exists("rasch_efrm", inherits = TRUE))
              get("rasch_efrm", inherits = TRUE) else
                getExportedValue("rasch", "rasch_efrm"), args),
            warning = function(w) {
              warnings <<- c(warnings, conditionMessage(w))
              invokeRestart("muffleWarning")
            })
          list(value = value, warnings = unique(warnings))
        },
        args = list(args = fit_args, progress_path = progress_file,
                    source_dir = .rasch_source_dir),
        libpath = .libPaths(), stdout = log_file, stderr = log_file,
        supervise = TRUE), error = function(e) e)
      if (inherits(worker, "error")) {
        progress_bar$close()
        unlink(c(progress_file, log_file), force = TRUE)
        complete_fit(worker, code_call, character(0), src_line,
                     if (!is.null(sim_data())) sim_gen() else NULL,
                     run_source)
        return(invisible(NULL))
      }
      efrm_job(list(
        process = worker, progress = progress_bar,
        progress_file = progress_file, log_file = log_file,
        code_call = code_call, src_line = src_line,
        se_method = input$ef_se %||% "hybrid",
        context = isolate(analysis_context()), data = df,
        simulation_stamp = if (!is.null(sim_data())) sim_gen() else NULL,
        source_state = run_source))
      return(invisible(NULL))
    }

    withProgress(message = "Estimating (pairwise conditional ML)…", value = 0.3, {
      fit <- tryCatch({
        if (identical(input$model_type, "btl")) {
          # a polytomous response column overrides the winner column (and the
          # ties rule: polytomous ties belong in a middle category); otherwise a
          # margin column combines with the winner into the polytomous response
          # (winner values matching neither object = ties = middle category)
          bt_graded <- !is.null(input$bt_response) && nzchar(input$bt_response)
          bt_marg <- !bt_graded && !is.null(input$bt_margin) &&
            nzchar(input$bt_margin)
          bt_thr <- input$bt_thr %||% "free"
          bt_exp <- identical(input$rasch_calibration %||% "free",
                              "explanatory")
          if (bt_exp && !is.null(input$bt_anchor_file))
            stop("paired-comparison anchors cannot be combined with an ",
                 "explanatory object calibration; remove the anchor file ",
                 "or choose free calibration")
          # the judgment-order column enables the within-judge dependence
          # analysis (exposure and carry-over); it only exists with a judge
          bt_ord <- if (!is.null(input$bt_order) && nzchar(input$bt_order) &&
                        !is.null(input$bt_judge) && input$bt_judge != NONE)
            input$bt_order else NULL
          bt_count <- if (!is.null(input$bt_count) && input$bt_count != NONE)
            input$bt_count else NULL
          # first-position advantage (object A is the first-presented of the
          # pair) and equating anchors (a named location per object) both feed
          # btl() directly; anchors come from a two-column CSV
          bt_pos <- isTRUE(input$bt_position)
          bt_anc_df <- if (bt_exp) NULL else bt_anchors_in()
          bt_anchor_vec <- if (!is.null(bt_anc_df))
            setNames(bt_anc_df$location, bt_anc_df$object) else NULL
          run_source$resources["bt_anchors"] <- list(bt_anc_df)
          if (!is.null(bt_anc_df)) {
            code_notes <- c(code_notes,
              if (!is.null(input$bt_anchor_file)) paste0(
                "bt_anchors <- read.csv(", qstr(input$bt_anchor_file$name),
                ", check.names = FALSE, stringsAsFactors = FALSE)") else
                  "bt_anchors <- project$resources$bt_anchors",
              paste0("bt_anchors$object <- ifelse(is.na(bt_anchors$object), ",
                     "NA_character_, trimws(as.character(bt_anchors$object)))"),
              "bt_anchor_values <- with(bt_anchors, setNames(location, object))")
          }
          if (any(c(input$bt_a, input$bt_b) == NONE) ||
              (!bt_graded && identical(input$bt_win, NONE)))
            stop("nominate the object A, object B, and winner (or polytomous response) columns")
          if (bt_exp) {
            ep <- exp_predictors_in()
            ef <- exp_formula()
            run_source$resources["predictors"] <- list(ep)
            code_notes <- c(
              if (is.null(sim_predictors_val()) && !is.null(input$exp_predictors))
                paste0("predictors <- read.csv(", qstr(input$exp_predictors$name),
                       ", check.names = FALSE, stringsAsFactors = FALSE)")
              else if (!is.null(restored_project_resources()$predictors))
                "predictors <- project$resources$predictors"
              else 'predictors <- attr(dat, "predictors")',
              exp_predictor_code(),
              paste0("explanatory_formula <- ",
                     paste(deparse(ef), collapse = " ")))
          }
          code_call <- paste0("bt <- ",
            if (bt_exp) "btl_explanatory" else "btl", "(dat,\n  ", paste(c(
            if (bt_exp) "predictors = predictors",
            if (bt_exp) "formula = explanatory_formula",
            paste0("object_a = ", qstr(input$bt_a)),
            paste0("object_b = ", qstr(input$bt_b)),
            if (bt_graded) paste0("response = ", qstr(input$bt_response))
            else paste0("winner = ", qstr(input$bt_win)),
            if (bt_marg) paste0("margin = ", qstr(input$bt_margin)),
            if (!is.null(input$bt_judge) && input$bt_judge != NONE)
              paste0("judge = ", qstr(input$bt_judge)),
            if (!is.null(bt_ord)) paste0("order = ", qstr(bt_ord)),
            if (bt_pos) "position = TRUE",
            if (!is.null(bt_anchor_vec)) "anchors = bt_anchor_values",
            if (!is.null(bt_count)) paste0("count = ", qstr(bt_count)),
            if (!bt_graded && !bt_marg)
              paste0("ties = ", qstr(input$bt_ties %||% "drop")),
            if ((bt_graded || bt_marg) && identical(bt_thr, "pc"))
              'thresholds = "pc"',
            code_est), collapse = ",\n  "), ")")
          # one shared argument list; the entry path (polytomous response,
          # winner + margin, winner only) contributes its own arguments
          bt_args <- c(
            list(df),
            if (bt_exp) list(predictors = ep, formula = ef) else list(),
            list(object_a = input$bt_a, object_b = input$bt_b,
                 judge = if (!is.null(input$bt_judge) &&
                             input$bt_judge != NONE) input$bt_judge else NULL,
                 order = bt_ord,
                 position = bt_pos,
                 count = bt_count,
                 maxit = eo$maxit, tol = eo$tol),
            if (!bt_exp) list(anchors = bt_anchor_vec) else list(),
            if (bt_graded)
              list(response = input$bt_response, thresholds = bt_thr)
            else if (bt_marg)
              list(winner = input$bt_win, margin = input$bt_margin,
                   thresholds = bt_thr)
            else
              list(winner = input$bt_win, ties = input$bt_ties %||% "drop"))
          do.call(if (bt_exp) btl_explanatory else btl, bt_args)
        } else if (identical(input$model_type, "mfrm")) {
          # the interaction is passed only in interactive facet mode, and
          # only when the interacting facet is one of the chosen facets
          lp_int <- if (identical(input$lp_structure %||% "additive",
                                  "interactive") &&
                        !is.null(input$lp_interaction) &&
                        input$lp_interaction %in% input$lp_facets)
            input$lp_interaction else NULL
          if (identical(input$lp_layout %||% "wide", "wide")) {
            if (identical(input$lp_person, NONE) || !length(input$lp_facets) ||
                !length(input$lp_items_wide))
              stop("nominate the person column, the item columns, and at least one facet column")
            code_call <- paste0("fit <- rasch_mfrm(dat,\n  ", paste(c(
              paste0("person = ", qstr(input$lp_person)),
              paste0("facets = ", qvec(input$lp_facets)),
              paste0("items = ", qvec(input$lp_items_wide)),
              code_args_common,
              if (!is.null(lp_int)) paste0("interaction = ", qstr(lp_int)),
              code_est), collapse = ",\n  "), ")")
            rasch_mfrm(df, person = input$lp_person,
                       facets = input$lp_facets, items = input$lp_items_wide,
                       n_groups = ng, interaction = lp_int,
                       maxit = eo$maxit, tol = eo$tol)
          } else {
            if (any(c(input$lp_person, input$lp_item, input$lp_score) == NONE) ||
                !length(input$lp_facets))
              stop("nominate the person, item, score, and at least one facet column")
            code_call <- paste0("fit <- rasch_mfrm(dat,\n  ", paste(c(
              paste0("person = ", qstr(input$lp_person)),
              paste0("item = ", qstr(input$lp_item)),
              paste0("score = ", qstr(input$lp_score)),
              paste0("facets = ", qvec(input$lp_facets)),
              code_args_common,
              if (!is.null(lp_int)) paste0("interaction = ", qstr(lp_int)),
              code_est), collapse = ",\n  "), ")")
            rasch_mfrm(df, person = input$lp_person, item = input$lp_item,
                       score = input$lp_score, facets = input$lp_facets,
                       n_groups = ng, interaction = lp_int,
                       maxit = eo$maxit, tol = eo$tol)
          }
        } else {
          idc <- if (!is.null(input$id_col) && input$id_col != NONE) input$id_col else NULL
          fac <- if (length(input$factor_cols)) input$factor_cols else NULL
          its <- if (length(input$item_cols)) input$item_cols else NULL
          # multiple-choice key: an uploaded CSV, the copy embedded in a
          # reopened project, or the bundled key for the dichotomous example
          mc_key <- key_in()
          exp_on <- identical(input$rasch_calibration %||% "free",
                              "explanatory")
          if (exp_on) {
            ep <- exp_predictors_in()
            ef <- exp_formula()
            run_source$resources["predictors"] <- list(ep)
            level <- input$exp_level %||% "item"
            code_notes <- c(
              if (is.null(sim_predictors_val()) && !is.null(input$exp_predictors))
                paste0("predictors <- read.csv(",
                       qstr(input$exp_predictors$name),
                       ", check.names = FALSE, stringsAsFactors = FALSE)")
              else if (!is.null(restored_project_resources()$predictors))
                "predictors <- project$resources$predictors"
              else 'predictors <- attr(dat, "predictors")',
              exp_predictor_code(),
              paste0("explanatory_formula <- ",
                     paste(deparse(ef), collapse = " ")))
            code_call <- paste0("fit <- rasch_explanatory(dat,\n  ", paste(c(
              "predictors = predictors",
              "formula = explanatory_formula",
              paste0("level = ", qstr(level)),
              if (!is.null(idc)) paste0("id = ", qstr(idc)),
              if (!is.null(fac)) paste0("factors = ", qvec(fac)),
              if (!is.null(its)) paste0("items = ", qvec(its)),
              code_args_common,
              if (!is.null(mc_key) && !is.null(input$key_file))
                paste0("key = read.csv(", qstr(input$key_file$name), ")")
              else if (!is.null(restored_project_resources()$key))
                "key = project$resources$key"
              else if (!is.null(mc_key))
                'key = setNames(rep("A", 15), sprintf("I%02d", 1:15))',
              code_est), collapse = ",\n  "), ")")
            rasch_explanatory(
              df, predictors = ep, formula = ef, level = level,
              id = idc, factors = fac, items = its,
              n_groups = ng, key = mc_key,
              maxit = eo$maxit, tol = eo$tol)
          } else {
          # anchors match by item name; a supplied map must apply in full
          anc <- anchors_in()
          run_source$resources["anchors"] <- list(anc)
          if (!is.null(anc)) {
            cand <- if (!is.null(its)) its else setdiff(names(df), c(idc, fac))
            present <- as.character(anc$item) %in% cand
            if (!all(present))
              stop("anchor item(s) are not in the selected analysis: ",
                   paste(unique(as.character(anc$item[!present])),
                         collapse = ", "))
            # location anchoring: collapse to one mean-location anchor per
            # item; average anchoring of the set adds the flag pcml() reads
            anchor_type <- input$anchor_type %||% "individual"
            if (!is.null(anc) && anchor_type %in% c("average", "set_average")) {
              mu <- tapply(anc$tau, as.character(anc$item), mean)
              anc <- data.frame(item = names(mu), k = NA, tau = as.numeric(mu))
              if (anchor_type == "set_average") anc$average <- TRUE
            }
          }
          # threshold structure: partial credit (item-specific) or rating
          # (common); dichotomous items are the one-threshold special case
          rsm_on <- identical(input$thr_structure %||% "pcm", "rsm")
          # principal-components (Andrich) threshold estimation; partial
          # credit route only, and not combinable with anchors
          pcc <- if (!rsm_on && identical(input$thr_mode, "pc"))
            as.integer(input$pc_rank %||% "4") else NULL
          if (!is.null(pcc) && !is.null(anc)) {
            stop("anchors cannot be combined with principal-components ",
                 "threshold estimation; remove the anchor file or choose ",
                 "ordinary partial-credit estimation")
          }
          code_notes <- if (!is.null(anc)) c(
            if (!is.null(input$anchor_file))
              paste0("anchors <- read.csv(", qstr(input$anchor_file$name),
                     ", check.names = FALSE, stringsAsFactors = FALSE)")
            else "anchors <- project$resources$anchors",
            paste0("anchors$item <- ifelse(is.na(anchors$item), NA_character_, ",
                   "trimws(as.character(anchors$item)))"),
            paste0("anchors <- anchors[as.character(anchors$item) %in% ",
                   qvec(cand), ", , drop = FALSE]"),
            if (anchor_type %in% c("average", "set_average")) c(
              "anchor_mean <- tapply(anchors$tau, as.character(anchors$item), mean)",
              "anchors <- data.frame(item = names(anchor_mean), k = NA, tau = as.numeric(anchor_mean))"),
            if (anchor_type == "set_average")
              "anchors$average <- TRUE")
          else character(0)
          code_call <- paste0("fit <- rasch(dat,\n  ", paste(c(
            paste0("model = ", qstr(if (rsm_on) "RSM" else "PCM")),
            if (!is.null(idc)) paste0("id = ", qstr(idc)),
            if (!is.null(fac)) paste0("factors = ", qvec(fac)),
            if (!is.null(its)) paste0("items = ", qvec(its)),
            code_args_common,
            if (!is.null(anc)) "anchors = anchors",
            if (!is.null(mc_key) && !is.null(input$key_file))
              paste0("key = read.csv(", qstr(input$key_file$name), ")")
            else if (!is.null(restored_project_resources()$key))
              "key = project$resources$key"
            else if (!is.null(mc_key))
              'key = setNames(rep("A", 15), sprintf("I%02d", 1:15))',
            if (!is.null(pcc)) paste0("pc_components = ", pcc),
            code_est), collapse = ",\n  "), ")")
          rasch(df, model = if (rsm_on) "RSM" else "PCM",
                id = idc, factors = fac, items = its,
                n_groups = ng, anchors = anc,
                key = mc_key, pc_components = pcc,
                maxit = eo$maxit, tol = eo$tol)
          }
        }
      }, error = function(e) e)
    })
    # Scoring metadata is part of either ordinary Rasch calibration. Store
    # the parsed object, not a temporary upload path, even when no key was
    # supplied (NULL is an exact part of the run specification).
    if (identical(input$model_type, "rasch"))
      run_source$resources["key"] <- list(
        tryCatch(key_in(), error = function(e) NULL))
    complete_fit(fit, code_call, code_notes, src_line,
                 if (!is.null(sim_data())) sim_gen() else NULL,
                 run_source)
  })
  fit <- reactive({
    s <- active_step()
    f <- if (is.null(s)) analysis() else s$fit
    req(f); f
  })
  # The same active-state resolution without req(): NULL before any run,
  # for UI that must render quietly in that state (navbar chips, report)
  fit_or_null <- function() {
    s <- active_step()
    if (!is.null(s)) return(s$fit)
    tryCatch(analysis(), error = function(e) NULL)
  }

  # A stable automatic logit range shared by the plot controls. Extreme-score
  # person estimates are omitted because their finite WLE values can otherwise
  # flatten the informative part of a plot. Thresholds are always retained.
  fitted_scale_range <- function(fallback = c(-5, 5), pad = .5) {
    f <- fit()
    th <- f$person$theta
    if (!is.null(f$person$extreme))
      th <- th[!f$person$extreme]
    th <- th[is.finite(th)]
    if (length(th) > 20L)
      th <- as.numeric(stats::quantile(th, c(.01, .99), na.rm = TRUE,
                                       names = FALSE))
    tau <- f$thresholds$tau
    z <- c(th, tau[is.finite(tau)])
    if (length(z) < 2L) return(fallback)
    r <- range(z) + c(-pad, pad)
    r <- c(floor(r[1] * 2) / 2, ceiling(r[2] * 2) / 2)
    if (!all(is.finite(r)) || r[1] >= r[2]) return(fallback)
    if (diff(r) < 4) {
      m <- mean(r)
      r <- c(floor((m - 2) * 2) / 2, ceiling((m + 2) * 2) / 2)
    }
    r
  }
  resolve_axis_range <- function(id, standard, wide, automatic) {
    mode <- input[[paste0(id, "_mode")]] %||% "auto"
    if (identical(mode, "auto")) return(automatic())
    if (identical(mode, "wide")) return(wide)
    if (!identical(mode, "custom")) return(standard)
    r <- suppressWarnings(as.numeric(c(input[[paste0(id, "_min")]],
                                       input[[paste0(id, "_max")]])))
    if (length(r) != 2L || any(!is.finite(r)) || r[1] >= r[2]) standard else r
  }

  output$has_mc <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    !is.null(f) && !is.null(f$mc)
  })
  outputOptions(output, "has_mc", suspendWhenHidden = FALSE)

  output$sel_item_title <- renderUI(span(class = "fw-semibold",
    tryCatch(sel_item(), error = function(e) "Selected item")))

  # only offer the pages that apply to the current analysis: Facets needs a
  # many-facet fit, Frames an extended-frames fit, and Guessing a
  # dichotomous one. Summary and Items serve Rasch and paired-comparison
  # (BTL) fits alike (each page shows the matching variant); Persons needs
  # person estimates (Rasch) or a judge column (BTL). Everything else stays.
  # Every nav_panel and nav_menu carries an explicit value, and visibility is
  # driven by those values through the rasch-nav-vis handler (shiny::hideTab,
  # which backs bslib::nav_hide, cannot reach entries inside a nav_menu).
  observe({
    f <- tryCatch(fit(), error = function(e) NULL)
    bf <- btl_fit()
    active_bf <- if (is.null(bf)) NULL else {
      s <- active_btl_step()
      if (is.null(s)) bf else s$fit
    }
    show <- function(value, on)
      session$sendCustomMessage("rasch-nav-vis",
                                list(value = value, show = isTRUE(on)))
    show("p_facets", inherits(f, "rasch_mfrm"))
    # the Frames page also carries the paired-comparison (BTL) extended frame
    # of reference branch, so it shows for either fit family
    show("p_frames", inherits(f, "rasch_efrm") || !is.null(bf))
    show("p_guess", !is.null(f) && !inherits(f, "rasch_mfrm") &&
           !inherits(f, "rasch_efrm") && max(f$m) == 1L)
    rasch_on <- !is.null(f)
    btl_on <- !is.null(bf)
    ordinary_rasch_on <- rasch_on &&
      !inherits(f, c("rasch_mfrm", "rasch_efrm"))
    equating_rasch_on <- rasch_on && .app_equating_fit(f)
    # judge-factor DIF applies to a paired-comparison fit once judge
    # factors are nominated in the Data roles
    btl_dif_on <- btl_on &&
      !inherits(active_bf, c("rasch_btl_efrm", "rasch_btl_explanatory")) &&
      length(input$bt_jfactors) > 0
    rasch_dif_factors <- if (rasch_on)
      setdiff(names(f$factors), f$frame_group %||% character(0)) else
        character(0)
    show("p_summary", rasch_on || btl_on)
    show("p_items", rasch_on || btl_on)
    show("p_explanatory", inherits(f, "rasch_explanatory") ||
           inherits(bf, "rasch_btl_explanatory"))
    show("p_persons", rasch_on || (btl_on && !is.null(bf$judges)))
    # Targeting is model aware. Common-item equating is defined only for an
    # ordinary person-by-item calibration; paired comparisons have their own
    # common-object method.
    show("p_targeting", rasch_on || btl_on)
    show("p_equating", equating_rasch_on || btl_on)
    # the Trait tab now carries paired-comparison dimensionality too
    # (transitivity loops + the residual bimension swirl)
    show("p_dim", ordinary_rasch_on || btl_on)
    # Local dependence: the Rasch Q3 suite, or the within-judge dependence
    # analysis of a paired-comparison fit
    show("p_ld", rasch_on || btl_on)
    # DIF needs at least one person factor in the fit (Rasch), or a
    # paired-comparison fit with judge factors nominated
    show("p_dif", length(rasch_dif_factors) > 0L || btl_dif_on)
    # Compare applies to either family: Rasch fits and paired-comparison
    # (BTL) fits can each be kept and compared among themselves (their
    # likelihoods are over different data, so the table groups kept fits by
    # family around whichever one is chosen as the reference)
    show("p_compare", rasch_on || btl_on)
    # menu headers hide too when everything inside them is hidden
    show("menu_independence", rasch_on || btl_on)
    show("menu_invariance", rasch_on || btl_dif_on)
    # More always shows: Simulate is the entry point when no data is loaded
    show("menu_more", TRUE)
    show("p_simulate", TRUE)
  })

  # ------------------------------------------------ UI visibility flags --
  # "nothing to show" areas hide instead of rendering empty; each flag pairs
  # with a conditionalPanel in the UI (same pattern as has_mc)
  output$has_interaction <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    inherits(f, "rasch_mfrm") && !is.null(f$interaction)
  })
  outputOptions(output, "has_interaction", suspendWhenHidden = FALSE)
  output$is_btl <- reactive(!is.null(btl_fit()))
  outputOptions(output, "is_btl", suspendWhenHidden = FALSE)
  output$btl_history_model <- reactive({
    b <- tryCatch(bfit(), error = function(e) NULL)
    is.data.frame(b$dependence) && any(b$dependence$effect != "position")
  })
  outputOptions(output, "btl_history_model", suspendWhenHidden = FALSE)
  output$has_expl_relaxations <- reactive({
    f <- tryCatch({
      b <- btl_fit()
      if (!is.null(b)) {
        s <- active_btl_step(); if (is.null(s)) b else s$fit
      } else fit()
    }, error = function(e) NULL)
    inherits(f, c("rasch_explanatory", "rasch_btl_explanatory")) &&
      nrow(f$explanatory$relaxations) > 0L
  })
  outputOptions(output, "has_expl_relaxations", suspendWhenHidden = FALSE)
  output$dep_magnitude_available <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    !is.null(f) && !inherits(f, "rasch_efrm")
  })
  outputOptions(output, "dep_magnitude_available", suspendWhenHidden = FALSE)
  output$subtest_available <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    !is.null(f) && !inherits(f, c("rasch_mfrm", "rasch_efrm"))
  })
  outputOptions(output, "subtest_available", suspendWhenHidden = FALSE)
  output$has_structural_items <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    inherits(f, c("rasch_mfrm", "rasch_efrm"))
  })
  outputOptions(output, "has_structural_items", suspendWhenHidden = FALSE)
  output$anchor_download_available <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    !is.null(f) && !inherits(f, c("rasch_mfrm", "rasch_efrm"))
  })
  outputOptions(output, "anchor_download_available", suspendWhenHidden = FALSE)
  output$spread_available <- reactive({
    f <- tryCatch(fit(), error = function(e) NULL)
    !is.null(f) && length(f$subtest_map) > 0L
  })
  outputOptions(output, "spread_available", suspendWhenHidden = FALSE)
  # a paired-comparison fit with a judge column: the Persons page shows the
  # judge fit table (and is offered in the navbar) only then
  output$has_judges <- reactive({
    b <- btl_fit()
    !is.null(b) && !is.null(b$judges)
  })
  outputOptions(output, "has_judges", suspendWhenHidden = FALSE)
  # empty-state flags for the run-on-demand cards: before a run the card
  # shows only its controls and one muted line (no table, plot, or download)
  output$has_dep <- reactive(!is.null(dep_res()))
  outputOptions(output, "has_dep", suspendWhenHidden = FALSE)
  output$has_spread <- reactive(!is.null(spread_res()))
  outputOptions(output, "has_spread", suspendWhenHidden = FALSE)
  output$has_contr <- reactive(!is.null(contr_res()))
  outputOptions(output, "has_contr", suspendWhenHidden = FALSE)
  output$has_dm <- reactive(!is.null(dm_res()))
  outputOptions(output, "has_dm", suspendWhenHidden = FALSE)
  output$has_rescore <- reactive(!is.null(rescore_res()))
  outputOptions(output, "has_rescore", suspendWhenHidden = FALSE)
  output$has_cmp <- reactive(length(kept_fits()) >= 2)
  outputOptions(output, "has_cmp", suspendWhenHidden = FALSE)
  # equating results exist once a reference is available (upload or kept fit)
  output$has_eq <- reactive({
    if (identical(input$eq_source, "kept"))
      !is.null(input$eq_kept) && nzchar(input$eq_kept) &&
        input$eq_kept %in% names(kept_fits())
    else !is.null(input$eq_file) ||
      !is.null(restored_project_resources()[["eq_reference"]])
  })
  outputOptions(output, "has_eq", suspendWhenHidden = FALSE)

  observeEvent(fit(), {
    its <- fit()$items$item
    updateSelectizeInput(session, "icc_compare_items", choices = its,
                         selected = character(0), server = TRUE)
    updateSelectizeInput(session, "subtest_items", choices = its, selected = character(0))
    fac <- setdiff(names(fit()$factors), fit()$frame_group %||% character(0))
    # MFRM DIF is pooled to underlying items by default. Planned contrasts
    # must therefore offer those item names, not the internal item-by-facet
    # response cells displayed in fit$items.
    dif_items <- .app_dif_item_choices(fit())
    updateSelectizeInput(session, "pc_items", choices = dif_items,
                         selected = character(0))
    fs <- if (inherits(fit(), "rasch_mfrm")) fit()$facet_spec else character(0)
    updateSelectizeInput(session, "facet_sel", choices = fs,
                         selected = if (length(fs)) fs[1] else character(0))
    fi <- if (inherits(fit(), "rasch_efrm"))
      unique(fit()$virtual_map$item) else character(0)
    updateSelectizeInput(session, "frame_item", choices = fi,
                         selected = if (length(fi)) fi[1] else character(0))
    f <- fit()
    pf <- names(f$factors) %||% character(0)
    person_choices <- c("One panel" = "")
    if (inherits(f, "rasch_efrm"))
      person_choices <- c(person_choices,
                          "Fitted person groups" = "groups")
    if (length(pf))
      person_choices <- c(person_choices,
        stats::setNames(pf, paste("Person factor:", pf)))
    person_choices <- person_choices[!duplicated(unname(person_choices))]
    old_person <- input$wright_person_panels %||% ""
    updateSelectInput(session, "wright_person_panels",
      choices = person_choices,
      selected = if (old_person %in% unname(person_choices)) old_person else "")

    item_choices <- c("One panel" = "")
    if (.wrightmap_item_panels && inherits(f, "rasch_efrm")) {
      item_choices <- c(item_choices, "Item sets" = "sets",
                        "Person groups" = "groups",
                        "Frames (set x group)" = "sets_groups")
    } else if (.wrightmap_item_panels && inherits(f, "rasch_mfrm") &&
               length(f$facet_spec)) {
      item_choices <- c(item_choices,
        stats::setNames(f$facet_spec, paste("Facet:", f$facet_spec)))
    }
    if (.wrightmap_item_panels)
      item_choices <- c(item_choices, "Uploaded item map" = "uploaded")
    old_item <- input$wright_item_panels %||% ""
    updateSelectInput(session, "wright_item_panels",
      choices = item_choices,
      selected = if (old_item %in% unname(item_choices)) old_item else "")
    updateSelectizeInput(session, "dim_pos", choices = its, selected = character(0))
    updateSelectizeInput(session, "dim_neg", choices = its, selected = character(0))
    # residual principal components available for the t-test default split:
    # 1..min(10, n_items - 1)
    updateSelectInput(session, "pca_component",
                      choices = seq_len(max(1L, min(10L, length(its) - 1L))),
                      selected = 1)
    updateSelectInput(session, "dep_item", choices = its,
                      selected = its[min(2L, length(its))])
    updateSelectInput(session, "ind_item", choices = its, selected = its[1])
    updateSelectizeInput(session, "guess_anchors", choices = its,
                         selected = character(0))
    # explorer class intervals start at the fit's own rule
    updateSelectInput(session, "ex_ng", selected = fit()$n_groups)
    # Results computed on request belong to the fit they came from. A project
    # restore is the exception: those results were serialised with that exact
    # fit and are reinstated later in the same flush cycle.
    if (!isTRUE(restoring_project())) {
      # Allow a 10% refit-failure margin above the Holm resolution threshold.
      updateNumericInput(session, "guess_boot_reps",
                         value = max(999L, ceiling(40 * length(its) / 0.9)))
      lr_res(NULL); dep_res(NULL); spread_res(NULL); dm_res(NULL)
      guess_res(NULL); contr_res(NULL); rescore_res(NULL)
      person_weight_state(NULL)
      restored_dimensionality(NULL); restored_subtest(NULL)
      dim_computed(NULL)
      restored_invariance(NULL)
      # An automatic resolution sets the override fit itself, so its trace
      # must survive its own refit; a fresh run or another override clears it.
      if (!identical(active_step_type(), "dif_auto")) resolve_res(NULL)
      # Manual dimensionality subsets name items of the previous fit.
      dim_subsets(NULL)
    }
  })

  observeEvent(input$make_subtest, {
    req(length(input$subtest_items) >= 2)
    res <- tryCatch(combine_items(fit(), list(input$subtest_items)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Subtest failed:", conditionMessage(res)), type = "error")
    } else {
      push_analysis_step(
        "superitem",
        paste("Superitem:", paste(input$subtest_items, collapse = " + ")),
        res, details = list(items = input$subtest_items),
        code = sprintf("fit <- combine_items(fit, list(%s))",
                       qvec(input$subtest_items))
      )
      showNotification("Re-analysed with the subtest in place. Use Reset to original data (or run again) to return to the base fit.",
                       type = "message", duration = 8)
    }
  })
  output$subtest_status <- renderUI({
    s <- active_step()
    if (is.null(s) || !identical(s$type, "superitem")) return(NULL)
    p(class = "text-success small mt-2", paste("Active:", s$label))
  })

  observeEvent(input$make_split, {
    f <- fit()
    it <- dif_sel_item()
    vars <- dif_sel_vars()
    selected <- dif_res()$summary[dif_sel_row(), , drop = FALSE]
    if (isTRUE(selected$nonuniform_DIF)) {
      showNotification(
        paste("This item has non-uniform DIF. A location split cannot model",
              "a change in discrimination; review or remove the item instead."),
        type = "warning", duration = 10)
      return()
    }
    req(it %in% f$items$item,
        !is.null(f$factors), length(vars) >= 1, all(vars %in% names(f$factors)))
    # one factor -> split by the factor name; an interaction row -> split by
    # the factor-combination cells. This is a uniform-DIF location resolution.
    by <- if (length(vars) == 1L) vars
          else app_factor_cells(f$factors[vars], sep = ":")
    res <- tryCatch(split_items(f, it, by = by), error = function(e) e)
    lab <- paste(vars, collapse = ":")
    if (inherits(res, "error")) {
      showNotification(paste("Split failed:", conditionMessage(res)), type = "error")
    } else {
      push_analysis_step(
        "dif_split", sprintf("DIF split: %s by %s", it, lab), res,
        details = list(item = it, factors = vars),
        code = if (length(vars) == 1L)
          sprintf("fit <- split_items(fit, %s, by = %s)",
                  qstr(it), qstr(vars))
        else paste0(
          "dif_group <- rasch:::.factor_cells(fit$factors[", qvec(vars),
          "], sep = \":\")\n",
          "fit <- split_items(fit, ", qstr(it), ", by = dif_group)")
      )
      showNotification(
        sprintf("Re-analysed with %s split by %s. Use Reset to original data (or run again) to return to the base fit.",
                it, lab),
        type = "message", duration = 8)
    }
  })

  # automatic iterative DIF resolution: split the largest-effect item, refit,
  # repeat until no item flags DIF (or the anchor set would fall too low). The
  # resolved fit becomes the active override; the trace is kept for the panel.
  resolve_res <- reactiveVal(NULL)
  observeEvent(input$resolve_all, {
    run_factors <- names(fit()$factors)
    run_alpha <- dif_alpha()
    run_adjust <- "holm"
    run_effects <- input$dif_effects %||% "main"
    rr <- tryCatch(
      resolve_dif(fit(), factors = run_factors, alpha = run_alpha,
                  p_adjust = run_adjust, effects = run_effects),
      error = function(e) e)
    if (inherits(rr, "error")) {
      showNotification(paste("Automatic resolution failed:", conditionMessage(rr)),
                       type = "error")
    } else {
      rr$run_factors <- run_factors
      rr$run_alpha <- run_alpha
      rr$run_p_adjust <- run_adjust
      rr$run_effects <- run_effects
      resolve_res(rr)
      push_analysis_step(
        "dif_auto", sprintf("Automatic DIF: %d split(s)", rr$n_splits),
        rr$fit, details = list(n_splits = rr$n_splits, stopped = rr$stopped),
        code = paste0(
          "dif_resolution <- resolve_dif(fit, factors = ", qvec(run_factors),
          ", alpha = ", run_alpha, ", p_adjust = ", qstr(run_adjust),
          ", effects = ", qstr(run_effects), ")\n",
          "fit <- dif_resolution$fit")
      )
      showNotification(
        sprintf("Automatic DIF resolution: %d split(s); %s. Use Reset to original data (or run again) to return to the base fit.",
                rr$n_splits, rr$stopped),
        type = "message", duration = 8)
    }
  })
  output$has_resolve <- reactive(!is.null(resolve_res()))
  outputOptions(output, "has_resolve", suspendWhenHidden = FALSE)
  output$resolve_summary <- renderUI({
    rr <- resolve_res(); req(!is.null(rr))
    p(class = "text-muted small mb-2",
      sprintf("%d split(s); %s; %d item(s) still flag DIF.",
              rr$n_splits, rr$stopped, rr$n_remaining_dif))
  })
  output$resolve_tbl <- DT::renderDT({
    rr <- resolve_res(); req(!is.null(rr))
    d <- rr$splits
    if (nrow(d)) { d$eta2 <- round(d$eta2, 3); d$magnitude <- round(d$magnitude, 3) }
    num_dt(d)
  })
  output$resolve_tbl_csv <- downloadHandler(
    filename = function() "dif_resolution.csv",
    content = function(file) {
      rr <- resolve_res(); req(!is.null(rr))
      write_csv_plain(rr$splits, file)
    })

  sel_item <- reactive({
    f <- fit()
    i <- input$items_tbl_rows_selected
    f$items$item[.app_selected_row(i, nrow(f$items))]
  })
  # A frame fit displays virtual item names. Resolve them through the fitted
  # map; delimiters are not parsed because item and factor labels may contain
  # colons or other punctuation.
  sel_source_item <- reactive({
    f <- fit(); it <- sel_item()
    src <- names(f$set_of)
    if (is.null(src) || it %in% src) return(it)
    vm <- f$virtual_map
    if (is.null(vm)) return(it)
    hit <- match(it, vm$vkey)
    if (is.na(hit)) it else vm$item[hit]
  })

  # Dropping the selected item and refitting is a step in the analysis, not
  # housekeeping: for frame models, removing an item changes both the
  # within-frame calibration and the person-side link between sets. The step
  # is recorded like any other, so it can be undone and is saved with the
  # project.
  output$drop_item_ui <- renderUI({
    f <- tryCatch(fit(), error = function(e) NULL)
    if (is.null(f) || inherits(f, "rasch_mfrm")) return(NULL)
    drop_btn <- span(class = "d-inline-flex align-items-center",
      actionButton("drop_item", "Drop item and refit", icon = bs_icon("trash"),
                   class = "btn-outline-secondary btn-xs rasch-control-button"),
      info_icon(paste("Removes the selected item and refits the same model.",
                      "The active item and person estimates are replaced by",
                      "the refitted estimates."), "About dropping an item"))
    # for a frame model the milder remedy belongs beside the blunt one: the
    # item keeps measuring inside each frame and only stops linking them
    if (!inherits(f, "rasch_efrm")) return(drop_btn)
    tagList(drop_btn,
      span(class = "d-inline-flex align-items-center",
        actionButton("resolve_item", "Resolve by frame and refit",
                     icon = bs_icon("arrows-expand"),
                     class = "btn-outline-secondary btn-xs rasch-control-button"),
        info_icon(paste("Estimates a separate location in each frame while",
                        "retaining the item in person measurement. The refit",
                        "uses the remaining links between frames. This addresses",
                        "a location difference, not a discrimination difference."),
                  "About resolving an item")))
  })

  observeEvent(input$resolve_item, {
    f <- fit()
    it <- sel_source_item()
    req(length(it) == 1L)
    if (identical(inv_se(), "bootstrap")) {
      inv <- efrm_invariance()
      if (!is.character(inv)) {
        loc <- any(inv$locations$item == it & inv$locations$flagged %in% TRUE)
        dsc <- any(inv$discrimination$item == it &
                     inv$discrimination$flagged %in% TRUE)
        if (dsc && !loc) {
          showNotification(
            paste("This item is flagged for discrimination, not location.",
                  "Resolving its location by frame would not fit that departure;",
                  "review or remove the item instead."),
            type = "warning", duration = 10)
          return()
        }
      }
    }
    res <- tryCatch(resolve_frames(f, it), error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Could not resolve", it, "--",
                             conditionMessage(res)),
                       type = "error", duration = 10)
    } else {
      push_analysis_step(
        "resolve_item", sprintf("Resolved by frame: %s", it), res,
        details = list(item = it),
        code = sprintf("fit <- resolve_frames(fit, %s)", qstr(it)))
      showNotification(
        sprintf(paste("%s now has a location in each frame and no longer",
                      "links them. Undo to restore the shared location."), it),
        type = "message", duration = 8)
    }
  })

  observeEvent(input$drop_item, {
    f <- fit()
    it <- sel_source_item()
    req(length(it) == 1L)
    res <- tryCatch(drop_items(f, it), error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Could not drop", it, "--", conditionMessage(res)),
                       type = "error", duration = 10)
    } else {
      push_analysis_step(
        "drop_item", sprintf("Dropped item: %s", it), res,
        details = list(item = it),
        code = sprintf("fit <- drop_items(fit, %s)", qstr(it)))
      showNotification(
        sprintf("Re-analysed without %s. Undo to restore it, or use Reset to original data.", it),
        type = "message", duration = 8)
    }
  })

  icc_items <- reactive({
    primary <- sel_item()
    comparison <- setdiff(input$icc_compare_items %||% character(0), primary)
    unique(c(primary, head(comparison, 7L)))
  })

  # ------------------------------------------------------- plot plumbing --
  # per-output "R code" disclosure: `code` is a function returning the exact
  # rasch call reproducing the output (it may read reactives, so the snippet
  # follows the current selections). Rendering is never suspended: the text
  # must be ready when the collapsed <details> footer is opened.
  register_code <- function(id, code) {
    cid <- paste0(id, "_code")
    output[[cid]] <- renderText(code())
    outputOptions(output, cid, suspendWhenHidden = FALSE)
  }
  register_code("resolve_tbl", function() "dif_resolution$splits")
  register_code("sim_recovery", function() {
    obj <- if (!is.null(tryCatch(btl_fit(), error = function(e) NULL))) "bt" else "fit"
    paste0("recovery <- sim_recovery(", obj, ", dat)\n",
           "recovery$summary\nplot_recovery(recovery)")
  })
  register_code("data_strip", function() paste(
    data_source_code(),
    "list(",
    "  rows = nrow(dat), columns = ncol(dat),",
    "  missing_percent = 100 * mean(is.na(as.matrix(dat)) |",
    "    trimws(as.matrix(dat)) == \"\", na.rm = FALSE)",
    ")", sep = "\n"))
  register_code("preview", function()
    paste(data_source_code(), "head(dat, 200)", sep = "\n"))
  register_code("sim_truth", function()
    paste(data_source_code(), 'attr(dat, "truth")', sep = "\n"))
  register_code("sim_preview", function()
    paste(data_source_code(), "head(dat, 12)", sep = "\n"))
  # `px` optionally sets a reactive on-screen height (a function returning
  # pixels), for plots whose natural size grows with the data (e.g. a matrix
  # heatmap that should track its item count and the table beside it)
  register_plot <- function(id, fun, w = 9, h = 6, code = NULL, px = NULL) {
    output[[id]] <- if (is.null(px)) renderPlot(soft(fun()), res = 96)
                    else renderPlot(soft(fun()), res = 96, height = px)
    if (!is.null(code)) register_code(id, code)
    for (fmt in c("png", "pdf")) local({
      fmt_ <- fmt
      output[[paste0(id, "_", fmt_)]] <- downloadHandler(
        filename = function() paste0("rasch_", id, ".", fmt_),
        content = function(file) {
          # 300 dpi PNG (and vector PDF) for publication. The device is
          # closed however the drawing ends: a plot that errors would
          # otherwise leave the device open for the rest of the session.
          before <- grDevices::dev.cur()
          if (fmt_ == "png") png(file, width = w, height = h, units = "in", res = 300)
          else pdf(file, width = w, height = h)
          on.exit(if (!identical(grDevices::dev.cur(), before))
            try(dev.off(), silent = TRUE), add = TRUE)
          fun()
        })
    })
  }
  # -------------------------------------------------- hover identification --
  # shared server-side plumbing for `plotCard(id, hover = TRUE)` (see its
  # definition for the client-side half: the hoverOpts plotOutput and the
  # floating uiOutput tip). Two shapes cover every hover-enabled plot in the
  # app: `register_hover_tip()` for ordinary scatter plots (nearPoints()
  # against a reproduction of the exact rows the plot draws), and
  # `register_hover_cormat()` for the symmetric item-by-item heatmaps, where
  # the hovered cell is looked up directly rather than nearest-matched.
  #
  # `data_fun()` must return a data frame with an `xvar`/`yvar` column pair in
  # exactly the same units and after exactly the same filtering as the plot
  # function itself (dropped NAs, excluded extremes, min-n screens, ...), so
  # the point nearest the cursor is always the point actually drawn there.
  # `label_fun(np)` formats the one-row match into the tooltip text.
  register_hover_tip <- function(id, data_fun, xvar, yvar, label_fun) {
    output[[paste0(id, "_tip")]] <- renderUI({
      hov <- input[[paste0(id, "_hover")]]
      req(hov)
      d <- data_fun()
      req(!is.null(d), nrow(d) > 0)
      np <- nearPoints(d, hov, xvar = xvar, yvar = yvar,
                       threshold = 12, maxpoints = 1)
      if (!nrow(np)) return(NULL)
      div(class = "rasch-hover-tip",
          style = sprintf("left:%.0fpx; top:%.0fpx;",
                          hov$coords_css$x, hov$coords_css$y - 12),
          label_fun(np))
    })
  }
  # `matrix_fun()` must return the same square, item-named matrix the plot
  # colours (before capping -- the tooltip shows the true value), so the cell
  # lookup mirrors plot_resid_cor()'s own layout: x is the column drawn
  # left-to-right, y is flipped top-to-bottom (image()'s row/column
  # convention), and the diagonal and lower triangle are masked out exactly
  # as the heatmap leaves them blank.
  register_hover_cormat <- function(id, matrix_fun, stat_label) {
    output[[paste0(id, "_tip")]] <- renderUI({
      hov <- input[[paste0(id, "_hover")]]
      req(hov)
      S <- matrix_fun()
      req(!is.null(S))
      L <- ncol(S)
      row <- round(hov$x); col <- L + 1L - round(hov$y)
      req(row >= 1, row <= L, col >= 1, col <= L)
      diag(S) <- NA
      S[lower.tri(S)] <- NA
      val <- S[row, col]
      if (is.na(val)) return(NULL)
      nms <- colnames(S)
      div(class = "rasch-hover-tip",
          style = sprintf("left:%.0fpx; top:%.0fpx;",
                          hov$coords_css$x, hov$coords_css$y - 12),
          sprintf("%s · %s · %s %.2f", nms[row], nms[col],
                  stat_label, val))
    })
  }
  # `csv_name` overrides the conventional rasch_<id>.csv download filename;
  # the CSV content is always the full table from `fun` (never the curated
  # on-screen display)
  register_table <- function(id, fun, dt_fun, code = NULL, csv_name = NULL) {
    output[[id]] <- renderDT(soft(dt_fun()))
    if (!is.null(code)) register_code(id, code)
    output[[paste0(id, "_csv")]] <- downloadHandler(
      filename = function() csv_name %||% paste0("rasch_", id, ".csv"),
      content = function(file) write_csv_plain(fun(), file))
  }
  # the curated stat boxes (test of fit, targeting, BTL test of fit) share
  # one registration: `ui_fun` builds the on-screen label-value rows, while
  # the CSV chip always downloads the COMPLETE table from `csv_fun`
  register_stat_box <- function(id, csv_fun, csv_name, ui_fun, code = NULL) {
    output[[id]] <- renderUI(ui_fun())
    if (!is.null(code)) register_code(id, code)
    output[[paste0(id, "_csv")]] <- downloadHandler(
      filename = function() csv_name,
      content = function(file) write_csv_plain(csv_fun(), file))
  }
  # paired PDF/ZIP batch downloads (all-items explorer plots, all-persons
  # kidmaps): one handler per extension, same content function for both
  # (the writer picks its format from the download file's extension)
  register_batch_download <- function(prefix, base, content) {
    for (ext in c("pdf", "zip")) local({
      ext_ <- ext
      output[[paste0(prefix, "_", ext_)]] <- downloadHandler(
        filename = function() paste0(base(), ".", ext_),
        content = content)
    })
  }
  # APA-leaning DT wrapper: Bootstrap 5 skin, right-aligned numerics, paging
  # controls only when the table needs them. Any fit_resid / infit_ms /
  # outfit_ms column is auto-flagged in the theme danger colour; `p_bold`
  # bolds p-values < .05.
  # curated display columns: fit objects carry every statistic, but the
  # tables show a readable core; the per-table "detailed columns" switch
  # reveals the rest (CSV downloads always contain everything)
  CORE <- list(
    items = c("item", "location", "se", "fit_resid", "infit_ms", "outfit_ms",
              "chisq", "df", "p_adj"),
    structural_item = c("item", "set", "location", "se",
                        "n", "fit_resid", "infit_ms", "outfit_ms"),
    person = c("id", "raw", "max_raw", "theta", "se", "extreme", "fit_resid"),
    dif_fact = c("item", "term", "F_uniform", "p_uniform_adj", "eta2_uniform",
                 "F_nonuniform", "p_nonuniform_adj", "eta2_nonuniform"),
    facet = c("level", "severity", "se", "n", "fit_resid"),
    btl_obj = c("object", "location", "se", "comparisons", "wins", "score",
                "infit_ms", "outfit_ms", "fit_resid", "extreme"),
    btl_judge = c("judge", "n", "infit_ms", "outfit_ms", "fit_resid", "misfit"),
    equate = c("item", "location_1", "location_2", "adj_difference", "t",
               "z", "df", "p_adj", "drift"),
    contrast = c("item", "contrast", "estimate", "se", "statistic", "p_adj"),
    frames = c("set", "group", "rho", "se_log_rho", "origin", "fit_resid",
               "n_responses"),
    compare = c("label", "model", "persons", "items", "judges", "objects",
                "eff_params", "cl_aic", "cl_bic", "two_delta_ll",
                "chisq_per_df", "item_fit_sd", "person_fit_sd", "PSI",
                "alpha", "OSI"),
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
    chisq_p_boot = "Boot p", chisq_p_boot_adj = "Boot adj. p",
    n_boot_chisq = "Chi-sq boot n",
    fit_resid_p_boot = "Fit resid boot p",
    fit_resid_p_boot_adj = "Fit resid boot adj. p",
    n_boot_fit_resid = "Fit resid boot n",
    p_anova_adj = "ANOVA adj. p", p_anova_bonf = "ANOVA Bonf. p",
    F_anova = "ANOVA F", F_uniform = "Uniform F",
    F_nonuniform = "Non-uniform F", p_uniform = "Uniform p",
    p_nonuniform = "Non-uniform p", p_uniform_adj = "Uniform adj. p",
    p_nonuniform_adj = "Non-uniform adj. p",
    p_uniform_boot = "Uniform boot p",
    p_uniform_boot_adj = "Uniform boot adj. p",
    p_nonuniform_boot = "Non-uniform boot p",
    p_nonuniform_boot_adj = "Non-uniform boot adj. p",
    eta2_uniform = "Uniform η²",
    eta2_nonuniform = "Non-uniform η²",
    eta2_partial = "Partial η²",
    uniform_DIF = "Uniform DIF", nonuniform_DIF = "Non-uniform DIF",
    min_judges = "Min judges",
    min_effective_judges = "Min effective judges",
    uniform_inference = "Uniform inference",
    nonuniform_inference = "Non-uniform inference",
    superseded = "Superseded", extreme = "Extreme", sum_sq = "Sum Sq",
    mean_sq = "Mean Sq",
    F_value = "F",
    mean_location = "Mean location", point_biserial = "Point-biserial",
    se_location = "SE", z_sep = "Separation z",
    alpha_drop = "α if deleted", item_total = "Item-total r",
    item_rest = "Item-rest r", di = "Discrimination", cum_pct = "Cum. %",
    exp_prop = "Expected", obs_prop = "Observed", obs_mean = "Observed mean",
    exp_mean = "Expected mean",
    exp_value = "Expected value", theta_mean = "Mean location",
    theta_max = "Max location", chisq_per_df = "Chi-sq/df",
    two_delta_ll = "2Δ log-lik", eff_params = "Eff. params",
    cl_aic = "CL-AIC", cl_bic = "CL-BIC", se_log_phi = "SE (log φ)",
    se_log_alpha = "SE (log α)", se_log_rho = "SE (log ρ)",
    mu = "Origin", comparisons = "Comparisons",
    obs_p = "Observed", est_p = "Expected", obs_t = "Threshold prop.",
    item_a = "Item A", item_b = "Item B", q3 = "Q3", q3_star = "Q3*",
    flagged = "Flagged", estimate = "Estimate (logits)",
    statistic = "Statistic", significant = "Significant",
    practical = "Practical", effect = "Effect",
    n_informative = "Informative comparisons")
  # p-value columns render as "<0.001" / 3 dp on the client, so sorting
  # still uses the raw value; detection runs on the ORIGINAL column names
  P_COL_RE <- "^p$|^p_|_p$|^prob$|p_anova|p_adj|p_bonf|p_uniform|p_nonuniform|p_boot"
  P_RENDER <- DT::JS("function(data,type,row){ if(type==='display'){ if(data===null||data==='') return ''; var x=Number(data); return x<0.001 ? '&lt;0.001' : x.toFixed(3);} return data; }")
  # styleInterval() includes the cut point in its lower interval. That is the
  # wrong boundary for p < alpha and makes symmetric fit bands asymmetric at
  # their negative cut. Use the same strict decisions as the result objects.
  colour_if <- function(condition) DT::JS(paste0(
    "function(value){if(value===null||value==='')return 'inherit';",
    "var x=Number(value);if(!Number.isFinite(x))return 'inherit';return ",
    condition, "?'var(--bs-danger)':'inherit';}"))
  weight_if <- function(condition) DT::JS(paste0(
    "function(value){if(value===null||value==='')return 'normal';",
    "var x=Number(value);if(!Number.isFinite(x))return 'normal';return ",
    condition, "?'bold':'normal';}"))
  # fit flags, consistent across every model table: a fit residual beyond
  # |2.5|, an outfit mean square outside 0.7-1.3, and an infit mean square
  # outside the tighter 0.8-1.2 (infit is information-weighted, so it varies
  # less under fit; conventional working bands). `orig` is the ORIGINAL
  # column-name vector of the displayed frame (positions match the widget).
  flag_fit_cols <- function(dt, orig) {
    for (fc in which(orig == "fit_resid"))
      dt <- formatStyle(dt, fc, color = colour_if("Math.abs(x)>2.5"))
    for (oc in which(orig == "outfit_ms"))
      dt <- formatStyle(dt, oc, color = colour_if("x<0.7||x>1.3"))
    for (ic in which(orig == "infit_ms"))
      dt <- formatStyle(dt, ic, color = colour_if("x<0.8||x>1.2"))
    dt
  }
  num_dt <- function(d, digits = 3, p_bold = NULL,
                     page_len = 15, paging = NULL, one_dp = NULL, ...) {
    orig <- names(d)
    # unname: which() on named logicals yields named position vectors, and
    # jsonlite warns whenever one reaches the widget payload
    num <- unname(vapply(d, is.numeric, TRUE))
    # integer-valued columns (counts, whole-number df) show no decimals;
    # fractional df columns fail the test and keep the 3-dp rounding
    intcol <- unname(vapply(d, function(v)
      is.numeric(v) && all(is.na(v) | v == round(v)), TRUE))
    pcol <- num & grepl(P_COL_RE, orig)
    # pager and info line appear only when the table overflows one page;
    # `paging` can force either way per table
    if (is.null(paging)) paging <- nrow(d) > page_len
    opts <- list(pageLength = if (paging) page_len else max(nrow(d), 1L),
                 scrollX = TRUE,
                 dom = if (paging) "tip" else "t")
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
    # unname: a named names-vector rides into the widget payload and makes
    # jsonlite warn on every table render
    names(d)[hit] <- unname(DISPLAY_NAMES[orig[hit]])
    dt <- datatable(d, rownames = FALSE, style = "bootstrap5",
                    class = "table-sm compact hover order-column", ...,
                    options = opts)
    rnd <- which(num & !intcol & !pcol)
    # a named subset (e.g. the composite-likelihood ICs) can round to one
    # decimal place instead of the table's default; each column is targeted
    # by exactly one formatRound() call, so the two never compete
    one <- orig %in% one_dp
    rnd1 <- intersect(rnd, which(one))
    rnd3 <- setdiff(rnd, rnd1)
    if (length(rnd3)) dt <- formatRound(dt, rnd3, digits)
    if (length(rnd1)) dt <- formatRound(dt, rnd1, 1)
    dt <- flag_fit_cols(dt, orig)
    for (pc in which(orig %in% p_bold))
      dt <- formatStyle(dt, pc,
                        fontWeight = weight_if("x<0.05"))
    dt
  }
  # Lower-triangular display of a symmetric matrix: the redundant strictly-upper
  # cells are blanked so each item pair is read once, the first column carries
  # the item name, and NA cells (e.g. the Q3* diagonal) show empty. Values are
  # pre-formatted, so the columns are right-aligned by column definition rather
  # than by numeric class.
  tri_dt <- function(M, digits = 2, flagged = NULL) {
    # tiny magnitudes (including negative zero) print as a clean 0.00
    M[!is.na(M) & abs(M) < 0.5 * 10^(-digits)] <- 0
    disp <- formatC(M, format = "f", digits = digits)
    dim(disp) <- dim(M); dimnames(disp) <- dimnames(M)
    # flagged pairs (Q3* above the threshold) are shown red in place of the
    # heatmap's flag mark; wrap before the triangle is blanked so only kept
    # cells carry the span
    if (!is.null(flagged)) {
      hot <- flagged & !is.na(M) & lower.tri(M, diag = FALSE)
      hot[is.na(hot)] <- FALSE
      disp[hot] <- sprintf(
        '<span style="color:var(--bs-danger);font-weight:600">%s</span>',
        disp[hot])
    }
    disp[upper.tri(M)] <- ""                 # redundant upper triangle
    disp[is.na(M)] <- ""                     # empty diagonal / missing pairs
    df <- data.frame(item = rownames(M), disp, check.names = FALSE,
                     stringsAsFactors = FALSE)
    # escape only the item column (user-controlled labels); the numeric cells
    # intentionally carry the red-flag <span>, so they must not be escaped
    dt <- datatable(df, rownames = FALSE, escape = "item", style = "bootstrap5",
              class = "table hover order-column",
              options = list(paging = FALSE, ordering = FALSE,
                             scrollX = TRUE, dom = "t",
                             columnDefs = list(list(className = "dt-right",
                               targets = seq_len(ncol(M))))))
    # a roomier grid than the compact default, so the matrix reads at a weight
    # closer to its heatmap alongside it
    formatStyle(dt, names(df), fontSize = "1rem",
                paddingTop = "7px", paddingBottom = "7px")
  }
  # Red-highlight a triggering value in place of a boolean flag column (the
  # items table pattern): colour the cell with the theme danger colour when
  # it crosses the threshold, leaving in-range values at the inherited text
  # colour. `d` is the DISPLAYED data frame, so column positions match the
  # rendered table; each helper is a no-op when the column is absent.
  style_lo_red <- function(dt, d, col, cut)     # red when value < cut (e.g. p)
    Reduce(function(x, j) formatStyle(x, j,
      color = colour_if(sprintf("x<%.17g", cut))),
      which(names(d) == col), dt)
  style_mag_red <- function(dt, d, col, mag)    # red when |value| >= mag
    Reduce(function(x, j) formatStyle(x, j,
      color = colour_if(sprintf("Math.abs(x)>=%.17g", mag))),
      which(names(d) == col), dt)

  # ------------------------------------------------ navbar status summary --
  # A single quiet capsule replaces the row of coloured bubbles. On narrower
  # screens the sample-size details drop away, leaving model and reliability.
  item_count_app <- function(f) {
    if (inherits(f, "rasch_mfrm")) nrow(f$item_effects)
    else if (inherits(f, "rasch_efrm")) nrow(f$item_arbitrary)
    else ncol(f$X)
  }
  alpha_design_applicable <- function(f)
    .classical_design_applicable(f)
  output$nav_status <- renderUI({
    sep <- span(class = "rasch-nav-sep", "·")
    b <- tryCatch(bfit(), error = function(e) NULL)
    if (!is.null(b)) {
      osi <- b$osi$PSI
      return(div(class = "rasch-nav-summary",
        span(class = "rasch-nav-model",
             if (inherits(b, "rasch_btl_efrm")) "CJ Extended Frames"
             else if (inherits(b, "rasch_btl_explanatory")) "Explanatory CJ"
             else "CJ"),
        span(class = "rasch-nav-secondary", sep,
             paste(nrow(b$objects), "objects"), sep,
             sprintf("%.0f comparisons", b$n_comparisons)),
        sep,
        span(class = if (finite1(osi) && osi >= 0.7)
          "rasch-nav-good" else "rasch-nav-warn",
          if (finite1(osi)) sprintf("OSI %.2f", osi) else "OSI —")))
    }
    f <- fit_or_null()
    if (is.null(f)) return(NULL)
    psi <- f$psi$PSI
    # the chip states what was actually fitted: PCM/RSM only make sense for
    # polytomous items, so an all-dichotomous fit reads "Dichotomous"
    model_lab <- if (inherits(f, "rasch_mfrm")) "MFRM"
      else if (inherits(f, "rasch_efrm")) "EFRM"
      else if (inherits(f, "rasch_explanatory")) f$explanatory_model
      else if (max(f$m) == 1L) "Dichotomous"
      else f$model
    div(class = "rasch-nav-summary",
      span(class = "rasch-nav-model", model_lab),
      span(class = "rasch-nav-secondary", sep,
           paste(nrow(f$X), "persons"), sep,
           paste(item_count_app(f), "items")),
      sep,
      span(class = if (finite1(psi) && psi >= 0.7)
        "rasch-nav-good" else "rasch-nav-warn",
        if (finite1(psi)) sprintf("PSI %.2f", psi) else "PSI —"))
  })

  # -------------------------------------------------------------- summary --
  output$vboxes <- renderUI({
    f <- fit()
    metric_grid(
      metric_tile("metric_persons", "Persons", nrow(f$X),
                  icon = "distribution", status = "person"),
      metric_tile("metric_items", "Items", item_count_app(f),
                  icon = "ruler", status = "item"),
      metric_tile("metric_psi", "PSI",
                  if (finite1(f$psi$PSI)) sprintf("%.3f", f$psi$PSI) else "—",
                  if (finite1(f$psi_noext$PSI))
                    sprintf("%.3f without extremes", f$psi_noext$PSI)
                  else "Without extremes —",
                  icon = "separation",
                  status = if (!finite1(f$psi$PSI)) "neutral"
                    else if (f$psi$PSI >= 0.7) "good" else "bad"),
      metric_tile("metric_alpha", "Alpha",
                  if (!alpha_design_applicable(f)) "—"
                  else if (finite1(f$alpha$alpha)) sprintf("%.3f", f$alpha$alpha)
                  else "—",
                  if (!alpha_design_applicable(f))
                    "Not applicable when items have several response cells"
                  else if (isFALSE(f$alpha$applicable))
                    sprintf("Complete cases only · n = %d", f$alpha$n)
                  else sprintf("Complete cases · n = %d", f$alpha$n),
                  icon = "alpha",
                  status = if (!alpha_design_applicable(f) ||
                               !finite1(f$alpha$alpha))
                    "neutral" else if (f$alpha$alpha >= 0.7)
                    "good" else "bad"),
      metric_tile("metric_item_trait",
                  if (inherits(f, c("rasch_mfrm", "rasch_efrm")))
                    "Approx. response-cell-trait p" else
                      "Approx. item-trait p",
                  if (!.has_repeated_residual_units(f) &&
                      finite1(f$total_chisq_p)) fmt_p(f$total_chisq_p) else "—",
                  icon = "chisq", status = "neutral"),
      metric_tile("metric_power", "Separation quality",
                  f$separation_quality %||% f$power_of_fit,
                  icon = "power", status = "neutral")
    )
  })
  register_code("vboxes", function() paste(
    "list(",
    paste0("  persons = nrow(fit$X), items = ",
      if (inherits(fit(), "rasch_mfrm")) "nrow(fit$item_effects)," else
      if (inherits(fit(), "rasch_efrm")) "nrow(fit$item_arbitrary)," else
        "ncol(fit$X),"),
    "  PSI = fit$psi$PSI, PSI_without_extremes = fit$psi_noext$PSI,",
    "  alpha = if (inherits(fit, c('rasch_mfrm', 'rasch_efrm')) && !isTRUE(fit$alpha$design_applicable)) NA_real_ else fit$alpha$alpha,",
    "  item_trait_p = fit$total_chisq_p,",
    "  separation_quality = if (!is.null(fit$separation_quality)) fit$separation_quality else fit$power_of_fit",
    ")", sep = "\n"))

  # test-of-fit and targeting/reliability summaries as curated stat boxes
  # (values read off the fit object directly); the CSV chips download the
  # COMPLETE fit_summary_table / targeting_table with raw column names.
  # The score-to-measure table stays available via score_table(fit) and the
  # model-results ZIP; thresholds live on the Items explorer Thresholds tab.
  register_stat_box("fitsum_tbl",
    csv_fun = function() fit_summary_table(fit()),
    csv_name = "fit_summary.csv",
    ui_fun = function() {
      f <- fit()
      est <- f$est %||% f
      dis <- names(which(vapply(f$thresholds_diag, function(d)
        !d$ordered && length(d$thresholds) > 1, TRUE)))
      link_failed <- inherits(f, "rasch_efrm") &&
        isTRUE(est$stage1_converged) && !isTRUE(est$converged)
      conv <- if (link_failed)
        span(class = "text-danger",
             "conditional stage converged; set link did not")
      else if (isTRUE(est$converged))
        sprintf("converged in %d iterations", est$iterations)
      else span(class = "text-danger",
                sprintf("did not converge in %d iterations", est$iterations))
      method <- .estimation_label(f)
      repeated_ids <- .has_repeated_residual_units(f)
      item_inference <- inference_count(if (repeated_ids)
        rep(NA_real_, nrow(f$items)) else f$items$p_adj)
      tagList(
        div(class = "stat-head",
            if (inherits(f, "rasch_explanatory")) f$explanatory_model else f$model,
            " · ", method, " · ", conv),
        stat_rows(
          stat_row(if (inherits(f, c("rasch_mfrm", "rasch_efrm")))
                     "Approx. response-cell-trait chi-square" else
                       "Approx. item-trait chi-square",
                   sprintf("%.2f on %d df, %s", f$total_chisq, f$total_df,
                           p_lab(if (repeated_ids) NA_real_
                                 else f$total_chisq_p))),
          stat_row(if (inherits(f, c("rasch_mfrm", "rasch_efrm")))
                     "Response-cell fit residual" else "Item fit residual",
                   sprintf("mean %.2f, SD %.2f", f$item_fit_summary$mean,
                           f$item_fit_summary$sd)),
          stat_row("Person fit residual",
                   sprintf("mean %.2f, SD %.2f", f$person_fit_summary$mean,
                           f$person_fit_summary$sd)),
          stat_row(if (inherits(f, c("rasch_mfrm", "rasch_efrm")))
                     "Response-cell location-residual correlation" else
                       "Item location-residual correlation",
                   { v <- f$summary_stats$cor_item_fit_location
                     if (finite1(v)) sprintf("%.3f", v) else "NA" }),
          stat_row("Person location-residual correlation",
                   { v <- f$summary_stats$cor_person_fit_location
                     if (finite1(v)) sprintf("%.3f", v) else "NA" }),
          stat_row(if (inherits(f, c("rasch_mfrm", "rasch_efrm")))
                     "Response cells with approx. Holm p < .05" else
                       "Items with approx. Holm p < .05",
                   item_inference$text),
          stat_row(if (inherits(f, c("rasch_mfrm", "rasch_efrm")))
                     "Disordered response-cell thresholds" else
                       "Disordered thresholds",
                   if (length(dis)) paste(dis, collapse = ", ") else "none")))
    },
    code = function() "fit_summary_table(fit)")
  # routine handling notes (the old text panel printed fit$notes)
  output$fitsum_notes <- renderUI({
    f <- fit()
    if (!length(f$notes)) return(NULL)
    sprintf("Note. %s.", paste(f$notes, collapse = "; "))
  })
  register_stat_box("targeting_tbl",
    csv_fun = function() targeting_table(fit()),
    csv_name = "targeting.csv",
    ui_fun = function() {
      f <- fit(); ss <- f$summary_stats; tg <- f$targeting
      psi_txt <- if (!finite1(f$psi$PSI)) "—"
      else if (finite1(f$psi_noext$PSI))
        sprintf("%.3f (%.3f without extremes)", f$psi$PSI, f$psi_noext$PSI)
      else sprintf("%.3f", f$psi$PSI)
      alpha_txt <- if (!alpha_design_applicable(f))
        "not applicable when items have several response cells"
      else if (finite1(f$alpha$alpha)) sprintf("%.3f", f$alpha$alpha)
      else sprintf("not applicable (%d complete cases)", f$alpha$n %||% 0L)
      stat_rows(
        stat_row("Person location",
                 sprintf("mean %.2f, SD %.2f", ss$person_location$mean,
                         ss$person_location$sd)),
        stat_row(if (inherits(f, c("rasch_mfrm", "rasch_efrm")))
                   "Response-cell location" else "Item location",
                 sprintf("SD %.2f (mean constrained to 0)",
                         ss$item_location$sd)),
        stat_row(if (inherits(f, c("rasch_mfrm", "rasch_efrm")))
                   "Calibration threshold range" else "Threshold range",
                 sprintf("%.2f to %.2f", tg$threshold_range[1],
                         tg$threshold_range[2])),
        stat_row("Persons beyond thresholds",
                 sprintf("%.1f%% below · %.1f%% above",
                         100 * tg$prop_below, 100 * tg$prop_above)),
        stat_row("PSI", psi_txt),
        stat_row(if (inherits(f, c("rasch_mfrm", "rasch_efrm")))
                   "Response-cell separation" else "Item reliability",
                 { v <- f$isi$PSI
                   if (finite1(v)) sprintf("%.3f", v) else "NA" }),
        stat_row("Person strata",
                 { v <- f$psi$strata
                   if (finite1(v)) sprintf("%.1f", v) else "NA" }),
        stat_row("Coefficient alpha", alpha_txt))
    },
    code = function() "targeting_table(fit)")

  # traditional (CTT) statistics: shown on the Items page (last accordion
  # panel), with the header line, CSV, and code footer registered here
  ctt_res <- reactive(tryCatch(ctt_table(fit()), error = function(e) e))
  register_code("ctt_tbl", function() "ctt_table(fit)$table")
  output$ctt_head <- renderUI({
    ct <- ctt_res()
    if (inherits(ct, "error"))
      return(p(class = "text-muted",
               paste("Traditional statistics unavailable:", conditionMessage(ct))))
    # total-score summaries need complete responders; with structural
    # missingness (linked forms) there may be none, so the header falls back
    # to the available-case framing instead of printing NA values
    if (!finite1(ct$mean))
      return(p(class = "small mb-2", HTML(sprintf(
        "Per-item statistics use available cases (item n %d&ndash;%d). Too few complete responders (n = %d) for the total-score mean, SD, SEM%s.",
        ct$n_range[1], ct$n_range[2], ct$n,
        if (finite1(ct$alpha)) sprintf("; alpha <b>%.3f</b>", ct$alpha)
        else ", and alpha"))))
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
    filename = function() "rasch_ctt_tbl.csv",
    content = function(file) {
      ct <- ctt_res(); req(!inherits(ct, "error"))
      write_csv_plain(ct$table, file)
    })

  # likelihood-ratio test of PCM against the rating parameterisation; only
  # meaningful for a PCM fit whose items share a common maximum score > 1.
  # The bottom row is server-built so the card hides when it does not apply.
  lr_res <- reactiveVal(NULL)
  output$summary_bottom <- renderUI({
    f <- fit()
    lr_applies <- identical(f$model, "PCM") && length(unique(f$m)) == 1L &&
      max(f$m) >= 2L
    if (!lr_applies) return(NULL)
    layout_columns(col_widths = breakpoints(sm = 12, xl = 6), .lr_card())
  })
  register_code("lr", function() "lr_test(fit)")
  observeEvent(input$run_lr, {
    f <- fit()
    r <- withProgress(message = "Refitting with the rating parameterisation…",
                      value = 0.4,
                      tryCatch(lr_test(f), error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("LR test failed:", conditionMessage(r)),
                       type = "error", duration = 10)
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
  # ---------------------------------------------------------------- items --
  structural_items <- reactive({
    f <- fit()
    if (inherits(f, "rasch_mfrm")) return(f$item_effects)
    if (inherits(f, "rasch_efrm")) return(f$item_arbitrary)
    NULL
  })
  output$items_table_title <- renderUI({
    if (inherits(fit(), c("rasch_mfrm", "rasch_efrm")))
      "Response-cell statistics" else "Item statistics"
  })

  # Explanatory calibration: the same fit object feeds the comparison,
  # coefficient, diagnostic and downstream pages. Diagnostics are cached until
  # the active calibration changes because each row requires a constrained
  # refit.
  expl_fit <- reactive({
    b <- tryCatch(btl_fit(), error = function(e) NULL)
    if (!is.null(b)) {
      s <- active_btl_step()
      f <- if (is.null(s)) b else s$fit
    } else f <- fit()
    validate(need(inherits(f, c("rasch_explanatory",
                                "rasch_btl_explanatory")),
                  "Fit an explanatory model to use this page."))
    f
  })
  expl_is_cj <- reactive(inherits(expl_fit(), "rasch_btl_explanatory"))
  expl_coef <- reactive({
    f <- expl_fit()
    if (inherits(f, "rasch_btl_explanatory")) f$object_coefficients
    else f$est$coefficients
  })
  expl_diag <- reactive(explanatory_diagnostics(expl_fit(), p_adjust = "holm"))
  output$expl_boxes <- renderUI({
    f <- expl_fit(); tst <- explanatory_test(f)
    p <- tst$p_kent[1L]
    metric_grid(
      metric_tile("metric_expl_model", "Model",
                  if (inherits(f, "rasch_btl_explanatory"))
                    "Explanatory CJ" else f$explanatory_model,
                  icon = "distribution", status = "accent"),
      metric_tile("metric_expl_terms", "Predictor effects",
                  nrow(expl_coef()), icon = "ruler"),
      metric_tile("metric_expl_relax", "Fixed departures",
                  nrow(f$explanatory$relaxations), icon = "separation"),
      metric_tile("metric_expl_r2", "Calibration R²",
                  if (is.finite(tst$r_squared[1L]))
                    sprintf("%.2f", tst$r_squared[1L]) else "—",
                  icon = "graph-up"),
      metric_tile("metric_expl_test", "Against free calibration", p_lab(p),
                  icon = "chisq",
                  status = if (!is.finite(p)) "neutral"
                    else if (p < .05) "bad" else "good"))
  })
  expl_object_name <- reactive(if (expl_is_cj()) "bt" else "fit")
  register_code("expl_boxes", function() sprintf(
    "list(formula = %s$explanatory$formula_text, test = explanatory_test(%s))",
    expl_object_name(), expl_object_name()))
  register_table("expl_test_tbl", function() explanatory_test(expl_fit()),
    function() {
      z <- explanatory_test(expl_fit())
      z <- z[, c("model", "parameters", "free_parameters",
                 "r_squared", "chisq_kent", "df", "p"), drop = FALSE]
      num_dt(z, p_bold = "p")
    },
    code = function() sprintf("explanatory_test(%s)", expl_object_name()))
  register_table("expl_coef_tbl", function() expl_coef(),
    function() num_dt(expl_coef(),
                       p_bold = "p_adj"),
    code = function() paste0(expl_object_name(),
      if (expl_is_cj()) "$object_coefficients" else "$est$coefficients"))
  register_table("expl_diag_tbl", function() expl_diag(), function() {
    num_dt(expl_diag(), page_len = 15, selection = "single",
           p_bold = "p_adj")
  }, code = function() sprintf(
    "explanatory_diagnostics(%s, p_adjust = \"holm\")", expl_object_name()))
  register_table("expl_relax_tbl", function() expl_fit()$explanatory$relaxations,
    function() num_dt(expl_fit()$explanatory$relaxations),
    code = function() paste0(expl_object_name(), "$explanatory$relaxations"))

  expl_selected <- reactive({
    i <- input$expl_diag_tbl_rows_selected
    d <- expl_diag()
    validate(need(length(i) == 1L && i >= 1L && i <= nrow(d),
                  "Select one diagnostic row."))
    d[i, , drop = FALSE]
  })
  output$expl_selected_note <- renderUI({
    z <- tryCatch(expl_selected(), error = function(e) NULL)
    if (is.null(z))
      return(p(class = "text-muted small mt-2", "No departure selected."))
    name <- if ("object" %in% names(z)) z$object else z$item
    p(class = "small mt-2",
      sprintf("Selected: %s, %s; adjusted p = %s.", name,
              tolower(z$component), fmt_p(z$p_adj)))
  })
  observeEvent(input$expl_relax, {
    z <- expl_selected()
    if (expl_is_cj()) {
      f <- tryCatch(relax_btl_explanatory(expl_fit(), z$object),
                    error = function(e) e)
      if (inherits(f, "error")) {
        showNotification(paste("Refit failed:", conditionMessage(f)),
                         type = "error", duration = 10)
        return()
      }
      lab <- sprintf("fixed explanatory departure for %s", z$object)
      push_btl_analysis_step(
        type = "explanatory_relax", label = lab, fitted = f,
        details = as.list(z),
        code = sprintf("bt <- relax_btl_explanatory(bt, %s)",
                       qstr(z$object)))
      showNotification(paste("Added", lab, "and repeated the calibration."),
                       type = "message", duration = 6)
      return()
    }
    component <- if (identical(z$component, "Item location"))
      "location" else "thresholds"
    f <- tryCatch(relax_explanatory(expl_fit(), z$item, component),
                  error = function(e) e)
    if (inherits(f, "error")) {
      showNotification(paste("Refit failed:", conditionMessage(f)),
                       type = "error", duration = 10)
      return()
    }
    lab <- sprintf("fixed explanatory departure for %s (%s)",
                   z$item, tolower(z$component))
    push_analysis_step(
      type = "explanatory_relax", label = lab, fitted = f,
      details = as.list(z),
      code = sprintf("fit <- relax_explanatory(fit, %s, component = %s)",
                     qstr(z$item), qstr(component)))
    showNotification(paste("Added", lab, "and repeated the calibration."),
                     type = "message", duration = 6)
  })
  expl_plot_data <- reactive({
    f <- expl_fit()
    if (inherits(f, "rasch_btl_explanatory")) {
      a <- f$objects; b <- f$reference_fit$objects
      j <- match(a$object, b$object)
      return(data.frame(name = a$object, explanatory = a$location,
                        free = b$location[j], label = a$object,
                        stringsAsFactors = FALSE))
    }
    a <- f$thresholds[, c("item", "k", "tau"), drop = FALSE]
    b <- f$reference_fit$thresholds[, c("item", "k", "tau"), drop = FALSE]
    a$item_name <- colnames(f$X)[a$item]
    b$item_name <- colnames(f$X)[b$item]
    a$key <- paste(a$item_name, a$k, sep = "\r")
    b$key <- paste(b$item_name, b$k, sep = "\r")
    j <- match(a$key, b$key)
    data.frame(item = a$item_name, threshold = a$k,
               explanatory = a$tau, free = b$tau[j],
               label = paste0(a$item_name, " [", a$k, "]"),
               stringsAsFactors = FALSE)
  })
  draw_expl_calibration <- function() {
    d <- expl_plot_data()
    lim <- range(c(d$explanatory, d$free), finite = TRUE)
    pad <- max(diff(lim) * .06, .15); lim <- lim + c(-pad, pad)
    unit <- if (expl_is_cj()) "object location" else "threshold"
    plot(d$free, d$explanatory, pch = 19, col = "#276FBF",
         xlab = paste("Free", unit, "(logits)"),
         ylab = paste("Explanatory", unit, "(logits)"),
         xlim = lim, ylim = lim)
    abline(0, 1, lty = 2, col = "#69727D")
    nlab <- min(5L, nrow(d))
    take <- order(abs(d$explanatory - d$free), decreasing = TRUE)[seq_len(nlab)]
    text(d$free[take], d$explanatory[take], d$label[take],
         pos = 4, cex = .75, xpd = NA)
  }
  register_plot("expl_calibration", draw_expl_calibration,
    code = function() if (expl_is_cj()) paste(
      "active <- bt$objects", "free <- bt$reference_fit$objects",
      "plot(free$location, active$location,",
      "     xlab = 'Free object location (logits)',",
      "     ylab = 'Explanatory object location (logits)')",
      "abline(0, 1, lty = 2)", sep = "\n") else paste(
        "active <- fit$thresholds", "free <- fit$reference_fit$thresholds",
        "plot(free$tau, active$tau, xlab = 'Free threshold (logits)',",
        "     ylab = 'Explanatory threshold (logits)')",
        "abline(0, 1, lty = 2)", sep = "\n"))
  register_table("structural_items_tbl", function() structural_items(),
    function() {
      d <- structural_items(); req(!is.null(d))
      num_dt(curate(d, "structural_item",
                    full = isTRUE(input$structural_items_full)), page_len = 10)
    }, code = function() {
      if (inherits(fit(), "rasch_mfrm")) "fit$item_effects" else
        "fit$item_arbitrary"
    })
  # ---- bootstrap null for the fit statistics (post-estimation, on request)
  boot_val <- reactiveVal(NULL)
  dif_boot_val <- reactiveVal(NULL)
  observeEvent(fit(), {
    if (!isTRUE(restoring_project())) {
      boot_val(NULL)
      dif_boot_val(NULL)
    }
  }, ignoreInit = TRUE)
  observeEvent(list(btl_fit(), active_btl_step()), {
    if (!isTRUE(restoring_project())) {
      boot_val(NULL)
      dif_boot_val(NULL)
      # Judge-group DIF belongs to the active BTL calibration, not merely the
      # base fit. Explanatory relaxations and undo operations change bfit()
      # through active_btl_step() while leaving btl_fit() untouched.
      bdif_res(NULL)
      bdif_meta(NULL)
    }
  }, ignoreInit = TRUE)
  output$can_boot <- reactive({
    f <- fit_or_null()
    !is.null(f) && inherits(f, "rasch") &&
      !inherits(f, c("rasch_mfrm", "rasch_efrm", "rasch_explanatory")) &&
      is.null(f$refit_spec$pc_components) &&
      (is.null(f$disc) || length(unique(f$disc)) == 1L)
  })
  outputOptions(output, "can_boot", suspendWhenHidden = FALSE)
  output$can_btl_boot <- reactive({
    f <- tryCatch(bfit(), error = function(e) NULL)
    !is.null(f) && inherits(f, "rasch_btl") &&
      !inherits(f, "rasch_btl_efrm") && !is.null(f$refit_spec)
  })
  outputOptions(output, "can_btl_boot", suspendWhenHidden = FALSE)
  output$can_dif_boot <- reactive({
    f <- fit_or_null()
    testable <- if (is.null(f) || is.null(f$factors)) character(0) else
      setdiff(names(f$factors), f$frame_group %||% character(0))
    !is.null(f) && length(testable) > 0L && inherits(f, "rasch") &&
      is.null(f$refit_spec$pc_components) &&
      (inherits(f, "rasch_efrm") || is.null(f$disc) ||
         length(unique(f$disc)) == 1L)
  })
  outputOptions(output, "can_dif_boot", suspendWhenHidden = FALSE)
  output$dif_boot_explainer <- renderUI({
    f <- fit_or_null()
    txt <- if (inherits(f, "rasch_efrm")) paste(
      "Generates responses conditional on each person's observed subtotal",
      "within each item set, then repeats the complete Extended Frames",
      "calibration and DIF analysis.") else paste(
        "Generates responses conditional on each person's observed score",
        "and missingness pattern, then repeats the calibration and complete",
        "DIF analysis.")
    accordion_info(paste(txt,
      "The adjusted residual ANOVA remains the primary result. Bootstrap",
      "familywise probabilities refer to the fitted global invariant null;",
      "they do not provide strong control when another item has DIF."))
  })

  boot_job <- reactiveVal(NULL)
  close_boot_job <- function(st) {
    if (!is.null(st$progress)) try(st$progress$close(), silent = TRUE)
    unlink(st$log_file, force = TRUE)
    invisible(NULL)
  }
  cancel_boot_job <- function() {
    st <- isolate(boot_job())
    if (is.null(st)) return(invisible(FALSE))
    if (st$process$is_alive()) stop_efrm_process(st$process)
    close_boot_job(st); boot_job(NULL)
    invisible(TRUE)
  }
  output$boot_job_controls <- renderUI({
    st <- boot_job()
    if (is.null(st) || !identical(st$kind, "rasch")) return(NULL)
    actionButton("cancel_boot", "Cancel", icon = bs_icon("x-circle"),
                 class = "btn-outline-danger btn-sm mb-2")
  })
  output$btl_boot_job_controls <- renderUI({
    st <- boot_job()
    if (is.null(st) || !identical(st$kind, "btl")) return(NULL)
    actionButton("cancel_btl_boot", "Cancel", icon = bs_icon("x-circle"),
                 class = "btn-outline-danger btn-sm mb-2")
  })
  output$dif_boot_job_controls <- renderUI({
    st <- boot_job()
    if (is.null(st) || !identical(st$kind, "dif")) return(NULL)
    actionButton("cancel_dif_boot", "Cancel", icon = bs_icon("x-circle"),
                 class = "btn-outline-danger btn-sm mb-2")
  })
  output$bdif_boot_job_controls <- renderUI({
    st <- boot_job()
    if (is.null(st) || !identical(st$kind, "bdif")) return(NULL)
    actionButton("cancel_bdif_boot", "Cancel", icon = bs_icon("x-circle"),
                 class = "btn-outline-danger btn-sm mb-2")
  })
  outputOptions(output, "boot_job_controls", suspendWhenHidden = FALSE)
  outputOptions(output, "btl_boot_job_controls", suspendWhenHidden = FALSE)
  outputOptions(output, "dif_boot_job_controls", suspendWhenHidden = FALSE)
  outputOptions(output, "bdif_boot_job_controls", suspendWhenHidden = FALSE)
  observeEvent(input$cancel_boot, {
    if (cancel_boot_job())
      showNotification("Fit bootstrap cancelled.", type = "message", duration = 5)
  })
  observeEvent(input$cancel_btl_boot, {
    if (cancel_boot_job())
      showNotification("Fit bootstrap cancelled.", type = "message", duration = 5)
  })
  observeEvent(input$cancel_dif_boot, {
    if (cancel_boot_job())
      showNotification("DIF bootstrap cancelled.", type = "message", duration = 5)
  })
  observeEvent(input$cancel_bdif_boot, {
    if (cancel_boot_job())
      showNotification("DIF bootstrap cancelled.", type = "message", duration = 5)
  })
  session$onSessionEnded(function() {
    st <- isolate(boot_job())
    if (!is.null(st) && st$process$is_alive()) stop_efrm_process(st$process)
    if (!is.null(st)) close_boot_job(st)
  })

  start_bootstrap <- function(f, B_raw, seed_raw, kind, dif = NULL) {
    if (!is.null(efrm_job()) || !is.null(btlef_job())) {
      showNotification("Another background analysis is running. Cancel it before starting a bootstrap.",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    if (!is.null(boot_job())) {
      showNotification("A bootstrap is already running. Cancel it before starting another.",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    if (identical(kind, "rasch") && !inherits(f, "rasch_btl") &&
        .has_repeated_residual_units(f)) {
      showNotification(paste(
        "The item-fit bootstrap is unavailable when person IDs repeat because",
        "its reference distribution assumes independent response rows.",
        "Analyse occasions separately for item-fit inference."),
        type = "warning", duration = 10)
      return(invisible(NULL))
    }
    B_raw <- suppressWarnings(as.numeric(B_raw))
    if (length(B_raw) != 1L || !is.finite(B_raw) || B_raw != floor(B_raw) ||
        B_raw < 99 || B_raw > 9999) {
      showNotification("Replicates must be a whole number from 99 to 9999",
                       type = "warning")
      return(invisible(NULL))
    }
    B <- as.integer(B_raw)
    seed_raw <- suppressWarnings(as.numeric(seed_raw))
    if (length(seed_raw) != 1L || !is.finite(seed_raw) || seed_raw < 0 ||
        seed_raw != floor(seed_raw) || seed_raw > .Machine$integer.max) {
      showNotification("Seed must be one non-negative whole number within the integer range",
                       type = "warning")
      return(invisible(NULL))
    }
    seed <- as.integer(seed_raw)
    log_file <- tempfile("rasch-bootstrap-", fileext = ".log")
    progress <- shiny::Progress$new(session, min = 0, max = 1)
    progress$set(message = sprintf("Bootstrapping %s (B = %d)",
                                   if (kind %in% c("dif", "bdif"))
                                     "DIF" else "fit", B),
                 detail = "Refitting the simulated datasets", value = .5)
    process <- tryCatch(callr::r_bg(
      function(fit, dif, B, seed, kind, source_dir) {
        if (!is.null(source_dir)) {
          if (!requireNamespace("pkgload", quietly = TRUE))
            stop("pkgload is needed for a source-tree background analysis")
          pkgload::load_all(dirname(source_dir), quiet = TRUE)
        }
        warnings <- character(0)
        value <- tryCatch(withCallingHandlers(
          do.call(if (kind %in% c("dif", "bdif")) {
              if (exists("dif_bootstrap", inherits = TRUE))
                get("dif_bootstrap", inherits = TRUE) else
                  getExportedValue("rasch", "dif_bootstrap")
            } else if (exists("fit_bootstrap", inherits = TRUE))
              get("fit_bootstrap", inherits = TRUE) else
                getExportedValue("rasch", "fit_bootstrap"),
            c(list(fit = fit), if (kind %in% c("dif", "bdif")) list(dif = dif)
              else list(), list(B = B, workers = 4L, seed = seed))),
          warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }), error = function(e) list(.error = conditionMessage(e),
                                       .refusal = inherits(e, "rasch_refusal")))
        list(value = value, warnings = unique(warnings))
      }, args = list(fit = f, dif = dif, B = B, seed = seed, kind = kind,
                     source_dir = .rasch_source_dir),
      libpath = .libPaths(), stdout = log_file, stderr = log_file,
      supervise = TRUE), error = function(e) e)
    if (inherits(process, "error")) {
      progress$close(); unlink(log_file, force = TRUE)
      showNotification(paste("Bootstrap failed:", conditionMessage(process)),
                       type = "error", duration = 10)
      return(invisible(NULL))
    }
    boot_job(list(process = process, progress = progress, log_file = log_file,
                  kind = kind, base_fit = f, base_dif = dif,
                  context = isolate(analysis_context()), B = B, seed = seed))
    invisible(NULL)
  }
  observeEvent(input$boot_run,
    start_bootstrap(fit(), input$boot_B, input$boot_seed, "rasch"))
  observeEvent(input$btl_boot_run,
    start_bootstrap(bfit(), input$btl_boot_B, input$btl_boot_seed, "btl"))
  observeEvent(input$dif_boot_run, {
    f <- fit()
    d <- tryCatch(dif_res(), error = function(e) e)
    if (inherits(d, "error")) {
      showNotification(paste("DIF analysis failed:", conditionMessage(d)),
                       type = "error", duration = 10)
      return()
    }
    start_bootstrap(f, input$dif_boot_B, input$dif_boot_seed, "dif", d)
  })

  observe({
    st <- boot_job()
    if (is.null(st)) return()
    invalidateLater(250, session)
    if (st$process$is_alive()) return()
    result <- tryCatch(st$process$get_result(), error = function(e) e)
    close_boot_job(st); boot_job(NULL)
    current <- if (st$kind %in% c("btl", "bdif"))
      tryCatch(bfit(), error = function(e) NULL) else fit_or_null()
    if (!identical(st$context, isolate(analysis_context())) ||
        !identical(current, st$base_fit) ||
        (st$kind %in% c("dif", "bdif") &&
         !identical(tryCatch(isolate(if (identical(st$kind, "bdif"))
           bdif_res() else dif_res()), error = function(e) NULL),
           st$base_dif))) {
      showNotification("The completed bootstrap was not used because the active analysis changed while it was running.",
                       type = "warning", duration = 8)
      return()
    }
    if (inherits(result, "error")) {
      showNotification(paste("Bootstrap failed:", conditionMessage(result)),
                       type = "error", duration = 10)
      return()
    }
    if (length(result$warnings))
      showNotification(tags$ul(class = "mb-0 ps-3",
        lapply(result$warnings, function(w) tags$li(w))),
        type = "warning", duration = 10)
    if (!is.null(result$value$.error)) {
      showNotification(result$value$.error,
        type = if (isTRUE(result$value$.refusal)) "warning" else "error",
        duration = 10)
      return()
    }
    if (st$kind %in% c("dif", "bdif")) {
      problem <- tryCatch({
        .validate_dif_bootstrap(result$value, st$base_fit, st$base_dif)
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(problem)) {
        showNotification(paste("DIF bootstrap result was not retained:",
                               problem),
                         type = "error", duration = 10)
        return()
      }
      dif_boot_val(list(db = result$value, B = st$B, seed = st$seed,
                        kind = st$kind))
    } else
      boot_val(list(bs = result$value, B = st$B, seed = st$seed,
                    kind = st$kind))
  })
  output$has_dif_boot <- reactive({
    z <- dif_boot_val(); !is.null(z) && !identical(z$kind, "bdif")
  })
  outputOptions(output, "has_dif_boot", suspendWhenHidden = FALSE)
  output$dif_boot_state <- renderUI({
    bv <- dif_boot_val()
    if (is.null(bv) || identical(bv$kind, "bdif")) return(NULL)
    span(class = "badge bg-primary-subtle text-primary-emphasis pb-2 mb-2",
         sprintf(paste0("Sensitivity analysis active: %d of %d used; %d did ",
                        "not converge; %d otherwise failed; seed %d"),
                 bv$db$B_used, bv$B,
                 bv$db$B_nonconverged %||% NA_integer_,
                 bv$db$B_errors %||% NA_integer_, bv$seed))
  })
  output$boot_state <- renderUI({
    bv <- boot_val()
    if (is.null(bv) || identical(bv$kind %||% "rasch", "btl")) return(NULL)
    span(class = "badge bg-primary-subtle text-primary-emphasis pb-2 mb-2",
         sprintf(paste0("Bootstrap null active: %d of %d used; %d did not ",
                        "converge; %d otherwise failed; seed %d"),
                 bv$bs$B_used, bv$B,
                 bv$bs$B_nonconverged %||% NA_integer_,
                 bv$bs$B_errors %||% NA_integer_, bv$seed))
  })
  output$btl_boot_state <- renderUI({
    bv <- boot_val()
    kind <- bv$kind %||% if (!is.null(bv$bs$model) &&
                              identical(bv$bs$model, "paired comparisons"))
      "btl" else "rasch"
    if (is.null(bv) || !identical(kind, "btl")) return(NULL)
    span(class = "badge bg-primary-subtle text-primary-emphasis pb-2 mb-2",
         sprintf(paste0("Bootstrap null active: %d of %d used; %d did not ",
                        "converge; %d otherwise failed; seed %d"),
                 bv$bs$B_used, bv$B,
                 bv$bs$B_nonconverged %||% NA_integer_,
                 bv$bs$B_errors %||% NA_integer_, bv$seed))
  })
  boot_kind <- function(bv) {
    if (is.null(bv)) return(NULL)
    bv$kind %||% if (identical(bv$bs$model %||% "", "paired comparisons"))
      "btl" else "rasch"
  }
  rasch_boot_val <- function() {
    bv <- boot_val()
    if (identical(boot_kind(bv), "rasch")) bv else NULL
  }
  btl_boot_val <- function() {
    bv <- boot_val()
    if (identical(boot_kind(bv), "btl")) bv else NULL
  }
  # fit$items with the bootstrap probabilities alongside, for the table and
  # its CSV; a NULL boot leaves the asymptotic table untouched
  items_with_boot <- function() {
    d <- fit()$items
    if (.has_repeated_residual_units(fit()))
      for (nm in intersect(c("p", "p_adj", "p_bonf", "p_anova",
                             "p_anova_adj", "p_anova_bonf"), names(d)))
        d[[nm]][] <- NA_real_
    bv <- rasch_boot_val()
    if (is.null(bv)) return(d)
    b <- bv$bs$items
    idx <- match(d$item, b$item)
    d$chisq_p_boot <- b$chisq_p_boot[idx]
    d$chisq_p_boot_adj <- b$chisq_p_boot_adj[idx]
    d$n_boot_chisq <- b$n_boot_chisq[idx]
    d$fit_resid_p_boot <- b$fit_resid_p_boot[idx]
    d$fit_resid_p_boot_adj <- b$fit_resid_p_boot_adj[idx]
    d$n_boot_fit_resid <- b$n_boot_fit_resid[idx]
    d
  }

  output$items_vboxes <- renderUI({
    f <- fit()
    bv <- rasch_boot_val()
    bp <- if (is.null(bv)) NULL else bv$bs$items$chisq_p_boot_adj
    inf <- inference_count(if (is.null(bp)) items_with_boot()$p_adj else bp)
    mis <- inf$flagged
    dis <- sum(vapply(f$thresholds_diag, function(d)
      !d$ordered && length(d$thresholds) > 1, TRUE))
    cell_fit <- inherits(f, c("rasch_mfrm", "rasch_efrm"))
    metric_grid(
      metric_tile("metric_cells", if (cell_fit) "Response cells" else "Items",
                  nrow(f$items), icon = "ruler", status = "item"),
      metric_tile("metric_item_misfit",
                  if (is.null(bv)) "Approx. Holm p < .05" else
                    "Bootstrap Holm p < .05",
                  if (is.na(mis)) "Unavailable" else mis, icon = "chisq",
                  status = if (is.na(mis) || is.null(bv)) "neutral"
                    else if (mis > 0) "bad" else "good"),
      metric_tile("metric_disordered", "Disordered thresholds", dis,
                  icon = "disorder", status = if (dis > 0) "bad" else "good"))
  })
  register_code("items_vboxes", function() {
    bv <- rasch_boot_val()
    flag <- if (is.null(bv))
      paste0("  adjusted_p_below_05 = if (any(is.finite(fit$items$p_adj))) ",
             "sum(fit$items$p_adj < .05, na.rm = TRUE) else NA_integer_,")
    else
      paste0("  bootstrap_p_below_05 = if (any(is.finite(",
             "bs$items$chisq_p_boot_adj))) sum(",
             "bs$items$chisq_p_boot_adj < .05, na.rm = TRUE) else NA_integer_,")
    paste(c(
      if (!is.null(bv))
        sprintf("bs <- fit_bootstrap(fit, B = %d, seed = %d)", bv$B, bv$seed),
      "list(",
      "  items = nrow(fit$items),",
      flag,
      "  disordered_thresholds = sum(vapply(fit$thresholds_diag, function(x)",
      "    !x$ordered && length(x$thresholds) > 1, logical(1)))",
      ")"), collapse = "\n")
  })
  output$items_note <- renderUI({
    f <- fit(); d <- items_with_boot()
    dis <- names(which(vapply(f$thresholds_diag, function(x)
      !x$ordered && length(x$thresholds) > 1, TRUE)))
    unit <- if (inherits(f, c("rasch_mfrm", "rasch_efrm")))
      "response cells" else "items"
    bv <- rasch_boot_val()
    chi <- if (is.null(bv)) {
      z <- inference_count(d$p_adj)
      if (!z$tested) "approximate asymptotic Holm inference unavailable"
      else paste(z$text, "with approximate asymptotic Holm p < .05")
    } else {
      z <- inference_count(bv$bs$items$chisq_p_boot_adj)
      if (!z$tested) "bootstrap chi-square inference unavailable"
      else sprintf("%s with Holm-adjusted bootstrap chi-square p < .05 (B = %d)",
                   z$text, bv$B)
    }
    sprintf("Note. %d of %d %s beyond |fit residual| 2.5; %s; disordered thresholds: %s.",
            sum(abs(d$fit_resid) > 2.5, na.rm = TRUE), nrow(d),
            unit, chi,
            if (length(dis)) paste(dis, collapse = ", ") else "none")
  })
  register_table("items_tbl", items_with_boot, function() {
    bv <- rasch_boot_val()
    # with a bootstrap null active the asymptotic probabilities give way to
    # the calibrated ones, in the display and in the CSV alike
    extra <- if (is.null(bv)) NULL
             else c("chisq_p_boot_adj", "n_boot_chisq",
                    "fit_resid_p_boot_adj", "n_boot_fit_resid")
    d <- curate(items_with_boot(), "items", full = isTRUE(input$items_full),
                extra = extra)
    if (!is.null(bv) && !isTRUE(input$items_full))
      d$p_adj <- NULL
    operative_p <- if (is.null(bv)) character(0) else
      c("chisq_p_boot_adj", "fit_resid_p_boot_adj")
    dt <- num_dt(d, page_len = 25, selection = "single",
                 p_bold = operative_p)
    # fit residual, infit and outfit are flagged by num_dt; the operative
    # chi-square probability turns red at .05 as well (no single flag column)
    style_lo_red(dt, d, "chisq_p_boot_adj", .05)
  }, code = function() {
    bv <- rasch_boot_val()
    if (is.null(bv)) return("fit$items")
    paste(sprintf("bs <- fit_bootstrap(fit, B = %d, seed = %d)", bv$B, bv$seed),
          "b <- bs$items[c(\"item\", \"chisq_p_boot\",",
          "  \"chisq_p_boot_adj\", \"n_boot_chisq\",",
          "  \"fit_resid_p_boot\", \"fit_resid_p_boot_adj\", \"n_boot_fit_resid\")]",
          "j <- match(fit$items$item, b$item)",
          "cbind(fit$items, b[j, setdiff(names(b), \"item\"), drop = FALSE])",
          sep = "\n")
  })

  # per-class-interval breakdown of the selected item's chi-square
  chisq_res <- reactive(chisq_detail(fit(), sel_item()))
  output$chisq_caption <- renderUI({
    cd <- chisq_res()
    p(class = "small mb-2",
      tags$b(cd$item),
      sprintf(" (location %.3f): total chi-square ", cd$location),
      tags$b(sprintf("%.3f", cd$chisq)),
      sprintf(paste0(" on %d df, p = %s; whole-sample mean = %.3f. ",
                     "Intervals with fewer than 2 responders carry no ",
                     "chi-square contribution."),
              cd$df, fmt_p(cd$p), cd$ave))
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
    filename = function() paste0("rasch_chisq_intervals_", sel_item(), ".csv"),
    content = function(file)
      write_csv_plain(chisq_res()$intervals, file))
  output$chisq_cat_csv <- downloadHandler(
    filename = function() paste0("rasch_chisq_categories_", sel_item(), ".csv"),
    content = function(file)
      write_csv_plain(chisq_res()$categories, file))
  register_code("chisq", function()
    sprintf('chisq_detail(fit, %s)', qstr(sel_item())))

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
  }, code = function() "fit$est$components")

  # explorer display settings (inline on the tab-strip controls row):
  # class intervals and scale
  # range, resolved with fallbacks; the code footers carry the displayed grid
  # explicitly so that automatic and preset ranges reproduce the app plot
  ex_ng <- reactive({
    ng <- input$ex_ng
    if (is.null(ng) || is.na(ng)) fit()$n_groups else as.integer(ng)
  })
  ex_rng <- reactive({
    resolve_axis_range("ex_axis", c(-5, 5), c(-8, 8),
                       function() fitted_scale_range(c(-5, 5), pad = .75))
  })
  ex_grid <- reactive(seq(ex_rng()[1], ex_rng()[2], 0.05))
  ex_code_args <- reactive(paste0(c(
    if (ex_ng() != fit()$n_groups) sprintf(", n_groups = %d", ex_ng()),
    sprintf(", grid = seq(%g, %g, 0.05)", ex_rng()[1], ex_rng()[2])),
    collapse = ""))
  # second code-footer line pointing at the matching all-items batch export
  ex_batch_line <- function(what) {
    grid <- ex_rng()
    sprintf(paste0('\n# all items: save_item_plots(fit, %s, ',
                   '%s, n_groups = %d, grid = seq(%g, %g, 0.05), ',
                   'observed = %s)'),
            qstr(what), qstr(paste0(what, "_all_items.pdf")), ex_ng(),
            grid[1], grid[2], if (isTRUE(input$show_obs)) "TRUE" else "FALSE")
  }
  register_plot("icc",  function()
    plot_icc(fit(), icc_items(), n_groups = ex_ng(), grid = ex_grid(),
             observed = isTRUE(input$show_obs)),
    code = function() paste0(sprintf('plot_icc(fit, %s, observed = %s%s)',
                                     qvec(icc_items()),
                                     isTRUE(input$show_obs), ex_code_args()),
                             ex_batch_line("icc")))
  register_plot("ccc",  function()
    plot_ccc(fit(), sel_item(), observed = isTRUE(input$show_obs),
             n_groups = ex_ng(), grid = ex_grid()),
    code = function() paste0(sprintf('plot_ccc(fit, %s, observed = %s%s)',
                                     qstr(sel_item()), isTRUE(input$show_obs),
                                     ex_code_args()),
                             ex_batch_line("ccc")))
  register_plot("tpc",  function()
    plot_threshold_prob(fit(), sel_item(), observed = isTRUE(input$show_obs),
                        n_groups = ex_ng(), grid = ex_grid()),
    code = function() paste0(
      sprintf('plot_threshold_prob(fit, %s, observed = %s%s)',
              qstr(sel_item()), isTRUE(input$show_obs), ex_code_args()),
      ex_batch_line("tpc")))
  register_plot("cfreq", function() plot_catfreq(fit(), sel_item()),
                code = function() paste0(
                  sprintf('plot_catfreq(fit, %s)', qstr(sel_item())),
                  ex_batch_line("cfreq")))
  # all-items batch downloads for the active explorer tab, honouring the
  # class-interval, range, and observed-points controls
  ex_what <- reactive(switch(input$items_nav %||% "ICC",
    ICC = "icc", Categories = "ccc", Thresholds = "tpc",
    Frequencies = "cfreq", "icc"))
  register_batch_download("items_all",
    base = function() paste0(ex_what(), "_all_items"),
    content = function(file)
      withProgress(message = "Drawing every item…", value = 0.4,
        save_item_plots(fit(), ex_what(), file, n_groups = ex_ng(),
                        grid = ex_grid(),
                        observed = isTRUE(input$show_obs))))
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
  }, code = function() "distractor_analysis(fit)")
  # the plotted item: the selected one if it is multiple-choice, else the
  # first MC item (the code disclosure mirrors this resolution)
  distractor_item <- reactive({
    f <- fit()
    req(!is.null(f$mc))
    if (sel_item() %in% colnames(f$mc$raw)) sel_item() else
      colnames(f$mc$raw)[1]
  })
  register_plot("distractor_plot", function() {
    f <- fit()
    validate(need(!is.null(f$mc),
                  "Provide a multiple-choice key (CSV: item,key) to see option curves."))
    plot_distractors(f, distractor_item())
  }, code = function()
    sprintf('plot_distractors(fit, %s)', qstr(distractor_item())))

  # polytomous option-scoring proposal (Andrich & Styles 2011)
  rescore_res <- reactiveVal(NULL)
  observeEvent(input$rescore_go, {
    f <- fit()
    if (is.null(f$mc)) {
      showNotification("Provide a multiple-choice key first.", type = "warning")
      return()
    }
    run_min_n <- max(2, input$rescore_min_n %||% 20)
    run_z <- max(0.1, input$rescore_z %||% 1.96)
    res <- tryCatch(distractor_rescore(f, min_n = run_min_n, z = run_z),
                    error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Rescore proposal failed:", conditionMessage(res)),
                       type = "error")
      return()
    }
    res$run_min_n <- run_min_n
    res$run_z <- run_z
    rescore_res(res)
    score_key <- app_factor_keys(res$option_scores[c("item", "option")])
    evidence_key <- app_factor_keys(res$evidence[c("item", "option")])
    n_cred <- sum(res$option_scores$score > 0 &
                  !res$evidence$keyed[match(score_key, evidence_key)])
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
  register_code("rescore_tbl", function()
    {
      r <- rescore_res(); req(!is.null(r))
      sprintf('distractor_rescore(fit, min_n = %s, z = %s)$evidence',
              r$run_min_n, r$run_z)
    })
  output$dl_rescore <- downloadHandler(
    filename = function() "option_scores.csv",
    content = function(file) {
      res <- rescore_res()
      if (is.null(res)) stop("run the proposal first")
      write_csv_plain(res$option_scores, file)
    })

  # -------------------------------------------------------------- persons --
  persons_with_boot <- function() {
    d <- fit()$person
    bv <- rasch_boot_val()
    if (is.null(bv) || is.null(bv$bs$persons)) return(d)
    b <- bv$bs$persons
    idx <- if (nrow(d) == nrow(b)) seq_len(nrow(d)) else match(d$id, b$id)
    add <- grep("(_p_boot(_adj)?$|^n_boot_)", names(b), value = TRUE)
    for (nm in add) d[[nm]] <- b[[nm]][idx]
    d
  }
  output$persons_vboxes <- renderUI({
    f <- fit(); d <- persons_with_boot(); bv <- rasch_boot_val()
    has_person_boot <- !is.null(bv) && !is.null(bv$bs$persons)
    pp <- if (has_person_boot) d$fit_resid_p_boot_adj else NULL
    n_testable <- if (is.null(pp)) NA_integer_ else sum(is.finite(pp))
    mis <- if (!has_person_boot)
      sum(abs(d$fit_resid) > 2.5, na.rm = TRUE) else
      if (n_testable) sum(pp < .05, na.rm = TRUE) else NA_integer_
    metric_grid(
      metric_tile("metric_persons", "Persons", nrow(d), icon = "distribution",
                  status = "person"),
      metric_tile("metric_extreme", "Extreme scores",
                  sum(d$extreme, na.rm = TRUE), icon = "range",
                  status = "person"),
      metric_tile("metric_person_misfit",
                  if (!has_person_boot) "|Fit residual| > 2.5" else
                    "Bootstrap adjusted p < .05",
                  if (is.na(mis)) "Unavailable" else mis,
                  icon = "outlier", status = if (is.na(mis)) "neutral" else
                    if (mis > 0) "bad" else "good"))
  })
  register_code("persons_vboxes", function() {
    bv <- rasch_boot_val()
    has_person_boot <- !is.null(bv) && !is.null(bv$bs$persons)
    paste(c(
      if (has_person_boot)
        sprintf("bs <- fit_bootstrap(fit, B = %d, seed = %d)", bv$B, bv$seed),
      "list(",
      "  persons = nrow(fit$person),",
      "  extreme_scores = sum(fit$person$extreme, na.rm = TRUE),",
      if (!has_person_boot)
        "  misfitting_persons = sum(abs(fit$person$fit_resid) > 2.5, na.rm = TRUE)" else
        paste0("  misfitting_persons = if (any(is.finite(",
               "bs$persons$fit_resid_p_boot_adj))) sum(",
               "bs$persons$fit_resid_p_boot_adj < .05, na.rm = TRUE) else NA_integer_"),
      ")"), collapse = "\n")
  })
  register_table("person_tbl", persons_with_boot, function() {
    fac <- names(fit()$factors)
    bv <- rasch_boot_val()
    extra <- c(fac, if (!is.null(bv) && !is.null(bv$bs$persons))
      c("fit_resid_p_boot_adj", "n_boot_fit_resid"))
    d <- curate(persons_with_boot(), "person",
                full = isTRUE(input$persons_full), extra = extra)
    dt <- datatable(d, rownames = FALSE, filter = "top", selection = "single",
                    style = "bootstrap5",
                    class = "table-sm compact hover order-column",
                    options = list(pageLength = 15, scrollX = TRUE, dom = "tip")) |>
      formatRound(names(d)[vapply(d, is.numeric, TRUE) &
                           !names(d) %in% c("raw", "max_raw", "n_items", "class_interval")], 3)
    # the same fit flags as every other table (fit residual, and the infit /
    # outfit mean squares revealed by the detailed-columns switch)
    dt <- flag_fit_cols(dt, names(d))
    if (!is.null(bv) && !is.null(bv$bs$persons))
      dt <- style_lo_red(dt, d, "fit_resid_p_boot_adj", .05)
    dt
  }, code = function() {
    bv <- rasch_boot_val()
    if (is.null(bv) || is.null(bv$bs$persons)) return("fit$person")
    paste(sprintf("bs <- fit_bootstrap(fit, B = %d, seed = %d)", bv$B, bv$seed),
          "cols <- grep('(_p_boot(_adj)?$|^n_boot_)', names(bs$persons), value = TRUE)",
          "cbind(fit$person, bs$persons[cols])",
          sep = "\n")
  })

  observeEvent(input$person_weights_go, {
    f <- fit()
    upload <- input$person_weights_file
    if (is.null(upload)) {
      showNotification("Choose a weights CSV first.", type = "warning")
      return()
    }
    parsed <- tryCatch({
      d <- read.csv(upload$datapath, check.names = FALSE,
                    stringsAsFactors = FALSE)
      if (anyDuplicated(names(d)))
        stop("the weights CSV has duplicate column names")
      by <- input$person_weight_level %||% "item"
      sets <- NULL
      if (identical(by, "item")) {
        if (!all(c("item", "weight") %in% names(d)))
          stop("item weights need columns item and weight")
        if (!nrow(d) || anyNA(d$item) || any(!nzchar(as.character(d$item))) ||
            anyDuplicated(as.character(d$item)))
          stop("item weights need one non-empty row per item")
        weights <- stats::setNames(as.numeric(d$weight), as.character(d$item))
      } else if ("item" %in% names(d)) {
        if (!all(c("item", "set", "weight") %in% names(d)))
          stop("set weights need columns item, set and weight")
        if (!nrow(d) || anyNA(d[, c("item", "set", "weight")]) ||
            any(!nzchar(as.character(d$item))) ||
            any(!nzchar(as.character(d$set))) ||
            anyDuplicated(as.character(d$item)))
          stop("set weights need one complete row per item")
        split_w <- split(as.numeric(d$weight), as.character(d$set))
        if (any(vapply(split_w, function(z)
          length(unique(z)) != 1L, logical(1))))
          stop("every row in an item set must carry the same set weight")
        weights <- vapply(split_w, `[`, numeric(1), 1L)
        sets <- stats::setNames(as.character(d$set), as.character(d$item))
      } else {
        if (!inherits(f, "rasch_efrm"))
          stop("set weights for this model need columns item, set and weight")
        if (!all(c("set", "weight") %in% names(d)) || !nrow(d) ||
            anyNA(d[, c("set", "weight")]) ||
            any(!nzchar(as.character(d$set))) ||
            anyDuplicated(as.character(d$set)))
          stop("Extended Frames set weights need one row per set with columns set and weight")
        weights <- stats::setNames(as.numeric(d$weight), as.character(d$set))
      }
      if (anyNA(weights)) stop("weights must be numeric")
      tab <- weighted_person_estimates(f, weights, by = by, sets = sets)
      state <- list(table = tab, weights = weights, by = by, sets = sets,
                    filename = upload$name,
                    fit_signature = .fit_boot_signature(f))
      state$result_signature <- .fit_boot_md5(state)
      state
    }, error = function(e) e)
    if (inherits(parsed, "error")) {
      showNotification(paste("Weighted estimates failed:",
                             conditionMessage(parsed)),
                       type = "error", duration = 10)
      return()
    }
    person_weight_state(parsed)
    showNotification("Weighted person estimates calculated.",
                     type = "message", duration = 5)
  })
  named_vector_code <- function(x) paste0(
    "c(", paste(sprintf("%s = %s", qstr(names(x)),
                         format(unname(x), digits = 15, trim = TRUE,
                                scientific = FALSE)), collapse = ", "), ")")
  person_weight_code <- function() {
    z <- person_weight_state(); req(!is.null(z))
    lines <- paste0("weights <- ", named_vector_code(z$weights))
    set_arg <- ""
    if (!is.null(z$sets)) {
      lines <- c(lines, paste0("item_sets <- ",
        "c(", paste(sprintf("%s = %s", qstr(names(z$sets)),
                             qstr(unname(z$sets))), collapse = ", "), ")"))
      set_arg <- ", sets = item_sets"
    }
    c(lines, sprintf("weighted_person_estimates(fit, weights, by = %s%s)",
                     qstr(z$by), set_arg)) |>
      paste(collapse = "\n")
  }
  register_table("person_weight_tbl", function() {
    z <- person_weight_state(); req(!is.null(z)); z$table
  }, function() {
    z <- person_weight_state()
    validate(need(!is.null(z),
                  "Upload the external weights and calculate the estimates."))
    d <- z$table
    datatable(d, rownames = FALSE, filter = "top", style = "bootstrap5",
      class = "table-sm compact hover order-column",
      options = list(pageLength = 15, scrollX = TRUE, dom = "tip")) |>
      formatRound(names(d)[vapply(d, is.numeric, TRUE) &
        !names(d) %in% c("raw", "max_raw", "n_items")], 3)
  }, code = person_weight_code)
  # the person drawn on the kidmap: the selected table row, defaulting to 1
  sel_person <- reactive({
    .app_selected_row(input$person_tbl_rows_selected,
                      nrow(fit()$person))
  })
  # confidence level for the kidmap band (header select; default 95%);
  # the code footers mention `level` only when it differs from the default
  kid_level <- reactive({
    lv <- suppressWarnings(as.numeric(input$kid_level %||% "0.95"))
    if (!isTRUE(is.finite(lv)) || lv <= 0 || lv >= 1) 0.95 else lv
  })
  kid_level_arg <- reactive(
    if (isTRUE(all.equal(kid_level(), 0.95))) ""
    else sprintf(", level = %g", kid_level()))
  register_plot("kidmap", function()
    plot_kidmap(fit(), person = sel_person(), level = kid_level()),
    code = function() paste0(
      sprintf('plot_kidmap(fit, person = %d%s)', sel_person(),
              kid_level_arg()),
      sprintf('\n# all persons: save_person_plots(fit, "kidmaps.pdf"%s)',
              kid_level_arg())))
  # batch kidmaps: multi-page PDF or ZIP of PNGs, chosen by the extension
  # (the download temp file carries the filename's extension)
  register_batch_download("kidmap_all",
    base = function() "kidmaps",
    content = function(file)
      withProgress(message = "Drawing a kidmap for every person…",
                   value = 0.4,
                   save_person_plots(fit(), file, level = kid_level())))
  register_plot("rdist_p", function() plot_resid_dist(fit(), "persons"),
                code = function() 'plot_resid_dist(fit, "persons")')
  register_plot("pfit",  function() plot_person_fit(fit()),
                code = function() "plot_person_fit(fit)")
  # hover identification for the person-fit scatter: nearPoints() against
  # exactly the rows plot_person_fit() draws (theta and fit_resid both
  # non-NA -- see R/plots.R), so the point under the cursor is always the
  # one actually plotted there
  register_hover_tip("pfit", function() {
    p <- fit()$person
    ok <- !is.na(p$theta) & !is.na(p$fit_resid)
    data.frame(id = p$id[ok], theta = p$theta[ok], fit_resid = p$fit_resid[ok])
  }, "theta", "fit_resid", function(np)
    sprintf("%s · location %.2f · fit residual %.2f",
            np$id[1], np$theta[1], np$fit_resid[1]))

  # ------------------------------------------------------- targeting plots --
  # sidebar-controlled bins and scale range shared by both targeting plots
  tg_bins <- reactive({
    b <- input$tg_bins
    if (is.null(b) || is.na(b)) 35L else as.integer(b)
  })
  tg_rng <- reactive({
    resolve_axis_range("tg_axis", c(-5, 5), c(-8, 8),
                       function() fitted_scale_range(c(-5, 5), pad = .5))
  })
  tg_code <- function(fun, information = FALSE) function() {
    r <- tg_rng() %||% c(-5, 5)
    info_arg <- if (isTRUE(information) && isTRUE(input$tg_information))
      ", information = TRUE" else ""
    # deparse() rather than a quoted %s: a level or set name carrying a
    # quote or backslash must still yield runnable code
    sel <- paste0(
      if (!is.null(tg_group())) paste0(", group = ", deparse(tg_group())) else "",
      if (!is.null(tg_items())) paste0(", items = ", deparse(tg_items())) else "")
    sprintf("%s(fit, bins = %d, xlim = c(%g, %g)%s%s)",
            fun, tg_bins(), r[1], r[2], info_arg, sel)
  }
  # the person groups a fit actually carries, and its item sets if any
  observeEvent(fit(), {
    f <- fit()
    # every choice is the qualified 'factor: level' address, which is what
    # plot_pimap() disambiguates on -- two factors may share a level name,
    # and a bare level would silently take the first
    lev <- if (!is.null(f$factors) && ncol(f$factors))
      sort(unique(unlist(lapply(names(f$factors), function(nm)
        sprintf("%s: %s", nm,
                sort(as.character(unique(f$factors[[nm]])))))))) else character(0)
    updateSelectInput(session, "tg_group",
                      choices = c("All persons" = "", stats::setNames(lev, lev)),
                      selected = "")
    sets <- if (!is.null(f$set_of)) sort(unique(as.character(f$set_of)))
            else character(0)
    updateSelectInput(session, "tg_items",
                      choices = c("All items" = "",
                                  stats::setNames(sets, sets)),
                      selected = "")
  }, ignoreNULL = TRUE)
  tg_group <- reactive({ z <- input$tg_group %||% ""; if (nzchar(z)) z else NULL })
  tg_items <- reactive({ z <- input$tg_items %||% ""; if (nzchar(z)) z else NULL })
  register_plot("pim_p", function()
    plot_pimap(fit(), bins = tg_bins(), xlim = tg_rng(),
               information = isTRUE(input$tg_information),
               group = tg_group(), items = tg_items()),
    code = tg_code("plot_pimap", information = TRUE))
  wright_package_on <- reactive(
    identical(input$wright_renderer, "wrightmap"))
  wright_person_arg <- reactive({
    z <- input$wright_person_panels %||% ""
    if (nzchar(z)) z else NULL
  })
  wright_item_map_data <- reactive({
    d <- if (!is.null(input$wright_item_map))
      tryCatch(read.csv(input$wright_item_map$datapath,
                        check.names = FALSE, stringsAsFactors = FALSE),
               error = function(e) NULL)
    else restored_project_resources()[["wright_item_map"]]
    validate(need(!is.null(d) && is.data.frame(d) &&
                    all(c("item", "panel") %in% names(d)),
                  "The item panel map needs columns named item and panel."))
    validate(need(!anyDuplicated(names(d)),
                  "The item panel map cannot have duplicate column names."))
    validate(need(!anyNA(d$item) && !anyNA(d$panel) &&
                    all(nzchar(trimws(as.character(d$item)))) &&
                    all(nzchar(trimws(as.character(d$panel)))),
                  "The item panel map cannot contain blank values."))
    d$item <- trimws(as.character(d$item))
    d$panel <- trimws(as.character(d$panel))
    validate(need(!anyDuplicated(d$item),
                  paste0("The item panel map names item(s) more than once: ",
                         paste(unique(d$item[duplicated(d$item)]),
                               collapse = ", "), ".")))
    d[, c("item", "panel"), drop = FALSE]
  })
  wright_item_arg <- reactive({
    z <- input$wright_item_panels %||% ""
    if (!nzchar(z)) return(NULL)
    if (identical(z, "uploaded")) {
      d <- wright_item_map_data()
      return(stats::setNames(d$panel, d$item))
    }
    if (identical(z, "sets_groups")) c("sets", "groups") else z
  })
  draw_wright <- function() {
    if (!wright_package_on())
      return(plot_wright(fit(), bins = tg_bins(), xlim = tg_rng()))
    validate(need(.wrightmap_available,
      paste("WrightMap is not installed. Install it with",
            "install.packages(\"WrightMap\"), restart the app, and try again.")))
    r <- tg_rng()
    side <- if (identical(input$wright_person_style, "density"))
      WrightMap::personDens else WrightMap::personHist
    wright_map(fit(), type = input$wright_type %||% "thresholds",
               person_panels = wright_person_arg(),
               item_panels = wright_item_arg(), person.side = side,
               min.l = r[1], max.l = r[2], breaks = tg_bins(),
               main.title = NULL)
  }
  wright_code <- function() {
    if (!wright_package_on()) return(tg_code("plot_wright")())
    r <- tg_rng()
    args <- c("fit",
      sprintf('type = "%s"', input$wright_type %||% "thresholds"),
      if (!is.null(wright_person_arg()))
        paste0("person_panels = ", qstr(wright_person_arg())),
      if (!is.null(wright_item_arg())) {
        if (identical(input$wright_item_panels, "uploaded"))
          "item_panels = setNames(item_panel_map$panel, item_panel_map$item)"
        else if (length(wright_item_arg()) == 1L)
          paste0("item_panels = ", qstr(wright_item_arg()))
        else 'item_panels = c("sets", "groups")'
      },
      if (identical(input$wright_person_style, "density"))
        "person.side = WrightMap::personDens",
      sprintf("min.l = %g", r[1]), sprintf("max.l = %g", r[2]),
      sprintf("breaks = %d", tg_bins()), "main.title = NULL")
    call <- paste0("wright_map(", paste(args, collapse = ", "), ")")
    if (identical(input$wright_item_panels, "uploaded")) {
      source <- if (!is.null(input$wright_item_map)) paste0(
        "item_panel_map <- read.csv(",
        qstr(input$wright_item_map$name %||% "item_panels.csv"),
        ", check.names = FALSE, stringsAsFactors = FALSE)") else
          'item_panel_map <- project$resources[["wright_item_map"]]'
      paste(source,
            "item_panel_map$item <- trimws(as.character(item_panel_map$item))",
            "item_panel_map$panel <- trimws(as.character(item_panel_map$panel))",
            call, sep = "\n")
    } else call
  }
  register_plot("wright", draw_wright, h = 7.5, code = wright_code)

  # ---------------------------------------------- test & item map plots --
  # thrmap and imap render on the Items page; tcc/tif/guttman on Test
  register_plot("thrmap", function() plot_threshold_map(fit()), h = 7,
                code = function() "plot_threshold_map(fit)")
  register_plot("imap",   function() plot_item_map(fit()),
                code = function() "plot_item_map(fit)")
  # hover identification for the item fit map: nearPoints() against
  # fit()$items (plot_item_map() draws every item -- no filtering)
  register_hover_tip("imap", function() fit()$items, "location", "fit_resid",
    function(np) sprintf("%s · location %.3f · fit residual %.2f",
                         np$item[1], np$location[1], np$fit_resid[1]))
  register_plot("rdist_i", function() plot_resid_dist(fit(), "items"),
                code = function() 'plot_resid_dist(fit, "items")')
  # Test-page scale range. The code footer carries the displayed grid
  # explicitly so that automatic and preset ranges reproduce the app plot.
  ts_rng <- reactive({
    resolve_axis_range("ts_axis", c(-6, 6), c(-8, 8),
                       function() fitted_scale_range(c(-6, 6), pad = 1))
  })
  ts_grid <- reactive(seq(ts_rng()[1], ts_rng()[2], 0.05))
  ts_code_arg <- reactive(
    sprintf(", grid = seq(%g, %g, 0.05)", ts_rng()[1], ts_rng()[2]))
  register_plot("tcc",    function() plot_tcc(fit(), grid = ts_grid()),
                code = function() paste0("plot_tcc(fit", ts_code_arg(), ")"))
  register_plot("tif",    function() plot_tif(fit(), grid = ts_grid()),
                code = function() paste0("plot_tif(fit", ts_code_arg(), ")"))
  register_plot("guttman", function() {
    validate(need(alpha_design_applicable(fit()),
                  paste("The whole-item scalogram is not defined when an",
                        "item has several frame or facet response cells.")))
    validate(need(all(fit()$m == 1L),
                  "The whole-item Guttman scalogram is available for dichotomous items only."))
    plot_guttman(fit())
  }, h = 7,
                code = function() "plot_guttman(fit)")
  # hover identification for the Guttman scalogram: neither axis is labelled
  # (persons are thinned to at most 80 rows and never named; items are
  # labelled but only along the x-axis), so a cell's person and item are
  # otherwise unreadable. Reproduces plot_guttman()'s own thinning of
  # guttman_table(fit)$matrix (max_persons = 80, its default and the app's
  # only call), then maps the hovered data-space cell to (row, col) exactly
  # as image() lays the matrix out: x is the item column left-to-right, y is
  # the person row flipped top-to-bottom -- see R/guttman.R.
  guttman_res <- reactive({
    validate(need(alpha_design_applicable(fit()),
                  paste("The whole-item scalogram is not defined when an",
                        "item has several frame or facet response cells.")))
    validate(need(all(fit()$m == 1L),
                  "The whole-item Guttman scalogram is available for dichotomous items only."))
    g <- guttman_table(fit())
    G <- g$matrix
    N <- nrow(G)
    if (N > 80) G <- G[round(seq(1, N, length.out = 80)), , drop = FALSE]
    G
  })
  output$guttman_tip <- renderUI({
    hov <- input$guttman_hover
    req(hov)
    G <- guttman_res()
    req(!is.null(G))
    L <- ncol(G); rN <- nrow(G)
    item_idx <- round(hov$x); person_row <- round(hov$y)
    req(item_idx >= 1, item_idx <= L, person_row >= 1, person_row <= rN)
    person_idx <- rN + 1L - person_row
    val <- G[person_idx, item_idx]
    div(class = "rasch-hover-tip",
        style = sprintf("left:%.0fpx; top:%.0fpx;",
                        hov$coords_css$x, hov$coords_css$y - 12),
        sprintf("%s · %s · score %s", rownames(G)[person_idx],
                colnames(G)[item_idx], if (is.na(val)) "missing" else val))
  })

  # ------------------------------------------------------------------ DIF --
  dif_alpha <- reactive(clamp01(input$dif_alpha, 0.05))
  # the Model toggle (main effects vs interactions) is immaterial with a
  # single nominated factor, where the model is always the one-way ANOVA
  dif_multi <- reactive({
    f <- fit()
    available <- setdiff(names(f$factors), f$frame_group %||% character(0))
    length(available) > 1L
  })
  output$dif_multifactor <- reactive(dif_multi())
  outputOptions(output, "dif_multifactor", suspendWhenHidden = FALSE)
  output$dif_refit_available <- reactive(
    !inherits(fit(), c("rasch_efrm", "rasch_mfrm")))
  outputOptions(output, "dif_refit_available", suspendWhenHidden = FALSE)
  output$dif_followup_available <- reactive(!inherits(fit(), "rasch_efrm"))
  outputOptions(output, "dif_followup_available", suspendWhenHidden = FALSE)
  output$dif_refit_note <- renderUI({
    f <- fit()
    txt <- if (inherits(f, "rasch_mfrm"))
      "Resolve Multiple Ratings DIF in the long-format data and refit the facet model."
    else
      "Resolve Extended Frames DIF in the source data; an ordinary split would discard the fitted frame units."
    p(class = "text-muted small mt-3", txt)
  })
  # one merged reactive: one factor -> one-way ANOVA (one row per item); several
  # factors -> joint model, main effects by default or factor-by-factor
  # interactions when requested (effects is ignored with a single factor)
  dif_res <- reactive({
    f <- fit(); req(!is.null(f$factors))
    soft(dif_anova(f, effects = input$dif_effects %||% "main",
                   p_adjust = "holm", alpha = dif_alpha()))
  })
  observeEvent(list(input$dif_effects, input$dif_alpha), {
    if (!isTRUE(restoring_project())) dif_boot_val(NULL)
  }, ignoreInit = TRUE)
  # code footer: omit the effects argument when there is only one factor
  dif_effects_arg <- function()
    if (dif_multi()) sprintf('effects = %s, ', qstr(input$dif_effects %||% "main"))
    else ""
  register_table("dif_tbl", function() dif_res()$summary, function() {
    d <- curate(dif_res()$summary, "dif_fact", full = isTRUE(input$dif_full))
    if ("superseded" %in% names(d))
      d$superseded <- ifelse(d$superseded, "(superseded)", "")
    dt <- num_dt(d, selection = "single")
    # no boolean DIF flag columns: the adjusted probability turns red when it
    # crosses alpha (uniform = factor effect, non-uniform = factor x interval)
    dt <- style_lo_red(dt, d, "p_uniform_adj", dif_alpha())
    style_lo_red(dt, d, "p_nonuniform_adj", dif_alpha())
  }, code = function()
    sprintf('dif_anova(fit, %sp_adjust = %s, alpha = %s)$summary',
            dif_effects_arg(), qstr("holm"), dif_alpha()))
  # full per-item ANOVA table: every model term, computed lazily when its
  # disclosure is first switched on (the DT renders only once visible)
  register_table("dif_full_tbl", function() dif_res()$terms, function() {
    d <- dif_res()$terms
    d$significant <- NULL                     # red adjusted p replaces the flag
    d$superseded <- ifelse(d$superseded, "(superseded)", "")
    style_lo_red(num_dt(d), d, "p_adj", dif_alpha())
  }, code = function()
    sprintf('dif_anova(fit, %sp_adjust = %s, alpha = %s)$terms',
            dif_effects_arg(), qstr("holm"), dif_alpha()))
  register_table("dif_boot_tbl", function() {
    bv <- dif_boot_val(); req(!is.null(bv)); bv$db$summary
  }, function() {
    bv <- dif_boot_val(); req(!is.null(bv))
    d <- bv$db$summary
    if (!isTRUE(input$dif_boot_full))
      d <- d[, intersect(c(
        "item", "term", "p_uniform_adj", "p_uniform_boot_adj",
        "p_nonuniform_adj", "p_nonuniform_boot_adj"), names(d)), drop = FALSE]
    dt <- num_dt(d)
    dt <- style_lo_red(dt, d, "p_uniform_boot_adj", dif_alpha())
    style_lo_red(dt, d, "p_nonuniform_boot_adj", dif_alpha())
  }, code = function() {
    bv <- dif_boot_val(); req(!is.null(bv))
    sprintf(paste0(
      "d <- dif_anova(fit, %sp_adjust = %s, alpha = %s)\n",
      "dif_bootstrap(fit, d, B = %d, workers = 4, seed = %d)$summary"),
      dif_effects_arg(), qstr("holm"), dif_alpha(), bv$B, bv$seed)
  })
  output$dif_note <- renderUI({
    r <- dif_res(); d <- r$summary
    pp <- c(d$p_uniform_adj, d$p_nonuniform_adj)
    tested <- sum(is.finite(pp))
    unavailable <- length(pp) - tested
    status <- if (tested) sprintf(
      "%d of %d available DIF effects significant after adjustment%s",
      sum(pp[is.finite(pp)] < dif_alpha()), tested,
      if (unavailable) sprintf(" (%d unavailable)", unavailable) else "")
    else sprintf("no DIF effect supports inference (%d unavailable)",
                 unavailable)
    base <- sprintf("Note. %s. Class intervals: %s (from the smallest cell).",
                    status, if (is.null(r$n_groups)) "NA" else r$n_groups)
    if (identical(r$effects, "factorial")) {
      sup <- sum(d$superseded, na.rm = TRUE)
      if (sup)
        base <- paste0(base,
          sprintf(" %d main-effect term(s) superseded by an interaction.", sup))
    }
    within <- r$within[!is.na(r$within)]
    if (length(within))
      base <- paste0(base,
        sprintf(" Within-subject factor(s) tested by repeated-measures ANOVA: %s.",
                paste(within, collapse = ", ")))
    base
  })
  # the items of the DIF summary in rendered row order (curate only drops
  # columns, so the order is preserved). The selected row drives the group-ICC
  # and the pairwise-comparisons panel; nothing selected defaults to the top row.
  dif_tbl_items <- reactive(dif_res()$summary$item)
  # the item and term of the currently selected summary row (top row by
  # default); the exact factor variables stored by dif_anova() feed dif_size
  # and the group ICC, without parsing their display label
  dif_sel_row <- reactive({
    r <- input$dif_tbl_rows_selected
    its <- dif_tbl_items()
    if (length(r) && !is.na(r[1]) && r[1] >= 1 && r[1] <= length(its))
      r[1] else 1L
  })
  dif_sel_item <- reactive({
    its <- dif_tbl_items(); req(length(its) >= 1); its[dif_sel_row()]
  })
  dif_sel_term <- reactive({
    tm <- dif_res()$summary$term; req(length(tm) >= 1); tm[dif_sel_row()]
  })
  dif_sel_vars <- reactive({
    r <- dif_res(); i <- dif_sel_row()
    if (!is.null(r$summary_factors) && length(r$summary_factors) >= i)
      r$summary_factors[[i]]
    else strsplit(dif_sel_term(), ":", fixed = TRUE)[[1]]
  })
  output$dif_selected_interaction <- reactive(length(dif_sel_vars()) > 1L)
  outputOptions(output, "dif_selected_interaction", suspendWhenHidden = FALSE)
  # the group-ICC uses the selected term's factor(s); plot_icc accepts several
  # factor names, so an interaction row overlays the factor-combination cells
  register_plot("dif_icc", function() {
    f <- fit()
    req(!is.null(f$factors))
    if (inherits(f, "rasch_efrm")) {
      req(dif_sel_item() %in% f$virtual_map$item)
      plot_icc_frames(f, dif_sel_item(), group = dif_sel_vars())
    } else {
      req(dif_sel_item() %in% if (inherits(f, "rasch_mfrm"))
        f$virtual_map$item else f$items$item)
      plot_icc(f, dif_sel_item(), group = dif_sel_vars())
    }
  }, code = function() sprintf(
    if (inherits(fit(), "rasch_efrm"))
      'plot_icc_frames(fit, %s, group = %s)' else
        'plot_icc(fit, %s, group = %s)',
    qstr(dif_sel_item() %||% ""), qvec(dif_sel_vars())))

  # Primary post-hoc: main effects use pairwise marginal resolved logits;
  # interactions use tensor contrasts (difference-in-differences for two
  # factors), so the displayed magnitude belongs to the interaction itself.
  # Holm controls precisely the selected post-hoc family across multifactor
  # and repeated-measures structures.
  dif_posthoc_res <- reactive({
    f <- fit(); req(!is.null(f$factors))
    flg <- max(0.05, input$dif_size_flag %||% 0.5)
    mn <- max(2, input$dif_size_minn %||% 20)
    dr <- dif_res()
    tryCatch(dif_posthoc(
      f, dif_sel_item(), term = dif_sel_term(), factors = dr$factor_names,
      within = dr$within, flag_logits = flg, min_n = mn,
      p_adjust = "holm", alpha = dif_alpha()), error = function(e) e)
  })
  output$dif_posthoc_heading <- renderUI({
    intr <- length(dif_sel_vars()) > 1L
    title <- if (intr) "Interaction magnitudes" else "Pairwise marginal differences"
    help <- if (intr)
      paste("Difference-in-differences in resolved item locations; with more",
            "than two levels, every pair of level differences is compared.",
            "Holm adjustment controls the selected family.")
    else paste("Pairwise differences between resolved item locations, averaged",
               "equally over complete cells of the other person factors.",
               "Holm adjustment controls the selected family.")
    h6(span(title, info_icon(help)))
  })
  output$dif_posthoc_note <- renderUI({
    ph <- dif_posthoc_res()
    if (inherits(ph, "error"))
      return(p(class = "text-muted small mb-2", conditionMessage(ph)))
    s <- dif_res()$summary[dif_sel_row(), , drop = FALSE]
    flagged <- isTRUE(s$uniform_DIF) || isTRUE(s$nonuniform_DIF)
    base <- p(
      class = paste("small mb-2", if (flagged) "text-body" else "text-muted"),
      if (flagged)
        "The omnibus term is significant; use these adjusted comparisons to locate it."
      else
        "The omnibus term is not significant; treat individual comparisons as exploratory."
    )
    extra <- if (length(ph$notes))
      p(class = "text-muted small mb-2", paste(ph$notes, collapse = " "))
    else NULL
    tagList(base, extra)
  })
  output$dif_posthoc_tbl <- DT::renderDT({
    ph <- dif_posthoc_res(); req(!inherits(ph, "error"))
    d <- ph$table[, c("contrast", "estimate", "se", "statistic", "df",
                      "lower", "upper", "p_adj")]
    names(d)[names(d) == "estimate"] <- "magnitude"
    dt <- style_mag_red(num_dt(d), d, "magnitude", ph$flag_logits)
    style_lo_red(dt, d, "p_adj", ph$alpha)
  })

  # Secondary interaction-cell comparisons locate which cells generate the
  # interaction. They are deliberately not labelled as its magnitude.
  dif_size_res <- reactive({
    f <- fit(); req(!is.null(f$factors))
    vars <- dif_sel_vars(); req(length(vars) >= 1, all(vars %in% names(f$factors)))
    valid_items <- if (inherits(f, "rasch_mfrm") && !is.null(f$virtual_map))
      unique(f$virtual_map$item) else f$items$item
    req(dif_sel_item() %in% valid_items)
    flg <- max(0.05, input$dif_size_flag %||% 0.5)
    mn <- max(2, input$dif_size_minn %||% 20)
    tryCatch(dif_size(f, dif_sel_item(), by = vars,
                      p_adjust = "holm", flag_logits = flg, min_n = mn,
                      alpha = dif_alpha()),
             error = function(e) e)
  })
  # the resolved-location summary above the pairwise table; a muted line
  # appears only when dif_size genuinely errors (e.g. too few responders)
  output$dif_levels_note <- renderUI({
    ds <- dif_size_res()
    if (inherits(ds, "error"))
      return(p(class = "text-muted small mb-0", conditionMessage(ds)))
    p(class = "text-muted small mb-2",
      sprintf("%s resolved by %s: %s.", ds$item, ds$by,
              paste(sprintf("%s %.2f", ds$levels$level, ds$levels$location),
                    collapse = ", ")))
  })
  # always the pairwise magnitude table when dif_size succeeds: one row for a
  # two-level factor (the single effect size), more for > 2 levels
  output$dif_size_tbl <- DT::renderDT({
    ds <- dif_size_res()
    req(!inherits(ds, "error"))
    keep <- c("level_a", "level_b", "difference", "se",
              .app_wald_columns(ds$pairs), "lower", "upper", "p_adj",
              "ets", "signed_area")
    d <- ds$pairs[, intersect(keep, names(ds$pairs))]
    for (nm in intersect(c("ets", "signed_area"), names(d)))
      if (all(is.na(d[[nm]]))) d[[nm]] <- NULL
    dt <- style_mag_red(num_dt(d), d, "difference", ds$flag_logits)
    style_lo_red(dt, d, "p_adj", ds$alpha)
  })
  register_code("dif_posthoc_tbl", function() {
    flg <- max(0.05, input$dif_size_flag %||% 0.5)
    mn <- max(2, input$dif_size_minn %||% 20)
    extra <- paste0(
      if (flg != 0.5) sprintf(", flag_logits = %s", flg) else "",
      if (mn != 20) sprintf(", min_n = %s", mn) else "",
      if (dif_alpha() != 0.05) sprintf(", alpha = %s", dif_alpha()) else "")
    within <- dif_res()$within
    sprintf('dif_posthoc(fit, %s, term = %s, factors = %s%s, p_adjust = "holm"%s)',
            qstr(dif_sel_item() %||% ""), qstr(dif_sel_term() %||% ""),
            qvec(dif_res()$factor_names),
            if (length(within)) paste0(", within = ", qvec(within)) else "",
            extra)
  })
  register_code("dif_size_tbl", function() {
    flg <- max(0.05, input$dif_size_flag %||% 0.5)
    mn <- max(2, input$dif_size_minn %||% 20)
    extra <- paste0(
      if (flg != 0.5) sprintf(", flag_logits = %s", flg) else "",
      if (mn != 20) sprintf(", min_n = %s", mn) else "",
      if (dif_alpha() != 0.05) sprintf(", alpha = %s", dif_alpha()) else "")
    sprintf('dif_size(fit, %s, by = %s, p_adjust = "holm"%s)$pairs',
            qstr(dif_sel_item() %||% ""), qvec(dif_sel_vars()), extra)
  })
  output$dl_dif_posthoc <- downloadHandler(
    filename = function() "dif_posthoc.csv",
    content = function(file) {
      ph <- dif_posthoc_res()
      if (inherits(ph, "error")) stop(conditionMessage(ph))
      write_csv_plain(ph$table, file)
    })

  # planned DIF contrasts: the family is derived from the factor structure
  # and every question tested at once; repeated rows use the person identifier
  # already stored by rasch()
  contr_res <- reactiveVal(NULL)
  observeEvent(input$pc_run, {
    f <- fit()
    res <- tryCatch({
      if (is.null(f$factors) || !length(names(f$factors)))
        stop("nominate at least one person factor on the Data page")
      fac <- names(f$factors)
      if (!length(fac))
        stop("nominate at least one person factor")
      its <- if (length(input$pc_items)) input$pc_items else NULL
      withProgress(message = "Resolving items…",
                   detail = "one refit per item", value = 0.15,
        dif_contrasts(f, factors = fac, items = its, p_adjust = "holm"))
    }, error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Planned contrasts:", conditionMessage(res)),
                       type = "warning", duration = 8)
      contr_res(NULL)
    } else {
      res$run_factors <- fac
      res$run_items <- its
      # the person-identifier control was retired with the repeated-measures
      # redesign; the code footer adds an id argument only when one is set
      res$run_id <- NULL
      contr_res(res)
    }
  })
  # the derived family in words, shown once a run has completed
  output$contr_family <- renderUI({
    r <- contr_res()
    if (is.null(r)) return(NULL)
    tagList(
      h6("The planned family"),
      div(class = "small font-monospace mb-2",
          lapply(seq_len(nrow(r$family)), function(i)
            div(paste0(r$family$contrast[i],
                       if (r$family$within[i]) "  [within subjects]" else "")))),
      if (isTRUE(r$paired))
        p(class = "text-muted small mb-2",
          "Stacked design: tests use person-level residual scores."),
      if (length(r$notes))
        p(class = "text-muted small mb-2", paste(r$notes, collapse = "; ")))
  })
  register_table("contr_tbl", function() {
    r <- contr_res(); req(!is.null(r)); r$table
  }, function() {
    r <- contr_res()
    validate(need(!is.null(r),
                  "Press the button to derive the planned family from the factor structure and test it."))
    d <- curate(r$table, "contrast", full = isTRUE(input$contr_full))
    if ("within" %in% names(d)) d$within <- ifelse(d$within, "*", "")
    # no significant/practical flags: the estimate turns red at the practical
    # criterion, the adjusted p at alpha
    dt <- style_mag_red(num_dt(d), d, "estimate", r$flag_logits %||% 0.5)
    style_lo_red(dt, d, "p_adj", r$alpha %||% 0.05)
  }, code = function() {
    r <- contr_res(); req(!is.null(r))
    sprintf('dif_contrasts(fit, factors = %s%s%s, p_adjust = "holm")$table',
            qvec(r$run_factors),
            if (length(r$run_items))
              paste0(', items = ', qvec(r$run_items)) else "",
            if (!is.null(r$run_id))
              paste0(', id = ', qstr(r$run_id)) else "")
  }, csv_name = "dif_contrasts.csv")

  # ------------------------------------------------------------------- BTL --
  bfit <- reactive({
    validate(need(!is.null(btl_fit()),
                  "Run a Comparative Judgement analysis from the Data page to see results here."))
    s <- active_btl_step()
    if (is.null(s)) btl_fit() else s$fit
  })
  output$active_btlef <- reactive({
    b <- tryCatch(bfit(), error = function(e) NULL)
    inherits(b, "rasch_btl_efrm")
  })
  outputOptions(output, "active_btlef", suspendWhenHidden = FALSE)
  app_item_estimates <- function(f) {
    d <- if (inherits(f, "rasch_efrm")) f$item_arbitrary else
      if (inherits(f, "rasch_mfrm")) f$item_effects else f$items
    d[, c("item", "location", "se"), drop = FALSE]
  }
  change_estimates <- reactive({
    if (!is.null(btl_fit())) {
      before <- btl_fit()$objects[, c("object", "location", "se")]
      after <- bfit()$objects[, c("object", "location", "se")]
      names(before) <- c("parameter", "original", "original_se")
      names(after) <- c("parameter", "active", "active_se")
      kind <- "object"
    } else {
      before <- app_item_estimates(analysis())
      after <- app_item_estimates(fit())
      names(before) <- c("parameter", "original", "original_se")
      names(after) <- c("parameter", "active", "active_se")
      kind <- "item"
    }
    d <- merge(before, after, by = "parameter", all = TRUE, sort = FALSE)
    d$type <- kind
    d$status <- ifelse(is.na(d$original), "added",
                       ifelse(is.na(d$active), "removed", "retained"))
    d$change <- d$active - d$original
    d[, c("type", "parameter", "status", "original", "active", "change",
          "original_se", "active_se")]
  })
  register_table("change_est_tbl", function() change_estimates(),
                 function() num_dt(change_estimates()),
                 code = function() {
                   bt <- !is.null(tryCatch(btl_fit(), error = function(e) NULL))
                   paste(
                     if (bt)
                       'before <- original_bt$objects[, c("object", "location", "se")]'
                     else paste(
                       'item_estimates <- function(x) {',
                       '  d <- if (inherits(x, "rasch_efrm")) x$item_arbitrary else',
                       '    if (inherits(x, "rasch_mfrm")) x$item_effects else x$items',
                       '  d[, c("item", "location", "se"), drop = FALSE]',
                       '}',
                       'before <- item_estimates(original_fit)', sep = "\n"),
                     if (bt)
                       'after <- bt$objects[, c("object", "location", "se")]'
                     else 'after <- item_estimates(fit)',
                     'names(before) <- c("parameter", "original", "original_se")',
                     'names(after) <- c("parameter", "active", "active_se")',
                     'd <- merge(before, after, by = "parameter", all = TRUE, sort = FALSE)',
                     sprintf('d$type <- "%s"', if (bt) "object" else "item"),
                     'd$status <- ifelse(is.na(d$original), "added", ifelse(is.na(d$active), "removed", "retained"))',
                     'd$change <- d$active - d$original',
                     'd[, c("type", "parameter", "status", "original", "active", "change", "original_se", "active_se")]',
                     sep = "\n")
                 })

  change_persons <- reactive({
    b <- analysis(); a <- fit()
    nb <- nrow(b$person); na <- nrow(a$person)
    before <- data.frame(
      row = seq_len(nb),
      id = if (is.null(b$person$id)) seq_len(nb) else as.character(b$person$id),
      original = b$person$theta, original_se = b$person$se)
    after <- data.frame(row = seq_len(na), active = a$person$theta,
                        active_se = a$person$se)
    d <- merge(before, after, by = "row", all = TRUE, sort = FALSE)
    d$change <- d$active - d$original
    d[, c("row", "id", "original", "active", "change",
          "original_se", "active_se")]
  })
  register_table("change_person_tbl", function() change_persons(),
                 function() num_dt(change_persons(), page_len = 25),
                 code = function() paste(
                   'before <- data.frame(row = seq_len(nrow(original_fit$person)),',
                   '  id = if (is.null(original_fit$person$id)) seq_len(nrow(original_fit$person)) else as.character(original_fit$person$id),',
                   '  original = original_fit$person$theta, original_se = original_fit$person$se)',
                   'after <- data.frame(row = seq_len(nrow(fit$person)),',
                   '  active = fit$person$theta, active_se = fit$person$se)',
                   'd <- merge(before, after, by = "row", all = TRUE, sort = FALSE)',
                   'd$change <- d$active - d$original',
                   'd[, c("row", "id", "original", "active", "change", "original_se", "active_se")]',
                   sep = "\n"))
  btl_boot_table <- function(which) {
    bv <- btl_boot_val()
    if (is.null(bv) || is.null(bv$bs[[which]])) bfit()[[which]]
    else bv$bs[[which]]
  }
  output$btl_boxes <- renderUI({
    f <- bfit(); bv <- btl_boot_val()
    framed <- inherits(f, "rasch_btl_efrm")
    # BTL--EFRM requires repeated comparisons by judges. Its pooled Pearson
    # statistic is descriptive, and no row-based or stale ordinary-BTL
    # probability may be displayed for it.
    pair_p <- if (framed) NA_real_ else if (is.null(bv)) f$total_p
      else bv$bs$total$chisq_p_boot
    metric_grid(
      metric_tile("metric_objects", "Objects", nrow(f$objects), icon = "podium"),
      metric_tile("metric_comparisons", "Comparisons",
                  sprintf("%.0f", f$n_comparisons), icon = "pair"),
      if (!is.null(f$judges))
        metric_tile("metric_judges", "Judges", nrow(f$judges), icon = "balance"),
      metric_tile("metric_osi", "Object separation",
                  if (finite1(f$osi$PSI)) sprintf("%.3f", f$osi$PSI) else "—",
                  icon = "separation",
                  status = if (!finite1(f$osi$PSI)) "neutral"
                    else if (f$osi$PSI >= 0.7) "good" else "bad"),
      metric_tile("metric_pair_fit", "Pairwise fit p",
                  if (finite1(pair_p)) fmt_p(pair_p) else "Unavailable",
                  icon = "chisq",
                  status = if (!finite1(pair_p)) "neutral" else
                    if (pair_p >= 0.05) "good" else "bad"))
  })
  register_code("btl_boxes", function() {
    f <- bfit()
    bv <- btl_boot_val()
    paste(c(
      if (!inherits(f, "rasch_btl_efrm") && !is.null(bv))
        sprintf("bs <- fit_bootstrap(bt, B = %d, seed = %d)", bv$B, bv$seed),
      "list(",
      "  objects = nrow(bt$objects), comparisons = bt$n_comparisons,",
      "  judges = if (is.null(bt$judges)) NULL else nrow(bt$judges),",
      if (inherits(f, "rasch_btl_efrm"))
        "  object_separation = bt$osi$PSI, pairwise_fit_p = NA_real_  # unavailable with judge clustering" else if (is.null(bv))
        "  object_separation = bt$osi$PSI, pairwise_fit_p = bt$total_p" else
        "  object_separation = bt$osi$PSI, pairwise_fit_p = bs$total$chisq_p_boot",
      ")"), collapse = "\n")
  })
  # test-of-fit stat box (Summary page): the paired-comparison headline set
  # read off the fit; the CSV chip downloads the COMPLETE table from
  # fit_summary_table()'s rasch_btl method
  register_stat_box("btl_fitsum_tbl",
    csv_fun = function() {
      d <- fit_summary_table(bfit())
      bv <- btl_boot_val()
      if (is.null(bv)) return(d)
      rbind(d, data.frame(
        statistic = "Bootstrap pairwise chi-square probability",
        value = fmt_p(bv$bs$total$chisq_p_boot)))
    },
    csv_name = "fit_summary.csv",
    ui_fun = function() {
      f <- bfit(); bv <- btl_boot_val()
      framed <- inherits(f, "rasch_btl_efrm")
      polytomous <- !is.null(f$m) && f$m > 1L
      model_lab <- if (inherits(f, "rasch_btl_efrm"))
        "Comparative Judgement with frame-dependent units"
      else if (polytomous)
        sprintf("Polytomous paired comparisons (%d categories)", f$m + 1L)
      else "Comparative Judgement"
      conv <- if (isTRUE(f$converged) && inherits(f, "rasch_btl_efrm"))
        "two-stage estimation converged"
      else if (isTRUE(f$converged))
        sprintf("converged in %d iterations", f$iterations)
      else span(class = "text-danger",
                sprintf("did not converge in %d iterations", f$iterations))
      design <- paste(c(
        sprintf("%d objects", nrow(f$objects)),
        sprintf("%.0f comparisons", f$n_comparisons),
        if (!is.null(f$judges)) sprintf("%d judges", nrow(f$judges))),
        collapse = " · ")
      dep_rows <- if (!is.null(f$dependence)) {
        lapply(seq_len(nrow(f$dependence)), function(r) {
          has_adjusted <- !is.null(f$dependence$p_adj) &&
            length(f$dependence$p_adj) >= r &&
            is.finite(f$dependence$p_adj[r])
          inference <- if (has_adjusted)
            paste("Holm-adjusted", p_lab(f$dependence$p_adj[r])) else
              if (is.null(f$dependence$p_adj))
                "adjusted p unavailable; refit" else
                  "adjusted p unavailable; see fit notes"
          effect_label <- if (identical(f$dependence$effect[r], "position"))
            "First-position advantage" else
              sprintf("Within-judge %s",
                      gsub("_", "-", f$dependence$effect[r]))
          stat_row(effect_label,
                   sprintf("%.2f logits (%s)", f$dependence$estimate[r],
                           inference))
        })
      }
      tagList(
        div(class = "stat-head",
            model_lab, " · maximum likelihood · ", conv,
            if (isTRUE(f$clustered)) " · SEs clustered by judge"),
        stat_rows(
          stat_row("Pairwise chi-square",
                   {
                     pair_p <- if (framed) NA_real_ else if (is.null(bv))
                       f$total_p else bv$bs$total$chisq_p_boot
                     if (!finite1(pair_p))
                       sprintf("%.2f on %d df; probability unavailable%s",
                               f$total_chisq, f$total_df,
                               if (isTRUE(f$clustered)) " (judge clustering)"
                               else "")
                     else sprintf("%.2f on %d df, %s", f$total_chisq,
                                  f$total_df, p_lab(pair_p))
                   }),
          stat_row("Design", design),
          stat_row("Object separation index",
                   if (finite1(f$osi$PSI)) sprintf("%.3f", f$osi$PSI)
                   else "—"),
          if (polytomous && !is.null(f$thr_structure))
            stat_row("Threshold structure",
                     if (identical(f$thr_structure, "pc"))
                       "principal components (spread)"
                     else "free symmetric"),
          dep_rows))
    },
    code = function() {
      f <- bfit()
      bv <- btl_boot_val()
      if (inherits(f, "rasch_btl_efrm") || is.null(bv))
        "fit_summary_table(bt)" else
        paste(sprintf("bs <- fit_bootstrap(bt, B = %d, seed = %d)", bv$B, bv$seed),
              "fit_summary_table(bt)",
              "bs$total  # bootstrap pairwise fit probability", sep = "\n")
    })
  # routine handling notes (the old text panel printed bt$notes)
  output$btl_fitsum_notes <- renderUI({
    f <- bfit()
    if (!length(f$notes)) return(NULL)
    sprintf("Note. %s.", paste(f$notes, collapse = "; "))
  })
  register_table("btl_obj_tbl", function() btl_boot_table("objects"),
                 function() {
    bv <- btl_boot_val()
    d <- curate(btl_boot_table("objects"), "btl_obj",
                full = isTRUE(input$btl_full),
                extra = if (!is.null(bv))
                  c("fit_resid_p_boot_adj", "n_boot_fit_resid"))
    # fit residual, infit and outfit are flagged by num_dt, as on every table
    dt <- num_dt(d, selection = "single")
    if (!is.null(bv)) dt <- style_lo_red(dt, d, "fit_resid_p_boot_adj", .05)
    dt
  }, code = function() {
    bv <- btl_boot_val()
    if (is.null(bv)) return("# bt from the Data page\nbt$objects")
    paste(sprintf("bs <- fit_bootstrap(bt, B = %d, seed = %d)", bv$B, bv$seed),
          "bs$objects", sep = "\n")
  })
  # the object selected by clicking a row of the table drives the object
  # characteristic curve on the right (master-detail, as the item table does)
  sel_object <- reactive({
    b <- bfit(); i <- input$btl_obj_tbl_rows_selected
    b$objects$object[.app_selected_row(i, nrow(b$objects))]
  })
  register_table("btl_pairs_tbl", function() btl_boot_table("pairs"),
                 function() {
    d <- btl_boot_table("pairs")
    d$chisq <- NULL   # residual^2; redundant on screen, kept in the CSV
    num_dt(d, p_bold = if (!is.null(btl_boot_val())) "chisq_p_boot_adj")
  }, code = function() {
    bv <- btl_boot_val()
    if (is.null(bv)) "bt$pairs" else
      paste(sprintf("bs <- fit_bootstrap(bt, B = %d, seed = %d)", bv$B, bv$seed),
            "bs$pairs", sep = "\n")
  })
  register_table("btl_judges_tbl", function() {
    validate(need(!is.null(bfit()$judges), "No judge column was nominated."))
    btl_boot_table("judges")
  }, function() {
    validate(need(!is.null(bfit()$judges), "No judge column was nominated."))
    bv <- btl_boot_val(); d <- btl_boot_table("judges")
    d$misfit <- if (!is.null(bv))
      ifelse(!is.na(d$fit_resid_p_boot_adj) &
               d$fit_resid_p_boot_adj < .05, "*", "") else
      ifelse(!is.na(d$fit_resid) & abs(d$fit_resid) > 2.5, "*", "")
    d <- curate(d, "btl_judge", full = isTRUE(input$btl_judges_full),
                extra = if (!is.null(bv))
                  c("fit_resid_p_boot_adj", "n_boot_fit_resid"))
    dt <- num_dt(d, selection = "single")
    if (!is.null(bv)) dt <- style_lo_red(dt, d, "fit_resid_p_boot_adj", .05)
    dt
  }, code = function() {
    bv <- btl_boot_val()
    if (is.null(bv)) "bt$judges" else
      paste(sprintf("bs <- fit_bootstrap(bt, B = %d, seed = %d)", bv$B, bv$seed),
            "bs$judges", sep = "\n")
  })
  # the judge whose surprise map is drawn: the selected row, defaulting to 1
  # (master-detail, exactly as the person table drives the kidmap)
  sel_judge <- reactive({
    j <- bfit()$judges; req(!is.null(j))
    n <- input$btl_judges_tbl_rows_selected
    as.character(j$judge[.app_selected_row(n, nrow(j))])
  })
  register_plot("btl_judge_map",
                function() plot_btl_judge_map(bfit(), sel_judge()),
                w = 7, h = 5.5, code = function()
                  sprintf('plot_btl_judge_map(bt, %s)', qstr(sel_judge())))
  # hover identification for the unexpected-judgements map: only the
  # surprising matchups are text-labelled on the plot (see R/btl-independence.R),
  # so most of the two points per matchup (the stronger and weaker object's
  # ends of each segment) are otherwise anonymous. judge_pair_surprise() is
  # the same call plot_btl_judge_map() makes internally, with its own
  # defaults (min_n = 1, flag_z = 1.96), which the app never overrides.
  btl_judge_pairs_res <- reactive(judge_pair_surprise(bfit(), sel_judge()))
  register_hover_tip("btl_judge_map", function() {
    p <- btl_judge_pairs_res()$pairs
    req(nrow(p) > 0)
    rbind(
      data.frame(z = p$z, y = p$loc_hi, p_adj = p$p_adj, tied = p$tied,
                object = p$object_hi,
                opponent = p$object_lo, stringsAsFactors = FALSE),
      data.frame(z = p$z, y = p$loc_lo, p_adj = p$p_adj, tied = p$tied,
                object = p$object_lo,
                opponent = p$object_hi, stringsAsFactors = FALSE))
  }, "z", "y", function(np)
    sprintf("%s vs %s · residual %+.2f · Holm p %s · location %.2f%s",
            np$object[1], np$opponent[1], np$z[1],
            if (is.finite(np$p_adj[1]))
              format.pval(np$p_adj[1], digits = 3, eps = 0.001) else "NA",
            np$y[1], if (isTRUE(np$tied[1])) " · tied locations" else ""))
  register_plot("btl_plot", function() plot_btl(bfit()),
                code = function() "# bt from the Data page\nplot_btl(bt)")
  # object characteristic curve: model expected response against opponent
  # location with per-opponent observed means (dichotomous and polytomous fits)
  observeEvent(btl_fit(), {
    b <- btl_fit()
    if (!is.null(b)) {
      jf <- input$bt_jfactors
      updateSelectizeInput(session, "bdif_factors",
                           choices = if (length(jf)) jf else character(0),
                           selected = jf)
    }
    # Results computed on request belong to the fit they came from. A project
    # restore is the exception: those results were saved with this exact fit
    # and are reinstated later in the same flush cycle.
    if (!isTRUE(restoring_project())) {
      clear_btl_fit_results()
      restored_dimensionality(NULL)
    }
  }, ignoreNULL = FALSE)
  register_plot("btl_occ", function() {
    o <- sel_object(); req(o %in% bfit()$objects$object)
    plot_btl_icc(bfit(), o)
  }, w = 8, h = 5.5, code = function()
    paste0(sprintf('plot_btl_icc(bt, %s)', qstr(sel_object() %||% "")),
           "\n# all objects: one page each in the PDF download"))
  output$btl_occ_all_pdf <- downloadHandler(
    filename = function() "occ_all_objects.pdf",
    content = function(file) {
      b <- bfit()
      pdf(file, width = 8, height = 5.5, onefile = TRUE)
      on.exit(dev.off(), add = TRUE)
      # a failed object still gets its page: a placeholder naming the
      # error, so no object silently vanishes from the batch
      for (o in b$objects$object)
        tryCatch(plot_btl_icc(b, o), error = function(e) {
          plot.new()
          text(0.5, 0.5, sprintf("%s: %s", o, conditionMessage(e)))
        })
    })
  # polytomous (ordinal) fits carry thresholds and category curves; the flag
  # hides both cards entirely for dichotomous fits
  output$btl_graded <- reactive({
    b <- btl_fit()
    !is.null(b) && !is.null(b$m) && b$m >= 2
  })
  outputOptions(output, "btl_graded", suspendWhenHidden = FALSE)
  register_table("btl_thr_tbl", function() bfit()$thresholds, function() {
    validate(need(!is.null(bfit()$thresholds),
                  "Dichotomous fit: no threshold structure to show."))
    num_dt(bfit()$thresholds)
  }, code = function() "bt$thresholds")
  register_table("btl_comp_tbl", function() bfit()$components, function() {
    validate(need(!is.null(bfit()$components),
                  "Dichotomous fit: no threshold components to show."))
    num_dt(bfit()$components)
  }, code = function() "bt$components")
  register_plot("btl_cats", function() plot_btl_categories(bfit()),
                code = function() "plot_btl_categories(bt)")

  # within-judge dependence (Independence > Local): exposure and carry-over
  # effects estimated when a judgment-order column was nominated
  # keyed on the comparison-level data, not the effects table: when every
  # effect is dropped (no informative comparisons, or separation) the panel
  # must still show the diagnostic table that explains why
  output$has_btl_dep <- reactive({
    b <- tryCatch(bfit(), error = function(e) NULL)
    # This page concerns the exposure and carry-over history covariates. A
    # position-only fit still retains row-level position data, but has no
    # judgment sequence; report it in the fit summary instead of opening the
    # history panel.
    !is.null(b) && .app_has_btl_history(b)
  })
  outputOptions(output, "has_btl_dep", suspendWhenHidden = FALSE)
  register_table("btl_dep_tbl", function() {
    validate(need(!is.null(bfit()$dependence_data),
                  "Nominate a judgment-order column in the Data roles to estimate within-judge dependence."))
    d <- bfit()$dependence
    if (!is.null(d) && !"p_adj" %in% names(d)) {
      d$p_adj <- NA_real_
      if ("significant" %in% names(d)) d$significant <- NA
    }
    d
  }, function() {
    validate(need(!is.null(bfit()$dependence_data),
                  "Nominate a judgment-order column in the Data roles to estimate within-judge dependence."))
    validate(need(!is.null(bfit()$dependence),
                  paste("No dependence effect was estimable:",
                        paste(bfit()$notes, collapse = "; "))))
    d <- bfit()$dependence
    if (!"p_adj" %in% names(d)) {
      d$p_adj <- NA_real_
      if ("significant" %in% names(d)) d$significant <- NA
    }
    num_dt(d, p_bold = "p_adj")
  }, code = function() "bt$dependence")
  # the graphical display of the selected dependence effect
  register_plot("btl_dep_plot", function() {
    b <- bfit()
    e <- input$btl_dep_effect %||% "exposure"
    # the first-position advantage is a single constant added to every
    # comparison, not a binned history covariate, so the departure-vs-covariate
    # display does not apply; it (and a fit with no order history at all) is
    # read from the table instead
    validate(need(!is.null(b$dependence_data) && e != "position",
                  "The first-position effect is a constant, not a history covariate; read it from the table."))
    validate(need(!is.null(b$dependence),
                  paste("No dependence effect was estimable:",
                        paste(b$notes, collapse = "; "))))
    validate(need(e %in% b$dependence$effect,
                  "This effect was not estimated for this fit (no order column, no informative comparisons, or separated; see the notes). Estimated effects are read from the table."))
    plot_btl_dependence(b, e)
  }, w = 8, h = 5.5, code = function()
    sprintf('plot_btl_dependence(bt, %s)',
            qstr(input$btl_dep_effect %||% "exposure")))
  # every comparison with its history covariates, for deep interrogation
  register_table("btl_dep_comps", function() bfit()$dependence_data,
                 function() {
    validate(need(!is.null(bfit()$dependence_data),
                  "Nominate a judgment-order column in the Data roles."))
    num_dt(bfit()$dependence_data, page_len = 20)
  }, code = function() "bt$dependence_data")

  # ---------------------- BTL trait dimensionality (loops + swirl) --------
  # the pair-structure analogues of transitivity (one order?) and residual
  # PCA (a second attribute steering contests?); both cached per fit
  btl_trans <- reactive({ b <- bfit(); req(!is.null(b)); btl_transitivity(b) })
  btl_dim <- reactive({
    b <- bfit(); req(!is.null(b))
    saved <- restored_dimensionality()
    if (!is.null(saved) &&
        !inherits(tryCatch(.validate_btl_dimensionality(saved, b),
                           error = function(e) e), "error")) return(saved)
    btl_dimensionality(b, seed = 1L)
  })
  register_plot("btl_scree", function() plot_btl_scree(btl_dim()),
                w = 7, h = 5, code = function()
                  "plot_btl_scree(btl_dimensionality(bt, seed = 1))")
  register_plot("btl_dim_map", function() plot_btl_dim_map(btl_dim()),
                w = 7, h = 5.5, code = function()
                  "plot_btl_dim_map(btl_dimensionality(bt, seed = 1))")
  register_table("btl_bimensions_tbl", function() btl_dim()$bimensions,
                 function() num_dt(btl_dim()$bimensions),
                 code = function() "btl_dimensionality(bt, seed = 1)$bimensions")
  register_table("btl_trans_tbl", function() btl_trans()$summary,
                 function() num_dt(btl_trans()$summary),
                 code = function() "btl_transitivity(bt)$summary")
  # structural lens (Trait tab): which objects sit in the most loops
  register_plot("btl_involve_plot",
                function() plot_btl_transitivity(btl_trans(), by = "object"),
                w = 7, h = 5, code = function()
                  'plot_btl_transitivity(btl_transitivity(bt), by = "object")')
  # per-judge lens (Persons tab): each judge's consistency
  register_plot("btl_judge_consist", function() {
    tr <- btl_trans()
    validate(need(!is.null(tr$judges),
      "No judge reached enough compared triples for a consistency estimate."))
    plot_btl_transitivity(tr, by = "judge")
  }, w = 7, h = 5, code = function()
    'plot_btl_transitivity(btl_transitivity(bt), by = "judge")')
  register_table("btl_trans_judges_tbl",
                 function() btl_trans()$judges,
                 function() {
                   j <- btl_trans()$judges
                   validate(need(!is.null(j),
                     "Per-judge consistency needs judges with enough compared triples."))
                   num_dt(j)
                 }, code = function() "btl_transitivity(bt)$judges")

  # ------------------------------------------------ BTL targeting (design --
  # information, the targeting plot, and the adaptive next-pairs recommender)
  btl_info <- reactive({ b <- bfit(); btl_information(b) })
  register_table("btl_info_tbl", function() btl_info()$objects,
                 function() num_dt(btl_info()$objects),
                 code = function() "btl_information(bt)$objects")
  register_plot("btl_targeting_plot", function() plot_btl_targeting(bfit()),
                w = 9, h = 6, code = function() "plot_btl_targeting(bt)")
  # the adaptive next-pairs recommender (numeric count + the SE-weighting
  # switch drive both the table and its reproducible code)
  btl_next_args <- reactive({
    nn <- input$btl_next_n %||% 10
    list(n = if (is.null(nn) || is.na(nn) || nn < 3) 10L else as.integer(nn),
         weight_se = isTRUE(input$btl_next_wse %||% TRUE))
  })
  register_table("btl_next_tbl",
                 function() { a <- btl_next_args()
                   btl_next_pairs(bfit(), n = a$n, weight_se = a$weight_se) },
                 function() { a <- btl_next_args()
                   num_dt(btl_next_pairs(bfit(), n = a$n,
                                         weight_se = a$weight_se)) },
                 code = function() { a <- btl_next_args()
                   sprintf("btl_next_pairs(bt, n = %d, weight_se = %s)",
                           a$n, if (a$weight_se) "TRUE" else "FALSE") })

  # ------------------------------------------------- BTL common-object equating --
  # a banked reference calibration (object, location, optional se), parsed defensively;
  # btl_equate() is carried as a value or an error so <2-common-object and
  # other failures surface as a message rather than a red crash
  bt_eq_bank <- reactive({
    a <- if (!is.null(input$bt_eq_file))
      tryCatch(read.csv(input$bt_eq_file$datapath, check.names = FALSE,
                        stringsAsFactors = FALSE), error = function(e) NULL)
    else restored_project_resources()[["bt_eq_bank"]]
    if (is.null(a)) return(NULL)
    if (!all(c("object", "location") %in% names(a))) {
      showNotification("Reference CSV needs columns object and location - ignored.",
                       type = "warning")
      return(NULL)
    }
    if (anyDuplicated(names(a))) {
      showNotification("Reference CSV has duplicate column names - ignored.",
                       type = "warning")
      return(NULL)
    }
    missing_object <- is.na(a$object)
    a$object <- trimws(as.character(a$object))
    a$object[missing_object] <- NA_character_
    if (!"se" %in% names(a)) a$se <- NA_real_
    if (is.logical(a$se) && all(is.na(a$se))) a$se <- NA_real_
    if ("m" %in% names(a)) {
      m <- suppressWarnings(as.numeric(as.character(a$m)))
      m <- unique(m[is.finite(m)])
      if (length(m) == 1L && m >= 1L && m == floor(m) &&
          m <= .Machine$integer.max)
        attr(a, "m") <- as.integer(m)
    }
    a
  })
  bt_eq_reference_code <- function() {
    source <- if (!is.null(input$bt_eq_file)) sprintf(
      "bank <- read.csv(%s, check.names = FALSE, stringsAsFactors = FALSE)",
      qstr(input$bt_eq_file$name %||% "reference.csv")) else
        'bank <- project$resources[["bt_eq_bank"]]'
    paste(source,
          paste0("bank$object <- ifelse(is.na(bank$object), NA_character_, ",
                 "trimws(as.character(bank$object)))"),
          'if (!"se" %in% names(bank)) bank$se <- NA_real_',
          'if ("m" %in% names(bank)) {',
          '  bank_m <- suppressWarnings(as.numeric(as.character(bank$m)))',
          '  bank_m <- unique(bank_m[is.finite(bank_m)])',
          paste0('  if (length(bank_m) == 1L && bank_m >= 1L && ',
                 'bank_m <= .Machine$integer.max && bank_m == floor(bank_m)) ',
                 'attr(bank, "m") <- as.integer(bank_m)'),
          '}', sep = "\n")
  }
  bt_equate <- reactive({
    bank <- bt_eq_bank()
    validate(need(!is.null(bank),
                  "Upload a reference calibration to test drift against it."))
    tryCatch(btl_equate(bfit(), bank,
                        independent = isTRUE(input$bt_eq_independent)),
             error = function(e) e)
  })
  output$btl_eq_summary <- renderUI({
    req(!is.null(bt_eq_bank()))
    r <- bt_equate(); req(!inherits(r, "error"))
    label <- if (identical(r$shift_method, "unweighted"))
      sprintf("Unweighted descriptive shift %.3f logits", r$shift)
    else sprintf("Precision-weighted shift %.3f ± %s logits", r$shift,
                 if (is.finite(r$shift_se)) sprintf("%.3f", r$shift_se)
                 else "withheld")
    p(class = "text-muted small mb-0 mt-2",
      sprintf("%s over %d common object%s.", label, r$n_common,
              if (r$n_common == 1L) "" else "s"))
  })
  # surface a failed btl_equate() as its own message. need()'s `message` is
  # forced eagerly (shiny::need calls force(message)), so conditionMessage(r)
  # must sit behind the inherits() guard -- inside need() it would also be
  # evaluated on a SUCCESSFUL equate and error on the non-condition object.
  bt_eq_ok <- function() {
    r <- bt_equate()
    if (inherits(r, "error")) validate(need(FALSE, conditionMessage(r)))
    r
  }
  register_table("btl_eq_tbl", function() {
    bt_eq_ok()$table
  }, function() {
    d <- bt_eq_ok()$table
    style_lo_red(num_dt(d), d, "p_adj", 0.05)
  }, code = function() paste0(
    bt_eq_reference_code(), "\n",
    sprintf('btl_equate(bt, bank, independent = %s)$table',
            if (isTRUE(input$bt_eq_independent)) "TRUE" else "FALSE")))
  register_plot("btl_eq_plot", function() {
    bt_eq_ok()
    plot_btl_equate(bfit(), bt_eq_bank(),
                    independent = isTRUE(input$bt_eq_independent))
  }, w = 8, h = 6, code = function() paste0(
    bt_eq_reference_code(), "\n",
    sprintf('plot_btl_equate(bt, bank, independent = %s)',
            if (isTRUE(input$bt_eq_independent)) "TRUE" else "FALSE")))
  # hover identification for the equating plot: only drifting objects are
  # text-labelled (see plot_btl_equate(), R/btl-equating.R); bt_equate()$table
  # is the exact table plot_btl_equate() draws from (same fit1/fit2, no extra
  # args either call), so reusing it guarantees identical rows.
  register_hover_tip("btl_eq_plot", function() {
    r <- bt_equate(); req(!inherits(r, "error"))
    r$table
  }, "location_2", "location_1", function(np)
    sprintf("%s · calibration 2 %.3f · calibration 1 %.3f",
            np$object[1], np$location_2[1], np$location_1[1]))

  # -------------------------------------------- BTL DIF by judge group --
  # The nominated judge factors as named judge -> level maps. Build the maps
  # only after checking the defining requirement: a judge has at most one
  # observed value of each factor. The observed judges come from the fitted
  # comparisons, so rows omitted from the analysis cannot add spurious names.
  bdif_fit_source <- function(fit = NULL) {
    f <- fit %||% tryCatch(bfit(), error = function(e) NULL)
    source <- if (!is.null(f)) .app_fit_source(f) else NULL
    # Structural CJ refits may not carry attributes through their returned
    # object. Their comparison design is still downstream of the original
    # calibration, whose authenticated source metadata remains authoritative.
    if (is.null(source)) {
      base <- tryCatch(btl_fit(), error = function(e) NULL)
      source <- if (!is.null(base)) .app_fit_source(base) else NULL
    }
    source
  }
  bdif_judge_role <- function(fit, df, jc) {
    .app_btl_dif_judge_role(bdif_fit_source(fit), df, jc)
  }
  bdif_factor_maps <- function() {
    if (is.null(btl_fit()))
      stop("run a Comparative Judgement analysis first")
    fcs <- input$bdif_factors
    if (is.null(fcs) || !length(fcs))
      stop("choose one or more judge factors in the sidebar")
    df <- raw_data()
    jc <- input$bt_judge
    if (is.null(jc) || identical(jc, NONE) || !jc %in% names(df))
      stop("judge-group DIF needs the judge column nominated on the Data page")
    role <- bdif_judge_role(btl_fit(), df, jc)
    jd <- role$ids
    judges <- unique(as.character(btl_fit()$comparisons$judge))
    judges <- judges[!is.na(judges)]
    # The saved calibration data may predate a correction to judge metadata.
    # Reuse the authenticated downstream maps after reopening, but only for
    # that unchanged data, judge role and calibration. A fresh data source
    # clears restored resources; changed metadata must take precedence.
    restored <- restored_project_resources()$btl_dif_factors
    saved_maps <- if (!is.null(restored) &&
        identical(df, restored$data) && identical(jc, restored$judge_col) &&
        .fit_boot_signature_matches(restored$fit_signature, btl_fit()))
      restored$factors else list()
    setNames(lapply(fcs, function(fc) {
      if (fc %in% names(saved_maps)) return(saved_maps[[fc]])
      if (!fc %in% names(df)) stop("column not found: ", fc)
      value <- as.character(df[[fc]])
      by_judge <- lapply(judges, function(j) {
        z <- value[!is.na(jd) & jd == j]
        unique(z[!is.na(z)])
      })
      varying <- lengths(by_judge) > 1L
      if (any(varying))
        stop("factor '", fc, "' varies within judge(s) ",
             paste(judges[varying], collapse = ", "),
             ": judge-group DIF needs judge-constant factors")
      setNames(vapply(by_judge, function(z)
        if (length(z)) z[1L] else NA_character_, ""), judges)
    }), fcs)
  }
  # the factors of one displayed ANOVA term, matched against a run's factor
  # maps: a factor name that itself contains ":" must match the whole term
  # before the term is split into pieces
  bdif_term_vars <- function(term, maps, row = NULL) {
    r <- bdif_res()
    if (!is.null(row) && !is.null(r$summary_factors) &&
        length(r$summary_factors) >= row)
      return(r$summary_factors[[row]])
    if (term %in% names(maps)) return(term)
    intersect(strsplit(term, ":", fixed = TRUE)[[1]], names(maps))
  }
  # the judge -> cell map for one ANOVA term of the DISPLAYED run (frozen at
  # run time, so later sidebar edits cannot silently regroup the overlay)
  bdif_term_group <- function(term, row = NULL) {
    r <- bdif_res()
    maps <- if (!is.null(r)) r$bootstrap_design$factors else NULL
    if (is.null(r) || !is.list(maps) || !length(maps))
      stop("run the DIF analysis first")
    vars <- bdif_term_vars(term, maps, row)
    if (!length(vars))
      stop("the selected term no longer matches the run's factors; re-run")
    js <- names(maps[[1]])
    setNames(as.character(app_factor_cells(as.data.frame(
      lapply(vars, function(v) maps[[v]][js]), check.names = FALSE),
      sep = ":")), js)
  }
  # the per-opponent, per-group means the *grouped* branch of plot_btl_icc()
  # draws for the DIF overlay (R/btl.R) -- unlike the ungrouped OCC, that
  # branch never text-labels its points, so hover is the only way to read an
  # opponent off the plot. Mirrors that branch's aggregation (including its
  # min_n = 10 default, which the app's plot_btl_icc() call never overrides)
  # exactly, so the reproduced rows match what is actually drawn.
  bdif_occ_points <- function(fit, object, group, min_n = 10) {
    ob <- fit$objects
    m <- if (is.null(fit$m)) 1L else fit$m
    cm <- fit$comparisons
    gv <- if (!is.null(names(group)))
      unname(as.character(group)[match(cm$judge, names(group))])
    else as.character(group)
    sel_a <- cm$object_a == object; sel_b <- cm$object_b == object
    opp <- c(cm$object_b[sel_a], cm$object_a[sel_b])
    resp <- c(cm$response[sel_a], m - cm$response[sel_b])
    wt <- c(cm$weight[sel_a], cm$weight[sel_b])
    gg <- c(gv[sel_a], gv[sel_b])
    keep <- opp %in% ob$object & !is.na(gg)
    opp <- opp[keep]; resp <- resp[keep]; wt <- wt[keep]; gg <- gg[keep]
    levs <- sort(unique(gg))
    do.call(rbind, lapply(levs, function(lv) {
      sel <- gg == lv
      nn <- tapply(wt[sel], opp[sel], sum)
      om <- tapply(wt[sel] * resp[sel], opp[sel], sum) / nn
      om <- om[nn >= min_n]
      data.frame(opponent = names(om),
                location = ob$location[match(names(om), ob$object)],
                mean = as.numeric(om), group = lv, stringsAsFactors = FALSE)
    }))
  }
  # quote source columns as R symbols in emitted code
  bq <- .app_quote_name
  # A completed run owns its judge-factor maps. The comparison source is
  # frozen at calibration, but judge metadata may have been edited since
  # then; rebuilding maps from `dat` need not reproduce the displayed DIF.
  # Before a run, show the same constancy check used by the app.
  bdif_code_grp <- function() {
    r <- bdif_res()
    if (!is.null(r)) {
      maps <- r$bootstrap_design$factors
      if (!is.list(maps) || !length(maps))
        stop("the saved DIF analysis has no judge-factor maps; rerun it before reproducing its R code")
      return(paste0("# Judge factors used in this DIF analysis\n",
                    "factors <- ",
                    paste(capture.output(dput(maps)), collapse = "\n")))
    }
    fcs <- input$bdif_factors %||% "factor"
    jc <- input$bt_judge %||% "judge"
    source <- bdif_fit_source()
    if (is.null(source))
      stop("the fitted Comparative Judgement calibration has no saved source metadata; refit it before generating judge-factor code")
    if (!is.null(source) && is.list(source$settings) &&
        is.character(source$settings$bt_judge) &&
        length(source$settings$bt_judge) == 1L &&
        !is.na(source$settings$bt_judge) &&
        nzchar(source$settings$bt_judge))
      # Before a run, emit code for the fitted role even if the live sidebar
      # has been edited. A completed run emits its own frozen maps above.
      jc <- source$settings$bt_judge
    one <- function(fc) sprintf("%s = judge_factor(dat$%s, %s)",
                                bq(fc), bq(fc), qstr(fc))
    paste0(
      "judge_ids <- unique(as.character(bt$comparisons$judge))\n",
      "judge_ids <- judge_ids[!is.na(judge_ids)]\n",
      "judge_factor <- function(x, label) {\n",
      sprintf("  j <- as.character(dat$%s)\n", bq(jc)),
      "  u <- lapply(judge_ids, function(id) {\n",
      "    z <- as.character(x)[!is.na(j) & j == id]\n",
      "    unique(z[!is.na(z)])\n",
      "  })\n",
      "  if (any(lengths(u) > 1L))\n",
      "    stop(sprintf(\"factor '%s' varies within judge\", label))\n",
      "  setNames(vapply(u, function(z) if (length(z)) z[1L] else NA_character_,\n",
      "                  character(1)), judge_ids)\n",
      "}\n",
      "factors <- list(\n  ",
      paste(vapply(fcs, one, ""), collapse = ",\n  "), ")")
  }
  bdif_alpha <- reactive(clamp01(input$bdif_alpha, 0.05))
  bdif_effects <- reactive(input$bdif_effects %||% "main")
  # the displayed run's own settings, for styling and snippets (falling back
  # to the live sidebar before any run)
  bdif_shown_alpha <- function() {
    r <- bdif_res(); if (!is.null(r)) r$alpha else bdif_alpha()
  }
  bdif_shown_effects <- function() {
    r <- bdif_res(); if (!is.null(r)) r$effects else bdif_effects()
  }
  output$bdif_multifactor <- reactive(length(input$bdif_factors) > 1L)
  outputOptions(output, "bdif_multifactor", suspendWhenHidden = FALSE)
  # judge factors nominated (or changed) after the run still reach the factor
  # select; the results themselves are computed on request only
  observeEvent(input$bt_jfactors, {
    jf <- input$bt_jfactors
    keep <- intersect(input$bdif_factors, jf)
    updateSelectizeInput(session, "bdif_factors",
                         choices = if (length(jf)) jf else character(0),
                         selected = if (length(keep)) keep else jf)
  }, ignoreNULL = FALSE)
  bdif_res <- reactiveVal(NULL)
  # Display-only run state is kept outside the signed btl_dif() result. The
  # judge-factor maps themselves already live in bootstrap_design; only the
  # source judge-column name is additional app metadata.
  bdif_meta <- reactiveVal(NULL)
  observeEvent(input$bdif_run, {
    dif_boot_val(NULL)
    active_bt <- tryCatch(bfit(), error = function(e) e)
    if (inherits(active_bt, "rasch_btl_efrm")) {
      showNotification(paste(
        "Judge-group DIF resolution is an equal-unit Comparative Judgement",
        "procedure and cannot be applied to the frame-adjusted fit. Undo the",
        "frame adjustment before running DIF."), type = "error", duration = 10)
      bdif_res(NULL)
      bdif_meta(NULL)
      return()
    }
    maps <- tryCatch(bdif_factor_maps(), error = function(e) e)
    if (inherits(maps, "error")) {
      showNotification(conditionMessage(maps), type = "error", duration = 10)
      bdif_res(NULL); bdif_meta(NULL); return()
    }
    r <- withProgress(message = "Resolving objects by judge group…",
                      value = 0.4,
                      tryCatch(btl_dif(active_bt,
                                       factors = maps,
                                       effects = bdif_effects(),
                                       p_adjust = "holm",
                                       alpha = bdif_alpha()),
                               error = function(e) e))
    if (inherits(r, "error")) {
      showNotification(paste("DIF analysis failed:", conditionMessage(r)),
                       type = "error", duration = 10)
      bdif_res(NULL)
      bdif_meta(NULL)
    } else {
      bdif_res(r)
      role <- bdif_judge_role(active_bt, raw_data(), input$bt_judge)
      bdif_meta(list(judge_col = input$bt_judge,
                     fitted_judge_col = role$fitted_col))
    }
  })
  observeEvent(input$bdif_boot_run, {
    f <- tryCatch(bfit(), error = function(e) e)
    d <- bdif_res()
    if (inherits(f, "error") || is.null(d)) {
      showNotification("Run the Comparative Judgement DIF analysis first.",
                       type = "warning", duration = 7)
      return()
    }
    start_bootstrap(f, input$bdif_boot_B, input$bdif_boot_seed,
                    "bdif", d)
  })
  output$has_bdif_boot <- reactive({
    z <- dif_boot_val(); !is.null(z) && identical(z$kind, "bdif")
  })
  outputOptions(output, "has_bdif_boot", suspendWhenHidden = FALSE)
  output$bdif_boot_state <- renderUI({
    bv <- dif_boot_val()
    if (is.null(bv) || !identical(bv$kind, "bdif")) return(NULL)
    span(class = "badge bg-primary-subtle text-primary-emphasis pb-2 mb-2",
         sprintf(paste0("Sensitivity analysis active: %d of %d used; %d did ",
                        "not converge; %d otherwise failed; seed %d"),
                 bv$db$B_used, bv$B,
                 bv$db$B_nonconverged %||% NA_integer_,
                 bv$db$B_errors %||% NA_integer_, bv$seed))
  })
  register_table("bdif_anova_tbl", function() {
    r <- bdif_res(); req(!is.null(r)); r$summary
  }, function() {
    r <- bdif_res()
    validate(need(!is.null(r),
                  "Choose one or more judge factors in the sidebar and run the DIF analysis."))
    d <- r$summary[, intersect(c("object", "term", "F_uniform",
                                 "p_uniform_adj", "min_judges",
                                 "min_effective_judges", "uniform_inference",
                                 "F_nonuniform", "p_nonuniform_adj",
                                 "nonuniform_inference", "superseded"),
                               names(r$summary)), drop = FALSE]
    # a superseded row's flags are read on its higher-order term instead
    d$superseded <- ifelse(d$superseded, "(superseded)", "")
    # no boolean DIF flags: the adjusted probabilities turn red at the RUN's
    # alpha (frozen with the result, not the live sidebar)
    dt <- style_lo_red(num_dt(d, selection = "single"), d,
                       "p_uniform_adj", bdif_shown_alpha())
    style_lo_red(dt, d, "p_nonuniform_adj", bdif_shown_alpha())
  }, code = function()
    paste0("# dat: the comparison data; bt: the fit from the Data page\n",
           bdif_code_grp(), "\n",
           sprintf(paste0('btl_dif(bt, factors, effects = %s, ',
                          'p_adjust = %s, alpha = %s)$summary'),
                   qstr(bdif_shown_effects()),
                   qstr(bdif_res()$p_adjust), bdif_shown_alpha())))
  register_table("bdif_boot_tbl", function() {
    bv <- dif_boot_val(); req(!is.null(bv), identical(bv$kind, "bdif"))
    bv$db$summary
  }, function() {
    bv <- dif_boot_val()
    validate(need(!is.null(bv) && identical(bv$kind, "bdif"),
                  "Run the bootstrap sensitivity analysis first."))
    d <- bv$db$summary
    if (!isTRUE(input$bdif_boot_full))
      d <- d[, intersect(c(
        "object", "term", "p_uniform_adj", "p_uniform_boot_adj",
        "p_nonuniform_adj", "p_nonuniform_boot_adj"), names(d)), drop = FALSE]
    dt <- num_dt(d)
    dt <- style_lo_red(dt, d, "p_uniform_boot_adj", bdif_shown_alpha())
    style_lo_red(dt, d, "p_nonuniform_boot_adj", bdif_shown_alpha())
  }, code = function() {
    bv <- dif_boot_val(); req(!is.null(bv), identical(bv$kind, "bdif"))
    paste0(bdif_code_grp(), "\n",
      sprintf(paste0(
        "d <- btl_dif(bt, factors, effects = %s, p_adjust = %s, alpha = %s)\n",
        "dif_bootstrap(bt, d, B = %d, workers = 4, seed = %d)$summary"),
        qstr(bdif_shown_effects()), qstr(bdif_res()$p_adjust),
        bdif_shown_alpha(), bv$B, bv$seed))
  })
  register_table("bdif_sizes_tbl", function() {
    r <- bdif_res(); req(!is.null(r), !is.null(r$sizes)); r$sizes
  }, function() {
    r <- bdif_res()
    validate(need(!is.null(r),
                  "Choose one or more judge factors in the sidebar and run the DIF analysis."))
    validate(need(!is.null(r$sizes),
                  "No object could be resolved (see the notes on the analysis-of-variance panel)."))
    keep <- c("object", "term", "level_a", "level_b", "difference",
              "n_judges_a", "n_judges_b", "effective_judges_a",
              "effective_judges_b", "se", .app_wald_columns(r$sizes),
              "p_adj")
    d <- r$sizes[, intersect(keep, names(r$sizes)), drop = FALSE]
    # no significant/practical flags: the difference turns red at 0.5 logits,
    # the adjusted p at the run's alpha
    dt <- style_mag_red(num_dt(d), d, "difference", 0.5)
    style_lo_red(dt, d, "p_adj", bdif_shown_alpha())
  }, code = function()
    paste0(bdif_code_grp(), "\n",
           sprintf(paste0('btl_dif(bt, factors, effects = %s, ',
                          'p_adjust = "holm", alpha = %s)$sizes'),
                   qstr(bdif_shown_effects()), bdif_shown_alpha())))
  output$bdif_notes <- renderUI({
    r <- bdif_res()
    if (is.null(r) || !length(r$notes)) return(NULL)
    sprintf("Note. %s.", paste(r$notes, collapse = "; "))
  })
  # the (object, term) whose by-group curve is shown: the selected
  # analysis-of-variance row (master-detail, as the Rasch DIF page drives its
  # ICC); the overlay groups judges by that term's factors
  sel_bdif <- reactive({
    r <- bdif_res(); req(!is.null(r))
    i <- input$bdif_anova_tbl_rows_selected
    selected <- .app_selected_row(i, nrow(r$summary))
    row <- r$summary[selected, ]
    list(object = row$object, term = row$term,
         row = selected)
  })
  register_plot("bdif_occ", function() {
    b <- bfit(); sb <- sel_bdif(); req(sb$object %in% b$objects$object)
    # surface the builder's own message (missing judge column, unchosen
    # factor, …) instead of one generic hint
    grp <- tryCatch(bdif_term_group(sb$term, sb$row), error = function(e) e)
    if (inherits(grp, "error"))
      validate(need(FALSE, conditionMessage(grp)))
    plot_btl_icc(b, sb$object, group = grp)
  }, w = 8, h = 5.5, code = function() {
    sb <- sel_bdif(); r <- bdif_res()
    maps <- if (!is.null(r)) r$bootstrap_design$factors else NULL
    vars <- if (is.list(maps) && length(maps))
      bdif_term_vars(sb$term, maps, sb$row) else
      strsplit(sb$term, ":", fixed = TRUE)[[1]]
    # setNames keeps the judge names, so the emitted grp is a judge -> cell map
    # (as bdif_term_group builds it) rather than an unnamed vector
    paste0(bdif_code_grp(), "\n",
           sprintf('grp <- setNames(as.character(rasch:::.factor_cells(as.data.frame(factors[%s], check.names = FALSE), sep = ":")), names(factors[[1]]))\n',
                   qvec(vars)),
           sprintf('plot_btl_icc(bt, %s, group = grp)',
                   qstr(sb$object %||% "")))
  })
  register_hover_tip("bdif_occ", function() {
    b <- bfit(); sb <- sel_bdif(); req(sb$object %in% b$objects$object)
    grp <- tryCatch(bdif_term_group(sb$term, sb$row), error = function(e) NULL)
    req(!is.null(grp))
    bdif_occ_points(b, sb$object, grp)
  }, "location", "mean", function(np)
    sprintf("%s (%s) · location %.3f · mean %.2f",
            np$opponent[1], np$group[1], np$location[1], np$mean[1]))

  # ---------------------------------------------------------------- facets --
  facet_dat <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm"),
                  "Run a Multiple Ratings (MFRM) analysis to see facet results."))
    req(input$facet_sel %in% f$facet_spec)
    f$facet_effects[[input$facet_sel]]
  })
  # additive fits carry no interaction card; the footer says so (and how to
  # get one), so toggling the facet structure visibly changes this page
  output$facet_structure_note <- renderUI({
    f <- tryCatch(fit(), error = function(e) NULL)
    if (!inherits(f, "rasch_mfrm") || !is.null(f$interaction)) return(NULL)
    "Additive structure: no item-by-facet terms estimated (choose Interactive in the data roles to test them)."
  })
  register_table("facet_tbl", function() facet_dat(), function()
    # num_dt flags the facet fit residual beyond |2.5|, as on every model table
    num_dt(curate(facet_dat(), "facet", full = isTRUE(input$facets_full))),
    code = function()
    sprintf('fit$facet_effects[[%s]]', qstr(input$facet_sel %||% "")))
  register_plot("facet_plot", function() {
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm"),
                  "Run a Multiple Ratings (MFRM) analysis to see facet results."))
    req(input$facet_sel %in% f$facet_spec)
    plot_facets(f, input$facet_sel)
  }, code = function()
    sprintf('plot_facets(fit, %s)', qstr(input$facet_sel %||% "")))
  facet_int <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm") && !is.null(f$interaction),
                  "Choose the interactive facet structure in the data roles and re-estimate to test item-by-facet interactions."))
    f$interaction_effects
  })
  register_table("facet_int_tbl", function() facet_int(), function() {
    d <- facet_int()
    d$significant <- ifelse(d$significant %in% TRUE, "*", "")
    num_dt(d)
  }, code = function() "fit$interaction_effects")
  register_table("facet_int_omnibus", function() {
    f <- fit()
    validate(need(inherits(f, "rasch_mfrm") && !is.null(f$interaction_test),
                  "No interaction omnibus test is available."))
    f$interaction_test
  }, function() num_dt(fit()$interaction_test),
  code = function() "fit$interaction_test")

  # ----------------------------------------------------------------- equating --
  eq_ref <- reactive({
    a <- if (!is.null(input$eq_file))
      tryCatch(read.csv(input$eq_file$datapath, check.names = FALSE,
                        stringsAsFactors = FALSE), error = function(e) NULL)
    else restored_project_resources()[["eq_reference"]]
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
      validate(need(!is.null(input$eq_file) ||
                      !is.null(restored_project_resources()[["eq_reference"]]),
                    "Upload a reference calibration (item, location, se) to equate against."))
      eq_ref()
    }
  })
  eq_independent <- reactive(if (identical(input$eq_source, "kept"))
    isTRUE(input$eq_kept_independent) else isTRUE(input$eq_csv_independent))
  eq_reference_code <- function() {
    if (!identical(input$eq_source, "kept")) {
      if (!is.null(input$eq_file))
        return(sprintf(paste0("reference <- read.csv(%s, check.names = FALSE, ",
                              "stringsAsFactors = FALSE)"),
                       qstr(input$eq_file$name %||% "reference.csv")))
      return('reference <- project$resources[["eq_reference"]]')
    }
    lab <- input$eq_kept
    z <- kept_fit_code()[[lab]]
    if (is.null(z)) return("reference <- kept_fit  # fit retained in this app session")
    fit_code_block(z$code, z$value, "reference")
  }
  eq_res <- reactive(equate_tests(fit(), eq_reference(), shift = input$eq_shift,
                                  independent = eq_independent()))
  register_table("eq_tbl", function() eq_res()$table, function() {
    d <- curate(eq_res()$table, "equate", full = isTRUE(input$eq_full))
    if ("drift" %in% names(d)) d$drift <- ifelse(d$drift, "*", "")
    num_dt(d)
  }, code = function()
    paste0(eq_reference_code(), "\n",
      sprintf('eq <- equate_tests(fit, reference, shift = %s, independent = %s)\neq$table',
              qstr(input$eq_shift %||% "mean"),
              if (eq_independent()) "TRUE" else "FALSE")))
  register_plot("eq_plot", function()
    plot_equate(fit(), eq_reference(), shift = input$eq_shift,
                independent = eq_independent()),
    code = function()
      paste0(eq_reference_code(), "\n",
        sprintf('plot_equate(fit, reference, shift = %s, independent = %s)',
                qstr(input$eq_shift %||% "mean"),
                if (eq_independent()) "TRUE" else "FALSE")))
  # hover identification for the equating plot: only drifting items are
  # text-labelled (see plot_equate(), R/equating.R); eq_res()$table is the
  # exact table plot_equate() draws from (same fit/reference/shift), so
  # reusing it guarantees identical rows.
  register_hover_tip("eq_plot", function() eq_res()$table,
    "location_2", "location_1", function(np)
    sprintf("%s · reference %.3f · current %.3f",
            np$item[1], np$location_2[1], np$location_1[1]))
  output$dl_anchors <- downloadHandler(
    filename = function() format(Sys.time(), "rasch_anchors_%Y%m%d_%H%M.csv"),
    content = function(file) {
      f <- fit()
      thr <- f$thresholds
      write_csv_plain(data.frame(item = f$items$item[thr$item], k = thr$k,
                                 tau = thr$tau), file)
    })

  output$dl_calib <- downloadHandler(
    filename = function() format(Sys.time(), "rasch_calibration_%Y%m%d_%H%M.csv"),
    content = function(file) {
      f <- fit()
      write_csv_plain(data.frame(item = f$items$item, location = f$items$location,
                                 se = f$items$se, max = f$items$max), file)
    })

  # ---------------------------------------------------------------- frames --
  efrm_fit <- reactive({
    f <- fit()
    validate(need(inherits(f, "rasch_efrm"),
                  "Run an Extended Frames (EFRM) analysis to see results here."))
    f
  })
  register_table("frame_tbl", function() efrm_fit()$frames,
                 function() num_dt(curate(efrm_fit()$frames, "frames",
                                          full = isTRUE(input$frames_full))),
                 code = function() "fit$frames")
  efrm_phi_tbl <- reactive({
    f <- efrm_fit(); d <- f$phi_table
    ut <- f$efrm_vs_rasch$unit_tests
    if (!is.null(ut)) {
      ii <- match(paste0("log phi[", d$group, "]"), ut$parameter)
      d$log_unit <- log(d$phi)
      d$z <- ut$z[ii]; d$p <- ut$p[ii]; d$p_adj <- ut$p_adj[ii]
    }
    d
  })
  # The frame model constrains each item to one location across frames, so
  # the fit cannot test that constraint. frame_invariance() calibrates each
  # frame separately and compares -- the check belongs beside the units it
  # validates, since a couple of items behaving differently across frames
  # move a unit by more than its standard error.
  # Conditional inference covers item locations; bootstrap inference covers
  # both location and discrimination comparisons. Decisions use one combined
  # Holm family in either case.
  inv_se <- reactive(input$inv_se %||% "conditional")
  inv_reps <- reactive({
    x <- suppressWarnings(as.numeric(input$inv_boot %||% NA_real_))
    validate(need(length(x) == 1L && is.finite(x) && x == floor(x) &&
                    x >= 30L && x <= 1000L,
                  "Replicates must be one whole number from 30 to 1000."))
    as.integer(x)
  })
  inv_seed <- reactive({
    x <- suppressWarnings(as.numeric(input$inv_seed %||% NA_real_))
    validate(need(length(x) == 1L && is.finite(x) && x == floor(x) &&
                    x >= 0L && x <= .Machine$integer.max,
                  paste("Seed must be one non-negative whole number within",
                        "the integer range.")))
    as.integer(x)
  })
  efrm_invariance_base <- reactive({
    f <- efrm_fit()
    saved <- restored_invariance()
    if (!is.null(saved) &&
        !inherits(tryCatch(.validate_frame_invariance(saved, f),
                           error = function(e) e), "error") &&
        identical(saved$se_method, inv_se()) &&
        (inv_se() != "bootstrap" ||
          (identical(saved$boot_reps, inv_reps()) &&
           identical(saved$seed, inv_seed())))) return(saved)
    tryCatch(withProgress(
      message = if (inv_se() == "bootstrap")
        "Resampling persons within person groups..." else "Calibrating frames...",
      value = 0.5,
      frame_invariance(f, adjust = "holm",
                       se_method = inv_se(), boot_reps = inv_reps(),
                       seed = inv_seed())),
             error = function(e) conditionMessage(e))
  })
  efrm_invariance <- reactive({
    z <- efrm_invariance_base()
    if (is.character(z)) return(z)
    z
  })
  inv_pcol <- function() "p_adj"
  inv_or_note <- function(part) {
    z <- efrm_invariance()
    if (is.character(z)) return(data.frame(note = z, stringsAsFactors = FALSE))
    d <- z[[part]]
    keep <- if (part == "summary") names(d) else if (part == "locations")
      c("set", "frame_1", "frame_2", "item", "difference", "se", "statistic",
        "p", "p_adj", "flagged")
    else if (identical(z$se_method, "bootstrap"))
      c("set", "frame_1", "frame_2", "item", "log_disc_ratio",
        "se_log_disc_ratio", "statistic", "p", "p_adj", "flagged",
        "disc_1", "disc_2", "disc_ratio", "disc_boundary")
    else c("set", "frame_1", "frame_2", "item", "infit_1", "infit_2",
           "infit_z", "p", "p_adj", "flagged", "disc_1", "disc_2",
           "disc_ratio", "disc_boundary")
    d[, intersect(keep, names(d)), drop = FALSE]
  }
  inv_dt <- function(part) {
    d <- inv_or_note(part)
    pc <- inv_pcol()
    if (pc %in% names(d)) style_lo_red(num_dt(d), d, pc, 0.05) else num_dt(d)
  }
  inv_code <- function(part) sprintf(
    paste0("inv <- frame_invariance(fit, adjust = \"%s\", ",
           "se_method = \"%s\"%s)\ninv$%s"),
    "holm", inv_se(),
    if (inv_se() == "bootstrap")
      sprintf(", boot_reps = %d, seed = %d",
              inv_reps(), inv_seed())
    else "", part)
  # the download carries the full table, per the register_table contract: the
  # curated screen version drops location_1 and location_2, which are the two
  # columns the card's own explainer points at
  inv_full <- function(part) {
    z <- efrm_invariance()
    if (is.character(z)) data.frame(note = z, stringsAsFactors = FALSE)
    else as.data.frame(z[[part]])
  }
  register_table("frame_inv_summary_tbl", function() inv_full("summary"),
    function() num_dt(inv_full("summary")),
    code = function() inv_code("summary"))
  register_table("frame_inv_loc_tbl", function() inv_full("locations"),
    function() inv_dt("locations"), code = function() inv_code("locations"))
  register_table("frame_inv_disc_tbl", function() inv_full("discrimination"),
    function() inv_dt("discrimination"), code = function() inv_code("discrimination"))

  register_table("phi_tbl", function() efrm_phi_tbl(), function() {
    d <- efrm_phi_tbl()
    if ("p_adj" %in% names(d)) style_lo_red(num_dt(d), d, "p_adj", 0.05)
    else num_dt(d)
  }, code = function() paste(
    "d <- fit$phi_table",
    "ut <- fit$efrm_vs_rasch$unit_tests",
    "if (!is.null(ut)) {",
    "  i <- match(paste0(\"log phi[\", d$group, \"]\"), ut$parameter)",
    "  d$log_unit <- log(d$phi)",
    "  d$z <- ut$z[i]; d$p <- ut$p[i]; d$p_adj <- ut$p_adj[i]",
    "}", "d", sep = "\n"))
  # keep the fit's original set order: merge() sorts by the key, so it is
  # re-matched to fit$set_table$set
  efrm_alpha_tbl <- reactive({
    f <- efrm_fit()
    d <- merge(f$alpha_table, f$set_table[, c("set", "mu", "n_items")],
               by = "set", sort = FALSE)
    d <- d[stats::na.omit(match(f$set_table$set, d$set)), , drop = FALSE]
    ut <- f$efrm_vs_rasch$unit_tests
    if (!is.null(ut)) {
      ii <- match(paste0("log alpha[", d$set, "]"), ut$parameter)
      d$log_unit <- log(d$alpha)
      d$z <- ut$z[ii]; d$p <- ut$p[ii]; d$p_adj <- ut$p_adj[ii]
    }
    rownames(d) <- NULL
    d
  })
  register_table("alpha_tbl", function() efrm_alpha_tbl(),
                 function() {
                   d <- efrm_alpha_tbl()
                   if ("p_adj" %in% names(d))
                     style_lo_red(num_dt(d), d, "p_adj", 0.05)
                   else num_dt(d)
                 },
                 code = function() paste(
                   "d <- merge(fit$alpha_table,",
                   "  fit$set_table[, c(\"set\", \"mu\", \"n_items\")],",
                   "  by = \"set\", sort = FALSE)",
                   "d <- d[na.omit(match(fit$set_table$set, d$set)), , drop = FALSE]",
                   "ut <- fit$efrm_vs_rasch$unit_tests",
                   "if (!is.null(ut)) {",
                   "  i <- match(paste0(\"log alpha[\", d$set, \"]\"), ut$parameter)",
                   "  d$log_unit <- log(d$alpha)",
                   "  d$z <- ut$z[i]; d$p <- ut$p[i]; d$p_adj <- ut$p_adj[i]",
                   "}", "rownames(d) <- NULL", "d", sep = "\n"))
  register_plot("frame_plot", function() plot_frames(efrm_fit()),
                code = function() "plot_frames(fit)")
  register_plot("frame_icc", function() {
    f <- efrm_fit()
    req(input$frame_item %in% f$virtual_map$item)
    plot_icc_frames(f, input$frame_item)
  }, code = function()
    sprintf('plot_icc_frames(fit, %s)', qstr(input$frame_item %||% "")))
  efrm_cmp_tbl <- reactive({
    f <- efrm_fit(); x <- f$efrm_vs_rasch
    data.frame(
      model = c("Equal group units", "Group-dependent units"),
      loglik = c(x$ll_equal, x$ll_efrm),
      unit_parameters = c(0L, x$extra_parameters),
      two_delta_loglik = c(NA_real_, x$two_delta_ll),
      item_fit_residual_sd = c(NA_real_, f$item_fit_summary$sd),
      stringsAsFactors = FALSE)
  })
  register_table("efrm_cmp_tbl", function() efrm_cmp_tbl(),
                 function() num_dt(efrm_cmp_tbl()),
                 code = function() paste(
                   "x <- fit$efrm_vs_rasch",
                   "data.frame(",
                   "  model = c(\"Equal group units\", \"Group-dependent units\"),",
                   "  loglik = c(x$ll_equal, x$ll_efrm),",
                   "  unit_parameters = c(0L, x$extra_parameters),",
                   "  two_delta_loglik = c(NA_real_, x$two_delta_ll),",
                   "  item_fit_residual_sd = c(NA_real_, fit$item_fit_summary$sd)",
                   ")", sep = "\n"))
  register_table("efrm_omnibus_tbl", function() {
    x <- efrm_fit()$efrm_vs_rasch$unit_omnibus
    if (is.null(x)) data.frame() else x
  }, function() {
    x <- efrm_fit()$efrm_vs_rasch$unit_omnibus
    validate(need(!is.null(x) && nrow(x),
                  "No unit family has more than one level to test."))
    style_lo_red(num_dt(x), x, "p_adj", 0.05)
  }, code = function() "fit$efrm_vs_rasch$unit_omnibus")

  # ------------------------------------ BTL frames (paired-comparison EFRM) --
  # object -> set map for btl_efrm(): an uploaded CSV wins, otherwise infer
  # from the object-name prefix (the part before trailing digits/separators),
  # exactly as ef_setmap infers item sets for the Rasch EFRM above. Returns a
  # named list (set -> object names, the shape btl_efrm() wants) plus a short
  # description of the source for the frozen-run notes and code snippet.
  btlef_build_sets <- function(objs) {
    if (!is.null(input$btlef_sets_file)) {
      out <- read_frame_map(input$btlef_sets_file, "object", objs, "object")
      miss <- setdiff(objs, names(out))
      if (length(miss)) out[miss] <- "(rest)"
      return(list(sets = split(objs, out[objs]),
                  source = "the uploaded object-set CSV"))
    }
    saved <- restored_project_resources()$btlef_sets
    if (!is.null(saved)) {
      assigned <- unlist(saved, use.names = FALSE)
      if (setequal(assigned, objs) && !anyDuplicated(assigned))
        return(list(sets = saved, source = "the map embedded in the saved analysis"))
    }
    pref <- sub("[_. -]*[0-9]+$", "", objs)
    pref[pref == ""] <- "(rest)"
    list(sets = split(objs, pref), source = "object-name prefixes (trailing digits stripped)")
  }
  btlef_res <- reactiveVal(NULL)
  btlef_job <- reactiveVal(NULL)
  output$btlef_job_controls <- renderUI({
    if (is.null(btlef_job())) return(NULL)
    actionButton("cancel_btlef", "Cancel frame estimation",
                 icon = bs_icon("x-circle"),
                 class = "btn-outline-danger w-100 mt-2")
  })
  outputOptions(output, "btlef_job_controls", suspendWhenHidden = FALSE)

  store_btlef_result <- function(r, st) {
    r$run_panel_col <- st$panel_col
    r$run_judge_col <- st$judge_col
    r$run_set_source <- st$set_source
    r$run_se_method <- st$se_method
    r$run_boot_reps <- st$boot_reps
    r$run_workers <- st$workers
    r$run_seed <- st$seed
    r$run_maxit <- st$maxit
    r$run_tol <- st$tol
    r$run_object_a <- st$object_a
    r$run_object_b <- st$object_b
    r$run_winner <- st$winner
    btlef_res(r)
    clear_btl_analysis_steps()
    push_btl_analysis_step(
      "btl_frames",
      sprintf("Extended frames: %d set%s × %d panel%s",
              length(r$sets), if (length(r$sets) == 1L) "" else "s",
              length(r$panels), if (length(r$panels) == 1L) "" else "s"),
      r,
      details = list(panel = st$panel_col, set_source = st$set_source,
                     se_method = st$se_method, boot_reps = st$boot_reps,
                     workers = st$workers, seed = st$seed,
                     maxit = st$maxit, tol = st$tol),
      code = paste0(btlef_code_setup(), "\nbt <- frm"))
    bdif_res(NULL)
    bdif_meta(NULL)
    invisible(r)
  }

  cancel_btlef_job <- function() {
    st <- isolate(btlef_job())
    if (is.null(st)) return(invisible(FALSE))
    if (st$process$is_alive()) stop_efrm_process(st$process)
    close_efrm_job(st)
    btlef_job(NULL)
    invisible(TRUE)
  }

  observeEvent(input$cancel_btlef, {
    if (!cancel_btlef_job()) return(invisible(NULL))
    showNotification("Frame estimation cancelled.", type = "message",
                     duration = 5)
  })

  session$onSessionEnded(function() {
    st <- isolate(btlef_job())
    if (!is.null(st) && st$process$is_alive())
      stop_efrm_process(st$process)
    if (!is.null(st)) close_efrm_job(st)
  })

  observe({
    st <- btlef_job()
    if (is.null(st)) return()
    invalidateLater(250, session)
    if (file.exists(st$progress_file)) {
      z <- tryCatch(strsplit(readLines(st$progress_file, warn = FALSE)[1],
                             "\t", fixed = TRUE)[[1]],
                    error = function(e) character(0))
      if (length(z) == 3L) {
        current <- suppressWarnings(as.numeric(z[2]))
        total <- suppressWarnings(as.numeric(z[3]))
        stage <- z[1]
        fraction <- if (is.finite(current) && is.finite(total) && total > 0)
          pmin(pmax(current / total, 0), 1) else 0
        value <- switch(stage,
          "two-stage fit" = 0.05 + 0.10 * fraction,
          "judge bootstrap" = 0.15 + 0.80 * fraction,
          "parametric bootstrap" = 0.15 + 0.80 * fraction,
          "finalising" = 0.99, 0.02)
        detail <- if (stage %in% c("judge bootstrap",
                                   "parametric bootstrap") && total > 0)
          sprintf("%s: %d of %d", stage, round(current), round(total)) else stage
        st$progress$set(value = value, detail = detail)
      }
    }
    if (st$process$is_alive()) return()

    result <- tryCatch(st$process$get_result(), error = function(e) e)
    close_efrm_job(st)
    btlef_job(NULL)
    if (!background_job_is_current(st)) {
      showNotification(paste(
        "The completed frame result was not used because the active data or",
        "Comparative Judgement analysis changed while it was running."),
        type = "warning", duration = 8)
      return()
    }
    if (inherits(result, "error")) {
      showNotification(paste("Frame estimation failed:",
                             conditionMessage(result)),
                       type = "error", duration = 10)
      return()
    }
    if (length(result$warnings))
      showNotification(paste(result$warnings, collapse = "\n"),
                       type = "warning", duration = 10)
    store_btlef_result(result$value, st)
  })

  observeEvent(input$btlef_run, {
    if (!is.null(boot_job())) {
      showNotification("A fit bootstrap is already running. Cancel it before estimating frames.",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    if (!is.null(efrm_job())) {
      showNotification("An EFRM estimation is already running. Cancel it before estimating frames.",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    if (!is.null(btlef_job())) {
      showNotification("A frame estimation is already running. Cancel it before starting another.",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    req(btl_fit())
    base <- btl_fit()
    frame_source <- tryCatch(.app_btl_frame_source(base, raw_data()),
                             error = function(e) e)
    if (inherits(frame_source, "error")) {
      showNotification(conditionMessage(frame_source), type = "error",
                       duration = 10)
      return(invisible(NULL))
    }
    if (inherits(base, "rasch_btl_explanatory")) {
      showNotification(paste(
        "Explanatory object effects and frame units are not currently",
        "estimated in one model. Fit the free comparative judgement",
        "calibration before adding frames."),
        type = "error", duration = 10)
      return()
    }
    if (!is.null(base$anchors)) {
      showNotification(paste(
        "The current comparison fit uses external anchors.",
        "The frame model does not yet accept anchored object locations;",
        "refit without anchors before adding frames."),
        type = "error", duration = 10)
      return()
    }
    if (!is.null(frame_source$count)) {
      showNotification(paste(
        "The current comparison fit uses aggregated count weights.",
        "The frame model requires one comparison per row; expand the counts",
        "before adding frames."), type = "error", duration = 10)
      return()
    }
    graded_btl <- !is.null(frame_source$response) ||
      !is.null(frame_source$margin)
    if (graded_btl) {
      showNotification(paste(
        "The current comparison fit has ordered response categories.",
        "The frame model currently fits dichotomous winners only."),
        type = "error", duration = 10)
      return()
    }
    if (identical(frame_source$ties, "half")) {
      showNotification(paste(
        "The current comparison fit assigns half a win to ties.",
        "The frame model can only drop ties; choose Drop before adding frames."),
        type = "error", duration = 10)
      return()
    }
    if (!is.null(base$dependence)) {
      showNotification(paste(
        "The current comparison fit includes order or position effects.",
        "The frame model does not yet estimate those effects jointly; remove",
        "the order/position terms before adding frames."),
        type = "error", duration = 10)
      return()
    }
    df <- frame_source$data
    jc <- frame_source$judge
    if (is.null(jc) || identical(jc, NONE) || !jc %in% names(df)) {
      showNotification("Paired-comparison frames need the judge column nominated on the Data page.",
                       type = "error", duration = 10)
      btlef_res(NULL); return()
    }
    pcol <- input$btlef_panel
    if (is.null(pcol) || !length(pcol) || !nzchar(pcol) || !pcol %in% names(df)) {
      showNotification("Choose a judge-panel column in the sidebar.",
                       type = "error", duration = 8)
      btlef_res(NULL); return()
    }
    if (is.null(frame_source$winner)) {
      showNotification(paste("Paired-comparison frames fit dichotomous winner data only;",
                             "nominate a Winner column on the Data page (not a polytomous response)."),
                       type = "error", duration = 10)
      btlef_res(NULL); return()
    }
    jd <- as.character(df[[jc]]); first <- !duplicated(jd)
    # the panel column must be constant within judge (a panel attribute, not
    # a per-comparison one) -- the same check bdif_factor_maps makes for a
    # single judge factor
    if (!all(tapply(as.character(df[[pcol]]), jd,
                    function(v) length(unique(v)) == 1L))) {
      showNotification(paste0('"', pcol, '" is not constant within judge; choose a ',
                              "judge-panel (not per-comparison) column."),
                       type = "error", duration = 10)
      btlef_res(NULL); return()
    }
    panel_map <- setNames(as.character(df[[pcol]])[first], jd[first])
    objs <- base$objects$object
    sm <- btlef_build_sets(objs)
    se_method <- input$btlef_se %||% "judge_bootstrap"
    boot_reps_raw <- suppressWarnings(as.numeric(input$btlef_boot))
    if (length(boot_reps_raw) != 1L || !is.finite(boot_reps_raw) ||
        boot_reps_raw != floor(boot_reps_raw) || boot_reps_raw < 50L ||
        boot_reps_raw > 999L) {
      showNotification(paste("Frame bootstrap replicates must be one whole",
                             "number from 50 to 999"),
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    boot_reps <- as.integer(boot_reps_raw)
    seed_raw <- suppressWarnings(as.numeric(input$btlef_seed))
    if (length(seed_raw) != 1L || !is.finite(seed_raw) || seed_raw < 0 ||
        seed_raw != floor(seed_raw) || seed_raw > .Machine$integer.max) {
      showNotification(paste("Frame bootstrap seed must be one non-negative",
                             "whole number within the integer range"),
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    seed <- as.integer(seed_raw)
    workers_raw <- suppressWarnings(as.numeric(input$btlef_workers))
    if (length(workers_raw) != 1L || !is.finite(workers_raw) ||
        workers_raw != floor(workers_raw) ||
        !workers_raw %in% .efrm_worker_values) {
      showNotification("Frame workers must be one of the available whole-number choices",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    workers <- if (identical(se_method, "judge_bootstrap"))
      as.integer(workers_raw) else 1L
    eo <- est_opts()

    fit_args <- list(
      data = df, object_a = frame_source$object_a,
      object_b = frame_source$object_b,
      winner = frame_source$winner, judge = jc, panels = panel_map,
      object_sets = sm$sets, se_method = se_method, boot_reps = boot_reps,
      workers = workers, seed = seed, maxit = eo$maxit, tol = eo$tol)
    progress_file <- tempfile("rasch-btlef-", fileext = ".progress")
    log_file <- tempfile("rasch-btlef-", fileext = ".log")
    progress_bar <- shiny::Progress$new(session, min = 0, max = 1)
    progress_bar$set(message = "Estimating paired-comparison frames",
                     detail = "two-stage fit", value = 0.02)
    process <- tryCatch(callr::r_bg(
      function(args, progress_path, source_dir) {
        if (!is.null(source_dir)) {
          if (!requireNamespace("pkgload", quietly = TRUE))
            stop("pkgload is needed for a source-tree background analysis")
          pkgload::load_all(dirname(source_dir), quiet = TRUE)
        }
        progress_fun <- function(stage, current, total) {
          writeLines(paste(stage, current, total, sep = "\t"), progress_path)
        }
        args$progress <- progress_fun
        warnings <- character(0)
        value <- withCallingHandlers(
          do.call(if (exists("btl_efrm", inherits = TRUE))
            get("btl_efrm", inherits = TRUE) else
              getExportedValue("rasch", "btl_efrm"), args),
          warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          })
        list(value = value, warnings = unique(warnings))
      },
      args = list(args = fit_args, progress_path = progress_file,
                  source_dir = .rasch_source_dir),
      libpath = .libPaths(), stdout = log_file, stderr = log_file,
      supervise = TRUE), error = function(e) e)
    if (inherits(process, "error")) {
      progress_bar$close()
      unlink(c(progress_file, log_file), force = TRUE)
      showNotification(paste("Frame estimation failed:",
                             conditionMessage(process)),
                       type = "error", duration = 10)
      return(invisible(NULL))
    }
    btlef_job(list(
      process = process, progress = progress_bar,
      progress_file = progress_file, log_file = log_file,
      context = isolate(analysis_context()), data = df,
      check_btl_fit = TRUE, base_fit = isolate(base),
      panel_col = pcol, judge_col = jc, set_source = sm$source,
      se_method = se_method, boot_reps = boot_reps,
      workers = workers, seed = seed,
      maxit = eo$maxit, tol = eo$tol,
      object_a = frame_source$object_a, object_b = frame_source$object_b,
      winner = frame_source$winner))
    invisible(NULL)
  })
  output$has_btlef <- reactive(!is.null(btlef_res()))
  outputOptions(output, "has_btlef", suspendWhenHidden = FALSE)
  output$btlef_multiset <- reactive({
    r <- btlef_res(); !is.null(r) && nrow(r$alpha_table) > 1L
  })
  outputOptions(output, "btlef_multiset", suspendWhenHidden = FALSE)
  # the DISPLAYED run's setup lines, for the reproducible-code footers: panels
  # built with setNames (as the bdif snippet builds its factors) and the
  # object sets recovered from the run's own objects table (frozen at run
  # time, so later sidebar edits cannot silently change what the snippet
  # reproduces); falls back to the live sidebar before any run
  btlef_code_setup <- function() {
    r <- btlef_res()
    pcol <- (if (!is.null(r)) r$run_panel_col else input$btlef_panel) %||% "panel"
    jc <- (if (!is.null(r)) r$run_judge_col else input$bt_judge) %||% "judge"
    se_m <- (if (!is.null(r)) r$run_se_method else input$btlef_se) %||%
      "judge_bootstrap"
    reps <- (if (!is.null(r)) r$run_boot_reps else input$btlef_boot) %||% 60
    workers <- (if (!is.null(r)) r$run_workers else input$btlef_workers) %||% 1
    seed <- (if (!is.null(r)) r$run_seed else input$btlef_seed) %||% 1
    maxit <- (if (!is.null(r)) r$run_maxit else est_opts()$maxit) %||% 60
    tol <- (if (!is.null(r)) r$run_tol else est_opts()$tol) %||% 1e-8
    oa <- (if (!is.null(r)) r$run_object_a else input$bt_a) %||% "object_a"
    ob <- (if (!is.null(r)) r$run_object_b else input$bt_b) %||% "object_b"
    wn <- (if (!is.null(r)) r$run_winner else input$bt_win) %||% "winner"
    sets_txt <- if (!is.null(r))
      paste(deparse(split(r$objects$object, r$objects$set)), collapse = "\n    ")
    else "list(...)  # inferred from object-name prefixes, or the uploaded CSV"
    paste0(
      "# dat contains the comparison data\n",
      sprintf("panels <- setNames(as.character(dat$%s), dat$%s)[!duplicated(dat$%s)]\n",
              bq(pcol), bq(jc), bq(jc)),
      sprintf("object_sets <- %s\n", sets_txt),
      sprintf(paste0('frm <- btl_efrm(dat, object_a = %s, object_b = %s, ',
                     'winner = %s,\n                 judge = %s, panels = panels, ',
                     'object_sets = object_sets,\n                 se_method = %s, ',
                     "boot_reps = %s, workers = %s, seed = %s, ",
                     "maxit = %s, tol = %s)"),
              qstr(oa), qstr(ob), qstr(wn), qstr(jc), qstr(se_m), reps,
              workers, seed, maxit, format(tol)))
  }
  btlef_code_call <- function(field) paste0(btlef_code_setup(), "\n", field)
  register_table("btlef_phi_tbl", function() {
    r <- btlef_res(); req(!is.null(r)); r$phi_table
  }, function() {
    r <- btlef_res()
    validate(need(!is.null(r),
                  "Choose the judge-panel column in the sidebar and press Estimate frame units."))
    d <- r$phi_table
    style_lo_red(num_dt(d), d, "p_adj", 0.05)
  }, code = function() btlef_code_call("frm$phi_table"))
  # the alpha/kappa tables share the set column: merged for one readable
  # table, keeping the fit's own set order (merge() sorts by key otherwise)
  btlef_units_tbl <- reactive({
    r <- btlef_res(); req(!is.null(r))
    d <- merge(
      r$alpha_table[, c("set", "alpha", "se_log_alpha", "p_adj",
                        "significant")],
      r$kappa_table[, c("set", "kappa", "se_kappa", "p_adj",
                        "significant")], by = "set", sort = FALSE)
    names(d)[names(d) == "p_adj.x"] <- "p_adj_alpha"
    names(d)[names(d) == "significant.x"] <- "significant_alpha"
    names(d)[names(d) == "p_adj.y"] <- "p_adj_kappa"
    names(d)[names(d) == "significant.y"] <- "significant_kappa"
    d <- d[stats::na.omit(match(r$alpha_table$set, d$set)), , drop = FALSE]
    rownames(d) <- NULL
    d
  })
  register_table("btlef_units_tbl", function() btlef_units_tbl(), function() {
    r <- btlef_res()
    validate(need(!is.null(r) && nrow(r$alpha_table) > 1L,
                  "A single object set: set units are not estimated (panel-units model)."))
    d <- btlef_units_tbl()
    dt <- style_lo_red(num_dt(d), d, "p_adj_alpha", 0.05)
    style_lo_red(dt, d, "p_adj_kappa", 0.05)
  }, code = function() btlef_code_call(paste0(
    'd <- merge(frm$alpha_table[, c("set", "alpha", "se_log_alpha", "p_adj", "significant")],\n',
    '           frm$kappa_table[, c("set", "kappa", "se_kappa", "p_adj", "significant")],\n',
    '           by = "set", sort = FALSE)\n',
    'names(d)[names(d) == "p_adj.x"] <- "p_adj_alpha"\n',
    'names(d)[names(d) == "significant.x"] <- "significant_alpha"\n',
    'names(d)[names(d) == "p_adj.y"] <- "p_adj_kappa"\n',
    'names(d)[names(d) == "significant.y"] <- "significant_kappa"\n',
    'd <- d[na.omit(match(frm$alpha_table$set, d$set)), , drop = FALSE]\n',
    'rownames(d) <- NULL\nd')))
  register_plot("btlef_units_plot", function() {
    r <- btlef_res(); req(!is.null(r))
    plot_btl_units(r)
  }, code = function() btlef_code_call("plot_btl_units(frm)"))
  register_table("btlef_frames_tbl", function() {
    r <- btlef_res(); req(!is.null(r)); r$frames
  }, function() {
    r <- btlef_res()
    validate(need(!is.null(r),
                  "Choose the judge-panel column in the sidebar and press Estimate frame units."))
    num_dt(r$frames)
  }, code = function() btlef_code_call("frm$frames"))
  btlef_cmp_tbl <- reactive({
    r <- btlef_res(); req(!is.null(r)); x <- r$equal_unit
    data.frame(
      model = c("Equal units", "Frame-dependent units"),
      loglik = c(x$loglik_single, x$loglik_frames),
      parameters = c(x$parameters_single, x$parameters_frames),
      two_delta_loglik = c(NA_real_, x$two_delta_ll),
      stringsAsFactors = FALSE)
  })
  register_table("btlef_cmp_tbl", function() btlef_cmp_tbl(),
                 function() num_dt(btlef_cmp_tbl()),
                 code = function() btlef_code_call(paste(
                   "x <- frm$equal_unit",
                   "data.frame(",
                   "  model = c(\"Equal units\", \"Frame-dependent units\"),",
                   "  loglik = c(x$loglik_single, x$loglik_frames),",
                   "  parameters = c(x$parameters_single, x$parameters_frames),",
                   "  two_delta_loglik = c(NA_real_, x$two_delta_ll)",
                   ")", sep = "\n")))
  register_table("btlef_omnibus_tbl", function() {
    r <- btlef_res(); req(!is.null(r)); r$unit_omnibus
  }, function() {
    r <- btlef_res(); req(!is.null(r)); d <- r$unit_omnibus
    validate(need(!is.null(d) && nrow(d),
                  "No unit family has more than one level to test."))
    style_lo_red(num_dt(d), d, "p_adj", 0.05)
  }, code = function() btlef_code_call("frm$unit_omnibus"))
  output$btlef_note <- renderUI({
    r <- btlef_res()
    if (is.null(r)) return(NULL)
    bits <- character(0)
    bits <- c(bits, r$se_note)
    if (length(r$notes)) bits <- c(bits, paste(r$notes, collapse = "; "))
    p(class = "text-muted small mb-0", paste0("Note. ", paste(bits, collapse = "; "), "."))
  })

  # --------------------------------------------------------- dimensionality --
  dim_subsets <- reactiveVal(NULL)
  dim_computed <- reactiveVal(NULL)
  observeEvent(input$dim_apply, {
    restored_subtest(NULL)
    dim_computed(NULL)
    dim_subsets(NULL)
    dm_res(NULL)
    s <- NULL
    if (length(input$dim_pos) >= 2 && length(input$dim_neg) >= 2) {
      if (length(intersect(input$dim_pos, input$dim_neg))) {
        showNotification("The two subsets must be disjoint.", type = "error")
        return()
      }
      s <- list(pos = input$dim_pos, neg = input$dim_neg)
    } else if (!length(input$dim_pos) && !length(input$dim_neg)) {
      s <- NULL
    } else {
      showNotification("Nominate at least two items in each subset (or leave both empty).",
                       type = "warning")
      return()
    }
    dim_subsets(s)
    B_raw <- suppressWarnings(as.numeric(input$dim_boot_B %||% NA_real_))
    workers_raw <- suppressWarnings(as.numeric(input$dim_workers %||% NA_real_))
    seed_raw <- suppressWarnings(as.numeric(input$dim_boot_seed %||% NA_real_))
    if (length(B_raw) != 1L || !is.finite(B_raw) || B_raw < 0L ||
        B_raw != floor(B_raw) || B_raw > .Machine$integer.max) {
      showNotification("Dimensionality replicates must be one non-negative whole number",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    if (length(workers_raw) != 1L || !is.finite(workers_raw) ||
        workers_raw != floor(workers_raw) ||
        !workers_raw %in% .efrm_worker_values) {
      showNotification("Dimensionality workers must be one of the available whole-number choices",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    if (length(seed_raw) != 1L || !is.finite(seed_raw) || seed_raw < 0L ||
        seed_raw != floor(seed_raw) || seed_raw > .Machine$integer.max) {
      showNotification("Dimensionality seed must be one non-negative whole number within the integer range",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    B <- as.integer(B_raw)
    workers <- as.integer(workers_raw)
    seed <- as.integer(seed_raw)
    f <- fit()
    value <- withProgress(
      message = if (is.finite(B) && B > 0L)
        "Calibrating the dimensionality test…" else
          "Running the dimensionality test…", value = 0.4,
      soft(if (is.null(s)) dimensionality_test(
        f, component = pca_k(), B = B, workers = workers, seed = seed)
      else dimensionality_test(f, items_positive = s$pos,
                               items_negative = s$neg, B = B,
                               workers = workers, seed = seed)))
    dim_computed(value)
    showNotification(if (is.null(s)) sprintf(
      "Ran the t-test on the automatic split (residual component %d).", pca_k())
      else "Ran the t-test on the nominated item subsets.", type = "message")
    # the magnitude table is computed from the subsets in force at ITS run;
    # a changed split makes it stale
    dm_res(NULL)
  })
  observeEvent(input$pca_component, {
    dm_res(NULL)
    if (!isTRUE(restoring_project())) {
      dim_computed(NULL)
      restored_subtest(NULL)
    }
  }, ignoreInit = TRUE)
  observeEvent(c(input$dim_boot_B, input$dim_workers, input$dim_boot_seed), {
    if (!isTRUE(restoring_project())) {
      dim_computed(NULL)
      restored_subtest(NULL)
    }
  }, ignoreInit = TRUE)
  # the residual principal component that, when no manual subsets are named,
  # defines the t-test default split
  pca_k <- reactive({
    k <- suppressWarnings(as.integer(input$pca_component %||% 1))
    if (is.na(k) || k < 1L) 1L else k
  })
  dim_res <- reactive({
    f <- fit()
    saved <- restored_subtest()
    if (!is.null(saved) &&
        !inherits(tryCatch(.validate_dimensionality_test(saved, f),
                           error = function(e) e), "error")) return(saved)
    value <- dim_computed()
    validate(need(!is.null(value),
                  "Choose item subsets or an automatic component, then press Run t-test."))
    value
  })
  output$dim_txt <- renderPrint({
    dt <- dim_res()
    if (!is.null(dt$note)) { cat(dt$note); return(invisible()) }
    cat(sprintf("Item split: %s\n", dt$split))
    cat(sprintf("First residual eigenvalue: %.3f\n", dt$first_eigenvalue))
    cat(sprintf("Significant person t-tests: %.1f%%  (Clopper-Pearson 95%% CI %.1f%% to %.1f%%, n = %d)\n",
                100 * dt$prop_significant, 100 * dt$ci[1], 100 * dt$ci[2], dt$n))
    cat(sprintf("Persons excluded (extreme on a subset): %d\n", dt$n_excluded_extreme))
    resolution_limited <- grepl(
      "resolution insufficient", dt$verdict_method %||% "", fixed = TRUE)
    verdict <- if (isTRUE(dt$multidimensional))
      "evidence against unidimensionality" else
      if (identical(dt$multidimensional, FALSE))
        "consistent with unidimensionality" else
          if (resolution_limited)
            "withheld because the bootstrap resolution is insufficient" else
            "withheld for the data-driven split"
    cat(sprintf("Verdict: %s\n", verdict))
    if (!is.null(dt$p_boot))
      cat(sprintf("Bootstrap p: %s (%d of %d replicates used)\n",
                  fmt_p(dt$p_boot), dt$bootstrap$B_used, dt$bootstrap$B))
    if (!is.null(dt$caution)) cat("Caution:", dt$caution, "\n")
    # a split chosen from the residuals is chosen to disagree: the binomial
    # rule is interpretable only for a split fixed in advance
    if (dt$split != "manual" && is.null(dt$p_boot))
      cat("Note: the split was chosen from the residuals; name the subsets by",
          "content or use bootstrap calibration for an inferential verdict.\n")
    if (isTRUE(resolution_limited))
      cat(sprintf(paste0(
        "Note: the smallest attainable bootstrap p is %.3f; increase B for ",
        "an inferential verdict at alpha %.3f.\n"),
        dt$bootstrap_resolution, dt$alpha))
    if (!is.null(dt$paired_t))
      cat(sprintf("Paired t-test of subset means: mean difference %.3f, t = %.2f (df %.0f), p = %s\n",
                  dt$paired_t$mean_difference, dt$paired_t$t,
                  dt$paired_t$df, fmt_p(dt$paired_t$p)))
    cat("\nSubset A items:\n ", paste(dt$items_positive, collapse = ", "), "\n")
    cat("Subset B items:\n ", paste(dt$items_negative, collapse = ", "), "\n")
  })
  register_code("dim", function() {
    s <- dim_subsets()
    B <- input$dim_boot_B %||% 0L
    workers <- suppressWarnings(as.numeric(input$dim_workers %||% 1L))
    seed <- input$dim_boot_seed %||% 1L
    extra <- if (length(B) == 1L && is.finite(B) && B > 0L) sprintf(
      ", B = %s, workers = %s, seed = %s", format(B, scientific = FALSE),
      format(workers, scientific = FALSE), format(seed, scientific = FALSE)) else ""
    if (is.null(s)) {
      k <- pca_k()
      return(if (k == 1L) sprintf("dimensionality_test(fit%s)", extra)
             else sprintf("dimensionality_test(fit, component = %d%s)", k, extra))
    }
    sprintf("dimensionality_test(fit, items_positive = %s,\n  items_negative = %s%s)",
            qvec(s$pos), qvec(s$neg), extra)
  })

  # magnitude of multidimensionality (Andrich 2016): needs every item in a
  # subset; the component split satisfies this by construction, manual subsets may not
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
        "Adjust subsets A and B (or leave both empty for the selected component's split) and re-run the t-test first."),
        type = "warning", duration = 10)
      return()
    }
    r <- withProgress(message = "Subtest re-analysis…", value = 0.4,
                      tryCatch(dimensionality_magnitude(f, list(s$pos, s$neg)),
                               error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Magnitude estimate failed:", conditionMessage(r)),
                       type = "error", duration = 10)
    else {
      r$run_subtests <- list(s$pos, s$neg)
      dm_res(r)
    }
  })
  output$dm_tbl <- renderDT({
    r <- dm_res()
    validate(need(!is.null(r),
                  "Press the button; the current subsets (manual, or the selected component's split) are combined into super-items and the two reliability calculations compared."))
    num_dt(r$table)
  })
  output$dm_tbl_csv <- downloadHandler(
    filename = function() "rasch_dimensionality_magnitude.csv",
    content = function(file) {
      r <- dm_res(); req(!is.null(r))
      write_csv_plain(r$table, file)
    })
  register_code("dm_tbl", function() {
    r <- dm_res(); req(!is.null(r))
    sprintf("dimensionality_magnitude(fit, list(%s, %s))$table",
            qvec(r$run_subtests[[1]]), qvec(r$run_subtests[[2]]))
  })
  scree_diag <- reactive({
    f <- fit()
    saved <- restored_dimensionality()
    if (!is.null(saved) &&
        !inherits(tryCatch(.validate_scree_result(saved, f),
                           error = function(e) e), "error")) return(saved)
    .scree_analysis(f,
      parallel = !inherits(f, c("rasch_efrm", "rasch_mfrm")), seed = 1L)
  })
  register_plot("scree", function() {
    plot_scree(fit(), result = scree_diag())
  }, code = function() {
    f <- fit()
    if (inherits(f, c("rasch_efrm", "rasch_mfrm")))
      "plot_scree(fit, parallel = FALSE, seed = 1)"
    else "plot_scree(fit, seed = 1)"
  })
  register_table("loadings_tbl", function() residual_pca(fit())$loadings_matrix,
                 function() {
    d <- residual_pca(fit())$loadings_matrix   # up to the first 10 components
    # roomier rows than the compact default, so the table reads at a weight
    # closer to the biplot beside it (matching the local dependence page)
    num_dt(d, paging = FALSE) |>
      formatStyle(names(d), fontSize = "1rem",
                  paddingTop = "7px", paddingBottom = "7px")
  }, code = function() "residual_pca(fit)$loadings_matrix")
  register_table("eigen_tbl", function() residual_pca(fit())$eigen_table,
                 function() {
    d <- residual_pca(fit())$eigen_table
    d$proportion <- 100 * d$proportion
    d$cumulative <- 100 * d$cumulative
    names(d)[match(c("proportion", "cumulative"), names(d))] <-
      c("Proportion %", "Cumulative %")
    num_dt(d) |> formatRound(c("Proportion %", "Cumulative %"), 1)
  }, code = function() "residual_pca(fit)$eigen_table")
  # biplot of the first two residual components; its card grows with the item
  # count so it stays level with the loadings table beside it (as on the local
  # dependence page)
  # an item-panel height that grows with the item count so a table and its
  # plot stay level, clamped so it never runs tiny or huge (~520px at 10 items)
  item_panel_px <- function(n) as.integer(max(400L, min(820L, 260L + 26L * n)))
  biplot_px <- function()
    item_panel_px(tryCatch(nrow(fit()$items), error = function(e) 10L))
  register_plot("pca_biplot", function() plot_pca_biplot(fit()),
                w = 7, h = 7, px = biplot_px,
                code = function() "plot_pca_biplot(fit)")

  # -------------------------------------------------------- local dependence --
  # the residual correlations as two paired matrices, each with its heatmap:
  # the Yen Q3 correlation (diagonal 1.00) and the adjusted Q3* (each Q3 less
  # the average off-diagonal, diagonal empty). Both tables show the lower
  # triangle only, so each item pair is read once, matching the heatmaps.
  rc_all <- reactive(residual_correlations(fit()))
  q3_mat <- reactive(rc_all()$matrix)
  q3s_mat <- reactive(rc_all()$star_matrix)
  # each matrix is flagged by its own rule at the shared threshold: the raw
  # matrix by |Q3| (Yen 1993), the adjusted matrix by Q3* (Christensen,
  # Makransky & Horton 2017). Flagging the raw matrix by the Q3* rule would
  # redden pairs whose own Q3 is nowhere near the cut.
  q3_flag <- reactive(abs(q3_mat()) >= ld_flag())
  q3s_flag <- reactive(q3s_mat() >= ld_flag())
  register_table("cormat_q3_tbl",
                 function() data.frame(item = rownames(q3_mat()), q3_mat(),
                                       check.names = FALSE),
                 function() tri_dt(q3_mat(), flagged = q3_flag()),
                 code = function() "residual_correlations(fit)$matrix")
  register_table("cormat_q3s_tbl",
                 function() data.frame(item = rownames(q3s_mat()), q3s_mat(),
                                       check.names = FALSE),
                 function() tri_dt(q3s_mat(), flagged = q3s_flag()),
                 code = function() "residual_correlations(fit)$star_matrix")
  # on-screen height grows with the item count (item_panel_px), so the triangle
  # keeps readable cells and stays level with the table beside it
  cormat_px <- function()
    item_panel_px(tryCatch(ncol(q3_mat()), error = function(e) 10L))
  register_plot("rcor_q3", function() plot_resid_cor(fit(), stat = "q3"),
                w = 8, h = 7, px = cormat_px,
                code = function() 'plot_resid_cor(fit, stat = "q3")')
  register_plot("rcor_q3s", function() plot_resid_cor(fit(), stat = "q3star"),
                w = 8, h = 7, px = cormat_px,
                code = function() 'plot_resid_cor(fit, stat = "q3star")')
  # hover identification: q3_mat()/q3s_mat() are the same matrices
  # plot_resid_cor() colours (R/plots.R), so the cell lookup is exact.
  register_hover_cormat("rcor_q3", q3_mat, "Q3")
  register_hover_cormat("rcor_q3s", q3s_mat, "Q3*")
  # the flag threshold (shared by both matrices' red highlighting); defaults to
  # a heuristic starting value; there is no universal Q3* critical value
  ld_flag <- reactive({
    fl <- input$ld_flag
    if (is.null(fl) || is.na(fl) || fl <= 0) 0.2 else fl
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
                       type = "error", duration = 10)
    else dep_res(r)
  })
  output$dep_txt <- renderPrint({
    r <- dep_res()
    validate(need(!is.null(r),
                  "Choose the dependent and independent items and press the button."))
    ref <- .app_wald_reference(r)
    cat(sprintf("Dependence of %s on %s: d = %.3f logits (se %.3f), %s = %.2f, p = %s\n",
                r$dependent, r$independent, r$d, r$se,
                ref$label, ref$value, fmt_p(r$p)))
  })
  output$dep_tbl <- renderDT({
    r <- dep_res()
    validate(need(!is.null(r), ""))
    num_dt(r$thresholds)
  })
  register_code("dep_tbl", function()
    {
      r <- dep_res(); req(!is.null(r))
      sprintf('dependence_magnitude(fit, dependent = %s, independent = %s)$thresholds',
              qstr(r$dependent), qstr(r$independent))
    })
  output$dep_tbl_csv <- downloadHandler(
    filename = function() "rasch_dependence_thresholds.csv",
    content = function(file) {
      r <- dep_res(); req(!is.null(r))
      write_csv_plain(r$thresholds, file)
    })

  # spread test against Andrich's least upper bounds
  spread_res <- reactiveVal(NULL)
  observeEvent(input$run_spread, {
    run_alpha <- input$spread_alpha %||% 0.05
    run_adjust <- "holm"
    r <- withProgress(message = "Principal-components refit…", value = 0.4,
                      tryCatch(spread_test(fit(), alpha = run_alpha,
                                           p_adjust = run_adjust),
                               error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Spread test failed:", conditionMessage(r)),
                       type = "error", duration = 10)
    else spread_res(r)
  })
  output$spread_tbl <- renderDT({
    r <- spread_res()
    validate(need(!is.null(r),
                  "Press the button after combining items into a superitem."))
    d <- r
    d$below_bound <- ifelse(is.na(d$below_bound), "",
                            ifelse(d$below_bound, "*", ""))
    d$dependent <- ifelse(is.na(d$dependent), "",
                          ifelse(d$dependent, "*", ""))
    style_lo_red(num_dt(d), d, "p_adj", attr(r, "alpha") %||% 0.05)
  })
  register_code("spread_tbl", function() {
    r <- spread_res()
    sprintf("spread_test(fit, alpha = %s, p_adjust = %s)",
            attr(r, "alpha") %||% 0.05,
            qstr(attr(r, "p_adjust") %||% "holm"))
  })
  output$spread_tbl_csv <- downloadHandler(
    filename = function() "rasch_spread_test.csv",
    content = function(file) {
      r <- spread_res(); req(!is.null(r))
      write_csv_plain(r, file)
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
    ch <- clamp01(input$guess_chance, 0.25)
    anc <- if (length(input$guess_anchors)) input$guess_anchors else NULL
    do_boot <- isTRUE(input$guess_bootstrap)
    boot_reps_raw <- suppressWarnings(as.numeric(
      input$guess_boot_reps %||% NA_real_))
    seed_raw <- suppressWarnings(as.numeric(input$guess_seed %||% NA_real_))
    if (do_boot &&
        (length(boot_reps_raw) != 1L || !is.finite(boot_reps_raw) ||
         boot_reps_raw < 50L || boot_reps_raw != floor(boot_reps_raw) ||
         boot_reps_raw > .Machine$integer.max)) {
      showNotification("Tailored bootstrap replicates must be one whole number of at least 50",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    if (do_boot &&
        (length(seed_raw) != 1L || !is.finite(seed_raw) || seed_raw < 0L ||
         seed_raw != floor(seed_raw) || seed_raw > .Machine$integer.max)) {
      showNotification("Tailored bootstrap seed must be one non-negative whole number within the integer range",
                       type = "warning", duration = 7)
      return(invisible(NULL))
    }
    boot_reps <- if (do_boot) as.integer(boot_reps_raw) else 999L
    boot_seed <- if (do_boot) as.integer(seed_raw) else NULL
    r <- withProgress(message = "Tailored analysis (three re-analyses)…",
                      value = 0.3,
                      tryCatch(withCallingHandlers(tailored_analysis(
                        f, chance = ch, anchor_items = anc,
                        se_method = if (do_boot) "bootstrap" else "none",
                        boot_reps = boot_reps, seed = boot_seed),
                        warning = function(w) {
                          showNotification(conditionMessage(w), type = "warning",
                                           duration = 15)
                          invokeRestart("muffleWarning")
                        }),
                               error = function(e) e))
    if (inherits(r, "error"))
      showNotification(paste("Tailored analysis failed:", conditionMessage(r)),
                       type = "error", duration = 10)
    else {
      guess_res(r)
    }
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
    if (identical(r$se_method, "bootstrap")) {
      up <- sum(r$table$significant %in% TRUE & r$table$shift > 0)
      cat(sprintf("Holm-adjusted positive shifts: %d of %d (%d usable bootstrap draws)\n",
                  up, nrow(r$table), r$boot_reps_used))
      if (2 * nrow(r$table) / (r$boot_reps_used + 1) >= 0.05)
        cat("Too few usable draws to detect significance after Holm adjustment; increase the bootstrap replicates.\n")
    } else {
      cat("Item shifts are descriptive; enable the person bootstrap for uncertainty.\n")
    }
  })
  guess_code_call <- function(result) {
    r <- guess_res(); req(!is.null(r))
    seed_arg <- if (identical(r$se_method, "bootstrap"))
      paste0(", seed = ", r$seed) else ""
    anchor_arg <- if (is.null(r$anchor_items_requested)) "NULL" else
      qvec(r$anchor_items_requested)
    paste0(sprintf(paste0("ta <- tailored_analysis(fit, chance = %s, ",
                          "anchor_items = %s, se_method = %s, boot_reps = %d%s)"),
                   r$chance, anchor_arg, qstr(r$se_method),
                   if (identical(r$se_method, "bootstrap")) r$boot_reps else
                     999L, seed_arg), "\n", result)
  }
  register_table("guess_tbl", function() {
    r <- guess_res(); req(!is.null(r)); r$table
  }, function() {
    validate(need(max(fit()$m) == 1L,
                  "Run a dichotomous (multiple-choice) analysis to use the tailored guessing procedure."))
    r <- guess_res()
    validate(need(!is.null(r), "Run the tailored analysis to see the comparison."))
    num_dt(r$table)
  }, code = function() guess_code_call("ta$table"))
  register_plot("guess_plot", function() {
    validate(need(max(fit()$m) == 1L,
                  "Run a dichotomous (multiple-choice) analysis to use the tailored guessing procedure."))
    r <- guess_res()
    validate(need(!is.null(r), "Run the tailored analysis to see the equating plot."))
    plot_equate(r$tailored, r$origin_equated, shift = "none",
                independent = FALSE)
  }, code = function() guess_code_call(
    'plot_equate(ta$tailored, ta$origin_equated, shift = "none", independent = FALSE)'))
  # hover identification: same equate_tests() table plot_equate() computes
  # internally from ta$tailored vs ta$origin_equated (shift = "none"); only
  # drifting items are text-labelled on the plot.
  guess_eq_res <- reactive({
    r <- guess_res(); req(!is.null(r))
    equate_tests(r$tailored, r$origin_equated, shift = "none",
                 independent = FALSE)
  })
  register_hover_tip("guess_plot", function() guess_eq_res()$table,
    "location_2", "location_1", function(np)
    sprintf("%s · reference %.3f · current %.3f",
            np$item[1], np$location_2[1], np$location_1[1]))

  # ---------------------------------------------------------------- compare --
  # kept fits hold either family (Rasch or paired-comparison/BTL); "keep
  # current fit" keeps whichever is active (btl_fit() when a BTL analysis is
  # current, the Rasch fit() otherwise -- the two are mutually exclusive, see
  # the run observer above)
  kept_fits <- reactiveVal(list())
  observeEvent(input$keep_fit, {
    bf <- if (!is.null(btl_fit())) bfit() else NULL
    f <- if (!is.null(bf)) bf else fit()
    k <- kept_fits()
    lab <- sprintf("%d_%s", length(k) + 1L,
                   if (inherits(f, "rasch_btl"))
                     if (inherits(f, "rasch_btl_efrm")) "CJ_frames"
                     else paste0("CJ_", f$thr_structure) else f$model)
    k[[lab]] <- f
    kept_fits(k)
    kc <- kept_fit_code()
    kc[[lab]] <- list(code = current_rcode(),
                      value = if (inherits(f, "rasch_btl")) "bt" else "fit")
    kept_fit_code(kc)
    updateSelectizeInput(session, "cmp_ref", choices = names(k),
                         selected = if (!is.null(input$cmp_ref) &&
                                        input$cmp_ref %in% names(k))
                           input$cmp_ref else names(k)[1])
    # the Equating "kept fit" reference is the Rasch common-item route only
    # (paired-comparison equating uses its own banked-CSV path on the
    # Equating page); offer only the Rasch-family kept fits there
    rasch_k <- names(k)[vapply(k, .app_equating_fit, logical(1))]
    updateSelectizeInput(session, "eq_kept", choices = rasch_k,
                         selected = if (!is.null(input$eq_kept) &&
                                        input$eq_kept %in% rasch_k)
                           input$eq_kept else rasch_k[length(rasch_k)])
    showNotification(sprintf("Kept '%s' (%d fit(s) held).", lab, length(k)),
                     type = "message", duration = 5)
  })
  observeEvent(input$clear_fits, {
    kept_fits(list())
    kept_fit_code(list())
    updateSelectizeInput(session, "cmp_ref", choices = character(0),
                         selected = character(0))
    updateSelectizeInput(session, "eq_kept", choices = character(0),
                         selected = character(0))
    showNotification("Cleared kept fits.", type = "message", duration = 4)
  })
  # compare_fits() refuses a mixture of Rasch and BTL fits ("not a
  # mixture": their likelihoods are over different data), so the kept fits
  # are grouped by family around whichever one the reference belongs to;
  # kept fits of the other family are set aside (n_other) rather than shown
  kept_fit_code <- reactiveVal(list())

  # A new source is not merely a pending set of controls: it changes the
  # observations to which every fit and downstream result refers.  Drop the
  # active analysis before the new source can be displayed, cancel workers,
  # and advance the context token so a worker that races the change cannot
  # repopulate an old result. Kept fits are deliberate comparison snapshots,
  # so they remain available even though the active fit is cleared. Project
  # restoration is the one exception: its
  # source inputs are assigned before the authenticated fit and results are
  # reinstated later in the same flush cycle.
  invalidate_source_results <- function(reason) {
    if (isTRUE(restoring_project())) return(invisible(FALSE))
    try(cancel_efrm_job(), silent = TRUE)
    try(cancel_btlef_job(), silent = TRUE)
    try(cancel_boot_job(), silent = TRUE)
    advance_analysis_context()
    clear_analysis_steps()
    clear_btl_analysis_steps()
    fit_val(NULL)
    btl_fit(NULL)
    rcode_str(NULL)
    fitted_sim_gen(NULL)
    person_weight_state(NULL)
    restored_dimensionality(NULL)
    restored_subtest(NULL)
    restored_invariance(NULL)
    resolve_res(NULL)
    lr_res(NULL)
    rescore_res(NULL)
    contr_res(NULL)
    bdif_res(NULL)
    bdif_meta(NULL)
    btlef_res(NULL)
    dim_subsets(NULL)
    dim_computed(NULL)
    dm_res(NULL)
    dep_res(NULL)
    spread_res(NULL)
    guess_res(NULL)
    boot_val(NULL)
    dif_boot_val(NULL)
    showNotification(
      paste0(reason, ". The previous fit and its results were cleared; ",
             "run Estimate again."),
      type = "warning", duration = 8)
    invisible(TRUE)
  }
  fit_code_block <- function(code, value, target) {
    indented <- paste0("  ", gsub("\n", "\n  ", code %||% ""))
    paste0(target, " <- local({\n", indented, "\n  ", value, "\n})")
  }
  cmp_group <- reactive({
    k <- kept_fits()
    if (!length(k)) return(list(kk = k, ref = 1, n_other = 0L))
    ref <- if (!is.null(input$cmp_ref) && input$cmp_ref %in% names(k))
      input$cmp_ref else 1
    is_btl <- vapply(k, inherits, TRUE, what = "rasch_btl")
    ref_is_btl <- unname(is_btl[[ref]])
    list(kk = k[is_btl == ref_is_btl], ref = ref,
         n_other = sum(is_btl != ref_is_btl))
  })
  cmp_res <- reactive({
    g <- cmp_group()
    validate(need(length(g$kk) >= 2,
                  "Keep at least two fits of the same family (Rasch or paired-comparison) as the reference to compare."))
    as.data.frame(do.call(compare_fits, c(g$kk, list(reference = g$ref))))
  })
  output$cmp_family_note <- renderUI({
    g <- cmp_group()
    if (g$n_other == 0L) return(NULL)
    p(class = "text-muted small mb-2",
      sprintf("%d kept fit(s) of the other family are not shown (likelihoods are over different data).",
              g$n_other))
  })
  outputOptions(output, "cmp_family_note", suspendWhenHidden = FALSE)
  register_table("cmp_tbl", function() cmp_res(), function() {
    d <- cmp_res()
    if ("same_data" %in% names(d))
      d$same_data <- ifelse(d$same_data, "yes", "no")
    num_dt(curate(d, "compare", full = isTRUE(input$cmp_full)),
          one_dp = c("eff_params", "cl_aic", "cl_bic"))
  }, code = function() {
    g <- cmp_group(); kc <- kept_fit_code()
    labs <- names(g$kk)
    blocks <- vapply(seq_along(labs), function(i) {
      z <- kc[[labs[i]]]
      if (is.null(z))
        return(sprintf("fit_%d <- kept_fit_%d", i, i))
      fit_code_block(z$code, z$value, sprintf("fit_%d", i))
    }, character(1))
    args <- paste(sprintf("%s = fit_%d", qstr(labs), seq_along(labs)),
                  collapse = ", ")
    paste(c(blocks,
            sprintf("do.call(compare_fits, c(list(%s), list(reference = %s)))",
                    args, if (is.character(g$ref)) qstr(g$ref)
                          else as.character(g$ref))), collapse = "\n\n")
  })

  # ----------------------------------------------------------------- export --
  project_resources <- function(base_source = NULL) {
    current_model <- if (!is.null(btl_fit())) "btl" else
      if (inherits(fit_or_null(), "rasch_efrm")) "efrm" else
      if (inherits(fit_or_null(), "rasch_mfrm")) "mfrm" else "rasch"
    out <- list(
      anchors = tryCatch(anchors_in(), error = function(e) NULL),
      bt_anchors = tryCatch(bt_anchors_in(), error = function(e) NULL),
      key = tryCatch(key_in(), error = function(e) NULL),
      predictors = tryCatch(exp_predictors_raw(), error = function(e) NULL),
      eq_reference = if (!identical(current_model, "btl") &&
                           identical(input$eq_source %||% "csv", "csv"))
        tryCatch(eq_ref(), error = function(e) NULL) else NULL,
      bt_eq_bank = if (identical(current_model, "btl"))
        tryCatch(bt_eq_bank(), error = function(e) NULL) else NULL,
      wright_item_map = if (!identical(current_model, "btl") &&
                              identical(input$wright_item_panels,
                                        "uploaded"))
        tryCatch(wright_item_map_data(), error = function(e) NULL) else NULL,
      ef_setmap = if (identical(current_model, "efrm"))
        tryCatch(ef_setmap(), error = function(e) NULL) else NULL,
      btlef_sets = if (identical(current_model, "btl"))
        tryCatch(btlef_build_sets(as.character(bfit()$objects$object))$sets,
                 error = function(e) NULL) else NULL)
    # Base-calibration resources are frozen with the fit. File inputs are
    # live controls and may now point at a different upload; never let those
    # replace the metadata that the retained fit actually used.
    if (!is.null(base_source))
      for (nm in names(base_source$resources))
        out[nm] <- base_source$resources[nm]
    # A completed CJ frame model carries its own realised set allocation.
    # Preserve that allocation rather than a map selected after estimation.
    frm <- tryCatch(btlef_res(), error = function(e) NULL)
    if (identical(current_model, "btl") &&
        inherits(frm, "rasch_btl_efrm") &&
        all(c("object", "set") %in% names(frm$objects)))
      out["btlef_sets"] <- list(split(as.character(frm$objects$object),
                                      as.character(frm$objects$set)))
    out
  }

  project_state <- function() {
    base <- if (!is.null(btl_fit())) btl_fit() else analysis()
    source <- .app_fit_source(base)
    settings <- .collect_app_settings(input)
    if (!is.null(source))
      for (nm in names(source$settings)) settings[nm] <- source$settings[nm]
    simulation <- if (!is.null(source)) source$simulation else list(
      data = sim_data(), truth = sim_truth_val(), code = sim_code_val(),
      predictors = sim_predictors_val(), interactions = sim_interactions_val(),
      generation = sim_gen(), fitted_generation = fitted_sim_gen())
    .seal_app_project(list(
      format = "rasch-shiny-project", schema = 2L,
      package_version = tryCatch(as.character(utils::packageVersion("rasch")),
                                 error = function(e) NA_character_),
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      data = if (!is.null(source)) source$data else
        as.data.frame(raw_data(), check.names = FALSE),
      model_type = if (!is.null(btl_fit())) "btl" else
        if (inherits(fit_or_null(), "rasch_mfrm")) "mfrm" else
        if (inherits(fit_or_null(), "rasch_efrm")) "efrm" else "rasch",
      base_fit = base,
      rasch_steps = analysis_steps(), btl_steps = btl_analysis_steps(),
      rcode = rcode_str(), kept_fits = kept_fits(),
      kept_fit_code = kept_fit_code(), settings = settings,
      resources = project_resources(source), simulation = simulation,
      results = list(
      resolve = resolve_res(), lr = lr_res(), rescore = rescore_res(),
      contrasts = contr_res(),
      dif = if (is.null(btl_fit())) app_dif_res() else NULL,
      btl_dif = bdif_res(), btl_dif_meta = bdif_meta(),
      btl_frames = btlef_res(),
      dimension_subsets = dim_subsets(), dimension_magnitude = dm_res(),
      dependence = dep_res(), spread = spread_res(), guessing = guess_res(),
      person_weights = person_weight_state(),
      subtest = tryCatch(if (is.null(btl_fit())) dim_res() else NULL,
                         error = function(e) NULL),
      dimensionality = tryCatch(if (!is.null(btl_fit())) btl_dim()
        else scree_diag(), error = function(e) NULL),
      frame_invariance = tryCatch({
        z <- if (inherits(fit_or_null(), "rasch_efrm"))
          efrm_invariance() else NULL
        if (is.character(z)) NULL else z
      }, error = function(e) NULL),
      bootstrap = boot_val(), dif_bootstrap = {
        bv <- dif_boot_val()
        if (is.null(bv)) NULL else {
          f <- if (!is.null(btl_fit()))
            tryCatch(bfit(), error = function(e) NULL) else fit_or_null()
          d <- app_dif_res()
          if (is.null(f) || is.null(d))
            stop("the saved DIF bootstrap no longer has its active fit and primary DIF analysis; rerun it before saving")
          .validate_dif_bootstrap(bv$db, f, d)
          bv
        }
      })))
  }

  # exports carry the analysis as run in this session: the DIF model the
  # analyst configured and any bootstrap null, not default recomputations
  app_dif_res <- function(strict = FALSE) {
    f <- if (!is.null(btl_fit())) tryCatch(bfit(), error = function(e) NULL)
      else fit_or_null()
    if (is.null(f)) return(NULL)
    if (inherits(f, "rasch_btl")) return(bdif_res())
    # A frame-defining factor is not an eligible DIF factor. Its presence
    # alone must not prevent an otherwise valid frame-model export.
    eligible <- setdiff(names(f$factors), f$frame_group %||% character(0))
    if (!length(eligible)) return(NULL)
    result <- tryCatch(dif_res(), error = function(e) e)
    if (inherits(result, "error")) {
      if (isTRUE(strict))
        stop("The selected DIF analysis is unavailable: ",
             conditionMessage(result),
             ". Revise the DIF specification before exporting; the report ",
             "cannot substitute a different analysis.", call. = FALSE)
      # Project files may retain the chosen settings without a computed
      # result. Report writers interpret NULL as a request for defaults, so
      # exports must use the strict path above.
      return(NULL)
    }
    result
  }
  app_boot_res <- function() {
    bv <- boot_val()
    if (is.null(bv)) NULL else bv$bs
  }
  app_dif_boot_res <- function() {
    bv <- dif_boot_val()
    f <- if (!is.null(btl_fit())) tryCatch(bfit(), error = function(e) NULL)
      else fit_or_null()
    if (is.null(bv) || is.null(f)) return(NULL)
    d <- app_dif_res()
    if (is.null(d))
      stop("the DIF bootstrap no longer has its primary DIF analysis; rerun it before exporting")
    .validate_dif_bootstrap(bv$db, f, d)
    bv$db
  }
  app_dim_res <- function() {
    tryCatch(if (!is.null(btl_fit())) btl_dim() else scree_diag(),
             error = function(e) NULL)
  }
  app_subtest_res <- function(strict = FALSE) {
    if (!is.null(btl_fit())) return(NULL)
    f <- fit_or_null()
    # These designs have no app person-subset analysis. The report records
    # the package's deliberate refusal rather than blocking other outputs.
    if (is.null(f) || inherits(f, c("rasch_mfrm", "rasch_efrm")) ||
        .has_repeated_residual_units(f)) return(NULL)
    value <- tryCatch(dim_res(), error = function(e) e)
    if (!inherits(value, "error")) return(value)
    B <- suppressWarnings(as.numeric(input$dim_boot_B %||% 0))
    selected <- !is.null(dim_subsets()) || length(input$dim_pos) > 0L ||
      length(input$dim_neg) > 0L || pca_k() != 1L ||
      length(B) != 1L || !is.finite(B) || B != 0
    if (isTRUE(strict) && selected)
      stop("The selected dimensionality t-test is unavailable: ",
           conditionMessage(value),
           " Run a supported t-test before exporting; the report cannot ",
           "substitute its default item split or bootstrap settings.",
           call. = FALSE)
    # No optional t-test has been selected: the writer may use its usual
    # descriptive default. Project saving also retains unfinished settings.
    NULL
  }
  app_inv_res <- function() {
    f <- fit_or_null()
    if (inherits(f, "rasch_efrm")) {
      z <- efrm_invariance()
      if (!is.character(z)) z else NULL
    } else NULL
  }

  output$dl_project <- downloadHandler(
    filename = function()
      format(Sys.time(), "rasch_analysis_%Y%m%d_%H%M.rasch"),
    content = function(file)
      .save_app_project(project_state(), file))

  observeEvent(input$project_file, {
    p <- tryCatch(.read_app_project(input$project_file$datapath),
                  error = function(e) e)
    if (inherits(p, "error")) {
      showNotification(paste("The analysis could not be opened:",
                             conditionMessage(p)),
                       type = "error", duration = 10)
      return()
    }
    if (isTRUE(attr(p, "rasch_project_legacy"))) {
      dropped <- attr(p, "rasch_project_legacy_dropped") %||% character(0)
      showNotification(paste(
        "This analysis was saved in the earlier project format.",
        "It passed structural checks; save it again to use the current format.",
        if (length(dropped)) paste0(
          "Results that could not be authenticated were omitted: ",
          paste(dropped, collapse = ", "), ".") else ""),
        type = "warning", duration = 10)
    } else if (length(attr(p, "rasch_project_legacy_dropped") %||%
                      character(0))) {
      dropped <- attr(p, "rasch_project_legacy_dropped")
      showNotification(paste0(
        "Saved derived result(s) using superseded inference were omitted: ",
        paste(dropped, collapse = ", "),
        ". Recompute them before reporting those results."),
        type = "warning", duration = 10)
    }
    cancelled_job <- cancel_efrm_job() | cancel_btlef_job() |
      cancel_boot_job()
    advance_analysis_context()
    if (cancelled_job)
      showNotification(paste(
        "The running background analysis was cancelled before the saved",
        "analysis was opened."), type = "message", duration = 7)
    restoring_project(TRUE)
    withProgress(message = "Opening saved analysis…", value = 0.5, {
      restored_project_settings(p$settings %||% list())
      resources <- p$resources %||% list()
      # This runtime resource is reconstructed from the signed result, not
      # accepted as an independent persisted override. It leaves all source
      # columns intact, including factors also used in a comparison role.
      saved_bdif <- p$results[["btl_dif"]]
      saved_bdif_meta <- p$results[["btl_dif_meta"]]
      resources$btl_dif_factors <- if (!is.null(saved_bdif) &&
          !is.null(saved_bdif_meta)) list(
        factors = saved_bdif$bootstrap_design$factors,
        # The source calibration's role is restored on the Data page. The
        # DIF run may have used another column containing the same judge IDs;
        # the stored maps are already keyed to the fitted judges themselves.
        judge_col = saved_bdif_meta$fitted_judge_col %||%
          p$settings$bt_judge %||% saved_bdif_meta$judge_col,
        data = as.data.frame(p$data, check.names = FALSE),
        fit_signature = .fit_boot_signature(p$base_fit)) else NULL
      restored_project_resources(resources)
      sim_data(as.data.frame(p$data, check.names = FALSE))
      sim_truth_val(p$simulation$truth %||% NULL)
      sim_code_val(p$simulation$code %||% NULL)
      restored_project_name(input$project_file$name %||% "analysis.rasch")
      sim_predictors_val(p$simulation$predictors %||% NULL)
      sim_interactions_val(p$simulation$interactions %||% character(0))
      sim_gen(p$simulation$generation %||% 0L)
      fitted_sim_gen(p$simulation$fitted_generation %||% NULL)
      updateSelectInput(session, "demo_choice", selected = "none")
      updateRadioButtons(session, "model_type", selected = p$model_type)

      if (inherits(p$base_fit, "rasch_btl")) {
        fit_val(NULL); btl_fit(p$base_fit)
        analysis_steps(list()); btl_analysis_steps(p$btl_steps %||% list())
      } else {
        btl_fit(NULL); fit_val(p$base_fit)
        btl_analysis_steps(list()); analysis_steps(p$rasch_steps %||% list())
      }
      rcode_str(p$rcode %||% NULL)
      k <- p$kept_fits %||% list(); kept_fits(k)
      kept_fit_code(p$kept_fit_code %||% list())
      updateSelectizeInput(session, "cmp_ref", choices = names(k),
                           selected = if (length(k)) names(k)[1] else character(0))
      rasch_k <- names(k)[vapply(k, .app_equating_fit, logical(1))]
      saved_eq <- p$settings$eq_kept %||% character(0)
      updateSelectizeInput(session, "eq_kept", choices = rasch_k,
                           selected = if (length(saved_eq) == 1L &&
                                           !is.na(saved_eq) &&
                                           saved_eq %in% rasch_k) saved_eq
                             else if (length(rasch_k)) tail(rasch_k, 1)
                             else character(0))

      rr <- p$results %||% list()
      resolve_res(rr$resolve %||% NULL); lr_res(rr$lr %||% NULL)
      rescore_res(rr$rescore %||% NULL); contr_res(rr$contrasts %||% NULL)
      bdif_res(rr$btl_dif %||% NULL)
      bdif_meta(rr$btl_dif_meta %||% NULL)
      btlef_res(rr$btl_frames %||% NULL)
      dim_subsets(rr$dimension_subsets %||% NULL)
      restored_subtest(rr$subtest %||% NULL)
      dm_res(rr$dimension_magnitude %||% NULL)
      dep_res(rr$dependence %||% NULL); spread_res(rr$spread %||% NULL)
      guess_res(rr$guessing %||% NULL)
      person_weight_state(rr$person_weights %||% NULL)
      boot_val(rr$bootstrap %||% NULL)
      dif_boot_val(rr$dif_bootstrap %||% NULL)
      restored_dimensionality(rr$dimensionality %||% NULL)
      restored_invariance(rr$frame_invariance %||% NULL)
    })
    # raw_data() first refreshes the choices available to every role control;
    # applying the stored selections after that flush prevents the automatic
    # name guesses from overwriting the saved analysis configuration.
    session$onFlushed(function() {
      .restore_app_settings(session, restored_project_settings())
      # Projects saved before the subset selectors became ordinary saved
      # settings still carry the nominated split with their result. Reflect
      # that split in the controls as well as in the restored calculation.
      s <- dim_subsets()
      if (is.list(s)) {
        updateSelectizeInput(session, "dim_pos", selected = s$pos %||% character(0))
        updateSelectizeInput(session, "dim_neg", selected = s$neg %||% character(0))
      }
      restoring_project(FALSE)
    }, once = TRUE)
    showNotification(paste("Saved analysis opened. The active fit, its history,",
                           "data roles and estimation settings have been restored."),
                     type = "message", duration = 7)
    try(nav_select("nav", "p_summary", session = session), silent = TRUE)
  })

  # The Export page offers HTML, Word and PDF. The navbar shortcut remains a
  # one-click HTML report. Comparative Judgement uses the R Markdown report;
  # the established self-contained HTML writer remains the richer Rasch path.
  fit_with_app_results <- function(f) {
    z <- person_weight_state()
    if (!is.null(z)) attr(f, "report_person_weights") <- z
    f
  }
  report_content <- function(file, format = input$report_format %||% "html") {
    f <- if (!is.null(btl_fit())) bfit() else fit_or_null()
    if (is.null(f)) {
      showNotification("Run an analysis first, then download the report.",
                       type = "warning", duration = 8)
      stop("no fit to report")
    }
    f <- fit_with_app_results(f)
    tailored <- if (inherits(f, "rasch_btl")) NULL else guess_res()
    withProgress(message = paste("Building the", toupper(format), "report…"),
                 value = 0.4, {
      if (identical(format, "html") && !inherits(f, "rasch_btl"))
        report_html(f, file, dif = app_dif_res(strict = TRUE), bootstrap = app_boot_res(),
                    dif_bootstrap = app_dif_boot_res(),
                    dimensionality = app_dim_res(),
                    subtest = app_subtest_res(strict = TRUE), invariance = app_inv_res(),
                    tailored = tailored)
      else report_document(f, file, format = format,
                           dif = app_dif_res(strict = TRUE), bootstrap = app_boot_res(),
                           dif_bootstrap = app_dif_boot_res(),
                           dimensionality = app_dim_res(),
                           subtest = app_subtest_res(strict = TRUE),
                           invariance = app_inv_res(), tailored = tailored)
    })
  }
  output$dl_report <- downloadHandler(
    filename = function() paste0("rasch_report.", input$report_format %||% "html"),
    content = function(file) report_content(file, input$report_format %||% "html"))
  output$dl_report_nav <- downloadHandler(
    filename = function() "rasch_report.html",
    content = function(file) report_content(file, "html"))

  output$dl_zip <- downloadHandler(
    filename = function() format(Sys.time(), "rasch_results_%Y%m%d_%H%M.zip"),
    content = function(file) {
      f <- if (!is.null(btl_fit())) bfit() else fit()
      f <- fit_with_app_results(f)
      tailored <- if (inherits(f, "rasch_btl")) NULL else guess_res()
      tmp <- tempfile("rasch-results-")
      dir.create(tmp)
      on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
      withProgress(message = "Writing all tables and plots…", value = 0.4, {
        save_outputs(f, tmp,
                     formats = if (length(input$exp_formats)) input$exp_formats else "png",
                     item_plots = isTRUE(input$exp_items),
                     dif = app_dif_res(strict = TRUE), bootstrap = app_boot_res(),
                     dif_bootstrap = app_dif_boot_res(),
                     dimensionality = app_dim_res(),
                     subtest = app_subtest_res(strict = TRUE),
                     invariance = app_inv_res(), tailored = tailored)
      })
      writeLines(c(
        "rasch model results archive",
        "",
        paste("Included: fitted-model tables and plots, plus any compatible",
              "DIF, fit-bootstrap, DIF-bootstrap, dimensionality, person-subset",
              "dimensionality, frame-invariance, tailored guessing, and externally weighted secondary",
              "person measures passed to save_outputs()."),
        paste("Not included: app-only comparison, equating, rescoring, planned",
              "contrast, dependence-magnitude, spread, and",
              "other displays for which save_outputs() has no serialization contract."),
        paste("Use the table-specific CSV downloads for those displays. Save the",
              ".rasch analysis to retain the app state and reopen it later.")),
        file.path(tmp, "README.txt"))
      local({
        old <- setwd(tmp)
        on.exit(setwd(old))
        files <- list.files(".", recursive = TRUE, all.files = TRUE,
                            no.. = TRUE)
        unlink(file)
        status <- utils::zip(zipfile = file, files = files, flags = "-r9Xq")
        if (!isTRUE(as.integer(status) == 0L) || !file.exists(file)) {
          unlink(file)
          stop("the model-results archive could not be created")
        }
      })
    })
}

shinyApp(ui, server)
