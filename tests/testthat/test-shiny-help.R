test_that("every Shiny result card has one canonical explainer", {
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  help <- testthat::test_path("..", "..", "inst", "shiny", "help.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  if (!file.exists(help)) help <- system.file("shiny", "help.R", package = "rasch")

  h <- new.env(parent = baseenv())
  sys.source(help, envir = h)
  # duplicate names silently shadow later entries (APP_HELP["x"] returns
  # only the first match), so the registry must be unique by construction
  expect_false(anyDuplicated(names(h$APP_HELP)) > 0,
               info = paste("duplicate APP_HELP names:",
                            paste(unique(names(h$APP_HELP)[
                              duplicated(names(h$APP_HELP))]), collapse = ", ")))

  src <- paste(readLines(app, warn = FALSE), collapse = "\n")
  card_ids <- unique(unlist(lapply(c("tableCard", "plotCard", "statCard"),
    function(fun) {
      hit <- regmatches(src, gregexpr(
        sprintf("%s\\(\\\"[^\\\"]+", fun), src, perl = TRUE))[[1]]
      if (!length(hit) || identical(hit, character(0))) return(character(0))
      sub(sprintf("^%s\\(\\\"", fun), "", hit)
    })))

  expect_true(all(card_ids %in% names(h$APP_HELP)),
              info = paste("missing:",
                           paste(setdiff(card_ids, names(h$APP_HELP)),
                                 collapse = ", ")))

  hand_built_tables <- c(
    "chisq_int_tbl", "chisq_cat_tbl", "ctt_tbl", "rescore_tbl",
    "dif_posthoc_tbl", "dif_size_tbl", "resolve_tbl", "contr_tbl",
    "eq_tbl", "dm_tbl", "dep_tbl", "spread_tbl", "cmp_tbl",
    "sim_recovery_tbl", "sim_preview", "preview"
  )
  expect_true(all(hand_built_tables %in% names(h$APP_HELP)))
  expect_true(all(vapply(hand_built_tables, function(id)
    grepl(sprintf('app_help\\("%s"\\)', id), src), logical(1))),
    info = "a hand-built table is missing its information control")
  expect_true(all(nzchar(h$APP_HELP)))
})

test_that("Shiny explainers remain succinct", {
  help <- testthat::test_path("..", "..", "inst", "shiny", "help.R")
  if (!file.exists(help)) help <- system.file("shiny", "help.R", package = "rasch")
  h <- new.env(parent = baseenv())
  sys.source(help, envir = h)
  words <- lengths(strsplit(trimws(h$APP_HELP), "[[:space:]]+"))
  expect_true(max(words) <= 60,
              info = paste(names(which.max(words)), "has", max(words), "words"))
})

test_that("the app retains the agreed model labels and frame safeguards", {
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  src <- paste(readLines(app, warn = FALSE), collapse = "\n")
  expect_match(src, "output.active_btlef != true && output.btl_history_model != true",
               fixed = TRUE)
  expect_match(src, "Pair recommendations require a specified judge and comparison history",
               fixed = TRUE)

  expect_match(src, '"Multiple Ratings \\(MFRM\\)" = "mfrm"')
  expect_match(src, '"Extended Frames \\(EFRM\\)" = "efrm"')
  expect_false(grepl('"Many-facet \\(MFRM\\)" = "mfrm"', src))
  expect_false(grepl('"Rated / many-facet \\(MFRM\\)" = "mfrm"', src))

  # Adding BTL frames must not silently discard specifications unsupported by
  # btl_efrm(): each case is stopped before the frame fit is called.
  for (text in c("uses external anchors", "uses aggregated count weights",
                 "has ordered response categories", "assigns half a win to ties",
                 "includes order or position effects"))
    expect_match(src, text, fixed = TRUE)

  expect_match(src, 'if (inherits(f, "rasch_efrm")) f$item_arbitrary',
               fixed = TRUE)
  expect_match(src, "before <- app_item_estimates(analysis())", fixed = TRUE)
  expect_match(src, "after <- app_item_estimates(fit())", fixed = TRUE)
  expect_match(src, "rho = phi / alpha", fixed = TRUE)
  expect_false(grepl("rho = phi x alpha", src, fixed = TRUE))
  expect_match(src, 'style_lo_red(num_dt(d), d, "p_adj", 0.05)',
               fixed = TRUE)
  expect_match(src, '!all(c("object", "location") %in% names(a))',
               fixed = TRUE)
  expect_match(src, 'if (!"se" %in% names(a)) a$se <- NA_real_',
               fixed = TRUE)
  expect_match(src,
               'anchors$item <- ifelse(is.na(anchors$item), NA_character_',
               fixed = TRUE)
  expect_match(src,
               'bt_anchors$object <- ifelse(is.na(bt_anchors$object)',
               fixed = TRUE)
  expect_match(src, 'tmp <- tempfile("rasch-results-")', fixed = TRUE)
  expect_match(src,
               'on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)',
               fixed = TRUE)
  expect_false(grepl('paste0("rasch_", as.integer(Sys.time()))', src,
                     fixed = TRUE))
  expect_match(src, "save_outputs(f, tmp", fixed = TRUE)
  expect_match(src, "for which save_outputs() has no", fixed = TRUE)
  expect_match(src, "tailored = tailored", fixed = TRUE)
  expect_false(grepl("r$run_boot_reps, seed_arg", src, fixed = TRUE))
  expect_false(grepl("r$run_seed) else", src, fixed = TRUE))
  expect_match(src, 'attr(f, "report_person_weights") <- z', fixed = TRUE)
  expect_match(src, 'is.null(r$anchor_items_requested)', fixed = TRUE)
  expect_match(src, "externally weighted secondary", fixed = TRUE)
  expect_match(src, 'identical(r$shift_method, "unweighted")',
               fixed = TRUE)
  expect_match(src,
               'operative_p <- if (is.null(bv)) character(0) else',
               fixed = TRUE)
  expect_match(src,
               'c("chisq_p_boot_adj", "fit_resid_p_boot_adj")',
               fixed = TRUE)
  expect_false(grepl('p_bold = c("p", "p_adj")', src, fixed = TRUE))
  expect_match(src, ': judge-group DIF needs judge-constant factors',
               fixed = TRUE)
  expect_match(src, 'judge_factor <- function(x, label)', fixed = TRUE)
  expect_match(src, '"p_adj_kappa"', fixed = TRUE)
  expect_match(src, "maxit = eo$maxit, tol = eo$tol", fixed = TRUE)
  expect_match(src, 'input$spread_alpha %||% 0.05', fixed = TRUE)
  expect_match(src, 'run_adjust <- "holm"', fixed = TRUE)
  expect_match(src, "Holm-adjusted bootstrap chi-square p < .05",
               fixed = TRUE)
  expect_match(src, '"adjusted p unavailable; refit"', fixed = TRUE)
  expect_false(grepl("else f$dependence$p[r]", src, fixed = TRUE))
  expect_match(src, 'num_dt(d, p_bold = "p_adj")', fixed = TRUE)
  expect_false(grepl("dif_padj|spread_padj|inv_screen", src))
  expect_match(src, 'spread_test(fit, alpha = %s, p_adjust = %s)',
               fixed = TRUE)
  expect_match(src, 'callr::r_bg', fixed = TRUE)
  expect_match(src, 'input_task_button("boot_run", "Bootstrap the fit statistics"',
               fixed = TRUE)
  expect_match(src, 'identical(kind, "rasch")', fixed = TRUE)
  expect_match(src, 'item-fit bootstrap is unavailable when person IDs repeat',
               fixed = TRUE)
  expect_match(src, 'input_task_button("btl_boot_run", "Bootstrap the fit statistics"',
               fixed = TRUE)
  expect_match(src, 'getExportedValue("rasch", "fit_bootstrap")',
               fixed = TRUE)
  expect_match(src, '.classical_design_applicable <-', fixed = TRUE)
  expect_false(grepl('rasch:::.classical_design_applicable', src,
                     fixed = TRUE))
  expect_match(src, 'equating_rasch_on <- rasch_on && .app_equating_fit(f)',
               fixed = TRUE)
  expect_match(src, 'show("p_equating", equating_rasch_on || btl_on)',
               fixed = TRUE)
  expect_match(src, 'pkgload::load_all(dirname(source_dir), quiet = TRUE)',
               fixed = TRUE)
  expect_match(src, 'persons_with_boot <- function()', fixed = TRUE)
  expect_match(src, 'btl_boot_table <- function(which)', fixed = TRUE)
  expect_match(src, 'if (is.na(mis)) "Unavailable"', fixed = TRUE)
  expect_match(src, 'j <- match(fit$items$item, b$item)', fixed = TRUE)
  expect_false(grepl('merge(fit$items, bs$items', src, fixed = TRUE))
  expect_match(src, 'bootstrap = boot_val()', fixed = TRUE)
  expect_match(src, 'cancel_boot_job()', fixed = TRUE)
  expect_match(src, 'Cancel EFRM estimation', fixed = TRUE)
  expect_match(src, 'input$ef_workers', fixed = TRUE)
  expect_match(src, 'selected = max(.efrm_worker_values)', fixed = TRUE)
  expect_match(src, 'input$ef_seed', fixed = TRUE)
  expect_match(src, 'process$kill_tree()', fixed = TRUE)
  expect_match(src, 'paste0("workers = ", workers)', fixed = TRUE)
  expect_match(src, 'paste0("seed = ", seed)', fixed = TRUE)
  expect_match(src, '!workers_raw %in% .efrm_worker_values', fixed = TRUE)
  expect_false(grepl('as.integer(round(seed_raw))', src, fixed = TRUE))
  expect_match(src, 'Dimensionality seed must be one non-negative whole number',
               fixed = TRUE)
  expect_match(src, '0.10 + 0.84 * fraction', fixed = TRUE)
  expect_match(src, '"full person bootstrap" = 0.50 + 0.44 * fraction',
               fixed = TRUE)
  expect_false(grepl('"Hybrid (fast)"', src, fixed = TRUE))

  # Inference defaults and structural remedies are deliberately conservative.
  expect_match(src, '"Judge bootstrap (recommended)" = "judge_bootstrap"',
               fixed = TRUE)
  expect_match(src, '"Parametric bootstrap" = "bootstrap"', fixed = TRUE)
  expect_match(src, 'input$btlef_se %||% "judge_bootstrap"', fixed = TRUE)
  expect_match(src, 'input$btlef_workers', fixed = TRUE)
  expect_match(src, 'input$btlef_seed', fixed = TRUE)
  expect_match(src, 'identical(se_method, "judge_bootstrap")', fixed = TRUE)
  expect_match(src, 'Cancel frame estimation', fixed = TRUE)
  expect_match(src, 'btlef_job <- reactiveVal(NULL)', fixed = TRUE)
  expect_match(src, 'process$kill_tree()', fixed = TRUE)
  expect_match(src, '"boot_reps = %s, workers = %s, seed = %s, "',
               fixed = TRUE)
  expect_false(grepl('input$pc_id', src, fixed = TRUE))
  expect_match(src,
               '"exp_interactions", "bdif_factors", "dim_pos", "dim_neg"',
               fixed = TRUE)
  expect_match(src,
               'updateSelectizeInput(session, "dim_pos", selected = s$pos',
               fixed = TRUE)
  expect_match(src,
               'updateSelectizeInput(session, "dim_neg", selected = s$neg',
               fixed = TRUE)
  expect_match(src, 'unique(f$virtual_map$item)', fixed = TRUE)
  expect_match(src, "A location split cannot model", fixed = TRUE)
  expect_match(src,
    '!inherits(active_bf, c("rasch_btl_efrm", "rasch_btl_explanatory"))',
    fixed = TRUE)
  # A structural BTL refit changes bfit() through active_btl_step(), without
  # changing btl_fit(). Its fit-dependent judge-group DIF cache must therefore
  # be cleared by the active-step observer as well.
  active_clear <- regmatches(src, regexpr(
    'observeEvent\\(list\\(btl_fit\\(\\), active_btl_step\\(\\)\\), \\{[^}]+\\}',
    src, perl = TRUE))
  expect_length(active_clear, 1L)
  expect_match(active_clear, "bdif_res(NULL)", fixed = TRUE)
  expect_match(active_clear, "bdif_meta(NULL)", fixed = TRUE)
})

test_that("the app uses structurally stable responsive control layouts", {
  skip_if_not_installed("xml2")
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  old_cache <- Sys.getenv("R_USER_CACHE_DIR", unset = NA_character_)
  Sys.setenv(R_USER_CACHE_DIR = file.path(tempdir(), "r-cache"))
  on.exit(if (is.na(old_cache)) Sys.unsetenv("R_USER_CACHE_DIR") else
    Sys.setenv(R_USER_CACHE_DIR = old_cache), add = TRUE)
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  doc <- xml2::read_html(htmltools::renderTags(e$ui)$html)
  ids <- xml2::xml_attr(xml2::xml_find_all(doc, "//*[@id]"), "id")
  expect_false(anyDuplicated(ids) > 0L,
               info = paste("duplicate DOM ids:",
                            paste(unique(ids[duplicated(ids)]),
                                  collapse = ", ")))
  expect_false(any(grepl("^bslib-card-[0-9]+$", ids)),
               info = "a full-screen card still relies on bslib's random id")
  expect_false(any(grepl("^bslib-accordion-panel-[0-9]+$", ids)),
               info = "an accordion panel still relies on bslib's random id")
  expect_false(any(grepl("^bslib-accordion-[0-9]+$", ids)),
               info = "an accordion still relies on bslib's random id")
  bad <- xml2::xml_find_all(doc,
    "//*[contains(concat(' ',normalize-space(@class),' '),' bslib-sidebar-layout ')]")
  expect_length(bad, 0L)
  top <- xml2::xml_find_all(doc,
    "//body/div[contains(@class,'container-fluid')]/div[contains(@class,'tab-content')]")
  expect_length(top, 1L)
  # Accordion headings are themselves buttons. An information button inside
  # one is invalid HTML and causes Chromium to close the tab container early.
  expect_length(xml2::xml_find_all(doc, "//button//button"), 0L)
  values <- xml2::xml_attr(xml2::xml_find_all(top,
    "./div[contains(concat(' ',normalize-space(@class),' '),' tab-pane ')]"),
    "data-value")
  expect_true(all(c("p_data", "p_targeting", "p_dif", "p_frames", "p_export")
                  %in% values))
})

test_that("every analytical result has an R-code disclosure", {
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  tree <- parse(app)

  calls_named <- function(x, names) {
    out <- list()
    walk <- function(z) {
      if (is.call(z)) {
        fun <- tryCatch(as.character(z[[1]])[1], error = function(e) "")
        if (fun %in% names) out[[length(out) + 1L]] <<- z
        if (length(z) > 1L)
          for (i in 2:length(z)) if (length(z[[i]])) walk(z[[i]])
      } else if (is.expression(z) || is.pairlist(z)) {
        for (i in seq_along(z)) if (length(z[[i]])) walk(z[[i]])
      }
    }
    walk(x)
    out
  }
  first_string <- function(z)
    if (length(z) >= 2L && is.character(z[[2]]) && length(z[[2]]) == 1L)
      z[[2]] else NA_character_

  ui_calls <- calls_named(tree,
    c("plotCard", "tableCard", "statCard", "rcode_details"))
  server_calls <- calls_named(tree,
    c("register_plot", "register_table", "register_stat_box", "register_code"))
  ui_ids <- unique(na.omit(vapply(ui_calls, first_string, character(1))))
  server_ids <- unique(na.omit(vapply(server_calls, first_string, character(1))))
  expect_setequal(ui_ids, server_ids)

  managed <- calls_named(tree,
    c("register_plot", "register_table", "register_stat_box"))
  has_code <- vapply(managed, function(z) "code" %in% names(as.list(z)),
                     logical(1))
  expect_true(all(has_code),
    info = paste("missing code argument:",
      paste(vapply(managed[!has_code], first_string, character(1)),
            collapse = ", ")))
})

test_that("bundled app examples are exactly reconstructible", {
  examples <- testthat::test_path("..", "..", "inst", "shiny", "examples.R")
  if (!file.exists(examples))
    examples <- system.file("shiny", "examples.R", package = "rasch")
  e <- new.env(parent = asNamespace("rasch"))
  sys.source(examples, envir = e)
  expect_identical(.app_example_data("pcm"), e$.demo_data())
  workflow_data <- simulate_rasch(
    n_persons = 600, n_items = 12, model = "PCM", n_categories = 4,
    difficulty = c(-1.5, 1.5), disordered = "I04",
    dependence = list(pairs = list(c("I10", "I11")), strength = 1.3),
    dif = list(items = "I08", uniform = 0.8), n_groups = 3, seed = 17)
  expect_identical(e$.demo_data(), workflow_data)
  expect_identical(.app_example_data("dich"), e$.demo_dich())
  expect_identical(.app_example_data("rsm"), e$.demo_rsm())
  expect_identical(.app_example_data("mfrm"), e$.demo_mfrm())
  expect_identical(.app_example_data("efrm"), e$.demo_efrm())
  expect_identical(.app_example_data("btl"), e$.demo_btl())
  expect_identical(.app_example_data("pl"), e$.demo_pl())
  expect_identical(.app_example_data("cj"), e$.demo_cj())
  cj_data <- simulate_btl(
    n_objects = 8, n_judges = 48, reps_per_pair = 84,
    erratic_judges = 2 / 48,
    dependence = list(exposure = 0.7, carry_over = 0), seed = 2)
  cj_number <- as.integer(sub("^J", "", cj_data$judge))
  cj_data$panel <- factor(ifelse(cj_number %% 2L,
                                 "panel A", "panel B"))
  cj_data$experience <- factor(ifelse(cj_number <= 24,
                                      "experienced", "novice"))
  expect_identical(e$.demo_btl(), cj_data)

  expect_equal(as.vector(table(cj_data$panel, cj_data$experience)),
               rep(588, 4))
  cj_fit <- btl(cj_data, "object_a", "object_b", winner = "winner",
                judge = "judge", order = "order")
  dep <- setNames(seq_len(nrow(cj_fit$dependence)),
                  cj_fit$dependence$effect)
  expect_gt(cj_fit$dependence$estimate[dep["exposure"]], 0.4)
  expect_lt(cj_fit$dependence$p_adj[dep["exposure"]], 0.05)
  expect_lt(abs(cj_fit$dependence$estimate[dep["carry_over"]]), 0.2)
  expect_gt(cj_fit$dependence$p_adj[dep["carry_over"]], 0.05)
  cj_dif <- btl_dif(cj_fit, cj_data[c("panel", "experience")])
  expect_false(any(cj_dif$summary$uniform_DIF, na.rm = TRUE))
  expect_false(any(cj_dif$summary$nonuniform_DIF, na.rm = TRUE))
})

test_that("tailored bootstrap controls are validated and reproducible", {
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  txt <- paste(readLines(app, warn = FALSE), collapse = "\n")
  expect_match(txt, 'numericInput\\("guess_seed", "Random seed"',
               perl = TRUE)
  expect_match(txt, "boot_reps_raw != floor\\(boot_reps_raw\\)",
               perl = TRUE)
  expect_match(txt, "seed_raw != floor\\(seed_raw\\)", perl = TRUE)
  expect_match(txt, "boot_reps = boot_reps, seed = boot_seed", fixed = TRUE)
  expect_match(txt, 'paste0\\(", seed = ", r\\$seed\\)', perl = TRUE)
  expect_false(grepl("as.integer(input$guess_boot_reps", txt, fixed = TRUE))
  expect_match(txt, "else 999L", fixed = TRUE)
  expect_match(txt, "Too few usable draws to detect significance after Holm adjustment",
               fixed = TRUE)
  expect_match(txt, "showNotification(conditionMessage(w)", fixed = TRUE)
})
