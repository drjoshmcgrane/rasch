# The Joint calibration page: judgements of the items (paired comparisons,
# rankings, or both) calibrated jointly with the responses by rasch_cj(),
# the joint calibration anchoring the response analysis on every item but
# the last, and the page that shows the calibration frame by frame.

.joint_app <- function(notes = NULL) {
  for (package in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(package)
  env <- new.env(parent = globalenv())
  app_path <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app_path))
    app_path <- system.file("shiny", "app.R", package = "rasch")
  suppressWarnings(sys.source(app_path, envir = env))
  if (!is.null(notes))
    env$showNotification <- function(ui, ..., type = "default", duration = NULL) {
      notes$msgs <- c(notes$msgs, paste0("[", type, "] ",
                                         paste(as.character(ui),
                                               collapse = " ")))
      invisible("id")
    }
  env
}

.joint_notes <- function() {
  n <- new.env(parent = emptyenv())
  n$msgs <- character(0)
  n
}

.joint_demo <- function() {
  path <- testthat::test_path("..", "..", "inst", "shiny", "examples.R")
  if (!file.exists(path))
    path <- system.file("shiny", "examples.R", package = "rasch")
  e <- new.env(parent = asNamespace("rasch"))
  sys.source(path, envir = e)
  e$.demo_cj()
}

.joint_items <- sprintf("I%02d", 1:10)

# the sidebar state of a joint calibration run on the bundled example
.joint_inputs <- function(session, ...) {
  session$setInputs(demo_choice = "cj", model_type = "rasch",
                    rasch_calibration = "free", thr_structure = "pcm",
                    thr_mode = "free", anchor_type = "individual",
                    item_cols = .joint_items, id_col = "person_id",
                    factor_cols = character(0), cj_fix_comp = FALSE,
                    cj_fix_rank = FALSE, maxit = 200, tol = 1e-8,
                    ng_auto = TRUE, ...)
  session$flushReact()
}

test_that("the bundled joint calibration example carries its judgements", {
  d <- .joint_demo()
  expect_identical(names(d), c("person_id", .joint_items))
  expect_equal(nrow(d), 250L)
  comp <- attr(d, "comparisons", exact = TRUE)
  rank <- attr(d, "rankings", exact = TRUE)
  expect_identical(names(comp), c("object_a", "object_b", "winner"))
  expect_equal(nrow(comp), 300L)
  expect_true(all(comp$winner == comp$object_a | comp$winner == comp$object_b))
  expect_identical(names(rank), c("ranking", "item", "rank"))
  expect_equal(nrow(rank), 240L)
  expect_true(all(table(rank$ranking) == 4L))
  # the attributes survive the frame conversion the app applies at launch
  expect_identical(attr(as.data.frame(d, check.names = FALSE), "comparisons"),
                   comp)
  cj <- rasch_cj(d, comparisons = comp, rankings = rank, items = .joint_items)
  expect_s3_class(cj, "rasch_cj")
  expect_true(cj$converged)
  expect_identical(unname(cj$n), c(232L, 300L, 60L))
  expect_true(all(cj$units$estimated[cj$units$frame != "responses"]))
  # both judgement frames place I07 harder than the responses do
  tab <- cj$invariance$items
  hit <- tab$item[!is.na(tab$p_adj) & tab$p_adj < 0.05]
  expect_identical(sort(unique(hit)), "I07")
  expect_identical(sort(unique(tab$frame[tab$item %in% hit])),
                   c("comparisons", "rankings"))
  expect_lt(cj$invariance$lr$p, 0.05)
  # the hand-off leaves the last item free
  expect_identical(.joint_app()$.app_joint_free_item(cj), "I10")
})

test_that("the joint calibration helpers read a run the same way twice", {
  e <- .joint_app()
  d <- .joint_demo()
  expect_null(e$.app_joint_items(d, NULL, NULL, NULL))
  expect_identical(e$.app_joint_items(d, .joint_items[1:3], NULL, NULL),
                   .joint_items[1:3])
  expect_identical(e$.app_joint_items(d, NULL, "person_id", NULL),
                   .joint_items)
  expect_identical(e$.app_joint_units(TRUE, TRUE, NULL, NULL),
                   c(comparisons = NA_real_, rankings = NA_real_))
  expect_identical(e$.app_joint_units(TRUE, FALSE, data.frame(x = 1),
                                      data.frame(x = 1)),
                   c(comparisons = 1, rankings = NA_real_))
  cj <- rasch_cj(d, comparisons = attr(d, "comparisons"),
                 rankings = attr(d, "rankings"), items = .joint_items)
  anc <- e$.app_joint_anchors(cj)
  expect_identical(names(anc), c("item", "k", "tau"))
  expect_false("I10" %in% anc$item)
  expect_equal(nrow(anc), 9L)
  # a fit without app metadata, or one not anchored on judgements, has no
  # joint calibration to refit
  f <- rasch(d, id = "person_id", anchors = anc)
  expect_null(e$.app_joint_refit(f))
  attr(f, "rasch_app_source") <- list(
    data = d, settings = list(model_type = "rasch"),
    resources = list(anchors = anc), simulation = list())
  expect_null(e$.app_joint_refit(f))
})

test_that("a joint calibration runs from the sidebar, anchors the analysis and is reproduced by its code", {
  skip_on_cran()
  notes <- .joint_notes()
  e <- .joint_app(notes)
  d <- .joint_demo()
  path <- tempfile("joint", fileext = ".rasch")
  anchors <- tempfile("anchors", fileext = ".csv")
  on.exit(unlink(c(path, anchors)), add = TRUE)
  write.csv(data.frame(item = "I01", k = 1, tau = -1.5), anchors,
            row.names = FALSE)
  shiny::testServer(e$server, {
    session$flushReact()
    .joint_inputs(session)
    expect_identical(raw_data(), d)
    s <- cj_sources()
    expect_identical(s$origin, "example")
    expect_identical(s$comparisons, attr(d, "comparisons"))
    expect_identical(s$rankings, attr(d, "rankings"))
    expect_false(is.null(cj_sources_quiet()$comparisons))
    expect_false(is.null(cj_sources_quiet()$rankings))
    expect_match(output$cj_sources_note$html,
                 "from the example dataset: 300 paired comparisons and 60 rankings",
                 fixed = TRUE)
    expect_match(output$data_info$html,
                 "Judgements of the items are loaded as well (300 paired comparisons and 60 rankings)",
                 fixed = TRUE)
    expect_null(cj_fit())

    notes$msgs <- character(0)
    session$setInputs(run = 1L); session$flushReact()
    expect_false(any(grepl("^\\[error\\]", notes$msgs)), info = notes$msgs)
    cj <- cj_fit()
    expect_s3_class(cj, "rasch_cj")
    expect_true(cj$converged)
    f <- fit_val()
    expect_s3_class(f, "rasch")
    expect_null(btl_fit()); expect_null(pl_fit())
    # every item but the last is anchored at its joint location; the last
    # is estimated free
    anchored <- f$items$item != "I10"
    expect_equal(f$items$location[anchored], cj$items$location[anchored],
                 tolerance = 1e-6)
    expect_true(all(f$items$se[anchored] == 0))
    expect_gt(f$items$se[!anchored], 0)
    src <- attr(f, "rasch_app_source", exact = TRUE)
    expect_identical(src$resources$cj_comparisons, attr(d, "comparisons"))
    expect_identical(src$resources$cj_rankings, attr(d, "rankings"))
    expect_identical(src$resources$anchors, e$.app_joint_anchors(cj))
    expect_false(src$settings$cj_fix_comp)

    # the disclosed code reproduces the joint calibration and the fit
    code <- current_rcode()
    expect_match(code, 'dat <- rasch:::.app_example_data("cj")', fixed = TRUE)
    expect_match(code, 'cj_comparisons <- attr(dat, "comparisons")', fixed = TRUE)
    expect_match(code, 'cj_rankings <- attr(dat, "rankings")', fixed = TRUE)
    expect_match(code, "cj <- rasch_cj(dat,", fixed = TRUE)
    expect_match(code, "comparisons = cj_comparisons", fixed = TRUE)
    expect_match(code, "rankings = cj_rankings", fixed = TRUE)
    expect_false(grepl("units =", code, fixed = TRUE))
    expect_match(code,
                 'anchors <- cj$anchors[cj$anchors$item != "I10", , drop = FALSE]',
                 fixed = TRUE)
    expect_match(code, "anchors = anchors", fixed = TRUE)
    env <- new.env(parent = globalenv())
    eval(parse(text = code), envir = env)
    expect_equal(env$cj$loglik, cj$loglik)
    expect_equal(env$cj$items, cj$items)
    expect_equal(env$fit$items$location, f$items$location)

    # every card renders and names the object its code reads
    for (id in c("joint_boxes", "joint_fitsum_tbl", "joint_map",
                 "joint_items_tbl", "joint_invariance_tbl",
                 "joint_anchors_tbl"))
      expect_no_error(output[[id]])
    expect_identical(output$joint_boxes_code, "cj")
    expect_identical(output$joint_items_tbl_code, "cj$items")
    expect_identical(output$joint_map_code, "plot_cj(cj)")
    expect_identical(output$joint_invariance_tbl_code, "cj$invariance$items")
    expect_identical(output$joint_anchors_tbl_code, "anchors")
    expect_identical(output$joint_fitsum_tbl_code, "print(cj)")
    html <- output$joint_fitsum_tbl$html
    expect_match(html, "Unit of comparisons relative to responses", fixed = TRUE)
    expect_match(html, "Unit of rankings relative to responses", fixed = TRUE)
    expect_match(html, "every item but I10, which is left free", fixed = TRUE)
    expect_match(html, "I07", fixed = TRUE)
    expect_match(output$joint_boxes$html, "300 comparisons, unit se", fixed = TRUE)
    expect_match(output$joint_boxes$html, "60 rankings, unit se", fixed = TRUE)
    expect_match(output$joint_invariance_note$html,
                 "Objects differing between frames (Holm p &lt; .05): I07",
                 fixed = TRUE)
    expect_match(output$joint_anchors_note$html, "I10 is left free", fixed = TRUE)
    expect_true(!is.null(output$joint_fitsum_tbl_csv))
    expect_true(!is.null(output$joint_anchors_tbl_csv))
    expect_match(output$nav_status$html, "Dichotomous", fixed = TRUE)

    # the saved analysis keeps the judgements beside the anchors they gave
    saved <- project_state()
    expect_identical(saved$resources$cj_comparisons, attr(d, "comparisons"))
    expect_identical(saved$resources$cj_rankings, attr(d, "rankings"))
    expect_identical(saved$resources$anchors, e$.app_joint_anchors(cj))
    expect_false(saved$settings$cj_fix_comp)
    expect_false(saved$settings$cj_fix_rank)
    expect_no_error(.validate_app_project(saved))
    expect_no_error(.save_app_project(saved, path))

    # a structure the joint calibration cannot anchor stops the run and
    # leaves the analysis in place
    session$setInputs(thr_structure = "rsm"); session$flushReact()
    notes$msgs <- character(0)
    session$setInputs(run = 2L); session$flushReact()
    expect_match(notes$msgs, "choose the partial credit structure", all = FALSE)
    expect_identical(cj_fit(), cj)
    expect_identical(fit_val()$items, f$items)
    session$setInputs(thr_structure = "pcm", thr_mode = "pc", pc_rank = "4")
    session$flushReact()
    notes$msgs <- character(0)
    session$setInputs(run = 3L); session$flushReact()
    expect_match(notes$msgs, "principal-components threshold estimation",
                 all = FALSE)
    session$setInputs(thr_mode = "free"); session$flushReact()

    # anchors for equating and judgements are alternatives
    session$setInputs(anchor_file = list(name = "anchors.csv",
                                         datapath = anchors))
    session$flushReact()
    notes$msgs <- character(0)
    session$setInputs(run = 4L); session$flushReact()
    expect_match(notes$msgs, "choose either anchors for equating or judgements",
                 all = FALSE)
    expect_identical(cj_fit(), cj)
    session$setInputs(anchor_file = NULL); session$flushReact()

    # a fixed unit reaches rasch_cj() and the code
    session$setInputs(cj_fix_comp = TRUE); session$flushReact()
    notes$msgs <- character(0)
    session$setInputs(run = 5L); session$flushReact()
    expect_false(any(grepl("^\\[error\\]", notes$msgs)), info = notes$msgs)
    cj2 <- cj_fit()
    expect_false(cj2$units$estimated[cj2$units$frame == "comparisons"])
    expect_true(cj2$units$estimated[cj2$units$frame == "rankings"])
    expect_match(current_rcode(), "units = c(comparisons = 1, rankings = NA)",
                 fixed = TRUE)
    expect_match(output$joint_boxes$html, "300 comparisons, unit fixed",
                 fixed = TRUE)
    expect_true(attr(fit_val(), "rasch_app_source")$settings$cj_fix_comp)

    # a rank analysis takes the joint calibration with it
    session$setInputs(demo_choice = "pl", model_type = "btl",
                      bt_layout = "ranks", pl_ranking = "ranking",
                      pl_object = "object", pl_rank = "rank",
                      pl_judge = "judge", pl_se = "sandwich",
                      pl_split = "first", pl_ties = "drop")
    session$flushReact()
    expect_null(cj_fit())
    expect_null(cj_sources())
  })

  # reopening the saved analysis restores the anchored fit and refits the
  # joint calibration from the judgements it carries
  shiny::testServer(e$server, {
    session$flushReact()
    session$setInputs(project_file = list(
      datapath = path, name = "joint.rasch", size = file.info(path)$size,
      type = "application/octet-stream"))
    session$flushReact()
    f <- fit_val()
    expect_s3_class(f, "rasch")
    cj <- cj_fit()
    expect_s3_class(cj, "rasch_cj")
    expect_equal(cj$items$location[1:9], f$items$location[1:9],
                 tolerance = 1e-6)
    expect_identical(restored_project_resources()$cj_comparisons,
                     attr(d, "comparisons"))
    s <- cj_sources()
    expect_identical(s$origin, "project")
    expect_identical(s$rankings, attr(d, "rankings"))
    expect_match(output$cj_sources_note$html, "from the reopened analysis",
                 fixed = TRUE)
    expect_no_error(output$joint_boxes)
    expect_no_error(output$joint_map)
    restored <- e$.restored_input_values(restored_project_settings())
    expect_false(restored$cj_fix_comp)
    expect_false(restored$cj_fix_rank)
    expect_no_error(.validate_app_project(project_state()))

    # estimating again reads the judgements from the reopened analysis
    session$setInputs(model_type = "rasch", rasch_calibration = "free",
                      thr_structure = "pcm", thr_mode = "free",
                      anchor_type = "individual", item_cols = .joint_items,
                      id_col = "person_id", factor_cols = character(0),
                      cj_fix_comp = FALSE, cj_fix_rank = FALSE,
                      maxit = 200, tol = 1e-8, ng_auto = TRUE)
    session$flushReact()
    session$setInputs(run = 1L); session$flushReact()
    expect_s3_class(cj_fit(), "rasch_cj")
    expect_equal(cj_fit()$loglik, cj$loglik)
    code <- current_rcode()
    expect_match(code, 'project <- readRDS("joint.rasch")', fixed = TRUE)
    expect_match(code, "cj_comparisons <- project$resources$cj_comparisons",
                 fixed = TRUE)
    expect_match(code, "cj_rankings <- project$resources$cj_rankings",
                 fixed = TRUE)
  })
})

test_that("uploaded judgements are checked, reach rasch_cj() and are reproduced by the code", {
  skip_on_cran()
  notes <- .joint_notes()
  e <- .joint_app(notes)
  d <- .joint_demo()
  data_csv <- tempfile("responses", fileext = ".csv")
  comp_csv <- tempfile("comparisons", fileext = ".csv")
  rank_csv <- tempfile("rankings", fileext = ".csv")
  bad_csv <- tempfile("bad", fileext = ".csv")
  on.exit(unlink(c(data_csv, comp_csv, rank_csv, bad_csv)), add = TRUE)
  write.csv(as.data.frame(unclass(d), check.names = FALSE,
                          stringsAsFactors = FALSE),
            data_csv, row.names = FALSE)
  write.csv(attr(d, "comparisons"), comp_csv, row.names = FALSE)
  write.csv(attr(d, "rankings"), rank_csv, row.names = FALSE)
  write.csv(data.frame(a = "I01", b = "I02"), bad_csv, row.names = FALSE)
  shiny::testServer(e$server, {
    session$flushReact()
    session$setInputs(file = list(name = "responses.csv", datapath = data_csv,
                                  size = file.info(data_csv)$size,
                                  type = "text/csv"))
    session$flushReact()
    session$setInputs(demo_choice = "none", model_type = "rasch",
                      rasch_calibration = "free", thr_structure = "pcm",
                      thr_mode = "free", anchor_type = "individual",
                      item_cols = .joint_items, id_col = "person_id",
                      factor_cols = character(0), cj_fix_comp = FALSE,
                      cj_fix_rank = FALSE, maxit = 200, tol = 1e-8,
                      ng_auto = TRUE)
    session$flushReact()
    expect_null(attr(raw_data(), "comparisons", exact = TRUE))
    expect_null(cj_sources())
    expect_null(cj_sources_quiet())
    expect_null(output$cj_sources_note)

    # a CSV without the named columns is refused before it is used
    session$setInputs(cj_comp_file = list(name = "bad.csv", datapath = bad_csv))
    session$flushReact()
    expect_error(cj_sources(), "needs columns object_a, object_b, winner")
    expect_match(output$cj_sources_note$html, "needs columns object_a",
                 fixed = TRUE)
    expect_null(cj_sources_quiet())
    notes$msgs <- character(0)
    session$setInputs(run = 1L); session$flushReact()
    expect_match(notes$msgs, "needs columns object_a, object_b, winner",
                 all = FALSE)
    expect_null(cj_fit()); expect_null(fit_val())

    # rankings alone are a joint calibration too
    session$setInputs(cj_comp_file = NULL,
                      cj_rank_file = list(name = "rankings.csv",
                                          datapath = rank_csv))
    session$flushReact()
    s <- cj_sources()
    expect_identical(s$origin, "upload")
    expect_null(s$comparisons)
    expect_equal(nrow(s$rankings), 240L)
    expect_null(cj_sources_quiet()$comparisons)
    expect_false(is.null(cj_sources_quiet()$rankings))
    expect_match(output$cj_sources_note$html, "uploaded: 60 rankings.",
                 fixed = TRUE)
    notes$msgs <- character(0)
    session$setInputs(run = 2L); session$flushReact()
    expect_false(any(grepl("^\\[error\\]", notes$msgs)), info = notes$msgs)
    cj <- cj_fit()
    expect_s3_class(cj, "rasch_cj")
    expect_identical(cj$units$frame, c("responses", "rankings"))
    expect_equal(unname(cj$n[["comparisons"]]), 0L)
    expect_match(output$joint_boxes$html, "60 rankings", fixed = TRUE)
    expect_false(grepl("Unit of comparisons", output$joint_boxes$html,
                       fixed = TRUE))
    code <- current_rcode()
    expect_match(code, 'cj_rankings <- read.csv("rankings.csv"', fixed = TRUE)
    expect_false(grepl("cj_comparisons", code, fixed = TRUE))
    expect_match(code, "rankings = cj_rankings", fixed = TRUE)

    # both sources: the code reads each upload by name and reproduces the fit
    session$setInputs(cj_comp_file = list(name = "comparisons.csv",
                                          datapath = comp_csv))
    session$flushReact()
    notes$msgs <- character(0)
    session$setInputs(run = 3L); session$flushReact()
    expect_false(any(grepl("^\\[error\\]", notes$msgs)), info = notes$msgs)
    cj <- cj_fit()
    expect_identical(unname(cj$n), c(232L, 300L, 60L))
    code <- current_rcode()
    expect_match(code, 'cj_comparisons <- read.csv("comparisons.csv"', fixed = TRUE)
    code <- sub('read.csv("responses.csv"', paste0('read.csv("', data_csv, '"'),
                code, fixed = TRUE)
    code <- sub('read.csv("comparisons.csv"', paste0('read.csv("', comp_csv, '"'),
                code, fixed = TRUE)
    code <- sub('read.csv("rankings.csv"', paste0('read.csv("', rank_csv, '"'),
                code, fixed = TRUE)
    env <- new.env(parent = globalenv())
    eval(parse(text = code), envir = env)
    expect_equal(env$cj$loglik, cj$loglik)
    expect_equal(env$fit$items$location, fit_val()$items$location)
    saved <- project_state()
    expect_equal(nrow(saved$resources$cj_comparisons), 300L)
    expect_no_error(.validate_app_project(saved))
  })
})
