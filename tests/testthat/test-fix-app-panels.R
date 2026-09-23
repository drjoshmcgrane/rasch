.panels_app_path <- function(file = "app.R") {
  path <- testthat::test_path("..", "..", "inst", "shiny", file)
  if (!file.exists(path)) path <- system.file("shiny", file, package = "rasch")
  path
}

.panels_app_env <- function() {
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.panels_app_path(), envir = e))
  e
}

.panels_csv <- function(txt) {
  p <- tempfile(fileext = ".csv")
  writeLines(txt, p)
  list(datapath = p, name = "map.csv", size = nchar(txt), type = "text/csv")
}

.panels_strip <- function(h) gsub("\\s+", " ", gsub("<[^>]+>", " ", as.character(h)))

test_that("a header-only set map is refused instead of becoming one '(rest)' set", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- .panels_app_env()
  d <- as.data.frame(simulate_rasch(n_persons = 120, n_items = 6, seed = 81))
  its <- sprintf("I%02d", 1:6)

  shiny::testServer(e$server, {
    sim_data(d)
    session$setInputs(model_type = "efrm", ef_items = its,
                      ef_sets = .panels_csv("item,set"))
    session$flushReact()
    # the missing-id fill-in below read a header-only map as "every item is
    # in (rest)", which is a design choice the analyst never made
    err <- tryCatch(ef_setmap(), error = function(e) e)
    expect_s3_class(err, "shiny.silent.error")
    expect_match(conditionMessage(err), "item-set CSV has no rows", fixed = TRUE)

    # the same reader serves the paired-comparison frames panel
    obj <- tryCatch(read_frame_map(.panels_csv("object,set"), "object",
                                   c("A", "B"), "object"),
                    error = function(e) e)
    expect_s3_class(obj, "shiny.silent.error")
    expect_match(conditionMessage(obj), "object-set CSV has no rows", fixed = TRUE)

    # a map with rows still assigns, and still fills the items it omits
    session$setInputs(ef_sets = .panels_csv(c("item,set", "I01,A", "I02,A")))
    session$flushReact()
    expect_equal(unname(ef_setmap()[c("I01", "I03")]), c("A", "(rest)"))
  })
})

test_that("the explanatory tile says no comparison exists at 0 df", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- .panels_app_env()
  # a four-level item factor on four items spans every free item parameter,
  # so the nested comparison has no degrees of freedom left
  set.seed(1)
  q <- data.frame(item = paste0("I", 1:4), grp = factor(c("a", "b", "c", "d")))
  X <- matrix(rbinom(300 * 4, 1, stats::plogis(outer(rnorm(300),
                                                     c(-1, -.3, .3, 1), "-"))),
              300, 4)
  colnames(X) <- q$item
  f <- rasch_explanatory(X, predictors = q, formula = ~ grp)
  tst <- explanatory_test(f)
  expect_equal(tst$df[1L], 0)
  expect_true(is.na(tst$p_kent[1L]))

  shiny::testServer(e$server, {
    session$flushReact()
    fit_val(f); session$flushReact()
    box <- .panels_strip(output$expl_boxes$html)
    expect_false(grepl("p = NA", box, fixed = TRUE))
    expect_match(box, "No comparison", fixed = TRUE)
    expect_match(box, "span every free item parameter (0 df)", fixed = TRUE)
  })
})

test_that("several DIF resolution notes are listed, not run into one sentence", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- .panels_app_env()
  set.seed(11)
  d <- simulate_rasch(n_persons = 400, n_items = 8, model = "PCM",
                      n_categories = 4,
                      dif = list(items = c("I02", "I05"), uniform = 1.5),
                      n_groups = 2, seed = 12)
  X <- as.data.frame(d)
  items <- grep("^I", names(X), value = TRUE)
  grp <- as.character(X$group)
  # a thin second group refuses one split, and an anchored item refuses another
  keep <- c(which(grp == "g1"), head(which(grp == "g2"), 12))
  gf <- data.frame(grp = grp[keep])
  f0 <- rasch(X[keep, items], factors = gf)
  anc <- data.frame(item = "I05", k = 1L,
                    tau = f0$items$location[f0$items$item == "I05"])
  f <- rasch(X[keep, items], factors = gf, anchors = anc)

  shiny::testServer(e$server, {
    session$flushReact()
    fit_val(f); session$flushReact()
    session$setInputs(resolve_all = 1L); session$flushReact()
    rr <- resolve_res()
    expect_true(length(rr$notes) > 1L)
    html <- output$resolve_notes$html
    expect_match(html, "<li>", fixed = TRUE)
    # the boundary between two refusals is visible
    expect_false(grepl(paste(rr$notes, collapse = " "), html, fixed = TRUE))
    for (n in rr$notes)
      expect_match(html, n, fixed = TRUE)
  })
})

test_that("the DIF footer carries every dif_anova note", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- .panels_app_env()
  set.seed(3)
  X <- as.data.frame(simulate_rasch(n_persons = 600, n_items = 10, seed = 3))
  X <- X[, sprintf("I%02d", 1:10)]
  id <- rep(1:300, 2)
  occ <- rep(c("t1", "t2"), each = 300)
  sex <- rep(rep(c("m", "f"), 150), 2)
  # 40 persons miss I03 at the second occasion: those panels leave the
  # within-person tests and the remaining between-person F is approximate
  X[occ == "t2" & id <= 40, "I03"] <- NA
  f <- rasch(X, factors = data.frame(occ = occ, sex = sex,
                                     stringsAsFactors = FALSE), id = id)

  shiny::testServer(e$server, {
    fit_val(f); session$flushReact()
    notes <- dif_res()$notes
    expect_true(length(notes) > 1L)
    foot <- .panels_strip(output$dif_note$html)
    expect_match(foot, "available DIF effects significant after adjustment",
                 fixed = TRUE)
    for (n in notes)
      expect_match(foot, .panels_strip(n), fixed = TRUE)
    # the count line is still given once
    expect_equal(length(gregexpr("Class intervals", foot, fixed = TRUE)[[1]]), 1L)
  })
})
