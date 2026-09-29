.audit_app <- function() {
  for (pkg in c("shiny", "bslib", "DT", "bsicons")) skip_if_not_installed(pkg)
  path <- test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(path)) path <- system.file("shiny", "app.R", package = "rasch")
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(path, e))
  e
}

.audit_tile <- function(html, label) {
  tiles <- strsplit(html, '<div class="metric-tile ', fixed = TRUE)[[1L]]
  tiles[vapply(tiles, function(x)
    grepl(paste0('<span class="metric-label">', label, '</span>'), x, fixed = TRUE),
    TRUE)]
}

test_that("invariance summaries separate withheld and nonsignificant tests", {
  x <- .invariance_summary(NULL, NULL)
  expect_true(is.na(x$n_moved))
  expect_identical(x$status, "neutral")
  expect_match(x$text, "unavailable")
  tab <- data.frame(p_adj = c(NA_real_, 0.5, 0.7))
  x <- .invariance_summary(tab, c("a", "b", "c"))
  expect_equal(x$tested, 2L)
  expect_equal(x$withheld, 1L)
  expect_equal(x$n_moved, 0L)
  expect_identical(x$status, "neutral")
  expect_match(x$text, "2 contrasts tested; 1 withheld", fixed = TRUE)
  tab$p_adj[1] <- .6
  expect_identical(.invariance_summary(tab, letters[1:3])$status, "good")
  tab$p_adj[2:3] <- .01
  x <- .invariance_summary(tab, c("a", "b", "b"))
  expect_identical(x$status, "bad")
  expect_equal(x$n_moved, 1L) # one object, two tested contrasts
  # A failed frame can contribute no rows, rather than rows with NA p.
  x <- .invariance_summary(data.frame(p_adj = .8), "a", n_expected = 3L)
  expect_equal(x$withheld, 2L)
  expect_identical(x$status, "neutral")
  expect_match(x$text, "2 withheld", fixed = TRUE)
})

test_that("duplicate bundles are refused, and generated bundle code round-trips", {
  expect_error(.parse_dif_bundles("domain: I1,I2\n domain : I3,I4"),
               "line 2 repeats the name 'domain'")
  skip_on_cran()
  e <- .audit_app()
  shiny::testServer(e$server, {
    b <- setNames(rep(list(c("I1", "I2")), 6),
                  c("domain\\q", "a`b", 'a"b', "a\nb", "if", "simple"))
    expect_identical(eval(parse(text = bundles_code(b))), b)
    expect_identical(bundles_code(list(simple = c("I1", "I2"))),
                     'list(simple = c("I1", "I2"))')
  })
})

test_that("rankings tiles retain inference availability and reversal direction", {
  skip_on_cran()
  set.seed(915)
  d <- do.call(rbind, lapply(1:200, function(i)
    data.frame(ranking = i, object = sample(LETTERS[1:4], prob = exp(c(-1,-.4,.2,1.2))),
               rank = 1:4, judge = paste0("J", (i-1) %% 5 + 1))))
  k <- pl(d, judge = "judge")
  expect_false(k$se_available)
  expect_true(all(is.na(k$invariance$objects$p_adj)))
  expect_output(print(k), "objects moving.*unavailable")
  e <- .audit_app()
  shiny::testServer(e$server, {
    session$flushReact()
    pl_fit(k); session$flushReact()
    tile <- .audit_tile(output$pl_boxes$html, "Invariance")
    expect_length(tile, 1L)
    expect_match(tile, '^metric-neutral"')
    expect_match(tile, 'class="metric-value">Unavailable', fixed = TRUE)
    expect_match(output$pl_invariance_note$html, "unavailable")
    # Hold the test result fixed to exercise the display's direction rule.
    a <- k; a$reversal$p <- .001
    a$reversal$loglik_forward <- -100; a$reversal$loglik_reversed <- -120
    pl_fit(a); session$flushReact()
    expect_match(.audit_tile(output$pl_boxes$html, "Reversal check"), '^metric-good"')
    a$reversal$loglik_reversed <- -80
    pl_fit(a); session$flushReact()
    expect_match(.audit_tile(output$pl_boxes$html, "Reversal check"), '^metric-bad"')
    a$reversal$p <- NA_real_
    pl_fit(a); session$flushReact()
    expect_match(.audit_tile(output$pl_boxes$html, "Reversal check"), '^metric-neutral"')
    a$invariance$objects$p_adj[1] <- .8
    pl_fit(a); session$flushReact()
    expect_match(.audit_tile(output$pl_boxes$html, "Invariance"), '^metric-neutral"')
    expect_match(output$pl_invariance_note$html, "1 contrasts tested; 3 withheld", fixed = TRUE)
  })
})

test_that("joint calibration summaries do not count untested objects as passing", {
  skip_on_cran()
  d <- simulate_rasch(180, 4, seed = 521)
  k <- rasch_cj(d, items = sprintf("I%02d", 1:4),
     comparisons = data.frame(object_a = rep("I01", 20),
                              object_b = "I02", winner = "I02"),
     units = c(comparisons = 1))
  expect_true(k$converged)
  expect_null(k$invariance$items)
  expect_equal(k$invariance$n_contrasts, 2L)
  expect_output(print(k), "Objects differing.*unavailable")
  expect_true(all(is.na(subset(k$frame_locations, frame == "comparisons")$location)))
  e <- .audit_app()
  shiny::testServer(e$server, {
    session$flushReact()
    cj_fit(k); session$flushReact()
    tile <- .audit_tile(output$joint_boxes$html, "Objects moving")
    expect_length(tile, 1L)
    expect_match(tile, '^metric-neutral"')
    expect_match(tile, 'class="metric-value">Unavailable', fixed = TRUE)
    expect_match(output$joint_fitsum_tbl$html, "unavailable")
    expect_match(output$joint_invariance_note$html, "unavailable")
    tab <- joint_summary_table(k)
    expect_true(is.na(tab$value[tab$statistic == "Objects moving (Holm p < .05)"]))
    expect_identical(tab$value[tab$statistic == "Object contrasts tested"], "0")
    expect_identical(tab$value[tab$statistic == "Object contrasts withheld"], "2")
    a <- k
    a$invariance$items <- data.frame(item = c("I01", "I02"), p_adj = c(.7, NA_real_))
    cj_fit(a); session$flushReact()
    expect_match(.audit_tile(output$joint_boxes$html, "Objects moving"), '^metric-neutral"')
    expect_match(output$joint_invariance_note$html, "1 contrasts tested; 1 withheld", fixed = TRUE)
  })
  # With two judgement sources the failed comparisons must still count as
  # withheld even though the valid rankings contribute a finite table.
  set.seed(44)
  ranks <- do.call(rbind, lapply(1:180, function(i)
    data.frame(ranking = i, item = sample(sprintf("I%02d", 1:4),
      prob = exp(c(-2.4, -.7, .8, 2.3))), rank = 1:4)))
  mixed <- rasch_cj(d, items = sprintf("I%02d", 1:4),
    comparisons = data.frame(object_a = rep("I01", 20), object_b = "I02", winner = "I02"),
    rankings = ranks, units = c(comparisons = 1, rankings = 1))
  expect_true(mixed$converged)
  expect_true(all(is.finite(mixed$invariance$items$p_adj)))
  expect_equal(nrow(mixed$invariance$items), 4L)
  expect_equal(mixed$invariance$n_contrasts, 6L)
  expect_output(print(mixed), "4 contrasts tested; 2 withheld")
})

test_that("DTF help describes test-level cancellation and pointwise intervals", {
  path <- test_path("..", "..", "inst", "shiny", "help.R")
  if (!file.exists(path)) path <- system.file("shiny", "help.R", package = "rasch")
  h <- new.env(parent = baseenv()); sys.source(path, h)
  expect_match(h$APP_HELP[["dtf_test_tbl"]], "Opposing item effects can cancel", fixed = TRUE)
  expect_match(h$APP_HELP[["dtf_plot"]], "pointwise 95 per cent", fixed = TRUE)
  expect_match(h$APP_HELP[["dtf_plot"]], "Excluding zero", fixed = TRUE)
  a <- list(tau = list(-1, 1), ids = list(1L, 2L), n_par = 4)
  b <- list(tau = list(1, -1), ids = list(3L, 4L), n_par = 4)
  s <- .dtf_summaries(seq(-2, 2, length.out = 7), a, b, c(-10, 10), diag(.01, 4), 2, NULL)
  expect_equal(unname(unlist(s[c("sDBF_score", "uDBF_score")])), c(0, 0),
               tolerance = 1e-12)
})
