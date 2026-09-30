# Regressions for the score-to-measure table of a split fit: each person
# answers one copy of a split item, so the conversion is given per test form
# rather than over a raw score no person can reach.

# 400 persons in two groups, item I3 harder for group b, split by group
sim_split_fit <- function(seed = 3, n = 400L) {
  set.seed(seed)
  d <- seq(-1.5, 1.5, length.out = 6)
  g <- rep(c("a", "b"), each = n / 2)
  sh <- matrix(0, n, 6)
  sh[g == "b", 3] <- 1.2
  X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, grp = g), factors = "grp")
  list(fit = fit, split = split_items(fit, "I3", by = "grp"), X = X, g = g)
}

test_that("a split fit converts each test form's own raw score", {
  s <- sim_split_fit()
  st <- score_table(s$split)
  expect_identical(names(st), c("form", "score", "theta", "se",
                                "extrapolated", "freq", "cum_pct"))
  expect_identical(sort(unique(st$form)), c("I3 (a)", "I3 (b)"))
  # six items answered, six the highest score: never the seven of the
  # calibrated columns
  expect_identical(st$score[st$form == "I3 (a)"], 0:6)
  expect_identical(st$score[st$form == "I3 (b)"], 0:6)
  expect_equal(max(st$score), 6L)
  expect_equal(nrow(st), 14L)
  expect_true(all(is.finite(st$theta)))
  expect_false(any(st$extrapolated))

  # each form is the Warm conversion over the unsplit items and its own copy
  items <- s$split$items$item
  for (lv in c("a", "b")) {
    own <- paste0("I3 (", lv, ")")
    idx <- match(c(setdiff(items, c("I3 (a)", "I3 (b)")), own), items)
    pe <- person_wle(s$split$tau_list[idx])
    rows <- st$form == own
    expect_equal(st$theta[rows], unname(pe$theta))
    expect_equal(st$se[rows], unname(pe$se))
    # its frequencies are the group's own raw-score distribution
    raw <- rowSums(s$X[s$g == lv, ])
    expect_identical(st$freq[rows],
                     as.integer(table(factor(raw, levels = 0:6))))
    expect_equal(st$cum_pct[rows][7L], 100)
  }
  expect_equal(sum(st$freq), nrow(s$X))
  # the two copies differ, so the conversions differ where the copies bite
  expect_false(isTRUE(all.equal(st$theta[st$form == "I3 (a)"],
                                st$theta[st$form == "I3 (b)"])))
})

test_that("the estimator and extreme-score options apply per form", {
  s <- sim_split_fit()
  mle <- score_table(s$split, method = "mle")
  expect_true(all(is.na(mle$theta[mle$score %in% c(0L, 6L)])))
  expect_true(all(is.finite(mle$theta[!mle$score %in% c(0L, 6L)])))
  ext <- score_table(s$split, method = "mle", extremes = "extrapolated")
  expect_identical(ext$extrapolated, ext$score %in% c(0L, 6L))
  expect_true(all(is.finite(ext$theta)))
  for (lv in c("I3 (a)", "I3 (b)")) {
    th <- ext$theta[ext$form == lv]
    expect_true(all(diff(th) > 0))
    # the geometric continuation: the outer difference is b^2 / a
    dd <- diff(th)
    expect_equal(dd[6L], dd[5L]^2 / dd[4L])
    expect_equal(dd[1L], dd[2L]^2 / dd[3L])
  }
})

test_that("persons who answered no copy of a split item belong to no form", {
  s <- sim_split_fit()
  X2 <- s$X
  X2[1:5, "I3"] <- NA          # group a persons without the split item
  X2[6:8, "I1"] <- NA          # group a persons with an unsplit item missing
  fit <- split_items(rasch(data.frame(X2, grp = s$g), factors = "grp"),
                     "I3", by = "grp")
  forms <- rasch:::.split_forms(fit)
  expect_length(forms, 2L)
  expect_identical(vapply(forms, `[[`, "", "label"), c("I3 (a)", "I3 (b)"))
  # the five without any copy are in no form; the three with a missing
  # unsplit item are in group a's form but are not complete responders
  expect_identical(sort(unlist(lapply(forms, `[[`, "rows"))),
                    setdiff(seq_len(nrow(X2)), 1:5))
  st <- score_table(fit)
  expect_equal(sum(st$freq[st$form == "I3 (a)"]), 200L - 8L)
  expect_equal(sum(st$freq[st$form == "I3 (b)"]), 200L)
  expect_equal(max(st$score), 6L)

  # a form nobody completed keeps its conversion, with empty counts
  none <- fit
  none$X[none$X[, "I3 (b)"] %in% 0:1 & !is.na(none$X[, "I3 (b)"]), "I2"] <- NA
  st2 <- score_table(none)
  expect_identical(st2$freq[st2$form == "I3 (b)"], rep(0L, 7L))
  expect_true(all(is.na(st2$cum_pct[st2$form == "I3 (b)"])))
  expect_equal(st2$theta[st2$form == "I3 (b)"], st$theta[st$form == "I3 (b)"])

  # and when no person answered a whole form, there is no table
  gone <- fit
  gone$X[, c("I3 (a)", "I3 (b)")] <- NA
  expect_identical(rasch:::.split_forms(gone), list())
  expect_null(score_table(gone))
})

test_that("an unsplit fit's table is unchanged and a resolved fit has forms", {
  s <- sim_split_fit()
  st <- score_table(s$fit)
  expect_identical(names(st), c("score", "theta", "se", "extrapolated",
                                "freq", "cum_pct"))
  expect_identical(st$score, 0:6)
  expect_equal(st$theta, unname(person_wle(s$fit$tau_list)$theta))
  expect_null(rasch:::.split_forms(s$fit))
  skip_on_cran()
  r <- resolve_dif(s$fit)
  expect_true(any(grepl("^I3 \\(", r$fit$items$item)))
  rt <- score_table(r$fit)
  expect_true("form" %in% names(rt))
  expect_equal(max(rt$score), 6L)
  expect_equal(sum(rt$freq), nrow(s$X))
})

test_that("the reports and the CSV export carry the per-form conversion", {
  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available("2.8"))
  s <- sim_split_fit()
  md <- tempfile(fileext = ".md")
  on.exit(unlink(md), add = TRUE)
  report_document(s$split, md, format = "md")
  lines <- readLines(md, warn = FALSE)
  at <- grep("^#+ Score to measure", lines)
  expect_length(at, 1L)
  section <- lines[at:length(lines)]
  nxt <- grep("^#+ ", section)[-1L]
  if (length(nxt)) section <- section[seq_len(nxt[1L] - 1L)]
  expect_true(any(grepl("no raw score is common to every person", section,
                        fixed = TRUE)))
  expect_true(any(grepl("I3 (a)", section, fixed = TRUE)))
  expect_true(any(grepl("I3 (b)", section, fixed = TRUE)))
  # rows are the per-form scores: none reaches the seven calibrated columns
  scores <- sub("^\\|[^|]*\\|[[:space:]]*([0-9]+)[[:space:]]*\\|.*$", "\\1",
                grep("^\\|[[:space:]]*I3 \\(", section, value = TRUE))
  expect_true(length(scores) >= 14L)
  expect_true(all(as.integer(scores) <= 6L))

  out <- tempfile("rasch-forms-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  # the split copies share no respondents, so the residual PCA plot is
  # skipped with a warning as before; the tables are still written
  expect_warning(save_outputs(s$split, out, formats = "png",
                              item_plots = FALSE), "pca_loadings")
  csv <- file.path(out, "tables", "score_to_measure.csv")
  expect_true(file.exists(csv))
  got <- read.csv(csv, stringsAsFactors = FALSE)
  expect_identical(names(got)[1:2], c("form", "score"))
  expect_equal(max(got$score), 6L)

  # a split fit nobody completed writes no conversion and says why
  gone <- s$split
  gone$X[, c("I3 (a)", "I3 (b)")] <- NA
  md2 <- tempfile(fileext = ".md")
  out2 <- tempfile("rasch-forms-")
  on.exit(unlink(c(md2, out2), recursive = TRUE), add = TRUE)
  report_document(gone, md2, format = "md")
  lines2 <- readLines(md2, warn = FALSE)
  expect_true(any(grepl("Not available: no person answered a whole form",
                        lines2, fixed = TRUE)))
  expect_warning(save_outputs(gone, out2, formats = "png",
                              item_plots = FALSE), "pca_loadings")
  expect_false(file.exists(file.path(out2, "tables", "score_to_measure.csv")))
})
