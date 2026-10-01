# Regressions for the designs of a split fit: the copies of a split item
# are answered by different persons, so the information, the characteristic
# curve and the score-to-measure conversion describe each test form, never
# the calibrated columns together; and a copy left alone by drop_items()
# is still a copy, whose source the other persons answered no copy of.

# 400 persons in two groups, I3 harder and I5 easier for group b
sim_two_dif <- function(seed = 3, n = 400L) {
  set.seed(seed)
  d <- seq(-1.5, 1.5, length.out = 6)
  g <- rep(c("a", "b"), each = n / 2)
  sh <- matrix(0, n, 6)
  sh[g == "b", 3] <- 1.2
  sh[g == "b", 5] <- -0.9
  X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
  colnames(X) <- paste0("I", 1:6)
  list(fit = rasch(data.frame(X, grp = g), factors = "grp"), X = X, g = g)
}

form_info <- function(fit, items, theta = 0) {
  idx <- match(items, fit$items$item)
  sum(vapply(idx, function(i) item_moments(theta, fit$tau_list[[i]])$V, 0))
}

test_that("a split fit's information and curves are per test form", {
  s <- sim_two_dif()
  sp <- split_items(s$fit, c("I3", "I5"), by = "grp")
  blocks <- rasch:::.design_blocks(sp)
  expect_identical(names(blocks), c("I3 (a), I5 (a)", "I3 (b), I5 (b)"))
  unsplit <- c("I1", "I2", "I4", "I6")
  expect_identical(sp$items$item[blocks[[1L]]],
                   c(unsplit, "I3 (a)", "I5 (a)"))
  expect_identical(sp$items$item[blocks[[2L]]],
                   c(unsplit, "I3 (b)", "I5 (b)"))
  # six items per form, never the eight calibrated columns
  expect_identical(unname(vapply(blocks, function(ii) sum(sp$m[ii]), 0)),
                   c(6, 6))

  ti <- test_information(sp, c(-1, 0, 1))
  expect_identical(names(ti), c("theta", "info", "sem", "design"))
  expect_identical(unique(ti$design), names(blocks))
  for (j in 1:2) {
    z <- ti$design == names(blocks)[j]
    expect_equal(ti$info[z], vapply(c(-1, 0, 1), function(th)
      form_info(sp, sp$items$item[blocks[[j]]], th), 0))
  }
  # the eight-column sum claimed the information of both copies at once
  eight <- form_info(sp, sp$items$item)
  expect_true(all(ti$sem[ti$theta == 0] > 1 / sqrt(eight)))
  expect_true(all(ti$sem[ti$theta == 0] > 0.85))
  expect_true(1 / sqrt(eight) < 0.8)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot_tcc(sp))
  expect_no_error(plot_tif(sp))
  expect_no_error(plot_pimap(sp, information = TRUE))
  expect_no_error(plot_pimap(sp, information = TRUE, group = "b"))
  expect_no_error(plot_pimap(sp, information = TRUE,
                             items = c("I1", "I2", "I3 (a)")))
})

test_that("a person group's map carries only the forms its persons answered", {
  s <- sim_two_dif()
  sp <- split_items(s$fit, c("I3", "I5"), by = "grp")
  one <- sp
  one$X <- sp$X[s$g == "b", , drop = FALSE]
  ti <- test_information(one, 0)
  expect_identical(ti$design, "I3 (b), I5 (b)")
  expect_equal(ti$info, form_info(sp, c("I1", "I2", "I4", "I6",
                                        "I3 (b)", "I5 (b)")))
})

test_that("an item selection names its designs within the selection", {
  s <- sim_two_dif()
  sp <- split_items(s$fit, c("I3", "I5"), by = "grp")
  # no copy selected: the forms agree and the curve is the test's
  ti <- test_information(sp, 0, items = c("I1", "I2"))
  expect_identical(names(ti), c("theta", "info", "sem"))
  expect_equal(ti$info, form_info(sp, c("I1", "I2")))
  # one copy selected: the form with it and the form without it
  ti <- test_information(sp, 0, items = c("I1", "I2", "I3 (a)"))
  expect_identical(ti$design, c("I3 (a)", "without I3"))
  expect_equal(ti$info, c(form_info(sp, c("I1", "I2", "I3 (a)")),
                          form_info(sp, c("I1", "I2"))))
  ti <- test_information(sp, 0, items = c("I3 (a)", "I3 (b)"))
  expect_identical(ti$design, c("I3 (a)", "I3 (b)"))
  expect_equal(ti$info, c(form_info(sp, "I3 (a)"), form_info(sp, "I3 (b)")))
})

test_that("a split copy left alone by drop_items keeps its source split", {
  s <- sim_two_dif()
  sp <- split_items(s$fit, "I1", by = "grp")
  dr <- drop_items(sp, "I1 (b)")
  expect_identical(dr$split_map[["I1 (a)"]], "I1")
  forms <- rasch:::.split_forms(dr)
  expect_length(forms, 2L)
  expect_identical(vapply(forms, `[[`, "", "label"), c("I1 (a)", "without I1"))
  five <- paste0("I", 2:6)
  expect_identical(forms[[1L]]$items, c(five, "I1 (a)"))
  expect_identical(forms[[2L]]$items, five)
  expect_identical(forms[[1L]]$rows, which(s$g == "a"))
  expect_identical(forms[[2L]]$rows, which(s$g == "b"))

  # group b's conversion is over the five items it answered
  st <- score_table(dr)
  expect_identical(names(st)[1:2], c("form", "score"))
  expect_identical(st$score[st$form == "I1 (a)"], 0:6)
  expect_identical(st$score[st$form == "without I1"], 0:5)
  idx <- match(five, dr$items$item)
  pe <- person_wle(dr$tau_list[idx])
  rows <- st$form == "without I1"
  expect_equal(st$theta[rows], unname(pe$theta))
  expect_equal(st$se[rows], unname(pe$se))
  raw_b <- rowSums(s$X[s$g == "b", five])
  expect_identical(st$freq[rows],
                   as.integer(table(factor(raw_b, levels = 0:5))))
  expect_equal(sum(st$freq), nrow(s$X))
  # the six-item conversion the group was given differs by over a logit
  six <- person_wle(dr$tau_list)$theta
  expect_true(max(abs(six[1:6] - pe$theta)) > 1)

  # the information and curves follow the same two forms
  ti <- test_information(dr, 0)
  expect_identical(ti$design, c("I1 (a)", "without I1"))
  expect_equal(ti$info, c(form_info(dr, c(five, "I1 (a)")),
                          form_info(dr, five)))
  blocks <- rasch:::.design_blocks(dr)
  expect_identical(unname(vapply(blocks, function(ii) sum(dr$m[ii]), 0)),
                   c(6, 5))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot_tcc(dr))
  expect_no_error(plot_tif(dr))
  expect_no_error(plot_pimap(dr, information = TRUE))
})

test_that("a level the item was not split for takes the form without it", {
  s <- sim_two_dif()
  g3 <- rep(c("a", "b", "c"), length.out = nrow(s$X))
  X3 <- s$X
  X3[g3 == "c", "I3"] <- 0L     # level c never scores 1 on I3
  fit <- rasch(data.frame(X3, grp = g3), factors = "grp")
  sp <- split_items(fit, "I3", by = "grp")
  expect_true(any(grepl("not split for unavailable level(s): c", sp$notes,
                        fixed = TRUE)))
  forms <- rasch:::.split_forms(sp)
  expect_identical(vapply(forms, `[[`, "", "label"),
                   c("I3 (a)", "I3 (b)", "without I3"))
  expect_identical(forms[[3L]]$rows, which(g3 == "c"))
  expect_identical(forms[[3L]]$items, paste0("I", c(1, 2, 4, 5, 6)))
  st <- score_table(sp)
  expect_equal(sum(st$freq[st$form == "without I3"]), sum(g3 == "c"))
  expect_identical(st$score[st$form == "without I3"], 0:5)
  ti <- test_information(sp, 0)
  expect_identical(ti$design, c("I3 (a)", "I3 (b)", "without I3"))
})

test_that("missing responses to copies give the sparser forms", {
  s <- sim_two_dif()
  X <- s$X
  X[1:3, "I3"] <- NA
  X[4, c("I3", "I5")] <- NA
  X[201:202, "I5"] <- NA
  X[5, "I2"] <- NA              # an unsplit item: still the full form
  fit <- rasch(data.frame(X, grp = s$g), factors = "grp")
  sp <- split_items(fit, c("I3", "I5"), by = "grp")
  forms <- rasch:::.split_forms(sp)
  labels <- vapply(forms, `[[`, "", "label")
  expect_identical(labels, c("I3 (a), I5 (a)", "I3 (b), I5 (b)",
                             "I3 (b), without I5", "I5 (a), without I3",
                             "without I3 and I5"))
  expect_identical(forms[[4L]]$rows, 1:3)
  expect_identical(forms[[5L]]$rows, 4L)
  expect_identical(forms[[3L]]$rows, 201:202)
  expect_true(5L %in% forms[[1L]]$rows)
  expect_identical(forms[[5L]]$items, c("I1", "I2", "I4", "I6"))
  st <- score_table(sp)
  expect_equal(sum(st$freq[st$form == "I3 (a), I5 (a)"]), 200L - 5L)
  expect_equal(sum(st$freq[st$form == "without I3 and I5"]), 1L)
  expect_identical(st$score[st$form == "without I3 and I5"], 0:4)
  # a person who answered several copies of an item took no form
  both <- sp
  both$X[7, "I3 (b)"] <- 1L
  forms2 <- rasch:::.split_forms(both)
  expect_false(7L %in% unlist(lapply(forms2, `[[`, "rows")))
})

test_that("form labels name the copies held and the sources without one", {
  map <- c(I1 = "I1", "I3 (a)" = "I3", "I3 (b)" = "I3", "I5 (a)" = "I5",
           "I2 (b)" = "I2")
  lab <- rasch:::.split_form_label
  expect_identical(lab(map, c("I1", "I3 (a)", "I5 (a)", "I2 (b)")),
                   "I3 (a), I5 (a), I2 (b)")
  expect_identical(lab(map, c("I1", "I3 (b)")), "I3 (b), without I5 and I2")
  expect_identical(lab(map, "I1"), "without I3, I5 and I2")
  expect_identical(rasch:::.split_copies(map),
                   c("I3 (a)", "I3 (b)", "I5 (a)", "I2 (b)"))
})

test_that("an unsplit fit keeps one whole-test design", {
  s <- sim_two_dif()
  expect_identical(rasch:::.design_blocks(s$fit), list(test = 1:6))
  expect_identical(rasch:::.design_label(s$fit, 1:3), "test")
  ti <- test_information(s$fit, 0)
  expect_identical(names(ti), c("theta", "info", "sem"))
  expect_equal(ti$info, form_info(s$fit, s$fit$items$item))
  expect_null(rasch:::.split_forms(s$fit))
})

test_that("the reports carry the form a dropped copy leaves", {
  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available("2.8"))
  s <- sim_two_dif()
  dr <- drop_items(split_items(s$fit, "I1", by = "grp"), "I1 (b)")
  md <- tempfile(fileext = ".md")
  on.exit(unlink(md), add = TRUE)
  report_document(dr, md, format = "md")
  lines <- readLines(md, warn = FALSE)
  expect_true(any(grepl("^\\|[[:space:]]*without I1", lines)))
  expect_true(any(grepl("split items they answered no copy of", lines,
                        fixed = TRUE)))
  html <- tempfile(fileext = ".html")
  on.exit(unlink(html), add = TRUE)
  report_html(dr, html)
  page <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_true(grepl("without I1", page, fixed = TRUE))
  expect_true(grepl("split items they answered no copy of", page,
                    fixed = TRUE))
  # and the HTML report says why a split fit with no form has no table
  gone <- split_items(s$fit, "I1", by = "grp")
  gone$X[, "I1 (a)"] <- 1L
  gone$X[, "I1 (b)"] <- 0L
  expect_null(score_table(gone))
  html2 <- tempfile(fileext = ".html")
  on.exit(unlink(html2), add = TRUE)
  # its residuals are still the split fit's, so the scree plot is skipped
  expect_warning(report_html(gone, html2), "scree")
  page2 <- paste(readLines(html2, warn = FALSE), collapse = "\n")
  expect_true(grepl("Not available: no person answered a form", page2,
                    fixed = TRUE))
})
