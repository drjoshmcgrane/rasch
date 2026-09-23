# Round-three repairs: an empty equating reference, a corrupt saved
# frame-invariance result, the legacy person note for a layout that pairs by
# identifier only, and the shape of the B = 0 verdict note.

test_that("an equating reference with no items is refused by name", {
  empty <- data.frame(item = character(), location = numeric(),
                      se = numeric())
  expect_error(.equate_ref(empty), "reference has no items")
  # the same refusal whichever optional column the caller supplied: the
  # length-one defaults below the check never see a zero-row frame
  expect_error(.equate_ref(data.frame(item = character(),
                                      location = numeric())),
               "reference has no items")
  expect_error(.equate_ref(data.frame(item = character(),
                                      location = numeric(), se = numeric(),
                                      max = integer())),
               "reference has no items")
  f <- rasch(simulate_rasch(120, 6, seed = 81), id = "id")
  msg <- tryCatch(equate_tests(f, empty), error = conditionMessage)
  expect_match(msg, "reference has no items")
  expect_false(grepl("replacement has 1 row", msg, fixed = TRUE))
  expect_error(plot_equate(f, empty), "reference has no items")
})

test_that("a non-list frame-invariance result is refused in the package's words", {
  fit <- structure(list(), class = c("rasch_efrm", "rasch", "list"))
  msg <- tryCatch(.validate_frame_invariance("not a list", fit),
                  error = conditionMessage)
  expect_match(msg, "must be a frame_invariance\\(\\) result")
  expect_false(grepl("atomic", msg, fixed = TRUE))
  # a class attribute pasted onto an atomic value is still not a result
  msg2 <- tryCatch(.validate_frame_invariance(
    structure(1:3, class = "rasch_frame_invariance"), fit),
    error = conditionMessage)
  expect_match(msg2, "must be a frame_invariance\\(\\) result")
  expect_false(grepl("atomic", msg2, fixed = TRUE))
})

test_that("a legacy many-facet truth is told what was actually tried", {
  d <- simulate_mfrm(60, 4, 5, seed = 5)
  attr(d, "truth")$person_id <- NULL
  names(attr(d, "truth")$theta) <- NULL
  fit <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                    facets = "rater")
  z <- sim_recovery(fit, d)
  expect_null(z$pieces[["person ability"]])
  expect_match(z$note, "no response pattern")
  # no pattern pairing runs for this layout, and rasch_mfrm() has no id=
  expect_false(grepl("uniquely pair", z$note, fixed = TRUE))
  expect_false(grepl("id=", z$note, fixed = TRUE))
  expect_match(z$note, "Re-simulate with the current simulator")
})

test_that("a legacy wide truth still reports the response-pattern pairing", {
  d <- simulate_rasch(200, 5, seed = 15)
  attr(d, "truth")$person_id <- NULL
  X <- as.matrix(d[rev(seq_len(nrow(d))),
                   names(attr(d, "truth")$difficulty)])
  rownames(X) <- NULL
  z <- sim_recovery(rasch(X), d)
  expect_null(z$pieces[["person ability"]])
  expect_match(z$note, "response patterns do not uniquely pair")
})

test_that("the B = 0 note reads as sentences whichever refusal it quotes", {
  # a refusal that closes with its own remedy sentence
  scoring <- structure(list(refit_spec = list(fixed_calibration = TRUE)),
                       class = c("rasch", "list"))
  note <- .dim_reference_advice(scoring)
  expect_identical(note, paste(
    "No bootstrap reference is available for this fit (this fully anchored",
    "calibration is a scoring fit; downstream recalibration and bootstrap",
    "procedures are not supported), so this comparison stays descriptive.",
    "Use the original calibration for a new analysis."))
  # no sentence closes inside the parentheses
  expect_false(grepl(".)", note, fixed = TRUE))
  # a refusal that is one sentence is quoted whole and closes the note
  stale <- structure(list(est = list(anchors = data.frame(
    item = "I1", location = 0))), class = c("rasch", "list"))
  note2 <- .dim_reference_advice(stale)
  expect_identical(note2, paste(
    "No bootstrap reference is available for this fit (the fitted calibration",
    "has anchors but its saved anchor settings are unavailable; refit from",
    "the source data with the original anchors before recalibration or",
    "bootstrapping), so this comparison stays descriptive."))
  expect_false(grepl(".)", note2, fixed = TRUE))
})

test_that("the verdict note of a scoring fit carries the repaired shape", {
  f <- rasch(simulate_rasch(300, 12, seed = 926))
  f$refit_spec$fixed_calibration <- TRUE
  z <- dimensionality_test(f, items_positive = sprintf("I%02d", 1:6),
                           items_negative = sprintf("I%02d", 7:12))
  expect_true(is.na(z$multidimensional))
  expect_match(z$verdict_note,
               "Use the original calibration for a new analysis\\.$")
  expect_false(grepl(".)", z$verdict_note, fixed = TRUE))
  expect_output(print(z), "Use the original calibration", fixed = TRUE)
})
