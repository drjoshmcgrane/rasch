# Round-four repair: a paired-comparison bank with no objects, the twin of
# the empty equating reference round three refused.

test_that("a bank with no objects is refused by name", {
  expect_error(.btl_equate_ref(data.frame(object = character(),
                                          location = numeric())),
               "a bank has no objects")
  # the same refusal whichever optional column the caller supplied: the
  # length-one default below the check never sees a zero-row frame
  expect_error(.btl_equate_ref(data.frame(object = character(),
                                          location = numeric(),
                                          se = numeric())),
               "a bank has no objects")
  # a header-only CSV is the shape that reaches this in practice
  expect_error(.btl_equate_ref(utils::read.csv(text = "object,location\n")),
               "a bank has no objects")
})

test_that("btl_equate and the project loader report the empty bank", {
  d <- simulate_btl(8, 20, reps_per_pair = 10, seed = 9681)
  fit <- btl(d, "object_a", "object_b", winner = "winner")
  empty <- data.frame(object = character(), location = numeric())
  msg <- tryCatch(btl_equate(fit, empty), error = conditionMessage)
  expect_match(msg, "a bank has no objects")
  expect_false(grepl("replacement has 1 row", msg, fixed = TRUE))
  # the saved bank records the same refusal, not the implementation message
  prob <- .app_project_resource_problem("bt_eq_bank", empty)
  expect_match(prob, "a bank has no objects")
  expect_false(grepl("replacement has 1 row", prob, fixed = TRUE))
})

test_that("the resolve_dif return documentation wraps at the file's width", {
  src_file <- testthat::test_path("..", "..", "R", "subtests.R")
  skip_if_not(file.exists(src_file), "sources are not installed")
  src <- readLines(src_file)
  roxygen <- grep("^#'", src, value = TRUE)
  expect_equal(roxygen[nchar(roxygen) > 100L], character(0))
})
