# Input honesty: misspelled column names, fractional scores, missing
# identifiers, unusable equating items, and HTML escaping. Each of these
# used to fail silently (wrong-data analyses, truncated responses, phantom
# "NA" facet levels, all-NA equating, markup injection in reports).

mkX <- function(n = 150, L = 6, seed = 1) {
  set.seed(seed)
  X <- matrix(rbinom(n * L, 1, plogis(outer(rnorm(n), seq(-1, 1, length.out = L), "-"))),
              n, L, dimnames = list(NULL, paste0("I", 1:L)))
  X
}

test_that("misspelled id, factor, and item columns are errors, not fallbacks", {
  d <- as.data.frame(mkX())
  d$sex <- rep(c("m", "f"), length.out = nrow(d))
  d$pid <- sprintf("P%03d", seq_len(nrow(d)))
  expect_error(rasch(d, id = "person_id", items = paste0("I", 1:6)),
               "id column 'person_id' not found")
  expect_error(rasch(d, factors = c("sex", "agee"), items = paste0("I", 1:6)),
               "factor column\\(s\\) not found.*agee")
  expect_error(rasch(d, items = c("I1", "I2", "Item3")),
               "item column\\(s\\) not found.*Item3")
  # correct names still work, including numeric item indices
  f <- rasch(d, id = "pid", factors = "sex", items = paste0("I", 1:6))
  expect_equal(ncol(f$X), 6L)
  f2 <- rasch(as.data.frame(mkX()), items = 1:6)
  expect_equal(ncol(f2$X), 6L)
})

test_that("fractional scores error instead of silently truncating", {
  X <- mkX()
  Xf <- as.data.frame(X)
  Xf$I3[5] <- 1.9
  expect_error(rasch(Xf), "non-integer score\\(s\\) in: I3.*1\\.9")
  # integer-valued doubles ("2.0") are fine
  Xd <- as.data.frame(X * 1.0)
  expect_s3_class(rasch(Xd), "rasch")
})

test_that("MFRM rows with missing identifiers are dropped with a note", {
  set.seed(3)
  long <- data.frame(person = rep(sprintf("P%03d", 1:80), each = 4),
                     rater = rep(c("R1", "R2"), 160),
                     item = rep(rep(c("A", "B"), each = 2), 80),
                     score = rbinom(320, 2, 0.5))
  long$rater[c(5, 9)] <- NA
  long$person[17] <- NA
  f <- rasch_mfrm(long, person = "person", item = "item", score = "score",
                  facets = "rater")
  expect_true(any(grepl("3 row\\(s\\) dropped: missing person, item, or facet",
                        f$notes)))
  # no phantom "NA" rater level
  expect_false("NA" %in% f$facet_effects$level)
})

test_that("equating excludes unusable items instead of returning all NA", {
  f <- rasch(as.data.frame(mkX(400, 8)))
  ref <- data.frame(item = paste0("I", 1:8),
                    location = f$items$location + 0.3,
                    se = f$items$se)
  ref$se[2] <- NA                       # e.g. a weakly determined item
  eq <- equate_tests(f, ref)
  expect_true(is.finite(eq$shift) && is.finite(eq$rmsd))
  expect_equal(eq$n, 7L)
  expect_true(is.na(eq$table$t[eq$table$item == "I2"]))
  expect_equal(sum(is.finite(eq$table$t)), 7L)
  expect_match(eq$note, "I2")
  expect_no_error(plot_equate(f, ref))
  ref$se[1:7] <- NA                     # fewer than two usable -> error
  expect_error(equate_tests(f, ref), "fewer than two common items")
})

test_that("report_html escapes data-derived text", {
  X <- mkX(200, 6)
  colnames(X) <- c(paste0("I", 1:5), "A<b>&x")
  f <- rasch(as.data.frame(X))
  out <- tempfile(fileext = ".html")
  report_html(f, out, title = "T <script>alert(1)</script>")
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_false(grepl("<script>alert", html, fixed = TRUE))
  expect_true(grepl("&lt;script&gt;alert", html, fixed = TRUE))
  expect_false(grepl("A<b>&x", html, fixed = TRUE))
  expect_true(grepl("A&lt;b&gt;&amp;x", html, fixed = TRUE) ||
              grepl("A&lt;b>&amp;x", html, fixed = TRUE))
  unlink(out)
})
