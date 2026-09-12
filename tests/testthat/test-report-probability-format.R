test_that("frame unit report probabilities are formatted before relabelling", {
  skip_if_not_installed("knitr")
  template <- testthat::test_path("..", "..", "inst", "rmarkdown",
                                 "rasch-report.Rmd")
  if (!file.exists(template))
    template <- system.file("rmarkdown", "rasch-report.Rmd", package = "rasch")
  lines <- readLines(template, warn = FALSE)
  start <- which(lines == "```{r setup, include=FALSE}")
  end <- which(seq_along(lines) > start & lines == "```")[1L]
  setup <- as.list(parse(text = lines[(start + 1L):(end - 1L)]))
  assignment <- Filter(function(e)
    is.call(e) && identical(e[[1L]], as.name("<-")) &&
      identical(e[[2L]], as.name("show_table")), setup)
  expect_length(assignment, 1L)
  env <- new.env(parent = asNamespace("rasch"))
  eval(assignment[[1L]], env)
  # Exercise the actual template formatter, not a duplicate implementation.
  tab <- data.frame(set = c("A", "B"), p_alpha_adj = c(1e-8, .04),
                     p_kappa_adj = c(2e-9, NA_real_))
  rendered <- paste(env$show_table(tab), collapse = "\n")
  expect_match(rendered, "adjusted p (alpha)", fixed = TRUE)
  expect_match(rendered, "adjusted p (kappa)", fixed = TRUE)
  expect_equal(lengths(regmatches(rendered,
    gregexpr("< 0.001", rendered, fixed = TRUE))), 2L)
  expect_match(rendered, "0.040", fixed = TRUE)
  text <- paste(lines, collapse = "\n")
  expect_match(text, 'names(au)[names(au) == "p_adj.x"] <- "p_alpha_adj"',
               fixed = TRUE)
  expect_match(text, 'names(au)[names(au) == "p_adj.y"] <- "p_kappa_adj"',
               fixed = TRUE)
})

test_that("report text escapes brackets for the dialect it renders in", {
  skip_if_not_installed("knitr")
  template <- testthat::test_path("..", "..", "inst", "rmarkdown",
                                 "rasch-report.Rmd")
  if (!file.exists(template))
    template <- system.file("rmarkdown", "rasch-report.Rmd", package = "rasch")
  lines <- readLines(template, warn = FALSE)
  start <- which(lines == "```{r setup, include=FALSE}")
  end <- which(seq_along(lines) > start & lines == "```")[1L]
  setup <- as.list(parse(text = lines[(start + 1L):(end - 1L)]))
  assignment <- Filter(function(e)
    is.call(e) && identical(e[[1L]], as.name("<-")) &&
      identical(e[[2L]], as.name("report_text")), setup)
  expect_length(assignment, 1L)
  env <- new.env(parent = asNamespace("rasch"))
  eval(assignment[[1L]], env)
  # Exercise the actual template escape, not a duplicate implementation.
  out <- env$report_text(c(
    "withheld because set-link uncertainty was omitted (boot_reps = 0)",
    "a [bracket] and a literal &#40; and <em>markup</em> and *emphasis*"))
  # The report renders under pandoc's tex_math_single_backslash dialect, in
  # which a backslash before a bracket opens mathematics instead of escaping
  # it: the note would be set as an equation with its brackets deleted.
  expect_false(any(grepl("\\(", out, fixed = TRUE)))
  expect_false(any(grepl("\\[", out, fixed = TRUE)))
  expect_match(out[1], "&#40;boot\\_reps = 0&#41;", fixed = TRUE)
  expect_match(out[2], "&#91;bracket&#93;", fixed = TRUE)
  # a character reference in the data stays text rather than becoming one
  expect_match(out[2], "&amp;\\#40;", fixed = TRUE)
  expect_match(out[2], "&lt;em&gt;markup&lt;/em&gt;", fixed = TRUE)
  expect_match(out[2], "\\*emphasis\\*", fixed = TRUE)
})
