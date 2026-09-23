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
  text_assignment <- Filter(function(e)
    is.call(e) && identical(e[[1L]], as.name("<-")) &&
      identical(e[[2L]], as.name("report_text")), setup)
  eval(text_assignment[[1L]], env)
  eval(assignment[[1L]], env)
  # Exercise the actual template formatter, not a duplicate implementation.
  tab <- data.frame(set = c("A", "B"), p_alpha_adj = c(1e-8, .04),
                     p_kappa_adj = c(2e-9, NA_real_))
  rendered <- paste(env$show_table(tab), collapse = "\n")
  expect_match(rendered, "adjusted p &#40;alpha&#41;", fixed = TRUE)
  expect_match(rendered, "adjusted p &#40;kappa&#41;", fixed = TRUE)
  expect_equal(lengths(regmatches(rendered,
    gregexpr("&lt; 0.001", rendered, fixed = TRUE))), 2L)
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
  # pandoc's default dialect reads ^x^ as superscript, ~x~ as subscript and
  # ~~x~~ as strikeout, so an item name that closes one of those spans is
  # escaped like the other inline markers, once
  sup <- env$report_text(c("x^2^", "H~2~O", "~~old~~"))
  expect_identical(sup, c("x\\^2\\^", "H\\~2\\~O", "\\~\\~old\\~\\~"))
})

test_that("report table labels remain literal after Pandoc conversion", {
  skip_if_not_installed("knitr"); skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  template <- testthat::test_path("..", "..", "inst", "rmarkdown", "rasch-report.Rmd")
  if (!file.exists(template))
    template <- system.file("rmarkdown", "rasch-report.Rmd", package = "rasch")
  lines <- readLines(template, warn = FALSE)
  start <- which(lines == "```{r setup, include=FALSE}")
  end <- which(seq_along(lines) > start & lines == "```")[1L]
  expr <- as.list(parse(text = lines[(start + 1L):(end - 1L)]))
  env <- new.env(parent = asNamespace("rasch"))
  for (e in expr) if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
    as.character(e[[2L]]) %in% c("report_text", "show_table")) eval(e, env)
  tab <- data.frame(label = c("a$b", "c$d", "item|set", "x^2^", "H~2~O",
                              "~~old~~"),
                    value = c(.4, .5, .6, .7, .8, .9))
  names(tab)[1] <- "item$label"
  md <- tempfile(fileext = ".md"); html <- tempfile(fileext = ".html")
  on.exit(unlink(c(md, html)), add = TRUE)
  writeLines(as.character(env$show_table(tab)), md)
  rmarkdown::pandoc_convert(md, to = "html", output = html,
    from = "markdown+tex_math_dollars+tex_math_single_backslash")
  out <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(out, "a$b", fixed = TRUE)
  expect_match(out, "c$d", fixed = TRUE)
  expect_match(out, "item$label", fixed = TRUE)
  expect_match(out, "item|set", fixed = TRUE)
  expect_match(out, "x^2^", fixed = TRUE)
  expect_match(out, "H~2~O", fixed = TRUE)
  expect_match(out, "~~old~~", fixed = TRUE)
  expect_false(grepl("<sup>", out, fixed = TRUE))
  expect_false(grepl("<sub>", out, fixed = TRUE))
  expect_false(grepl("<del>", out, fixed = TRUE))
  expect_false(grepl('class="math', out, fixed = TRUE))
  tex <- tempfile(fileext = ".tex")
  docx <- tempfile(fileext = ".docx")
  xml_dir <- tempfile(); dir.create(xml_dir)
  on.exit(unlink(c(tex, docx, xml_dir), recursive = TRUE), add = TRUE)
  rmarkdown::pandoc_convert(md, to = "latex", output = tex,
    from = "markdown+tex_math_dollars+tex_math_single_backslash")
  expect_match(paste(readLines(tex, warn = FALSE), collapse = "\n"),
               "a\\$b", fixed = TRUE)
  rmarkdown::pandoc_convert(md, to = "docx", output = docx,
    from = "markdown+tex_math_dollars+tex_math_single_backslash")
  utils::unzip(docx, files = "word/document.xml", exdir = xml_dir)
  xml <- paste(readLines(file.path(xml_dir, "word/document.xml"),
                         warn = FALSE), collapse = "\n")
  expect_match(xml, "a$b", fixed = TRUE)
  expect_match(xml, "x^2^", fixed = TRUE)
  expect_match(xml, "H~2~O", fixed = TRUE)
  expect_false(grepl("vertAlign", xml, fixed = TRUE))
  expect_false(grepl("<w:strike", xml, fixed = TRUE))
  expect_false(grepl("<m:oMath", xml, fixed = TRUE))
})
