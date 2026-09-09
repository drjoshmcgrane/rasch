test_that("save_outputs refuses non-empty destinations before writing", {
  fit <- rasch(simulate_rasch(60, 4, seed = 981))
  nonempty <- tempfile("rasch-export-nonempty-")
  dir.create(nonempty)
  writeLines("keep", file.path(nonempty, ".keep"))
  nested <- file.path(nonempty, "nested")
  dir.create(nested)
  on.exit(unlink(nonempty, recursive = TRUE), add = TRUE)

  expect_error(
    save_outputs(fit, nonempty, formats = "png", item_plots = FALSE),
    "new or empty directory")
  expect_true(file.exists(file.path(nonempty, ".keep")))
  expect_false(file.exists(file.path(nonempty, "tables")))

  existing_file <- tempfile("rasch-export-file-")
  writeLines("not a directory", existing_file)
  on.exit(unlink(existing_file), add = TRUE)
  expect_error(save_outputs(fit, existing_file, formats = "png",
                            item_plots = FALSE), "existing file")

  bt_data <- simulate_btl(4, 8, 5, seed = 983)
  bt <- btl(bt_data, "object_a", "object_b", winner = "winner",
            judge = "judge")
  btl_nonempty <- tempfile("rasch-btl-export-nonempty-")
  dir.create(btl_nonempty)
  writeLines("keep", file.path(btl_nonempty, ".keep"))
  on.exit(unlink(btl_nonempty, recursive = TRUE), add = TRUE)
  expect_error(save_outputs(bt, btl_nonempty, formats = "png",
                            item_plots = FALSE), "new or empty directory")
  expect_false(file.exists(file.path(btl_nonempty, "tables")))

  empty <- tempfile("rasch-export-empty-")
  dir.create(empty)
  on.exit(unlink(empty, recursive = TRUE), add = TRUE)
  expect_no_error(suppressWarnings(
    save_outputs(fit, empty, formats = "png", item_plots = FALSE)))
  expect_true(file.exists(file.path(empty, "tables", "fit_summary.csv")))

  exported <- tempfile("rasch-export-reuse-")
  on.exit(unlink(exported, recursive = TRUE), add = TRUE)
  expect_no_error(suppressWarnings(
    save_outputs(fit, exported, formats = "png", item_plots = FALSE)))
  before <- sort(list.files(exported, recursive = TRUE, all.files = TRUE,
                            no.. = TRUE, full.names = TRUE))
  before_hash <- unname(tools::md5sum(before))
  expect_error(
    save_outputs(fit, exported, formats = "png", item_plots = FALSE),
    "new or empty directory")
  after <- sort(list.files(exported, recursive = TRUE, all.files = TRUE,
                           no.. = TRUE, full.names = TRUE))
  expect_identical(after, before)
  expect_identical(unname(tools::md5sum(after)), before_hash)
})

test_that("supplied external-factor DIF is included in exports", {
  set.seed(982)
  X <- matrix(rbinom(120 * 5, 1, .5), 120, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X)
  expect_null(fit$factors)
  factors <- data.frame(group = factor(rep(c("A", "B"), each = 60)))
  dif <- dif_anova(fit, factors = factors)
  out <- tempfile("rasch-external-dif-")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_no_error(suppressWarnings(
    save_outputs(fit, out, formats = "png", item_plots = FALSE, dif = dif)))
  expect_true(file.exists(file.path(out, "tables", "dif_anova.csv")))
  expect_true(file.exists(file.path(out, "tables", "dif_anova_terms.csv")))
  csv_summary <- read.csv(file.path(out, "tables", "dif_anova.csv"),
                          check.names = FALSE, stringsAsFactors = FALSE)
  csv_terms <- read.csv(file.path(out, "tables", "dif_anova_terms.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
  compare_export_table <- function(actual, expected) {
    expected <- as.data.frame(expected)
    expect_identical(names(actual), names(expected))
    expect_identical(nrow(actual), nrow(expected))
    for (nm in names(expected)) {
      if (is.numeric(actual[[nm]]) || is.numeric(expected[[nm]])) {
        expect_equal(as.numeric(actual[[nm]]), as.numeric(expected[[nm]]),
                     tolerance = 1e-12)
      } else {
        expect_equal(as.character(actual[[nm]]), as.character(expected[[nm]]))
      }
    }
  }
  compare_export_table(csv_summary, dif$summary)
  compare_export_table(csv_terms, dif$terms)

  html <- tempfile(fileext = ".html")
  on.exit(unlink(html), add = TRUE)
  expect_no_error(suppressWarnings(report_html(fit, html, dif = dif, dpi = 40)))
  html_text <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(html_text, "Differential item functioning", fixed = TRUE)

  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "Pandoc is required for report_document")
  dif_boot <- suppressWarnings(
    dif_bootstrap(fit, dif, B = 1L, workers = 1L, seed = 983))
  report <- tempfile(fileext = ".html")
  on.exit(unlink(report), add = TRUE)
  expect_no_error(suppressWarnings(
    report_document(fit, report, format = "html", dif = dif,
                    dif_bootstrap = dif_boot)))
  report_text <- paste(readLines(report, warn = FALSE), collapse = "\n")
  expect_match(report_text, "Differential item functioning", fixed = TRUE)
  expect_match(report_text, "Bootstrap sensitivity analysis", fixed = TRUE)
})
