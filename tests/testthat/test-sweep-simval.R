# The simulation-validation battery lives in tools/, which is not installed
# with the package, so every check here needs a source tree to run against.
.simval_path <- function(...) testthat::test_path("..", "..", "tools",
                                                  "simval", ...)

# A study that loads a prefix of another file pins that file's content hash.
# Every such pin must name the committed file, or the study aborts at its
# first guard for everyone: check them all, not only the one last repaired.
test_that("every study pins the committed hash of the file it loads", {
  root <- testthat::test_path("..", "..")
  dir <- file.path(root, "tools", "simval", "studies")
  skip_if_not(dir.exists(dir))
  pins <- 0L
  for (f in list.files(dir, pattern = "[.]R$", full.names = TRUE)) {
    code <- readLines(f)
    for (i in grep('^[A-Za-z._]+_md5 <- "[0-9a-f]{32}"', code)) {
      pinned <- sub('^[A-Za-z._]+_md5 <- "([0-9a-f]{32})".*$', "\\1", code[i])
      paths <- grep('^[A-Za-z._]+ <- "tools/simval/[^"]+[.]R"$',
                    code[seq_len(i)], value = TRUE)
      expect_length(paths, 1L)
      target <- file.path(root, sub('^[^"]*"([^"]+)".*$', "\\1",
                                    paths[length(paths)]))
      expect_true(file.exists(target))
      expect_identical(pinned, unname(tools::md5sum(target)))
      pins <- pins + 1L
    }
  }
  expect_gt(pins, 0L)
})

test_that("the loaded prefixes still define what their callers use", {
  helper <- .simval_path("studies", "dif-conditional-bootstrap-extended.R")
  screen <- .simval_path("studies", "dif-score-purification.R")
  skip_if_not(file.exists(helper) && file.exists(screen))
  assigned <- function(file, n) vapply(parse(file)[seq_len(n)], function(e)
    if (is.call(e) && identical(as.character(e[[1L]]), "<-"))
      as.character(e[[2L]]) else "", character(1))
  # The confirmation evaluates the helper's first nineteen expressions, and
  # the refined purification study the screen's first fifteen; those must
  # still be the ones that define the simulation and its scenarios.
  expect_true(all(c("scenario", "one_replicate") %in% assigned(helper, 19L)))
  expect_true("anchor_state" %in% assigned(screen, 15L))
})

test_that("the interval checker reports the R-source hash, not gates on it", {
  checker <- .simval_path("check-item-fit-intervals.R")
  results <- .simval_path("results", paste0("item-fit-bootstrap-intervals",
    c("", "-attempts", "-items"), ".csv"))
  root <- testthat::test_path("..", "..")
  skip_if_not(all(file.exists(c(checker, results))) &&
                dir.exists(file.path(root, "R")))
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  out <- capture.output(source("tools/simval/check-item-fit-intervals.R",
                               local = new.env(parent = globalenv())))
  expect_true(any(grepl("^All 200 datasets:", out)))
  expect_true(any(grepl("^R source tree ", out)))
})

test_that("the scree category-guard results come from the committed script", {
  study <- .simval_path("studies", "scree-conditional-reference.R")
  guard <- .simval_path("results",
                        "scree-conditional-reference-category-guard.csv")
  readme <- .simval_path("README.md")
  skip_if_not(all(file.exists(c(study, guard, readme))))
  g <- read.csv(guard)
  # Pin the script the rows name, not its hash: the rows come from a
  # 200-replicate run that no test can redo, so hashing the study would fail
  # the suite on the next edit to it. Check instead that the study is the
  # committed one and still classifies the way the rows were classified.
  expect_identical(unique(g$script),
                   "tools/simval/studies/scree-conditional-reference.R")
  src <- readLines(study)
  expect_true(any(grepl("conditional-reference support guard", src,
                        fixed = TRUE)))
  expect_true(any(grepl("rasch_fit_bootstrap_refusal", src, fixed = TRUE)))
  expect_true(any(grepl("empty_result(refused = 1, reference = reference)",
                        src, fixed = TRUE)))
  # A replicate the support guard rejects is a refusal, not an error, and a
  # wholly refused cell reports the guard rather than a rejection rate.
  s <- g[g$scenario == "sparse PCM, null", ]
  expect_identical(nrow(s), 1L)
  expect_identical(s$quantity, "conditional-reference support guard")
  expect_identical(s$n_reps, 0L)
  expect_gt(s$n_refused, 0L)
  expect_true(any(grepl("category-guard", readLines(readme), fixed = TRUE)))
})
