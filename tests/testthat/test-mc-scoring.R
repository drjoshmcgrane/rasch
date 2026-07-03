# Multiple-choice scoring: double keying and polytomous option scoring
# (Andrich & Styles 2011), with the distractor_rescore proposal.

sim_mc_partial <- function(N = 800, L = 6, seed = 4) {
  # trait-driven 3-category items relabelled as options: A = full credit,
  # B = informative distractor, D = uninformative wrong, (C = rare noise)
  set.seed(seed)
  th <- rnorm(N)
  d0 <- seq(-0.8, 0.8, length.out = L)
  raw <- sapply(d0, function(d) {
    x <- vapply(th, function(b) {
      p <- item_moments(b, c(d - 0.8, d + 0.8))$P
      sample(0:2, 1, prob = p)
    }, 0L)
    o <- c("D", "B", "A")[x + 1]
    swap <- runif(N) < 0.05 & o == "D"
    o[swap] <- "C"
    o
  })
  colnames(raw) <- paste0("M", seq_len(L))
  list(raw = raw, th = th)
}

test_that("double keying credits every listed option", {
  set.seed(9); N <- 400
  th <- rnorm(N)
  raw <- sapply(seq(-1, 1, length.out = 5), function(d) {
    ok <- rbinom(N, 1, plogis(th - d))
    ifelse(ok == 1, sample(c("A", "C"), N, replace = TRUE),
           sample(c("B", "D"), N, replace = TRUE))
  })
  colnames(raw) <- paste0("M", 1:5)
  fit <- rasch(raw, key = setNames(rep("A/C", 5), colnames(raw)))
  expect_true(all(fit$m == 1))
  # both A and C score 1
  expect_equal(unname(fit$X[, 1]),
               as.integer(raw[, 1] %in% c("A", "C")))
  da <- distractor_analysis(fit)
  expect_true(all(da$keyed[da$option %in% c("A", "C")]))
  expect_true(all(!da$keyed[da$option %in% c("B", "D")]))
  expect_equal(unname(fit$mc$key[1]), "A/C")
})

test_that("polytomous option scoring fits credited distractors as categories", {
  s <- sim_mc_partial()
  os <- do.call(rbind, lapply(colnames(s$raw), function(it)
    data.frame(item = it, option = c("A", "B"), score = c(2, 1))))
  fit <- rasch(s$raw, key = os)
  expect_true(all(fit$m == 2))
  expect_equal(unname(fit$X[, 3]),
               unname(c(A = 2L, B = 1L, C = 0L, D = 0L)[s$raw[, 3]]))
  # the polytomous scoring recovers the trait better than binary scoring
  bin <- rasch(s$raw, key = setNames(rep("A", 6), colnames(s$raw)))
  expect_gt(cor(fit$person$theta, s$th, use = "complete.obs"),
            cor(bin$person$theta, s$th, use = "complete.obs"))
  # display forms
  expect_equal(unname(fit$mc$key[1]), "A=2, B=1")
  da <- distractor_analysis(fit)
  expect_equal(da$score[da$item == "M1" & da$option == "B"], 1L)
  expect_true(da$keyed[da$item == "M1" & da$option == "A"])
})

test_that("distractor_rescore proposes credit for the informative distractor", {
  s <- sim_mc_partial()
  bin <- rasch(s$raw, key = setNames(rep("A", 6), colnames(s$raw)))
  pr <- distractor_rescore(bin)
  os <- pr$option_scores
  # B credited on (nearly) all items; C and D never; A keeps top score
  b <- os[os$option == "B", "score"]
  expect_gte(mean(b > 0), 5 / 6)
  expect_true(all(os[os$option %in% c("C", "D"), "score"] == 0))
  for (it in unique(os$item)) {
    d <- os[os$item == it, ]
    expect_equal(d$score[d$option == "A"], max(d$score))
  }
  # the proposal covers every observed option and feeds rasch directly
  refit <- rasch(s$raw, key = os)
  expect_s3_class(refit, "rasch")
  expect_true(all(refit$m >= 1))
  expect_true(is.finite(refit$psi$PSI))
  # evidence table carries the separation statistic
  expect_true(all(c("z_sep", "proposed", "se_location") %in% names(pr$evidence)))
})

test_that("key validation guards remain informative", {
  raw <- matrix(sample(c("A", "B"), 60, TRUE), 20, 3,
                dimnames = list(NULL, paste0("M", 1:3)))
  expect_error(rasch(raw, key = data.frame(item = "M1", option = "A",
                                           score = -1)),
               "non-negative")
  expect_error(rasch(raw, key = data.frame(item = "M1", option = c("A", "A"),
                                           score = c(1, 2))),
               "duplicate")
  expect_error(rasch(raw, key = data.frame(item = "M1", option = "A",
                                           score = 0)),
               "credits no option")
  expect_error(rasch(raw, key = data.frame(item = "M9", key = "A")),
               "no key item matches")
})
