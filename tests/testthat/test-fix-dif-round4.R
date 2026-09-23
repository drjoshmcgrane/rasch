# Regressions for the round-four DIF fixes: every withheld row names the
# cause that actually applies to it, and the class-interval nuisance row is
# counted and worded as the nuisance term it is.

# 120 persons in three groups, the third seen on one occasion only, so the
# panels are unequal and the joint CR3 branch takes over. Level "c" occupies
# the top class interval alone, which aliases one of the ci contrasts while
# leaving the ci term two Type II degrees of freedom: the ci row is withheld.
sim_withheld_ci <- function() {
  set.seed(77)
  np <- 120L
  a0 <- factor(rep(c("a", "b", "c"), each = np / 3),
               levels = c("a", "b", "c"))
  pid <- seq_len(np)
  ids <- c(pid, pid[1:100])
  occ <- factor(c(rep("pre", np), rep("post", 100)),
                levels = c("pre", "post"))
  x <- matrix(rbinom(length(ids) * 4, 1, 0.5), length(ids), 4,
              dimnames = list(NULL, paste0("I", 1:4)))
  fit <- rasch(x, id = ids,
               factors = data.frame(A = a0[ids], occasion = occ))
  th <- ifelse(a0 == "c", 5, seq(-2, 2, length.out = np))
  fit$person$theta <- th[match(fit$person$id, pid)]
  fit$person$extreme <- FALSE
  fit
}

test_that("the two within-stratum refusals name different causes", {
  set.seed(3)
  # occupied between cells are (a, p), (b, p) and (c, q) only: every factor
  # keeps all its levels, but the empty cells alias the balanced grand mean
  f1 <- factor(rep(c("a", "b", "c"), each = 8))
  f2 <- factor(ifelse(f1 == "c", "q", "p"))
  n <- length(f1)
  y <- matrix(rnorm(n * 2), n, 2,
              dimnames = list(paste0("p", seq_len(n)), c("t1", "t2")))
  pd <- data.frame(pid = paste0("p", seq_len(n)), f1 = f1, f2 = f2)
  expect_equal(nlevels(droplevels(f1)), 3L)
  expect_equal(nlevels(droplevels(f2)), 2L)
  xb <- model.matrix(~ f1 + f2, pd,
                     contrasts.arg = list(f1 = "contr.sum", f2 = "contr.sum"))
  expect_equal(qr(xb)$rank, ncol(xb) - 1L)

  aliased <- rasch:::.dif_within_tests(y, pd, "occ", list(occ = 2L),
                                       "occ", c("f1", "f2"))
  why <- aliased$unavailable_reason[aliased$term == "occ"]
  expect_true(is.na(aliased$F_value[aliased$term == "occ"]))
  expect_match(why, "empty cell of the between-person design")
  expect_match(why, "balanced grand mean")
  expect_false(grepl("two levels of every factor", why, fixed = TRUE))

  # the other refusal: g2 keeps a single level in the complete panels, so a
  # term built on it is not estimable at all
  pd2 <- data.frame(pid = pd$pid, g1 = rep(c("A", "B"), 12),
                    g2 = factor("only"))
  lost <- rasch:::.dif_within_tests(y, pd2, "occ", list(occ = 2L),
                                    "occ:g2", c("g1", "g2", "g1:g2"))
  why2 <- lost$unavailable_reason[lost$term == "occ:g2"]
  expect_match(why2, "two levels of every factor in this term")
  expect_false(identical(why, why2))
})

test_that("a fully fitted person cluster is not blamed on the term", {
  set.seed(4)
  n <- 60L
  # level "c" holds one person, whose single row the design fits exactly:
  # its delete-cluster residual, not the term's covariance, is what fails
  f1 <- factor(c(rep("a", 30), rep("b", 29), "c"))
  d <- data.frame(z = rnorm(n), f1 = f1)
  got <- rasch:::.dif_type2(d, "f1", variance = "cr3",
                            cluster = factor(seq_len(n)),
                            weights = rep(1, n))
  why <- got$unavailable_reason[got$term == "f1"]
  expect_true(is.na(got$F_value[got$term == "f1"]))
  expect_match(why, "delete-cluster residual")
  expect_false(grepl("covariance of the tested term", why, fixed = TRUE))
})

test_that("a withheld class-interval row is not counted as a DIF test", {
  skip_on_cran()
  da <- dif_anova(sim_withheld_ci(), factors = c("A", "occasion"),
                  within = "occasion", n_groups = 4, effects = "main")
  tab <- as.data.frame(da$terms)
  withheld_ci <- tab$term == "ci" & !is.finite(tab$F_value)
  expect_gt(sum(withheld_ci), 0L)               # the row this test is about

  rows <- grep("unavailable because", da$notes, value = TRUE)
  ci_rows <- grep("[ci]:", rows, fixed = TRUE, value = TRUE)
  expect_equal(length(ci_rows), sum(withheld_ci))
  # the class-interval row is told what it is, not that the family keeps it
  expect_true(all(grepl("not a requested DIF test", ci_rows, fixed = TRUE)))
  expect_false(any(grepl("adjustment family", ci_rows, fixed = TRUE) &
                     !grepl("joins no adjustment family", ci_rows,
                            fixed = TRUE)))

  # the two integers on the Notes lines match the per-row notes they cover
  req <- grep("requested item-term test\\(s\\) were not estimable", da$notes,
              value = TRUE)
  nui <- grep("class-interval main-effect row\\(s\\) carry no F", da$notes,
              value = TRUE)
  expect_length(req, 1L)
  expect_length(nui, 1L)
  n_req <- as.integer(sub("^([0-9]+) .*$", "\\1", req))
  n_nui <- as.integer(sub("^([0-9]+) .*$", "\\1", nui))
  expect_equal(n_nui, length(ci_rows))
  expect_equal(n_req, length(rows) - length(ci_rows))
  expect_equal(n_req, sum(grepl("adjustment family", rows, fixed = TRUE) &
                            !grepl("joins no adjustment family", rows,
                                   fixed = TRUE)))
  # and the aggregate no longer promises the family keeps every withheld row
  expect_match(req, "all of them remain in the adjusted-probability family")
})

test_that("dif_anova documents the class interval's place in the family", {
  rd <- test_path("..", "..", "man", "dif_anova.Rd")
  skip_if_not(file.exists(rd), "man/ is not in the checked tree")
  txt <- gsub("[[:space:]]+", " ", paste(readLines(rd, warn = FALSE),
                                         collapse = " "))
  expect_false(grepl(paste("Every withheld test is named in \\code{notes}",
                           "with its reason and stays in the multiplicity",
                           "family"), txt, fixed = TRUE))
  expect_match(txt, paste("withheld class-interval row is a nuisance term,",
                          "never a member of it"), fixed = TRUE)
})
