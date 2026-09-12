# pkgload::load_all() does not regenerate NAMESPACE, so the freshly added
# print method is not yet in the S3 table; register it for the test session
# exactly as the roxygen @export tag will once NAMESPACE is regenerated.
registerS3method("print", "rasch_btl_equate", print.rasch_btl_equate)

test_that("btl_equate recovers the origin shift and flags nothing when no object drifts", {
  set.seed(20)
  objs <- paste0("O", 1:10)
  beta <- setNames(seq(-2, 2, length.out = 10), objs)

  # two panels judging overlapping object sets: 10 objects in all, 7 common.
  # each panel's sum-zero origin is the mean of its OWN set, so the two
  # origins differ and equating must recover the difference.
  set1 <- paste0("O", 1:9)                 # unique: O8, O9
  set2 <- c(paste0("O", 1:7), "O10")       # unique: O10
  common <- intersect(set1, set2)
  expect_equal(length(common), 7L)

  sim_panel <- function(objset, betas, reps = 60) {
    pr <- t(utils::combn(objset, 2))
    d <- data.frame(a = rep(pr[, 1], each = reps),
                    b = rep(pr[, 2], each = reps), stringsAsFactors = FALSE)
    p <- plogis(betas[d$a] - betas[d$b])
    d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
    btl(d, "a", "b", "win")
  }

  f1 <- sim_panel(set1, beta)
  f2 <- sim_panel(set2, beta)
  bad_convergence <- f1
  bad_convergence$converged <- NA
  expect_error(btl_equate(bad_convergence, f2), "did not converge")
  withheld <- btl_equate(f1, f2)
  expect_false(withheld$inferential)
  expect_true(all(is.na(withheld$table$p)))
  expect_match(paste(withheld$notes, collapse = " "), "independence")
  eq <- btl_equate(f1, f2, independent = TRUE)
  expect_error(btl_equate(f1, f2, independent = matrix(TRUE)),
               "NULL, TRUE, or FALSE")
  expect_error(btl_equate(f1, f2, independent = TRUE,
                          p_adjust = matrix("holm")),
               "p_adjust")

  expect_s3_class(eq, "rasch_btl_equate")
  expect_equal(eq$n_common, 7L)

  # sum-zero origins differ by mean(beta over set2) - mean(beta over set1);
  # loc2 + shift lands fit2 on fit1's scale, so shift ~ that difference
  expected_shift <- mean(beta[set2]) - mean(beta[set1])
  expect_lt(abs(eq$shift - expected_shift), 0.35)

  # essentially nothing drifts (at most the ~alpha rate)
  expect_lte(sum(eq$table$drifting), 1L)

  # fit2's equated locations track the truth (a pure shift of it) ~ perfectly
  truth <- beta[eq$equated$object]
  expect_gt(cor(eq$equated$location, truth), 0.97)

  bad_covariance <- f1
  bad_covariance$cov_beta[1, 1] <- -1e6
  guarded <- btl_equate(bad_covariance, f2, independent = TRUE)
  expect_false(guarded$inferential)
  expect_true(all(is.na(guarded$table$p_adj)))
  expect_match(paste(guarded$notes, collapse = " "),
               "not positive semidefinite")
})

test_that("btl_equate flags a planted drift and essentially only that object", {
  set.seed(6)
  objs <- paste0("O", 1:10)
  beta <- setNames(seq(-2, 2, length.out = 10), objs)
  set1 <- paste0("O", 1:9)
  set2 <- c(paste0("O", 1:7), "O10")

  sim_panel <- function(objset, betas, reps = 60) {
    pr <- t(utils::combn(objset, 2))
    d <- data.frame(a = rep(pr[, 1], each = reps),
                    b = rep(pr[, 2], each = reps), stringsAsFactors = FALSE)
    p <- plogis(betas[d$a] - betas[d$b])
    d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
    btl(d, "a", "b", "win")
  }

  # panel 2 values one common object (O4) 1.2 logits higher than panel 1 does
  beta2 <- beta
  beta2["O4"] <- beta2["O4"] + 1.2

  f1 <- sim_panel(set1, beta)
  f2 <- sim_panel(set2, beta2)
  eq <- btl_equate(f1, f2, independent = TRUE)

  drift <- setNames(eq$table$drifting, eq$table$object)
  expect_true(drift["O4"])                       # the planted object flags
  expect_equal(sum(drift[names(drift) != "O4"]), 0L)   # and only that one
  # the planted object dominates the drift statistics
  expect_equal(eq$table$object[which.max(abs(eq$table$t))], "O4")
})

test_that("clustered equating uses contrast-specific finite degrees of freedom", {
  objs <- paste0("O", 1:5)
  C <- 0.04 * (diag(5) - matrix(1 / 5, 5, 5))
  make_fit <- function(loc, prefix) structure(list(
    objects = data.frame(object = objs, location = loc,
                         se = sqrt(diag(C))),
    cov_beta = C, converged = TRUE, m = 1L, categories = 0:1,
    thr_structure = "none", clustered = TRUE,
    comparisons = data.frame(judge = rep(paste0(prefix, 1:12), each = 2))),
    class = "rasch_btl")
  f1 <- make_fit(c(-1, -0.4, 0, 0.5, 0.9), "A")
  f2 <- make_fit(c(-1, -0.4, 0.35, 0.5, 0.55), "B")
  eq <- btl_equate(f1, f2, independent = TRUE, p_adjust = "none")
  expect_true(all(is.finite(eq$table$df)))
  expect_equal(eq$table$df, rep(22, 5), tolerance = 1e-8)
  expect_equal(eq$table$p,
               2 * pt(-abs(eq$table$t), df = eq$table$df),
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(eq$table$p,
                                2 * pnorm(-abs(eq$table$t)))))
})

test_that("precision-weighted errors state their estimated weights", {
  objs <- paste0("O", 1:5)
  C <- 0.04 * (diag(5) - matrix(1 / 5, 5, 5))
  make_fit <- function(loc, prefix) structure(list(
    objects = data.frame(object = objs, location = loc,
                         se = sqrt(diag(C))),
    cov_beta = C, converged = TRUE, m = 1L, categories = 0:1,
    thr_structure = "none", clustered = TRUE,
    comparisons = data.frame(judge = rep(paste0(prefix, 1:12), each = 2))),
    class = "rasch_btl")
  f1 <- make_fit(c(-0.9, -0.3, 0.1, 0.4, 0.7), "A")
  f2 <- make_fit(c(-1.1, -0.5, 0.3, 0.5, 0.8), "B")
  eq <- btl_equate(f1, f2, independent = TRUE)
  expect_identical(eq$shift_method, "precision-weighted")
  expect_true(is.finite(eq$shift_se))
  # the weights are built from the calibrations' own estimated errors, so
  # every quadratic form in them conditions on estimated weights: the note
  # must name the drift and equated errors too, not the shift alone
  expect_true(any(is.finite(eq$table$p_adj)))
  expect_true(any(is.finite(eq$equated$se)))
  msg <- paste(eq$notes, collapse = " ")
  expect_match(msg, "Every error built from the estimated precision weights",
               fixed = TRUE)
  expect_match(msg, "the shift standard error", fixed = TRUE)
  expect_match(msg, "each drift contrast standard error", fixed = TRUE)
  expect_match(msg, "the equated location standard errors", fixed = TRUE)
  fixed_origin <- btl_equate(f1, f2, independent = TRUE, shift = "none")
  expect_identical(fixed_origin$shift_method, "none")
  expect_false(any(grepl("estimated precision weights", fixed_origin$notes,
                         fixed = TRUE)))
})

test_that("BTL-EFRM equating df follows its uncertainty method", {
  z <- structure(list(
    se_method = "bootstrap",
    comparisons = data.frame(judge = rep(paste0("J", 1:12), each = 2)),
    clustered = TRUE), class = c("rasch_btl_efrm", "rasch_btl"))
  expect_identical(.btl_equate_cov_df(z), Inf)
  z$se_method <- "conditional"
  expect_true(is.na(.btl_equate_cov_df(z)))
  z$se_method <- "judge_bootstrap"
  expect_equal(.btl_equate_cov_df(z), 11)
})

test_that("a bank link is descriptive without its joint covariance", {
  set.seed(22)
  objs <- paste0("O", 1:8)
  beta <- setNames(seq(-1.8, 1.8, length.out = 8), objs)
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 50),
                  b = rep(pr[, 2], each = 50), stringsAsFactors = FALSE)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  f1 <- btl(d, "a", "b", "win")

  bank <- data.frame(object = f1$objects$object,
                     location = f1$objects$location,
                     se = f1$objects$se, stringsAsFactors = FALSE)
  eq <- btl_equate(f1, bank)

  expect_lt(abs(eq$shift), 1e-8)
  expect_identical(eq$shift_method, "precision-weighted")
  expect_false(eq$inferential)
  expect_true(all(is.na(eq$table$drifting)))
  expect_match(paste(eq$notes, collapse = " "), "joint object-location covariance")
  expect_equal(eq$n_common, nrow(f1$objects))

  bank_no_se <- bank
  bank_no_se$se[] <- NA_real_
  eq_unweighted <- btl_equate(f1, bank_no_se)
  expect_identical(eq_unweighted$shift_method, "unweighted")
  expect_equal(eq_unweighted$n_inference, 0L)
  expect_equal(eq_unweighted$shift,
               mean(f1$objects$location - bank_no_se$location),
               tolerance = 1e-12)
  expect_true(all(is.na(eq_unweighted$table$p_adj)))
  expect_match(paste(eq_unweighted$notes, collapse = " "), "unweighted mean")
  expect_match(paste(eq_unweighted$notes, collapse = " "),
               "included in the descriptive shift")

  # A separately calibrated bank can carry the full covariance explicitly.
  attr(bank, "cov_location") <- f1$cov_beta
  eq_cov <- btl_equate(f1, bank, independent = TRUE)
  expect_true(eq_cov$inferential)
  expect_equal(sum(eq_cov$table$drifting), 0L)
})

test_that("btl_equate permits a two-object link but withholds drift tests", {
  set.seed(23)
  objs <- paste0("O", 1:6)
  beta <- setNames(seq(-1.5, 1.5, length.out = 6), objs)
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 40),
                  b = rep(pr[, 2], each = 40), stringsAsFactors = FALSE)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  f1 <- btl(d, "a", "b", "win")

  # a bank sharing only two objects with fit1
  bank2 <- data.frame(object = c("O1", "O2", "Z1", "Z2"),
                      location = c(-1, -0.5, 0.5, 1),
                      se = rep(0.2, 4), stringsAsFactors = FALSE)
  eq2 <- btl_equate(f1, bank2)
  expect_equal(eq2$n_common, 2L)
  expect_true(is.finite(eq2$shift))
  expect_false(eq2$inferential)
  expect_true(all(is.na(eq2$table$p_adj)))
  expect_match(paste(eq2$notes, collapse = " "), "at least three common")

  bank1 <- bank2[bank2$object != "O2", ]
  expect_error(btl_equate(f1, bank1), "at least two common")

  bank_dup <- rbind(
    data.frame(object = f1$objects$object, location = f1$objects$location,
               se = f1$objects$se),
    data.frame(object = f1$objects$object[1], location = 0, se = 0.2))
  expect_error(btl_equate(f1, bank_dup), "must be unique")
  expect_error(btl_equate(f1, bank_dup[!duplicated(bank_dup$object), ],
                           independent = 1), "NULL, TRUE, or FALSE")

  # non-btl fit1
  expect_error(btl_equate(42, bank2), "btl")
  # non-btl, non-bank fit2
  expect_error(btl_equate(f1, 42), "btl fit or a bank")
})

test_that("btl_equate refuses a malformed bank covariance degree of freedom", {
  set.seed(224)
  d <- simulate_btl(n_objects = 5, n_judges = 20,
                    reps_per_pair = 20, seed = 224)
  fit <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  bank <- fit$objects[, c("object", "location", "se")]
  attr(bank, "cov_location") <- fit$cov_beta
  attr(bank, "df_location") <- -1
  expect_error(btl_equate(fit, bank), "df_location.*positive")
  attr(bank, "cov_location") <- NULL
  expect_error(btl_equate(fit, bank), "df_location.*positive")
})

test_that("a polytomous bank declares its score scale as numeric", {
  set.seed(225)
  d <- simulate_btl(n_objects = 5, n_judges = 20,
                    reps_per_pair = 25, model = "polytomous",
                    n_categories = 3, seed = 225)
  fit <- btl(d, "object_a", "object_b", response = "response",
             judge = "judge")
  bank <- fit$objects[, c("object", "location", "se")]
  attr(bank, "m") <- factor(as.character(fit$m))
  expect_error(btl_equate(fit, bank), "attr\\(bank, 'm'\\)")
  attr(bank, "m") <- .Machine$integer.max + 1
  expect_error(btl_equate(fit, bank), "attr\\(bank, 'm'\\)")
})

test_that("bank covariance validation is invariant to uncertainty scale", {
  set.seed(226)
  d <- simulate_btl(n_objects = 5, n_judges = 20,
                    reps_per_pair = 10, seed = 226)
  fit <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  bank <- fit$objects[, c("object", "location", "se")]
  bank$se <- 1e-8
  attr(bank, "cov_location") <- diag(rep((2e-8)^2, nrow(bank)))
  expect_error(btl_equate(fit, bank), "standard errors must agree")
  attr(bank, "cov_location") <- diag(rep((1e-8)^2, nrow(bank)))
  expect_no_error(btl_equate(fit, bank))

  bank$se <- 1e-10
  Cbad <- diag(rep(1e-20, nrow(bank)))
  Cbad[1, 2] <- Cbad[2, 1] <- 2e-20
  attr(bank, "cov_location") <- Cbad
  expect_error(btl_equate(fit, bank), "positive semidefinite")

  Casym <- diag(rep(1e-20, nrow(bank)))
  Casym[1, 2] <- 5e-21
  attr(bank, "cov_location") <- Casym
  expect_error(btl_equate(fit, bank), "symmetric")
})

test_that("plot_btl_equate draws without error", {
  set.seed(24)
  objs <- paste0("O", 1:10)
  beta <- setNames(seq(-2, 2, length.out = 10), objs)
  set1 <- paste0("O", 1:9)
  set2 <- c(paste0("O", 1:7), "O10")

  sim_panel <- function(objset, betas, reps = 50) {
    pr <- t(utils::combn(objset, 2))
    d <- data.frame(a = rep(pr[, 1], each = reps),
                    b = rep(pr[, 2], each = reps), stringsAsFactors = FALSE)
    p <- plogis(betas[d$a] - betas[d$b])
    d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
    btl(d, "a", "b", "win")
  }
  f1 <- sim_panel(set1, beta)
  f2 <- sim_panel(set2, beta)

  pdf(NULL)
  on.exit(dev.off())
  res <- plot_btl_equate(f1, f2, independent = TRUE)
  expect_s3_class(res, "rasch_btl_equate")

  # print method exercises its formatting path too
  expect_output(print(res), "Common-object equating")
})
