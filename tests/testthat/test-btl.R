# Bradley-Terry-Luce paired comparisons: equivalence with the conditional
# Rasch model, parameter recovery, judge diagnostics, and guards.

sim_btl <- function(beta, n_per_pair, seed = 1, judges = NULL) {
  set.seed(seed)
  pr <- t(utils::combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = n_per_pair),
                  b = rep(pr[, 2], each = n_per_pair))
  p <- plogis(beta[d$a] - beta[d$b])
  d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
  if (!is.null(judges)) d$judge <- sample(judges, nrow(d), replace = TRUE)
  d
}

test_that("BTL is the conditional form of the dichotomous Rasch model", {
  # extract from Rasch data every pair response with exactly one correct:
  # 'the correct item beat the incorrect one'. The pairwise conditional
  # likelihood pcml() maximises is exactly the BTL likelihood of those
  # comparisons, so the estimates must agree up to solver tolerance
  # (with beta = -delta: winning means being easier).
  set.seed(4)
  d0 <- seq(-1.5, 1.5, length.out = 6)
  N <- 350
  X <- matrix(rbinom(N * 6, 1, plogis(outer(rnorm(N), d0, "-"))), N, 6)
  colnames(X) <- paste0("I", 1:6)
  cmp <- list()
  for (i in 1:5) for (j in (i + 1):6) {
    sel <- which(X[, i] + X[, j] == 1)
    if (!length(sel)) next
    cmp[[length(cmp) + 1L]] <- data.frame(
      person = sel,
      a = colnames(X)[i], b = colnames(X)[j],
      win = ifelse(X[sel, i] == 1, colnames(X)[i], colnames(X)[j]))
  }
  cmp <- do.call(rbind, cmp)
  bt <- btl(cmp, "a", "b", "win", judge = "person")
  pc <- pcml(X)
  delta <- vapply(1:6, function(i) mean(pc$thr$tau[pc$thr$item == i]), 0)
  loc <- bt$objects$location[match(colnames(X), bt$objects$object)]
  expect_equal(unname(loc), unname(-delta), tolerance = 1e-4)
  # person-clustered sandwich SEs agree with pcml's to a similar order
  se_pc <- vapply(1:6, function(i) {
    rows <- pc$thr$id[pc$thr$item == i]
    sqrt(mean(pc$cov_tau[rows, rows]))
  }, 0)
  se_bt <- bt$objects$se[match(colnames(X), bt$objects$object)]
  expect_equal(unname(se_bt), unname(se_pc), tolerance = 0.25)
})

test_that("BTL recovers simulated locations with calibrated fit", {
  beta <- c(A = -1.2, B = -0.4, C = 0, D = 0.5, E = 1.1)
  d <- sim_btl(beta, 60, seed = 7)
  ft <- btl(d, "a", "b", "win")
  expect_true(ft$converged)
  loc <- ft$objects$location[match(names(beta), ft$objects$object)]
  expect_gt(cor(loc, beta), 0.98)
  expect_equal(sum(loc), 0, tolerance = 1e-8)          # sum-zero
  expect_lt(abs(max(loc - (beta - mean(beta)))), 0.35) # centred recovery
  # model-true data: pairwise chi-square unremarkable, fit residuals tame
  expect_gt(ft$total_p, 0.01)
  expect_lt(max(abs(ft$objects$fit_resid), na.rm = TRUE), 3)
  expect_true(ft$osi$PSI > 0.7)
})

test_that("an erratic judge is flagged by the judge fit residual", {
  beta <- c(A = -1.2, B = -0.4, C = 0.2, D = 0.6, E = 0.8)
  d <- sim_btl(beta, 80, seed = 11, judges = paste0("J", 1:8))
  # judge J1 answers at random
  sel <- d$judge == "J1"
  d$win[sel] <- ifelse(runif(sum(sel)) < 0.5, d$a[sel], d$b[sel])
  ft <- btl(d, "a", "b", "win", judge = "judge")
  expect_true(ft$clustered)
  jt <- ft$judges
  expect_equal(which.max(jt$fit_resid), match("J1", jt$judge))
  expect_gt(jt$fit_resid[jt$judge == "J1"], 2.5)
  expect_lt(max(jt$fit_resid[jt$judge != "J1"]), 2.5)
})

test_that("ties, extremes, counts, and disconnection are handled", {
  beta <- c(A = -0.8, B = 0, C = 0.8)
  d <- sim_btl(beta, 40, seed = 3)
  d$win[1:5] <- "tie"
  expect_error(btl(d, "a", "b", "win", ties = "error"), "tie")
  ft_d <- btl(d, "a", "b", "win", ties = "drop")
  expect_true(any(grepl("tie", ft_d$notes)))
  ft_h <- btl(d, "a", "b", "win", ties = "half")
  expect_gt(ft_h$n_comparisons, ft_d$n_comparisons)
  expect_equal(sum(ft_h$observed_comparisons$weight), nrow(d))
  # The two half rows are one comparison for sandwich purposes. Giving every
  # original row its own judge reproduces that meat, apart from CR1.
  d$row_cluster <- sprintf("R%04d", seq_len(nrow(d)))
  ft_h_cluster <- btl(d, "a", "b", "win", judge = "row_cluster",
                      ties = "half")
  cr1 <- nrow(d) / (nrow(d) - 1)
  expect_equal(ft_h$objects$location, ft_h_cluster$objects$location,
               tolerance = 1e-10)
  expect_equal(ft_h$cov_beta, ft_h_cluster$cov_beta / cr1,
               tolerance = 1e-8)
  # counts replicate rows
  dc <- data.frame(a = c("A", "A", "B"), b = c("B", "C", "C"),
                   win = c("A", "C", "B"), k = c(30, 30, 30))
  dc2 <- rbind(data.frame(a = "A", b = "B", win = "B", k = 10),
               data.frame(a = "B", b = "C", win = "C", k = 10),
               data.frame(a = "A", b = "C", win = "A", k = 10), dc)
  ft_c <- btl(dc2, "a", "b", "win", count = "k")
  expect_equal(ft_c$n_comparisons, 120)
  # an undefeated object is removed with a note
  d2 <- sim_btl(beta, 30, seed = 9)
  d2$win[d2$a == "C" | d2$b == "C"] <- "C"
  expect_error(btl(d2, "a", "b", "win"), "three comparable")
  d3 <- sim_btl(c(beta, D = 0.2, E = -0.2), 30, seed = 9)
  d3$win[d3$a == "E" | d3$b == "E"] <- "E"
  ft3 <- btl(d3, "a", "b", "win")
  expect_true(any(grepl("response boundary", ft3$notes)))
  # the undefeated object is set aside from estimation but reported at an
  # extrapolated location with the extreme flag and no standard error
  e_row <- ft3$objects[ft3$objects$object == "E", ]
  expect_true(isTRUE(e_row$extreme))
  expect_true(is.na(e_row$se))
  expect_gt(e_row$location,
            max(ft3$objects$location[!ft3$objects$extreme]))
  expect_true("E" %in% c(ft3$observed_comparisons$object_a,
                         ft3$observed_comparisons$object_b))
  expect_false("E" %in% c(ft3$comparisons$object_a,
                          ft3$comparisons$object_b))
  # disconnected comparison graphs are refused with the components listed
  dd <- data.frame(a = c("A", "A", "C", "C"), b = c("B", "B", "D", "D"),
                   win = c("A", "B", "C", "D"))
  expect_error(btl(dd, "a", "b", "win"), "disconnected")
})

test_that("comparison histories key judges and objects structurally", {
  a <- c("A\rB", "B", "A\rB", "B")
  b <- c("C", "D", "D", "C")
  judge <- c("J", "J\rA", "J", "J\rA")
  x <- c(1, 0, 0, 1)
  got <- .btl_exposure(a, b, x, 1, judge, c(1, 1, 2, 2))
  a2 <- c("X", "Y", "X", "Y")
  judge2 <- c("U", "V", "U", "V")
  expected <- .btl_exposure(a2, b, x, 1, judge2, c(1, 1, 2, 2))
  expect_identical(got, expected)
})

test_that("plot_btl draws and print method runs", {
  beta <- c(A = -1, B = 0, C = 0.4, D = 0.6)
  ft <- btl(sim_btl(beta, 40, seed = 5), "a", "b", "win")
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_btl(ft))
  drawn <- NULL
  testthat::with_mocked_bindings(
    plot_btl(ft),
    .btl_equate_cov_df = function(...) 4,
    segments = function(x0, y0, x1, y1, ...) {
      drawn <<- list(lower = x0, upper = x1)
    },
    .package = "rasch")
  d <- ft$objects[order(ft$objects$location), ]
  critical <- qt(0.975, 4)
  expect_equal(drawn$lower, d$location - critical * d$se)
  expect_equal(drawn$upper, d$location + critical * d$se)
  expect_output(print(ft), "Bradley-Terry-Luce")
})

test_that("graded comparisons recover locations and symmetric thresholds", {
  set.seed(3)
  beta <- c(A = -1.2, B = -0.5, C = 0, D = 0.6, E = 1.1)
  beta <- beta - mean(beta)
  tau <- c(-1.4, 0, 1.4)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 60), b = rep(pr[, 2], each = 60),
                  judge = sample(sprintf("J%02d", 1:15), 600, TRUE))
  d$grade <- vapply(seq_len(nrow(d)), function(r) {
    p <- item_moments(beta[d$a[r]] - beta[d$b[r]], tau)$P
    sample(0:3, 1, prob = p)
  }, 0L)
  f <- btl(d, "a", "b", response = "grade", judge = "judge")
  expect_true(f$converged)
  expect_lt(abs(sum(f$objects$location)), 1e-8)
  expect_gt(cor(f$objects$location, beta[f$objects$object]), 0.99)
  # thresholds symmetric with mirrored SEs, middle fixed at zero
  expect_equal(f$thresholds$tau[1], -f$thresholds$tau[3])
  expect_equal(f$thresholds$tau[2], 0)
  expect_equal(f$thresholds$se[1], f$thresholds$se[3])
  expect_lt(abs(f$thresholds$tau[1] - (-1.4)), 3 * f$thresholds$se[1])
  expect_true(f$clustered)
  expect_true(all(c("obs_mean", "exp_mean") %in% names(f$pairs)))

  # presentation-order invariance: swap objects and reverse the grades
  d2 <- data.frame(a = d$b, b = d$a, grade = 3 - d$grade, judge = d$judge)
  f2 <- btl(d2, "a", "b", response = "grade", judge = "judge")
  expect_equal(f2$objects$location, f$objects$location, tolerance = 1e-6)
  expect_equal(f2$thresholds$tau, f$thresholds$tau, tolerance = 1e-6)

  # two categories reproduce the dichotomous path exactly
  set.seed(4)
  d$win01 <- rbinom(nrow(d), 1, plogis(beta[d$a] - beta[d$b]))
  d$winner <- ifelse(d$win01 == 1, d$a, d$b)
  fd <- btl(d, "a", "b", winner = "winner", judge = "judge")
  fg <- btl(d, "a", "b", response = "win01", judge = "judge")
  expect_equal(fg$objects$location, fd$objects$location, tolerance = 1e-6)
  expect_equal(fg$objects$se, fd$objects$se, tolerance = 1e-6)

  # the fitted point maximises the likelihood along feasible directions
  ll_of <- function(bv, tv) sum(vapply(seq_len(nrow(d)), function(r) {
    p <- item_moments(bv[d$a[r]] - bv[d$b[r]], tv)$P
    log(p[d$grade[r] + 1])
  }, 0))
  bhat <- setNames(f$objects$location, f$objects$object)
  that <- f$thresholds$tau
  ll0 <- ll_of(bhat, that)
  set.seed(9)
  worse <- 0L
  for (rep in 1:6) {
    db <- rnorm(5); db <- db - mean(db)
    dt1 <- rnorm(1)
    if (ll_of(bhat + 0.004 * db, that + 0.004 * c(dt1, 0, -dt1)) <=
        ll0 + 1e-9) worse <- worse + 1L
  }
  expect_equal(worse, 6L)
})

test_that("three graded categories give the Davidson ties structure", {
  set.seed(6)
  beta <- c(P = -0.8, Q = 0, R = 0.8)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 80), b = rep(pr[, 2], each = 80))
  d$grade <- vapply(seq_len(nrow(d)), function(r) {
    p <- item_moments(beta[d$a[r]] - beta[d$b[r]], c(-0.7, 0.7))$P
    sample(0:2, 1, prob = p)
  }, 0L)
  f <- btl(d, "a", "b", response = "grade")
  expect_equal(nrow(f$thresholds), 2L)
  expect_equal(f$thresholds$tau[1], -f$thresholds$tau[2])
  expect_equal(f$m, 2L)
  # ordered-factor input maps by level order and keeps the labels; a plain
  # factor is refused (its level order may be alphabetical accident)
  d$lab <- factor(c("worse", "tie", "better")[d$grade + 1],
                  levels = c("worse", "tie", "better"), ordered = TRUE)
  d$plain <- factor(c("worse", "tie", "better")[d$grade + 1])
  expect_error(btl(d, "a", "b", response = "plain"), "ORDERED")
  fl <- btl(d, "a", "b", response = "lab")
  expect_equal(fl$objects$location, f$objects$location, tolerance = 1e-8)
  expect_identical(fl$categories, c("worse", "tie", "better"))
  # unused category errors informatively
  d$bad <- ifelse(d$grade == 2, 3L, 0L)
  expect_error(btl(d, "a", "b", response = "bad"), "never used")
  # category curves plot renders
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_btl_categories(f))
})

test_that("the object characteristic curve renders and the fit keeps its comparisons", {
  set.seed(7)
  beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  f <- btl(d, "a", "b", winner = "win")
  expect_true(is.data.frame(f$comparisons))
  expect_identical(names(f$comparisons),
                   c("object_a", "object_b", "response", "weight", "judge"))
  d$grade <- vapply(seq_len(nrow(d)), function(r) {
    p <- item_moments(beta[d$a[r]] - beta[d$b[r]], c(-1, 0, 1))$P
    sample(0:3, 1, prob = p)
  }, 0L)
  fg <- btl(d, "a", "b", response = "grade")
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_btl_icc(f, "C"))
  expect_no_error(plot_btl_icc(fg, "C"))
  expect_error(plot_btl_icc(f, "Z"), "no such object")
})

test_that("winner + margin entry and PC thresholds behave as designed", {
  set.seed(3)
  beta <- c(A = -1.2, B = -0.5, C = 0, D = 0.6, E = 1.1)
  beta <- beta - mean(beta)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 60), b = rep(pr[, 2], each = 60))
  d$grade <- vapply(seq_len(nrow(d)), function(r) {
    p <- item_moments(beta[d$a[r]] - beta[d$b[r]], c(-1.4, 0, 1.4))$P
    sample(0:3, 1, prob = p)
  }, 0L)
  # winner + margin is exactly the single-column coding, orientation-free
  d$winner <- ifelse(d$grade >= 2, d$a, d$b)
  d$margin <- factor(c("much", "a little", "a little", "much")[d$grade + 1],
                     levels = c("a little", "much"), ordered = TRUE)
  f1 <- btl(d, "a", "b", response = "grade")
  f2 <- btl(d, "a", "b", winner = "winner", margin = "margin")
  expect_equal(f2$objects$location, f1$objects$location, tolerance = 1e-10)
  expect_equal(f2$thresholds$tau, f1$thresholds$tau, tolerance = 1e-10)
  expect_identical(f2$categories,
                   c("worse by much", "worse by a little",
                     "better by a little", "better by much"))
  # margin without winner is an error
  expect_error(btl(d, "a", "b", margin = "margin"), "winner")
  # ties in the winner column form the middle category (5 categories)
  set.seed(11)
  d$winner2 <- ifelse(runif(nrow(d)) < 0.15, "tie", d$winner)
  f3 <- btl(d, "a", "b", winner = "winner2", margin = "margin")
  expect_equal(f3$m, 4L)
  expect_identical(f3$categories[3], "tie")
  # components: spread + kurtosis for four thresholds, skewness nowhere
  expect_setequal(f3$components$component, c("spread", "kurtosis"))
  # PC thresholds are exactly linear in the threshold index, and their
  # spread agrees with the free-mode spread component
  f4 <- btl(d, "a", "b", winner = "winner2", margin = "margin",
            thresholds = "pc")
  tv <- f4$thresholds$tau
  k <- seq_along(tv) - (length(tv) + 1) / 2
  expect_lt(max(abs(tv - sum(tv * k) / sum(k^2) * k)), 1e-10)
  expect_lt(abs(f4$components$estimate[1] - f3$components$estimate[1]), 0.25)
  expect_equal(f4$thr_structure, "pc")
})

test_that("the identifiability guards distinguish interior from extreme sparseness", {
  set.seed(5)
  beta <- c(A = -0.8, B = 0, C = 0.8)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 100), b = rep(pr[, 2], each = 100))
  d$grade <- vapply(seq_len(nrow(d)), function(r) {
    p <- item_moments(beta[d$a[r]] - beta[d$b[r]], c(-1.5, -0.5, 0.5, 1.5))$P
    sample(0:4, 1, prob = p)
  }, 0L)
  # interior categories emptied: free errors toward pc; pc fits with a note
  d$g2 <- ifelse(d$grade %in% c(1L, 3L), 2L, d$grade)
  expect_error(btl(d, "a", "b", response = "g2"), "thresholds = 'pc'")
  f <- btl(d, "a", "b", response = "g2", thresholds = "pc")
  expect_true(f$converged)
  expect_true(any(grepl("pooled", f$notes)))
  # empty extremes have no finite estimate under either structure: seven
  # declared levels with only the middle five used leaves categories 0 and
  # 6 empty in both orientations
  d$g <- factor(paste0("L", d$grade + 1),
                levels = paste0("L", 0:6), ordered = TRUE)
  expect_error(btl(d, "a", "b", response = "g"), "extreme category")
  expect_error(btl(d, "a", "b", response = "g", thresholds = "pc"),
               "extreme category")
})

test_that("exposure and carry-over are recovered and null when absent", {
  set.seed(21)
  K <- 20; objs <- sprintf("S%02d", 1:K)
  beta <- setNames(rnorm(K, 0, 0.9), objs); beta <- beta - mean(beta)
  jids <- sprintf("J%02d", 1:30)
  sim <- function(phi, psi) {
    rows <- list()
    for (j in jids) {
      pool <- sample(objs, 12)
      hc <- setNames(numeric(K), objs); hs <- hc
      for (t in 1:14) {
        p2 <- sample(pool, 2); a <- p2[1]; b <- p2[2]
        Fa <- as.numeric(hc[a] > 0); Fb <- as.numeric(hc[b] > 0)
        Wa <- if (hc[a] > 0) hs[a] / hc[a] else 0
        Wb <- if (hc[b] > 0) hs[b] / hc[b] else 0
        y <- rbinom(1, 1, plogis(beta[a] - beta[b] + phi * (Fa - Fb) +
                                   psi * (Wa - Wb)))
        rows[[length(rows) + 1L]] <- data.frame(
          a = a, b = b, judge = j, t = t, winner = if (y == 1) a else b)
        hc[a] <- hc[a] + 1; hc[b] <- hc[b] + 1
        hs[a] <- hs[a] + (2 * y - 1); hs[b] <- hs[b] - (2 * y - 1)
      }
    }
    do.call(rbind, rows)
  }
  f1 <- btl(sim(0.8, 1.0), "a", "b", winner = "winner", judge = "judge",
            order = "t")
  dp <- f1$dependence
  expect_identical(dp$effect, c("exposure", "carry_over"))
  expect_true(all(dp$p < 0.05))
  expect_true(all(dp$estimate > 0))
  f0 <- btl(sim(0, 0), "a", "b", winner = "winner", judge = "judge",
            order = "t")
  expect_true(all(f0$dependence$p > 0.05))
  # order without judge, and with half-ties, are refused
  expect_error(btl(sim(0, 0), "a", "b", winner = "winner", order = "t"),
               "judge")
})

test_that("btl_dif finds a planted judge-group effect on the right object only", {
  set.seed(2)
  K <- 12; objs <- sprintf("S%02d", 1:K)
  beta <- setNames(seq(-1.4, 1.4, length.out = K), objs)
  jids <- sprintf("J%02d", 1:28)
  grp <- setNames(rep(c("g1", "g2"), each = 14), jids)
  pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 28), b = rep(pr[, 2], each = 28))
  d$judge <- sample(jids, nrow(d), TRUE)
  shift <- ifelse(grp[d$judge] == "g2" & d$a == "S06", 1,
           ifelse(grp[d$judge] == "g2" & d$b == "S06", -1, 0))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b] + shift),
                  d$a, d$b)
  f <- btl(d, "a", "b", winner = "win", judge = "judge")
  dif <- btl_dif(f, grp)
  expect_s3_class(dif, "rasch_btl_dif")
  tested <- dif$terms$term != "band" & is.finite(dif$terms$p)
  expect_equal(dif$terms$p_adj[tested],
               p.adjust(dif$terms$p[tested], method = "holm"))
  # summary route: a single factor gives one "group" term, only S06 flagged
  expect_true(dif$summary$uniform_DIF[dif$summary$object == "S06"])
  expect_equal(sum(dif$summary$uniform_DIF), 1L)
  # magnitude route: right size, right object, nothing else
  s6 <- dif$sizes[dif$sizes$object == "S06", ]
  expect_identical(dif$size_family_n, nrow(dif$sizes))
  expect_true(s6$significant && s6$practical)
  expect_lt(abs(abs(s6$difference) - 1), 3 * s6$se)
  expect_equal(sum(dif$sizes$significant), 1L)
  expect_output(print(dif), "Resolved locations")
  # grouped characteristic curve renders
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_btl_icc(f, "S06", group = grp))

  old_graded <- rasch:::.btl_graded
  invalid_covariance <- testthat::with_mocked_bindings(
    btl_dif(f, grp, objects = "S06"),
    .btl_graded = function(...) {
      z <- old_graded(...)
      copy <- grep("^S06 \\(", rownames(z$cov_beta))[1]
      z$cov_beta[copy, copy] <- -1e6
      z
    },
    .package = "rasch")
  expect_true(all(is.finite(invalid_covariance$levels$location)))
  expect_true(all(is.na(invalid_covariance$levels$se)))
  expect_true(all(is.finite(invalid_covariance$sizes$difference)))
  expect_true(all(is.na(invalid_covariance$sizes$se)))
  expect_true(all(is.na(invalid_covariance$sizes$p_adj)))
  expect_true(any(grepl("not positive semidefinite",
                        invalid_covariance$notes)))

  no_supported_cells <- btl_dif(f, grp, objects = "S06", min_n = 1000L)
  expect_null(no_supported_cells$sizes)
  expect_identical(no_supported_cells$size_family_n, 1L)
  expect_match(paste(no_supported_cells$notes, collapse = " "),
               "unavailable follow-up question(s)", fixed = TRUE)

  testthat::local_mocked_bindings(
    .btl_graded = function(...) {
      z <- old_graded(...)
      z$converged <- FALSE
      z
    },
    .package = "rasch")
  failed_resolution <- btl_dif(f, grp)
  expect_null(failed_resolution$sizes)
  expect_gt(failed_resolution$size_family_n, 0L)
  expect_true(any(grepl("resolved calibration did not converge",
                        failed_resolution$notes)))
})

test_that("btl_dif withholds pairwise inference for thin judge-factor levels", {
  set.seed(76001)
  K <- 6L; J <- 20L
  objs <- sprintf("O%d", seq_len(K))
  jids <- sprintf("J%d", seq_len(J))
  grp <- setNames(rep(sprintf("g%d", 1:5), each = 4), jids)
  pairs <- t(utils::combn(objs, 2))
  judge_shift <- matrix(rnorm(J * K, 0, 0.8), J, K,
                        dimnames = list(jids, objs))
  rows <- lapply(jids, function(j) {
    d <- data.frame(a = rep(pairs[, 1], each = 3),
                    b = rep(pairs[, 2], each = 3), judge = j)
    beta <- as.numeric(scale(seq_len(K)))
    planted <- if (grp[j] == "g1") 2.5 else 0
    eta <- beta[match(d$a, objs)] - beta[match(d$b, objs)] +
      judge_shift[j, d$a] - judge_shift[j, d$b] +
      planted * ((d$a == "O3") - (d$b == "O3"))
    d$winner <- ifelse(runif(nrow(d)) < plogis(eta), d$a, d$b)
    d
  })
  fit <- btl(do.call(rbind, rows), "a", "b", winner = "winner",
             judge = "judge")
  out <- btl_dif(fit, grp, objects = "O3")

  expect_true(is.na(out$summary$p_uniform))
  expect_false(out$summary$uniform_DIF)
  expect_false(out$terms$inference_available[out$terms$term == "group"])
  expect_null(out$levels)
  expect_null(out$sizes)
  expect_match(paste(out$notes, collapse = " "),
               "below eight judges or eight effective judges")
})

test_that("btl_dif uses effective rather than raw judges for concentrated levels", {
  set.seed(76002)
  K <- 6L; J <- 20L
  objs <- sprintf("O%d", seq_len(K))
  jids <- sprintf("J%d", seq_len(J))
  grp <- setNames(rep(c("g1", "g2"), each = 10), jids)
  workload <- setNames(c(5L, rep(1L, 19L)), jids)
  pairs <- t(utils::combn(objs, 2))
  judge_shift <- matrix(rnorm(J * K, 0, 0.8), J, K,
                        dimnames = list(jids, objs))
  rows <- lapply(jids, function(j) {
    d <- data.frame(a = rep(pairs[, 1], each = 3 * workload[j]),
                    b = rep(pairs[, 2], each = 3 * workload[j]),
                    judge = j)
    beta <- as.numeric(scale(seq_len(K)))
    planted <- if (grp[j] == "g1") 2.5 else 0
    eta <- beta[match(d$a, objs)] - beta[match(d$b, objs)] +
      judge_shift[j, d$a] - judge_shift[j, d$b] +
      planted * ((d$a == "O3") - (d$b == "O3"))
    d$winner <- ifelse(runif(nrow(d)) < plogis(eta), d$a, d$b)
    d
  })
  fit <- btl(do.call(rbind, rows), "a", "b", winner = "winner",
             judge = "judge")
  expect_true(fit$cl$inference_available)
  out <- btl_dif(fit, grp, objects = "O3")

  expect_true(is.na(out$summary$p_uniform))
  expect_false(out$summary$uniform_DIF)
  tr <- out$terms[out$terms$term == "group", ]
  expect_equal(tr$min_judges, 10)
  expect_lt(tr$min_effective_judges, 8)
  expect_false(tr$inference_available)
  expect_null(out$levels)
  expect_null(out$sizes)
  expect_match(paste(out$notes, collapse = " "), "eight effective judges")
})

test_that("btl_dif fits several judge factors jointly (main and factorial)", {
  set.seed(2)
  K <- 12; objs <- sprintf("S%02d", 1:K)
  beta <- setNames(seq(-1.4, 1.4, length.out = K), objs)
  jids <- sprintf("J%02d", 1:28)
  A <- setNames(rep(c("g1", "g2"), each = 14), jids)   # real DIF on S06
  B <- setNames(rep(c("h1", "h2"), 14), jids)          # null second factor
  pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 28), b = rep(pr[, 2], each = 28))
  d$judge <- sample(jids, nrow(d), TRUE)
  shift <- ifelse(A[d$judge] == "g2" & d$a == "S06", 2,
           ifelse(A[d$judge] == "g2" & d$b == "S06", -2, 0))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b] + shift),
                  d$a, d$b)
  f <- btl(d, "a", "b", winner = "win", judge = "judge")

  # main effects: both factors modelled together, one term each
  mn <- btl_dif(f, list(A = A, B = B), effects = "main")
  expect_setequal(unique(mn$summary$term), c("A", "B"))
  # DIF is found on factor A (object S06) and not on the null factor B
  expect_true(mn$summary$uniform_DIF[mn$summary$object == "S06" &
                                     mn$summary$term == "A"])
  expect_equal(sum(mn$summary$uniform_DIF[mn$summary$term == "B"]), 0L)
  expect_equal(sum(mn$summary$uniform_DIF), 1L)

  # factorial adds the factor-by-factor interaction term
  fac <- btl_dif(f, list(A = A, B = B), effects = "factorial")
  expect_true("A:B" %in% fac$summary$term)
})

test_that("fit_summary_table dispatches for paired-comparison fits", {
  set.seed(1)
  beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30),
                  judge = sample(sprintf("J%d", 1:8), 180, TRUE))
  d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  bt <- btl(d, "a", "b", winner = "win", judge = "judge", order = "t")
  ft <- fit_summary_table(bt)
  expect_identical(names(ft), c("statistic", "value"))
  expect_equal(ft$value[ft$statistic == "Objects"], "4")
  expect_equal(ft$value[ft$statistic == "Standard errors"],
               "sandwich, clustered by judge")
  expect_identical(ft$value[ft$statistic == "Pairwise fit probability"],
                   "unavailable (judge clustering)")
  printed_fit <- capture.output(print(bt))
  expect_true(any(grepl("probability unavailable", printed_fit,
                        fixed = TRUE)))
  expect_false(any(grepl("p = $", printed_fit)))
  expect_true(any(grepl("Within-judge exposure", ft$statistic)))
  expect_true(any(grepl("Within-judge carry-over", ft$statistic)))
  expect_true(all(grepl("Holm p", ft$value[
    grepl("Within-judge", ft$statistic)])))
  # An older or incomplete fit must not silently substitute an unadjusted
  # probability on a surface whose decisions use the Holm family.
  no_adjustment <- bt
  no_adjustment$dependence$p_adj <- NULL
  ft_unavailable <- fit_summary_table(no_adjustment)
  dependence_rows <- grepl("Within-judge", ft_unavailable$statistic)
  expect_true(all(grepl("Holm p unavailable",
                        ft_unavailable$value[dependence_rows], fixed = TRUE)))
  expect_false(any(grepl(" p = ", ft_unavailable$value[dependence_rows],
                         fixed = TRUE)))
  printed_unavailable <- capture.output(print(no_adjustment))
  printed_dependence <- printed_unavailable[
    grepl("Within-judge|First-position", printed_unavailable)]
  expect_true(any(grepl("Holm p unavailable", printed_dependence,
                        fixed = TRUE)))
  expect_false(any(grepl(" p = ", printed_dependence, fixed = TRUE)))
  expect_no_error(plot_btl_dependence(no_adjustment, "exposure"))
  # graded fit reports its category count and threshold structure
  d$grade <- vapply(seq_len(nrow(d)), function(r) {
    p <- item_moments(beta[d$a[r]] - beta[d$b[r]], c(-1, 0, 1))$P
    sample(0:3, 1, prob = p)
  }, 0L)
  fg <- fit_summary_table(btl(d, "a", "b", response = "grade"))
  expect_true(any(grepl("Graded paired comparisons \\(4 categories\\)",
                        fg$value)))
  expect_equal(fg$value[fg$statistic == "Threshold structure"],
               "free symmetric")
})

test_that("objects and judges carry infit and outfit mean squares", {
  set.seed(31)
  K <- 8; objs <- sprintf("E%02d", 1:K)
  beta <- setNames(seq(-1.4, 1.4, length.out = K), objs)
  pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 24), b = rep(pr[, 2], each = 24))
  d$judge <- sample(sprintf("J%02d", 1:6), nrow(d), TRUE)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  f <- btl(d, "a", "b", winner = "win", judge = "judge")
  expect_true(all(c("infit_ms", "outfit_ms") %in% names(f$objects)))
  expect_true(all(c("infit_ms", "outfit_ms") %in% names(f$judges)))
  # well-fitting data: mean squares near 1, finite and positive
  expect_true(all(f$objects$infit_ms > 0 & is.finite(f$objects$infit_ms)))
  expect_lt(abs(mean(f$objects$infit_ms) - 1), 0.2)
  expect_lt(abs(mean(f$objects$outfit_ms) - 1), 0.2)
})

test_that("plot_btl_icc omits opponents met too few times", {
  set.seed(32)
  objs <- sprintf("O%d", 1:6)
  beta <- setNames(seq(-1.2, 1.2, length.out = 6), objs)
  pr <- t(combn(objs, 2))
  # O1 vs O2 dense (30), every other pair sparse (2)
  reps <- ifelse(pr[, 1] == "O1" & pr[, 2] == "O2", 30, 2)
  d <- do.call(rbind, Map(function(a, b, n)
    data.frame(a = a, b = b, k = seq_len(n)), pr[, 1], pr[, 2], reps))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  f <- btl(d, "a", "b", winner = "win")
  pdf(NULL); on.exit(dev.off())
  shown <- plot_btl_icc(f, "O1", min_n = 10)
  # only O2 (met 30 times) survives; the 2-comparison opponents are dropped
  expect_setequal(shown, "O2")
  # a lower threshold keeps them all
  expect_gt(length(plot_btl_icc(f, "O1", min_n = 1)), 1L)
})

test_that("btl_dif resolves interactions by cells and supersedes lower terms", {
  set.seed(1)
  K <- 10; objs <- sprintf("S%02d", 1:K)
  beta <- setNames(seq(-1.0, 1.0, length.out = K), objs)
  jids <- sprintf("J%02d", 1:40)
  A <- setNames(rep(c("g1", "g2"), each = 20), jids)
  B <- setNames(rep(rep(c("h1", "h2"), each = 10), 2), jids)  # A, B crossed
  pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
  d$judge <- sample(jids, nrow(d), TRUE)
  # DIF on S05 concentrated in the single g2:h2 cell: a main effect AND an
  # interaction both surface, so the A main effect is superseded by A:B
  cell <- paste(A[d$judge], B[d$judge])
  sh <- ifelse(cell == "g2 h2" & d$a == "S05", 2.5,
        ifelse(cell == "g2 h2" & d$b == "S05", -2.5, 0))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b] + sh), d$a, d$b)
  f <- btl(d, "a", "b", winner = "win", judge = "judge")
  r <- btl_dif(f, list(A = A, B = B), effects = "factorial")

  s5 <- r$summary[r$summary$object == "S05", ]
  # the interaction is flagged and supersedes the significant A main effect
  expect_true(s5$uniform_DIF[s5$term == "A:B"])
  expect_true(s5$superseded[s5$term == "A"])
  # only the non-superseded interaction is resolved. The cell locations stay
  # visible, while its magnitude is the difference-in-differences contrast.
  s5sz <- r$sizes[r$sizes$object == "S05" & r$sizes$term == "A:B", ]
  expect_equal(nrow(s5sz), 1L)
  expect_match(s5sz$contrast, "g2 - g1 x h2 - h1", fixed = TRUE)
  expect_setequal(r$levels$level[r$levels$object == "S05" &
                                  r$levels$term == "A:B"],
                  c("g1:h1", "g1:h2", "g2:h1", "g2:h2"))
  expect_gt(abs(s5sz$difference), 2)
  expect_true(s5sz$significant)
})

test_that("btl_dif tolerates adversarial factor names (band, f1)", {
  set.seed(2)
  K <- 12; objs <- sprintf("S%02d", 1:K)
  beta <- setNames(seq(-1.4, 1.4, length.out = K), objs)
  jids <- sprintf("J%02d", 1:28)
  g <- setNames(rep(c("lo", "hi"), each = 14), jids)   # real DIF on S06
  h <- setNames(rep(c("p", "q"), 14), jids)            # null second factor
  pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 28), b = rep(pr[, 2], each = 28))
  d$judge <- sample(jids, nrow(d), TRUE)
  sh <- ifelse(g[d$judge] == "hi" & d$a == "S06", 2,
        ifelse(g[d$judge] == "hi" & d$b == "S06", -2, 0))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b] + sh), d$a, d$b)
  f <- btl(d, "a", "b", winner = "win", judge = "judge")
  # a factor named "band" must not collide with the opponent-band variable
  rb <- btl_dif(f, list(band = g))
  expect_identical(unique(rb$summary$term), "band")
  expect_true(rb$summary$uniform_DIF[rb$summary$object == "S06"])
  # a factor named like a stand-in ("f1") must not corrupt the other labels
  rf <- btl_dif(f, list(x = g, f1 = h), effects = "main")
  expect_setequal(unique(rf$summary$term), c("x", "f1"))
  expect_true(rf$summary$uniform_DIF[rf$summary$object == "S06" &
                                     rf$summary$term == "x"])
  expect_equal(sum(rf$summary$uniform_DIF[rf$summary$term == "f1"]), 0L)

  sparse <- setNames(ifelse(jids %in% jids[1:5], "rare", "common"), jids)
  rs <- btl_dif(f, list(cohort = sparse), objects = "S06")
  note <- paste(rs$notes, collapse = " ")
  expect_match(note, "cohort")
  expect_false(grepl("term\\(s\\) f1", note))
})

test_that("unavailable object DIF terms remain in the Holm family", {
  set.seed(901)
  objects <- LETTERS[1:5]
  judges <- sprintf("J%02d", 1:20)
  group <- setNames(rep(c("g1", "g2"), each = 10), judges)
  beta <- setNames(seq(-1, 1, length.out = 5), objects)
  pairs <- t(combn(objects, 2))
  rows <- list()
  for (i in seq_len(nrow(pairs))) {
    # A is observed only in g1. The other objects retain both groups, so the
    # requested family is usable but A's two questions are not estimable.
    js <- if ("A" %in% pairs[i, ]) judges[1:10] else judges
    for (r in 1:3) for (j in js)
      rows[[length(rows) + 1L]] <- data.frame(
        a = pairs[i, 1], b = pairs[i, 2], judge = j)
  }
  d <- do.call(rbind, rows)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]),
                  d$a, d$b)
  fit <- btl(d, "a", "b", "win", judge = "judge")
  out <- btl_dif(fit, list(group = group), min_n = 2)

  unavailable <- out$terms$object == "A" &
    out$terms$term %in% c("group", "group:band")
  expect_equal(sum(unavailable), 2L)
  expect_true(all(is.na(out$terms$p[unavailable])))
  family <- out$terms$term != "band"
  expect_equal(sum(family), 2L * length(objects))
  usable <- family & is.finite(out$terms$p)
  expect_equal(out$terms$p_adj[usable],
               p.adjust(out$terms$p[usable], "holm", n = sum(family)))
  expect_match(paste(out$notes, collapse = " "),
               "remain in the adjusted-probability family")
})

test_that("btl stores per-comparison dependence covariates and plots them", {
  set.seed(1)
  beta <- c(A = -0.8, B = -0.2, C = 0.4, D = 0.9)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
  d$judge <- sample(sprintf("J%02d", 1:8), nrow(d), TRUE)
  d <- d[order(d$judge), ]
  d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  f <- btl(d, "a", "b", winner = "win", judge = "judge", order = "t")

  dd <- f$dependence_data
  expect_false(is.null(dd))
  expect_setequal(names(dd), c("judge", "order", "object_a", "object_b",
                               "response", "weight", "exposure", "carry_over"))
  expect_equal(nrow(dd), f$n_comparisons)
  # the informative-comparison counts equal the non-zero covariate counts
  ni <- setNames(f$dependence$n_informative, f$dependence$effect)
  expect_equal(sum(dd$exposure != 0), ni[["exposure"]])
  expect_equal(sum(dd$carry_over != 0), ni[["carry_over"]])
  # a fit without an order column carries no dependence data
  expect_null(btl(d, "a", "b", winner = "win")$dependence_data)

  # the display runs and returns the binned observed and fitted departures
  pdf(NULL); on.exit(dev.off())
  be <- plot_btl_dependence(f, "exposure")
  expect_true(all(c("covariate", "observed", "fitted", "n") %in% names(be)))
  expect_setequal(be$covariate, c(-1, 0, 1))          # exposure's three levels
  expect_equal(sum(be$n), f$n_comparisons)
  expect_equal(plot_btl_dependence(f), be)
  expect_no_error(plot_btl_dependence(f, "carry_over"))
  expect_error(plot_btl_dependence(f, c("exposure", "carry_over")),
               "effect.*one of")
  expect_error(plot_btl_dependence(btl(d, "a", "b", winner = "win"),
                                   "exposure"), "no dependence data")
})

test_that("graded free-threshold fits with dependence estimate correctly (C1)", {
  # the threshold x dependence Hessian block: 5 categories (q = 2 free
  # thresholds) with both dependence effects once diverged from a transposed
  # block assignment; the fit must converge and recover the planted effects
  set.seed(11)
  objs <- sprintf("O%d", 1:6); beta <- setNames(seq(-1, 1, length.out = 6), objs)
  tau <- c(-1.6, -0.6, 0.6, 1.6)
  pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
  d$judge <- sample(sprintf("J%d", 1:10), nrow(d), TRUE)
  d <- d[order(d$judge), ]; d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  seen <- new.env(parent = emptyenv()); hs <- new.env(parent = emptyenv())
  hc <- new.env(parent = emptyenv())
  g0 <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
  resp <- integer(nrow(d))
  for (r in seq_len(nrow(d))) {
    j <- d$judge[r]; a <- d$a[r]; b <- d$b[r]
    ka <- paste(j, a); kb <- paste(j, b)
    Fa <- as.numeric(g0(seen, ka) > 0); Fb <- as.numeric(g0(seen, kb) > 0)
    Wa <- if (g0(hc, ka) > 0) g0(hs, ka) / g0(hc, ka) else 0
    Wb <- if (g0(hc, kb) > 0) g0(hs, kb) / g0(hc, kb) else 0
    dd <- beta[[a]] - beta[[b]] + 0.4 * (Fa - Fb) + 0.8 * (Wa - Wb)
    x <- sample(0:4, 1, prob = item_moments(dd, tau)$P); resp[r] <- x
    assign(ka, g0(seen, ka) + 1, seen); assign(kb, g0(seen, kb) + 1, seen)
    assign(ka, g0(hc, ka) + 1, hc); assign(kb, g0(hc, kb) + 1, hc)
    assign(ka, g0(hs, ka) + (2 * x / 4 - 1), hs)
    assign(kb, g0(hs, kb) + (2 * (4 - x) / 4 - 1), hs)
  }
  d$resp <- resp
  f <- btl(d, "a", "b", response = "resp", judge = "judge", order = "t")
  expect_true(f$converged)
  expect_lt(f$iterations, 20)
  dep <- setNames(f$dependence$estimate, f$dependence$effect)
  expect_lt(abs(dep[["exposure"]] - 0.4), 3 * f$dependence$se[1])
  expect_lt(abs(dep[["carry_over"]] - 0.8), 3 * f$dependence$se[2])
  expect_true(is.na(f$dependence$p[f$dependence$effect == "carry_over"]))
  use <- is.finite(f$dependence$p)
  expect_equal(f$dependence$p_adj[use],
               p.adjust(f$dependence$p[use], "holm",
                        n = nrow(f$dependence)))
  expect_match(paste(f$notes, collapse = " "), "fewer than 30 judges")
  # thresholds recovered too (the corrupted block used to distort them)
  expect_lt(max(abs(f$thresholds$tau - tau)), 0.25)
})

test_that("btl_dif holds fitted dependence effects fixed (H1)", {
  # strong planted carry-over, NO DIF: the dependence-adjusted moments must
  # not let a judge factor absorb the sequential structure as spurious DIF
  set.seed(107)
  objs <- sprintf("O%d", 1:6); beta <- setNames(seq(-1, 1, length.out = 6), objs)
  jids <- sprintf("J%d", 1:10); grp <- setNames(rep(c("x", "y"), each = 5), jids)
  pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 24), b = rep(pr[, 2], each = 24))
  d$judge <- sample(jids, nrow(d), TRUE)
  d <- d[order(d$judge), ]; d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  hs <- new.env(parent = emptyenv()); hc <- new.env(parent = emptyenv())
  g0 <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
  win <- character(nrow(d))
  for (r in seq_len(nrow(d))) {
    j <- d$judge[r]; a <- d$a[r]; b <- d$b[r]
    ka <- paste(j, a); kb <- paste(j, b)
    Wa <- if (g0(hc, ka) > 0) g0(hs, ka) / g0(hc, ka) else 0
    Wb <- if (g0(hc, kb) > 0) g0(hs, kb) / g0(hc, kb) else 0
    y <- as.integer(runif(1) < plogis(beta[[a]] - beta[[b]] + 1.5 * (Wa - Wb)))
    win[r] <- if (y == 1) a else b
    assign(ka, g0(hc, ka) + 1, hc); assign(kb, g0(hc, kb) + 1, hc)
    assign(ka, g0(hs, ka) + (2 * y - 1), hs)
    assign(kb, g0(hs, kb) - (2 * y - 1), hs)
  }
  d$win <- win
  f <- btl(d, "a", "b", winner = "win", judge = "judge", order = "t")
  # the fit stores the row-aligned covariates for downstream analyses
  expect_true(all(c("exposure", "carry_over") %in% names(f$comparisons)))
  adj <- btl_dif(f, grp)
  expect_equal(sum(adj$summary$uniform_DIF), 0L)
  # judge-level aggregation makes even the unadjusted screen robust to
  # within-judge dependence here (the dependence lives inside judges, and
  # judges are the analysis units): both screens stay clean. The fitted
  # covariates still matter for the resolved magnitudes, which refit
  # locations with the dependence effects held fixed.
  fU <- f; fU$comparisons$exposure <- NULL; fU$comparisons$carry_over <- NULL
  unadj <- btl_dif(fU, grp)
  expect_equal(sum(unadj$summary$uniform_DIF), 0L)
})

test_that("btl_dif weights count-aggregated comparisons correctly (H2)", {
  set.seed(5)
  objs <- sprintf("S%d", 1:6); beta <- setNames(seq(-1, 1, length.out = 6), objs)
  jids <- sprintf("J%d", 1:16); grp <- setNames(rep(c("g1", "g2"), each = 8), jids)
  pr <- t(combn(objs, 2)); rows <- list()
  for (i in seq_len(nrow(pr))) for (j in jids) {
    sh <- ifelse(grp[j] == "g2" & pr[i, 1] == "S3", 1.8,
          ifelse(grp[j] == "g2" & pr[i, 2] == "S3", -1.8, 0))
    wins <- rbinom(1, 5, plogis(beta[pr[i, 1]] - beta[pr[i, 2]] + sh))
    rows[[length(rows) + 1]] <- data.frame(a = pr[i, 1], b = pr[i, 2],
                                           judge = j, win = pr[i, 1], k = wins)
    rows[[length(rows) + 1]] <- data.frame(a = pr[i, 1], b = pr[i, 2],
                                           judge = j, win = pr[i, 2], k = 5 - wins)
  }
  agg <- do.call(rbind, rows); agg <- agg[agg$k > 0, ]
  expd <- agg[rep(seq_len(nrow(agg)), agg$k), ]; expd$k <- 1
  ra <- btl_dif(btl(agg, "a", "b", winner = "win", judge = "judge", count = "k"), grp)
  re <- btl_dif(btl(expd, "a", "b", winner = "win", judge = "judge"), grp)
  # judge-cell means and weighted-quantile bands are invariant to row
  # expansion, so the aggregated analysis IS the expanded one exactly
  Fa <- ra$summary$F_uniform[order(ra$summary$object)]
  Fe <- re$summary$F_uniform[order(re$summary$object)]
  expect_lt(max(abs(Fa - Fe) / pmax(Fe, 1)), 1e-8)
  expect_true(ra$summary$uniform_DIF[ra$summary$object == "S3"])
  expect_equal(ra$summary$uniform_DIF, re$summary$uniform_DIF)
})

test_that("NA judges are dropped with a note, not a crash (M2)", {
  set.seed(1)
  beta <- c(A = -1, B = 0, C = 1); pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
  d$judge <- sample(c("J1", "J2", "J3", NA), nrow(d), TRUE)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  f <- btl(d, "a", "b", winner = "win", judge = "judge")
  expect_s3_class(f, "rasch_btl")
  expect_true(any(grepl("dropped", f$notes)))
  expect_false(anyNA(f$comparisons$judge))
})

test_that("a separated dependence effect is dropped with a note", {
  # few informative comparisons all pointing one way: the coefficient would
  # run to a boundary with a maxit-dependent 'estimate'; it must be set
  # aside, not reported as decisively significant
  set.seed(3)
  objs <- c("A", "B", "C"); beta <- setNames(c(-1, 0, 1), objs)
  pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
  d$judge <- sample(sprintf("J%d", 1:6), nrow(d), TRUE)
  d <- d[order(d$judge), ]; d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  f <- btl(d, "a", "b", winner = "win", judge = "judge", order = "t")
  expect_true(any(grepl("separated", f$notes)))
  if (!is.null(f$dependence))
    expect_true(all(abs(f$dependence$estimate) < 10))
})

test_that("btl_transitivity and btl_dimensionality read one-D vs a swirl", {
  mk <- function(cyc, seed) {
    set.seed(seed); K <- 8; objs <- sprintf("O%d", 1:K)
    beta <- setNames(seq(-1.5, 1.5, length.out = K), objs)
    th <- setNames(2 * pi * (0:(K - 1)) / K, objs)
    pr <- t(combn(objs, 2))
    d <- data.frame(a = rep(pr[, 1], each = 25), b = rep(pr[, 2], each = 25))
    d$judge <- sample(sprintf("J%d", 1:10), nrow(d), TRUE)
    lp <- beta[d$a] - beta[d$b] + cyc * sin(th[d$a] - th[d$b])
    d$win <- ifelse(runif(nrow(d)) < plogis(lp), d$a, d$b)
    btl(d, "a", "b", "win", judge = "judge")
  }
  # one-dimensional: consistent, leading bimension within the noise band
  f1 <- mk(0, 2)
  t1 <- btl_transitivity(f1)
  d1 <- btl_dimensionality(f1, reps = 40, independent_comparisons = TRUE)
  expect_lt(t1$summary$circular_rate, 0.1)
  expect_gt(t1$summary$consistency, 0.6)
  expect_false(d1$leading_structured)
  expect_true(all(abs(d1$residual_matrix + t(d1$residual_matrix)) < 1e-8)) # skew
  expect_s3_class(t1, "rasch_btl_transitivity")
  expect_s3_class(d1, "rasch_btl_dim")
  expect_false(is.null(t1$judges))    # judges present -> per-judge table
  expect_equal(d1$reference$p,
               (1 + sum(d1$reference$draws >= d1$bimensions$strength[1])) /
                 (length(d1$reference$draws) + 1))
  expect_equal(d1$reference$p_adj, d1$reference$p)
  expect_identical(d1$leading_structured, d1$reference$p_adj <= 0.05)
  current_print <- capture.output(print(d1))
  expect_true(any(grepl(sprintf("adjusted p = %.3f", d1$reference$p_adj),
                        current_print, fixed = TRUE)))
  old_d1 <- d1
  old_d1$reference$p <- old_d1$reference$p_adj <- NULL
  old_print <- expect_no_error(capture.output(print(old_d1)))
  expect_true(any(grepl("adjusted p unavailable", old_print, fixed = TRUE)))

  # a cyclic swirl: leading bimension clears the reference, most of residual
  f2 <- mk(1.6, 2)
  d2 <- btl_dimensionality(f2, reps = 40, independent_comparisons = TRUE)
  expect_true(d2$leading_structured)
  expect_gt(d2$bimensions$strength[1], d2$reference$p95)
  expect_identical(d2$leading_structured, d2$reference$p_adj <= 0.05)
  expect_gt(d2$bimensions$prop_residual[1], 0.5)

  # genuine intransitivity (flat locations, strong cycle) -> loops above chance
  set.seed(3); K <- 7; objs <- sprintf("O%d", 1:K)
  th <- setNames(2 * pi * (0:(K - 1)) / K, objs); pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
  d$win <- ifelse(runif(nrow(d)) < plogis(1.4 * sin(th[d$a] - th[d$b])), d$a, d$b)
  t3 <- btl_transitivity(btl(d, "a", "b", "win"))
  expect_gt(t3$summary$circular_rate, 0.25)   # worse than chance
  expect_lt(t3$summary$consistency, 0)
})

test_that("judge_surprise flags a judge's systematic contrary judgements", {
  set.seed(1); K <- 8; objs <- sprintf("O%d", 1:K)
  beta <- setNames(seq(-1.5, 1.5, length.out = K), objs)
  jids <- sprintf("J%d", 1:8); pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 12), b = rep(pr[, 2], each = 12))
  d$judge <- sample(jids, nrow(d), TRUE)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  # J1 is contrarian on the extremes: always backs the weakest (O1), always
  # sinks the strongest (O8)
  b <- d$judge == "J1"
  d$win[b & (d$a == "O1" | d$b == "O1")] <- "O1"
  d$win[b & d$a == "O8"] <- d$b[b & d$a == "O8"]
  d$win[b & d$b == "O8"] <- d$a[b & d$b == "O8"]
  f <- btl(d, "a", "b", "win", judge = "judge")

  js <- judge_surprise(f, "J1")
  expect_s3_class(js, "rasch_btl_judge")
  eligible <- js$objects$n >= js$min_n
  expect_equal(js$objects$p_adj[eligible],
               p.adjust(js$objects$p[eligible], "holm"))
  expect_true(all(is.na(js$objects$p_adj[!eligible])))
  s <- js$objects[js$objects$surprise, ]
  expect_true(all(c("O1", "O8") %in% s$object))          # both extremes flagged
  # correct direction: O8 strong under-rated (z<0), O1 weak over-rated (z>0)
  expect_lt(js$objects$z[js$objects$object == "O8"], 0)
  expect_gt(js$objects$z[js$objects$object == "O1"], 0)
  expect_equal(js$objects$type[js$objects$object == "O8"],
               "strong object under-rated")

  # The scale origin is arbitrary. Anchoring two objects ten logits higher is
  # a pure translation of this fit and must not change which judgements are
  # surprising or how their direction is described.
  anchors <- setNames(
    f$objects$location[match(c("O1", "O8"), f$objects$object)] + 10,
    c("O1", "O8"))
  translated <- btl(d, "a", "b", "win", judge = "judge",
                    anchors = anchors)
  js_translated <- judge_surprise(translated, "J1")
  ref <- js$objects[order(js$objects$object), ]
  got <- js_translated$objects[order(js_translated$objects$object), ]
  expect_equal(got$z, ref$z, tolerance = 1e-10)
  expect_identical(got$surprise, ref$surprise)
  expect_identical(got$type, ref$type)
  # a model-conforming judge shows no systematic surprise
  expect_equal(sum(judge_surprise(f, "J3")$objects$surprise), 0L)

  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_btl_judge_map(f, "J1"))
  expect_error(judge_surprise(f, "nobody"), "no comparisons")
})

test_that("judge_pair_surprise flags the matchups a judge got against the grain", {
  set.seed(1); K <- 8; objs <- sprintf("O%d", 1:K)
  beta <- setNames(seq(-1.5, 1.5, length.out = K), objs)
  jids <- sprintf("J%d", 1:8); pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 12), b = rep(pr[, 2], each = 12))
  d$judge <- sample(jids, nrow(d), TRUE)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  b <- d$judge == "J1"
  d$win[b & (d$a == "O1" | d$b == "O1")] <- "O1"
  d$win[b & d$a == "O8"] <- d$b[b & d$a == "O8"]
  d$win[b & d$b == "O8"] <- d$a[b & d$b == "O8"]
  f <- btl(d, "a", "b", "win", judge = "judge")

  jp <- judge_pair_surprise(f, "J1")
  expect_s3_class(jp, "rasch_btl_judge_pairs")
  eligible <- jp$pairs$n >= jp$min_n
  expect_equal(jp$pairs$p_adj[eligible],
               p.adjust(jp$pairs$p[eligible], "holm"))
  expect_true(all(is.na(jp$pairs$p_adj[!eligible])))
  s <- jp$pairs[jp$pairs$surprise, ]
  expect_gt(nrow(s), 3L)
  # a flagged matchup is one where the stronger object under-performed
  expect_true(all(s$z < 0))
  expect_true(all(s$loc_hi >= s$loc_lo))          # orientation to the stronger
  # after familywise adjustment the retained surprises involve an extreme
  # whose judgements J1 deliberately distorted
  involves <- vapply(seq_len(nrow(s)), function(i)
    any(c("O1", "O8") %in% c(s$object_hi[i], s$object_lo[i])), TRUE)
  expect_true(all(involves))
  # a model-conforming judge produces no familywise flag in this fixture
  expect_equal(sum(judge_pair_surprise(f, "J3")$pairs$surprise), 0L)

  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_btl_judge_map(f, "J1"))
  expect_error(judge_pair_surprise(f, "nobody"), "no comparisons")
})

test_that("unavailable judge residuals remain in the Holm family", {
  d <- simulate_btl(6, 6, reps_per_pair = 4, seed = 10041)
  f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  judge <- f$comparisons$judge[1L]
  real_moments <- rasch:::.btl_fitted_moments
  unavailable <- function(fit, cmp) {
    z <- real_moments(fit, cmp)
    z$E[1L] <- NA_real_
    z$V[1L] <- NA_real_
    z
  }

  js <- testthat::with_mocked_bindings(
    judge_surprise(f, judge, min_n = 1L),
    .btl_fitted_moments = unavailable, .package = "rasch")
  planned_object <- js$objects$n >= js$min_n
  usable_object <- planned_object & is.finite(js$objects$p)
  expect_true(any(planned_object & !usable_object))
  expect_equal(
    js$objects$p_adj[usable_object],
    p.adjust(js$objects$p[usable_object], "holm", n = sum(planned_object)))

  jp <- testthat::with_mocked_bindings(
    judge_pair_surprise(f, judge, min_n = 1L),
    .btl_fitted_moments = unavailable, .package = "rasch")
  planned_pair <- jp$pairs$n >= jp$min_n
  usable_pair <- planned_pair & is.finite(jp$pairs$p)
  expect_true(any(planned_pair & !usable_pair))
  expect_equal(
    jp$pairs$p_adj[usable_pair],
    p.adjust(jp$pairs$p[usable_pair], "holm", n = sum(planned_pair)))
})

test_that("btl_dimensionality is calibrated and powered on non-cyclic 2-D data", {
  # two camps of judges ranking by two ORTHOGONAL attributes -- an
  # interpretable multidimensionality the bimension method is not tuned for
  K <- 8; njudge <- 12; objs <- sprintf("O%d", 1:K)
  jids <- sprintf("J%d", 1:njudge); pr <- t(combn(objs, 2))
  set.seed(99); u0 <- rnorm(K); v0 <- resid(lm(rnorm(K) ~ u0))
  u <- as.numeric(scale(u0)) * 1.3; v <- as.numeric(scale(v0)) * 1.3
  sim <- function(seed, twoD, s, nper) {
    set.seed(seed)
    d <- data.frame(a = rep(pr[, 1], each = nper), b = rep(pr[, 2], each = nper))
    d$judge <- sample(jids, nrow(d), TRUE)
    camp <- setNames(rep(c("u", "v"), length.out = njudge), jids)
    useu <- if (twoD) camp[d$judge] == "u" else rep(TRUE, nrow(d))
    ai <- ifelse(useu, u[match(d$a, objs)], v[match(d$a, objs)])
    aj <- ifelse(useu, u[match(d$b, objs)], v[match(d$b, objs)])
    d$win <- ifelse(runif(nrow(d)) < plogis(s * (ai - aj)), d$a, d$b)
    btl(d, "a", "b", "win", judge = "judge")
  }
  # genuine 2-D structure is flagged, with the leading bimension dominant
  d2 <- btl_dimensionality(sim(1, TRUE, 2.4, 70), reps = 80,
                          independent_comparisons = TRUE)
  expect_true(d2$leading_structured)
  expect_gt(d2$bimensions$strength[1], d2$reference$p95)
  expect_gt(d2$bimensions$prop_residual[1], 0.5)
  # a single-attribute (truly 1-D) fit is not flagged, even well separated
  d1 <- btl_dimensionality(sim(7, FALSE, 1.6, 40), reps = 80,
                          independent_comparisons = TRUE)
  expect_false(d1$leading_structured)
})

test_that("ordered count rows are refused before fitting", {
  d <- simulate_btl(5, 12, reps_per_pair = 3,
                    dependence = list(exposure = 0.3, carry_over = 0.2),
                    seed = 991)
  d$count <- 1L
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
             order = "order", count = "count")
  ignored <- d
  ignored$count[1L] <- 5L
  ignored$winner[1L] <- NA_character_
  expect_no_error(btl(ignored, "object_a", "object_b", winner = "winner",
                      judge = "judge", order = "order", count = "count"))
  d$count <- 5L
  expect_error(
    btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
        order = "order", count = "count"),
    "count-compressed rows cannot be combined")

  # Keep the downstream guard for a fit saved by an earlier package version.
  # Such a fit may still carry an ordered row with a replication weight.
  fit$comparisons$weight <- 5
  fit$dependence_data$weight <- 5
  out <- btl_dimensionality(fit, reps = 20)
  expect_true(is.finite(out$bimensions$strength[1L]))
  expect_true(is.na(out$leading_structured))
  expect_true(is.na(out$reference$mean))
  expect_identical(out$reference$reps, 0L)
  expect_length(out$reference$draws, 0L)
  expect_true(any(grepl("within-row comparison sequence", out$notes,
                        fixed = TRUE)))
  printed <- capture.output(print(out))
  expect_true(any(grepl("reference unavailable", printed, fixed = TRUE)))
  expect_false(any(grepl("reference 5% upper limit: NA", printed,
                         fixed = TRUE)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_btl_scree(out))

  # Whole count rows are estimable when position is the only fitted effect:
  # no within-row response history has to be invented.
  pos <- btl(d, "object_a", "object_b", winner = "winner",
             position = TRUE, count = "count")
  pos_out <- suppressWarnings(btl_dimensionality(pos, reps = 20))
  expect_identical(pos_out$reference$reps, 20L)
  expect_true(is.finite(pos_out$reference$mean))
})

test_that("the default (unanchored, no-position) path is unchanged", {
  # regression net: the sum-zero path must be bit-identical to the pre-anchor
  # engine. Locations and SEs are hard-coded from the code as it stood before
  # the position/anchor design matrix was generalised (captured via git stash).
  set.seed(42)
  beta <- c(A = -1.2, B = -0.4, C = 0.1, D = 0.6, E = 0.9)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 50), b = rep(pr[, 2], each = 50))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  f <- btl(d, "a", "b", winner = "win")
  loc0 <- c(-1.43966223417452, -0.174047810032907, 0.0484307842564718,
            0.72306904740734, 0.842210212543618)
  se0 <- c(0.15899493095724, 0.123599149854014, 0.130315397443216,
           0.129753126063131, 0.130417229533075)
  expect_identical(f$objects$object, c("A", "B", "C", "D", "E"))
  expect_equal(f$objects$location, loc0, tolerance = 1e-10)
  expect_equal(f$objects$se, se0, tolerance = 1e-10)
  expect_equal(sum(f$objects$location), 0, tolerance = 1e-12)  # sum-zero
  expect_null(f$anchors)
  expect_null(f$dependence)
})

test_that("position bias: a first-position advantage is recovered", {
  # planted +0.5 first-position advantage with randomised orientation
  # (object_a is the first-presented object of every comparison)
  set.seed(101)
  beta <- c(A = -1.0, B = -0.3, C = 0.2, D = 0.5, E = 0.9); beta <- beta - mean(beta)
  objs <- names(beta); pr <- t(combn(objs, 2)); n_per <- 80
  sim_pos <- function(gamma) {
    rows <- list()
    for (i in seq_len(nrow(pr))) for (k in seq_len(n_per)) {
      # randomise which object is presented first (object_a)
      if (runif(1) < 0.5) { aa <- pr[i, 1]; bb <- pr[i, 2] }
      else { aa <- pr[i, 2]; bb <- pr[i, 1] }
      p <- plogis(beta[aa] - beta[bb] + gamma)   # object_a first -> +gamma
      rows[[length(rows) + 1L]] <- data.frame(
        a = aa, b = bb, win = if (runif(1) < p) aa else bb)
    }
    do.call(rbind, rows)
  }
  f <- btl(sim_pos(0.5), "a", "b", winner = "win", position = TRUE)
  pos <- f$dependence[f$dependence$effect == "position", ]
  expect_identical(f$dependence$effect, "position")
  expect_lt(abs(pos$estimate - 0.5), 3 * pos$se)          # within 3 SE of 0.5
  expect_equal(pos$n_informative, f$n_comparisons)         # every row informative
  # locations still recover, and stay sum-zero
  loc <- f$objects$location[match(objs, f$objects$object)]
  expect_gt(cor(loc, beta), 0.95)
  expect_equal(sum(f$objects$location), 0, tolerance = 1e-8)

  # planted zero: the coefficient is near zero and not significant
  set.seed(202)
  f0 <- btl(sim_pos(0), "a", "b", winner = "win", position = TRUE)
  pos0 <- f0$dependence[f0$dependence$effect == "position", ]
  expect_lt(abs(pos0$t), 2)

  # identified through triangle closure even with a FIXED orientation (each
  # pair always presented in the same order); the separation guard leaves it in
  set.seed(707)
  d2 <- data.frame(a = rep(pr[, 1], each = 120), b = rep(pr[, 2], each = 120))
  d2$win <- ifelse(runif(nrow(d2)) < plogis(beta[d2$a] - beta[d2$b] + 0.5),
                   d2$a, d2$b)
  ff <- btl(d2, "a", "b", winner = "win", position = TRUE)
  fp <- ff$dependence[ff$dependence$effect == "position", ]
  expect_false(any(grepl("separated", ff$notes)))
  expect_lt(abs(fp$estimate), 10)
  expect_gt(fp$t, 2)
  # the print label is the positional one, not "Within-judge"
  expect_output(print(ff), "First-position advantage")
  position_summary <- fit_summary_table(ff)
  expect_true(any(grepl("First-position advantage",
                        position_summary$statistic, fixed = TRUE)))
  expect_false(any(grepl("Within-judge position",
                         position_summary$statistic, fixed = TRUE)))
})

test_that("position effects are refused when confounded with object locations", {
  d <- rbind(
    data.frame(a = "A", b = "B", win = rep(c("A", "B"), 20)),
    data.frame(a = "B", b = "C", win = rep(c("B", "C"), 20)))
  expect_no_error(btl(d, "a", "b", winner = "win"))
  expect_error(btl(d, "a", "b", winner = "win", position = TRUE),
               "not separately identified")
})

test_that("anchored estimation reproduces the free scale and equates panels", {
  set.seed(303)
  beta <- c(A = -1.2, B = -0.4, C = 0.1, D = 0.6, E = 0.9)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 60), b = rep(pr[, 2], each = 60))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  ff <- btl(d, "a", "b", winner = "win")
  loc_free <- setNames(ff$objects$location, ff$objects$object)

  # (a) anchoring two objects AT their free-fit values reproduces the free fit
  fa <- btl(d, "a", "b", winner = "win", anchors = loc_free[c("B", "D")])
  loc_anc <- setNames(fa$objects$location, fa$objects$object)
  expect_equal(unname(loc_anc[names(loc_free)]), unname(loc_free),
               tolerance = 1e-6)
  # (c) anchored rows are exact and report se 0
  expect_equal(unname(loc_anc["B"]), unname(loc_free["B"]), tolerance = 1e-10)
  expect_equal(unname(loc_anc["D"]), unname(loc_free["D"]), tolerance = 1e-10)
  expect_equal(fa$objects$se[fa$objects$object %in% c("B", "D")], c(0, 0))
  expect_setequal(names(fa$anchors), c("B", "D"))
  expect_output(print(fa), "Anchored at 2 object")

  fa_pad <- btl(d, "a", "b", winner = "win",
                anchors = stats::setNames(loc_free[c("B", "D")],
                                          c(" B ", "D ")))
  expect_setequal(names(fa_pad$anchors), c("B", "D"))
  expect_equal(fa_pad$objects$location, fa$objects$location,
               tolerance = 1e-8)
  expect_error(
    btl(d, "a", "b", winner = "win", anchors = c(B = 0, " B " = 1)),
    "after trimming")

  # (b) two-panel equating: overlapping objects anchor panel 2 onto panel 1's
  # scale; the non-common objects then land at their true spacing
  truth <- c(A = -1.5, B = -0.7, C = 0, D = 0.7, E = 1.5, F = 2.1)
  mkpanel <- function(objs, nper, seed) {
    set.seed(seed)
    pp <- t(combn(objs, 2))
    dd <- data.frame(a = rep(pp[, 1], each = nper), b = rep(pp[, 2], each = nper))
    dd$win <- ifelse(runif(nrow(dd)) < plogis(truth[dd$a] - truth[dd$b]),
                     dd$a, dd$b)
    dd
  }
  f1 <- btl(mkpanel(c("A", "B", "C", "D"), 80, 11), "a", "b", winner = "win")
  l1 <- setNames(f1$objects$location, f1$objects$object)
  f2 <- btl(mkpanel(c("C", "D", "E", "F"), 80, 22), "a", "b", winner = "win",
            anchors = l1[c("C", "D")])
  l2 <- setNames(f2$objects$location, f2$objects$object)
  expect_equal(unname(l2["C"]), unname(l1["C"]), tolerance = 1e-10)  # anchored exact
  expect_equal(unname(l2["D"]), unname(l1["D"]), tolerance = 1e-10)
  # E and F recover their true spacing relative to the anchored D
  expect_lt(abs((l2["E"] - l2["D"]) - (truth["E"] - truth["D"])), 0.35)
  expect_lt(abs((l2["F"] - l2["D"]) - (truth["F"] - truth["D"])), 0.35)

  # guard rails
  expect_error(btl(d, "a", "b", winner = "win", anchors = c(ZZ = 1)),
               "do not match any object")
  # a single misspelled name among valid ones must also error, not be
  # silently dropped leaving that object free
  expect_error(btl(d, "a", "b", winner = "win",
                   anchors = c(A = -1, Z_typo = 0.5)),
               "do not match any object")
  expect_error(btl(d, "a", "b", winner = "win", anchors = c(1, 2)),
               "named numeric")
  bad_name <- c(0, 1); names(bad_name) <- c("A", NA_character_)
  expect_error(btl(d, "a", "b", winner = "win", anchors = bad_name),
               "named numeric")
  expect_error(btl(d, "a", "b", winner = "win", anchors = c(" " = 0)),
               "named numeric")
  # an anchored boundary object is an error, not silent removal
  set.seed(9)
  db <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
  db$win <- ifelse(runif(nrow(db)) < plogis(beta[db$a] - beta[db$b]),
                   db$a, db$b)
  db$win[db$a == "A" | db$b == "A"] <- "A"           # A undefeated
  expect_error(btl(db, "a", "b", winner = "win", anchors = c(A = 0, C = 0.1)),
               "boundary")
  # the same boundary object, left free, is still removed with a note
  fb <- btl(db, "a", "b", winner = "win", anchors = c(C = 0.1))
  expect_true(isTRUE(fb$objects$extreme[fb$objects$object == "A"]))
  expect_true(is.na(fb$objects$se[fb$objects$object == "A"]))
  expect_true(any(grepl("boundary", fb$notes)))
})

test_that("an anchor is not silently dropped with its only comparison", {
  d <- data.frame(
    a = c(rep("A", 20), rep("B", 20), "D"),
    b = c(rep("B", 20), rep("C", 20), "A"),
    winner = c(rep(c("A", "B"), 10), rep(c("B", "C"), 10), "D"),
    count = c(rep(1L, 40), 0L))
  expect_error(
    btl(d, "a", "b", winner = "winner", count = "count",
        anchors = c(A = 0, D = 1)),
    "anchored object.*no usable comparisons.*D")
})

test_that("position and order covariates are estimated together", {
  set.seed(505)
  K <- 6; objs <- sprintf("O%d", 1:K)
  beta <- setNames(seq(-1, 1, length.out = K), objs); beta <- beta - mean(beta)
  jids <- sprintf("J%02d", 1:12); pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
  d$judge <- sample(jids, nrow(d), TRUE)
  d <- d[order(d$judge), ]; d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  # randomise orientation so position is well identified
  flip <- runif(nrow(d)) < 0.5
  aa <- ifelse(flip, d$b, d$a); bb <- ifelse(flip, d$a, d$b)
  d$a <- aa; d$b <- bb
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b] + 0.4), d$a, d$b)
  f <- btl(d, "a", "b", winner = "win", judge = "judge", order = "t",
           position = TRUE)
  expect_true(f$converged)
  expect_setequal(f$dependence$effect, c("exposure", "carry_over", "position"))
  # both covariate machineries are present in the row-aligned tables
  expect_true(all(c("exposure", "carry_over", "position") %in%
                    names(f$comparisons)))
  expect_true(all(c("exposure", "carry_over", "position") %in%
                    names(f$dependence_data)))
  # every row informs position; not every row informs exposure/carry-over
  pos <- f$dependence[f$dependence$effect == "position", ]
  expect_equal(pos$n_informative, f$n_comparisons)
  expect_gt(pos$t, 2)   # the planted first-position advantage is detected
  # Surprise diagnostics must use the complete fitted predictor, not only
  # the object-location difference.
  mo <- .btl_fitted_moments(f, f$comparisons)
  bl <- setNames(f$objects$location, f$objects$object)
  Z <- as.matrix(f$comparisons[, f$dependence$effect, drop = FALSE])
  expected_lp <- bl[f$comparisons$object_a] - bl[f$comparisons$object_b] +
    drop(Z %*% f$dependence$estimate)
  expect_equal(mo$lp, unname(expected_lp), tolerance = 1e-12)
  expect_gt(max(abs(mo$lp -
    (bl[f$comparisons$object_a] - bl[f$comparisons$object_b]))), 0.1)
})

test_that("count-weighted rows give the SAME standard errors as expanded rows", {
  # the sandwich meat must treat a count=w row as w independent comparisons
  # (w * (x-E)^2), not one cluster of weight w ((w*(x-E))^2): the latter
  # inflated every SE by ~sqrt(w)
  set.seed(1); beta <- c(A = -1, B = -0.2, C = 0.5, D = 0.7)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 12), b = rep(pr[, 2], each = 12))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  agg <- aggregate(cbind(k = rep(1, nrow(d))) ~ a + b + win, d, sum)
  fe <- btl(d, "a", "b", "win")
  fw <- btl(agg, "a", "b", "win", count = "k")
  expect_equal(fw$objects$location, fe$objects$location, tolerance = 1e-10)
  expect_equal(fw$objects$se, fe$objects$se, tolerance = 1e-10)
  expect_equal(fw$osi$PSI, fe$osi$PSI, tolerance = 1e-8)
})

test_that("anchoring converges regardless of the anchor origin", {
  set.seed(2); objs <- paste0("O", 1:5)
  b <- setNames(seq(-1, 1, length.out = 5), objs)
  pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
  d$win <- ifelse(runif(nrow(d)) < plogis(b[d$a] - b[d$b]), d$a, d$b)
  free <- btl(d, "a", "b", "win")
  for (delta in c(3, 10)) {
    anc <- setNames(free$objects$location[match(c("O1", "O5"),
                                                free$objects$object)] + delta,
                    c("O1", "O5"))
    fa <- btl(d, "a", "b", "win", anchors = anc)
    expect_true(fa$converged)
    # anchoring at translated values IS a pure translation of the free fit
    expect_equal(fa$objects$location, free$objects$location + delta,
                 tolerance = 1e-8)
  }
})

test_that("btl_next_pairs one-step priority beats the lowest-priority pair", {
  set.seed(9); objs <- paste0("E", 1:6)
  b <- setNames(seq(-1.2, 1.2, length.out = 6), objs)
  pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 15), b = rep(pr[, 2], each = 15))
  d$win <- ifelse(runif(nrow(d)) < plogis(b[d$a] - b[d$b]), d$a, d$b)
  d <- d[!(d$a == "E3" | d$b == "E3") | seq_len(nrow(d)) %% 5 == 0, ]  # starve E3
  f <- btl(d, "a", "b", "win")
  np <- btl_next_pairs(f, n = 15)
  addvar <- function(pa, pb) {
    set.seed(77)
    ex <- data.frame(a = rep(pa, 25), b = rep(pb, 25))
    ex$win <- ifelse(runif(25) < plogis(b[ex$a] - b[ex$b]), ex$a, ex$b)
    ff <- btl(rbind(d[, c("a", "b", "win")], ex), "a", "b", "win")
    sum(ff$objects$se^2)
  }
  expect_lt(addvar(np$object_a[1], np$object_b[1]),
            addvar(np$object_a[nrow(np)], np$object_b[nrow(np)]))
  expect_true("E3" %in% c(np$object_a[1], np$object_b[1]))
})

test_that("model-based BTL diagnostics refuse an unconverged calibration", {
  set.seed(87)
  objs <- LETTERS[1:5]
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 20),
                  b = rep(pr[, 2], each = 20),
                  judge = rep(sprintf("J%02d", 1:10), length.out = 200))
  d$win <- ifelse(stats::runif(nrow(d)) < .5, d$a, d$b)
  f <- btl(d, "a", "b", "win", judge = "judge")
  limited <- suppressWarnings(
    btl(d, "a", "b", "win", judge = "judge", maxit = 1L))
  expect_false(limited$converged)
  expect_true(all(is.na(limited$objects$se)))
  expect_true(is.na(limited$osi$PSI))
  expect_false(limited$cl$inference_available)
  expect_true(is.na(limited$cl$eff_params))
  expect_true(all(is.na(limited$cov_parameters)))
  f$converged <- FALSE
  expect_error(btl_information(f), "did not converge")
  expect_error(btl_next_pairs(f), "did not converge")
  expect_error(btl_dimensionality(f, reps = 20), "did not converge")
  expect_error(judge_surprise(f, "J01"), "did not converge")
  expect_error(judge_pair_surprise(f, "J01"), "did not converge")
  expect_error(plot_btl(f), "did not converge")
  expect_error(plot_btl_categories(f), "did not converge")
  expect_error(plot_btl_icc(f, "A"), "did not converge")
  expect_error(plot_btl_dependence(f), "did not converge")
})

test_that("BTL result plots refuse the wrong result class directly", {
  expect_error(plot_btl_transitivity(list()), "btl_transitivity")
  expect_error(plot_btl_scree(list()), "btl_dimensionality")
  expect_error(plot_btl_dim_map(list()), "btl_dimensionality")
})

test_that("BTL DIF does not redefine an externally anchored object", {
  set.seed(49)
  beta <- c(A = -1, B = -.5, C = 0, D = .5, E = 1)
  pr <- t(combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = 100),
                  b = rep(pr[, 2], each = 100))
  judges <- sprintf("J%02d", 1:20)
  d$judge <- rep(judges, length.out = nrow(d))
  grp <- setNames(rep(c("g1", "g2"), each = 10), judges)
  shift <- ifelse(grp[d$judge] == "g2" & d$a == "C", 2,
            ifelse(grp[d$judge] == "g2" & d$b == "C", -2, 0))
  p <- plogis(beta[d$a] - beta[d$b] + shift)
  d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
  f <- btl(d, "a", "b", "win", judge = "judge", anchors = c(C = 0))
  z <- btl_dif(f, grp, objects = "C")
  expect_true(z$summary$uniform_DIF)
  expect_null(z$sizes)
  expect_true(any(grepl("externally anchored", z$notes)))
})

test_that("named judge factors are matched by judge before row length", {
  set.seed(731)
  objects <- LETTERS[1:5]
  pairs <- t(combn(objects, 2))
  d <- data.frame(a = rep(pairs[, 1], length.out = 120),
                  b = rep(pairs[, 2], length.out = 120),
                  judge = sprintf("J%03d", 1:120))
  beta <- setNames(seq(-1, 1, length.out = 5), objects)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]),
                  d$a, d$b)
  fit <- btl(d, "a", "b", "win", judge = "judge")
  group <- setNames(rep(c("g1", "g2"), each = 60), d$judge)
  a <- btl_dif(fit, list(group = group))
  b <- btl_dif(fit, list(group = group[sample(names(group))]))
  expect_equal(a$terms, b$terms)
  expect_equal(a$summary, b$summary)
})

test_that("a boundary object is reported at an extrapolated location", {
  set.seed(701)
  d <- simulate_btl(n_objects = 8, n_judges = 30, reps_per_pair = 4)
  # Make O1 deterministically winless; relying on one seed to happen to
  # produce a boundary object makes this a simulator-RNG test instead.
  has_o1 <- d$object_a == "O1" | d$object_b == "O1"
  d$winner[has_o1] <- ifelse(d$object_a[has_o1] == "O1",
                              d$object_b[has_o1], d$object_a[has_o1])
  f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  ext <- f$objects[f$objects$object == "O1", ]
  cal <- f$objects[!f$objects$extreme, ]
  expect_identical(ext$object, "O1")
  expect_true(ext$extreme)
  # winless: below every calibrated location, finite, no SE, no fit
  expect_true(is.finite(ext$location))
  expect_lt(ext$location, min(cal$location))
  expect_true(is.na(ext$se))
  expect_true(is.na(ext$fit_resid))
  expect_match(paste(f$notes, collapse = " "), "extrapolated location")
  expect_true("O1" %in% c(f$observed_comparisons$object_a,
                          f$observed_comparisons$object_b))
  expect_setequal(f$observed_comparisons$judge, d$judge)
  trn <- btl_transitivity(f)
  expect_equal(trn$summary$n_objects, 8L)
  expect_equal(trn$summary$n_pairs, choose(8L, 2L))
  expect_true("O1" %in% trn$objects$object)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_btl_icc(f, " O2 ", min_n = 1))
  expect_error(plot_btl_icc(f, "O1", min_n = 1), "response boundary")
  judges <- sort(unique(f$comparisons$judge))
  group <- stats::setNames(rep(c("g1", "g2"), length.out = length(judges)),
                           judges)
  bd <- btl_dif(f, group, min_n = 2)
  expect_false("O1" %in% bd$summary$object)
  expect_match(paste(bd$notes, collapse = " "), "extrapolations")
  expect_error(btl_dif(f, group, objects = "O1", min_n = 2),
               "response boundary")
  one_judge <- f$comparisons$judge[1]
  expect_false("O1" %in% names(judge_surprise(
    f, paste0(" ", one_judge, " "), min_n = 1)$all_locations))
  expect_false("O1" %in% names(judge_pair_surprise(
    f, paste0(" ", one_judge, " "), min_n = 1)$all_locations))
  # the extrapolated row takes no part in equating
  set.seed(9000)
  d2 <- simulate_btl(n_objects = 8, n_judges = 30, reps_per_pair = 2)
  f2 <- btl(d2, "object_a", "object_b", winner = "winner", judge = "judge")
  eq <- btl_equate(f, f2, independent = TRUE)
  expect_false("O1" %in% eq$table$object)
})

test_that("the observed-ICC display survives a fully omitted comparator set", {
  set.seed(21)
  d <- simulate_btl(n_objects = 8, n_judges = 30, reps_per_pair = 2)
  f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  # every opponent falls below the default min_n: the model curve and the
  # omission note still draw, with no visible points
  expect_no_error(plot_btl_icc(f, "O2"))
})

test_that("the grouped BTL ICC refuses an empty grouping display", {
  d <- simulate_btl(n_objects = 8, n_judges = 30, reps_per_pair = 10,
                    seed = 22)
  f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_error(plot_btl_icc(f, "O2", group = rep(NA_character_,
                                                   nrow(f$comparisons))),
               "no observed judge-group values")
  blank <- stats::setNames(rep("  ", length(unique(f$comparisons$judge))),
                           unique(f$comparisons$judge))
  expect_error(plot_btl_icc(f, "O2", group = blank),
               "no observed judge-group values")

  group <- ifelse(f$comparisons$object_a == "O2" |
                    f$comparisons$object_b == "O2", NA_character_, "g1")
  expect_error(plot_btl_icc(f, "O2", group = group),
               "for object 'O2'")
})

test_that("the dependence plot's baseline carries a fitted position effect", {
  d <- simulate_btl(8, 14, reps_per_pair = 12,
                    dependence = list(exposure = 0.6), seed = 42)
  f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
           order = "order", position = TRUE)
  expect_true("position" %in% f$dependence$effect)
  # the partial residual under the FULL fitted model centres on zero; a
  # baseline that omits the fitted position coefficient displaces every
  # point through the nonlinear expected-score map
  dd <- f$dependence_data
  bl <- setNames(f$objects$location, f$objects$object)
  tau <- if (!is.null(f$thresholds)) f$thresholds$tau else numeric(1)
  dep <- setNames(f$dependence$estimate, f$dependence$effect)
  Em <- function(v) vapply(v, function(t) item_moments(t, tau)$E, 0)
  base <- unname(bl[dd$object_a] - bl[dd$object_b])
  lin <- base + dep[["exposure"]] * dd$exposure +
    dep[["carry_over"]] * dd$carry_over + dep[["position"]] * dd$position
  expect_lt(abs(mean(dd$response - Em(lin))), 0.005)
  expect_gt(abs(mean(dd$response - Em(lin - dep[["position"]] * dd$position))),
            0.02)
  # and the plot itself runs on a position-fitted model
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_btl_dependence(f, "exposure"))
  expect_no_error(plot_btl_dependence(f, "carry_over"))
})
