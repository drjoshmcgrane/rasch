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
  expect_true(any(grepl("no wins or no losses", ft3$notes)))
  expect_false("E" %in% ft3$objects$object)
  # disconnected comparison graphs are refused with the components listed
  dd <- data.frame(a = c("A", "A", "C", "C"), b = c("B", "B", "D", "D"),
                   win = c("A", "B", "C", "D"))
  expect_error(btl(dd, "a", "b", "win"), "disconnected")
})

test_that("plot_btl draws and print method runs", {
  beta <- c(A = -1, B = 0, C = 0.4, D = 0.6)
  ft <- btl(sim_btl(beta, 40, seed = 5), "a", "b", "win")
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_btl(ft))
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
  # ordered-factor input maps by level order and keeps the labels
  d$lab <- factor(c("worse", "tie", "better")[d$grade + 1],
                  levels = c("worse", "tie", "better"))
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
                     levels = c("a little", "much"))
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
                levels = paste0("L", 0:6))
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
  jids <- sprintf("J%02d", 1:20)
  grp <- setNames(rep(c("g1", "g2"), each = 10), jids)
  pr <- t(combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 14), b = rep(pr[, 2], each = 14))
  d$judge <- sample(jids, nrow(d), TRUE)
  shift <- ifelse(grp[d$judge] == "g2" & d$a == "S06", 1,
           ifelse(grp[d$judge] == "g2" & d$b == "S06", -1, 0))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b] + shift),
                  d$a, d$b)
  f <- btl(d, "a", "b", winner = "win", judge = "judge")
  dif <- btl_dif(f, groups = grp)
  expect_s3_class(dif, "rmt_btl_dif")
  # ANOVA route: only S06 flagged
  expect_true(dif$anova$uniform_DIF[dif$anova$object == "S06"])
  expect_equal(sum(dif$anova$uniform_DIF), 1L)
  # magnitude route: right size, right object, nothing else
  s6 <- dif$sizes[dif$sizes$object == "S06", ]
  expect_true(s6$significant && s6$practical)
  expect_lt(abs(abs(s6$difference) - 1), 3 * s6$se)
  expect_equal(sum(dif$sizes$significant), 1L)
  # grouped characteristic curve renders
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_btl_icc(f, "S06", group = grp))
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
  expect_true(any(grepl("Within-judge exposure", ft$statistic)))
  expect_true(any(grepl("Within-judge carry-over", ft$statistic)))
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
