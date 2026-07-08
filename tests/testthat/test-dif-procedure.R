# The full DIF procedure: single-factor two-way ANOVA with the
# class-interval row and effect sizes; x-way factorial with interactions;
# post-hoc pairwise comparisons for significant interactions and for
# multi-level main effects, familywise-corrected; and DIF magnitudes in
# logits for practical significance.

sim_dif <- function(n = 900, seed = 2, shifts) {
  # shifts: named list group-label -> per-item logit shift vector
  set.seed(seed)
  d <- seq(-1.5, 1.5, length.out = 8)
  g <- factor(rep(names(shifts), each = n / length(shifts)))
  th <- rnorm(n)
  X <- matrix(NA_integer_, n, 8)
  for (lv in names(shifts))
    X[g == lv, ] <- matrix(rbinom(sum(g == lv) * 8, 1,
      plogis(outer(th[g == lv], d + shifts[[lv]], "-"))), sum(g == lv), 8)
  colnames(X) <- paste0("I", 1:8)
  list(X = X, g = g)
}

test_that("dif_anova reports the full two-way table with effect sizes", {
  s <- sim_dif(shifts = list(a = rep(0, 8), b = c(0, 0, 1, rep(0, 5))))
  fit <- rasch(data.frame(s$X, grp = s$g), factors = "grp")
  da <- dif_anova(fit)
  su <- da$summary
  expect_true(all(c("eta2_uniform", "eta2_nonuniform", "uniform_DIF",
                    "nonuniform_DIF") %in% names(su)))
  # one row per item (single factor); the planted item has the largest effect
  expect_equal(which.max(su$eta2_uniform), 3L)
  expect_gt(su$eta2_uniform[3], max(su$eta2_uniform[-3]) * 3)
  expect_true(su$uniform_DIF[3])
  expect_true(all(su$eta2_uniform > 0 & su$eta2_uniform < 1, na.rm = TRUE))
  # the class-interval main effect lives in the full terms table
  expect_true(all(is.finite(da$terms$F_value[da$terms$term == "ci"])))
  # familywise option flows through
  db <- dif_anova(fit, p_adjust = "bonferroni")
  expect_true(all(db$summary$p_uniform_adj >= su$p_uniform_adj - 1e-12,
                  na.rm = TRUE))
})

test_that("dif_size recovers a planted uniform DIF in logits", {
  s <- sim_dif(shifts = list(a = rep(0, 8), b = c(0, 0, 0.8, rep(0, 5))))
  fit <- rasch(data.frame(s$X, grp = s$g), factors = "grp")
  ds <- dif_size(fit, "I3", by = "grp")
  expect_equal(nrow(ds$pairs), 1)
  expect_lt(abs(abs(ds$pairs$difference) - 0.8), 0.25)
  expect_true(ds$pairs$significant)
  expect_true(ds$pairs$practical)          # 0.8 > 0.5 logits
  expect_true(ds$pairs$lower < ds$pairs$difference &
              ds$pairs$difference < ds$pairs$upper)
  # a clean item shows neither statistical nor practical DIF
  ds0 <- dif_size(fit, "I6", by = "grp")
  expect_false(ds0$pairs$significant)
  expect_false(ds0$pairs$practical)
})

test_that("multi-level factors get familywise pairwise comparisons in logits", {
  s <- sim_dif(n = 1200, seed = 5,
               shifts = list(a = rep(0, 8),
                             b = c(0, 0.5, rep(0, 6)),
                             c = c(0, 1.0, rep(0, 6))))
  fit <- rasch(data.frame(s$X, grp = s$g), factors = "grp")
  ds <- dif_size(fit, "I2", by = "grp")
  expect_equal(nrow(ds$levels), 3)
  expect_equal(nrow(ds$pairs), 3)          # all pairs of three levels
  d_ac <- ds$pairs$difference[ds$pairs$level_a == "a" & ds$pairs$level_b == "c"]
  d_ab <- ds$pairs$difference[ds$pairs$level_a == "a" & ds$pairs$level_b == "b"]
  expect_lt(abs(abs(d_ac) - 1.0), 0.35)
  expect_lt(abs(abs(d_ab) - 0.5), 0.35)
  expect_gt(abs(d_ac), abs(d_ab))          # graded shifts recovered in order
  # Holm adjustment is monotone in the raw p over the family
  expect_true(all(ds$pairs$p_adj >= ds$pairs$p - 1e-12))
  # the extreme pair is flagged practical, and a-c also significant
  sel <- ds$pairs$level_a == "a" & ds$pairs$level_b == "c"
  expect_true(ds$pairs$significant[sel] && ds$pairs$practical[sel])
})

test_that("factorial procedure: interaction post-hocs and sizes for significant terms", {
  set.seed(7); n <- 1600
  d <- seq(-1.5, 1.5, length.out = 6)
  g1 <- factor(rep(c("a", "b"), each = n / 2))
  g2 <- factor(rep(c("x", "y"), times = n / 2))
  th <- rnorm(n)
  sh <- ifelse(g1 == "b" & g2 == "y", 1.1, 0)   # DIF in one cell only
  X <- matrix(rbinom(n * 6, 1,
    plogis(outer(th, d, "-") - outer(sh, c(0, 0, 1, 0, 0, 0)))), n, 6)
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))
  fa <- dif_anova(fit, sizes = TRUE, effects = "factorial")

  # the g1:g2 interaction is significant for the planted item and
  # supersedes the main effects it involves
  t3 <- fa$terms[fa$terms$item == "I3", ]
  expect_true(t3$significant[t3$term == "g1:g2"])
  expect_true(all(c("eta2_partial") %in% names(fa$terms)))
  expect_true(t3$eta2_partial[t3$term == "g1:g2"] > 0)
  sup <- t3$term[t3$superseded]
  expect_true(all(sup %in% c("g1", "g2")))

  # Tukey post-hocs exist for the interaction cells of the planted item
  tk3 <- fa$tukey[fa$tukey$item == "I3" & fa$tukey$term == "g1:g2", ]
  expect_equal(nrow(tk3), 6)               # 4 cells -> 6 pairs
  worst <- tk3$comparison[which.min(tk3$p_tukey)]
  expect_true(grepl("b:y", worst))

  # sizes: logit magnitudes for the significant term, the b:y cell apart
  sz <- fa$sizes[fa$sizes$item == "I3" & fa$sizes$term == "g1:g2", ]
  expect_gt(nrow(sz), 0)
  by_pairs <- sz[sz$level_a == "b:y" | sz$level_b == "b:y", ]
  other_pairs <- sz[!(sz$level_a == "b:y" | sz$level_b == "b:y"), ]
  expect_gt(min(abs(by_pairs$difference)), max(abs(other_pairs$difference)))
  expect_lt(abs(max(abs(by_pairs$difference)) - 1.1), 0.4)
  expect_true(any(by_pairs$practical))
  # clean items produce no size rows
  expect_false(any(fa$sizes$item == "I5"))
})

test_that("dif_size guards: thin levels dropped, unknown factor errors", {
  s <- sim_dif(shifts = list(a = rep(0, 8), b = rep(0, 8)))
  g3 <- as.character(s$g); g3[1:5] <- "tiny"
  fit <- rasch(data.frame(s$X, grp = factor(g3)), factors = "grp")
  ds <- dif_size(fit, "I1", by = "grp", min_n = 20)
  expect_true(any(grepl("tiny", ds$notes)))
  expect_equal(nrow(ds$levels), 2)
  expect_error(dif_size(fit, "I1", by = "nofactor"), "factor")
})

test_that("person_extrapolated continues the score table geometrically", {
  set.seed(23)
  d <- seq(-2, 2, length.out = 10)
  X <- matrix(rbinom(400 * 10, 1, plogis(outer(rnorm(400, 0, 2.2), d, "-"))),
              400, 10)
  colnames(X) <- paste0("I", 1:10)
  fit <- rasch(X)
  skip_if(!any(fit$person$extreme), "no extreme persons in this draw")
  pe <- person_extrapolated(fit)
  st <- score_table(fit, extremes = "extrapolated")
  top <- pe$extreme & pe$raw == pe$max_raw & pe$n_items == 10
  bot <- pe$extreme & pe$raw == 0 & pe$n_items == 10
  if (any(top)) {
    expect_equal(unique(pe$theta_extrapolated[top]),
                 st$theta[nrow(st)], tolerance = 1e-8)
    # the geometric continuation differs from Warm's own extreme value
    # (Warm's finite extreme step exceeds the continued interior steps)
    expect_false(isTRUE(all.equal(unique(pe$theta_extrapolated[top]),
                                  unique(pe$theta[top]))))
    # still beyond the last interior score's measure
    expect_true(all(pe$theta_extrapolated[top] >
                    st$theta[nrow(st) - 1]))
  }
  if (any(bot)) {
    expect_equal(unique(pe$theta_extrapolated[bot]), st$theta[1],
                 tolerance = 1e-8)
    expect_true(all(pe$theta_extrapolated[bot] < st$theta[2]))
  }
  # non-extreme persons unchanged
  ne <- !pe$extreme
  expect_equal(pe$theta_extrapolated[ne], pe$theta[ne])
  expect_equal(pe$se_extrapolated[ne], pe$se[ne])
})

test_that("MFRM facet fit reports margin and pooled statistics with df", {
  set.seed(31)
  simP <- function(th, tau) { x <- 0:length(tau)
    p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
  persons <- sprintf("P%03d", 1:150); raters <- paste0("R", 1:3)
  th <- setNames(rnorm(150, 0, 1.2), persons)
  rho <- setNames(c(-0.5, 0, 0.5), raters)
  tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4))
  dd <- expand.grid(person = persons, item = names(tau), rater = raters,
                    stringsAsFactors = FALSE)
  dd$score <- mapply(function(p, i, r)
    sample(0:2, 1, prob = simP(th[p], tau[[i]] + rho[r])),
    dd$person, dd$item, dd$rater)
  mf <- rasch_mfrm(dd, person = "person", item = "item", score = "score",
                   facets = "rater")
  fe <- mf$facet_effects$rater
  expect_true(all(c("fit_resid", "fit_resid_pooled", "df_fit") %in% names(fe)))
  expect_true(all(is.finite(fe$fit_resid)))
  # the margin statistic: the mean of the level's virtual items' fits
  vsel <- mf$virtual_map$rater == fe$level[1]
  expect_equal(fe$fit_resid[1],
               mean(mf$items$fit_resid[vsel], na.rm = TRUE),
               tolerance = 1e-10)
  # pooled df equals the cell df factor times the level's cells
  expect_equal(fe$df_fit, mf$summary_stats$df_factor * fe$n, tolerance = 1e-8)
  expect_true(all(is.finite(mf$item_effects$df_fit)))
  expect_true(all(is.finite(mf$item_effects$fit_resid_pooled)))
})

test_that("the factorial summary pivots to uniform/non-uniform per group term", {
  set.seed(1); n <- 600
  g1 <- rep(c("m", "f"), n / 2); g2 <- sample(c("young", "old"), n, TRUE)
  d <- seq(-1.5, 1.5, length.out = 6)
  sh <- matrix(0, n, 6); sh[g1 == "f", 3] <- 0.8
  X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, sex = g1, age = g2), factors = c("sex", "age"))
  fa <- dif_anova(fit, effects = "factorial")
  s <- fa$summary
  # one row per item and group term (no ci terms, no residual row)
  expect_setequal(unique(s$term), c("sex", "age", "sex:age"))
  expect_equal(nrow(s), 6 * 3)
  # the pivot agrees with the full table
  u <- fa$terms[fa$terms$item == "I3" & fa$terms$term == "sex", ]
  nu <- fa$terms[fa$terms$item == "I3" & fa$terms$term == "sex:ci", ]
  r <- s[s$item == "I3" & s$term == "sex", ]
  expect_equal(r$F_uniform, u$F_value)
  expect_equal(r$p_uniform_adj, u$p_adj)
  expect_equal(r$F_nonuniform, nu$F_value)
  # the planted uniform DIF is flagged as uniform, not non-uniform
  expect_true(r$uniform_DIF)
  # no misfit flag on the items table any more
  expect_false("misfit" %in% names(fit$items))
})

test_that("DIF class intervals adapt to the cells each analysis uses", {
  set.seed(1); n <- 800
  g1 <- rep(c("m", "f"), n / 2)
  g2 <- sample(c("a", "b", "c", "d"), n, TRUE)
  d <- seq(-1.5, 1.5, length.out = 8)
  X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-"))), n, 8)
  colnames(X) <- paste0("I", 1:8)
  f <- rasch(data.frame(X, sex = g1, age = g2), factors = c("sex", "age"))
  # the joint model requests one interval count from the smallest cell used
  # (a single factor from that factor's smallest group); the returned
  # n_groups is the count actually FORMED, which quantile ties can reduce
  ds <- dif_anova(f, factors = "sex")
  req1 <- max(2L, min(10L, min(table(g1[!is.na(f$person$theta)])) %/% 30L))
  expect_lte(ds$n_groups, req1)
  expect_gte(ds$n_groups, 2L)
  # several factors from the smallest factor-combination cell
  da <- dif_anova(f)
  cells <- interaction(g1, g2, drop = TRUE)
  req2 <- max(2L, min(10L, min(table(cells[!is.na(f$person$theta)])) %/% 30L))
  expect_lte(da$n_groups, req2)
  expect_gt(ds$n_groups, da$n_groups)   # combination cells are smaller
  # explicit n_groups still overrides (and is achievable here)
  expect_equal(dif_anova(f, n_groups = 5)$n_groups, 5)
})

test_that("dif_anova tests a within-subject factor against person clustering", {
  set.seed(7)
  N <- 300; d <- seq(-1.5, 1.5, length.out = 8); th <- rnorm(N)
  gen <- function(shift) {
    s <- matrix(0, N, 8); s[, 3] <- shift
    matrix(rbinom(N * 8, 1, plogis(outer(th, d, "-") - s)), N, 8)
  }
  X <- rbind(gen(0), gen(0.9)); colnames(X) <- paste0("I", 1:8)
  dat <- data.frame(X, occasion = rep(c("t1", "t2"), each = N))
  id <- rep(sprintf("P%03d", 1:N), 2)
  fit <- rasch(dat, factors = "occasion", id = id)

  da <- dif_anova(fit)
  # the repeated occasion factor is auto-detected as within-subject
  expect_identical(da$within, "occasion")
  su <- da$summary
  # the planted within-subject DIF on I3 is recovered, clean items are null
  expect_true(su$uniform_DIF[su$item == "I3"])
  expect_equal(sum(su$uniform_DIF), 1L)
  # the within uniform test equals the person-level paired test (its gold
  # standard) up to the class-interval filtering
  z <- fit$residuals[, 3]; g <- factor(fit$factors$occasion); pid <- factor(id)
  dp <- tapply(z[g == "t2"], pid[g == "t2"], mean) -
        tapply(z[g == "t1"], pid[g == "t1"], mean)
  expect_lt(abs(su$F_uniform[su$item == "I3"] -
                t.test(dp)$statistic^2), 3)

  # a cross-sectional design (unique ids) is unchanged: nothing within,
  # and the ordinary between-subjects group DIF still flags the shifted item
  set.seed(1); gg <- rep(c("a", "b"), each = 250)
  sh <- matrix(0, 500, 8); sh[gg == "b", 3] <- 0.8
  Xc <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-") - sh)),
               500, 8)
  colnames(Xc) <- paste0("I", 1:8)
  dc <- dif_anova(rasch(data.frame(Xc, grp = gg), factors = "grp"))
  expect_length(dc$within, 0L)
  expect_true(dc$summary$uniform_DIF[dc$summary$item == "I3"])
})

test_that("factorial DIF uses a mixed ANOVA when a factor is within-subject", {
  set.seed(11)
  N <- 350; d <- seq(-1.5, 1.5, length.out = 8); th <- rnorm(N)
  sex <- rep(c("m", "f"), length.out = N)
  gen <- function(occ_shift) {
    s <- matrix(0, N, 8)
    s[sex == "f", 4] <- 0.8      # between DIF on I4
    s[, 6] <- occ_shift          # within DIF on I6
    matrix(rbinom(N * 8, 1, plogis(outer(th, d, "-") - s)), N, 8)
  }
  X <- rbind(gen(0), gen(0.9)); colnames(X) <- paste0("I", 1:8)
  dat <- data.frame(X, sex = rep(sex, 2), occasion = rep(c("t1", "t2"), each = N))
  id <- rep(sprintf("P%03d", 1:N), 2)
  fit <- rasch(dat, factors = c("sex", "occasion"), id = id)

  fa <- dif_anova(fit)
  expect_identical(fa$within, "occasion")
  s <- fa$summary
  # the between factor's DIF lands on the between item, the within factor's
  # on the within item, and nowhere else
  expect_true(s$uniform_DIF[s$item == "I4" & s$term == "sex"])
  expect_false(s$uniform_DIF[s$item == "I4" & s$term == "occasion"])
  expect_true(s$uniform_DIF[s$item == "I6" & s$term == "occasion"])
  expect_false(s$uniform_DIF[s$item == "I6" & s$term == "sex"])
  expect_equal(sum(s$uniform_DIF), 2L)

  # forcing the between-subjects treatment reproduces the ordinary factorial
  fb <- dif_anova(fit, within = character(0))
  expect_length(fb$within, 0L)
  fc <- dif_anova(rasch(data.frame(X[1:N, ], sex = sex),
                                  factors = "sex"))
  expect_length(fc$within, 0L)
})

test_that("resolve_dif splits DIF items by effect size and protects anchors", {
  # one strong real DIF item: resolve_dif splits it and ends with no DIF
  set.seed(3); n <- 800
  d <- seq(-2, 2, length.out = 10); g <- rep(c("a", "b"), each = n / 2)
  sh <- matrix(0, n, 10); sh[g == "b", 3] <- 1.4
  X <- matrix(rbinom(n * 10, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 10)
  colnames(X) <- sprintf("I%02d", 1:10)
  fit <- rasch(data.frame(X, grp = g), factors = "grp")
  rr <- resolve_dif(fit)
  expect_s3_class(rr, "rmt_resolve_dif")
  expect_true("I03" %in% rr$splits$item)          # the planted item is resolved
  expect_equal(rr$n_remaining_dif, 0L)            # nothing left
  expect_gt(ncol(rr$fit$X), ncol(fit$X))          # the fit gained resolved copies
  expect_lt(abs(rr$splits$magnitude[rr$splits$item == "I03"] - 1.4), 0.3)

  # pervasive, bidirectional DIF on most items must not resolve everything:
  # the anchor floor is preserved (the invariant), whatever the stop reason
  set.seed(5)
  shp <- matrix(0, n, 10)
  shp[g == "b", 1:4] <- 1.0; shp[g == "b", 5:8] <- -1.0
  Xp <- matrix(rbinom(n * 10, 1, plogis(outer(rnorm(n), d, "-") - shp)), n, 10)
  colnames(Xp) <- sprintf("I%02d", 1:10)
  fp <- rasch(data.frame(Xp, grp = g), factors = "grp")
  rp <- resolve_dif(fp, min_anchors = 3)
  # never fewer than min_anchors original items left unsplit
  expect_gte(10L - length(unique(rp$splits$item)), 3L)
})

test_that("mixed-design DIF survives missing data (strata dedup)", {
  # missing responses unbalance the design and aov then projects the within
  # terms onto the between-person stratum as well as their own; the flatten
  # step must keep each term once (its own stratum) or the planted DIF is
  # silently un-flagged and the summary duplicates rows
  set.seed(4); Np <- 400; L <- 6
  d0 <- scale(seq(-1.5, 1.5, length.out = L), scale = FALSE)[, 1]
  th <- rnorm(Np, 0, 1.3)
  sh <- matrix(0, Np, L); sh[, 3] <- 0.8
  X1 <- matrix(rbinom(Np * L, 1, plogis(outer(th, d0, "-"))), Np, L)
  X2 <- matrix(rbinom(Np * L, 1, plogis(outer(th, d0, "-") - sh)), Np, L)
  colnames(X1) <- colnames(X2) <- sprintf("W%02d", 1:L)
  stk <- rbind(data.frame(pid = sprintf("P%03d", 1:Np), X1, time = "t1"),
               data.frame(pid = sprintf("P%03d", 1:Np), X2, time = "t2"))
  set.seed(9)
  for (j in 1:L) stk[runif(nrow(stk)) < 0.08, 1 + j] <- NA
  fit <- rasch(stk, id = "pid", factors = "time")
  # aov warns "Error() model is singular" for some unbalanced items; the
  # flatten step is what is under test here
  r <- suppressWarnings(dif_anova(fit))
  # one row per item-term, and the planted within DIF is flagged
  expect_equal(nrow(r$summary),
               length(unique(r$summary$item)) * length(unique(r$summary$term)))
  expect_true(r$summary$uniform_DIF[r$summary$item == "W03"])
})

test_that("dif_anova tolerates factors named like its internals (ci, f1)", {
  set.seed(2); Np <- 600; L <- 6
  d0 <- scale(seq(-1.5, 1.5, length.out = L), scale = FALSE)[, 1]
  th <- rnorm(Np)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, d0, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  g <- rep(c("u", "v"), each = Np / 2)
  X[, 4] <- rbinom(Np, 1, plogis(th - d0[4] - ifelse(g == "v", 0.9, 0)))
  fit <- rasch(data.frame(X, ci = g, f1 = sample(c("p", "q"), Np, TRUE)),
               factors = c("ci", "f1"))
  r <- dif_anova(fit)
  expect_setequal(unique(r$summary$term), c("ci", "f1"))
  expect_true(r$summary$uniform_DIF[r$summary$item == "I04" &
                                    r$summary$term == "ci"])
  expect_equal(sum(r$summary$uniform_DIF[r$summary$term == "f1"]), 0L)
})

test_that("an item's own ci-crossing does not supersede it", {
  # uniform AND non-uniform DIF on the same item: both are reported on one
  # row, so the ci interaction must not knock the group term out of the
  # follow-ups (it used to, silently excluding the most DIF-affected item
  # from resolve_dif)
  set.seed(6); Np <- 800; L <- 6
  d0 <- scale(seq(-1.5, 1.5, length.out = L), scale = FALSE)[, 1]
  th <- rnorm(Np, 0, 1.4); g <- rep(c("r", "f"), each = Np / 2)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, d0, "-"))), Np, L)
  colnames(X) <- sprintf("H%02d", 1:L)
  X[, 5] <- rbinom(Np, 1, plogis(th - d0[5] -
                                 ifelse(g == "f", 0.9 + 0.5 * th, 0)))
  fit <- rasch(data.frame(X, grp = g), factors = "grp")
  r <- dif_anova(fit)
  s5 <- r$summary[r$summary$item == "H05", ]
  expect_true(s5$uniform_DIF)
  expect_true(s5$nonuniform_DIF)
  expect_false(s5$superseded)
})
