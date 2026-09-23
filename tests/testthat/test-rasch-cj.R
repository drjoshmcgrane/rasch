# rasch_cj(): item responses and comparative judgements in one likelihood

cj_sim <- function(seed = 11, I = 10, N = 400, K = 500, alpha = 0.5,
                   shift = NULL, n_rank = 0, kappa = 0.8, rank_size = 4) {
  set.seed(seed)
  delta <- seq(-1.8, 1.8, length.out = I); delta <- delta - mean(delta)
  names(delta) <- sprintf("I%02d", seq_len(I))
  theta <- stats::rnorm(N)
  X <- sapply(delta, function(d) as.integer(stats::runif(N) < stats::plogis(theta - d)))
  d_j <- delta
  if (!is.null(shift)) d_j[names(shift)] <- d_j[names(shift)] + shift
  pairs <- t(utils::combn(names(delta), 2))[sample(choose(I, 2), K, TRUE), ]
  p_a <- stats::plogis(alpha * (d_j[pairs[, 1]] - d_j[pairs[, 2]]))
  cj <- data.frame(a = pairs[, 1], b = pairs[, 2],
                   winner = ifelse(stats::runif(K) < p_a, pairs[, 1], pairs[, 2]),
                   stringsAsFactors = FALSE)
  rk <- NULL
  if (n_rank > 0) {
    rk <- do.call(rbind, lapply(seq_len(n_rank), function(j) {
      rem <- sample(names(delta), rank_size); ord <- character(0)
      while (length(rem) > 1) {
        pick <- sample(rem, 1, prob = exp(kappa * d_j[rem]))
        ord <- c(ord, pick); rem <- setdiff(rem, pick)
      }
      data.frame(ranking = j, item = c(ord, rem), rank = seq_len(rank_size),
                 stringsAsFactors = FALSE)
    }))
  }
  list(X = X, cj = cj, rk = rk, delta = delta)
}

num_grad <- function(f, x, h = 1e-6)
  vapply(seq_along(x), function(k) { e <- x; e[k] <- e[k] + h; (f(e) - f(x)) / h }, 0)

test_that("each likelihood block has the gradient and Hessian it claims", {
  set.seed(4)
  I <- 6; delta <- seq(-1, 1, length.out = I)
  th <- c(delta, log(0.7))
  W <- matrix(stats::rpois(I * I, 5), I, I); diag(W) <- 0
  p <- rasch:::.cj_bt_parts(delta, 0.7, W)
  f <- function(t) rasch:::.cj_bt_ll(t[1:I], exp(t[I + 1]), W)
  expect_equal(p$g, num_grad(f, th), tolerance = 1e-4)
  Hn <- sapply(seq_along(th), function(k) { e <- th; e[k] <- e[k] + 1e-5
    (rasch:::.cj_bt_parts(e[1:I], exp(e[I + 1]), W)$g - p$g) / 1e-5 })
  expect_equal(p$H, Hn, tolerance = 1e-3)

  rk <- list(c(1L, 3L, 5L, 2L), c(6L, 4L, 1L), c(2L, 6L))
  p <- rasch:::.cj_pl_parts(delta, 0.7, rk)
  f <- function(t) rasch:::.cj_pl_ll(t[1:I], exp(t[I + 1]), rk)
  expect_equal(p$g, num_grad(f, th), tolerance = 1e-4)
  Hn <- sapply(seq_along(th), function(k) { e <- th; e[k] <- e[k] + 1e-5
    (rasch:::.cj_pl_parts(e[1:I], exp(e[I + 1]), rk)$g - p$g) / 1e-5 })
  expect_equal(p$H, Hn, tolerance = 1e-3)

  N <- 150; theta <- stats::rnorm(N)
  X <- sapply(delta, function(d) as.integer(stats::runif(N) < stats::plogis(theta - d)))
  X[sample(length(X), 30)] <- NA
  o <- rasch:::.cj_cml_prep(X)
  d <- delta + 0.1
  expect_equal(rasch:::.cj_cml_grad(d, o),
               num_grad(function(z) rasch:::.cj_cml_ll(z, o), d), tolerance = 1e-4)
  Hn <- sapply(seq_along(d), function(k) { e <- d; e[k] <- e[k] + 1e-5
    (rasch:::.cj_cml_grad(e, o) - rasch:::.cj_cml_grad(d, o)) / 1e-5 })
  expect_equal(rasch:::.cj_cml_hess(d, o), Hn, tolerance = 1e-3)
})

test_that("the pattern-collapsed CML equals the brute-force conditional likelihood", {
  set.seed(5)
  I <- 5; N <- 60; delta <- c(-1, -0.5, 0, 0.5, 1)
  X <- sapply(delta, function(d) as.integer(stats::runif(N) < stats::plogis(stats::rnorm(N) - d)))
  o <- rasch:::.cj_cml_prep(X)
  d <- delta + 0.2
  brute <- 0
  for (n in seq_len(N)) {
    r <- sum(X[n, ]); if (r == 0 || r == I) next
    pats <- as.matrix(expand.grid(rep(list(0:1), I)))
    pats <- pats[rowSums(pats) == r, , drop = FALSE]
    brute <- brute + (-sum(X[n, ] * d)) - log(sum(exp(-pats %*% d)))
  }
  expect_equal(rasch:::.cj_cml_ll(d, o), brute, tolerance = 1e-10)
})

test_that("locations and the unit are recovered from responses and comparisons", {
  s <- cj_sim()
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner")
  expect_s3_class(fit, "rasch_cj")
  expect_true(fit$converged)
  expect_equal(fit$items$item, names(s$delta))
  expect_lt(max(abs(fit$items$location - s$delta)), 0.35)
  expect_lt(mean(abs(fit$items$location - s$delta)), 0.15)
  u <- fit$units[fit$units$frame == "comparisons", ]
  expect_true(u$estimated)
  expect_lt(abs(u$unit - 0.5), 2.5 * u$se)
  expect_equal(sum(fit$items$location), 0, tolerance = 1e-8)
  expect_true(all(fit$items$se > 0))
  expect_gt(fit$invariance$lr$p, 0.01)
  expect_equal(fit$invariance$lr$df, 8L)
  expect_true(all(fit$invariance$items$p_adj > 0.05))
  expect_output(print(fit), "Unit of comparisons relative to responses")
})

test_that("the combined fit is more precise than the response frame alone", {
  s <- cj_sim(seed = 12, K = 2000)
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner")
  alone <- rasch:::.cj_fit_alone(
    function(d) rasch:::.cj_cml_ll(d, rasch:::.cj_cml_prep(s$X)),
    function(d) { o <- rasch:::.cj_cml_prep(s$X); pp <- rasch:::.cj_cml_p(d, o)
      list(g = rasch:::.cj_cml_grad(d, o, pp), H = rasch:::.cj_cml_hess(d, o, pp)) },
    rasch:::.cj_basis(ncol(s$X)))
  expect_true(all(fit$items$se < alone$se))
  expect_equal(fit$items$location_responses, alone$par, tolerance = 1e-6)
})

test_that("rankings of two items reproduce the comparison fit", {
  s <- cj_sim(seed = 13, K = 300)
  rk <- data.frame(ranking = rep(seq_len(nrow(s$cj)), each = 2),
                   item = as.vector(t(cbind(s$cj$winner,
                     ifelse(s$cj$winner == s$cj$a, s$cj$b, s$cj$a)))),
                   rank = rep(1:2, nrow(s$cj)), stringsAsFactors = FALSE)
  f1 <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                 winner = "winner")
  f2 <- rasch_cj(s$X, rankings = rk)
  expect_equal(f2$items$location, f1$items$location, tolerance = 1e-6)
  expect_equal(f2$units$unit[2], f1$units$unit[2], tolerance = 1e-6)
  expect_equal(f2$loglik, f1$loglik, tolerance = 1e-8)
  expect_equal(f2$units$frame[2], "rankings")
})

test_that("three frames fit together and each unit is recovered", {
  s <- cj_sim(seed = 14, n_rank = 200, kappa = 0.8)
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner", rankings = s$rk)
  expect_true(fit$converged)
  expect_equal(fit$units$frame, c("responses", "comparisons", "rankings"))
  expect_lt(abs(fit$units$unit[2] - 0.5), 2.5 * fit$units$se[2])
  expect_lt(abs(fit$units$unit[3] - 0.8), 2.5 * fit$units$se[3])
  expect_equal(fit$invariance$lr$df, 2L * 9L - 2L)
  expect_equal(unique(fit$invariance$items$frame), c("comparisons", "rankings"))
  expect_lt(max(abs(fit$items$location - s$delta)), 0.35)
})

test_that("a fixed unit is honoured and changes the test's degrees of freedom", {
  s <- cj_sim(seed = 15, alpha = 1)
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner", units = c(comparisons = 1))
  expect_equal(fit$units$unit[2], 1)
  expect_false(fit$units$estimated[2])
  expect_true(is.na(fit$units$se[2]))
  expect_equal(fit$invariance$lr$df, 9L)
  expect_output(print(fit), "1.000 \\(fixed\\)")
  free <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                   winner = "winner")
  expect_gte(free$loglik, fit$loglik - 1e-8)
  expect_lt(abs(free$units$unit[2] - 1), 2.5 * free$units$se[2])
})

test_that("an item the judges see differently fails the invariance test", {
  s <- cj_sim(seed = 16, K = 1500, shift = c(I03 = 2.5))
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner")
  expect_lt(fit$invariance$lr$p, 0.001)
  tab <- fit$invariance$items
  expect_lt(tab$p_adj[tab$item == "I03"], 0.05)
  expect_gt(tab$judgements[tab$item == "I03"] - tab$reference[tab$item == "I03"], 1)
  expect_output(print(fit), "I03")
})

test_that("an item with no response variation is located by the judgements", {
  s <- cj_sim(seed = 17, K = 800)
  X <- s$X; X[, "I01"] <- 1L
  # the judgement data were generated with I01 easiest, consistent with a
  # column every person answered correctly
  fit <- rasch_cj(X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner")
  expect_true(fit$converged)
  expect_true(any(grepl("no response variation", fit$notes)))
  expect_equal(fit$items$item, colnames(X))
  loc <- fit$items$location
  expect_equal(which.min(loc), 1L)
  expect_true(is.finite(fit$items$se[1]))
  expect_error(rasch_cj(X, comparisons = s$cj[s$cj$a != "I01" & s$cj$b != "I01", ],
                        object_a = "a", object_b = "b", winner = "winner"),
               "no response variation and no judgements")
})

test_that("the anchor table carries the calibration into rasch()", {
  s <- cj_sim(seed = 18)
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner")
  expect_equal(names(fit$anchors), c("item", "k", "tau"))
  r <- rasch(s$X, anchors = fit$anchors[-1, ])
  expect_equal(r$items$location[-1], fit$items$location[-1], tolerance = 1e-6)
  expect_equal(r$items$se[-1], rep(0, ncol(s$X) - 1L))
})

test_that("ties are dropped with a note and reversed orientation is flagged", {
  s <- cj_sim(seed = 19)
  cj <- s$cj; cj$winner[1:5] <- "tie"
  fit <- rasch_cj(s$X, comparisons = cj, object_a = "a", object_b = "b",
                  winner = "winner")
  expect_true(any(grepl("5 tied comparison", fit$notes)))
  expect_equal(fit$n[["comparisons"]], nrow(s$cj) - 5L)
  rev <- s$cj; rev$winner <- ifelse(rev$winner == rev$a, rev$b, rev$a)
  fit <- rasch_cj(s$X, comparisons = rev, object_a = "a", object_b = "b",
                  winner = "winner")
  expect_true(any(grepl("orientation is reversed", fit$notes)))
})

test_that("input errors are specific", {
  s <- cj_sim(seed = 20)
  expect_error(rasch_cj(s$X), "supply `comparisons`, `rankings`, or both")
  bad <- s$cj; bad$a[1] <- "Z9"
  expect_error(rasch_cj(s$X, comparisons = bad, object_a = "a", object_b = "b",
                        winner = "winner"), "not items: 'Z9'")
  bad <- s$cj; bad$winner[1] <- "I07"; bad$a[1] <- "I01"; bad$b[1] <- "I02"
  expect_error(rasch_cj(s$X, comparisons = bad, object_a = "a", object_b = "b",
                        winner = "winner"), "neither compared object")
  expect_error(rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                        winner = "result"), "lacks column")
  expect_error(rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                        winner = "winner", units = c(comparisons = 2)),
               "fixed unit must be 1")
  expect_error(rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                        winner = "winner", units = c(other = NA)), "names must be")
  # a gap in the categories is rescored, as rasch() does, and noted
  poly <- s$X; poly[1, 1] <- 3L
  gap <- rasch_cj(poly, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner")
  expect_true(any(grepl("I01 rescored", gap$notes)))
  expect_equal(nrow(gap$thresholds), ncol(s$X) + 1L)
  rk <- data.frame(ranking = c(1, 1, 2, 2), item = c("I01", "I01", "I02", "I03"),
                   rank = c(1, 2, 1, 2))
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner", rankings = rk)
  expect_true(any(grepl("1 ranking\\(s\\) dropped", fit$notes)))
  expect_equal(fit$n[["rankings"]], 1L)
  rk_only <- data.frame(ranking = 1, item = c("I01", "I01"), rank = 1:2)
  expect_error(rasch_cj(s$X, rankings = rk_only), "no usable rows")
})

test_that("a disconnected design is refused", {
  s <- cj_sim(seed = 21, I = 6, K = 200)
  X <- s$X
  X[1:200, 4:6] <- NA; X[201:400, 1:3] <- NA
  within <- with(s$cj, (a %in% c("I01", "I02", "I03")) == (b %in% c("I01", "I02", "I03")))
  expect_error(rasch_cj(X, comparisons = s$cj[within, ], object_a = "a",
                        object_b = "b", winner = "winner"), "not connected")
  fit <- rasch_cj(X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner")
  expect_true(fit$converged)
})

test_that("comparisons and rankings combine without responses", {
  s <- cj_sim(11, I = 8, K = 800, alpha = 0.5, n_rank = 300, kappa = 1.5)
  fit <- rasch_cj(NULL, comparisons = s$cj, object_a = "a", object_b = "b",
                  rankings = s$rk)
  expect_true(fit$converged)
  expect_equal(fit$reference, "comparisons")
  expect_equal(fit$n[["persons"]], 0L)
  # the rankings' unit is relative to the comparisons' scale
  expect_equal(fit$units$unit[fit$units$frame == "comparisons"], 1)
  expect_true(fit$units$estimated[fit$units$frame == "rankings"])
  expect_equal(fit$units$unit[fit$units$frame == "rankings"], 1.5 / 0.5,
               tolerance = 0.25)
  expect_true(cor(fit$items$location, s$delta) > 0.95)
  expect_equal(fit$invariance$lr$df, 6L)
  expect_output(print(fit), "Unit of rankings relative to comparisons")
  expect_true("location_comparisons" %in% names(fit$items))
  expect_error(rasch_cj(NULL, comparisons = s$cj, object_a = "a", object_b = "b"),
               "supply both")
  expect_error(rasch_cj(NULL, rankings = s$rk), "supply both")
})

# polytomous simulation: partial credit items with thresholds psi[[i]],
# judgements of item locations (mean threshold) or of single thresholds
pcm_sim <- function(seed, N = 500, m = c(1, 2, 3, 2, 3, 1), K = 800, alpha = 0.6,
                    thr_cj = FALSE, n_rank = 0, kappa = 1, rank_size = 4,
                    shift = NULL) {
  set.seed(seed)
  I <- length(m); items <- sprintf("I%02d", seq_len(I))
  psi <- lapply(seq_len(I), function(i) {
    centre <- (i - (I + 1) / 2) * 0.5
    if (m[i] == 1) centre else centre + seq(-0.8, 0.8, length.out = m[i])
  })
  psi_all <- unlist(psi); w <- rep(1 / m, m)
  psi_all <- psi_all - sum(w * psi_all) / sum(w)      # mean item location zero
  psi <- split(psi_all, rep(seq_len(I), m))
  theta <- stats::rnorm(N)
  X <- sapply(seq_len(I), function(i) {
    cs <- c(0, cumsum(psi[[i]]))
    vapply(theta, function(t) {
      p <- exp((0:m[i]) * t - cs); sample(0:m[i], 1, prob = p / sum(p))
    }, 0L)
  })
  colnames(X) <- items
  # objects: items, or every threshold
  if (thr_cj) {
    obj <- data.frame(item = rep(items, m), k = unlist(lapply(m, seq_len)),
                      loc = psi_all, stringsAsFactors = FALSE)
  } else {
    obj <- data.frame(item = items, k = NA, loc = vapply(psi, mean, 0),
                      stringsAsFactors = FALSE)
  }
  # shift = c(key, amount): the judges see one object elsewhere
  key <- ifelse(is.na(obj$k), obj$item, paste0(obj$item, ":", obj$k))
  if (!is.null(shift)) obj$loc[key == shift[1]] <- obj$loc[key == shift[1]] + as.numeric(shift[2])
  pairs <- t(utils::combn(nrow(obj), 2))[sample(choose(nrow(obj), 2), K, TRUE), ]
  p_a <- stats::plogis(alpha * (obj$loc[pairs[, 1]] - obj$loc[pairs[, 2]]))
  win_a <- stats::runif(K) < p_a
  cj <- data.frame(a = obj$item[pairs[, 1]], b = obj$item[pairs[, 2]],
                   threshold_a = obj$k[pairs[, 1]], threshold_b = obj$k[pairs[, 2]],
                   stringsAsFactors = FALSE)
  cj$winner <- ifelse(win_a, key[pairs[, 1]], key[pairs[, 2]])
  rk <- NULL
  if (n_rank > 0) {
    rk <- do.call(rbind, lapply(seq_len(n_rank), function(j) {
      rem <- sample(nrow(obj), rank_size); ord <- integer(0)
      while (length(rem) > 1) {
        pick <- sample(rem, 1, prob = exp(kappa * obj$loc[rem]))
        ord <- c(ord, pick); rem <- setdiff(rem, pick)
      }
      o <- c(ord, rem)
      data.frame(ranking = j, item = obj$item[o], threshold = obj$k[o],
                 rank = seq_len(rank_size), stringsAsFactors = FALSE)
    }))
  }
  list(X = X, cj = cj, rk = rk, psi = psi_all, m = m, obj = obj)
}

test_that("the polytomous CML block has the derivatives it claims", {
  s <- pcm_sim(21, N = 120)
  X <- s$X; X[sample(length(X), 40)] <- NA
  o <- rasch:::.cj_cml_prep(X)
  expect_equal(o$m, s$m)
  d <- s$psi + 0.1
  expect_equal(rasch:::.cj_cml_grad(d, o),
               num_grad(function(z) rasch:::.cj_cml_ll(z, o), d), tolerance = 1e-4)
  Hn <- sapply(seq_along(d), function(k) { e <- d; e[k] <- e[k] + 1e-5
    (rasch:::.cj_cml_grad(e, o) - rasch:::.cj_cml_grad(d, o)) / 1e-5 })
  expect_equal(rasch:::.cj_cml_hess(d, o), Hn, tolerance = 1e-3)
})

test_that("the polytomous CML equals the brute-force conditional likelihood", {
  set.seed(22)
  m <- c(2, 1, 3, 2); N <- 40
  s <- pcm_sim(22, N = N, m = m)
  o <- rasch:::.cj_cml_prep(s$X, m)
  psi <- s$psi + 0.2
  idx <- split(seq_along(psi), rep(seq_along(m), m))
  term <- function(pat) sum(vapply(seq_along(m), function(i)
    -sum(psi[idx[[i]]][seq_len(pat[i])]), 0))
  pats <- as.matrix(expand.grid(lapply(m, function(k) 0:k)))
  brute <- 0
  for (n in seq_len(N)) {
    r <- sum(s$X[n, ]); if (r == 0 || r == sum(m)) next
    same <- pats[rowSums(pats) == r, , drop = FALSE]
    brute <- brute + term(s$X[n, ]) - log(sum(exp(apply(same, 1, term))))
  }
  expect_equal(rasch:::.cj_cml_ll(psi, o), brute, tolerance = 1e-10)
})

test_that("polytomous items are located by judgements of their mean threshold", {
  s <- pcm_sim(23, N = 600, K = 1500)
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner")
  expect_true(fit$converged)
  expect_equal(nrow(fit$thresholds), sum(s$m))
  expect_equal(fit$thresholds$k, unlist(lapply(s$m, seq_len)))
  expect_true(all(is.na(fit$objects$threshold)))
  expect_equal(fit$thresholds$threshold, s$psi, tolerance = 0.35)
  expect_equal(fit$units$unit[2], 0.6, tolerance = 0.25)
  expect_gt(fit$invariance$lr$p, 0.01)
  # the response frame alone agrees with rasch()'s pairwise calibration
  r <- rasch(s$X)
  expect_equal(fit$items$location_responses, r$items$location, tolerance = 0.1)
  # anchoring every item but one reproduces the combined thresholds
  a <- rasch(s$X, anchors = fit$anchors[fit$anchors$item != "I01", ])
  expect_equal(a$items$location[-1], fit$items$location[-1], tolerance = 1e-6)
  expect_output(print(fit), "12 thresholds")
})

test_that("judgements at the threshold level inform single thresholds", {
  s <- pcm_sim(24, N = 400, K = 2500, thr_cj = TRUE, n_rank = 400, kappa = 1.2,
               rank_size = 5)
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner", rankings = s$rk)
  expect_true(fit$converged)
  expect_equal(nrow(fit$objects), sum(s$m))
  expect_true(all(!is.na(fit$objects$threshold)))
  expect_equal(fit$thresholds$threshold, s$psi, tolerance = 0.3)
  expect_equal(fit$units$unit[2:3], c(0.6, 1.2), tolerance = 0.3)
  # every threshold is judged in each frame: (P - 1) free per frame, 2 units
  expect_equal(fit$invariance$lr$df, 2L * (sum(s$m) - 1L) - 2L)
  tab <- fit$invariance$items
  expect_equal(nrow(tab), 2L * sum(s$m))
  expect_true(all(!is.na(tab$threshold)))
  expect_true(all(tab$p_adj > 0.001, na.rm = TRUE))
  # a threshold the judges see differently is flagged, and thresholds of one
  # item can be compared with each other when the winner names the key
  s2 <- pcm_sim(24, N = 400, K = 2500, thr_cj = TRUE, shift = c("I03:3", 2.5))
  fit2 <- rasch_cj(s$X, comparisons = s2$cj, object_a = "a", object_b = "b",
                   winner = "winner")
  tab2 <- fit2$invariance$items
  expect_lt(tab2$p_adj[tab2$item == "I03" & tab2$threshold == 3], 0.05)
  expect_output(print(fit2), "I03:3")
  # a threshold the item does not have, and thresholds without responses
  bad <- s$cj; bad$threshold_a[1] <- 4
  expect_error(rasch_cj(s$X, comparisons = bad, object_a = "a", object_b = "b",
                        winner = "winner"), "thresholds an item does not have")
  expect_error(rasch_cj(NULL, comparisons = s$cj, object_a = "a", object_b = "b",
                        winner = "winner", rankings = s$rk),
               "need response data")
})

test_that("a frame judging a subset of items is tested on those alone", {
  s <- cj_sim(seed = 25, K = 600)
  sub <- s$cj[s$cj$a %in% sprintf("I%02d", 1:5) & s$cj$b %in% sprintf("I%02d", 1:5), ]
  fit <- rasch_cj(s$X, comparisons = sub, object_a = "a", object_b = "b",
                  winner = "winner")
  expect_true(fit$converged)
  expect_equal(sum(fit$objects$comparisons), 5L)
  expect_equal(fit$invariance$lr$df, 4L - 1L)
  expect_equal(nrow(fit$invariance$items), 5L)
  expect_true(all(is.na(fit$items$location_comparisons[6:10])))
  # two disconnected blocks of comparisons: one constraint per block
  blocks <- s$cj[(s$cj$a %in% sprintf("I%02d", 1:5)) == (s$cj$b %in% sprintf("I%02d", 1:5)), ]
  fit2 <- rasch_cj(s$X, comparisons = blocks, object_a = "a", object_b = "b",
                   winner = "winner")
  expect_true(fit2$converged)
  expect_equal(fit2$invariance$lr$df, 10L - 2L - 1L)
  expect_gt(fit2$invariance$lr$p, 0.001)
})

test_that("judgement identifiers are trimmed and mismatches are named", {
  X <- pcm_sim(5, N = 150, m = rep(1, 5), K = 20)$X
  cmp <- data.frame(object_a = c(" I01", "I02 "), object_b = c("I03", "I04"),
                    winner = c("I03 ", " I02"), stringsAsFactors = FALSE)
  fit <- rasch_cj(X, comparisons = cmp)
  expect_equal(fit$n[["comparisons"]], 2L)
  cmp$object_a[1] <- "i01"
  expect_error(rasch_cj(X, comparisons = cmp), "not items: 'i01'")
})

