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


# person mode: responses to anchored items and judgements of the persons

cj_person_sim <- function(seed = 3, I = 12, N = 150, K = 600, alpha = 0.7,
                          n_rank = 0, kappa = 0.9, rank_size = 4, m = 1,
                          shift = NULL) {
  set.seed(seed)
  tau <- lapply(seq_len(I), function(i)
    sort(seq(-0.8, 0.8, length.out = m) + (i - (I + 1) / 2) * 3.6 / I))
  names(tau) <- sprintf("I%02d", seq_len(I))
  theta <- stats::rnorm(N)
  ids <- sprintf("S%03d", seq_len(N))
  X <- sapply(tau, function(t) {
    cs <- c(0, cumsum(t)); lp <- outer(theta, 0:m) - rep(cs, each = N)
    P <- exp(lp - apply(lp, 1, max)); P <- P / rowSums(P)
    apply(P, 1, function(p) sample(0:m, 1, prob = p))
  })
  rownames(X) <- ids
  th_j <- theta
  if (!is.null(shift)) th_j[match(names(shift), ids)] <- th_j[match(names(shift), ids)] + shift
  pa <- t(replicate(K, sample(N, 2)))
  p <- stats::plogis(alpha * (th_j[pa[, 1]] - th_j[pa[, 2]]))
  cj <- data.frame(a = ids[pa[, 1]], b = ids[pa[, 2]],
                   winner = ids[ifelse(stats::runif(K) < p, pa[, 1], pa[, 2])],
                   stringsAsFactors = FALSE)
  rk <- NULL
  if (n_rank > 0) {
    rk <- do.call(rbind, lapply(seq_len(n_rank), function(j) {
      s <- sample(N, rank_size)
      g <- -log(-log(stats::runif(rank_size)))
      data.frame(ranking = j, item = ids[s[order(-(kappa * th_j[s] + g))]],
                 rank = seq_len(rank_size), stringsAsFactors = FALSE)
    }))
  }
  anchors <- data.frame(item = rep(names(tau), each = m), k = rep(seq_len(m), I),
                        tau = unlist(tau), stringsAsFactors = FALSE)
  list(X = X, cj = cj, rk = rk, theta = theta, anchors = anchors, ids = ids)
}

test_that("the person mode recovers the persons, the units and its standard errors", {
  s <- cj_person_sim(n_rank = 80)
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner", rankings = s$rk, objects = "persons",
                  anchors = s$anchors)
  expect_s3_class(fit, "rasch_cj")
  expect_identical(fit$mode, "persons")
  expect_true(fit$converged)
  ps <- fit$persons
  expect_identical(ps$person, s$ids)
  expect_equal(fit$units$frame, c("responses", "comparisons", "rankings"))
  expect_true(abs(fit$units$unit[2] - 0.7) < 2.5 * fit$units$se[2])
  expect_true(abs(fit$units$unit[3] - 0.9) < 2.5 * fit$units$se[3])
  ok <- !ps$extreme & is.finite(ps$location_responses)
  rmse_joint <- sqrt(mean((ps$location - s$theta)[ok]^2))
  rmse_resp <- sqrt(mean((ps$location_responses - s$theta)[ok]^2))
  expect_lt(rmse_joint, rmse_resp)
  z <- (ps$location - s$theta)[!ps$extreme] / ps$se[!ps$extreme]
  expect_true(abs(stats::sd(z) - 1) < 0.15)
  # the response-only reference is the per-person maximum likelihood
  # location given the anchors
  j <- which(ok)[1]
  ml <- stats::optimize(function(t) rasch:::.cj_person_parts(t, s$X[j, , drop = FALSE],
                                                             split(s$anchors$tau, s$anchors$item)[colnames(s$X)])$ll,
                        c(-6, 6), maximum = TRUE)$maximum
  expect_equal(ps$location_responses[j], ml, tolerance = 1e-4)
  # the response block is the person likelihood: its gradient checks
  th <- c(-0.4, 1.1); Xs <- s$X[1:2, ]
  tl <- split(s$anchors$tau, s$anchors$item)[colnames(s$X)]
  p <- rasch:::.cj_person_parts(th, Xs, tl)
  expect_equal(p$g, num_grad(function(t) sum(rasch:::.cj_person_parts(t, Xs, tl)$ll), th),
               tolerance = 1e-4)
  inv <- fit$invariance$persons
  expect_null(fit$invariance$lr)
  expect_true(all(c("frame", "person", "reference", "judgements", "difference",
                    "se", "z", "p", "p_adj") %in% names(inv)))
  expect_lt(mean(inv$p < 0.05, na.rm = TRUE), 0.1)
  expect_output(print(fit), "Combined person measurement: 150 persons")
  expect_output(print(fit), "Unit of comparisons relative to responses")
  expect_output(print(fit), "Persons placed differently by the rankings")
  # fixed unit
  f1 <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                 winner = "winner", objects = "persons", anchors = s$anchors,
                 units = c(comparisons = 1))
  expect_equal(f1$units$unit[2], 1)
  expect_false(f1$units$estimated[2])
})

test_that("the person mode places extreme scorers and judged-only persons by the judgements", {
  s <- cj_person_sim(seed = 5, I = 10, N = 60, K = 300, alpha = 0.8)
  X <- s$X; X[1, ] <- 1L; X[2, ] <- 0L; X[3, ] <- 1L
  cj <- s$cj
  # S001 wins every comparison it is in: no downward evidence
  w1 <- cj$a == "S001" | cj$b == "S001"; cj$winner[w1] <- "S001"
  extra <- data.frame(a = c("Q1", "Q2", "S010"), b = c("S005", "Q1", "Q2"),
                      winner = c("Q1", "Q1", "S010"), stringsAsFactors = FALSE)
  rf <- rasch(X[4:60, ])
  fit <- rasch_cj(X, comparisons = rbind(cj, extra), object_a = "a", object_b = "b",
                  winner = "winner", objects = "persons", anchors = rf)
  ps <- fit$persons
  expect_equal(nrow(ps), 62)
  expect_true(ps$extreme[1])
  wle <- rasch:::.person_estimates(X[1, , drop = FALSE], rf$tau_list)
  expect_equal(ps$location[1], wle$theta)
  expect_true(is.infinite(ps$location_responses[1]))
  # S002 and S003 have extreme scores but wins and losses: finite joint locations
  expect_false(any(ps$extreme[2:3]))
  expect_true(all(is.finite(ps$location[2:3])))
  expect_lt(ps$location[2], ps$location[3])
  # Q1 lost to no one and Q2 won nothing: set aside with no location; they
  # have no responses to fall back on
  q <- match(c("Q1", "Q2"), ps$person)
  expect_true(all(ps$extreme[q]))
  expect_true(all(is.na(ps$location[q])))
  expect_true(all(is.na(ps$location_responses[q])))
  expect_true(any(grepl("judged person\\(s\\) without responses", fit$notes)))
  expect_true(any(grepl("extreme person\\(s\\) set aside", fit$notes)))
  # a judged-only person with evidence both ways is placed
  extra2 <- data.frame(a = c("Q1", "S010"), b = c("S005", "Q1"),
                       winner = c("Q1", "S010"), stringsAsFactors = FALSE)
  f2 <- rasch_cj(X, comparisons = rbind(cj, extra2), object_a = "a", object_b = "b",
                 winner = "winner", objects = "persons", anchors = rf)
  q1 <- match("Q1", f2$persons$person)
  expect_true(is.finite(f2$persons$location[q1]))
  expect_equal(f2$persons$n_items[q1], 0L)
  # a judged group with no link to a person with responses is refused
  bad <- rbind(cj, data.frame(a = c("Q3", "Q4"), b = c("Q4", "Q3"),
                              winner = c("Q3", "Q4"), stringsAsFactors = FALSE))
  expect_error(rasch_cj(X, comparisons = bad, object_a = "a", object_b = "b",
                        winner = "winner", objects = "persons", anchors = rf),
               "cannot be placed on the test scale: Q3, Q4")
})

test_that("the person mode takes anchors from a data frame, rasch() or rasch_cj(), and names ids", {
  s <- cj_person_sim(seed = 8, I = 8, N = 120, K = 400, m = 3)
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner", objects = "persons", anchors = s$anchors)
  expect_true(fit$converged)
  expect_equal(nrow(fit$anchors), 24)
  expect_true(all(fit$persons$max_raw == 24))
  # from a rasch_cj() item fit, an id column and an item subset
  cj_items <- rasch_cj(s$X, comparisons = data.frame(a = "I01", b = "I05", winner = "I05")[rep(1, 20), ],
                       object_a = "a", object_b = "b", winner = "winner")
  D <- data.frame(id = s$ids, s$X, check.names = FALSE)
  f2 <- rasch_cj(D, comparisons = s$cj, object_a = "a", object_b = "b",
                 winner = "winner", objects = "persons", anchors = cj_items,
                 id = "id", items = colnames(s$X)[1:6])
  expect_identical(f2$persons$person, s$ids)
  expect_true(any(grepl("2 anchored item\\(s\\) not in the data ignored", f2$notes)))
  expect_equal(f2$persons$n_items[1], 6L)
  # from a rasch() fit, with the id as a vector
  rf <- rasch(s$X)
  f3 <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                 winner = "winner", objects = "persons", anchors = rf, id = s$ids)
  expect_true(f3$converged)
  # refusals
  expect_error(rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                        winner = "winner", objects = "persons"),
               "needs `anchors`")
  expect_error(rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                        winner = "winner", anchors = rf),
               "`anchors` is for the person mode")
  expect_error(rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                        winner = "winner", objects = "persons", anchors = rf,
                        id = "nope"),
               "`id` must name a column")
  anc_bad <- s$anchors[s$anchors$k != 2, ]
  expect_error(rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                        winner = "winner", objects = "persons", anchors = anc_bad),
               "thresholds 1..m")
  Xd <- s$X[, 1:4]
  expect_error(rasch_cj(Xd, comparisons = s$cj, object_a = "a", object_b = "b",
                        winner = "winner", objects = "persons",
                        anchors = s$anchors[s$anchors$item != "I01", ]),
               "without anchored thresholds: I01")
  expect_error(rasch_cj(NULL, comparisons = s$cj, object_a = "a", object_b = "b",
                        winner = "winner", rankings = s$cj, objects = "persons",
                        anchors = rf),
               "needs response data")
})

test_that("the person mode flags a person whose work is judged unlike their responses", {
  s <- cj_person_sim(seed = 11, I = 8, N = 200, K = 800, alpha = 1.3, m = 3,
                     shift = c(S001 = 3))
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner", objects = "persons", anchors = s$anchors)
  inv <- fit$invariance$persons
  r <- inv$person == "S001"
  expect_gt(inv$difference[r], 1.5)
  expect_lt(inv$p[r], 0.05)
  expect_gt(fit$persons$location[1], fit$persons$location_responses[1])
})

test_that("the person mode reports a unit with no maximum when judgements are too sparse", {
  s <- cj_person_sim(seed = 2, I = 8, N = 100, K = 110, alpha = 0.8)
  fit <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                  winner = "winner", objects = "persons", anchors = s$anchors,
                  maxit = 60)
  expect_true(any(grepl("unit ran away", fit$notes)))
  f1 <- rasch_cj(s$X, comparisons = s$cj, object_a = "a", object_b = "b",
                 winner = "winner", objects = "persons", anchors = s$anchors,
                 units = c(comparisons = 1))
  expect_true(f1$converged)
  expect_true(all(is.finite(f1$persons$se[!f1$persons$extreme])))
})

# Two tests of one construct with no item in common, linked only by
# judgements. Test B is written in a unit rho times the reference: its
# item spread is rho times the true spread, and in the person mode its
# persons respond at rho theta - c.
cj_two_tests <- function(seed = 21, I = 10, N = 400, K = 800, alpha = 0.8,
                         rho = 1.5, shift = 0.6, persons = FALSE) {
  set.seed(seed)
  delta <- seq(-1.6, 1.6, length.out = I); delta <- delta - mean(delta)
  dA <- delta; names(dA) <- sprintf("A%02d", seq_len(I))
  dB <- delta; names(dB) <- sprintf("B%02d", seq_len(I))
  gen <- function(theta, d) {
    X <- sapply(d, function(v) as.integer(stats::runif(length(theta)) <
                                            stats::plogis(theta - v)))
    X
  }
  if (!persons) {
    XA <- gen(stats::rnorm(N), dA)
    XB <- gen(stats::rnorm(N), rho * dB)
    all <- c(dA, dB)
    pairs <- t(utils::combn(names(all), 2))[sample(choose(2 * I, 2), K, TRUE), ]
    p_a <- stats::plogis(alpha * (all[pairs[, 1]] - all[pairs[, 2]]))
    cj <- data.frame(a = pairs[, 1], b = pairs[, 2],
                     winner = ifelse(stats::runif(K) < p_a, pairs[, 1], pairs[, 2]),
                     stringsAsFactors = FALSE)
    return(list(XA = XA, XB = XB, cj = cj, dA = dA, dB = dB))
  }
  theta <- stats::rnorm(2 * N, 0, 1.2)
  ids <- sprintf("S%03d", seq_len(2 * N))
  XA <- gen(theta[seq_len(N)], dA); rownames(XA) <- ids[seq_len(N)]
  XB <- gen(rho * theta[N + seq_len(N)] - shift, dB); rownames(XB) <- ids[N + seq_len(N)]
  pa <- t(replicate(K, sample(2 * N, 2)))
  p <- stats::plogis(alpha * (theta[pa[, 1]] - theta[pa[, 2]]))
  cj <- data.frame(a = ids[pa[, 1]], b = ids[pa[, 2]],
                   winner = ids[ifelse(stats::runif(K) < p, pa[, 1], pa[, 2])],
                   stringsAsFactors = FALSE)
  anchors <- data.frame(item = c(names(dA), names(dB)), k = 1L, tau = c(dA, dB),
                        stringsAsFactors = FALSE)
  list(XA = XA, XB = XB, cj = cj, theta = theta, ids = ids, anchors = anchors)
}

test_that("two disconnected tests are calibrated together with a unit for the second", {
  s <- cj_two_tests()
  fit <- rasch_cj(list(A = s$XA, B = s$XB), comparisons = s$cj, object_a = "a",
                  object_b = "b", winner = "winner")
  expect_true(fit$converged)
  expect_identical(fit$tests, c("A", "B"))
  expect_identical(fit$reference, "A")
  expect_identical(fit$units$frame, c("A", "B", "comparisons"))
  uB <- fit$units[fit$units$frame == "B", ]
  expect_true(uB$estimated)
  expect_lt(abs(uB$unit - 1.5), 2.5 * uB$se)
  ua <- fit$units[fit$units$frame == "comparisons", ]
  expect_lt(abs(ua$unit - 0.8), 2.5 * ua$se)
  # items of both tests on one scale in the unit of A
  it <- fit$items
  expect_identical(it$item, c(names(s$dA), names(s$dB)))
  expect_lt(sqrt(mean((it$location - c(s$dA, s$dB))^2)), 0.2)
  expect_true(all(is.na(it$location_B[seq_len(10)])))
  expect_true(all(is.na(it$location_A[10 + seq_len(10)])))
  expect_true(all(is.finite(it$location_A[seq_len(10)])))
  # the separate calibration of B is in B's own unit; scaled by the unit
  # it matches the joint locations
  expect_lt(sqrt(mean((it$location_B[10 + 1:10] -
                         (it$location[10 + 1:10] - mean(it$location[10 + 1:10])))^2)), 0.15)
  # one free location per item in each frame (9 + 9 + 19), less the joint
  # parameters (19 locations and 2 units)
  expect_equal(fit$invariance$lr$df, 16L)
  expect_gt(fit$invariance$lr$p, 0.001)
  inv <- fit$invariance$items
  expect_setequal(unique(paste(inv$frame, inv$against)),
                  c("B comparisons", "comparisons A"))
  expect_equal(sum(inv$frame == "B"), 10L)
  expect_true(all(inv$p_adj > 0.05, na.rm = TRUE))
  expect_output(print(fit), "Unit of B relative to A")
  expect_output(print(fit), "Unit of comparisons relative to A")
})

test_that("a test's unit can be fixed and the tests are checked for consistency", {
  s <- cj_two_tests()
  fix <- rasch_cj(list(A = s$XA, B = s$XB), comparisons = s$cj, object_a = "a",
                  object_b = "b", winner = "winner", units = c(B = 1))
  expect_equal(fix$units$unit[fix$units$frame == "B"], 1)
  expect_false(fix$units$estimated[fix$units$frame == "B"])
  expect_equal(fix$invariance$lr$df, 17L)
  free <- rasch_cj(list(A = s$XA, B = s$XB), comparisons = s$cj, object_a = "a",
                   object_b = "b", winner = "winner")
  expect_gte(free$loglik, fix$loglik - 1e-8)
  expect_output(print(fix), "Unit of B relative to A: 1.000 \\(fixed\\)")
  # the reference test's unit is 1 by definition
  expect_error(rasch_cj(list(A = s$XA, B = s$XB), comparisons = s$cj, object_a = "a",
                        object_b = "b", winner = "winner", units = c(A = 1)),
               "reference")
  expect_error(rasch_cj(list(s$XA, s$XB), comparisons = s$cj, object_a = "a",
                        object_b = "b", winner = "winner"), "named")
  expect_error(rasch_cj(list(A = s$XA, comparisons = s$XB), comparisons = s$cj,
                        object_a = "a", object_b = "b", winner = "winner"),
               "comparisons")
  # an item in both tests must have the same categories
  XB2 <- cbind(s$XB, A01 = rep(0:2, length.out = nrow(s$XB)))
  expect_error(rasch_cj(list(A = s$XA, B = XB2), comparisons = s$cj, object_a = "a",
                        object_b = "b", winner = "winner"),
               "different categories across tests")
  # without a judgement between the tests there is no link
  cjA <- s$cj[s$cj$a %in% names(s$dA) & s$cj$b %in% names(s$dA), ]
  expect_error(rasch_cj(list(A = s$XA, B = s$XB), comparisons = cjA, object_a = "a",
                        object_b = "b", winner = "winner"), "connect|link|disconnected")
})

test_that("persons from two disconnected tests are measured on one scale", {
  s <- cj_two_tests(persons = TRUE, N = 150, K = 2000, alpha = 0.9, rho = 1.4)
  fit <- rasch_cj(list(A = s$XA, B = s$XB), comparisons = s$cj, object_a = "a",
                  object_b = "b", winner = "winner", objects = "persons",
                  anchors = s$anchors)
  expect_true(fit$converged)
  expect_identical(fit$mode, "persons")
  expect_identical(fit$reference, "A")
  expect_identical(fit$units$frame, c("A", "B", "comparisons"))
  tt <- fit$tests
  expect_identical(tt$test, c("A", "B"))
  expect_equal(tt$unit[1], 1); expect_equal(tt$shift[1], 0)
  # the unit carries the same upward bias as in the single-test person
  # mode, so the check is loose; the origin shift is recovered
  expect_lt(abs(tt$unit[2] - 1.4), 0.45)
  expect_lt(abs(tt$shift[2] - 0.6), 3 * tt$se_shift[2])
  expect_equal(tt$unit[2], fit$units$unit[2])
  ps <- fit$persons
  expect_identical(ps$person, s$ids)
  expect_identical(ps$test, rep(c("A", "B"), each = 150))
  expect_gt(stats::cor(ps$location, s$theta), 0.85)
  expect_lt(sqrt(mean((ps$location - s$theta)^2)), 0.75)
  # response-only locations of B are mapped to the common scale
  finite <- is.finite(ps$location_responses)
  expect_gt(stats::cor(ps$location_responses[finite], s$theta[finite]), 0.75)
  expect_true(all(c("test", "item", "k", "tau") %in% names(fit$anchors)))
  expect_output(print(fit), "Unit of B relative to A")
  expect_output(print(fit), "Unit of comparisons relative to A")
  # the anchors and ids can come per test
  fit2 <- rasch_cj(list(A = s$XA, B = s$XB), comparisons = s$cj, object_a = "a",
                   object_b = "b", winner = "winner", objects = "persons",
                   anchors = list(A = s$anchors[1:10, ], B = s$anchors[11:20, ]),
                   id = list(A = s$ids[1:150], B = s$ids[151:300]))
  expect_equal(fit2$loglik, fit$loglik)
  expect_error(rasch_cj(list(A = s$XA, B = s$XB), comparisons = s$cj, object_a = "a",
                        object_b = "b", winner = "winner", objects = "persons",
                        anchors = list(A = s$anchors[1:10, ])), "no element for test B")
  # a fixed unit leaves the shift free
  fix <- rasch_cj(list(A = s$XA, B = s$XB), comparisons = s$cj, object_a = "a",
                  object_b = "b", winner = "winner", objects = "persons",
                  anchors = s$anchors, units = c(B = 1))
  expect_equal(fix$tests$unit[2], 1)
  expect_true(is.na(fix$tests$se_unit[2]))
  expect_false(is.na(fix$tests$se_shift[2]))
  expect_gte(fit$loglik, fix$loglik - 1e-8)
})

test_that("the person mode refuses tests it cannot place", {
  s <- cj_two_tests(persons = TRUE, N = 60, K = 600, alpha = 0.9, rho = 1.4)
  idsA <- s$ids[1:60]; idsB <- s$ids[61:120]
  # no judgement across the tests
  cjA <- s$cj[s$cj$a %in% idsA & s$cj$b %in% idsA, ]
  expect_error(rasch_cj(list(A = s$XA, B = s$XB), comparisons = cjA, object_a = "a",
                        object_b = "b", winner = "winner", objects = "persons",
                        anchors = s$anchors), "no judgement links test B to test A")
  # one judged person in B: the origin is placed once the unit is fixed
  one <- idsB[which(rowSums(s$XB) > 0 & rowSums(s$XB) < 10)[1]]
  cj1 <- rbind(cjA, s$cj[(s$cj$a == one & s$cj$b %in% idsA) |
                           (s$cj$b == one & s$cj$a %in% idsA), ])
  expect_error(rasch_cj(list(A = s$XA, B = s$XB), comparisons = cj1, object_a = "a",
                        object_b = "b", winner = "winner", objects = "persons",
                        anchors = s$anchors), "units = c\\(B = 1\\)")
  fix <- rasch_cj(list(A = s$XA, B = s$XB), comparisons = cj1, object_a = "a",
                  object_b = "b", winner = "winner", objects = "persons",
                  anchors = s$anchors, units = c(B = 1))
  expect_true(fix$converged)
  expect_true(is.finite(fix$tests$shift[2]))
  # identifiers are unique across tests
  XB2 <- s$XB; rownames(XB2)[1] <- idsA[1]
  expect_error(rasch_cj(list(A = s$XA, B = XB2), comparisons = s$cj, object_a = "a",
                        object_b = "b", winner = "winner", objects = "persons",
                        anchors = s$anchors), "unique across tests: S001")
})
