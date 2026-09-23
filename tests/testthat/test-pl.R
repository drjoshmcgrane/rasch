# Plackett-Luce rank analysis: derivatives, agreement with btl() on pairs,
# recovery, extremes, partial rankings, judges and the reversal check

pl_sim <- function(seed, R, size = 5, judges = NULL,
                   beta = c(A = -1.5, B = -0.5, C = 0, D = 0.5, E = 1.5),
                   worst_first = FALSE) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(R), function(r) {
    rem <- sample(names(beta), size); ord <- character(0)
    sgn <- if (worst_first) -1 else 1
    while (length(rem) > 1) {
      pick <- sample(rem, 1, prob = exp(sgn * beta[rem]))
      ord <- c(ord, pick); rem <- setdiff(rem, pick)
    }
    ord <- c(ord, rem)
    if (worst_first) ord <- rev(ord)
    data.frame(ranking = r, object = ord, rank = seq_len(size),
               judge = if (is.null(judges)) NA_character_ else
                 paste0("J", (r - 1) %% judges + 1),
               stringsAsFactors = FALSE)
  }))
}

test_that("the Plackett-Luce likelihood has the derivatives it claims", {
  rk <- pl_sim(1, 30)
  rk$rank[rk$rank > 3] <- NA          # a top-3 design
  o <- rasch:::.pl_cells(match(rk$object, LETTERS[1:5]), rk$rank, rk$ranking, 5)
  expect_equal(o$n_stages, 90L)
  b <- c(-1.4, -0.4, 0.1, 0.4, 1.3)
  ng <- sapply(1:5, function(k) { e <- b; e[k] <- e[k] + 1e-6
    (rasch:::.pl_ll(e, o) - rasch:::.pl_ll(b, o)) / 1e-6 })
  expect_equal(rasch:::.pl_grad(b, o, 5), ng, tolerance = 1e-4)
  nh <- sapply(1:5, function(k) { e <- b; e[k] <- e[k] + 1e-6
    (rasch:::.pl_grad(e, o, 5) - rasch:::.pl_grad(b, o, 5)) / 1e-6 })
  expect_equal(rasch:::.pl_hess(b, o, 5), nh, tolerance = 1e-3)
  # brute force: the probability of a ranking is the product of the stage choices
  one <- rk[rk$ranking == 1, ]
  ord <- one$object[order(one$rank, na.last = TRUE)]
  names(b) <- LETTERS[1:5]
  brute <- 0
  rem <- ord
  for (s in 1:3) { brute <- brute + b[ord[s]] - log(sum(exp(b[rem]))); rem <- setdiff(rem, ord[s]) }
  o1 <- rasch:::.pl_cells(match(one$object, LETTERS[1:5]), one$rank, one$ranking, 5)
  expect_equal(rasch:::.pl_ll(unname(b), o1), unname(brute))
})

test_that("rankings of two objects reproduce btl()", {
  pr <- pl_sim(2, 300, size = 2)
  f <- pl(pr)
  cmp <- data.frame(a = pr$object[pr$rank == 1], b = pr$object[pr$rank == 2])
  cmp$w <- cmp$a
  fb <- btl(cmp, "a", "b", "w")
  expect_equal(f$objects$location, fb$objects$location, tolerance = 1e-6)
  expect_equal(f$objects$se, fb$objects$se, tolerance = 1e-5)
  expect_null(f$reversal)
})

test_that("locations and their errors are recovered from rankings", {
  beta <- c(A = -1.5, B = -0.5, C = 0, D = 0.5, E = 1.5)
  rk <- pl_sim(3, 400)
  f <- pl(rk)
  expect_true(f$converged)
  expect_equal(f$objects$location, unname(beta), tolerance = 0.2)
  expect_true(all(abs(f$objects$fit_resid) < 3))
  expect_gt(f$osi$PSI, 0.9)
  expect_equal(sum(f$objects$location), 0)
  # error calibration over replications
  zs <- unlist(lapply(1:15, function(i) {
    g <- pl(pl_sim(100 + i, 100))
    (g$objects$location - beta) / g$objects$se
  }))
  expect_lt(abs(mean(zs)), 0.3)
  expect_lt(abs(stats::sd(zs) - 1), 0.3)
  # model-based errors are close to the sandwich ones under the model
  fm <- pl(rk, se = "model")
  expect_equal(fm$objects$se, f$objects$se, tolerance = 0.25)
  expect_equal(fm$objects$location, f$objects$location)
  # the reversal check does not reject best-first for best-first data
  expect_gt(f$reversal$z, -2)
  expect_gt(f$reversal$correlation, 0.98)
  expect_output(print(f), "sandwich SEs clustered by ranking")
  expect_output(print(fm), "information-based SEs")
  expect_silent(plot_pl(f))
})

test_that("the reversal check recognises worst-first rankings", {
  rw <- pl_sim(4, 300, worst_first = TRUE)
  f <- pl(rw)
  expect_lt(f$reversal$z, -2)
  expect_lt(f$reversal$p, 0.05)
  expect_equal(nrow(f$reversal$objects), 5L)
  expect_output(print(f), "Reversal check on 300 complete rankings")
})

test_that("judges are clustered, fitted and flagged", {
  rk <- pl_sim(5, 240, judges = 12)
  # one judge ranks at random
  rnd <- rk$judge == "J12"
  rk$rank[rnd] <- ave(rk$rank[rnd], rk$ranking[rnd], FUN = sample)
  f <- pl(rk, judge = "judge")
  expect_true(f$clustered)
  expect_equal(nrow(f$judges), 12L)
  expect_equal(f$judges$rankings, rep(20L, 12))
  expect_equal(f$judges$stages, rep(80L, 12))
  fine <- f$judges$judge != "J12"
  expect_true(all(abs(f$judges$fit_resid[fine]) < 3))
  expect_gt(f$judges$outfit_ms[!fine], max(f$judges$outfit_ms[fine]))
  expect_gt(f$judges$fit_resid[!fine], 2.5)
  expect_gt(f$judges$mean_surprise_z[!fine], max(f$judges$mean_surprise_z[fine]))
  expect_output(print(f), "clustered by judge")
  expect_output(print(f), "J12")
  expect_equal(f$rankings$judge[1:2], c("J1", "J2"))
  # too few judges: sandwich withheld, model errors available
  few <- rk[rk$judge %in% paste0("J", 1:5), ]
  g <- pl(few, judge = "judge")
  expect_true(all(is.na(g$objects$se)))
  expect_true(any(grepl("withheld", g$notes)))
  expect_true(is.na(g$osi$PSI))
  gm <- pl(few, judge = "judge", se = "model")
  expect_true(all(is.finite(gm$objects$se)))
  expect_error(pl(transform(rk, judge = ifelse(seq_len(nrow(rk)) == 1, "X", judge)),
                  judge = "judge"), "more than one judge")
})

test_that("extreme objects are set aside and placed by extrapolation", {
  rk <- pl_sim(6, 60)
  # Z always wins, Y always loses, in the rankings that hold them
  add <- rk$ranking <= 30
  top <- rk[add & rk$rank == 1, ]; top$object <- "Z"; top$rank <- 0
  bot <- rk[add & rk$rank == 1, ]; bot$object <- "Y"; bot$rank <- 9
  f <- pl(rbind(rk, top, bot))
  expect_true(f$converged)
  ext <- f$objects[f$objects$extreme, ]
  expect_equal(ext$object, c("Y", "Z"))
  expect_true(all(is.na(ext$se)))
  expect_gt(ext$location[ext$object == "Z"], max(f$objects$location[!f$objects$extreme]))
  expect_lt(ext$location[ext$object == "Y"], min(f$objects$location[!f$objects$extreme]))
  expect_equal(f$objects$rankings[f$objects$object == "Z"], 30L)
  expect_true(any(grepl("extreme", f$notes)))
  # the calibrated objects are unchanged by the extremes
  g <- pl(rk)
  expect_equal(f$objects$location[!f$objects$extreme], g$objects$location,
               tolerance = 1e-6)
  expect_output(print(f), "\\*")
})

test_that("partial rankings and ties are handled", {
  rk <- pl_sim(7, 200)
  # top-2 of five
  tk <- rk; tk$rank[tk$rank > 2] <- NA
  f <- pl(tk)
  expect_true(f$converged)
  expect_equal(f$n_stages, 400L)
  expect_equal(f$rankings$n_ranked, rep(2L, 200))
  expect_null(f$reversal)
  expect_equal(f$objects$location, c(-1.5, -0.5, 0, 0.5, 1.5), tolerance = 0.35)
  # rankings over subsets of the objects
  sub <- rk[!(rk$ranking %% 2 == 0 & rk$object == "C"), ]
  sub$rank <- ave(sub$rank, sub$ranking, FUN = rank)
  g <- pl(sub)
  expect_equal(g$size, c(4L, 5L))
  expect_equal(g$objects$rankings[3], 100L)
  # ties
  tie <- rk; tie$rank[tie$ranking == 1 & tie$rank == 2] <- 1
  expect_error(pl(tie, ties = "error"), "tied ranks")
  h <- pl(tie)
  expect_equal(h$n_rankings, 199L)
  expect_true(any(grepl("tied ranks dropped", h$notes)))
  # an object listed twice
  dup <- rbind(rk, rk[1, ])
  expect_error(pl(dup), "more than once")
  # a disconnected design: A-B-C never ranked below D-E
  dis <- rk
  dis$rank <- ave(seq_len(nrow(dis)), dis$ranking, FUN = function(i)
    rank(10 * (dis$object[i] %in% c("D", "E")) + dis$rank[i]))
  expect_error(pl(dis), "never ranked below")
  # anchors
  a <- pl(rk, anchors = c(A = -1.5, E = 1.5))
  expect_equal(a$objects$location[c(1, 5)], c(-1.5, 1.5))
  expect_equal(a$objects$se[c(1, 5)], c(0, 0))
  expect_true(all(a$objects$se[2:4] > 0))
  expect_error(pl(rk, anchors = c(Q = 1)), "do not match")
})
