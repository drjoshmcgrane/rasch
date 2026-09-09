pair_tie_fixture <- function() {
  objects <- c("A", "B", "C", "D")
  pairs <- utils::combn(objects, 2L, simplify = FALSE)
  rows <- lapply(seq_along(pairs), function(k) {
    pair <- pairs[[k]]
    key <- paste(pair, collapse = "")
    first_wins <- c(AB = 20L, AC = 10L, AD = 30L,
                    BC = 10L, BD = 30L, CD = 35L)[[key]]
    winner <- c(rep(pair[1L], first_wins),
                rep(pair[2L], 40L - first_wins))
    loser <- ifelse(winner == pair[1L], pair[2L], pair[1L])
    judge <- paste0("J", 2L + ((seq_len(40L) + 7L * k) %% 29L))
    if (key == "AB") judge[seq_len(10L)] <- "J1"
    if (key == "CD") judge[36:40] <- "J1"
    data.frame(a = winner, b = loser, winner = winner, judge = judge)
  })
  do.call(rbind, rows)
}

relabel_comparisons <- function(data, map) {
  out <- data
  for (column in c("a", "b", "winner"))
    out[[column]] <- unname(map[as.character(data[[column]])])
  out
}

restore_pair_labels <- function(pairs, map) {
  back <- stats::setNames(names(map), unname(map))
  pairs$object_hi <- unname(back[pairs$object_hi])
  pairs$object_lo <- unname(back[pairs$object_lo])
  pairs$net_winner <- unname(back[pairs$net_winner])
  key <- vapply(seq_len(nrow(pairs)), function(i)
    paste(sort(c(pairs$object_hi[i], pairs$object_lo[i])), collapse = "|"), "")
  pairs[order(key), , drop = FALSE]
}

select_pair <- function(pairs, a, b) {
  pair <- vapply(seq_len(nrow(pairs)), function(i)
    paste(sort(c(pairs$object_hi[i], pairs$object_lo[i])), collapse = "|"), "")
  pairs[pair == paste(sort(c(a, b)), collapse = "|"), , drop = FALSE]
}

test_that("pair surprise is label invariant when fitted locations tie", {
  data <- pair_tie_fixture()
  map <- c(A = "Z", B = "A", C = "M", D = "N")
  fit <- btl(data, "a", "b", "winner", judge = "judge")
  renamed <- btl(relabel_comparisons(data, map), "a", "b", "winner",
                 judge = "judge")

  original <- restore_pair_labels(
    judge_pair_surprise(fit, "J1", min_n = 5L)$pairs,
    stats::setNames(names(map), names(map)))
  relabelled <- restore_pair_labels(
    judge_pair_surprise(renamed, "J1", min_n = 5L)$pairs, map)

  expect_equal(fit$objects$location[fit$objects$object == "A"],
               fit$objects$location[fit$objects$object == "B"],
               tolerance = 1e-12)
  expect_true(any(original$tied))
  expect_true(any(!original$tied))
  expect_true(all(is.finite(original$p_adj)))
  expect_false(any(original$surprise[original$tied]))
  expect_equal(relabelled, original, tolerance = 1e-12)
  expect_equal(original$p_adj,
               p.adjust(original$p, method = "holm", n = nrow(original)))
})

test_that("pair tie tolerance is fixed and independent of the fitted origin", {
  data <- pair_tie_fixture()
  fit <- btl(data, "a", "b", "winner", judge = "judge")
  ab <- match(c("A", "B"), fit$objects$object)
  midpoint <- mean(fit$objects$location[ab])

  inside <- fit
  inside$objects$location[ab] <- midpoint + c(0.4e-10, -0.4e-10)
  translated <- inside
  translated$objects$location <- translated$objects$location + 500
  outside <- fit
  outside$objects$location[ab] <- midpoint + c(0.6e-10, -0.6e-10)

  row_inside <- select_pair(judge_pair_surprise(inside, "J1")$pairs, "A", "B")
  row_translated <- select_pair(
    judge_pair_surprise(translated, "J1")$pairs, "A", "B")
  row_outside <- select_pair(judge_pair_surprise(outside, "J1")$pairs, "A", "B")

  expect_true(row_inside$tied)
  expect_true(row_translated$tied)
  expect_false(row_outside$tied)
  expect_false(row_inside$surprise)
  expect_equal(row_translated[c("gap", "tied", "z", "p", "p_adj", "surprise")],
               row_inside[c("gap", "tied", "z", "p", "p_adj", "surprise")],
               tolerance = 1e-9)
})

test_that("pair surprise retains direction under non-tied relabeling", {
  data <- simulate_btl(
    4, 12, reps_per_pair = 16,
    object_locations = c(O1 = -1.5, O2 = -0.4, O3 = 0.5, O4 = 1.4),
    seed = 772)
  map <- c(O1 = "Z", O2 = "A", O3 = "Y", O4 = "B")
  fit <- btl(data, "object_a", "object_b", "winner", judge = "judge")
  renamed_data <- data
  names(renamed_data)[names(renamed_data) == "object_a"] <- "a"
  names(renamed_data)[names(renamed_data) == "object_b"] <- "b"
  renamed_data <- relabel_comparisons(renamed_data, map)
  renamed <- btl(renamed_data, "a", "b", "winner", judge = "judge")

  original <- restore_pair_labels(
    judge_pair_surprise(fit, "J1")$pairs,
    stats::setNames(names(map), names(map)))
  relabelled <- restore_pair_labels(
    judge_pair_surprise(renamed, "J1")$pairs, map)

  expect_false(any(original$tied))
  expect_equal(relabelled, original, tolerance = 1e-9)
})
