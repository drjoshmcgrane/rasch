test_that("a third BTL fit cannot erase two fits' observed sequence mismatch", {
  d <- simulate_btl(6, 30, reps_per_pair = 80, seed = 21)
  d$order <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  changed <- d
  set.seed(33)
  changed$order <- ave(changed$order, changed$judge, FUN = sample)
  ref <- btl(d, "object_a", "object_b", "winner", judge = "judge",
             order = "order")
  other <- btl(changed, "object_a", "object_b", "winner", judge = "judge",
               order = "order")
  plain <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  cols <- c("same_data", "cl_aic", "cl_bic", "two_delta_ll")
  two <- compare_fits(ref = ref, other = other)
  three <- compare_fits(ref = ref, other = other, plain = plain)
  expect_false(two$same_data[2L])
  expect_equal(as.data.frame(two)[2L, cols], as.data.frame(three)[2L, cols])
  expect_true(all(is.na(unlist(three[2L, cols[-1L]]))))
  expect_true(three$same_data[3L])

  # Conversely, an ordered third model must not impose presentation
  # orientation on a comparison between two plain models.
  reversed <- d
  reversed$object_a <- d$object_b
  reversed$object_b <- d$object_a
  reverse_fit <- btl(reversed, "object_a", "object_b", "winner",
                     judge = "judge")
  two <- compare_fits(ref = plain, other = reverse_fit)
  three <- compare_fits(ref = plain, other = reverse_fit, ordered = ref)
  expect_true(two$same_data[2L])
  expect_equal(as.data.frame(two)[2L, cols], as.data.frame(three)[2L, cols])
})

test_that("BTL sequence identity ignores labels but retains ordering", {
  d <- simulate_btl(5, 18, reps_per_pair = 50, seed = 721)
  d$order <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  ref <- btl(d, "object_a", "object_b", "winner", judge = "judge",
             order = "order")

  # Relabel each judge's sequence monotonically, relabel the judges, and
  # shuffle rows.  None changes the histories seen by the likelihood.
  relabelled <- d
  judge_number <- match(relabelled$judge, unique(relabelled$judge))
  relabelled$order <- 1000 * judge_number + 7 * relabelled$order
  relabelled$judge <- paste0("anonymous-", 100 + judge_number)
  set.seed(722)
  relabelled <- relabelled[sample.int(nrow(relabelled)), , drop = FALSE]
  position <- btl(relabelled, "object_a", "object_b", "winner",
                  judge = "judge", order = "order", position = TRUE)

  out <- compare_fits(history = ref, position = position)
  expect_true(out$same_data[2L])
  expect_true(all(is.finite(unlist(out[2L, c("cl_aic", "cl_bic",
                                             "two_delta_ll")]))))

  # A decreasing relabelling changes which comparison is first and must
  # remain protected as genuinely different comparison data.
  reordered <- d
  reordered$order <- ave(reordered$order, reordered$judge,
                         FUN = function(x) max(x) + 1 - x)
  reordered_fit <- btl(reordered, "object_a", "object_b", "winner",
                       judge = "judge", order = "order")
  expect_false(compare_fits(history = ref,
                            reordered = reordered_fit)$same_data[2L])
})

test_that("BTL comparison signatures preserve sequence ties", {
  d <- simulate_btl(4, 12, reps_per_pair = 40, seed = 723)
  d$order <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)

  tied <- d
  first <- which(tied$judge == tied$judge[1L])[1:2]
  tied$order[first] <- tied$order[first[1L]]
  expect_error(
    btl(tied, "object_a", "object_b", "winner", judge = "judge",
        order = "order"),
    "repeats within judge"
  )

  # Legacy objects can predate that guard.  Equal tie patterns are invariant
  # to monotone labels, but are not confused with the untied sequence.
  fit <- btl(d, "object_a", "object_b", "winner", judge = "judge",
             order = "order")
  legacy_tied <- fit
  legacy_tied$comparisons$order[first] <-
    legacy_tied$comparisons$order[first[1L]]
  relabelled_tie <- legacy_tied
  relabelled_tie$comparisons$order <-
    13 + 11 * relabelled_tie$comparisons$order

  expect_true(compare_fits(original = legacy_tied,
                           relabelled = relabelled_tie)$same_data[2L])
  expect_false(compare_fits(original = fit,
                            tied = legacy_tied)$same_data[2L])
})
