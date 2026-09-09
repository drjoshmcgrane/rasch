test_that("within-set strong components follow indirect win paths", {
  cycle <- matrix(FALSE, 4, 4)
  cycle[cbind(1:4, c(2:4, 1))] <- TRUE
  expect_identical(.btlef_strong_components(cycle), rep(1L, 4))
  cycle[4, 1] <- FALSE
  expect_length(unique(.btlef_strong_components(cycle)), 4L)
  cycle[2, 1] <- TRUE
  comp <- .btlef_strong_components(cycle)
  expect_identical(comp[1], comp[2])
  expect_length(unique(comp), 3L)
})

test_that("BTL-EFRM reports within-set outcome separation explicitly", {
  core <- data.frame(
    a = rep(c("A", "A", "B"), each = 20),
    b = rep(c("B", "C", "C"), each = 20),
    stringsAsFactors = FALSE)
  core$winner <- c(rep(c("A", "B"), 10),
                   rep(c("A", "C"), 10),
                   rep(c("B", "C"), 10))
  boundary <- data.frame(
    a = rep("D", 60), b = rep(c("A", "B", "C"), each = 20),
    winner = "D", stringsAsFactors = FALSE)
  d <- rbind(core, boundary)
  d$judge <- rep(sprintf("J%02d", 1:20), length.out = nrow(d))
  d$panel <- "P1"

  ordinary <- btl(d, "a", "b", "winner", judge = "judge")
  expect_true(ordinary$converged)
  expect_true(ordinary$objects$extreme[ordinary$objects$object == "D"])

  err <- tryCatch({
    btl_efrm(d, "a", "b", "winner", "judge", panels = "panel",
      object_sets = list(S1 = c("A", "B", "C", "D")),
      se_method = "conditional", boot_reps = 0)
    NULL
  }, error = identity)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "within-set outcomes are separated")
  expect_match(conditionMessage(err), "set 'S1'")
  expect_match(conditionMessage(err), "not strongly connected")
  expect_false(grepl("within-set information .* singular",
                     conditionMessage(err)))
})

test_that("the within-set Ford check pools incomplete panel allocations", {
  make_rows <- function(a, b, panel, wins_a, wins_b) {
    out <- data.frame(a = a, b = b, panel = panel,
      winner = c(rep(a, wins_a), rep(b, wins_b)),
      stringsAsFactors = FALSE)
    judges <- if (panel == "P1") sprintf("J%02d", 1:10) else
      sprintf("J%02d", 11:20)
    out$judge <- rep(judges, length.out = nrow(out))
    out
  }
  # P1 observes only A--B; P2 observes B--C and A--C. Neither panel covers the
  # complete object graph, but their pooled directed graph is strongly connected
  # because every observed pair has outcomes in both directions.
  d <- rbind(
    make_rows("A", "B", "P1", 45, 15),
    make_rows("B", "C", "P2", 40, 20),
    make_rows("A", "C", "P2", 50, 10))

  fit <- btl_efrm(d, "a", "b", "winner", "judge", panels = "panel",
    object_sets = list(S1 = c("A", "B", "C")),
    se_method = "conditional", boot_reps = 0)
  expect_true(fit$converged)
  expect_equal(sort(fit$objects$object), c("A", "B", "C"))
  expect_equal(sort(fit$panels), c("P1", "P2"))
})
