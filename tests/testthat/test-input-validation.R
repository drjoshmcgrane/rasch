# Input honesty: misspelled column names, fractional scores, missing
# identifiers, unusable equating items, and HTML escaping. Each of these
# used to fail silently (wrong-data analyses, truncated responses, phantom
# "NA" facet levels, all-NA equating, markup injection in reports).

mkX <- function(n = 150, L = 6, seed = 1) {
  set.seed(seed)
  X <- matrix(rbinom(n * L, 1, plogis(outer(rnorm(n), seq(-1, 1, length.out = L), "-"))),
              n, L, dimnames = list(NULL, paste0("I", 1:L)))
  X
}

test_that("misspelled id, factor, and item columns are errors, not fallbacks", {
  d <- as.data.frame(mkX())
  d$sex <- rep(c("m", "f"), length.out = nrow(d))
  d$pid <- sprintf("P%03d", seq_len(nrow(d)))
  expect_error(rasch(d, id = "person_id", items = paste0("I", 1:6)),
               "id column 'person_id' not found")
  expect_error(rasch(d, factors = c("sex", "agee"), items = paste0("I", 1:6)),
               "factor column\\(s\\) not found.*agee")
  expect_error(rasch(d, items = c("I1", "I2", "Item3")),
               "item column\\(s\\) not found.*Item3")
  # correct names still work, including numeric item indices
  f <- rasch(d, id = "pid", factors = "sex", items = paste0("I", 1:6))
  expect_equal(ncol(f$X), 6L)
  f2 <- rasch(as.data.frame(mkX()), items = 1:6)
  expect_equal(ncol(f2$X), 6L)
})

test_that("fractional scores error instead of silently truncating", {
  X <- mkX()
  Xf <- as.data.frame(X)
  Xf$I3[5] <- 1.9
  expect_error(rasch(Xf), "non-integer score\\(s\\) in: I3.*1\\.9")
  # integer-valued doubles ("2.0") are fine
  Xd <- as.data.frame(X * 1.0)
  expect_s3_class(rasch(Xd), "rasch")
})

test_that("MFRM rows with missing identifiers are dropped with a note", {
  set.seed(3)
  long <- data.frame(person = rep(sprintf("P%03d", 1:80), each = 4),
                     rater = rep(c("R1", "R2"), 160),
                     item = rep(rep(c("A", "B"), each = 2), 80),
                     score = rbinom(320, 2, 0.5))
  long$rater[c(5, 9)] <- NA
  long$person[17] <- NA
  f <- rasch_mfrm(long, person = "person", item = "item", score = "score",
                  facets = "rater")
  expect_true(any(grepl("3 row\\(s\\) dropped: missing person, item, or facet",
                        f$notes)))
  # no phantom "NA" rater level
  expect_false("NA" %in% f$facet_effects$level)
})

test_that("equating excludes unusable items instead of returning all NA", {
  f <- rasch(as.data.frame(mkX(400, 8)))
  ref <- data.frame(item = paste0("I", 1:8),
                    location = f$items$location + 0.3,
                    se = f$items$se)
  ref$se[2] <- NA                       # e.g. a weakly determined item
  eq <- equate_tests(f, ref)
  expect_true(is.finite(eq$shift) && is.finite(eq$rmsd))
  expect_equal(eq$n, 7L)
  expect_true(is.na(eq$table$t[eq$table$item == "I2"]))
  expect_equal(sum(is.finite(eq$table$t)), 7L)
  expect_match(eq$note, "I2")
  expect_no_error(plot_equate(f, ref))
  ref$se[1:7] <- NA                     # fewer than two usable -> error
  expect_error(equate_tests(f, ref), "fewer than two common items")
})

test_that("report_html escapes data-derived text", {
  X <- mkX(200, 6)
  colnames(X) <- c(paste0("I", 1:5), "A<b>&x")
  f <- rasch(as.data.frame(X))
  out <- tempfile(fileext = ".html")
  report_html(f, out, title = "T <script>alert(1)</script>")
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_false(grepl("<script>alert", html, fixed = TRUE))
  expect_true(grepl("&lt;script&gt;alert", html, fixed = TRUE))
  expect_false(grepl("A<b>&x", html, fixed = TRUE))
  expect_true(grepl("A&lt;b&gt;&amp;x", html, fixed = TRUE) ||
              grepl("A&lt;b>&amp;x", html, fixed = TRUE))
  unlink(out)
})

test_that("every public estimator rejects fractional scores", {
  X <- matrix(c(0, 1, 1.9, 1, 0, 1, 0, 1, 1, 0, 1, 0), 6, 2,
              dimnames = list(NULL, c("A", "B")))
  expect_error(pcml(X), "non-integer")
  expect_error(pcml_pc(X), "non-integer")
  long <- data.frame(person = rep(sprintf("P%02d", 1:20), each = 2),
                     item = rep(c("A", "B"), 20),
                     score = c(1.9, rep(c(0, 1, 1, 0), 9), 0, 1, 1))
  expect_error(rasch_mfrm(long, "person", "item", "score", facets = NULL),
               "non-integer")
  d <- data.frame(a = rep("X", 30), b = rep("Y", 30),
                  resp = rep(c(0, 1, 1.5), 10))
  expect_error(btl(d, "a", "b", response = "resp"), "non-integer")
})

test_that("item_moments is overflow-stable and person_wle survives wide items", {
  im <- item_moments(8, seq(-3, 3, length.out = 30))
  expect_true(all(is.finite(unlist(im))))
  expect_equal(sum(im$P), 1, tolerance = 1e-12)
  w <- person_wle(list(seq(-3, 3, length.out = 30)))
  expect_true(all(is.finite(w$theta)))
})

test_that("secondary simulated trait keeps the requested mean and sd", {
  d <- simulate_rasch(n_persons = 20000, n_items = 6, theta_mean = 2,
                      theta_sd = 1, second_dim = list(items = 4:6, rho = 0.5),
                      seed = 1)
  t2 <- attr(d, "truth")$theta2
  expect_false(is.null(t2))
  expect_lt(abs(mean(t2) - 2), 0.05)
  expect_lt(abs(sd(t2) - 1), 0.05)
})

test_that("a judge in two panels is rejected by btl_efrm", {
  d <- simulate_btl_efrm(n_objects_per_set = 5, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 6, reps_within = 15,
                         reps_cross = 15, seed = 5)
  jj <- unique(d$judge)[1]
  d$panel[d$judge == jj][1] <- setdiff(unique(d$panel), d$panel[d$judge == jj][1])[1]
  expect_error(
    btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel",
             attr(d, "truth")$object_sets, se_method = "conditional"),
    "more than one panel")
})

test_that("OSI is withheld when the clustered covariance is rank-deficient", {
  set.seed(6)
  K <- 10; beta <- seq(-1.5, 1.5, length.out = K); n <- 600
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  win <- rbinom(n, 1, plogis(beta[ia] - beta[ib]))
  d <- data.frame(object_a = paste0("O", ia), object_b = paste0("O", ib),
                  winner = paste0("O", ifelse(win == 1, ia, ib)),
                  judge = sample(sprintf("J%d", 1:5), n, TRUE))
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  expect_true(is.na(f$osi$PSI))
  expect_true(any(grepl("OSI is withheld", f$notes)))
})

test_that("alpha is NA, not -Inf, when the total score is constant", {
  X <- cbind(I1 = rep(c(0L, 1L), 40), I2 = rep(c(1L, 0L), 40))
  f <- suppressWarnings(rasch(X, n_groups = 2))
  expect_true(is.na(f$alpha$alpha))
})

test_that("factor scores cannot bypass the integer guard", {
  long <- data.frame(person = rep(sprintf("P%02d", 1:20), each = 2),
                     item = rep(c("A", "B"), 20),
                     score = factor(c("1.9", rep(c("0", "1", "1", "0"), 9),
                                      "0", "1", "1")))
  expect_error(rasch_mfrm(long, "person", "item", "score", facets = NULL),
               "non-integer")
  expect_error(pcml(matrix(c(0, 1, Inf, 1, 0, 1, 0, 1), 4, 2)), "non-finite")
  expect_error(pcml(matrix(c("0", "1", "abc", "1", "0", "1", "0", "1"), 4, 2)),
               "non-numeric")
})

test_that("graded btl requires an ordered factor", {
  d <- data.frame(a = rep(c("X", "Y", "Z"), 40), b = rep(c("Y", "Z", "X"), 40))
  set.seed(1)
  resp <- sample(c("worse", "same", "better"), 120, TRUE)
  d$plain <- factor(resp)
  d$ord <- factor(resp, levels = c("worse", "same", "better"), ordered = TRUE)
  expect_error(btl(d, "a", "b", response = "plain"), "ORDERED")
  f <- btl(d, "a", "b", response = "ord")
  expect_identical(f$categories, c("worse", "same", "better"))
})

test_that("clustered dependence tests use a t reference with G - 1 df", {
  set.seed(3)
  K <- 6; b <- seq(-1, 1, length.out = K); n <- 800
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  d <- data.frame(object_a = paste0("O", ia), object_b = paste0("O", ib),
                  winner = paste0("O", ifelse(
                    rbinom(n, 1, plogis(b[ia] - b[ib] + 0.4)) == 1, ia, ib)),
                  judge = sample(sprintf("J%d", 1:6), n, TRUE))
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge",
           position = TRUE)
  dp <- f$dependence
  expect_equal(unique(dp$df), 5L)
  expect_true("t" %in% names(dp))                   # labelled for its reference
  expect_equal(dp$p, 2 * pt(-abs(dp$t), df = 5), tolerance = 1e-12)
  expect_true(all(dp$p >= 2 * pnorm(-abs(dp$t))))   # wider than normal theory
})

test_that("simulate_rasch validates the second-dimension specification", {
  expect_error(simulate_rasch(50, 6, second_dim = list(items = 4:6, rho = 1.2)),
               "correlation in")
  expect_error(simulate_rasch(50, 6, second_dim = list(items = "I99", rho = .5)),
               "unknown item")
})

test_that("margin ordering must be explicit (ordered factor or numeric)", {
  d <- data.frame(a = rep(c("X", "Y", "Z"), 40), b = rep(c("Y", "Z", "X"), 40))
  set.seed(1)
  d$win <- ifelse(runif(120) < .5, d$a, d$b)
  mg <- sample(c("small", "large"), 120, TRUE)
  d$m_plain <- factor(mg); d$m_chr <- mg
  d$m_ord <- factor(mg, levels = c("small", "large"), ordered = TRUE)
  expect_error(btl(d, "a", "b", winner = "win", margin = "m_plain"), "ORDERED")
  expect_error(btl(d, "a", "b", winner = "win", margin = "m_chr"), "character")
  expect_s3_class(btl(d, "a", "b", winner = "win", margin = "m_ord"),
                  "rasch_btl")
})

test_that("invalid graded numeric responses error instead of being dropped", {
  d <- data.frame(a = rep(c("X", "Y", "Z"), 40), b = rep(c("Y", "Z", "X"), 40))
  d$r_chr <- rep(c("0", "1", "abc"), 40)
  d$r_inf <- rep(c(0, 1, Inf), 40)
  expect_error(btl(d, "a", "b", response = "r_chr"), "non-numeric")
  expect_error(btl(d, "a", "b", response = "r_inf"), "non-finite")
})

test_that("simulator rejects malformed second-dimension specifications", {
  expect_error(simulate_rasch(50, 6, second_dim = list(items = 4.9, rho = .5)),
               "whole numbers")
  expect_error(simulate_rasch(50, 6, second_dim = list(items = integer(0),
                                                       rho = .5)),
               "at least one item")
  expect_error(simulate_rasch(50, 6, second_dim = list(items = 4:6,
                                                       rho = c(.5, .6))),
               "single correlation")
})

test_that("fits saved before the t rename still print their dependence", {
  set.seed(3)
  K <- 6; b <- seq(-1, 1, length.out = K); n <- 400
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  d <- data.frame(object_a = paste0("O", ia), object_b = paste0("O", ib),
                  winner = paste0("O", ifelse(
                    rbinom(n, 1, plogis(b[ia] - b[ib] + 0.4)) == 1, ia, ib)),
                  judge = sample(sprintf("J%d", 1:8), n, TRUE))
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge",
           position = TRUE)
  # 1.11.4 transitional schema: z name, df present, t-based p -> label t
  trans <- f
  names(trans$dependence)[names(trans$dependence) == "t"] <- "z"
  expect_output(print(trans), "t = ")
  # pre-1.11.4 schema: z, NO df, normal-reference p -> keep the z label,
  # since relabelling it t would misrepresent how its p was computed
  legacy <- f
  names(legacy$dependence)[names(legacy$dependence) == "t"] <- "z"
  legacy$dependence$df <- NULL
  legacy$dependence$p <- 2 * pnorm(-abs(legacy$dependence$z))
  expect_output(print(legacy), "z = ")
  expect_output(print(f), "t = ")
})

test_that("MFRM wide input carries person factors through the melt", {
  set.seed(21)
  N <- 80; I <- 4
  th <- rnorm(N); del <- seq(-1, 1, length.out = I)
  sev <- c(R1 = -0.3, R2 = 0.3)
  # every person is rated by BOTH raters (one wide row per person-rater
  # combination): persons link the raters, which the conditional
  # likelihood requires -- a fully nested rater design would be refused
  # by the connectivity check, and rightly so
  wide <- data.frame(pid = rep(seq_len(N), each = 2),
                     rater = rep(names(sev), N),
                     grp = rep(sample(c("boy", "girl"), N, replace = TRUE),
                               each = 2))
  for (i in seq_len(I))
    wide[[paste0("it", i)]] <- rbinom(2 * N, 1,
      plogis(th[wide$pid] - del[i] - sev[wide$rater]))
  fw <- rasch_mfrm(wide, person = "pid", items = paste0("it", 1:I),
                   facets = "rater", factors = "grp")
  long <- reshape(wide, direction = "long", varying = paste0("it", 1:I),
                  v.names = "score", timevar = "item",
                  times = paste0("it", 1:I), idvar = "..rid")
  fl <- rasch_mfrm(long, person = "pid", item = "item", score = "score",
                   facets = "rater", factors = "grp")
  d1 <- dif_anova(fw, factors = "grp")
  d2 <- dif_anova(fl, factors = "grp")
  expect_gt(nrow(d1$summary), 0)
  expect_equal(d1$summary$F_uniform, d2$summary$F_uniform, tolerance = 1e-8)
  expect_equal(d1$summary$p_uniform, d2$summary$p_uniform, tolerance = 1e-8)
  # a data-frame factor is replicated row-wise the same way
  fw2 <- rasch_mfrm(wide, person = "pid", items = paste0("it", 1:I),
                    facets = "rater",
                    factors = data.frame(grp = wide$grp))
  expect_equal(dif_anova(fw2, factors = "grp")$summary$F_uniform,
               d2$summary$F_uniform, tolerance = 1e-8)
  # misspelled wide factor column errors instead of silently dropping
  expect_error(rasch_mfrm(wide, person = "pid", items = paste0("it", 1:I),
                          facets = "rater", factors = "grpp"),
               "not found")
})

test_that("MFRM duplicate person-by-cell responses are an error", {
  set.seed(22)
  d <- data.frame(pid = rep(1:40, each = 3),
                  item = rep(c("A", "B", "C"), 40),
                  rater = rep(c("R1", "R2"), 60),
                  score = rbinom(120, 1, 0.5))
  f <- rasch_mfrm(d, person = "pid", item = "item", score = "score",
                  facets = "rater")
  expect_s3_class(f, "rasch_mfrm")
  # a repeated row would make the kept response depend on row order
  expect_error(rasch_mfrm(d[c(seq_len(nrow(d)), 1L), ], person = "pid",
                          item = "item", score = "score", facets = "rater"),
               "duplicate")
})

test_that("MFRM structurally confounded facet designs are an error", {
  set.seed(23)
  d <- expand.grid(pid = 1:60, item = c("A", "B", "C", "D"),
                   stringsAsFactors = FALSE)
  # each rater sees a disjoint half of the items: severity is confounded
  # with the item locations
  d$rater <- ifelse(d$item %in% c("A", "B"), "R1", "R2")
  d$score <- rbinom(nrow(d), 1, 0.5)
  expect_error(rasch_mfrm(d, person = "pid", item = "item", score = "score",
                          facets = "rater"),
               "unidentified")
})

test_that("BTL refuses directed separation of the win graph (Ford 1957)", {
  set.seed(24)
  # O1, O2 never lose to O3, O4: no finite ML locations exist, and the
  # optimiser's boundary values must not be presented as a converged fit
  d <- data.frame(a = c(rep("O1", 15), rep("O2", 15), rep("O1", 10),
                        rep("O3", 10)),
                  b = c(rep("O3", 15), rep("O4", 15), rep("O2", 10),
                        rep("O4", 10)))
  d$win <- c(rep("O1", 15), rep("O2", 15),
             ifelse(runif(10) < .5, "O1", "O2"),
             ifelse(runif(10) < .5, "O3", "O4"))
  expect_error(btl(d, "a", "b", winner = "win"), "not strongly connected")
  # a design with wins in both directions across the divide is untouched
  d2 <- d
  d2$win[1:2] <- c("O3", "O3")
  expect_silent(f2 <- btl(d2, "a", "b", winner = "win"))
  expect_true(f2$converged)
})

test_that("zero-information response pairs do not connect item blocks", {
  set.seed(1)
  N <- 120
  th <- rnorm(N)
  X <- matrix(NA_integer_, N + 1, 4,
              dimnames = list(NULL, c("A", "B", "C", "D")))
  X[1:60, 1:2] <- cbind(rbinom(60, 1, plogis(th[1:60])),
                        rbinom(60, 1, plogis(th[1:60] - 0.5)))
  X[61:120, 3:4] <- cbind(rbinom(60, 1, plogis(th[61:120] + 0.5)),
                          rbinom(60, 1, plogis(th[61:120])))
  # the only bridge respondent scores (0, 0) on the B-C pair: total zero
  # has a single feasible conditional allocation and carries no
  # information, so the blocks stay unlinked
  X[121, c("B", "C")] <- c(0L, 0L)
  expect_error(rasch(X), "not connected")
  # a SINGLE informative bridge is still perfect separation in the
  # conditional pair logit: the pair MLE runs to the boundary, the
  # information vanishes at the solution, and the projected-information
  # backstop refuses what the graph check alone cannot see
  X[121, c("B", "C")] <- c(1L, 0L)
  expect_error(rasch(X), "singular")
  # two bridges in opposite directions give an interior maximum: a real link
  X <- rbind(X, NA_integer_)
  X[122, c("B", "C")] <- c(0L, 1L)
  f <- rasch(X)
  expect_true(f$est$converged)
  expect_true(all(f$items$se > 0, na.rm = TRUE))
  # MFRM analogue: an extreme-total bridge does not join the blocks
  set.seed(2)
  d <- expand.grid(pid = 1:80, item = c("A", "B", "C", "D"),
                   rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d <- d[(d$pid <= 40 & d$item %in% c("A", "B")) |
         (d$pid >  40 & d$item %in% c("C", "D")), ]
  d$score <- rbinom(nrow(d), 1, 0.5)
  d <- rbind(d, data.frame(pid = 81L, item = c("B", "C"), rater = "R1",
                           score = c(0L, 0L)))
  expect_error(rasch_mfrm(d, person = "pid", item = "item", score = "score",
                          facets = "rater"),
               "does not bridge")
})

test_that("anchors relax the Ford condition to anchored recession bounds", {
  # two balanced components, each with its own anchor, joined only by
  # A always beating C: with both endpoints of the crossing edge fixed,
  # nothing diverges and every free object is tied to an anchor in both
  # directions -- the unrestricted Ford condition would wrongly refuse it
  d <- rbind(data.frame(a = "A", b = "B",
                        win = rep(c("A", "B"), each = 15)),
             data.frame(a = "C", b = "D",
                        win = rep(c("C", "D"), each = 15)),
             data.frame(a = "A", b = "C", win = rep("A", 10)))
  f <- btl(d, "a", "b", winner = "win", anchors = c(A = 1, C = -1))
  expect_true(f$converged)
  loc <- setNames(f$objects$location, f$objects$object)
  expect_equal(unname(loc["A"]), 1)
  expect_equal(unname(loc["C"]), -1)
  expect_lt(abs(loc["B"] - 1), 0.75)     # balanced against the anchor at 1
  expect_lt(abs(loc["D"] + 1), 0.75)
  # without the C anchor the C-D cluster can recede: still refused
  expect_error(btl(d, "a", "b", winner = "win", anchors = c(A = 1)),
               "not tied to an anchor")
})

test_that("MFRM refuses disconnected response blocks the facet map cannot bridge", {
  set.seed(43)
  d <- expand.grid(pid = 1:80, item = c("A", "B", "C", "D"),
                   rater = c("R1", "R2"), stringsAsFactors = FALSE)
  # persons 1-40 answer A/B only, persons 41-80 answer C/D only: B has
  # full algebraic rank, but no person compares the blocks, so their
  # relative locations are a flat direction of the conditional likelihood
  d <- d[(d$pid <= 40 & d$item %in% c("A", "B")) |
         (d$pid >  40 & d$item %in% c("C", "D")), ]
  d$score <- rbinom(nrow(d), 1, 0.5)
  expect_error(rasch_mfrm(d, person = "pid", item = "item", score = "score",
                          facets = "rater"),
               "does not bridge")
})

test_that("MFRM factors stay aligned when rows with missing identifiers drop", {
  set.seed(6)
  d <- data.frame(pid = rep(1:40, each = 3),
                  item = rep(c("A", "B", "C"), 40),
                  rater = rep(c("R1", "R2"), 60),
                  score = rbinom(120, 1, 0.5))
  d$sx <- sample(c("m", "f"), 40, TRUE)[d$pid]
  d$pid[5] <- NA
  f <- rasch_mfrm(d, person = "pid", item = "item", score = "score",
                  facets = "rater", factors = "sx")
  expect_true(any(grepl("dropped", f$notes)))
  expect_equal(nrow(f$factors), length(unique(d$pid[!is.na(d$pid)])))
  # data-frame factors with one row per original data row align the same way
  f2 <- rasch_mfrm(d, person = "pid", item = "item", score = "score",
                   facets = "rater", factors = data.frame(sx = d$sx))
  expect_equal(f2$factors$sx, f$factors$sx)
})

test_that("EFRM removes data-frame factor columns from the item matrix", {
  set.seed(5)
  N <- 300; g <- rep(c("g1", "g2"), each = 150); th <- rnorm(N, sd = 1.5)
  X <- as.data.frame(matrix(0L, N, 8)); names(X) <- paste0("v", 1:8)
  del <- rep(seq(-1.5, 1.5, length.out = 4), 2)
  for (i in 1:8) X[[i]] <- rbinom(N, 1, plogis(th - del[i]))
  X$grp <- g
  X$sex <- rep(c(1L, 2L), N / 2)
  # without items=, the numeric factor column must not become an item
  # (a single set keeps the test off the person-side linking machinery)
  f <- rasch_efrm(X, groups = "grp",
                  item_sets = list(A = paste0("v", 1:8)),
                  factors = data.frame(sex = X$sex))
  expect_setequal(unique(f$thresholds_arbitrary$item), paste0("v", 1:8))
})
