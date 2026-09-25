# Bundled example data used by the graphical interface. Kept in a separate
# source file so the exact datasets can also be reconstructed by the R-code
# disclosures without starting the app.

# Polytomous workflow example. These defaults are the simulation printed in
# the Rasch workflow vignette, so its code and application walkthrough analyse
# the same observations.
.demo_data <- function(seed = 17, Np = 600) {
  simulate_rasch(
    n_persons = Np,
    n_items = 12,
    model = "PCM",
    n_categories = 4,
    difficulty = c(-1.5, 1.5),
    disordered = "I04",
    dependence = list(pairs = list(c("I10", "I11")), strength = 1.3),
    dif = list(items = "I08", uniform = 0.8),
    n_groups = 3,
    seed = seed
  )
}

# dichotomous demo: 15 multiple-choice items (raw A-D responses), DIF planted
# on I05 by group, and I07 deliberately miskeyed (true correct C, key says A)
.demo_dich <- function(seed = 41, Np = 1000) {
  set.seed(seed)
  d <- seq(-2, 2, length.out = 15)
  grp <- rep(c("reference", "focal"), each = Np / 2)
  sex <- sample(c("female", "male"), Np, replace = TRUE)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(seq_along(d), function(i) {
    sft <- if (i == 5) ifelse(grp == "focal", 0.8, 0) else 0
    correct <- if (i == 7) "C" else "A"
    ok <- rbinom(Np, 1, plogis(th - d[i] - sft))
    ifelse(ok == 1, correct,
           sample(setdiff(c("A", "B", "C", "D"), correct), Np, replace = TRUE))
  })
  colnames(X) <- sprintf("I%02d", seq_along(d))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
             group = grp, sex = sex, check.names = FALSE)
}

# the demo key: all "A" (so I07 is the discoverable miskey)
.demo_dich_key <- function()
  setNames(rep("A", 15), sprintf("I%02d", 1:15))

# Comparative judgement workflow example. Exposure and erratic judging are
# planted; the two crossed judge factors are null and remain above the
# inference boundary. Keeping residual dimensionality out of this dataset
# prevents it from being mistaken for a history effect.
.demo_btl <- function(seed = 2, reps = 84) {
  d <- simulate_btl(
    n_objects = 8,
    n_judges = 48,
    reps_per_pair = reps,
    erratic_judges = 2 / 48,
    dependence = list(exposure = 0.7, carry_over = 0),
    seed = seed
  )
  judge_number <- as.integer(sub("^J", "", d$judge))
  d$panel <- factor(ifelse(judge_number %% 2L,
                           "panel A", "panel B"))
  d$experience <- factor(ifelse(judge_number <= 24,
                                "experienced", "novice"))
  d
}

# rankings demo: each judge ranks three random sets of five of eight
# objects, choosing successively from those left with probability
# proportional to exp(location); judge 7 ranks at random. One ranked
# object per row, the layout pl() reads.
.demo_pl <- function(seed = 5, n_judges = 24, sets = 3) {
  set.seed(seed)
  beta <- stats::setNames(seq(-1.75, 1.75, length.out = 8),
                          sprintf("O%d", 1:8))
  rows <- list()
  for (j in seq_len(n_judges)) {
    for (k in seq_len(sets)) {
      rem <- sample(names(beta), 5)
      ord <- character(0)
      while (length(rem) > 1) {
        pick <- if (j == 7L) sample(rem, 1)
          else sample(rem, 1, prob = exp(beta[rem]))
        ord <- c(ord, pick)
        rem <- setdiff(rem, pick)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        judge = sprintf("J%02d", j), ranking = (j - 1L) * sets + k,
        object = c(ord, rem), rank = 1:5, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# joint calibration demo: 10 dichotomous items answered by 250 persons and
# judged directly as well, in 300 paired comparisons of difficulty (unit 0.7
# of the response scale) and 60 rankings of four items (unit 1.1). Both
# judgement frames place I07 1.2 logits harder than the responses do. The
# judgements travel as attributes of the response table, so one dataset
# carries every frame of the joint calibration.
.demo_cj <- function(seed = 16, Np = 250) {
  set.seed(seed)
  delta <- stats::setNames(seq(-1.6, 1.6, length.out = 10),
                           sprintf("I%02d", 1:10))
  theta <- rnorm(Np, 0, 1.2)
  X <- sapply(delta, function(d) as.integer(runif(Np) < plogis(theta - d)))
  judged <- delta
  judged["I07"] <- judged["I07"] + 1.2
  pairs <- t(utils::combn(names(delta), 2))[sample(45, 300, replace = TRUE), ]
  p_a <- plogis(0.7 * (judged[pairs[, 1]] - judged[pairs[, 2]]))
  comparisons <- data.frame(
    object_a = pairs[, 1], object_b = pairs[, 2],
    winner = ifelse(runif(300) < p_a, pairs[, 1], pairs[, 2]),
    stringsAsFactors = FALSE)
  rankings <- do.call(rbind, lapply(1:60, function(j) {
    rem <- sample(names(delta), 4)
    ord <- character(0)
    while (length(rem) > 1) {
      pick <- sample(rem, 1, prob = exp(1.1 * judged[rem]))
      ord <- c(ord, pick)
      rem <- setdiff(rem, pick)
    }
    data.frame(ranking = j, item = c(ord, rem), rank = 1:4,
               stringsAsFactors = FALSE)
  }))
  out <- data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
                    check.names = FALSE)
  attr(out, "comparisons") <- comparisons
  attr(out, "rankings") <- rankings
  out
}

# rating scale demo: common step structure, item locations vary
.demo_rsm <- function(seed = 51, Np = 1000) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  loc <- seq(-1.2, 1.2, length.out = 8)
  step <- c(-0.9, 0.0, 0.9)
  grp <- rep(c("reference", "focal"), each = Np / 2)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(loc, function(b) sapply(th, function(t)
    sample(0:3, 1, prob = simP(t, b + step))))
  colnames(X) <- sprintf("R%02d", seq_along(loc))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X,
             group = grp, check.names = FALSE)
}

# rated (MFRM) demo, wide layout: 5 item columns, 6 raters (one erratic),
# incomplete design — one row per person-by-rater combination. The responses
# are simulated in long form (same structure and seed as always) and
# reshaped, so results are unchanged.
.demo_mfrm <- function(seed = 21, Np = 250) {
  set.seed(seed)
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  persons <- sprintf("P%04d", seq_len(Np)); raters <- paste0("Rater_", 1:6)
  th <- setNames(rnorm(Np, 0, 1.3), persons)
  rho <- setNames(c(-0.9, -0.4, -0.1, 0.1, 0.4, 0.9), raters)
  tau <- list(Essay = c(-1.2, 0.2, 1.1), Argument = c(-0.8, 0.5, 1.3),
              Evidence = c(-1.5, -0.2, 0.9), Style = c(-0.6, 0.4, 1.2),
              Mechanics = c(-1.0, 0.0, 1.0))
  d <- expand.grid(person = persons, item = names(tau), rater = raters,
                   stringsAsFactors = FALSE)
  seen <- unlist(lapply(persons, function(p) paste(p, sample(raters, 3))))
  d <- d[paste(d$person, d$rater) %in% seen, ]
  d$score <- mapply(function(p, i, r) {
    if (r == "Rater_6" && runif(1) < 0.2) return(sample(0:3, 1))  # erratic rater
    sample(0:3, 1, prob = simP(th[p], tau[[i]] + rho[r]))
  }, d$person, d$item, d$rater)
  # wide: one row per person-by-rater with one column per item
  w <- reshape(d, idvar = c("person", "rater"), timevar = "item",
               v.names = "score", direction = "wide")
  names(w) <- sub("^score\\.", "", names(w))
  w <- w[order(w$person, w$rater), c("person", "rater", names(tau))]
  rownames(w) <- NULL
  w
}

# frames demo: 2 person groups x 3 item sets with distinct units
.demo_efrm <- function(seed = 31, per_g = 350) {
  set.seed(seed)
  simP <- function(th, tau, r) { x <- 0:length(tau); p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p) }
  glev <- c("year5", "year7"); grp <- rep(glev, each = per_g); Np <- length(grp)
  phi <- c(year5 = 0.8, year7 = 1.25)
  sets <- rep(c("Number", "Algebra", "Space"), each = 6)
  alpha <- c(Number = 0.75, Algebra = 1.0, Space = 4 / 3)
  th <- rnorm(Np, 0, 1.3) + ifelse(grp == "year7", 0.5, 0)
  d <- as.numeric(sapply(c(-0.3, 0.1, 0.2), function(m) m + seq(-1.2, 1.2, length.out = 6)))
  X <- sapply(seq_along(sets), function(i) sapply(seq_len(Np), function(n)
    sample(0:2, 1, prob = simP(th[n], d[i] + c(-0.5, 0.5),
                               alpha[sets[i]] * phi[grp[n]]))))
  colnames(X) <- sprintf("%s_%02d", sets, seq_along(sets))
  data.frame(person_id = sprintf("P%04d", seq_len(Np)), X, year_group = grp,
             check.names = FALSE)
}
