# Validation of joint incomplete-panel DIF and multiple-maximum WLE scoring.
# Run from the package root. SV_PART=dif, public, or wle selects a component;
# the default runs all. Results are separate files with execution provenance.
options(rasch.max_workers = 1L, device = function(...) grDevices::pdf(NULL))
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
part <- Sys.getenv("SV_PART", "all")

design <- function(n, coverage, seed) {
  set.seed(seed)
  A <- factor(rep(c("a0", "a1"), each = n / 2))
  B <- factor(ifelse(runif(n) < ifelse(A == "a0", .15, .85), "b1", "b0"))
  ci <- factor(rep(1:3, length.out = n))
  follow <- if (coverage == "one-group") which(A == "a1") else
    which(runif(n) < ifelse(A == "a0", .25, .8))
  id <- c(seq_len(n), follow)
  data.frame(pid = id, A = A[id], B = B[id], ci = ci[id],
    occasion = factor(c(rep("pre", n), rep("post", length(follow)))))
}

if (part %in% c("all", "dif")) {
  scenarios <- data.frame(n = c(600, 300, 600, rep(600, 4)),
    coverage = c("one-group", "one-group", "differential", rep("one-group", 4)),
    effect = c(0, 0, 0, .4, .8, .4, .8),
    kind = c(rep("null", 3), "uniform", "uniform", "nonuniform", "nonuniform"),
    reps = c(3000, 300, 300, rep(200, 4)))
  rows <- list()
  terms <- c("A", "B", "ci", "occasion", "A:ci", "B:ci", "occasion:ci")
  for (s in seq_len(nrow(scenarios))) {
    sc <- scenarios[s, ]; d <- design(sc$n, sc$coverage, 92001 + s)
    w <- 1 / as.numeric(table(d$pid)[as.character(d$pid)])
    p <- matrix(NA_real_, sc$reps, 2, dimnames = list(NULL, c("A", "A:ci")))
    errors <- 0L
    set.seed(92100 + s)
    for (r in seq_len(sc$reps)) {
      person_error <- rnorm(sc$n, sd = .6)
      d$z <- 2 * (d$B == "b1") + .2 * (d$B == "b1") * as.numeric(d$ci) +
        .7 * (d$occasion == "post") * as.numeric(d$ci) +
        sc$effect * (d$A == "a1") *
          (if (sc$kind == "nonuniform") as.numeric(d$ci) - 2 else 1) +
        person_error[d$pid] + rnorm(nrow(d), sd = ifelse(d$A == "a0", .4, .8))
      ans <- tryCatch(.dif_type2(d, terms, variance = "cr3", cluster = d$pid,
        weights = w, report_terms = c("A", "A:ci")), error = identity)
      if (inherits(ans, "error")) { errors <- errors + 1L; next }
      if (!is.null(ans)) p[r, ] <- ans$p[match(colnames(p), ans$term)]
    }
    label <- paste(sc$n, sc$coverage, sc$kind, sc$effect)
    for (j in seq_len(ncol(p))) {
      ok <- is.finite(p[, j]); rate <- mean(p[ok, j] < .05)
      rows[[length(rows) + 1L]] <- sv_row("joint-dif-wle-maxima", label,
        paste("residual-model", colnames(p)[j], "rejection"), sum(ok),
        type1 = if (sc$kind == "null") rate else NA_real_,
        power = if (sc$kind != "null") rate else NA_real_, effect = sc$effect,
        n_attempted = sc$reps, n_refused = 0L, n_nonconv = 0L,
        n_error = errors, n_withheld = sum(!ok) - errors,
        n_metric_unavailable = 0L,
        notes = paste("Fixed factor design; correlated A/B; known Gaussian",
          "residual model with random person intercept and unequal residual",
          "variance; not an end-to-end Rasch calibration claim."))
    }
    cat(label, "rates", colMeans(p < .05, na.rm = TRUE), "errors", errors, "\n")
  }
  sv_write(do.call(rbind, rows), "joint-dif-residual-model")
}

if (part %in% c("all", "public")) {
  reps <- 200L; p <- rep(NA_real_, reps); errors <- refused <- nonconv <- 0L
  for (r in seq_len(reps)) {
    d <- design(600, "one-group", 93000 + r)
    set.seed(94000 + r); theta <- rnorm(600)
    X <- sapply(seq(-1.5, 1.5, length.out = 10), function(b)
      rbinom(nrow(d), 1, plogis(theta[d$pid] - b)))
    X[, 3] <- rbinom(nrow(d), 1,
      plogis(theta[d$pid] + .5 - 2 * (d$B == "b1")))
    colnames(X) <- paste0("I", 1:10)
    fit <- tryCatch(rasch(X, id = d$pid,
      factors = d[c("A", "B", "occasion")]), error = identity)
    if (inherits(fit, "error")) { errors <- errors + 1L; next }
    if (!isTRUE(fit$est$converged)) { nonconv <- nonconv + 1L; next }
    out <- tryCatch(dif_anova(fit, within = "occasion", n_groups = 3),
      error = identity)
    if (inherits(out, "error")) { errors <- errors + 1L; next }
    p[r] <- out$terms$p[out$terms$item == "I3" & out$terms$term == "A"]
    if (r %% 25 == 0) cat("Public Rasch replicates", r, "\n")
  }
  ok <- is.finite(p)
  rows <- list(sv_row("joint-dif-wle-maxima", "600 one-group public Rasch",
    "null A test on item with B DIF", sum(ok), type1 = mean(p[ok] < .05),
    n_attempted = reps, n_refused = refused, n_nonconv = nonconv,
    n_error = errors, n_withheld = sum(!ok) - errors - nonconv,
    n_metric_unavailable = 0L,
    notes = "End-to-end dichotomous calibration and estimated class intervals; A has no direct item effect given B. Limited 200-replicate screen."))
  sv_write(do.call(rbind, rows), "joint-dif-public-screen")
  cat("Public Rasch null rejection", mean(p[ok] < .05), "usable", sum(ok), "\n")
}

if (part %in% c("all", "wle")) {
  set.seed(95001); deficits <- numeric(100); withheld <- integer(100)
  for (r in 1:100) {
    gap <- runif(1, 0, 16)
    tau <- lapply(seq_len(8), function(j)
      sort(rnorm(if (r %% 2) 1L else 2L, sd = .7)) + if (j <= 3) -gap else gap)
    curve <- .person_wle_curve(tau, 1)
    grid <- seq(min(unlist(tau)) - 12, max(unlist(tau)) + 12, length.out = 4001)
    v <- curve$evaluate((grid - curve$origin) * curve$scale)
    for (score in 0:sum(lengths(tau))) {
      estimate <- .person_wle_maximum(curve, score)
      if (!is.finite(estimate)) { withheld[r] <- withheld[r] + 1L; next }
      u <- (estimate - curve$origin) * curve$scale
      at <- curve$evaluate(u)
      value <- score * u - at[, "psi"] + log(at[, "info"]) / 2
      reference <- score * (grid - curve$origin) * curve$scale -
        v[, "psi"] + log(v[, "info"]) / 2
      deficits[r] <- max(deficits[r], max(reference, na.rm = TRUE) - value)
    }
  }
  rows <- list(sv_row("joint-dif-wle-maxima", "100 dichotomous/PCM banks; gaps 0-16",
    "bank with dense-grid objective above selected WLE by >1e-7", 100,
    type1 = mean(deficits > 1e-7), n_attempted = 100L, n_refused = 0L,
    n_nonconv = 0L, n_error = 0L, n_withheld = sum(withheld > 0),
    n_metric_unavailable = 0L,
    notes = paste("Numerical conformance, not inferential Type I; all scores tested.",
      "Maximum positive dense-grid objective advantage:", max(deficits),
      "Unavailable individual scores:", sum(withheld))))
  sv_write(do.call(rbind, rows), "wle-separated-bank-maxima")
  cat("Max grid advantage", max(deficits), "withheld scores", sum(withheld), "\n")
}
