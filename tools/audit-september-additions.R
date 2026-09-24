# Targeted null screens for the September 2026 additions. Run from the root:
# Rscript tools/audit-september-additions.R 400 output/september-audit
# These screens supplement, not replace, the wider simulation programme.
args <- commandArgs(TRUE)
B <- if (length(args)) as.integer(args[1]) else 400L
out <- if (length(args) >= 2L) args[2] else "output/september-audit"
stopifnot(is.finite(B), B >= 10L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
files <- sort(c(list.files("R", full.names = TRUE), "tools/audit-september-additions.R"))
source_hashes <- tools::md5sum(files)
started <- Sys.time()
pkgload::load_all(quiet = TRUE)
options(rasch.efrm_workers = 1L)
rows <- list()
record <- function(label, values, n_attempted) {
  x <- do.call(rbind, values)
  ok <- is.finite(x[, "p"])
  p <- x[ok, "p"]; n <- length(p)
  rejected <- sum(p < .05)
  rate <- if (n) rejected / n else NA_real_
  ci <- if (n) binom.test(rejected, n)$conf.int else c(NA_real_, NA_real_)
  data.frame(study = label, attempted = n_attempted, analysed = n,
    unavailable = n_attempted - n, type1 = rate,
    mcse = sqrt(rate * (1 - rate) / n), lower = ci[1], upper = ci[2],
    mean_difference = mean(x[ok, "difference"]),
    sd_over_mean_se = sd(x[ok, "difference"]) / mean(x[ok, "se"]))
}
simulate_rankings <- function(beta, n, judges = 20L) do.call(rbind, lapply(seq_len(n), function(i) {
  ord <- sample(names(beta), prob = exp(beta))
  data.frame(ranking = i, object = ord, rank = seq_along(ord),
             judge = paste0("J", (i - 1L) %% judges + 1L))
}))
for (study in c("pl_model", "pl_judges", "pl_judges_dependent",
                "pl_judges_dependent_50", "cj_fixed", "cj_estimated", "dif_wald", "dtf")) {
  result <- vector("list", B)
  errors <- rep(NA_character_, B)
  for (b in seq_len(B)) {
    set.seed(924000L + b)
    result[[b]] <- tryCatch({
      delta <- setNames(seq(-1.5, 1.5, length.out = 6L), paste0("I", 1:6))
      value <- if (startsWith(study, "pl")) {
        d <- simulate_rankings(delta, 120L)
        if (startsWith(study, "pl_judges_dependent")) {
          # Independent judges, each submitting identical
          # rankings. The marginal distribution remains exactly PL while
          # dependence within judge is perfect, not merely nominal.
          J <- if (study == "pl_judges_dependent_50") 50L else 20L
          repeats <- if (J == 20L) 6L else 3L
          base <- simulate_rankings(delta, J, J)
          d <- do.call(rbind, lapply(seq_len(repeats), function(i) {
            z <- base; z$ranking <- z$ranking + J * (i - 1L); z
          }))
        }
        f <- if (study == "pl_model") pl(d, se = "model") else pl(d, judge = "judge")
        v <- f$invariance$objects
        v <- v[v$object == "I3", , drop = FALSE]
        c(p = v$p, difference = v$difference, se = v$se)
      } else {
        theta <- rnorm(240)
        X <- sapply(delta, function(d) rbinom(length(theta), 1, plogis(theta - d)))
        if (startsWith(study, "cj")) {
          pairs <- t(combn(names(delta), 2))
          pr <- pairs[sample(nrow(pairs), 400, TRUE), ]
          a <- pr[, 1]; bb <- pr[, 2]
          u <- if (study == "cj_fixed") 1 else .7
          cj <- data.frame(a = a, b = bb,
            w = ifelse(runif(length(a)) < plogis(u * (delta[a] - delta[bb])), a, bb))
          f <- rasch_cj(X, cj, "a", "b", "w",
                        units = c(comparisons = if (study == "cj_fixed") 1 else NA))
          v <- f$invariance$items
          v <- v[v$item == "I3", ]
          if (!f$converged || !nrow(v)) stop("unavailable")
          c(p = v$p, difference = v$difference, se = v$se)
        } else {
          f <- rasch(data.frame(X, group = rep(c("a", "b"), each = 120)),
                     factors = "group")
          if (study == "dif_wald") {
            v <- dif_wald(f, items = "I3")$summary
            c(p = v$p, difference = v$shift, se = v$se)
          } else {
            v <- dtf(f, "group", items = c("I2", "I3"), grid = 5L)$test
            c(p = v$p, difference = v$shift_mean, se = v$se)
          }
        }
      }
      if (length(value) != 3L || any(!is.finite(value)))
        stop("the prespecified contrast is unavailable")
      value
    }, error = function(e) {
      errors[b] <<- conditionMessage(e)
      cat(study, b, conditionMessage(e), "\n")
      c(p = NA_real_, difference = NA_real_, se = NA_real_)
    })
    if (b %% 50L == 0L) cat(study, b, "of", B, "\n")
  }
  rows[[study]] <- record(study, result, B)
  saveRDS(result, file.path(out, paste0(study, ".rds")))
  saveRDS(errors, file.path(out, paste0(study, "-unavailable.rds")))
  print(rows[[study]], row.names = FALSE)
}
summary <- do.call(rbind, rows)
saveRDS(list(summary = summary, source_hashes = source_hashes, started = started,
             session = sessionInfo(), B = B), file.path(out, "summary.rds"))
print(summary, row.names = FALSE)
