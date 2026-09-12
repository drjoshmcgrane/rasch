# Independent checks of the completed interval-policy validation outputs.
# Run from the package root after item-fit-bootstrap-intervals.R finishes.
stem <- "tools/simval/results/item-fit-bootstrap-intervals"
attempts <- read.csv(paste0(stem, "-attempts.csv"))
items <- read.csv(paste0(stem, "-items.csv"))
summary <- read.csv(paste0(stem, ".csv"))
scenarios <- c("linked booklets", "unequal exposure")
stopifnot(nrow(attempts) == 200L, nrow(summary) == 18L,
          setequal(attempts$scenario, scenarios))

for (x in list(attempts, items, summary)) {
  stopifnot(length(unique(x$script_md5)) == 1L,
    identical(unique(x$script_md5), unname(tools::md5sum(x$script[1]))),
    identical(unique(x$harness_md5), unname(tools::md5sum("tools/simval/harness.R"))),
    identical(unique(x$algorithm), "loo-maxt-2"))
}
rfiles <- sort(list.files("R", "[.]R$", full.names = TRUE))
tf <- tempfile()
writeLines(paste(unname(tools::md5sum(rfiles)), collapse = ""), tf)
rtree <- substr(unname(tools::md5sum(tf)), 1, 12)
unlink(tf)
# The recorded R-tree hash identifies the source the study ran against, which
# later commits move away from. Require the three files to agree on one hash,
# then report whether the working tree is still that source rather than
# failing the reconstruction over a state the recorded run cannot control.
recorded_rtree <- unique(c(attempts$r_tree_md5, items$r_tree_md5,
                           summary$r_tree_md5))
stopifnot(length(recorded_rtree) == 1L)
cat(if (identical(recorded_rtree, rtree))
      sprintf("R source tree is the one recorded for the run (%s).\n", rtree)
    else sprintf(paste0("R source tree has moved since the run: recorded %s, ",
                        "working tree %s. The rates below are reconstructed ",
                        "from the recorded outputs, not re-simulated against ",
                        "current code.\n"), recorded_rtree, rtree))

equal <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-12,
                                        check.attributes = FALSE))
for (scenario in scenarios) {
  a <- attempts[attempts$scenario == scenario, ]
  s <- summary[summary$scenario == scenario, ]
  stopifnot(nrow(a) == 100L, identical(sort(a$replicate), 1:100),
    all(a$n_boot_attempted == a$n_boot_used + a$n_boot_nonconv + a$n_boot_errors))
  stopifnot(all(is.finite(a$same_data_max_error[a$status == "ok"])),
    all(a$same_data_max_error[a$status == "ok"] < 1e-8))
  for (j in seq_len(nrow(s))) {
    metric <- s$quantity[j]
    ok <- a$status == "ok" & is.finite(a[[metric]])
    v <- a[[metric]][ok]
    field <- if (grepl("familywise$", metric)) "familywise" else "type1"
    stopifnot(s$n_attempted[j] == 100L, s$n_reps[j] == sum(ok),
      s$n_refused[j] == sum(a$status == "refused"),
      s$n_nonconv[j] == sum(a$status == "nonconverged"),
      s$n_error[j] == sum(a$status == "error"),
      s$n_metric_unavailable[j] == sum(a$status == "ok" & !is.finite(a[[metric]])),
      equal(s[[field]][j], if (length(v)) mean(v) else NA_real_),
      equal(s[[paste0("mc_se_", field)]][j],
            if (length(v) > 1L) sd(v) / sqrt(length(v)) else NA_real_))
    for (nm in c("n_boot_attempted", "n_boot_used", "n_boot_nonconv", "n_boot_errors"))
      stopifnot(s[[nm]][j] == sum(a[[nm]]))
  }
  for (r in a$replicate[a$status == "ok"]) {
    it <- items[items$scenario == scenario & items$replicate == r, ]
    row <- a[a$replicate == r, ]
    L <- if (scenario == "linked booklets") 15L else 8L
    stopifnot(nrow(it) == L, !anyDuplicated(it$item))
    for (st in c("chisq", "fit_resid", "infit_z", "outfit_z")) {
      p <- it[[paste0(st, "_p_boot")]]
      adj <- it[[paste0(st, "_p_boot_adj")]]
      stopifnot(equal(row[[paste0(st, "_marginal")]],
                       if (all(is.finite(p))) mean(p < .05) else NA_real_),
        equal(row[[paste0(st, "_familywise")]],
               if (all(is.finite(adj))) as.numeric(any(adj < .05)) else NA_real_))
      # Both one-sided chi-square and two-sided residual Holm tests must
      # have enough retained draws to cross the strict .05 threshold.
      multiplier <- if (st == "chisq") 1L else 2L
      stopifnot(all(multiplier * L / (it[[paste0("n_boot_", st)]] + 1) < .05))
    }
  }
}
cat("All 200 datasets: provenance, item probabilities, accounting, resolution and summaries verified.\n")
print(summary[, c("scenario", "quantity", "n_reps", "type1", "mc_se_type1",
                   "familywise", "mc_se_familywise")], row.names = FALSE)
for (scenario in scenarios) for (metric in c("chisq_familywise", "fit_resid_familywise")) {
  d <- attempts[attempts$scenario == scenario & attempts$status == "ok", ]
  x <- d[[metric]][is.finite(d[[metric]])]
  ci <- binom.test(sum(x), length(x))$conf.int
  cat(sprintf("%s %s: %g/%d; exact 95%% CI %.2f%% to %.2f%%\n",
                scenario, metric, sum(x), length(x), 100 * ci[1], 100 * ci[2]))
}
