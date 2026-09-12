# Refined purification study. The first screen rejected per-item leave-one-out
# matching and showed that recomputing residuals from a short anchor scale can
# lose substantial power. This study separates anchor-based class intervals
# from anchor-based residual recalibration, selects one strongest item per
# round, and checks the public split-and-refit resolution for uniform DIF.

# The pinned hash is the screen as committed. Its own hash moved when its
# stale helper pin was repaired, which postdates the recorded results, whose
# rows name the earlier screen they ran against
# (afd7a676ac870245ace7844f5e452ed7), so a rerun is distinguishable from them.

base_study <- "tools/simval/studies/dif-score-purification.R"
base_md5 <- "532cbc3095c6067cb3f09883e7ff3e99"
if (!identical(unname(tools::md5sum(base_study)), base_md5))
  stop("purification helper hash does not match the committed version")
ex2 <- parse(file = base_study)
for (i in seq_len(15L)) eval(ex2[[i]], envir = environment())

NREP <- as.integer(Sys.getenv("SV_REPS", "500"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "1")))

stats_anchor_mode <- function(fit, factors, model, anchors,
                              mode = c("ci-only", "recalibrated")) {
  mode <- match.arg(mode)
  st <- anchor_state(fit, factors, model, anchors)
  if (is.null(st)) return(NULL)
  if (mode == "ci-only") st$residuals <- fit$residuals
  stats_from_state(fit, factors, st)
}

strongest_pipeline <- function(fit, factors, model, mode,
                               max_rounds = 3L) {
  items <- colnames(fit$X)
  cur <- full_stats(fit)
  excluded <- first <- character(0)
  for (r in seq_len(max_rounds)) {
    cand <- cur[cur$p_adj < 0.05 & !cur$item %in% excluded, , drop = FALSE]
    if (!nrow(cand)) break
    item_p <- tapply(cand$p_adj, cand$item, min)
    pick <- names(which.min(item_p))[1L]
    if (!length(first)) first <- pick
    proposed <- c(excluded, pick)
    anchors <- setdiff(items, proposed)
    if (length(anchors) < MIN_ANCHORS) break
    nxt <- stats_anchor_mode(fit, factors, model, anchors, mode)
    if (is.null(nxt)) return(NULL)
    excluded <- proposed
    cur <- nxt
  }
  list(stats = cur, excluded = excluded, first = first,
       n_anchors = length(setdiff(items, excluded)))
}

resolve_summary <- function(fit, sc) {
  rr <- tryCatch(resolve_dif(fit), error = function(e) NULL)
  if (is.null(rr)) return(NULL)
  ds <- tryCatch(full_stats(rr$fit), error = function(e) NULL)
  if (is.null(ds)) return(NULL)
  base <- rr$fit$split_map %||% stats::setNames(ds$item, ds$item)
  item_base <- unname(base[ds$item])
  unaffected <- !item_base %in% "I3"
  c(first_target = nrow(rr$splits) > 0L && rr$splits$base_item[1L] == "I3",
    first_wrong = nrow(rr$splits) > 0L && rr$splits$base_item[1L] != "I3",
    target_split = "I3" %in% rr$splits$base_item,
    wrong_split = any(rr$splits$base_item != "I3"),
    unaffected_fwer = any(ds$p_adj[unaffected] < 0.05),
    any_final = any(ds$p_adj < 0.05), n_splits = rr$n_splits,
    n_nonuniform = rr$n_nonuniform)
}

one_refined <- function(seed, sc) {
  set.seed(seed)
  dat <- simulate_scenario(sc)
  fit <- fit_scenario(dat$X, dat$factors, sc$model)
  if (is.null(fit) || !isTRUE(fit$est$converged)) return(NULL)
  items <- colnames(fit$X)
  fs <- full_stats(fit)
  oci <- stats_anchor_mode(fit, dat$factors, sc$model,
                           setdiff(items, "I3"), "ci-only")
  ore <- stats_anchor_mode(fit, dat$factors, sc$model,
                           setdiff(items, "I3"), "recalibrated")
  pci <- strongest_pipeline(fit, dat$factors, sc$model, "ci-only")
  pre <- strongest_pipeline(fit, dat$factors, sc$model, "recalibrated")
  if (any(vapply(list(fs, oci, ore, pci, pre), is.null, TRUE))) return(NULL)
  dec <- c(full = method_decisions(fs, sc),
           oracle_ci = method_decisions(oci, sc),
           oracle_recal = method_decisions(ore, sc),
           pipeline_ci = method_decisions(pci$stats, sc),
           pipeline_recal = method_decisions(pre$stats, sc))
  sel <- c(
    pipeline_ci_first_target = identical(pci$first, "I3"),
    pipeline_ci_first_wrong = length(pci$first) && pci$first != "I3",
    pipeline_ci_target_excluded = "I3" %in% pci$excluded,
    pipeline_ci_wrong_excluded = any(pci$excluded != "I3"),
    pipeline_ci_n_anchors = pci$n_anchors,
    pipeline_recal_first_target = identical(pre$first, "I3"),
    pipeline_recal_first_wrong = length(pre$first) && pre$first != "I3",
    pipeline_recal_target_excluded = "I3" %in% pre$excluded,
    pipeline_recal_wrong_excluded = any(pre$excluded != "I3"),
    pipeline_recal_n_anchors = pre$n_anchors)
  rs <- if (sc$nonuniform == 0) resolve_summary(fit, sc) else
    setNames(rep(NA_real_, 8L), c("first_target", "first_wrong",
      "target_split", "wrong_split", "unaffected_fwer", "any_final",
      "n_splits", "n_nonuniform"))
  c(dec, sel, setNames(rs, paste0("resolve.", names(rs))))
}

run_refined <- function(sc) {
  seeds <- sample.int(.Machine$integer.max, NREP)
  z <- parallel::mclapply(seeds, one_refined, sc = sc, mc.cores = NCORE,
                          mc.preschedule = TRUE)
  refused <- sum(vapply(z, is.null, TRUE))
  z <- z[!vapply(z, is.null, TRUE)]
  if (!length(z)) stop("no replicate completed in ", sc$name)
  z <- do.call(rbind, z); n <- nrow(z)
  null <- sc$uniform == 0 && sc$nonuniform == 0
  effect <- max(abs(sc$uniform), abs(sc$nonuniform))
  mk <- function(quantity, col, field, note) {
    v <- z[, col]
    args <- list(study = "dif-score-purification-refined",
      scenario = sc$name, quantity = quantity, n_reps = n,
      n_attempted = NREP, n_refused = refused, n_nonconv = 0L,
      effect = effect, notes = paste0(sc$model, "; ", sc$design, "; ",
        note, "; base helper md5 ", base_md5))
    args[[field]] <- mean(v, na.rm = TRUE)
    do.call(sv_row, args)
  }
  rows <- list()
  methods <- c("full", "oracle_ci", "oracle_recal", "pipeline_ci",
               "pipeline_recal")
  for (method in methods) {
    label <- gsub("_", " ", method)
    if (null) rows[[length(rows) + 1L]] <- mk(
      paste(label, "Holm FWER"), paste0(method, ".fwer"), "familywise",
      "final reported family") else {
      rows[[length(rows) + 1L]] <- mk(
        paste(label, "target power"), paste0(method, ".target"), "power",
        paste("I3", sc$target_kind, "group term"))
      rows[[length(rows) + 1L]] <- mk(
        paste(label, "unaffected-item FWER"),
        paste0(method, ".unaffected_fwer"), "familywise", "I1-I2 and I4-I6")
      if (sc$design == "multifactor") rows[[length(rows) + 1L]] <- mk(
        paste(label, "non-target-factor FWER"),
        paste0(method, ".other_factor_fwer"), "familywise", "all region terms")
    }
  }
  for (mode in c("pipeline_ci", "pipeline_recal")) {
    label <- gsub("_", " ", mode)
    rows[[length(rows) + 1L]] <- mk(paste(label, "selected target first"),
      paste0(mode, "_first_target"), "power", "selection diagnostic")
    rows[[length(rows) + 1L]] <- mk(paste(label, "selected wrong item first"),
      paste0(mode, "_first_wrong"), "familywise", "selection diagnostic")
    rows[[length(rows) + 1L]] <- mk(paste(label, "excluded target"),
      paste0(mode, "_target_excluded"), "power", "selection diagnostic")
    rows[[length(rows) + 1L]] <- mk(paste(label, "excluded any wrong item"),
      paste0(mode, "_wrong_excluded"), "familywise", "selection diagnostic")
    rows[[length(rows) + 1L]] <- sv_row(
      study = "dif-score-purification-refined", scenario = sc$name,
      quantity = paste(label, "mean anchors retained"), n_reps = n,
      n_attempted = NREP, n_refused = refused, n_nonconv = 0L,
      effect = effect, mean_se = mean(z[, paste0(mode, "_n_anchors")]),
      se_ratio = NA_real_, notes = "mean_se stores the mean count")
  }
  if (sc$nonuniform == 0) {
    for (nm in c("first_target", "first_wrong", "target_split", "wrong_split",
                 "unaffected_fwer", "any_final"))
      rows[[length(rows) + 1L]] <- mk(paste("resolve", gsub("_", " ", nm)),
        paste0("resolve.", nm), if (grepl("target", nm)) "power" else
          "familywise", "public split-and-refit resolution")
    for (nm in c("n_splits", "n_nonuniform")) rows[[length(rows) + 1L]] <-
      sv_row(study = "dif-score-purification-refined", scenario = sc$name,
        quantity = paste("resolve mean", gsub("_", " ", nm)), n_reps = n,
        n_attempted = NREP, n_refused = refused, n_nonconv = 0L,
        effect = effect, mean_se = mean(z[, paste0("resolve.", nm)]),
        se_ratio = NA_real_, notes = "mean_se stores the mean count")
  }
  do.call(rbind, rows)
}

scenarios <- list(
  scenario("PCM", "imbalanced"),
  scenario("PCM", "imbalanced", "uniform", 0.6),
  scenario("PCM", "imbalanced", "nonuniform", 1.4),
  scenario("RSM", "imbalanced"),
  scenario("RSM", "imbalanced", "uniform", 0.6),
  scenario("RSM", "imbalanced", "nonuniform", 1.4),
  scenario("dichotomous", "multifactor"),
  scenario("dichotomous", "multifactor", "nonuniform", 1.4))

set.seed(8.71e7)
rows <- list()
for (i in seq_along(scenarios)) {
  message(sprintf("[%d/%d] %s", i, length(scenarios), scenarios[[i]]$name))
  rows[[i]] <- run_refined(scenarios[[i]])
  sv_write(do.call(rbind, rows), "dif-score-purification-refined")
}
