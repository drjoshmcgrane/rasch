# Validation of score purification after the extended bootstrap study found
# artificial DIF under partial alternatives. Each anchor method refits the
# calibration, recomputes person measures, residuals and class intervals, and
# then reconstructs the joint Type II DIF tests.

# The pinned hash is the helper as committed. Its reference-probability fix
# postdates the recorded results, whose rows name the earlier helper they ran
# against (813c7f6af9cbb06cac9a2b8c24ab07e5), so a rerun is distinguishable
# from them.

helper <- "tools/simval/studies/dif-conditional-bootstrap-extended.R"
helper_md5 <- "4624788f06cc401781b2cc774a246b8e"
if (!identical(unname(tools::md5sum(helper)), helper_md5))
  stop("extended-study helper hash does not match the committed version")
ex <- parse(file = helper)
for (i in seq_len(19L)) eval(ex[[i]], envir = environment())

NREP <- as.integer(Sys.getenv("SV_REPS", "300"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "1")))
MIN_ANCHORS <- 3L

anchor_state <- function(fit, factors, model, anchors) {
  X <- fit$X
  af <- fit_scenario(X[, anchors, drop = FALSE], factors, model)
  if (is.null(af) || !isTRUE(af$est$converged)) return(NULL)
  loc_full <- vapply(fit$tau_list[anchors], mean, 0)
  loc_anchor <- vapply(af$tau_list[anchors], mean, 0)
  shift <- mean(loc_anchor - loc_full)
  tau <- lapply(fit$tau_list, function(x) x + shift)
  tau[anchors] <- af$tau_list[anchors]
  pe <- .person_estimates(X[, anchors, drop = FALSE], af$tau_list[anchors])
  mo <- .moment_arrays(pe$theta, tau)
  Z <- (X - mo$E) / sqrt(mo$V)
  colnames(Z) <- colnames(X)
  list(tau = tau, person = pe, residuals = Z, anchor_fit = af)
}

stats_from_state <- function(fit, factors, state, only_item = NULL) {
  fnames <- names(factors)
  safe <- paste0("f", seq_along(fnames))
  cells <- .factor_cells(factors, sep = ".")
  ng <- .dif_n_groups(fit, cells)
  ci <- factor(.class_intervals(state$person$theta, state$person$extreme, ng))
  term_labels <- c(safe, "ci", paste0(safe, ":ci"))
  items <- if (is.null(only_item)) colnames(fit$X) else only_item
  out <- list()
  for (it in items) {
    d <- data.frame(z = state$residuals[, it], ci = ci)
    for (j in seq_along(safe)) d[[safe[j]]] <- factor(factors[[fnames[j]]])
    d <- d[stats::complete.cases(d), , drop = FALSE]
    ft <- .dif_type2(d, term_labels, variance = "hc3", robust_terms = safe)
    if (is.null(ft)) return(NULL)
    for (j in seq_along(safe)) {
      u <- ft[ft$term == safe[j], , drop = FALSE]
      n <- ft[ft$term == paste0(safe[j], ":ci"), , drop = FALSE]
      if (nrow(u) != 1L || nrow(n) != 1L) return(NULL)
      out[[length(out) + 1L]] <- data.frame(
        key = c(paste(it, fnames[j], "uniform", sep = "|"),
                paste(it, fnames[j], "nonuniform", sep = "|")),
        item = it, term = fnames[j], kind = c("uniform", "nonuniform"),
        F = c(u$F_value, n$F_value), p = c(u$p, n$p))
    }
  }
  ans <- do.call(rbind, out)
  ans$p_adj <- stats::p.adjust(ans$p, "holm")
  ans
}

stats_with_anchors <- function(fit, factors, model, anchors,
                               only_item = NULL) {
  st <- anchor_state(fit, factors, model, anchors)
  if (is.null(st)) return(NULL)
  stats_from_state(fit, factors, st, only_item)
}

full_stats <- function(fit) dif_statistics(fit, "main")

loo_stats <- function(fit, factors, model) {
  items <- colnames(fit$X)
  z <- lapply(items, function(it) stats_with_anchors(
    fit, factors, model, setdiff(items, it), only_item = it))
  if (any(vapply(z, is.null, TRUE))) return(NULL)
  out <- do.call(rbind, z)
  out$p_adj <- stats::p.adjust(out$p, "holm")
  out
}

iterative_stats <- function(fit, factors, model, max_rounds = 4L) {
  items <- colnames(fit$X)
  cur <- full_stats(fit)
  if (is.null(cur)) return(NULL)
  excluded <- character(0)
  first <- unique(cur$item[cur$p_adj < 0.05])
  for (r in seq_len(max_rounds)) {
    flagged <- unique(cur$item[cur$p_adj < 0.05])
    add <- setdiff(flagged, excluded)
    if (!length(add)) break
    proposed <- union(excluded, add)
    anchors <- setdiff(items, proposed)
    if (length(anchors) < MIN_ANCHORS) break
    nxt <- stats_with_anchors(fit, factors, model, anchors)
    if (is.null(nxt)) return(NULL)
    excluded <- proposed
    cur <- nxt
  }
  list(stats = cur, excluded = excluded, first = first,
       n_anchors = length(setdiff(items, excluded)))
}

method_decisions <- function(s, sc) {
  target <- s$key == paste("I3", "group", sc$target_kind, sep = "|")
  unaffected <- s$item != "I3"
  other <- if (sc$design == "multifactor") s$term == "region" else
    rep(FALSE, nrow(s))
  c(fwer = any(s$p_adj < 0.05),
    target = if (any(target)) any(s$p_adj[target] < 0.05) else NA,
    unaffected_fwer = any(s$p_adj[unaffected] < 0.05),
    other_factor_fwer = if (any(other)) any(s$p_adj[other] < 0.05) else NA,
    raw_uniform = mean(s$p[s$kind == "uniform"] < 0.05),
    raw_nonuniform = mean(s$p[s$kind == "nonuniform"] < 0.05))
}

one_purification <- function(seed, sc) {
  set.seed(seed)
  dat <- simulate_scenario(sc)
  fit <- fit_scenario(dat$X, dat$factors, sc$model)
  if (is.null(fit) || !isTRUE(fit$est$converged)) return(NULL)
  items <- colnames(fit$X)
  fs <- full_stats(fit)
  ls <- loo_stats(fit, dat$factors, sc$model)
  os <- stats_with_anchors(fit, dat$factors, sc$model, setdiff(items, "I3"))
  is <- iterative_stats(fit, dat$factors, sc$model)
  if (any(vapply(list(fs, ls, os, is), is.null, TRUE))) return(NULL)
  ans <- c(full = method_decisions(fs, sc),
           loo = method_decisions(ls, sc),
           oracle = method_decisions(os, sc),
           iterative = method_decisions(is$stats, sc),
           iterative_target_first = "I3" %in% is$first,
           iterative_unaffected_first = any(is$first != "I3"),
           iterative_target_excluded = "I3" %in% is$excluded,
           iterative_unaffected_excluded = any(is$excluded != "I3"),
           iterative_n_anchors = is$n_anchors)
  ans
}

run_purification <- function(sc) {
  seeds <- sample.int(.Machine$integer.max, NREP)
  z <- parallel::mclapply(seeds, one_purification, sc = sc, mc.cores = NCORE,
                          mc.preschedule = TRUE)
  refused <- sum(vapply(z, is.null, TRUE))
  z <- z[!vapply(z, is.null, TRUE)]
  if (!length(z)) stop("no replicate completed in ", sc$name)
  z <- do.call(rbind, z)
  n <- nrow(z)
  null <- sc$uniform == 0 && sc$nonuniform == 0
  effect <- max(abs(sc$uniform), abs(sc$nonuniform))
  mk <- function(quantity, col, field, note) {
    v <- z[, col]
    args <- list(study = "dif-score-purification", scenario = sc$name,
      quantity = quantity, n_reps = n, n_attempted = NREP,
      n_refused = refused, n_nonconv = 0L, effect = effect,
      notes = paste0(sc$model, "; ", sc$design, "; anchor calibration ",
        "and person measures refitted; helper md5 ", helper_md5, "; ", note))
    args[[field]] <- mean(v, na.rm = TRUE)
    if (field == "type1") args$mc_override <- list(
      type1 = stats::sd(v, na.rm = TRUE) / sqrt(sum(is.finite(v))))
    do.call(sv_row, args)
  }
  rows <- list()
  for (method in c("full", "loo", "oracle", "iterative")) {
    if (null) {
      rows[[length(rows) + 1L]] <- mk(
        paste(method, "Holm FWER"), paste0(method, ".fwer"), "familywise",
        if (method == "full") "public analysis" else paste(method, "matching"))
      rows[[length(rows) + 1L]] <- mk(
        paste(method, "uniform item-wise Type I"),
        paste0(method, ".raw_uniform"), "type1", paste(method, "matching"))
      rows[[length(rows) + 1L]] <- mk(
        paste(method, "non-uniform item-wise Type I"),
        paste0(method, ".raw_nonuniform"), "type1", paste(method, "matching"))
    } else {
      rows[[length(rows) + 1L]] <- mk(
        paste(method, "target power"), paste0(method, ".target"), "power",
        paste("I3", sc$target_kind, "group term"))
      rows[[length(rows) + 1L]] <- mk(
        paste(method, "unaffected-item FWER"),
        paste0(method, ".unaffected_fwer"), "familywise", "I1-I2 and I4-I6")
      if (sc$design == "multifactor") rows[[length(rows) + 1L]] <- mk(
        paste(method, "non-target-factor FWER"),
        paste0(method, ".other_factor_fwer"), "familywise", "all region terms")
    }
  }
  rows[[length(rows) + 1L]] <- mk("iterative target selected first",
    "iterative_target_first", "power", "selection diagnostic")
  rows[[length(rows) + 1L]] <- mk("iterative unaffected item selected first",
    "iterative_unaffected_first", "familywise", "selection diagnostic")
  rows[[length(rows) + 1L]] <- mk("iterative target ultimately excluded",
    "iterative_target_excluded", "power", "selection diagnostic")
  rows[[length(rows) + 1L]] <- mk("iterative unaffected item ultimately excluded",
    "iterative_unaffected_excluded", "familywise", "selection diagnostic")
  # The mean number of retained anchors is a continuous diagnostic. Store it
  # in mean_se so the common result schema remains rectangular.
  rows[[length(rows) + 1L]] <- sv_row(
    study = "dif-score-purification", scenario = sc$name,
    quantity = "iterative mean anchors retained", n_reps = n,
    n_attempted = NREP, n_refused = refused, n_nonconv = 0L,
    effect = effect, mean_se = mean(z[, "iterative_n_anchors"]),
    se_ratio = NA_real_, notes = paste0("mean_se stores the mean count; ",
      "helper md5 ", helper_md5))
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

set.seed(8.57e7)
rows <- list()
for (i in seq_along(scenarios)) {
  message(sprintf("[%d/%d] %s", i, length(scenarios), scenarios[[i]]$name))
  rows[[i]] <- run_purification(scenarios[[i]])
  sv_write(do.call(rbind, rows), "dif-score-purification")
}
