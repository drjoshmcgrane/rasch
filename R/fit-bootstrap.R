# rasch :: a null distribution for the item fit statistics
# ===========================================================================
# Every item fit statistic in this package is computed at ESTIMATED person
# locations, and referred to a distribution derived as though those locations
# were known. Each pays for that differently, and all of them pay more as the
# sample grows.
#
# The class-interval chi-square forms its intervals by selecting on the
# estimates, so within an interval the true locations regress toward the mean
# relative to the estimates that selected them and the expected interval mean
# is biased. The bias does not shrink as the sample grows while the interval
# sizes do, so the statistic grows about linearly in N against a fixed
# reference: at 8 items a correctly fitting item is rejected 9% of the time
# at 250 persons and 100% at 4,000, and the whole-test total is rejecting
# every correctly fitting instrument by 1,000.
#
# The fit residual and the Wilson-Hilferty standardisations fail in a
# different direction. Under a true model the fit residual's null SD runs
# 0.71 at 250 persons to 1.00 at 4,000 with a mean drifting from -0.12 to
# -0.37, so the conventional +/-2.5 cut is far too lenient at small samples
# and only reaches its nominal meaning at a few thousand persons. The
# standardised statistics are worse: infit z beyond 1.96 flags 12% of
# correctly fitting items at 250 persons and 69% at 4,000. Only the mean
# squares themselves are stable, sitting at 1.06 at every sample size -- it
# is the standardisation that drifts, because the spread it divides by
# shrinks as N grows.
#
# Rescaling the chi-square to a reference sample size, the remedy RUMM2030
# offers and this package carried until 1.12.1, is one constant applied
# across items: it can move how many items are flagged but not which, and at
# a reference of 500 it silenced a planted 2.5-slope item in every one of 120
# replicates. Analysing a subsample instead inherits whatever miscalibration
# holds at the subsample's size. Neither addresses the reference
# distribution, which is the thing that is wrong.
#
# The parametric bootstrap replaces it, once, for all of them. Data are
# generated from the fitted model and each replicate travels the same road as
# the observed data -- estimate the item side, estimate person locations from
# it, form class intervals on those estimates, take residuals against those
# same estimates -- so the same bias enters the null as entered the observed
# statistic and cancels in the comparison. One set of replicates serves every
# statistic, because a replicate has to refit the whole model anyway.
#
# What the correction cannot supply is power the statistic never had. A
# flatter-than-Rasch item varies less across the class intervals and so
# carries LESS of the selection bias than a fitting item does: its chi-square
# comes out smaller than a fitting item's, and no reference distribution can
# make a small statistic significant. The fit residual sees that departure
# plainly. They have complementary blind spots, which is the argument for
# calibrating both rather than choosing between them.
# ===========================================================================

.fit_boot_failure <- function(status = c("error", "nonconverged")) {
  status <- match.arg(status)
  structure(list(status = status), class = "rasch_fit_boot_failure")
}

.fit_boot_refuse <- function(..., B, B_used, B_nonconverged, B_errors) {
  stop(structure(
    class = c("rasch_fit_bootstrap_refusal", "rasch_refusal", "error",
              "condition"),
    list(message = paste0(...), call = NULL, B = B, B_used = B_used,
         B_nonconverged = B_nonconverged, B_errors = B_errors)))
}

.fit_boot_status <- function(x) {
  if (inherits(x, "rasch_fit_boot_failure")) return(x$status)
  if (is.null(x)) return("error")
  "ok"
}

.fit_boot_md5_raw <- function(x) {
  path <- tempfile("rasch-fit-signature-", fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  writeBin(x, path)
  unname(tools::md5sum(path))
}

# Current signatures hash R's platform-independent XDR serialisation. The
# prefix identifies the encoding so the validator can continue to read the
# text-dump signatures written before this format was introduced.
.fit_boot_canonical <- function(x) {
  # A formula's evaluation environment is not part of the fitted design, but
  # serialising it would hash every mutable binding in the caller's frame.
  # The formula expression and all model matrices remain in the object.
  if (inherits(x, "formula")) {
    environment(x) <- emptyenv()
    return(x)
  }
  if (is.environment(x)) return("<environment>")
  if (is.function(x)) return(list(formals = formals(x), body = body(x)))
  if (is.list(x)) {
    for (i in seq_along(x)) x[i] <- list(.fit_boot_canonical(x[[i]]))
  }
  x
}

.fit_boot_md5 <- function(x) {
  paste0("xdr3:", .fit_boot_md5_raw(serialize(
    .fit_boot_canonical(x), connection = NULL, ascii = FALSE,
    xdr = TRUE, version = 3)))
}

.fit_boot_md5_legacy_candidates <- function(x) {
  path <- tempfile("rasch-fit-signature-", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  dput(x, path, control = c("keepNA", "keepInteger", "showAttributes",
                            "hexNumeric"))
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  z <- as.integer(bytes)
  # Text connections use the host line ending. Reconstruct both historical
  # byte streams so a legacy result saved on Windows validates on Unix and
  # vice versa.
  drop_cr <- z == 13L & c(z[-1L], -1L) == 10L
  lf <- z[!drop_cr]
  crlf <- unlist(lapply(lf, function(b)
    if (b == 10L) c(13L, 10L) else b), use.names = FALSE)
  unique(c(.fit_boot_md5_raw(as.raw(lf)),
           .fit_boot_md5_raw(as.raw(crlf))))
}

.fit_boot_hash_matches <- function(stored, x) {
  if (!is.character(stored) || length(stored) != 1L || is.na(stored))
    return(FALSE)
  if (startsWith(stored, "xdr3:")) identical(stored, .fit_boot_md5(x)) else
    stored %in% .fit_boot_md5_legacy_candidates(x)
}

.fit_signature_target <- function(fit) {
  # report_document() carries already-validated results to its template as
  # attributes. They do not change the fitted model and must not make those
  # same results appear to belong to a different fit during rendering.
  for (nm in c("report_dif", "report_bootstrap", "report_dif_bootstrap",
               "report_dimensionality", "report_subtest",
               "report_invariance", "report_person_weights",
               "report_tailored"))
    attr(fit, nm) <- NULL
  fit
}

.fit_boot_signature <- function(fit) {
  fit <- .fit_signature_target(fit)
  if (inherits(fit, "rasch_btl")) {
    return(list(
      kind = "btl", m = fit$m, n_comparisons = fit$n_comparisons,
      objects = as.character(fit$objects$object),
      locations = unname(as.numeric(fit$objects$location)),
      pairs = .btl_pair_key(fit$pairs$object_a, fit$pairs$object_b),
      total_chisq = unname(as.numeric(fit$total_chisq)),
      # The result is tied to the complete fitted object, not only to its
      # displayed locations. In particular, fitted probabilities generate
      # an unordered bootstrap and thresholds plus history coefficients
      # generate an ordered one; omitting either accepted a bootstrap from a
      # different null model when its object locations happened to agree.
      fingerprint = .fit_boot_md5(fit)))
  }
  list(
    kind = "rasch", model = fit$model, dimensions = dim(fit$X),
    items = as.character(fit$items$item),
    item_chisq = unname(as.numeric(fit$items$chisq)),
    raw_scores = unname(as.numeric(fit$person$raw)),
    person_ids = as.character(fit$person$id),
    # DIF and fitted-model bootstraps use more than the response matrix and
    # threshold table: person identities determine repeated-person strata,
    # while frame maps, virtual items, discriminations and refit designs can
    # all change the analysis without changing those two objects. Hash the
    # fitted model as a whole so no active part is silently omitted again.
    fingerprint = .fit_boot_md5(fit))
}

.fit_boot_signature_matches <- function(signature, fit) {
  fit <- .fit_signature_target(fit)
  current <- .fit_boot_signature(fit)
  if (!is.list(signature) || !identical(names(signature), names(current)) ||
      !is.character(signature$fingerprint) ||
      length(signature$fingerprint) != 1L || is.na(signature$fingerprint))
    return(FALSE)
  stored_fingerprint <- signature$fingerprint
  signature$fingerprint <- NULL
  current$fingerprint <- NULL
  identical(signature, current) &&
    .fit_boot_hash_matches(stored_fingerprint, fit)
}

# The algorithm tag a stored result must carry to be read as current. The
# ability-sampling generators once re-derived an automatic class-interval
# count inside every replicate, so their nulls mixed interval counts, and
# degrees of freedom, the observed chi-square never used; those results are
# superseded. The conditional generator retains every raw score, so its
# non-extreme set, and hence its interval count, was already constant across
# replicates: those nulls are unchanged and stay current.
.fit_boot_algorithm <- function(kind, theta) {
  if (identical(kind, "btl")) "loo-maxt-1"
  else if (identical(theta, "conditional")) "loo-maxt-2" else "loo-maxt-3"
}

.new_fit_bootstrap <- function(x, fit, kind) {
  x$algorithm <- .fit_boot_algorithm(kind, x$theta)
  x$model_kind <- kind
  x$fit_signature <- .fit_boot_signature(fit)
  x <- .tag_tables(x)
  # The fitted-model fingerprint prevents a result being paired with another
  # calibration. This second fingerprint protects the result itself: reports
  # and restored app projects must not accept edited tables or accounting as
  # though they were the output of the retained replicate arrays.
  x$result_signature <- .fit_boot_md5(x)
  class(x) <- c("rasch_fit_bootstrap", "list")
  x
}

.validate_fit_bootstrap <- function(bootstrap, fit) {
  if (is.null(bootstrap)) return(invisible(NULL))
  .require_refittable_calibration(fit)
  fail <- function() stop(
    "`bootstrap` is incomplete or internally inconsistent; recompute it with fit_bootstrap()",
    call. = FALSE)
  if (!inherits(bootstrap, "rasch_fit_bootstrap"))
    stop("`bootstrap` must be a current fit_bootstrap() result", call. = FALSE)
  if (!is.character(bootstrap$result_signature) ||
      length(bootstrap$result_signature) != 1L ||
      is.na(bootstrap$result_signature)) fail()
  unsigned <- unclass(bootstrap)
  unsigned$result_signature <- NULL
  if (!.fit_boot_hash_matches(bootstrap$result_signature, unsigned)) fail()
  expected <- if (inherits(fit, "rasch_btl")) "btl" else "rasch"
  if (!identical(bootstrap$model_kind, expected))
    stop("`bootstrap` was computed for a different model family")
  if (!.fit_boot_signature_matches(bootstrap$fit_signature, fit))
    stop("`bootstrap` was computed from a different fitted model")

  required <- c("algorithm", "total", "replicates", "B", "B_used", "B_failed",
                "B_nonconverged", "B_errors", "minimum_usable", "theta")
  required <- c(required, if (expected == "btl")
    c("model", "pairs", "objects", "judges", "adjustment") else
    c("items", "persons", "person_adjustment"))
  if (!all(required %in% names(bootstrap)) ||
      !identical(bootstrap$algorithm,
                 .fit_boot_algorithm(expected, bootstrap$theta)) ||
      !is.list(bootstrap$total) || !is.list(bootstrap$replicates)) fail()

  whole <- function(x, lower = 0L)
    length(x) == 1L && is.numeric(x) && is.finite(x) && x == floor(x) &&
      x >= lower && x <= .Machine$integer.max
  count_names <- c("B", "B_used", "B_failed", "B_nonconverged",
                   "B_errors", "minimum_usable")
  count_lower <- c(B = 1L, B_used = 1L, B_failed = 0L,
                   B_nonconverged = 0L, B_errors = 0L,
                   minimum_usable = 1L)
  if (any(!vapply(count_names, function(nm)
    whole(bootstrap[[nm]], count_lower[[nm]]), logical(1))) ||
      bootstrap$B_used + bootstrap$B_failed != bootstrap$B ||
      bootstrap$B_nonconverged + bootstrap$B_errors != bootstrap$B_failed ||
      bootstrap$minimum_usable != .fit_min_boot_success(bootstrap$B) ||
      bootstrap$B_used < bootstrap$minimum_usable) fail()

  same_num <- function(x, y) is.numeric(x) && length(x) == length(y) &&
    !any(is.nan(x) | is.infinite(x)) &&
    isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 0,
                     check.attributes = FALSE))
  same_frame <- function(x, y) is.data.frame(x) && nrow(x) == nrow(y) &&
    all(names(y) %in% names(x)) && isTRUE(all.equal(
      as.data.frame(x[names(y)]), as.data.frame(y), tolerance = 0,
      check.attributes = FALSE))
  valid_matrix <- function(x, nr, nc)
    is.matrix(x) && is.numeric(x) && identical(dim(x), c(nr, nc)) &&
      !any(is.nan(x) | is.infinite(x))
  same_meta <- function(x, y) is.list(x) &&
    identical(names(x), c("family_n", "adjusted_n", "joint_replicates")) &&
    all(vapply(names(y), function(nm)
      same_num(x[[nm]], y[[nm]]), logical(1)))
  p_ok <- function(x) is.numeric(x) && !any(is.nan(x) | is.infinite(x)) &&
    all(is.na(x) | (x >= 0 & x <= 1))
  min_success <- as.integer(bootstrap$minimum_usable)
  B_used <- as.integer(bootstrap$B_used)

  validate_table <- function(tab, base, reps, stats, sides) {
    if (!same_frame(tab, base) || !is.list(reps) ||
        !identical(names(reps), stats)) fail()
    meta <- list()
    for (st in stats) {
      if (!st %in% names(base) || !valid_matrix(reps[[st]], B_used,
                                                nrow(base))) fail()
      mode <- if (st == "chisq") "raw" else if (st == "fit_resid")
        "centred" else "studentised"
      trans <- if (st %in% c("infit_ms", "outfit_ms")) log else identity
      z <- .boot_maxt(base[[st]], reps[[st]], sides[[st]], min_success,
                      mode = mode, transform = trans)
      fields <- c(paste0(st, "_p_boot"), paste0(st, "_p_boot_adj"),
                  paste0("n_boot_", st))
      if (!all(fields %in% names(tab)) ||
          !same_num(tab[[fields[1L]]], z$p) ||
          !same_num(tab[[fields[2L]]], z$p_adj) ||
          !same_num(tab[[fields[3L]]], z$n_boot) ||
          !p_ok(tab[[fields[1L]]]) || !p_ok(tab[[fields[2L]]])) fail()
      meta[[st]] <- list(family_n = z$family_n,
                         adjusted_n = z$adjusted_n,
                         joint_replicates = z$family_boot)
    }
    meta
  }

  if (expected == "btl") {
    if (!identical(bootstrap$model, "paired comparisons") ||
        !is.null(bootstrap$theta) || !is.list(bootstrap$adjustment) ||
        !all(c("total_chisq", "pairs", "objects", "judges") %in%
             names(bootstrap$replicates))) fail()
    pair_stats <- "chisq"
    object_stats <- c("fit_resid", "infit_ms", "outfit_ms")
    pair_meta <- validate_table(
      bootstrap$pairs, fit$pairs, bootstrap$replicates$pairs, pair_stats,
      c(chisq = "upper"))
    object_meta <- validate_table(
      bootstrap$objects, fit$objects, bootstrap$replicates$objects,
      object_stats,
      c(fit_resid = "two", infit_ms = "two", outfit_ms = "two"))
    if (is.null(fit$judges)) {
      if (!is.null(bootstrap$judges) ||
          !is.null(bootstrap$replicates$judges) ||
          !is.null(bootstrap$adjustment$judges)) fail()
      judge_meta <- NULL
    } else {
      judge_meta <- validate_table(
        bootstrap$judges, fit$judges, bootstrap$replicates$judges,
        object_stats,
        c(fit_resid = "two", infit_ms = "two", outfit_ms = "two"))
    }
    expected_method <- paste("single-step maximum-statistic bootstrap within each",
                             "statistic and parameter family")
    if (!identical(bootstrap$adjustment$method, expected_method) ||
        !identical(names(bootstrap$adjustment),
                   c("method", "pairs", "objects", "judges")) ||
        !identical(names(bootstrap$adjustment$pairs), names(pair_meta)) ||
        !identical(names(bootstrap$adjustment$objects), names(object_meta)) ||
        !all(vapply(names(pair_meta), function(st)
          same_meta(bootstrap$adjustment$pairs[[st]], pair_meta[[st]]),
          logical(1))) ||
        !all(vapply(names(object_meta), function(st)
          same_meta(bootstrap$adjustment$objects[[st]], object_meta[[st]]),
          logical(1))) ||
        (!is.null(judge_meta) &&
         (!identical(names(bootstrap$adjustment$judges), names(judge_meta)) ||
          !all(vapply(names(judge_meta), function(st)
            same_meta(bootstrap$adjustment$judges[[st]], judge_meta[[st]]),
            logical(1)))))) fail()
    tot <- bootstrap$replicates$total_chisq
    if (!is.numeric(tot) || length(tot) != B_used ||
        any(is.nan(tot) | is.infinite(tot))) fail()
    n_tot <- sum(is.finite(tot))
    p_tot <- if (n_tot >= min_success)
      .boot_p(fit$total_chisq, tot, "upper") else NA_real_
    expected_total <- list(chisq = fit$total_chisq, df = fit$total_df,
                           p = fit$total_p, chisq_p_boot = p_tot,
                           n_boot = n_tot)
    if (!identical(names(bootstrap$total), names(expected_total)) ||
        !all(vapply(names(expected_total), function(nm)
          same_num(bootstrap$total[[nm]], expected_total[[nm]]),
          logical(1))) || !p_ok(bootstrap$total$chisq_p_boot)) fail()
    return(invisible(NULL))
  }

  if (!is.data.frame(bootstrap$items) || !is.list(bootstrap$replicates) ||
      !identical(names(bootstrap$replicates), names(.FIT_STATS)) ||
      nrow(bootstrap$items) != nrow(fit$items) ||
      !identical(as.character(bootstrap$items$item),
                 as.character(fit$items$item)) ||
      !is.character(bootstrap$theta) || length(bootstrap$theta) != 1L ||
      is.na(bootstrap$theta) ||
      !bootstrap$theta %in% c("conditional", "resample", "fixed", "normal"))
    fail()
  L <- nrow(fit$items)
  for (st in names(.FIT_STATS)) {
    M <- bootstrap$replicates[[st]]
    fields <- c(st, paste0(st, "_p_boot"), paste0(st, "_p_boot_adj"),
                paste0("n_boot_", st))
    if (!st %in% names(fit$items) || !all(fields %in% names(bootstrap$items)) ||
        !valid_matrix(M, B_used, L) ||
        !same_num(bootstrap$items[[st]], fit$items[[st]])) fail()
    n_stat <- colSums(is.finite(M))
    p <- vapply(seq_len(L), function(i) {
      if (n_stat[i] < min_success) return(NA_real_)
      .boot_p(fit$items[[st]][i], M[, i], .FIT_STATS[[st]])
    }, numeric(1))
    adj <- rep(NA_real_, L)
    usable <- is.finite(p)
    adj[usable] <- stats::p.adjust(p[usable], method = "holm", n = L)
    if (!same_num(bootstrap$items[[fields[2L]]], p) ||
        !same_num(bootstrap$items[[fields[3L]]], adj) ||
        !same_num(bootstrap$items[[fields[4L]]], n_stat) ||
        !p_ok(bootstrap$items[[fields[2L]]]) ||
        !p_ok(bootstrap$items[[fields[3L]]])) fail()
  }
  if (!"n_boot" %in% names(bootstrap$items) ||
      !same_num(bootstrap$items$n_boot, bootstrap$items$n_boot_chisq)) fail()

  if (identical(bootstrap$theta, "conditional")) {
    keep_cols <- intersect(c("id", "raw", "theta", "se"), names(fit$person))
    if (!is.data.frame(bootstrap$persons) ||
        !same_frame(bootstrap$persons, fit$person[, keep_cols, drop = FALSE]) ||
        !is.list(bootstrap$person_adjustment) ||
        !identical(bootstrap$person_adjustment$method,
                   "single-step maximum-statistic bootstrap") ||
        !is.list(bootstrap$person_adjustment$statistics) ||
        !identical(names(bootstrap$person_adjustment$statistics),
                   names(.PERSON_FIT_STATS))) fail()
    for (st in names(.PERSON_FIT_STATS)) {
      # Person nulls are not part of the item matrices. They are retained in
      # the public result only through their derived tables. The result
      # fingerprint protects their probabilities; these checks also pin the
      # copied observations and the visible maxT accounting.
      fields <- c(st, paste0(st, "_p_boot"), paste0(st, "_p_boot_adj"),
                  paste0("n_boot_", st))
      if (!all(fields %in% names(bootstrap$persons)) ||
          !same_num(bootstrap$persons[[st]], fit$person[[st]]) ||
          !p_ok(bootstrap$persons[[fields[2L]]]) ||
          !p_ok(bootstrap$persons[[fields[3L]]]) ||
          !is.numeric(bootstrap$persons[[fields[4L]]]) ||
          any(bootstrap$persons[[fields[4L]]] < 0 |
              bootstrap$persons[[fields[4L]]] > B_used)) fail()
      meta <- bootstrap$person_adjustment$statistics[[st]]
      if (!is.list(meta) ||
          !all(c("family_n", "adjusted_n", "joint_replicates") %in%
               names(meta)) ||
          !same_num(meta$family_n, sum(is.finite(fit$person[[st]]))) ||
          !same_num(meta$adjusted_n,
                    sum(is.finite(bootstrap$persons[[fields[3L]]]))) ||
          !whole(meta$joint_replicates, 0L) ||
          meta$joint_replicates > B_used) fail()
    }
  } else if (!is.null(bootstrap$persons) ||
             !is.null(bootstrap$person_adjustment)) fail()

  M <- bootstrap$replicates
  ok <- is.finite(fit$items$chisq)
  tot_rep <- rep(NA_real_, B_used)
  if (any(ok)) {
    complete <- rowSums(is.finite(M$chisq[, ok, drop = FALSE])) == sum(ok)
    tot_rep[complete] <- rowSums(M$chisq[complete, ok, drop = FALSE])
  }
  fr_ok <- is.finite(fit$items$fit_resid)
  fr_mean <- fr_sd <- rep(NA_real_, B_used)
  if (any(fr_ok)) {
    complete <- rowSums(is.finite(M$fit_resid[, fr_ok, drop = FALSE])) ==
      sum(fr_ok)
    fr_mean[complete] <- rowMeans(M$fit_resid[complete, fr_ok, drop = FALSE])
    if (sum(fr_ok) >= 2L)
      fr_sd[complete] <- apply(M$fit_resid[complete, fr_ok, drop = FALSE],
                               1L, stats::sd)
  }
  total_p <- function(obs, null, side) if (sum(is.finite(null)) >= min_success)
    .boot_p(obs, null, side) else NA_real_
  expected_total <- list(
    chisq = fit$total_chisq, df = fit$total_df, p = fit$total_chisq_p,
    chisq_p_boot = total_p(fit$total_chisq, tot_rep, "upper"),
    fit_resid_mean = fit$item_fit_summary$mean,
    fit_resid_mean_p_boot = total_p(fit$item_fit_summary$mean, fr_mean, "two"),
    fit_resid_sd = fit$item_fit_summary$sd,
    fit_resid_sd_p_boot = total_p(fit$item_fit_summary$sd, fr_sd, "two"),
    n_boot = sum(is.finite(tot_rep)),
    n_boot_fit_resid_mean = sum(is.finite(fr_mean)),
    n_boot_fit_resid_sd = sum(is.finite(fr_sd)))
  if (!identical(names(bootstrap$total), names(expected_total)) ||
      !all(vapply(names(expected_total), function(nm)
        same_num(bootstrap$total[[nm]], expected_total[[nm]]), logical(1))) ||
      !all(vapply(c("chisq_p_boot", "fit_resid_mean_p_boot",
                    "fit_resid_sd_p_boot"), function(nm)
        p_ok(bootstrap$total[[nm]]), logical(1)))) fail()
  invisible(NULL)
}

# The estimation chain a replicate shares with the observed fit: the item
# calibration, the person estimates, and the standardised residuals. A
# failure object stands in for a replicate that could not be estimated.
.fit_refit_residuals <- function(X, model, anchors, expected_m, maxit, tol) {
  est <- tryCatch(pcml(X, model = model, anchors = anchors,
                       maxit = maxit, tol = tol),
                  error = function(e) .fit_boot_failure("error"))
  if (inherits(est, "rasch_fit_boot_failure")) return(est)
  if (!isTRUE(est$converged)) return(.fit_boot_failure("nonconverged"))
  L <- ncol(X)
  tau_list <- lapply(seq_len(L), function(i) est$thr$tau[est$thr$item == i])
  if (!all(vapply(tau_list, function(t) length(t) && all(is.finite(t)), TRUE)))
    return(.fit_boot_failure("error"))
  # A sparse replicate can leave its highest category unvisited. pcml() then
  # fits a shorter item without error, but that is a different model and its
  # residual statistic cannot enter the null reference for the observed
  # scale. Treat it like any other failed same-model refit.
  if (!identical(vapply(tau_list, length, integer(1)),
                 as.integer(expected_m)))
    return(.fit_boot_failure("error"))
  person <- .person_estimates(X, tau_list, disc = 1)
  if (all(is.na(person$theta))) return(.fit_boot_failure("error"))
  mo <- .moment_arrays(person$theta, tau_list, disc = rep(1, L))
  Z <- (X - mo$E) / sqrt(mo$V)
  colnames(Z) <- colnames(X)
  list(est = est, tau_list = tau_list, person = person, moments = mo, Z = Z)
}

# Everything the item fit statistics need from a response matrix and nothing
# else: the person tables, score curves and targeting that rasch() also builds
# are not used here and would multiply the cost by the number of replicates.
# Repeat the requested allocation rule: NULL selects intervals automatically
# (per item with missing responses), while a supplied count stays fixed.
# fixed_groups, when supplied, holds the interval counts the observed fit
# realised under the automatic rule, so a replicate is scored on the same
# intervals as the statistic it forms a null for.
.fit_refit <- function(X, model, n_groups, anchors, expected_m, maxit, tol,
                       fixed_groups = NULL) {
  r <- .fit_refit_residuals(X, model, anchors, expected_m, maxit, tol)
  if (inherits(r, "rasch_fit_boot_failure")) return(r)
  est <- r$est; tau_list <- r$tau_list; person <- r$person
  mo <- r$moments; Z <- r$Z
  L <- ncol(X)
  m_i <- vapply(tau_list, length, 1L)
  item_extreme <- vapply(seq_len(L), function(j) {
    col <- X[!person$extreme, j]
    tot <- sum(col, na.rm = TRUE); nn <- sum(!is.na(col))
    nn == 0L || tot == 0 || tot == nn * m_i[j]
  }, logical(1))
  ci <- .class_intervals(person$theta, person$extreme,
                         fixed_groups$common %||% n_groups)
  ng_item <- fixed_groups$item
  ci_list <- if (!anyNA(X)) NULL else if (is.null(ng_item))
    .class_intervals_by_item(X, person$theta, person$extreme, n_groups) else
      lapply(seq_len(L), function(i)
        .class_intervals_by_item(X[, i, drop = FALSE], person$theta,
                                 person$extreme, ng_item[i])[[1L]])
  it <- .item_trait(X, mo, ci, ci_list = ci_list)
  ifit <- .item_fit(X, Z, mo, extreme = person$extreme)
  pfit <- .person_fit(X, Z, mo, item_extreme = item_extreme)
  n_par <- if (is.null(est$n_parameters)) nrow(est$thr) - 1L else est$n_parameters
  rf <- .fitres(Z, mo, person$extreme, n_par)
  list(chisq = it$chisq, fit_resid = rf$items$fit_resid,
       infit_ms = ifit$infit_ms, outfit_ms = ifit$outfit_ms,
       infit_z = ifit$infit_z, outfit_z = ifit$outfit_z,
       persons = list(fit_resid = rf$persons$fit_resid,
                      infit_ms = pfit$infit_ms,
                      outfit_ms = pfit$outfit_ms,
                      infit_z = pfit$infit_z,
                      outfit_z = pfit$outfit_z))
}

# The interval counts the observed item-trait chi-squares were computed on.
# The automatic rule reads the non-extreme sample size, which the
# ability-sampling generators move from replicate to replicate: re-deriving
# the rule inside each replicate builds the null from chi-squares carrying a
# different number of intervals, and different degrees of freedom, from the
# observed statistic. Resolve the rule once, against the fit itself.
.fit_boot_realised_groups <- function(fit) {
  common <- fit$n_groups
  if (!(length(common) == 1L && is.finite(common) && common >= 2L))
    return(NULL)
  item <- if (length(fit$ci_item) != ncol(fit$X)) NULL else
    vapply(fit$ci_item, function(g) {
      g <- g[!is.na(g)]
      if (length(g)) as.integer(max(g)) else 1L
    }, 1L)
  list(common = as.integer(common), item = item)
}

# The statistics this bootstrap calibrates, and how each is read. The
# chi-square is a discrepancy and only its upper tail is misfit; the others
# depart in two directions -- above for a flatter item than the model
# predicts, below for a steeper one -- and both readings are misfit, so their
# p-values are equal-tailed.
.FIT_STATS <- c(chisq = "upper", fit_resid = "two", infit_ms = "two",
                outfit_ms = "two", infit_z = "two", outfit_z = "two")
.PERSON_FIT_STATS <- c(fit_resid = "two", infit_ms = "two",
                       outfit_ms = "two", infit_z = "two", outfit_z = "two")

# Bootstrap p-value for one statistic against its own replicated null.
# (1 + r)/(1 + B) never returns zero; the two-sided form doubles the smaller
# tail and is capped, so its resolution is half the one-sided form's.
.boot_p <- function(obs, null, side) {
  null <- null[is.finite(null)]
  n <- length(null)
  if (!is.finite(obs) || n < 1L) return(NA_real_)
  up <- (1 + sum(null >= obs)) / (1 + n)
  if (side == "upper") return(up)
  lo <- (1 + sum(null <= obs)) / (1 + n)
  min(1, 2 * min(up, lo))
}

# A requested exploratory run below 30 draws may still be useful at its very
# coarse Monte Carlo resolution, but it must retain a majority. At 30 or
# above, at least 30 and 90% of the requested refits must succeed. Failed
# refits need not be a random subset of the null, particularly when a sparse
# polytomous replicate loses an interior category, so a more depleted null is
# withheld rather than interpreted conditionally on being estimable.
.fit_min_boot_success <- function(B) {
  if (B < 30L) return(as.integer(floor(B / 2) + 1L))
  as.integer(max(30L, ceiling(0.9 * B)))
}

# Single-step maximum-statistic adjustment from the joint bootstrap null. The
# centring and scale are statistic-specific, while each replicate retains the
# dependence among the quantities. Unlike an empirical minimum-p reference,
# maxT does not collapse to the endpoint when the family is larger than B.
.boot_maxt <- function(obs, null, side, min_success,
                       mode = c("studentised", "centred", "raw"),
                       transform = identity) {
  mode <- match.arg(mode)
  if (!is.matrix(null)) null <- as.matrix(null)
  n <- ncol(null)
  n_stat <- colSums(is.finite(null))
  raw <- vapply(seq_len(n), function(i) {
    if (n_stat[i] < min_success) return(NA_real_)
    .boot_p(obs[i], null[, i], side)
  }, 0)
  family <- is.finite(obs)
  adj <- rep(NA_real_, n)
  family_n <- sum(family)
  if (!family_n)
    return(list(p = raw, p_adj = adj, n_boot = n_stat,
                family_n = 0L, family_boot = 0L, adjusted_n = 0L))

  V <- transform(null[, family, drop = FALSE])
  obs_e <- transform(obs[family])
  # complete.cases() accepts +/-Inf. That is unsuitable after transforming a
  # statistic: for example log(0) is -Inf and would poison every centring,
  # scale and row maximum in the joint reference. A usable maxT row must be
  # finite for every member of the declared family.
  complete <- rowSums(is.finite(V)) == ncol(V)
  family_boot <- sum(complete)
  # The family is declared by the finite observed statistics, not by which
  # null columns happen to survive. A partial maxT family would make every
  # remaining adjusted probability depend on an estimation failure.
  # Studentisation needs at least two training degrees of freedom. Tiny
  # exploratory runs still return marginal probabilities, but cannot support
  # an adjusted family.
  enough_for_mode <- mode != "studentised" || family_boot >= 3L
  if (all(is.finite(raw[family])) && all(is.finite(obs_e)) &&
      family_boot >= min_success && enough_for_mode) {
    V <- V[complete, , drop = FALSE]
    B <- nrow(V); K <- ncol(V)
    centre <- if (mode == "raw") rep(0, K) else colMeans(V)
    spread <- apply(V, 2L, stats::sd)
    scale <- if (mode == "studentised") spread else rep(1, K)
    stable <- is.finite(scale) & scale > sqrt(.Machine$double.eps)

    # An observed statistic is external to the bootstrap null. Each null row
    # must therefore also be standardised externally: using that row in its
    # own mean and SD shrinks its extremes and makes studentised maxT tests
    # anti-conservative. Raw statistics need no estimated nuisance quantities;
    # centred and studentised statistics use leave-one-out estimates.
    Z <- matrix(0, B, K)
    if (mode == "raw") {
      Z <- V
    } else {
      for (i in seq_len(B)) {
        training <- V[-i, , drop = FALSE]
        centre_i <- colMeans(training)
        if (mode == "centred") {
          Z[i, ] <- V[i, ] - centre_i
        } else {
          scale_i <- apply(training, 2L, stats::sd)
          stable_i <- is.finite(scale_i) &
            scale_i > sqrt(.Machine$double.eps)
          Z[i, stable_i] <- (V[i, stable_i] - centre_i[stable_i]) /
            scale_i[stable_i]
          if (any(!stable_i)) {
            tol_i <- sqrt(.Machine$double.eps) *
              pmax(1, abs(centre_i[!stable_i]))
            delta_i <- V[i, !stable_i] - centre_i[!stable_i]
            Z[i, !stable_i] <- ifelse(delta_i > tol_i, Inf,
              ifelse(delta_i < -tol_i, -Inf, 0))
          }
        }
      }
    }
    z_obs <- if (mode == "raw") obs_e else obs_e - centre
    if (mode == "studentised") {
      z_obs[stable] <- z_obs[stable] / scale[stable]
      if (any(!stable)) {
        tol <- sqrt(.Machine$double.eps) * pmax(1, abs(centre[!stable]))
        delta <- obs_e[!stable] - centre[!stable]
        z_obs[!stable] <- ifelse(delta > tol, Inf,
          ifelse(delta < -tol, -Inf, 0))
      }
    }
    if (side == "two") {
      Z <- abs(Z)
      z_obs <- abs(z_obs)
    }
    max_null <- apply(Z, 1L, max)
    a <- vapply(z_obs, function(z)
      (1 + sum(max_null >= z)) / (1 + length(max_null)), 0)
    pos <- which(family)
    adj[pos] <- pmax(raw[pos], a)
  }
  list(p = raw, p_adj = adj, n_boot = n_stat,
       family_n = family_n, family_boot = family_boot,
       adjusted_n = sum(is.finite(adj)))
}

.btl_pair_key <- function(a, b) {
  lo <- ifelse(a <= b, a, b)
  hi <- ifelse(a <= b, b, a)
  paste0(nchar(lo, type = "bytes"), ":", lo,
         nchar(hi, type = "bytes"), ":", hi)
}

.btl_boot_data <- function(fit) {
  cmp <- fit$comparisons
  spec <- fit$refit_spec %||% list()
  if (is.null(cmp) || !all(c("object_a", "object_b", "weight") %in% names(cmp)))
    .refuse("the paired-comparison fit does not carry the comparison design ",
            "needed for a fit bootstrap; refit it with the current package version")
  if (any(!is.finite(cmp$weight)) || any(cmp$weight <= 0) ||
      any(cmp$weight != floor(cmp$weight)))
    .refuse("the paired-comparison fit bootstrap requires whole positive ",
            "comparison counts; a fit using half-weighted ties cannot be generated ",
            "as independent comparisons. Code ties as an ordered middle category")
  if (any(cmp$weight > .Machine$integer.max))
    .refuse("the paired-comparison fit contains a compressed row whose count ",
            "exceeds the integer range supported by the multinomial bootstrap; ",
            "split that count over several otherwise identical rows before fitting")
  if (any(grepl("half a win", fit$notes %||% character(0), fixed = TRUE)))
    .refuse("the paired-comparison fit used half-weighted ties, which the fit ",
            "bootstrap cannot generate as independent comparisons. Code ties ",
            "as an ordered middle category")

  requested_effects <- intersect(c("exposure", "carry_over", "position"),
                                 names(cmp))
  fitted_effects <- if (is.null(fit$dependence)) character(0) else
    as.character(fit$dependence$effect)
  dropped_effects <- setdiff(requested_effects, fitted_effects)
  if (length(dropped_effects))
    .refuse("the paired-comparison fit dropped the following requested ",
            "history or position effect(s): ",
            paste(dropped_effects, collapse = ", "),
            ". The bootstrap cannot refit a reduced dependence model through ",
            "the public btl() specification")

  has_order <- isTRUE(spec$has_order)
  if (has_order && (any(cmp$weight != 1) || !"order" %in% names(cmp)))
    .refuse("an ordered paired-comparison bootstrap needs one row per judgement ",
            "with a unique sequence value; count-weighted ordered rows are not a ",
            "generative history")

  if (!has_order) {
    P <- fit$fitted_prob
    if (!is.matrix(P) || nrow(P) != nrow(cmp) || ncol(P) != fit$m + 1L ||
        any(!is.finite(P)) || any(abs(rowSums(P) - 1) > 1e-8))
      .refuse("the paired-comparison fit does not carry valid fitted category ",
              "probabilities; refit it with the current package version")
    rows <- vector("list", nrow(cmp))
    for (r in seq_len(nrow(cmp))) {
      z <- as.vector(stats::rmultinom(1L, size = cmp$weight[r], prob = P[r, ]))
      k <- which(z > 0L)
      rows[[r]] <- data.frame(
        object_a = cmp$object_a[r], object_b = cmp$object_b[r],
        response = k - 1L, count = z[k],
        judge = cmp$judge[r], stringsAsFactors = FALSE)
    }
    out <- do.call(rbind, rows)
    # Keep the fitted response scale even when a replicate happens not to
    # visit an extreme category. Numeric scores would make btl() infer a
    # shorter scale from that replicate and fit a different model.
    out$response <- ordered(out$response, levels = 0:fit$m)
    return(out)
  }

  # History-dependent effects must be generated in sequence. Holding the
  # observed carry-over covariate fixed would condition on the outcome being
  # tested and produce the wrong null.
  beta <- stats::setNames(fit$objects$location, fit$objects$object)
  if (anyNA(beta[cmp$object_a]) || anyNA(beta[cmp$object_b]))
    .refuse("the ordered comparison design names an object without a fitted location")
  tau <- if (is.null(fit$thresholds)) numeric(fit$m) else fit$thresholds$tau
  dep <- if (is.null(fit$dependence)) numeric(0) else
    stats::setNames(fit$dependence$estimate, fit$dependence$effect)
  out <- cmp[, c("object_a", "object_b", "judge", "order"), drop = FALSE]
  out$response <- NA_integer_
  out$count <- 1L
  cnt <- new.env(hash = TRUE, parent = emptyenv())
  tot <- new.env(hash = TRUE, parent = emptyenv())
  gets <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
  key_a <- .factor_keys(data.frame(judge = cmp$judge,
                                   object = cmp$object_a,
                                   stringsAsFactors = FALSE))
  key_b <- .factor_keys(data.frame(judge = cmp$judge,
                                   object = cmp$object_b,
                                   stringsAsFactors = FALSE))
  for (r in order(cmp$judge, cmp$order)) {
    ka <- key_a[r]; kb <- key_b[r]
    na <- gets(cnt, ka); nb <- gets(cnt, kb)
    exposure <- as.numeric(na > 0) - as.numeric(nb > 0)
    carry <- (if (na > 0) gets(tot, ka) / na else 0) -
      (if (nb > 0) gets(tot, kb) / nb else 0)
    d <- beta[[cmp$object_a[r]]] - beta[[cmp$object_b[r]]]
    if ("exposure" %in% names(dep)) d <- d + dep[["exposure"]] * exposure
    if ("carry_over" %in% names(dep)) d <- d + dep[["carry_over"]] * carry
    if ("position" %in% names(dep)) d <- d + dep[["position"]]
    eta <- d * (0:fit$m) - c(0, cumsum(tau))
    eta <- eta - max(eta)
    pr <- exp(eta) / sum(exp(eta))
    x <- sample.int(fit$m + 1L, 1L, prob = pr) - 1L
    out$response[r] <- x
    cnt[[ka]] <- na + 1; cnt[[kb]] <- nb + 1
    tot[[ka]] <- gets(tot, ka) + (2 * x / fit$m - 1)
    tot[[kb]] <- gets(tot, kb) + (2 * (fit$m - x) / fit$m - 1)
  }
  out$response <- ordered(out$response, levels = 0:fit$m)
  out
}

.btl_boot_refit <- function(fit) {
  d <- .btl_boot_data(fit)
  spec <- fit$refit_spec %||% list()
  args <- list(data = d, object_a = "object_a", object_b = "object_b",
               response = "response", count = "count",
               thresholds = spec$thresholds %||% "free",
               position = isTRUE(spec$position),
               anchors = spec$anchors,
               maxit = spec$maxit %||% 60L, tol = spec$tol %||% 1e-8,
               .object_design = spec$object_design)
  if (any(!is.na(d$judge))) args$judge <- "judge"
  if (isTRUE(spec$has_order)) args$order <- "order"
  nonconv <- FALSE
  z <- tryCatch(withCallingHandlers(
    do.call(btl, args), warning = function(w) {
      # A non-convergence warning makes that replicate unusable. Other
      # warnings are expected occasionally at sparse bootstrap boundaries;
      # the returned fit decides whether the replicate is estimable.
      if (grepl("did NOT converge", conditionMessage(w), fixed = TRUE))
        nonconv <<- TRUE
      invokeRestart("muffleWarning")
    }), error = function(e) .fit_boot_failure("error"))
  if (inherits(z, "rasch_fit_boot_failure")) return(z)
  if (nonconv || !isTRUE(z$converged))
    return(.fit_boot_failure("nonconverged"))
  wanted_effects <- if (is.null(fit$dependence)) character(0) else
    sort(as.character(fit$dependence$effect))
  fitted_effects <- if (is.null(z$dependence)) character(0) else
    sort(as.character(z$dependence$effect))
  if (!identical(wanted_effects, fitted_effects))
    return(.fit_boot_failure("error"))
  active <- if (is.null(fit$objects$extreme))
    rep(TRUE, nrow(fit$objects)) else !fit$objects$extreme
  active_objects <- fit$objects$object[active]
  if (!setequal(as.character(z$objects$object), as.character(active_objects)) ||
      any(z$objects$extreme %||% FALSE)) return(.fit_boot_failure("error"))

  align <- function(tab, key, target, stats) {
    if (is.null(tab)) return(NULL)
    j <- match(target, key)
    ans <- lapply(stats, function(st) {
      v <- rep(NA_real_, length(target)); ok <- !is.na(j)
      v[ok] <- tab[[st]][j[ok]]; v
    })
    names(ans) <- stats
    ans
  }
  pair_target <- .btl_pair_key(fit$pairs$object_a, fit$pairs$object_b)
  pair_key <- .btl_pair_key(z$pairs$object_a, z$pairs$object_b)
  list(
    total_chisq = z$total_chisq,
    pairs = align(z$pairs, pair_key, pair_target, c("chisq")),
    objects = align(z$objects, z$objects$object, fit$objects$object,
                    c("fit_resid", "infit_ms", "outfit_ms")),
    judges = if (is.null(fit$judges)) NULL else
      align(z$judges, z$judges$judge, fit$judges$judge,
            c("fit_resid", "infit_ms", "outfit_ms")))
}

.btl_boot_table <- function(base, reps, stats, side, id, min_success) {
  out <- base
  meta <- list()
  matrices <- list()
  for (st in stats) {
    M <- do.call(rbind, lapply(reps, function(z) z[[st]]))
    mode <- if (st == "chisq") "raw" else if (st == "fit_resid")
      "centred" else "studentised"
    trans <- if (st %in% c("infit_ms", "outfit_ms")) log else identity
    z <- .boot_maxt(base[[st]], M, side[[st]], min_success,
                    mode = mode, transform = trans)
    out[[paste0(st, "_p_boot")]] <- z$p
    out[[paste0(st, "_p_boot_adj")]] <- z$p_adj
    out[[paste0("n_boot_", st)]] <- z$n_boot
    matrices[[st]] <- M
    meta[[st]] <- list(family_n = z$family_n,
                       adjusted_n = z$adjusted_n,
                       joint_replicates = z$family_boot)
  }
  list(table = out, adjustment = meta, replicates = matrices)
}

.btl_fit_bootstrap <- function(fit, B, workers, seed) {
  if (inherits(fit, "rasch_btl_efrm"))
    .refuse("the paired-comparison fit bootstrap does not yet generate the ",
            "frame-dependent unit structure; use it with the single-frame ",
            "paired-comparison fit")
  if (!isTRUE(fit$converged))
    .refuse("the observed paired-comparison fit did not converge; refit it ",
            "successfully before bootstrapping fit")
  B <- .check_whole(B, "B", 1)
  workers <- .check_whole(workers, "workers", 1)
  if (!is.null(seed)) seed <- .check_whole(seed, "seed", 0)
  workers <- min(workers, .rasch_available_workers())
  .sim_seed_check()
  if (!is.null(seed)) {
    old <- .sim_seed_capture(); on.exit(.sim_seed_restore(old), add = TRUE)
    set.seed(seed)
  }
  # Check the design before starting the workers and before spending B refits.
  # The check generates one candidate dataset; restore the stream so that
  # validation cannot change the requested bootstrap draws.
  validation_stream <- .sim_seed_capture()
  tryCatch(invisible(.btl_boot_data(fit)),
           finally = .sim_seed_restore(validation_stream))
  seeds <- sample.int(.Machine$integer.max, B)
  one <- function(b) {
    old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old_stream), add = TRUE)
    set.seed(seeds[b])
    .btl_boot_refit(fit)
  }
  raw <- .rasch_boot_apply(B, one, workers = workers,
                           label = "paired-comparison fit bootstrap")
  status <- vapply(raw, .fit_boot_status, "")
  keep <- status == "ok"
  B_used <- sum(keep); min_success <- .fit_min_boot_success(B)
  if (B_used < min_success)
    .fit_boot_refuse(
      "only ", B_used, " of ", B, " paired-comparison bootstrap replicates ",
      "were usable (", sum(status == "nonconverged"),
      " did not converge; ", sum(status == "error"),
      " otherwise failed); at least ", min_success,
      " are required for the bootstrap null",
      B = B, B_used = B_used,
      B_nonconverged = sum(status == "nonconverged"),
      B_errors = sum(status == "error"))
  reps <- raw[keep]
  if (B_used < .9 * B)
    warning(B - B_used, " of ", B, " paired-comparison bootstrap replicates ",
            "were unusable (", sum(status == "nonconverged"),
            " did not converge; ", sum(status == "error"),
            " otherwise failed); the null uses the remaining ", B_used,
            call. = FALSE)

  pair_reps <- lapply(reps, `[[`, "pairs")
  object_reps <- lapply(reps, `[[`, "objects")
  judge_reps <- if (is.null(fit$judges)) NULL else lapply(reps, `[[`, "judges")
  pairs <- .btl_boot_table(fit$pairs, pair_reps, "chisq",
                           c(chisq = "upper"), "pair", min_success)
  objects <- .btl_boot_table(
    fit$objects, object_reps, c("fit_resid", "infit_ms", "outfit_ms"),
    c(fit_resid = "two", infit_ms = "two", outfit_ms = "two"),
    "object", min_success)
  judges <- if (is.null(fit$judges)) NULL else .btl_boot_table(
    fit$judges, judge_reps, c("fit_resid", "infit_ms", "outfit_ms"),
    c(fit_resid = "two", infit_ms = "two", outfit_ms = "two"),
    "judge", min_success)
  tot <- vapply(reps, `[[`, 0, "total_chisq")
  total <- list(chisq = fit$total_chisq, df = fit$total_df, p = fit$total_p,
                chisq_p_boot = if (sum(is.finite(tot)) >= min_success)
                  .boot_p(fit$total_chisq, tot, "upper") else NA_real_,
                n_boot = sum(is.finite(tot)))
  adjustment <- list(
    method = paste("single-step maximum-statistic bootstrap within each",
                   "statistic and parameter family"),
    pairs = pairs$adjustment, objects = objects$adjustment,
    judges = if (is.null(judges)) NULL else judges$adjustment)
  flat_meta <- c(pairs$adjustment, objects$adjustment,
                 if (is.null(judges)) list() else judges$adjustment)
  incomplete <- vapply(flat_meta, function(z)
    z$family_n > z$adjusted_n ||
      (z$family_n > 0L && z$joint_replicates < min_success), logical(1))
  if (any(incomplete))
    warning("at least one paired-comparison fit family had too few joint ",
            "replicates or no replicated variation for every requested ",
            "adjusted probability; affected values are NA", call. = FALSE)
  .new_fit_bootstrap(list(
    model = "paired comparisons", pairs = pairs$table,
    objects = objects$table,
    judges = if (is.null(judges)) NULL else judges$table,
    total = total, adjustment = adjustment,
    replicates = list(total_chisq = tot, pairs = pairs$replicates,
                      objects = objects$replicates,
                      judges = if (is.null(judges)) NULL else judges$replicates),
    B = B, B_used = B_used, B_failed = B - B_used,
    B_nonconverged = sum(status == "nonconverged"),
    B_errors = sum(status == "error"),
    minimum_usable = min_success, theta = NULL), fit, "btl")
}

# One replicate conditional on the observed scores: for each person, a
# response pattern drawn from the Rasch conditional distribution given their
# raw score over their own observed items. Sufficiency cancels the person
# parameter, so no ability has to be drawn at all -- which is what makes
# this the default. Generating at sampled abilities either inflates the
# person variance (resampled estimates carry their error variance, and the
# standardised fit statistics feel it as the sample grows) or breaks the tie
# between who a person is and what they answered (a booklet design's
# missingness is informative, and an independently drawn ability severs it).
# Conditioning keeps each person's score and missingness glued together and
# reproduces the observed score margins exactly.
#
# Sampling is sequential with suffix elementary-symmetric functions: for
# items j = 1..k and remaining score r, P(x_j = x) is proportional to
# w_j(x) G_{j+1}(r - x), where w_j(x) = exp(-cumulative tau_j(x)) and
# G_{j+1} is the ESF of the items after j. Persons sharing a missingness
# pattern share the suffix tables, so booklets cost one recursion each.
.fit_gen_conditional <- function(X, tau_list, na_mask) {
  N <- nrow(X); L <- length(tau_list)
  obs_mask <- if (is.null(na_mask)) matrix(TRUE, N, L) else !na_mask
  # Work throughout on the log scale. Merely scaling exp(-cumsum(tau)) after
  # exponentiation still produces Inf/Inf for a finite but extreme threshold.
  logW <- lapply(tau_list, function(tau) {
    z <- -c(0, cumsum(tau))
    z - max(z)
  })
  log_add <- function(a, b) {
    use <- is.finite(a) | is.finite(b)
    out <- rep(-Inf, length(a))
    hi <- pmax(a[use], b[use])
    out[use] <- hi + log1p(exp(-abs(a[use] - b[use])))
    out
  }
  out <- matrix(NA_integer_, N, L, dimnames = dimnames(X))
  key <- apply(obs_mask, 1, function(z) paste(which(z), collapse = ","))
  for (kk in unique(key)) {
    rows <- which(key == kk)
    items <- if (nzchar(kk))
      as.integer(strsplit(kk, ",", fixed = TRUE)[[1]]) else integer(0)
    k <- length(items)
    if (!k) next
    m <- vapply(items, function(j) length(logW[[j]]) - 1L, 0L)
    # suffix log-ESFs: G[[t]][s + 1] carries items t..k at score s. A scalar
    # shift at each recursion cancels from the subsequent conditional draw.
    G <- vector("list", k + 1L)
    G[[k + 1L]] <- 0
    for (t in k:1) {
      prev <- G[[t + 1L]]
      cur <- rep(-Inf, m[t] + length(prev))
      w <- logW[[items[t]]]
      for (x in 0:m[t]) {
        at <- (x + 1):(x + length(prev))
        cur[at] <- log_add(cur[at], w[x + 1L] + prev)
      }
      finite <- is.finite(cur)
      if (any(finite)) cur[finite] <- cur[finite] - max(cur[finite])
      G[[t]] <- cur
    }
    for (i in rows) {
      r <- sum(X[i, items])
      for (t in seq_len(k)) {
        w <- logW[[items[t]]]
        xs <- 0:m[t]
        rest <- G[[t + 1L]]
        ok <- (r - xs) >= 0 & (r - xs) <= (length(rest) - 1L)
        lp <- rep(-Inf, length(xs))
        lp[ok] <- w[xs[ok] + 1L] + rest[r - xs[ok] + 1L]
        if (!any(is.finite(lp)))
          stop("internal conditional-bootstrap score has no feasible pattern")
        pr <- exp(lp - max(lp))
        x <- if (t == k) r else sample.int(length(xs), 1L, prob = pr) - 1L
        out[i, items[t]] <- x
        r <- r - x
      }
    }
  }
  out
}

# One replicate of the fitted model: responses for the given locations under
# the estimated thresholds, carrying the observed missing-data pattern so
# that items are tested on the same number of responders they were fitted
# on.
.fit_gen <- function(theta, tau_list, na_mask, item_names) {
  X <- vapply(seq_along(tau_list),
              function(i) .sim_item(theta, tau_list[[i]]),
              integer(length(theta)))
  dim(X) <- c(length(theta), length(tau_list))
  colnames(X) <- item_names
  if (!is.null(na_mask)) X[na_mask] <- NA_integer_
  X
}

# Person locations for the non-default, ability-sampling schemes. Each has a
# documented cost, which is why the score-conditional generator above is the
# default. "resample" generates at the spread of the ESTIMATES -- the true
# variance plus their error variance -- and the standardised fit statistics
# feel that inflation as the sample grows (infit z reached 14% at 4,000
# persons); it also draws abilities independently of each person's
# missingness pattern, which severs a linked-booklet design's informative
# missingness. "normal" corrects the variance but leaves the extreme
# categories of easy and hard items unvisited often enough to lose a
# non-random share of replicates (a third, on a four-category rating scale
# at 400 persons). "fixed" keeps each person's own estimate and mask paired,
# at the estimates' inflated spread. All three remain for comparison and for
# designs where conditioning is not wanted.
.fit_theta_source <- function(theta, X, scheme) {
  if (length(theta) != nrow(X))
    .refuse("the fit does not contain one person location per response row; ",
            "the item fit bootstrap is unavailable")
  finite <- theta[is.finite(theta)]
  if (scheme == "conditional") return(NULL)
  if (scheme %in% c("normal", "resample")) {
    if (length(finite) < 2L)
      .refuse("the fit has fewer than two finite person locations to generate ",
              "from; the item fit bootstrap is unavailable")
    return(finite)
  }
  observed <- rowSums(!is.na(X)) > 0L
  if (any(!is.finite(theta) & observed))
    .refuse("`theta = \"fixed\"` requires a finite person location for ",
            "every row with an observed response")
  theta[!is.finite(theta)] <- if (length(finite)) mean(finite) else 0
  theta
}

.fit_theta <- function(scheme, theta, psi, n = length(theta)) {
  resample_theta <- function(x, size)
    x[sample.int(length(x), size = size, replace = TRUE)]
  switch(scheme,
    normal = {
      vt <- psi$var_theta; mse <- psi$mean_error_var
      if (!is.finite(vt) || !is.finite(mse))
        .refuse("the error-corrected person variance is unavailable; ",
                "`theta = \"normal\"` cannot be generated for this fit")
      sd_c <- sqrt(max(vt - mse, 0))
      if (sd_c <= sqrt(.Machine$double.eps)) rep(mean(theta), n)
      else stats::rnorm(n, mean(theta), sd_c)
    },
    resample = resample_theta(theta, n),
    fixed = {
      if (length(theta) != n)
        stop("fixed person locations must have one value per response row")
      theta
    })
}

#' Bootstrap fit statistics
#'
#' Refers fit statistics to replicated datasets fitted in the same way as the
#' observed data. For a person-by-item Rasch model, the default generator
#' conditions on each person's observed raw score and missingness pattern.
#' The person parameter then cancels by sufficiency. Item parameters and
#' person locations are re-estimated in every replicate.
#' The original class-interval rule is repeated, including automatic per-item
#' allocation with missing responses. An automatic count is resolved against
#' the observed fit and held across replicates, so every replicate chi-square
#' is scored on the intervals the observed one used. Tied locations remain in
#' one interval. A result stored by an earlier version under
#' \code{theta = "resample"}, \code{"fixed"} or \code{"normal"} was built
#' from a null that mixed interval counts, so it is refused and must be
#' recomputed; results from the default generator are unchanged.
#' The generator assumes independent response rows. A fit with repeated person
#' IDs is therefore refused because this bootstrap does not reproduce
#' within-person dependence.
#'
#' Item chi-squares use the upper tail. Fit residuals, infit and outfit use
#' equal-tailed probabilities. Holm adjustment is applied separately to each
#' predeclared item-statistic family; an unavailable item remains in the
#' multiplicity count. Under the conditional generator, the same replicates
#' also give person-specific null distributions. Person probabilities are
#' adjusted with a single-step maximum-statistic distribution across persons
#' for each statistic. A maximum-statistic adjustment is withheld for the
#' complete family if any testable member lacks a usable joint null.
#'
#' The adjustments describe the fitted global null. They do not guarantee
#' familywise error control among otherwise fitting items, persons, objects or
#' judges when another member departs from the model. Each fit statistic is a
#' separate family. The marginal probabilities are retained.
#'
#' For a \code{\link{btl}} fit, outcomes are generated on the fitted comparison
#' design and the model is refitted. The result covers the total pairwise
#' chi-square and pair, object and judge fit. Ordered response thresholds,
#' judge allocation, counts, anchors, position effects and explanatory object
#' restrictions are retained. History-dependent effects are generated in
#' sequence. Half-weighted ties and frame-dependent paired-comparison fits are
#' refused because the present generator does not reproduce those models.
#'
#' Bootstrap probabilities use \eqn{(1+r)/(1+B)}. Their one-sided resolution
#' is therefore \eqn{1/(1+B)} and their equal-tailed resolution is
#' \eqn{2/(1+B)}. For \eqn{L} items, Holm-adjusted inference at .05 requires
#' at least \eqn{20L} replicates for the chi-square and \eqn{40L} for the
#' two-sided statistics.
#' Runs of at least 30 replicates also require at least 30 and 90 percent of
#' the requested refits to be usable. Otherwise inference is withheld because
#' failed refits can select a non-random subset of the bootstrap null.
#'
#' The approach follows Wolfe (2013), who bootstrapped critical values for
#' Rasch fit statistics by generating from the estimated person and item
#' parameters and averaging replicate quantiles into cut points. This
#' implementation conditions on the observed scores instead --- generating
#' at estimated locations carries their estimation error into the null,
#' which his single 1,000-person design could not surface --- and returns
#' per-statistic probabilities under a declared multiplicity policy rather
#' than averaged cut points.
#'
#' The need for a replicated null is not a real-data artefact. Wu and Adams
#' (2013) derive the mean squares' null variance as about \eqn{2/N} but
#' conclude the standardised statistics have a sample-size-independent null,
#' reading contrary reports as flawed simulation or as genuine misfit
#' surfacing in large samples. Under data generated from the model ---
#' neither explanation applies --- infit z beyond 1.96 still flags 11.5
#' percent of correctly fitting items at 250 persons and 68.4 percent at
#' 4,000 with eight items: parameter estimation alone moves the null, and
#' reproducing that estimation in every replicate is what calibrates it.
#'
#' @param fit A fitted object from \code{\link{rasch}} or \code{\link{btl}}.
#'   Extended-frame, many-facet and explanatory person-by-item fits are not
#'   supported. Explanatory paired-comparison fits are supported.
#'   Fully anchored scoring fits are not supported by this refitting procedure.
#'   An older anchored fit must retain its original anchor settings.
#' @param B Positive whole number of bootstrap replicates.
#' @param theta Generator for a person-by-item fit. \code{"conditional"}
#'   retains each observed raw score and missingness pattern. \code{"resample"}
#'   resamples estimated person locations, \code{"fixed"} reuses each
#'   location with the same response row and missingness pattern, and
#'   \code{"normal"} draws from a normal distribution with error-corrected
#'   variance. Fixed generation requires a finite location for each row with
#'   an observed response. This argument does not apply to paired comparisons.
#' @param workers Number of parallel bootstrap workers. The default is four,
#'   reduced when fewer physical cores are available or the R process has a
#'   lower system limit. Per-replicate seeds and random-number generator
#'   settings are shared, so results do not depend on the worker count.
#' @param seed Optional non-negative whole-number seed within the integer
#'   range. The caller's random stream is restored on exit. This calculation
#'   does not support Box--Muller, including when \code{seed = NULL};
#'   see \code{\link{rasch_rng}}.
#' @return An object of class \code{rasch_fit_bootstrap}. For a person-by-item
#'   fit, it contains \code{items},
#'   \code{persons}, \code{total}, \code{replicates}, adjustment metadata and
#'   replicate counts, including separate non-convergence and other-failure
#'   counts. For a paired-comparison fit, the corresponding tables are
#'   \code{pairs}, \code{objects}, \code{judges} and \code{total}.
#'   Runs requesting at least 30 replicates are withheld unless at least 30
#'   and 90 percent of the requested refits are usable.
#' @references Andrich, D. and Marais, I. (2019) \emph{A Course in Rasch
#'   Measurement Theory}. Springer.
#'
#' Wolfe, E. W. (2013). A bootstrap approach to evaluating person and item
#'   fit to the Rasch model. \emph{Journal of Applied Measurement}, 14(1),
#'   1--9.
#'
#' Wu, M. and Adams, R. J. (2013). Properties of Rasch residual fit
#'   statistics. \emph{Journal of Applied Measurement}, 14(4), 339--355.
#'
#' Molenaar, I. W. and Hoijtink, H. (1996). Person-fit test statistics for the
#'   Rasch model. \emph{Applied Measurement in Education}, 9(1), 87--106.
#'
#' Westfall, P. H. and Young, S. S. (1993). \emph{Resampling-Based Multiple
#'   Testing}. Wiley.
#' @seealso \code{\link{chisq_detail}}, \code{\link{btl}},
#'   \code{\link{rasch_rng}}.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' # an exploratory run, kept small for speed: raw probabilities are usable,
#' # and the warning says what Holm-adjusted flagging at .05 would need
#' bs <- suppressWarnings(fit_bootstrap(rasch(X), B = 49, seed = 1))
#' bs$items[c("item", "chisq", "chisq_p_boot", "fit_resid", "fit_resid_p_boot")]
#' @export
fit_bootstrap <- function(fit, B = 200,
                          theta = c("conditional", "resample", "fixed", "normal"),
                          workers = 4L, seed = NULL) {
  if (inherits(fit, "rasch_btl")) {
    if (!missing(theta))
      stop("`theta` applies to person-by-item Rasch fits, not paired comparisons")
    return(.btl_fit_bootstrap(fit, B = B, workers = workers, seed = seed))
  }
  if (!inherits(fit, "rasch"))
    stop("`fit` must be a fitted model from rasch()")
  .require_refittable_calibration(fit)
  if (inherits(fit, c("rasch_efrm", "rasch_mfrm", "rasch_explanatory")))
    .refuse("the item fit bootstrap generates from a single-facet Rasch ",
            "model; an extended-frame, many-facet or explanatory fit has a ",
            "generating structure this function does not reproduce")
  if (.has_repeated_residual_units(fit))
    .refuse("the item fit bootstrap assumes independent response rows and ",
            "does not reproduce within-person dependence; fit the occasions ",
            "separately or use the repeated-measures DIF procedures")
  if (!is.null(fit$disc) && length(unique(fit$disc)) > 1L)
    .refuse("the item fit bootstrap generates under equal discriminations; ",
            "this fit carries frame units that differ across items")
  if (!is.null((fit$refit_spec %||% list())$pc_components))
    .refuse("the item fit bootstrap re-estimates each replicate the way the ",
            "fit was estimated; thresholds estimated through principal ",
            "components are not reproduced")
  if (!isTRUE(fit$est$converged))
    .refuse("the observed Rasch fit did not converge; refit it successfully ",
            "before bootstrapping fit")
  B <- .check_whole(B, "B", 1)
  theta <- match.arg(theta)
  workers <- .check_whole(workers, "workers", 1)
  workers <- min(workers, .rasch_available_workers())
  .sim_seed_check()

  if (!is.null(seed)) {
    seed <- .check_whole(seed, "seed", 0)
    old <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old), add = TRUE)
    set.seed(seed)
  }

  X <- fit$X
  L <- ncol(X)
  # the resolution floor: a two-sided bootstrap probability cannot fall
  # below 2/(B + 1), so Holm across the items cannot reach .05 unless
  # B is at least 40 times the item count (20 times for the one-sided
  # chi-square). Warn now rather than let a whole run produce probabilities
  # that could never have flagged anything.
  if (2 * L >= 0.05 * (B + 1))
    warning("B = ", B, " cannot reach Holm-adjusted significance at .05 ",
            "for ", L, " two-sided tests (floor 2L/(B+1) = ",
            signif(2 * L / (B + 1), 2), "); use B >= ", 40 * L,
            " for adjusted two-sided inference", call. = FALSE)
  spec <- fit$refit_spec %||% list()
  obs <- list(chisq = fit$items$chisq, fit_resid = fit$items$fit_resid,
              infit_ms = fit$items$infit_ms, outfit_ms = fit$items$outfit_ms,
              infit_z = fit$items$infit_z, outfit_z = fit$items$outfit_z)
  th_full <- fit$person$theta
  na_mask <- if (anyNA(X)) is.na(X) else NULL
  # A fixed-location replicate must keep row n paired with row n's observed
  # missingness pattern. Non-finite locations are harmless only for entirely
  # unobserved rows, whose generated values are all masked again.
  th <- .fit_theta_source(th_full, X, theta)

  # every random draw a replicate depends on is made here, in the parent:
  # its person locations (for the ability-sampling schemes), and the seed its
  # responses are generated under. A worker therefore reproduces its
  # replicate from what it is handed, and the worker count cannot move the
  # result.
  th_b <- if (theta == "conditional") NULL
          else lapply(seq_len(B), function(b)
            .fit_theta(theta, th, fit$psi, n = nrow(X)))
  seeds <- sample.int(.Machine$integer.max, B)
  tau_list <- fit$tau_list; item_names <- colnames(X)
  model <- fit$model; ng <- .refit_n_groups(fit); anchors <- spec$anchors
  maxit <- spec$maxit %||% 60L; tol <- spec$tol %||% 1e-8
  # An automatic count is resolved against the observed fit rather than
  # re-derived per replicate; an explicit request is already fixed.
  fixed_groups <- if (is.null(ng)) .fit_boot_realised_groups(fit) else NULL
  one <- function(b) {
    old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old_stream), add = TRUE)
    set.seed(seeds[b])
    Xb <- if (is.null(th_b)) .fit_gen_conditional(X, tau_list, na_mask)
          else .fit_gen(th_b[[b]], tau_list, na_mask, item_names)
    .fit_refit(Xb, model, ng, anchors, fit$m, maxit, tol,
               fixed_groups = fixed_groups)
  }
  reps <- .rasch_boot_apply(B, one, workers = workers,
                            label = "item fit bootstrap")

  status <- vapply(reps, .fit_boot_status, "")
  keep <- status == "ok"
  B_used <- sum(keep)
  min_success <- .fit_min_boot_success(B)
  if (B_used < min_success)
    .fit_boot_refuse(
      "only ", B_used, " of ", B, " bootstrap replicates were usable (",
      sum(status == "nonconverged"), " did not converge; ",
      sum(status == "error"), " otherwise failed); at least ", min_success,
      " are required for the bootstrap null. The fitted model generates data ",
      "this estimator cannot fit reliably",
      B = B, B_used = B_used,
      B_nonconverged = sum(status == "nonconverged"),
      B_errors = sum(status == "error"))
  reps <- reps[keep]
  # replicates are not lost at random -- an unvisited extreme category is
  # likelier in the less dispersed ones -- so a materially thinned null is
  # worth saying out loud rather than leaving to be read off B_used
  if (B_used < 0.9 * B)
    warning(B - B_used, " of ", B, " bootstrap replicates were unusable (",
            sum(status == "nonconverged"), " did not converge; ",
            sum(status == "error"), " otherwise failed); the null is formed ",
            "from the remaining ", B_used,
            "; replicates are not lost at random, so read the result with ",
            "that in mind", call. = FALSE)

  M <- lapply(names(.FIT_STATS), function(st)
    do.call(rbind, lapply(reps, `[[`, st)))
  names(M) <- names(.FIT_STATS)

  items <- data.frame(item = fit$items$item, stringsAsFactors = FALSE)
  resolution_lost <- character(0)
  insufficient <- character(0)
  for (st in names(.FIT_STATS)) {
    side <- .FIT_STATS[[st]]
    n_stat <- colSums(is.finite(M[[st]]))
    p <- vapply(seq_len(L), function(i) {
      if (n_stat[i] < min_success) return(NA_real_)
      .boot_p(obs[[st]][i], M[[st]][, i], side)
    }, 0)
    usable <- is.finite(p)
    # Keep the predeclared family intact when one item's bootstrap null is
    # unavailable. Dropping that item from Holm's family would make the
    # remaining adjusted probabilities less conservative because a failed
    # test happened to be unreportable.
    family_n <- L
    if (any(is.finite(obs[[st]]) & n_stat < min_success))
      insufficient <- c(insufficient, st)
    adj <- rep(NA_real_, L)
    adj[usable] <- stats::p.adjust(p[usable], method = "holm", n = family_n)
    items[[st]] <- obs[[st]]
    items[[paste0(st, "_p_boot")]] <- p
    items[[paste0(st, "_p_boot_adj")]] <- adj
    items[[paste0("n_boot_", st)]] <- n_stat

    # The pre-run warning uses requested B.  Recheck the resolution from the
    # usable null for each statistic, since failures or statistic-specific NAs
    # can make an initially attainable Holm threshold unattainable.
    mult <- if (side == "upper") 1 else 2
    requested_possible <- mult * family_n / (B + 1) < 0.05
    actual_possible <- any(usable) &&
      family_n * min(mult / (n_stat[usable] + 1)) < 0.05
    if (requested_possible && !actual_possible)
      resolution_lost <- c(resolution_lost, st)
  }
  # Backward-compatible alias for the original chi-square count.
  items$n_boot <- items$n_boot_chisq
  rownames(items) <- NULL
  if (length(insufficient))
    warning("too few usable replicated statistics for ",
            paste(unique(insufficient), collapse = ", "),
            " in at least one item; the affected bootstrap probabilities are NA",
            call. = FALSE)
  if (length(resolution_lost))
    warning("after unusable replicates were removed, Holm-adjusted significance ",
            "at .05 is no longer resolvable for ",
            paste(unique(resolution_lost), collapse = ", "),
            "; increase `B`", call. = FALSE)

  # Person fit has a different nuisance structure from item fit. Under the
  # score-conditional generator, row n keeps its observed raw score and
  # missingness pattern in every replicate, so its replicated statistics are
  # the relevant null for that same response pattern. A joint maximum-
  # statistic reference retains the dependence among people in every replicate.
  persons <- NULL
  person_adjustment <- NULL
  if (theta == "conditional") {
    keep_cols <- intersect(c("id", "raw", "theta", "se"), names(fit$person))
    persons <- fit$person[, keep_cols, drop = FALSE]
    person_adjustment <- list(method = "single-step maximum-statistic bootstrap",
                              statistics = list())
    person_warnings <- character(0)
    for (st in names(.PERSON_FIT_STATS)) {
      Mp <- do.call(rbind, lapply(reps, function(z) z$persons[[st]]))
      mode <- if (st %in% c("fit_resid", "infit_z", "outfit_z"))
        "centred" else "studentised"
      trans <- if (st %in% c("infit_ms", "outfit_ms")) log else identity
      z <- .boot_maxt(fit$person[[st]], Mp, .PERSON_FIT_STATS[[st]],
                      min_success, mode = mode, transform = trans)
      persons[[st]] <- fit$person[[st]]
      persons[[paste0(st, "_p_boot")]] <- z$p
      persons[[paste0(st, "_p_boot_adj")]] <- z$p_adj
      persons[[paste0("n_boot_", st)]] <- z$n_boot
      person_adjustment$statistics[[st]] <- list(
        family_n = z$family_n, adjusted_n = z$adjusted_n,
        joint_replicates = z$family_boot)
      if (z$family_n > z$adjusted_n ||
          (z$family_n > 0L && z$family_boot < min_success))
        person_warnings <- c(person_warnings, st)
    }
    rownames(persons) <- NULL
    if (length(person_warnings))
      warning("too few joint person-fit bootstrap replicates, or no replicated ",
              "variation, for the maxT adjustment of ",
              paste(unique(person_warnings), collapse = ", "),
              "; marginal probabilities remain available but affected ",
              "adjusted probabilities are NA", call. = FALSE)
  }

  # whole-test readings: the summed chi-square, and the mean and SD of the
  # item fit residuals that the fit summary reports. The SD is the reason
  # this one matters -- the convention reads it against 1, and under a true
  # model it is nowhere near 1 until the sample runs to a few thousand.
  ok <- is.finite(obs$chisq)
  tot_rep <- rep(NA_real_, B_used)
  if (any(ok)) {
    complete_chisq <- rowSums(is.finite(M$chisq[, ok, drop = FALSE])) == sum(ok)
    tot_rep[complete_chisq] <- rowSums(
      M$chisq[complete_chisq, ok, drop = FALSE])
  }
  # Use the same observed item set in every replicate. Allowing na.rm to
  # change the set from row to row compares means and SDs with different
  # composition and therefore does not form one bootstrap null statistic.
  fr_ok <- is.finite(obs$fit_resid)
  fr_mean <- fr_sd <- rep(NA_real_, B_used)
  if (any(fr_ok)) {
    complete_fr <- rowSums(
      is.finite(M$fit_resid[, fr_ok, drop = FALSE])) == sum(fr_ok)
    fr_mean[complete_fr] <- rowMeans(
      M$fit_resid[complete_fr, fr_ok, drop = FALSE])
    if (sum(fr_ok) >= 2L)
      fr_sd[complete_fr] <- apply(
        M$fit_resid[complete_fr, fr_ok, drop = FALSE], 1, stats::sd)
  }
  total_p <- function(obs_value, null, side) {
    if (sum(is.finite(null)) < min_success) return(NA_real_)
    .boot_p(obs_value, null, side)
  }
  total <- list(
    chisq = fit$total_chisq, df = fit$total_df, p = fit$total_chisq_p,
    chisq_p_boot = total_p(fit$total_chisq, tot_rep, "upper"),
    fit_resid_mean = fit$item_fit_summary$mean,
    fit_resid_mean_p_boot = total_p(
      fit$item_fit_summary$mean, fr_mean, "two"),
    fit_resid_sd = fit$item_fit_summary$sd,
    fit_resid_sd_p_boot = total_p(fit$item_fit_summary$sd, fr_sd, "two"),
    n_boot = sum(is.finite(tot_rep)),
    n_boot_fit_resid_mean = sum(is.finite(fr_mean)),
    n_boot_fit_resid_sd = sum(is.finite(fr_sd)))
  total_counts <- c(chisq = total$n_boot,
                    fit_resid_mean = total$n_boot_fit_resid_mean,
                    fit_resid_sd = total$n_boot_fit_resid_sd)
  if (any(total_counts < min_success))
    warning("too few usable replicates for whole-test ",
            paste(names(total_counts)[total_counts < min_success],
                  collapse = ", "),
            "; the affected bootstrap probabilities are NA", call. = FALSE)

  .new_fit_bootstrap(list(
    items = items, total = total, replicates = M,
    persons = persons, person_adjustment = person_adjustment,
    B = B, B_used = B_used, B_failed = B - B_used,
    B_nonconverged = sum(status == "nonconverged"),
    B_errors = sum(status == "error"),
    minimum_usable = min_success, theta = theta), fit, "rasch")
}

#' @export
print.rasch_fit_bootstrap <- function(x, ...) {
  kind <- if (identical(x$model_kind, "btl"))
    "Paired-comparison" else "Person-by-item"
  cat(kind, "fit bootstrap\n", sep = " ")
  cat("Usable replicates: ", x$B_used, " of ", x$B,
      if (isTRUE(x$B_failed > 0L))
        sprintf(" (%d failed)", x$B_failed) else "", "\n", sep = "")
  if (identical(x$model_kind, "btl")) {
    cat("Results: ", nrow(x$pairs), " pairs, ", nrow(x$objects), " objects",
        if (!is.null(x$judges)) sprintf(", %d judges", nrow(x$judges)) else "",
        "\n", sep = "")
  } else {
    cat("Results: ", nrow(x$items), " items",
        if (!is.null(x$persons)) sprintf(", %d persons", nrow(x$persons)) else "",
        "\n", sep = "")
  }
  p <- x$total$chisq_p_boot
  cat("Whole-test bootstrap p: ",
      if (length(p) == 1L && is.finite(p)) .fmt_p(p) else "unavailable",
      "\n", sep = "")
  invisible(x)
}
