# Conditional bootstrap sensitivity analysis for DIF
#
# The observed DIF ANOVA uses residuals at estimated person locations and
# class intervals formed from those estimates.  A single set of residuals
# cannot be bootstrapped after the fact: every null sample must be calibrated,
# scored, grouped and tested again.  This file deliberately uses the package's
# exact score-conditional Rasch generator and the public fitting and DIF paths
# so all of those nuisance operations recur in every replicate.

.dif_boot_binding <- function(dif) {
  list(
    fit_signature = dif$fit_signature,
    design = dif$bootstrap_design,
    term_ids = dif$term_ids,
    summary_term_ids = dif$summary_term_ids,
    terms = dif$terms[, intersect(
      c("item", "object", "term", "df", "df_denom", "gg_epsilon",
        "F_value", "p"),
      names(dif$terms)), drop = FALSE])
}

.dif_boot_signature <- function(dif) {
  .fit_boot_md5(.dif_boot_binding(dif))
}

.validate_btl_dif_result <- function(dif, fit) {
  if (is.null(dif)) return(invisible(NULL))
  if (!inherits(dif, "rasch_btl_dif") || !is.data.frame(dif$summary) ||
      !is.data.frame(dif$terms))
    stop("`dif` must be a current btl_dif() result")
  if (!inherits(fit, "rasch_btl") ||
      inherits(fit, c("rasch_btl_efrm", "rasch_btl_explanatory")))
    stop("a btl_dif() result needs its ordinary paired-comparison fit")
  if (is.null(dif$fit_signature))
    stop("`dif` predates fitted-model provenance; recompute it from this fit")
  if (!.fit_boot_signature_matches(dif$fit_signature, fit))
    stop("`dif` was computed from a different fitted model")
  .validate_primary_dif_tables(dif, "object", "band", ":band")
  invisible(NULL)
}

.validate_boot_dif_result <- function(dif, fit) {
  if (inherits(fit, "rasch_btl")) .validate_btl_dif_result(dif, fit)
  else .validate_dif_result(dif, fit)
}

.validate_dif_bootstrap <- function(bootstrap, fit, dif = NULL) {
  if (is.null(bootstrap)) return(invisible(NULL))
  .require_refittable_calibration(fit)
  fail <- function() stop(
    "`dif_bootstrap` is incomplete or internally inconsistent; recompute it with dif_bootstrap()",
    call. = FALSE)
  required <- c(
    "summary", "terms", "replicates", "algorithm", "adjustment", "family_n", "B",
    "B_used", "B_failed", "B_nonconverged", "B_errors",
    "minimum_usable", "alpha", "null_method", "model_kind", "unit",
    "fit_signature", "dif_signature", "result_signature")
  if (!inherits(bootstrap, "rasch_dif_bootstrap") ||
      !all(required %in% names(bootstrap)) ||
      !is.data.frame(bootstrap$terms) ||
      !is.data.frame(bootstrap$summary) ||
      !is.list(bootstrap$replicates) ||
      !is.character(bootstrap$result_signature) ||
      length(bootstrap$result_signature) != 1L ||
      is.na(bootstrap$result_signature))
    fail()
  if (!identical(bootstrap$algorithm, "reference-minp-1")) fail()
  unsigned <- unclass(bootstrap)
  unsigned$result_signature <- NULL
  if (!.fit_boot_hash_matches(bootstrap$result_signature, unsigned)) fail()

  is_whole <- function(x, lower = 0L)
    length(x) == 1L && is.numeric(x) && is.finite(x) &&
      x == floor(x) && x >= lower && x <= .Machine$integer.max
  counts <- c("B", "B_used", "B_failed", "B_nonconverged", "B_errors",
              "minimum_usable", "family_n")
  lower <- c(B = 1L, B_used = 1L, B_failed = 0L,
             B_nonconverged = 0L, B_errors = 0L,
             minimum_usable = 1L, family_n = 1L)
  if (any(!vapply(counts, function(nm)
    is_whole(bootstrap[[nm]], lower[[nm]]), logical(1)))) fail()
  if (bootstrap$B_used + bootstrap$B_failed != bootstrap$B ||
      bootstrap$B_nonconverged + bootstrap$B_errors != bootstrap$B_failed ||
      bootstrap$minimum_usable != .fit_min_boot_success(bootstrap$B) ||
      bootstrap$B_used < bootstrap$minimum_usable) fail()
  if (length(bootstrap$alpha) != 1L || !is.numeric(bootstrap$alpha) ||
      !is.finite(bootstrap$alpha) || bootstrap$alpha <= 0 ||
      bootstrap$alpha >= 1 ||
      length(bootstrap$adjustment) != 1L ||
      !is.character(bootstrap$adjustment) || is.na(bootstrap$adjustment) ||
      length(bootstrap$null_method) != 1L ||
      !is.character(bootstrap$null_method) || is.na(bootstrap$null_method) ||
      length(bootstrap$model_kind) != 1L ||
      !is.character(bootstrap$model_kind) || is.na(bootstrap$model_kind) ||
      length(bootstrap$unit) != 1L || !is.character(bootstrap$unit) ||
      is.na(bootstrap$unit)) fail()

  is_btl <- inherits(fit, "rasch_btl")
  expected_kind <- if (is_btl) "btl" else if (inherits(fit, "rasch_efrm"))
    "efrm" else if (inherits(fit, "rasch_mfrm")) "mfrm" else if (
      inherits(fit, "rasch_explanatory")) "explanatory" else "rasch"
  expected_unit <- if (is_btl) "object" else "item"
  expected_null <- if (is_btl) "fitted-outcome" else if (
    inherits(fit, "rasch_efrm")) "item-set-score conditional" else
      "person-score conditional"
  if (!identical(bootstrap$model_kind, expected_kind) ||
      !identical(bootstrap$unit, expected_unit) ||
      !identical(bootstrap$null_method, expected_null)) fail()
  expected_adjustment <- paste(
    "single-step minimum reference-p over the complete",
    paste0(expected_unit, "-by-DIF-term family"))
  if (!identical(bootstrap$adjustment, expected_adjustment)) fail()

  reps <- bootstrap$replicates
  if (!all(c("F", "p", "min_p") %in% names(reps)) ||
      !is.matrix(reps$F) || !is.numeric(reps$F) ||
      !is.matrix(reps$p) || !is.numeric(reps$p) ||
      !is.numeric(reps$min_p) ||
      !identical(dim(reps$F), c(as.integer(bootstrap$B_used),
                                as.integer(bootstrap$family_n))) ||
      !identical(dim(reps$p), dim(reps$F)) ||
      length(reps$min_p) != bootstrap$B_used ||
      any(!is.finite(reps$F)) || any(!is.finite(reps$p)) ||
      any(reps$p < 0 | reps$p > 1) || any(!is.finite(reps$min_p)) ||
      any(reps$min_p < 0 | reps$min_p > 1) ||
      !identical(colnames(reps$F), colnames(reps$p)) ||
      !isTRUE(all.equal(reps$min_p, apply(reps$p, 1L, min),
                        tolerance = 64 * .Machine$double.eps,
                        check.attributes = FALSE))) fail()

  if (!.fit_boot_signature_matches(bootstrap$fit_signature, fit))
    stop("`dif_bootstrap` was computed from a different fitted model")
  if (!is.null(dif)) {
    .validate_boot_dif_result(dif, fit)
    if (is.null(bootstrap$dif_signature) ||
        !.fit_boot_hash_matches(bootstrap$dif_signature,
                                .dif_boot_binding(dif)))
      stop("`dif_bootstrap` was computed from a different DIF analysis")

    # The observed tables and the bootstrap family must be the exact objects
    # identified by the DIF signature. A class and a matching signature alone
    # do not make a truncated saved result safe to export or restore.
    if (!all(names(dif$terms) %in% names(bootstrap$terms)) ||
        !all(names(dif$summary) %in% names(bootstrap$summary)) ||
        !isTRUE(all.equal(
          as.data.frame(bootstrap$terms[names(dif$terms)]),
          as.data.frame(dif$terms), check.attributes = FALSE,
          tolerance = 0)) ||
        !isTRUE(all.equal(
          as.data.frame(bootstrap$summary[names(dif$summary)]),
          as.data.frame(dif$summary), check.attributes = FALSE,
          tolerance = 0)) ||
        !isTRUE(all.equal(bootstrap$alpha, dif$alpha, tolerance = 0))) fail()

    family <- if (is_btl) dif$term_ids != "band" else
      !dif$term_ids %in% c("Residuals", "ci")
    id_name <- expected_unit
    obs <- dif$terms[family, , drop = FALSE]
    obs_key <- .factor_keys(data.frame(
      unit = as.character(obs[[id_name]]), term = dif$term_ids[family],
      stringsAsFactors = FALSE))
    if (bootstrap$family_n != length(obs_key) ||
        !identical(colnames(reps$F), obs_key)) fail()

    term_cols <- c("p_boot", "p_boot_adj", "n_boot", "significant_boot")
    summary_cols <- c(
      "p_uniform_boot", "p_uniform_boot_adj", "p_nonuniform_boot",
      "p_nonuniform_boot_adj", "n_boot_uniform", "n_boot_nonuniform",
      "uniform_DIF_boot", "nonuniform_DIF_boot")
    if (!all(term_cols %in% names(bootstrap$terms)) ||
        !all(summary_cols %in% names(bootstrap$summary))) fail()

    p_boot <- vapply(seq_along(obs_key), function(j)
      (1 + sum(reps$p[, j] <= obs$p[j])) /
        (bootstrap$B_used + 1), numeric(1))
    p_min <- vapply(seq_along(obs_key), function(j)
      (1 + sum(reps$min_p <= obs$p[j])) /
        (bootstrap$B_used + 1), numeric(1))
    p_adj <- p_min
    same_num <- function(x, y) is.numeric(x) && length(x) == length(y) &&
      isTRUE(all.equal(as.numeric(x), as.numeric(y),
                       tolerance = 64 * .Machine$double.eps,
                       check.attributes = FALSE))
    if (!same_num(bootstrap$terms$p_boot[family], p_boot) ||
        !same_num(bootstrap$terms$p_boot_adj[family], p_adj) ||
        !same_num(bootstrap$terms$n_boot[family],
                  rep(bootstrap$B_used, length(obs_key))) ||
        !is.logical(bootstrap$terms$significant_boot) ||
        anyNA(bootstrap$terms$significant_boot) ||
        !identical(bootstrap$terms$significant_boot[family],
                   p_adj < bootstrap$alpha) ||
        any(!is.na(bootstrap$terms$p_boot[!family])) ||
        any(!is.na(bootstrap$terms$p_boot_adj[!family])) ||
        !same_num(bootstrap$terms$n_boot[!family], rep(0, sum(!family))) ||
        !identical(bootstrap$terms$significant_boot[!family],
                   rep(FALSE, sum(!family)))) fail()

    sm <- bootstrap$summary
    exp_sm <- data.frame(
      p_uniform_boot = NA_real_, p_uniform_boot_adj = NA_real_,
      p_nonuniform_boot = NA_real_, p_nonuniform_boot_adj = NA_real_,
      n_boot_uniform = 0, n_boot_nonuniform = 0)
    exp_sm <- exp_sm[rep(1L, nrow(sm)), , drop = FALSE]
    for (j in seq_len(nrow(sm))) {
      base <- dif$summary_term_ids[j]
      unit_j <- as.character(sm[[id_name]][j])
      u <- which(as.character(obs[[id_name]]) == unit_j &
                 dif$term_ids[family] == base)
      suffix <- if (is_btl) ":band" else ":ci"
      nu <- which(as.character(obs[[id_name]]) == unit_j &
                  dif$term_ids[family] == paste0(base, suffix))
      if (length(u) == 1L) {
        exp_sm$p_uniform_boot[j] <- p_boot[u]
        exp_sm$p_uniform_boot_adj[j] <- p_adj[u]
        exp_sm$n_boot_uniform[j] <- bootstrap$B_used
      }
      if (length(nu) == 1L) {
        exp_sm$p_nonuniform_boot[j] <- p_boot[nu]
        exp_sm$p_nonuniform_boot_adj[j] <- p_adj[nu]
        exp_sm$n_boot_nonuniform[j] <- bootstrap$B_used
      }
    }
    for (nm in names(exp_sm))
      if (!same_num(sm[[nm]], exp_sm[[nm]])) fail()
    if (!identical(as.logical(sm$uniform_DIF_boot),
                   is.finite(exp_sm$p_uniform_boot_adj) &
                     exp_sm$p_uniform_boot_adj < bootstrap$alpha) ||
        !identical(as.logical(sm$nonuniform_DIF_boot),
                   is.finite(exp_sm$p_nonuniform_boot_adj) &
                     exp_sm$p_nonuniform_boot_adj < bootstrap$alpha)) fail()
  }
  invisible(NULL)
}

.dif_boot_refit_ordinary <- function(X, fit, design) {
  .require_refittable_calibration(fit)
  spec <- fit$refit_spec %||% list()
  bf <- tryCatch(
    rasch(
      X, model = fit$model, id = design$id, factors = design$factors,
      n_groups = .refit_n_groups(fit),
      anchors = spec$anchors, maxit = spec$maxit %||% 60L,
      tol = spec$tol %||% 1e-8),
    error = function(e) .fit_boot_failure("error"))
  if (inherits(bf, "rasch_fit_boot_failure")) return(bf)
  if (!isTRUE(bf$est$converged)) return(.fit_boot_failure("nonconverged"))

  # A sparse polytomous replicate can fail to visit a category.  If fitting
  # then shortens or rescales an item, it is no longer the model that formed
  # the declared DIF family and the complete replicate is unusable.
  if (!identical(colnames(bf$X), colnames(fit$X)) ||
      !identical(as.integer(bf$m), as.integer(fit$m)))
    return(.fit_boot_failure("error"))

  bd <- tryCatch(
    dif_anova(
      bf, factors = design$factors, n_groups = design$n_groups,
      p_adjust = design$p_adjust, alpha = design$alpha,
      effects = design$effects,
      sizes = FALSE, id = design$id, within = design$within,
      pool_facets = design$pool_facets),
    error = function(e) .fit_boot_failure("error"))
  if (inherits(bd, "rasch_fit_boot_failure")) return(bd)
  list(dif = bd)
}

.dif_boot_refit_explanatory <- function(X, fit, design) {
  bf <- tryCatch(.explanatory_refit_modified(fit, X, inherit_mc = FALSE),
                 error = function(e) .fit_boot_failure("error"))
  if (inherits(bf, "rasch_fit_boot_failure")) return(bf)
  if (!isTRUE(bf$est$converged)) return(.fit_boot_failure("nonconverged"))
  if (!identical(colnames(bf$X), colnames(fit$X)) ||
      !identical(as.integer(bf$m), as.integer(fit$m)) ||
      ncol(bf$est$B) != ncol(fit$est$B))
    return(.fit_boot_failure("error"))
  bd <- tryCatch(dif_anova(
    bf, factors = design$factors, n_groups = design$n_groups,
    p_adjust = design$p_adjust, alpha = design$alpha,
    effects = design$effects, sizes = FALSE, id = design$id,
    within = design$within, pool_facets = design$pool_facets),
    error = function(e) .fit_boot_failure("error"))
  if (inherits(bd, "rasch_fit_boot_failure")) return(bd)
  list(dif = bd)
}

.dif_boot_refit_mfrm <- function(X, fit, design) {
  vm <- fit$virtual_map
  if (is.null(vm) || nrow(vm) != ncol(X))
    return(.fit_boot_failure("error"))
  spec <- fit$refit_spec %||% list()
  taken <- unique(c(spec$facets %||% fit$facet_spec,
                    names(fit$factors %||% data.frame())))
  fresh <- function(x) {
    while (x %in% taken) x <- paste0(x, ".")
    taken <<- c(taken, x); x
  }
  pn <- fresh("..person"); itn <- fresh("..item"); scn <- fresh("..score")
  N <- nrow(X); L <- ncol(X)
  d <- data.frame(
    person = rep(as.character(fit$person$id), times = L),
    item = rep(as.character(vm$item), each = N),
    score = as.vector(X), stringsAsFactors = FALSE,
    check.names = FALSE)
  names(d)[seq_len(3L)] <- c(pn, itn, scn)
  facets <- spec$facets %||% fit$facet_spec
  for (f in facets) d[[f]] <- rep(as.character(vm[[f]]), each = N)
  bf <- tryCatch(rasch_mfrm(
    d, person = pn, item = itn, score = scn, facets = facets,
    n_groups = .refit_n_groups(fit),
    interaction = spec$interaction %||% fit$interaction,
    factors = fit$factors,
    maxit = spec$maxit %||% 60L, tol = spec$tol %||% 1e-8),
    error = function(e) .fit_boot_failure("error"))
  if (inherits(bf, "rasch_fit_boot_failure")) return(bf)
  if (!isTRUE(bf$est$converged)) return(.fit_boot_failure("nonconverged"))
  if (!identical(as.character(bf$virtual_map$vkey),
                 as.character(fit$virtual_map$vkey)) ||
      !identical(as.integer(bf$m), as.integer(fit$m)))
    return(.fit_boot_failure("error"))
  bd <- tryCatch(dif_anova(
    bf, factors = design$factors, n_groups = design$n_groups,
    p_adjust = design$p_adjust, alpha = design$alpha,
    effects = design$effects, sizes = FALSE, id = design$id,
    within = design$within, pool_facets = design$pool_facets),
    error = function(e) .fit_boot_failure("error"))
  if (inherits(bd, "rasch_fit_boot_failure")) return(bd)
  list(dif = bd)
}

.dif_efrm_matrix <- function(X, fit) {
  vm <- fit$virtual_map
  items <- names(fit$set_of)
  out <- matrix(NA_integer_, nrow(X), length(items),
                dimnames = list(NULL, items))
  for (it in items) {
    cols <- which(vm$item == it)
    if (!length(cols)) next
    z <- X[, cols, drop = FALSE]
    nr <- rowSums(!is.na(z))
    if (any(nr > 1L))
      stop("an EFRM row contains several frame copies of one item")
    hit <- which(nr == 1L)
    if (length(hit)) out[hit, it] <- z[cbind(hit, max.col(!is.na(z[hit, ,
      drop = FALSE]), ties.method = "first"))]
  }
  out
}

.dif_gen_efrm <- function(fit) {
  X <- fit$X; vm <- fit$virtual_map
  out <- matrix(NA_integer_, nrow(X), ncol(X), dimnames = dimnames(X))
  block <- .factor_keys(vm[, c("set", "group"), drop = FALSE])
  for (bb in unique(block)) {
    cols <- which(block == bb)
    rr <- unique(as.numeric(fit$disc[cols]))
    if (length(rr) != 1L || !is.finite(rr) || rr <= 0)
      stop("an EFRM item-set frame does not have one positive unit")
    tl <- lapply(fit$tau_list[cols], `*`, rr)
    z <- .fit_gen_conditional(X[, cols, drop = FALSE], tl,
                              is.na(X[, cols, drop = FALSE]))
    out[, cols] <- z
  }
  out
}

.dif_boot_refit_efrm <- function(X, fit, design) {
  spec <- fit$refit_spec %||% list()
  base <- tryCatch(.dif_efrm_matrix(X, fit),
                   error = function(e) .fit_boot_failure("error"))
  if (inherits(base, "rasch_fit_boot_failure")) return(base)

  # A by-value frame gets a default name, which must not replace the original
  # role names when DIF is tested. Retain the frame components and ordinary
  # factors independently of the response matrix's item names.
  frame_names <- spec$groups
  if (is.null(frame_names)) {
    frame_names <- if (length(fit$frame_group) > 1L)
      fit$frame_group[-1L] else fit$frame_group[1L]
  }
  frame_names <- as.character(frame_names)
  frame_group <- as.character(fit$frame_group %||% character(0))
  factor_names <- spec$factors
  if (is.null(factor_names))
    factor_names <- setdiff(names(fit$factors %||% data.frame()), frame_group)
  factor_names <- as.character(factor_names)
  role_names <- unique(c(frame_names, factor_names))
  fac <- fit$factors
  primary_frame <- frame_group[1L]
  if (is.null(fac) || !is.data.frame(fac) ||
      !length(frame_names) || anyNA(role_names) ||
      !length(primary_frame) || !primary_frame %in% names(fac) ||
      any(!frame_names %in% names(fac)) ||
      any(!factor_names %in% names(fac)) ||
      any(vapply(fac[unique(c(frame_group, factor_names))], length,
                 integer(1)) != nrow(base)))
    return(.fit_boot_failure("error"))
  item_names <- colnames(base)
  if (is.null(item_names))
    return(.fit_boot_failure("error"))
  # Keep the reduced response matrix's item names intact.  Supplying the
  # combined frame by value avoids a collision when an external ordinary
  # factor happens to use an item name; explicit items= still makes that
  # distinction unambiguous.  The original role metadata is restored below.
  factor_data <- if (length(factor_names)) fac[factor_names] else NULL
  bf <- tryCatch(rasch_efrm(
    base, item_sets = fit$set_of, groups = fac[[primary_frame]],
    id = fit$person$id,
    factors = factor_data,
    items = item_names,
    n_groups = .refit_n_groups(fit),
    maxit = spec$maxit %||% 50L, tol = spec$tol %||% 1e-7,
    min_link_persons = spec$min_link_persons %||% 30L,
    se_method = "hybrid", boot_reps = 0L, workers = 1L),
    error = function(e) .fit_boot_failure("error"))
  if (inherits(bf, "rasch_fit_boot_failure")) return(bf)
  if (!isTRUE(bf$est$converged)) return(.fit_boot_failure("nonconverged"))
  bf$frame_group <- fit$frame_group
  bf$factors <- fac[unique(c(frame_group, factor_names)), , drop = FALSE]
  bf$refit_spec$groups <- frame_names
  bf$refit_spec$factors <- factor_names
  if (!identical(as.character(bf$virtual_map$vkey),
                 as.character(fit$virtual_map$vkey)) ||
      !identical(as.integer(bf$m), as.integer(fit$m)))
    return(.fit_boot_failure("error"))
  bd <- tryCatch(dif_anova(
    bf, factors = design$factors, n_groups = design$n_groups,
    p_adjust = design$p_adjust, alpha = design$alpha,
    effects = design$effects, sizes = FALSE, id = design$id,
    within = design$within, pool_facets = design$pool_facets),
    error = function(e) .fit_boot_failure("error"))
  if (inherits(bd, "rasch_fit_boot_failure")) return(bd)
  list(dif = bd)
}

.dif_boot_refit_rasch <- function(X, fit, design) {
  if (inherits(fit, "rasch_efrm"))
    .dif_boot_refit_efrm(X, fit, design)
  else if (inherits(fit, "rasch_mfrm"))
    .dif_boot_refit_mfrm(X, fit, design)
  else if (inherits(fit, "rasch_explanatory"))
    .dif_boot_refit_explanatory(X, fit, design)
  else .dif_boot_refit_ordinary(X, fit, design)
}

.dif_boot_refit_btl <- function(fit, design) {
  d <- tryCatch(.btl_boot_data(fit),
                error = function(e) .fit_boot_failure("error"))
  if (inherits(d, "rasch_fit_boot_failure")) return(d)
  spec <- fit$refit_spec %||% list()
  args <- list(data = d, object_a = "object_a", object_b = "object_b",
               response = "response", count = "count",
               thresholds = spec$thresholds %||% "free",
               position = isTRUE(spec$position), anchors = spec$anchors,
               maxit = spec$maxit %||% 60L, tol = spec$tol %||% 1e-8,
               .object_design = spec$object_design)
  if (any(!is.na(d$judge))) args$judge <- "judge"
  if (isTRUE(spec$has_order)) args$order <- "order"
  nonconv <- FALSE
  bf <- tryCatch(withCallingHandlers(do.call(btl, args), warning = function(w) {
    if (grepl("did NOT converge", conditionMessage(w), fixed = TRUE))
      nonconv <<- TRUE
    invokeRestart("muffleWarning")
  }), error = function(e) .fit_boot_failure("error"))
  if (inherits(bf, "rasch_fit_boot_failure")) return(bf)
  if (nonconv || !isTRUE(bf$converged))
    return(.fit_boot_failure("nonconverged"))
  wanted_effects <- if (is.null(fit$dependence)) character(0) else
    sort(as.character(fit$dependence$effect))
  fitted_effects <- if (is.null(bf$dependence)) character(0) else
    sort(as.character(bf$dependence$effect))
  active <- if (is.null(fit$objects$extreme))
    rep(TRUE, nrow(fit$objects)) else !fit$objects$extreme
  if (!identical(wanted_effects, fitted_effects) ||
      !setequal(as.character(bf$objects$object),
                as.character(fit$objects$object[active])) ||
      any(bf$objects$extreme %||% FALSE))
    return(.fit_boot_failure("error"))
  bd <- tryCatch(do.call(btl_dif, c(list(fit = bf), design)),
                 error = function(e) .fit_boot_failure("error"))
  if (inherits(bd, "rasch_fit_boot_failure")) return(bd)
  list(dif = bd)
}

#' Bootstrap sensitivity analysis for DIF
#'
#' Repeats a residual DIF analysis under the fitted invariant model. Rasch,
#' explanatory Rasch and Multiple Ratings fits condition on each person's
#' raw score. Extended Frames fits condition on the person's subtotal within
#' each observed item set. Comparative Judgement fits draw outcomes from the
#' fitted comparison model. Every replicate refits the model and repeats the
#' complete DIF design. The ordinary adjusted analysis remains the primary
#' DIF result.
#'
#' @details
#' Responses are drawn from the conditional distribution
#' \deqn{P(\mathbf X_p=\mathbf x\mid R_p=r_p,\boldsymbol\delta),}
#' where \eqn{R_p} is person \eqn{p}'s observed raw score over that person's
#' observed items. The person parameter cancels by sufficiency. Thus the
#' generator retains the observed score, booklet or missingness pattern and
#' repeated-person row structure without drawing an ability distribution.
#' For Extended Frames, the same calculation is applied within each item set:
#' the frame unit is common to the set, so conditioning on the set subtotal
#' cancels the person parameter. This conditions on more information than a
#' single total when a person sees several sets, but remains an exact null
#' conditional distribution. For paired comparisons, responses are drawn
#' from the fitted category probabilities (and generated sequentially when
#' history effects were fitted), retaining judges and the comparison design.
#' Half-weighted ties and undefeated or winless objects are refused because
#' they do not supply fitted outcome probabilities for the required null.
#'
#' A replicate contributes only when its calibration converges and every
#' member of the declared item- or object-by-DIF-term family is estimable.
#' Marginal probabilities compare each observed term's F-reference probability
#' with the corresponding replicated probabilities. This puts refits whose
#' numerator or denominator degrees of freedom differ on the same tail scale.
#' Familywise probabilities use the single-step distribution of the smallest
#' term-wise F-reference probability in each replicate. The same transformation
#' is applied to the observed and replicated statistics. These
#' adjustments describe the fitted global invariant null. They do not guarantee
#' familywise error control among otherwise invariant items or objects when
#' another member has DIF, because that departure can affect the fitted
#' calibration and matching scores. For
#' \code{B >= 30}, at least 90 per cent of the requested replicates and no
#' fewer than 30 must be usable; a smaller exploratory run must retain a
#' majority.
#'
#' BTL-EFRM and explanatory Comparative Judgement fits are refused because
#' \code{\link{btl_dif}} does not define judge-group DIF after those structural
#' restrictions. A bootstrap cannot supply an estimand that the fitted model
#' does not define.
#'
#' @param fit A fitted Rasch, Multiple Ratings, Extended Frames,
#'   explanatory Rasch or ordinary Comparative Judgement model.
#'   Fully anchored scoring fits are not supported by this refitting procedure.
#' @param dif A current \code{\link{dif_anova}} or \code{\link{btl_dif}}
#'   result from \code{fit}. For person-by-item models it may be omitted and
#'   the default DIF analysis is then computed. Comparative Judgement requires
#'   an explicit result because its judge factors have no default.
#' @param B Number of bootstrap replicates.
#' @param workers Number of parallel workers. The default is four, subject to
#'   the limits reported by the operating system and job scheduler.
#' @param seed Optional non-negative whole-number seed.
#' @return An object of class \code{"rasch_dif_bootstrap"}. Its
#'   \code{summary} and \code{terms} tables add marginal and familywise
#'   bootstrap probabilities to the observed DIF analysis. \code{replicates}
#'   contains the complete replicated F statistics and their F-reference
#'   probabilities; the requested, usable, non-converged and
#'   failed counts are recorded separately.
#' @seealso \code{\link{dif_anova}}, \code{\link{btl_dif}},
#'   \code{\link{fit_bootstrap}}, \code{\link{rasch_rng}}
#' @references
#' Andrich, D. and Marais, I. (2019). \emph{A Course in Rasch Measurement
#' Theory}. Springer.
#'
#' Westfall, P. H. and Young, S. S. (1993). \emph{Resampling-Based Multiple
#' Testing: Examples and Methods for p-Value Adjustment}. Wiley.
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(300 * 6, 1, .5), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' g <- factor(rep(c("A", "B"), each = 150))
#' fit <- rasch(X, factors = data.frame(group = g))
#' d <- dif_anova(fit)
#' # Small only to keep the example quick; use substantially more replicates
#' # for inferential work.
#' db <- suppressWarnings(dif_bootstrap(fit, d, B = 3, seed = 1))
#' head(db$summary)
#' @export
dif_bootstrap <- function(fit, dif = NULL, B = 999, workers = 4L,
                          seed = NULL) {
  is_btl <- inherits(fit, "rasch_btl")
  if (!inherits(fit, "rasch") && !is_btl)
    stop("`fit` must be a fitted model from rasch(), rasch_mfrm(), ",
         "rasch_efrm(), rasch_explanatory(), or btl()")
  .require_refittable_calibration(fit)
  if (inherits(fit, "rasch_btl_efrm"))
    .refuse("judge-group DIF is not defined after a BTL-EFRM frame ",
            "adjustment, so there is no corresponding DIF null to bootstrap")
  if (inherits(fit, "rasch_btl_explanatory"))
    .refuse("judge-group DIF is not defined for an explanatory Comparative ",
            "Judgement fit, so there is no corresponding DIF null to bootstrap")
  if (is_btl && any(fit$objects$extreme %||% FALSE))
    .refuse("the Comparative Judgement fit contains an undefeated or winless ",
            "object. Its displayed location is an extrapolation rather than ",
            "a fitted null parameter, so DIF cannot be bootstrapped from it")
  converged <- if (is_btl) fit$converged else fit$est$converged
  if (!isTRUE(converged))
    .refuse("the observed calibration did not converge; refit it successfully ",
            "before bootstrapping DIF")
  if (!is_btl && !inherits(fit, "rasch_efrm") &&
      !is.null(fit$disc) && length(unique(fit$disc)) > 1L)
    .refuse("the score-conditional DIF bootstrap requires equal item ",
            "discriminations outside an Extended Frames fit")
  if (!is_btl && !is.null((fit$refit_spec %||% list())$pc_components))
    .refuse("the conditional DIF bootstrap does not reproduce threshold ",
            "estimation through principal components")
  B <- .check_whole(B, "B", 1)
  workers <- .check_whole(workers, "workers", 1)
  if (!is.null(seed)) seed <- .check_whole(seed, "seed", 0)
  workers <- min(workers, .rasch_available_workers())

  if (is.null(dif)) {
    if (is_btl)
      .refuse("Comparative Judgement needs an explicit btl_dif() result ",
              "because its judge factors have no default")
    dif <- dif_anova(fit)
  }
  if (is.null(dif$bootstrap_design) || is.null(dif$term_ids) ||
      is.null(dif$summary_term_ids))
    .refuse("the DIF result predates conditional-bootstrap design ",
            "provenance; recompute it with dif_anova()")
  .validate_boot_dif_result(dif, fit)
  if (length(dif$term_ids) != nrow(dif$terms) ||
      length(dif$summary_term_ids) != nrow(dif$summary))
    stop("the DIF result carries malformed bootstrap design provenance")

  family <- if (is_btl) dif$term_ids != "band" else
    !dif$term_ids %in% c("Residuals", "ci")
  unit_label <- if (is_btl) "object" else "item"
  if (!any(family))
    .refuse("the DIF result contains no ", unit_label,
            "-by-factor tests to bootstrap")
  if (any(!is.finite(dif$terms$F_value[family]) |
          !is.finite(dif$terms$p[family])))
    .refuse("the observed DIF family contains an unavailable term; the ",
            "complete family must be estimable before it can be bootstrapped")
  obs <- dif$terms[family, , drop = FALSE]
  id_name <- if (is_btl) "object" else "item"
  obs_ids <- data.frame(unit = as.character(obs[[id_name]]),
                        term = dif$term_ids[family], stringsAsFactors = FALSE)
  obs_key <- .factor_keys(obs_ids)
  if (anyDuplicated(obs_key))
    stop("the DIF result contains duplicate ", unit_label,
         "-by-term identifiers")

  # A minimum-p probability has resolution 1/(B+1), irrespective of family
  # size.  Ninety-nine draws are still only a one-percentage-point grid, so
  # label smaller runs as exploratory rather than pretending to precision.
  .sim_seed_check()
  K <- nrow(obs)
  if (B < 99L)
    warning("B = ", B, " gives a coarse minimum-p reference for ", K,
            " DIF tests; use at least 99, and preferably more, replicates ",
            "for familywise inference", call. = FALSE)

  # Refuse a design the null generator cannot reproduce before paying for B
  # failed refits.  The check makes one draw but restores the caller's random
  # stream, so it cannot alter the requested bootstrap samples.
  validation_stream <- .sim_seed_capture()
  tryCatch({
    if (is_btl) invisible(.btl_boot_data(fit))
    else if (inherits(fit, "rasch_efrm")) invisible(.dif_gen_efrm(fit))
    else invisible(.fit_gen_conditional(
      fit$X, fit$tau_list, if (anyNA(fit$X)) is.na(fit$X) else NULL))
  }, finally = .sim_seed_restore(validation_stream))

  if (!is.null(seed)) {
    seed <- as.integer(seed)
    old <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old), add = TRUE)
    set.seed(seed)
  }
  seeds <- sample.int(.Machine$integer.max, B)
  X <- if (is_btl) NULL else fit$X
  tau_list <- if (is_btl) NULL else fit$tau_list
  na_mask <- if (!is_btl && anyNA(X)) is.na(X) else NULL
  design <- dif$bootstrap_design

  one <- function(b) {
    old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old_stream), add = TRUE)
    set.seed(seeds[b])
    z <- if (is_btl) {
      .dif_boot_refit_btl(fit, design)
    } else {
      Xb <- if (inherits(fit, "rasch_efrm")) .dif_gen_efrm(fit) else
        .fit_gen_conditional(X, tau_list, na_mask)
      .dif_boot_refit_rasch(Xb, fit, design)
    }
    if (inherits(z, "rasch_fit_boot_failure")) return(z)
    bd <- z$dif
    fam_b <- if (is_btl) bd$term_ids != "band" else
      !bd$term_ids %in% c("Residuals", "ci")
    ids_b <- data.frame(unit = as.character(bd$terms[[id_name]][fam_b]),
                        term = bd$term_ids[fam_b], stringsAsFactors = FALSE)
    key_b <- .factor_keys(ids_b)
    at <- match(obs_key, key_b)
    if (length(key_b) != length(obs_key) || anyNA(at) ||
        anyDuplicated(key_b) ||
        any(!is.finite(bd$terms$F_value[fam_b][at])) ||
        any(!is.finite(bd$terms$p[fam_b][at])))
      return(.fit_boot_failure("error"))
    list(F = bd$terms$F_value[fam_b][at], p = bd$terms$p[fam_b][at])
  }

  null_method <- if (is_btl) "fitted-outcome" else if (
    inherits(fit, "rasch_efrm")) "item-set-score conditional" else
      "person-score conditional"
  reps <- .rasch_boot_apply(B, one, workers = workers,
                            label = paste(null_method, "DIF bootstrap"))
  status <- vapply(reps, .fit_boot_status, "")
  keep <- status == "ok"
  B_used <- sum(keep)
  min_success <- .fit_min_boot_success(B)
  if (B_used < min_success)
    .fit_boot_refuse(
      "only ", B_used, " of ", B,
      " DIF bootstrap replicates retained the complete family (",
      sum(status == "nonconverged"), " did not converge; ",
      sum(status == "error"), " otherwise failed); at least ", min_success,
      " are required", B = B, B_used = B_used,
      B_nonconverged = sum(status == "nonconverged"),
      B_errors = sum(status == "error"))
  if (B_used < 0.9 * B)
    warning(B - B_used, " of ", B,
            " DIF bootstrap replicates were unusable; the null ",
            "uses the remaining complete-family replicates", call. = FALSE)

  reps <- reps[keep]
  Fmat <- do.call(rbind, lapply(reps, `[[`, "F"))
  Pmat <- do.call(rbind, lapply(reps, `[[`, "p"))
  if (!is.matrix(Fmat)) Fmat <- matrix(Fmat, nrow = B_used)
  if (!is.matrix(Pmat)) Pmat <- matrix(Pmat, nrow = B_used)
  colnames(Fmat) <- colnames(Pmat) <- obs_key
  # Sparse class intervals and repeated-measures corrections can change the
  # reference degrees of freedom between refits. Raw F values are then not
  # exchangeable across replicates. Compare the corresponding upper-tail
  # reference probabilities, as the minimum-p family adjustment below does.
  p_boot <- vapply(seq_len(K), function(j)
    (1 + sum(Pmat[, j] <= obs$p[j])) / (B_used + 1), 0)
  # Each term's own F-reference probability puts unlike numerator and
  # denominator degrees of freedom onto a common monotone scale. The
  # bootstrap distribution of their minimum calibrates that transformation;
  # replacing it with ranks computed from these same B rows would collapse
  # toward the 1/B endpoint as the family approaches B and destroy power.
  min_p <- apply(Pmat, 1L, min)
  p_min <- vapply(seq_len(K), function(j)
    (1 + sum(min_p <= obs$p[j])) / (B_used + 1), 0)
  p_boot_adj <- p_min

  terms <- dif$terms
  terms$p_boot <- terms$p_boot_adj <- NA_real_
  terms$n_boot <- 0L
  terms$p_boot[family] <- p_boot
  terms$p_boot_adj[family] <- p_boot_adj
  terms$n_boot[family] <- B_used
  terms$significant_boot <- is.finite(terms$p_boot_adj) &
    terms$p_boot_adj < dif$alpha

  sm <- dif$summary
  sm$p_uniform_boot <- sm$p_uniform_boot_adj <- NA_real_
  sm$p_nonuniform_boot <- sm$p_nonuniform_boot_adj <- NA_real_
  sm$n_boot_uniform <- sm$n_boot_nonuniform <- 0L
  for (j in seq_len(nrow(sm))) {
    base <- dif$summary_term_ids[j]
    unit_j <- as.character(sm[[id_name]][j])
    u <- which(as.character(obs[[id_name]]) == unit_j &
               dif$term_ids[family] == base)
    suffix <- if (is_btl) ":band" else ":ci"
    nu <- which(as.character(obs[[id_name]]) == unit_j &
                dif$term_ids[family] == paste0(base, suffix))
    if (length(u) == 1L) {
      sm$p_uniform_boot[j] <- p_boot[u]
      sm$p_uniform_boot_adj[j] <- p_boot_adj[u]
      sm$n_boot_uniform[j] <- B_used
    }
    if (length(nu) == 1L) {
      sm$p_nonuniform_boot[j] <- p_boot[nu]
      sm$p_nonuniform_boot_adj[j] <- p_boot_adj[nu]
      sm$n_boot_nonuniform[j] <- B_used
    }
  }
  sm$uniform_DIF_boot <- is.finite(sm$p_uniform_boot_adj) &
    sm$p_uniform_boot_adj < dif$alpha
  sm$nonuniform_DIF_boot <- is.finite(sm$p_nonuniform_boot_adj) &
    sm$p_nonuniform_boot_adj < dif$alpha

  out <- list(
    summary = sm, terms = terms,
    replicates = list(F = Fmat, p = Pmat, min_p = min_p),
    algorithm = "reference-minp-1",
    adjustment = paste("single-step minimum reference-p over the complete",
                       paste0(unit_label, "-by-DIF-term family")),
    family_n = K, B = B, B_used = B_used, B_failed = B - B_used,
    B_nonconverged = sum(status == "nonconverged"),
    B_errors = sum(status == "error"), minimum_usable = min_success,
    alpha = dif$alpha, null_method = null_method,
    model_kind = if (is_btl) "btl" else if (inherits(fit, "rasch_efrm"))
      "efrm" else if (inherits(fit, "rasch_mfrm")) "mfrm" else if (
        inherits(fit, "rasch_explanatory")) "explanatory" else "rasch",
    unit = unit_label,
    fit_signature = .fit_boot_signature(fit),
    dif_signature = .dif_boot_signature(dif))
  out <- .tag_tables(out)
  out$result_signature <- .fit_boot_md5(out)
  class(out) <- c("rasch_dif_bootstrap", "list")
  out
}

#' @export
print.rasch_dif_bootstrap <- function(x, ...) {
  cat("DIF bootstrap sensitivity analysis\n")
  cat("Null generator: ", x$null_method %||% "conditional", "\n", sep = "")
  cat("Usable replicates: ", x$B_used, " of ", x$B,
      if (isTRUE(x$B_failed > 0L)) sprintf(" (%d failed)", x$B_failed) else "",
      "\n", sep = "")
  cat("Family: ", x$family_n, " ", x$unit %||% "item",
      "-by-term tests\n", sep = "")
  nflag <- sum(x$terms$significant_boot, na.rm = TRUE)
  cat("Global-null familywise bootstrap flags: ", nflag, "\n", sep = "")
  invisible(x)
}
