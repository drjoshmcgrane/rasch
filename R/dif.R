# rasch :: differential item functioning
# ===========================================================================
# DIF by analysis of variance of the standardised residuals (Hagquist &
# Andrich 2017). For each item, residuals are analysed by person factor and
# trait class interval: a factor main effect indicates uniform DIF and a
# factor-by-interval interaction indicates non-uniform DIF. With several
# person factors they are modelled jointly by dif_anova (main effects by
# default, factor-by-factor interactions optional), with covariance-aware
# logit contrasts on significant group terms and the convention that a
# significant interaction supersedes the main effects of the factors
# involved. Multiplicity across items is handled by Holm familywise
# adjustment.
# ===========================================================================

.dif_factors <- function(fit, factors) {
  n <- nrow(fit$X)
  if (is.null(factors)) factors <- fit$factors
  if (is.null(factors)) stop("no person factors supplied or stored in the fit")
  stored_names <- if (is.null(fit$factors)) character(0) else
    names(fit$factors)
  if (.role_columns(factors, stored_names, n)) {
    # a short character vector names fitted factors; anything else risks
    # recycling a fabricated grouping over the persons
    unknown <- if (is.null(fit$factors)) factors else
      setdiff(factors, names(fit$factors))
    if (length(unknown))
      stop("factor name(s) not stored in the fit: ",
           paste(unknown, collapse = ", "),
           "; name fitted factors or supply one value per person")
    if (anyDuplicated(factors))
      stop("factor(s) named more than once: ",
           paste(unique(factors[duplicated(factors)]), collapse = ", "))
    factors <- fit$factors[, factors, drop = FALSE]
  }
  if (!is.data.frame(factors)) {
    if (!is.atomic(factors) || !is.null(dim(factors)))
      stop("`factors` must be a data frame or an ordinary vector with one ",
           "value per person", call. = FALSE)
    if (length(factors) != n)
      stop("`factors` must give one value per person (", n, "); recycling ",
           "would analyse a fabricated grouping")
    factors <- data.frame(group = factors)
  } else if (nrow(factors) != n)
    stop("`factors` has ", nrow(factors), " rows but the fit has ", n,
         " persons")
  if (is.null(names(factors)) || anyNA(names(factors)) ||
      any(!nzchar(trimws(names(factors)))))
    stop("every DIF factor needs a non-empty name")
  if (anyDuplicated(names(factors)))
    stop("duplicate factor name(s): ",
         paste(unique(names(factors)[duplicated(names(factors))]),
               collapse = ", "))
  # Empty text is missing factor metadata, not a substantive group. Normalise
  # it here so every DIF entry point uses the same respondents and labels.
  # Trimming also prevents visually identical levels such as "A" and " A "
  # from defining different hypotheses.
  factors[] <- lapply(factors, function(v) {
    if (!is.character(v) && !is.factor(v)) return(v)
    z <- .role_text_values(v)
    z[!is.na(z) & !nzchar(z)] <- NA_character_
    if (!is.factor(v)) return(z)
    lev <- unique(.role_text_values(levels(v)))
    lev <- lev[!is.na(lev) & nzchar(lev)]
    factor(z, levels = lev, ordered = is.ordered(v))
  })
  factors
}

.check_dif_factor_levels <- function(factors) {
  n_observed <- vapply(factors, function(v)
    length(unique(v[!is.na(v)])), integer(1))
  if (any(n_observed < 2L))
    stop("DIF factor(s) have fewer than two observed levels: ",
         paste(names(factors)[n_observed < 2L], collapse = ", "),
         call. = FALSE)
  invisible(factors)
}

.check_dif_within <- function(within) {
  if (is.null(within)) return(NULL)
  if (!(is.character(within) || is.factor(within)) ||
      !is.null(dim(within)) || anyNA(within) ||
      any(!nzchar(trimws(as.character(within)))) ||
      anyDuplicated(as.character(within)))
    stop("`within` must be an ordinary vector of unique, non-missing factor names",
         call. = FALSE)
  as.character(within)
}

# A missing person identifier denotes an unknown person, not one shared
# identity. Give every missing entry its own collision-free internal key so
# complete-case filtering, tapply(), and split() retain the response row while
# never treating two unknown identifiers as repeated observations.
.dif_ids <- function(id) {
  if (is.null(id)) return(NULL)
  z <- .role_text_values(id)
  missing <- is.na(z) | !nzchar(z)
  z[missing] <- NA_character_
  known <- unique(z[!missing])
  out <- match(z, known)
  if (any(missing))
    out[missing] <- length(known) + seq_len(sum(missing))
  as.character(out)
}

# Rows that can contribute residual evidence to a DIF analysis. Duplicate
# identifiers on rows with no finite fitted residual must not turn an
# independent design into a repeated-person one or create within-person
# factor variation.
.dif_residual_support <- function(residuals) {
  if (!is.matrix(residuals) || !is.numeric(residuals))
    stop("internal DIF residuals must be a numeric matrix", call. = FALSE)
  rowSums(is.finite(residuals)) > 0L
}

.dif_repeated_support <- function(id, usable) {
  if (is.null(id)) return(FALSE)
  if (length(id) != length(usable))
    stop("internal DIF identifier and residual support lengths differ",
         call. = FALSE)
  # This decides the sampling units for a residual analysis, not for the
  # conditional calibration.  A response row can carry a finite residual
  # without contributing an informative conditional item pair, so the
  # calibration support is not authoritative here.  Restrict the fitted (or
  # supplied) identifiers to the rows that can enter the residual statistic.
  keep <- usable & !is.na(id)
  anyDuplicated(id[keep]) > 0L
}

.dif_varies_within <- function(x, id, usable) {
  keep <- usable & !is.na(id) & !is.na(x)
  if (!any(keep)) return(FALSE)
  by_id <- split(as.character(x[keep]), id[keep])
  any(vapply(by_id, function(v) length(unique(v)) > 1L, logical(1)))
}

# Number of independent respondents contributing to each resolved DIF cell.
# A stacked design can contain several response rows from one person, possibly
# in more than one within-person cell. Count that person once in every cell in
# which they answered the item; counting rows would let replication satisfy a
# minimum-sample rule without adding independent information.
.dif_cell_n <- function(grp, observed, id = NULL) {
  if (length(observed) != length(grp))
    stop("internal DIF support mask has the wrong length", call. = FALSE)
  lev <- if (is.factor(grp)) levels(grp) else
    levels(factor(grp[!is.na(grp)]))
  out <- stats::setNames(integer(length(lev)), lev)
  keep <- !is.na(observed) & observed & !is.na(grp)
  if (!any(keep)) return(out)
  g <- as.character(grp[keep])
  zid <- if (is.null(id)) as.character(which(keep)) else .dif_ids(id)[keep]
  out[] <- vapply(lev, function(lv)
    length(unique(zid[g == lv])), integer(1))
  out
}

# Class intervals for a DIF analysis are set from the cells the analysis
# actually uses: the residual ANOVA crosses trait intervals with group
# levels (or with the factor-combination cells in the factorial), so the
# interval count is chosen to keep the smallest group's expected cell size
# adequate -- independently of the interval count of the overall fit.
.dif_n_groups <- function(fit, grp, cell_min = 30L, id = NULL) {
  ok <- !is.na(grp) & !is.na(fit$person$theta)
  if (!any(ok)) return(2L)
  # with repeated ids the cells are counted in PERSONS, not rows: stacked
  # or duplicated observations must not widen the interval rule
  n_min <- if (!is.null(id)) {
    id <- .dif_ids(id)
    min(tapply(id[ok], droplevels(factor(grp[ok])),
               function(v) length(unique(v))))
  } else min(table(droplevels(factor(grp[ok]))))
  max(2L, min(10L, as.integer(n_min) %/% as.integer(cell_min)))
}

.dif_class_intervals <- function(fit, n_groups) {
  ci <- fit$person$class_interval
  if (is.null(ci) || !identical(n_groups, fit$n_groups))
    ci <- .class_intervals(fit$person$theta, fit$person$extreme, n_groups)
  factor(ci)
}

# Person-level class intervals for a within-subjects analysis: each person
# gets one interval from their mean location, so the interval is a clean
# whole-plot factor. Returned aligned to the rows of the fit.
.dif_person_ci <- function(fit, id, n_groups) {
  id <- .dif_ids(id)
  th <- fit$person$theta; ex <- fit$person$extreme
  pth <- tapply(th, id, mean, na.rm = TRUE)
  pex <- tapply(ex, id, function(v) all(v, na.rm = TRUE))
  pci <- .class_intervals(as.numeric(pth), as.logical(pex), n_groups)
  factor(pci[match(id, names(pth))])
}

# variables of an ANOVA term label, e.g. "g1:ci" -> c("g1", "ci")
.term_vars <- function(term) strsplit(term, ":", fixed = TRUE)[[1]]

# Formula terms are built from safe internal names. For display, quote a
# literal factor name containing the interaction separator so a main effect
# called "age:band" cannot be mistaken for the age-by-band interaction.
.dif_term_label <- function(vars) {
  z <- vapply(vars, function(v) {
    if (grepl("[:`]", v)) paste0("`", gsub("`", "``", v, fixed = TRUE), "`")
    else v
  }, "")
  paste(z, collapse = ":")
}

# ---------------------------------------------------------------------------
# Order-invariant person-level DIF tests. Between-person terms get Type II
# sums of squares on person-level residual means: each term is adjusted for
# every term NOT containing it (so the class interval is always adjusted
# out of every group test), with F against the full model's between-person
# residual -- sequential (Type I) tests let a group factor absorb trait or
# correlated-factor variance in unbalanced designs, flipping which factor
# flags with entry order. Within-person terms are tested on person-by-cell
# means through orthonormal contrasts with the Greenhouse-Geisser epsilon
# correction (Maxwell & Delaney 2004): classical split-plot strata assume
# sphericity, and a nonspherical 4-level null rejected at ~9% nominal 5%.
# ---------------------------------------------------------------------------
.dif_type2 <- function(d, term_labels, resp = "z",
                       variance = c("classical", "hc3", "cr3"),
                       robust_terms = NULL, cluster = NULL,
                       weights = NULL, report_terms = term_labels) {
  variance <- match.arg(variance)
  d$.dif_weights <- if (is.null(weights)) rep(1, nrow(d)) else weights
  mk <- function(tl) stats::as.formula(paste(
    resp, "~", if (length(tl)) paste(tl, collapse = " + ") else "1"))
  fit_model <- function(tl) stats::lm(mk(tl), data = d,
                                      weights = d$.dif_weights)
  full <- tryCatch(fit_model(term_labels),
                   error = function(e) NULL)
  if (is.null(full)) return(NULL)
  rss <- function(m) sum(stats::weighted.residuals(m)^2)
  rss_full <- rss(full)
  df_res <- if (variance == "cr3") length(unique(cluster)) - full$rank else
    stats::df.residual(full)
  if (df_res < 1 || rss_full <= 0) return(NULL)
  mse <- rss_full / df_res
  out <- list()
  for (tt in report_terms) {
    tv <- .term_vars(tt)
    not_cont <- term_labels[!vapply(term_labels, function(u)
      all(tv %in% .term_vars(u)), TRUE)]
    m0 <- fit_model(not_cont)
    m1 <- fit_model(c(not_cont, tt))
    df_t <- stats::df.residual(m0) - stats::df.residual(m1)
    if (df_t < 1) next
    ss_t <- max(rss(m0) - rss(m1), 0)
    Fv <- (ss_t / df_t) / mse
    p_t <- stats::pf(Fv, df_t, df_res, lower.tail = FALSE)
    df_denom <- df_res
    # Judge-level residual means can have very different precision when
    # comparison workloads differ. HC3 retains the equal-judge estimand but
    # does not impose a common residual variance. The robust Wald statistic is
    # reported as an F with the model residual denominator; the separate
    # cell-support guard below avoids presenting this small-sample
    # approximation where a factor level has too few independent judges.
    if (variance == "cr3" ||
        (variance == "hc3" && (is.null(robust_terms) || tt %in% robust_terms))) {
      Fv <- p_t <- NA_real_
      X <- stats::model.matrix(m1)
      asg <- attr(X, "assign")
      labs <- attr(stats::terms(m1), "term.labels")
      ti <- which(vapply(labs, function(lab)
        setequal(.term_vars(tt), .term_vars(lab)), TRUE))
      jj <- which(asg == ti)
      X <- X * sqrt(d$.dif_weights)
      qrX <- qr(X)
      if (length(jj) && qrX$rank == ncol(X)) {
        Xi <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
        if (!is.null(Xi)) {
          if (variance == "cr3") {
            # CR3 is the cluster analogue of HC3. Weighted rows sum to one
            # per person; delete-person leverage accounts for the jointly
            # estimated occasion adjustment as well as the group effects.
            er <- stats::weighted.residuals(m1)
            blocks <- split(seq_len(nrow(X)), cluster)
            scores <- lapply(blocks, function(ii) {
              Xg <- X[ii, , drop = FALSE]
              A <- diag(length(ii)) - Xg %*% Xi %*% t(Xg)
              if (!is.finite(rcond(A)) || rcond(A) < 1e-10) return(NULL)
              drop(crossprod(Xg, solve(A, er[ii])))
            })
            Vr <- if (any(vapply(scores, is.null, TRUE)))
              matrix(NA_real_, ncol(X), ncol(X)) else
                Xi %*% crossprod(do.call(rbind, scores)) %*% Xi
          } else {
            h <- pmin(stats::hatvalues(m1), 1 - 1e-8)
            ae <- stats::weighted.residuals(m1) / pmax(1 - h, 1e-8)
            Vr <- Xi %*% crossprod(X * ae) %*% Xi
          }
          Vt <- Vr[jj, jj, drop = FALSE]
          bt <- stats::coef(m1)[jj]
          Wr <- if (.covariance_is_psd(Vt))
            tryCatch(drop(t(bt) %*% solve(Vt, bt)),
                     error = function(e) NA_real_) else NA_real_
          if (is.finite(Wr) && Wr >= 0) {
            Fv <- Wr / length(jj)
            p_t <- if (is.finite(df_denom) && df_denom > 0)
              stats::pf(Fv, length(jj), df_denom, lower.tail = FALSE) else NA_real_
          }
        }
      }
    }
    out[[length(out) + 1L]] <- data.frame(
      term = tt, df = df_t, df_denom = df_denom, gg_epsilon = NA_real_,
      sum_sq = ss_t, mean_sq = ss_t / df_t,
      F_value = Fv, p = p_t,
      resid_ss = rss_full, stringsAsFactors = FALSE)
  }
  if (!length(out)) return(NULL)
  rbind(do.call(rbind, out),
        data.frame(term = "Residuals", df = df_res, df_denom = NA_real_,
                   gg_epsilon = NA_real_, sum_sq = rss_full,
                   mean_sq = mse, F_value = NA_real_, p = NA_real_,
                   resid_ss = NA_real_, stringsAsFactors = FALSE))
}

# Within-stratum tests on the person-by-within-cell mean matrix Y (complete
# cases over cells). For a term pairing the within subspace (Kronecker
# contrast matrix over the within factors) with a between portion, the
# contrast scores S = Y C are tested by Type II model comparison of each
# score column on the between design, pooled over columns, with df scaled
# by the Greenhouse-Geisser epsilon of the residual score covariance.
.dif_within_tests <- function(Y, pdat, wname, wlv, within_terms,
                              bterms_all) {
  n <- nrow(Y)
  contr_of <- function(k) {
    C <- stats::contr.helmert(k)
    sweep(C, 2, sqrt(colSums(C^2)), "/")
  }
  meanvec_of <- function(k) matrix(1 / sqrt(k), k, 1)
  mk <- function(tl, resp) stats::as.formula(paste(
    resp, "~", if (length(tl)) paste(tl, collapse = " + ") else "1"))
  out <- list(); resid_pool <- 0; resid_df <- 0
  na_row <- function(tt) data.frame(
    term = tt, df = NA_real_, df_denom = NA_real_, gg_epsilon = NA_real_,
    sum_sq = NA_real_, mean_sq = NA_real_, F_value = NA_real_, p = NA_real_,
    resid_ss = NA_real_, stringsAsFactors = FALSE)
  # between design for the scores: only terms whose factors survive the
  # complete-panel filtering with at least two levels (a group observed at
  # a single occasion pattern can lose every complete panel; its
  # interactions are then non-estimable and are reported NA rather than
  # crashing lm with a one-level factor)
  pdat <- droplevels(pdat)
  ok_var <- vapply(names(pdat), function(cn)
    !is.factor(pdat[[cn]]) || nlevels(pdat[[cn]]) >= 2L, TRUE)
  bad_vars <- names(pdat)[!ok_var]
  bt_full <- bterms_all[!vapply(bterms_all, function(u)
    any(.term_vars(u) %in% bad_vars), TRUE)]
  for (tt in within_terms) {
    tv <- .term_vars(tt)
    w_t <- intersect(tv, names(wlv))
    b_t <- setdiff(tv, w_t)
    if (any(b_t %in% bad_vars)) {         # non-estimable after filtering
      out[[length(out) + 1L]] <- na_row(tt)
      next
    }
    Cm <- matrix(1, 1, 1)
    for (wf in names(wlv)) {
      k <- wlv[[wf]]
      Cm <- Cm %x% (if (wf %in% w_t) contr_of(k) else meanvec_of(k))
    }
    S <- Y %*% Cm
    m <- ncol(S)
    fits_j <- lapply(seq_len(m), function(j) {
      dd <- pdat; dd$s_ <- S[, j]
      stats::lm(mk(bt_full, "s_"), data = dd)
    })
    rss_f <- sum(vapply(fits_j, function(f) sum(stats::resid(f)^2), 0))
    dfr1 <- stats::df.residual(fits_j[[1]])
    df_err <- m * dfr1
    if (df_err < 1 || rss_f <= 0) next
    # Greenhouse-Geisser epsilon from the residual score covariance
    E <- vapply(fits_j, stats::resid, numeric(n))
    Sg <- crossprod(as.matrix(E)) / dfr1
    lam <- eigen(Sg, symmetric = TRUE, only.values = TRUE)$values
    lam <- pmax(lam, 0)
    eps <- if (m == 1L || sum(lam^2) <= 0) 1 else
      max(min(sum(lam)^2 / (m * sum(lam^2)), 1), 1 / m)
    if (!length(b_t)) {
      # the within main effect is the grand mean of the contrast scores,
      # adjusted for the between design. Removing the intercept from a
      # FORMULA does nothing when factors are present (R re-parameterises
      # them to absorb the constant), so the design is built explicitly
      # with sum-to-zero factor coding, where the intercept column is the
      # balanced grand mean. Empty factorial cells can still alias this
      # intercept, so its estimability is checked below.
      Xb <- if (length(bt_full)) {
        fml <- stats::as.formula(paste("~", paste(bt_full, collapse = " + ")))
        used <- unique(unlist(lapply(bt_full, .term_vars)))
        fac_cols <- intersect(
          names(pdat)[vapply(pdat, is.factor, TRUE)], used)
        ctr <- stats::setNames(
          rep(list("contr.sum"), length(fac_cols)), fac_cols)
        stats::model.matrix(fml, data = pdat, contrasts.arg = ctr)
      } else matrix(1, n, 1)
      ss_t <- 0
      for (j in seq_len(m)) {
        f_full <- stats::lm.fit(Xb, S[, j])
        f_red <- stats::lm.fit(Xb[, -1, drop = FALSE], S[, j])
        df_t1 <- f_full$rank - f_red$rank
        ss_t <- ss_t + max(sum(f_red$residuals^2) -
                            sum(f_full$residuals^2), 0)
      }
      if (df_t1 < 1L) {
        out[[length(out) + 1L]] <- na_row(tt)
        next
      }
      df_t <- m * df_t1
    } else {
      not_cont <- bt_full[!vapply(bt_full, function(u)
        all(b_t %in% .term_vars(u)), TRUE)]
      rss0 <- rss1 <- 0; df_t1 <- NA_integer_
      for (j in seq_len(m)) {
        dd <- pdat; dd$s_ <- S[, j]
        f0 <- stats::lm(mk(not_cont, "s_"), data = dd)
        f1 <- stats::lm(mk(c(not_cont, paste(b_t, collapse = ":")), "s_"),
                        data = dd)
        rss0 <- rss0 + sum(stats::resid(f0)^2)
        rss1 <- rss1 + sum(stats::resid(f1)^2)
        df_t1 <- stats::df.residual(f0) - stats::df.residual(f1)
      }
      if (is.na(df_t1) || df_t1 < 1) next
      ss_t <- max(rss0 - rss1, 0)
      df_t <- m * df_t1
    }
    Fv <- (ss_t / df_t) / (rss_f / df_err)
    # the p-value is computed at the Greenhouse-Geisser corrected degrees
    # of freedom (eps * df, eps * df_denom); the nominal df, the
    # denominator df, and epsilon are all returned so the test is
    # reproducible and reportable
    out[[length(out) + 1L]] <- data.frame(
      term = tt, df = df_t, df_denom = df_err, gg_epsilon = eps,
      sum_sq = ss_t, mean_sq = ss_t / df_t,
      F_value = Fv,
      p = stats::pf(Fv, eps * df_t, eps * df_err, lower.tail = FALSE),
      resid_ss = rss_f, stringsAsFactors = FALSE)
    resid_pool <- rss_f; resid_df <- df_err
  }
  if (!length(out)) return(NULL)
  rbind(do.call(rbind, out),
        data.frame(term = "Residuals", df = resid_df, df_denom = NA_real_,
                   gg_epsilon = NA_real_, sum_sq = resid_pool,
                   mean_sq = if (resid_df > 0) resid_pool / resid_df else
                     NA_real_, F_value = NA_real_, p = NA_real_,
                   resid_ss = NA_real_, stringsAsFactors = FALSE))
}

# Shared validation of the DIF-family arguments; every public DIF entry
# point applies the same rules.
.check_dif_args <- function(alpha, p_adjust, flag_logits = NULL,
                            min_n = NULL, n_groups = NULL) {
  .check_prob(alpha, "alpha")
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      !is.null(dim(p_adjust)) || !is.null(oldClass(p_adjust)) ||
      !p_adjust %in% stats::p.adjust.methods)
    stop("`p_adjust` must name a method in stats::p.adjust.methods",
         call. = FALSE)
  if (!is.null(flag_logits) &&
      (length(flag_logits) != 1L || !is.numeric(flag_logits) ||
       is.complex(flag_logits) || !is.null(dim(flag_logits)) ||
       !is.null(oldClass(flag_logits)) ||
       !is.finite(flag_logits) || flag_logits <= 0))
    stop("`flag_logits` must be one positive finite practical threshold",
         call. = FALSE)
  if (!is.null(min_n)) .check_whole(min_n, "min_n", 1)
  if (!is.null(n_groups)) .check_whole(n_groups, "n_groups", 2)
  invisible(NULL)
}

# Automatic follow-ups must use the exact normalized factor design from the
# omnibus analysis.  This stamp is deliberately separate from the mixed-panel
# algorithm stamp: a result may have the current omnibus adjustment while its
# stored post-hoc estimates still come from the older name-only hand-off.
.dif_followup_algorithm <- "normalized-design-1"

.dif_followups_current <- function(dif) {
  if (!is.list(dif)) return(FALSE)
  present <- intersect(c("sizes", "posthoc"), names(dif))
  # This compatibility helper may also be used while restoring projects, so
  # malformed follow-up fields must never be mistaken for an empty analysis.
  if (length(present) &&
      any(!vapply(dif[present], is.data.frame, logical(1)))) return(FALSE)
  has_followups <- any(vapply(dif[present], nrow, integer(1)) > 0L)
  !has_followups ||
    identical(dif$followup_algorithm, .dif_followup_algorithm)
}

#' Differential item functioning by residual analysis of variance
#'
#' Tests uniform and non-uniform DIF by analysing each item's standardised
#' residuals over person factors and trait class intervals (Andrich and Marais
#' 2019, ch. 16). Several person factors are fitted jointly. The function also
#' supports designs containing both between-person and within-person factors.
#'
#' @details
#' With one factor \eqn{G} and class interval \eqn{C}, the residual model is
#' \deqn{z=\mu+G+C+G\mathbin{:}C+\varepsilon.}
#' The factor term tests uniform DIF and its interaction with class interval
#' tests non-uniform DIF. With several factors, \code{effects = "main"} fits
#' \code{(f1 + f2 + ...) * ci}; \code{effects = "factorial"} also includes
#' factor-by-factor interactions. Type II sums of squares are used. The
#' multiplicity adjustment covers all item-by-DIF-term tests, including both
#' uniform and non-uniform DIF; the class-interval main effect is a nuisance
#' term and is not included. A reported term remains in this family when its
#' probability is unavailable.
#' Effects that cannot be estimated from the retained design are reported
#' as \code{NA}, including within-person effects whose adjusted mean is
#' confounded with between-person terms in an incomplete factorial design.
#'
#' When identifiers repeat, the person is the unit of analysis. Between-person
#' terms use person means and the between-person error stratum. Within-person
#' terms use orthonormal contrasts of person-by-cell means. A
#' Greenhouse--Geisser correction is applied to within-person factors with
#' more than two levels. Persons missing a required cell are excluded from the
#' corresponding within-person test. Required cells include every combination
#' of the within-person factor levels, even when a combination or level has
#' no observations for an item. Uniform
#' between-person factor terms use HC3 covariance so unequal group sizes,
#' leverage, and differing precision of person means do not impose a common
#' residual variance. Class-interval interactions retain the residual-ANOVA
#' reference used to test non-uniform DIF.
#' In incomplete mixed designs, the between-person tests instead fit the
#' declared occasion and person-factor model jointly to person-by-cell means.
#' Each person has total weight one. All between-person terms then use
#' person-cluster CR3 covariance, including uncertainty in the occasion
#' adjustment, with an approximate F reference whose denominator degrees of
#' freedom are the number of persons minus the full model rank. This branch
#' does not use marginal occasion means to adjust the residuals.
#' For between-person design matrix \eqn{X}, residuals \eqn{e_i}, and leverages
#' \eqn{h_i},
#' \deqn{\widehat{V}_{\mathrm{HC3}}=(X^{\mathsf T}X)^{-1}X^{\mathsf T}
#' \operatorname{diag}\left\{\frac{e_i^2}{(1-h_i)^2}\right\}X
#' (X^{\mathsf T}X)^{-1}.}
#'
#' A significant higher-order factor term supersedes its component terms in
#' the summary. For EFRM fits, frame-defining factors are excluded because
#' they define the model rather than a separate DIF contrast; testing such a
#' factor means stepping outside the model, which is what
#' \code{\link{frame_invariance}} does. MFRM residuals are pooled to
#' underlying items unless \code{pool_facets = FALSE}. EFRM response cells
#' are always pooled by item; the frame-defining factors remain excluded.
#' Inference is available only from a converged calibration.
#'
#' @param fit A fitted object from \code{\link{rasch}},
#'   \code{\link{rasch_mfrm}}, or \code{\link{rasch_efrm}}.
#' @param factors A vector (one factor), a data frame of person factors, or a
#'   character vector naming factor columns nominated in the fit. Defaults to
#'   every factor stored in the fit.
#' @param n_groups Number of trait class intervals. The default uses the
#'   smallest joint factor cell to retain about 30 expected responses per
#'   interval and cell, with between 2 and 10 intervals. The selected value is
#'   returned in \code{n_groups}.
#' @param p_adjust Multiplicity adjustment over all item-by-term tests;
#'   default \code{"holm"}. Use \code{"BH"} only for
#'   false-discovery-rate screening rather than familywise control.
#' @param alpha Significance level applied to the adjusted probabilities.
#' @param effects \code{"main"} (default) models several factors additively
#'   (each factor's main effect and its class-interval interaction, but no
#'   factor-by-factor terms); \code{"factorial"} also crosses the factors
#'   with each other. Immaterial with a single factor.
#' @param id Person identifier for stacked or repeated-measures data. It may
#'   be a column name stored in the fit or a vector with one value per row;
#'   by default the identifier carried by the fit is used.
#' @param within Names of within-person factors. With repeated identifiers,
#'   varying factors are detected automatically when this is omitted. See
#'   Details for the mixed-design analysis.
#' @param pool_facets For MFRM fits: pool residuals to the underlying
#'   items (the default), so DIF is tested per item rather than per
#'   item-by-facet cell; \code{FALSE} tests each cell as its own item.
#'   EFRM response cells are always pooled by item, so this argument does not
#'   alter EFRM fits. Ignored for ordinary fits.
#' @param sizes If \code{TRUE}, refit each flagged item-term and calculate
#'   marginal contrasts in logits using \code{\link{dif_posthoc}}. Their
#'   probabilities are adjusted together over the complete family opened by
#'   all flagged, non-superseded uniform terms.
#' @return A list with:
#' \describe{
#'   \item{\code{summary}}{One row per item and group term, containing the
#'   uniform and non-uniform tests, partial eta-squared, adjusted
#'   probabilities, DIF flags, and supersession flag.}
#'   \item{\code{terms}}{The complete item-wise analysis-of-variance tables.}
#'   \item{\code{sizes}}{When requested, marginal pairwise differences for
#'   main effects and difference-in-differences magnitudes for interactions,
#'   adjusted over the complete nominated factor design. This is retained as
#'   an alias of \code{posthoc}.}
#'   \item{\code{posthoc}}{When \code{sizes = TRUE}, marginal pairwise
#'   differences for main effects and difference-in-differences magnitudes
#'   for interactions, calculated by \code{\link{dif_posthoc}} and adjusted
#'   together over the opened follow-up family.}
#'   \item{\code{posthoc_family_n}}{When \code{sizes = TRUE}, the number of
#'   planned questions in that family, including unavailable comparisons.}
#'   \item{\code{followup_algorithm}}{When \code{sizes = TRUE}, records that
#'   stored contrasts used the same normalized factor values as the omnibus
#'   analysis.}
#'   \item{\code{between_covariance}}{The covariance reference used for
#'   uniform between-person terms.}
#' }
#' The remaining components record the factors, class intervals, adjustment,
#' significance level, and design settings.
#' @references
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' Scandinavian Journal of Statistics, 6(2), 65--70.
#'
#' Hagquist, C. and Andrich, D. (2017). Recent advances in analysis of
#' differential item functioning in health research using the Rasch model.
#' Health and Quality of Life Outcomes, 15, 181.
#'
#' MacKinnon, J. G. and White, H. (1985). Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' Journal of Econometrics, 29(3), 305--325.
#'
#' Maxwell, S. E. and Delaney, H. D. (2004). Designing Experiments and
#' Analyzing Data: A Model Comparison Perspective (2nd ed.). Lawrence Erlbaum.
#' @seealso \code{\link{dif_size}}, \code{\link{dif_contrasts}}, and
#'   \code{\link{resolve_dif}}; and \code{\link{frame_invariance}} for the
#'   frame-defining factor this function excludes.
#' @examples
#' set.seed(1); n <- 800
#' d <- seq(-1.5, 1.5, length.out = 6)
#' g1 <- rep(c("a", "b"), each = n / 2)
#' g2 <- rep(c("x", "y"), times = n / 2)
#' sh <- matrix(0, n, 6); sh[g1 == "b", 2] <- 0.8
#' X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
#' colnames(X) <- paste0("I", 1:6)
#' fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))
#' dif_anova(fit)$summary
#'
#' \donttest{
#' # Mixed design: group is between persons and occasion is within persons.
#' N <- 320; theta <- rnorm(N); group <- rep(c("A", "B"), each = N / 2)
#' make_wave <- function(occasion_shift) {
#'   shift <- matrix(0, N, 6)
#'   shift[group == "B", 2] <- 0.9
#'   shift[, 5] <- occasion_shift
#'   matrix(rbinom(N * 6, 1,
#'          plogis(outer(theta, d, "-") - shift)), N, 6)
#' }
#' Xm <- rbind(make_wave(0), make_wave(1.0))
#' colnames(Xm) <- paste0("I", 1:6)
#' repeated <- data.frame(Xm, group = rep(group, 2),
#'                        occasion = rep(c("T1", "T2"), each = N))
#' mixed_fit <- rasch(repeated, id = rep(seq_len(N), 2),
#'                    factors = c("group", "occasion"))
#' mixed_dif <- dif_anova(mixed_fit, within = "occasion")
#' subset(mixed_dif$summary, uniform_DIF | nonuniform_DIF)
#' }
#' @export
dif_anova <- function(fit, factors = NULL, n_groups = NULL,
                                p_adjust = "holm", alpha = 0.05,
                                effects = c("main", "factorial"),
                                sizes = FALSE, id = NULL, within = NULL,
                                pool_facets = TRUE) {
  .check_dif_args(alpha, p_adjust, n_groups = n_groups)
  if (!inherits(fit, "rasch"))
    stop("dif_anova needs a rasch fit")
  if (!isTRUE(sizes) && !isFALSE(sizes))
    stop("`sizes` must be TRUE or FALSE")
  if (!isTRUE(pool_facets) && !isFALSE(pool_facets))
    stop("`pool_facets` must be TRUE or FALSE")
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; DIF inference is unavailable")
  within <- .check_dif_within(within)
  effects <- match.arg(effects)
  Z <- fit$residuals; L <- ncol(Z)
  # Structural residuals pool to UNDERLYING items: users ask whether item A
  # shows DIF, not whether an internal item-by-frame or item-by-facet cell
  # does. In an EFRM a person contributes to at most one frame cell for an
  # underlying item, so the standardised sum is exactly that observed cell.
  # MFRM users may still request per-cell tests with pool_facets = FALSE.
  pooled_note <- NULL
  pooled_structural <- (inherits(fit, "rasch_mfrm") && isTRUE(pool_facets)) ||
    inherits(fit, "rasch_efrm")
  if (pooled_structural &&
      !is.null(fit$virtual_map)) {
    vm <- fit$virtual_map
    items_u <- unique(vm$item)
    Zp <- vapply(items_u, function(it) {
      zz <- Z[, vm$vkey[vm$item == it], drop = FALSE]
      nn <- rowSums(is.finite(zz))
      out <- rowSums(zz, na.rm = TRUE) / sqrt(pmax(nn, 1L))
      out[nn == 0L] <- NA_real_
      out
    }, numeric(nrow(Z)))
    Zp[!is.finite(Zp)] <- NA_real_
    colnames(Zp) <- items_u
    Z <- Zp; L <- ncol(Z)
    pooled_note <- if (inherits(fit, "rasch_efrm")) paste(
      "EFRM response-cell residuals pooled to the underlying items; each",
      "person contributes its observed frame cell and frame-defining factors",
      "remain excluded") else paste(
      "MFRM residuals pooled to the underlying items (standardised sum over",
      "each item's observed facet cells, so rows with different facet",
      "coverage are normalised by the square root of that count);",
      "pool_facets = FALSE tests",
      "each item-by-facet cell as its own item")
  }
  usable_id_rows <- .dif_residual_support(Z)
  factors <- .dif_factors(fit, factors)
  # the EFRM frame group IS the frame structure: each frame has its own
  # virtual items, so the group factor has a single level among any
  # virtual item's responders and cannot be tested as DIF
  drop_frame_note <- NULL
  if (!is.null(fit$frame_group) && any(names(factors) %in% fit$frame_group)) {
    hit <- intersect(names(factors), fit$frame_group)
    if (length(hit) == length(factors))
      .refuse("'", paste(hit, collapse = "', '"),
           "' define(s) the EFRM frame structure itself: the factor is ",
           "constant among the persons responding within any frame, so ",
           "no within-frame comparison remains to test as DIF; nominate ",
           "other person factors")
    factors <- factors[!names(factors) %in% fit$frame_group]
    drop_frame_note <- paste0("frame factor(s) '",
                              paste(hit, collapse = "', '"),
                              "' excluded: they are the frame structure, ",
                              "not testable DIF factors")
  }
  # Check the factors that will actually be analysed. A one-level EFRM frame
  # factor is structural and is removed above; it must not prevent another,
  # estimable person factor from being tested.
  .check_dif_factor_levels(factors)

  # within-subject factors (levels repeating within a person) turn the
  # analysis into a mixed (split-plot) one: the class interval is taken at
  # the person level so it is a clean whole-plot factor, and the within
  # factors carry a person error stratum.
  fitted_id <- is.null(id)
  if (!is.null(id) && (!is.atomic(id) || !is.null(dim(id))))
    stop("`id` must be an ordinary vector or one fitted factor name",
         call. = FALSE)
  if (is.character(id) && length(id) == 1L) {
    if (is.null(fit$factors) || !id %in% names(fit$factors))
      stop("id column '", id, "' not found among the fit's factors")
    id <- fit$factors[[id]]
  } else if (!is.null(id) && length(id) != nrow(factors)) {
    # a recycled identifier understates the sampling units: it would change
    # the person error stratum, and with it every test in the table
    stop("`id` has ", length(id), " entries but the fit has ", nrow(factors),
         " rows; the repeated-measures structure needs one identifier per row")
  }
  if (is.null(id) && !is.null(fit$person$id)) id <- as.character(fit$person$id)
  if (!is.null(id)) id <- .dif_ids(id)
  # a missing identifier is unknown, not shared: counting NAs as repeats
  # would declare a repeated-measures design and change every test
  repeated <- .dif_repeated_support(id, usable_id_rows)
  if (!is.null(within)) {
    unknown <- setdiff(within, names(factors))
    if (length(unknown))
      stop("within-subject factor(s) not among the nominated factors: ",
           paste(unknown, collapse = ", "))
    if (!repeated)
      stop("within-subject factors need repeated person ids (each id ",
           "observed more than once); no id repeats here")
    varies <- vapply(within, function(fn)
      .dif_varies_within(factors[[fn]], id, usable_id_rows), TRUE)
    if (any(!varies))
      stop("factor(s) declared within-subject never vary within any id: ",
           paste(within[!varies], collapse = ", "))
    # the converse is equally ill-defined: a factor that varies within
    # persons has no person-level value, so it cannot be treated as
    # between-subjects (the old row-level treatment pseudo-replicated)
    other <- setdiff(names(factors), within)
    ovaries <- vapply(other, function(fn)
      .dif_varies_within(factors[[fn]], id, usable_id_rows), TRUE)
    if (length(other) && any(ovaries))
      stop("factor(s) vary within persons but are not declared in ",
           "`within`: ", paste(other[ovaries], collapse = ", "),
           "; declare them within-subject (or aggregate the data to one ",
           "row per person) -- treating repeated observations as ",
           "independent between-person rows manufactures information")
  }
  if (is.null(within) && repeated) {
    within <- names(factors)[vapply(names(factors), function(fn)
      .dif_varies_within(factors[[fn]], id, usable_id_rows), TRUE)]
  }
  if (is.null(within)) within <- character(0)
  within <- intersect(within, names(factors))
  mixed <- length(within) > 0L

  if (is.null(n_groups)) {
    cells <- .factor_cells(factors, sep = ".")
    n_groups <- .dif_n_groups(fit, cells,
                              id = if (repeated) id else NULL)
  }
  # PERSONS are the units whenever ids repeat (stacked or duplicated rows):
  # observation-level tests would let copied observations manufacture
  # information. The class interval is taken at the person level so it is a
  # clean whole-plot covariate.
  ci <- if (repeated) .dif_person_ci(fit, id, n_groups) else
    .dif_class_intervals(fit, n_groups)

  fnames <- names(factors)
  safe <- paste0("f", seq_along(fnames))           # syntactic stand-ins
  wsafe <- safe[match(within, fnames)]
  bsafe <- setdiff(safe, wsafe)
  op <- if (effects == "factorial") " * " else " + "
  form_all <- stats::as.formula(
    paste("z ~ (", paste(safe, collapse = op), ") * ci"))
  all_terms <- attr(stats::terms(form_all), "term.labels")
  planned_dif_terms <- setdiff(all_terms, "ci")
  bterms <- all_terms[!vapply(all_terms, function(tt)
    any(.term_vars(tt) %in% wsafe), TRUE)]
  wterms <- setdiff(all_terms, bterms)

  rows <- list()
  incomplete_note <- 0L
  joint_between_items <- character(0)
  for (i in seq_len(L)) {
    d <- data.frame(z = Z[, i], ci = ci)
    d$pid <- if (is.null(id)) sprintf("p%06d", seq_len(nrow(d))) else id
    for (j in seq_along(fnames)) d[[safe[j]]] <- factor(factors[[fnames[j]]])
    d <- d[stats::complete.cases(d), ]
    if (nrow(d) < 10 || any(vapply(safe, function(s)
      length(unique(d[[s]])) < 2, TRUE))) next

    # aggregate to one mean residual per person per within-cell (persons
    # without repeats aggregate to themselves). Within cells are ordered
    # with the LAST within factor varying fastest, matching the Kronecker
    # construction of the contrast matrices. Reversing the structural factor
    # order gives that mixed-radix order without joining factor labels, which
    # could collide when a level contains the separator.
    if (mixed) {
      wcell <- .factor_cells(d[rev(wsafe)], sep = "\r")
    } else wcell <- factor(rep("all", nrow(d)))
    key <- .factor_cells(data.frame(pid = d$pid, wcell = wcell), sep = "\r")
    agz <- tapply(d$z, key, mean)
    firsts <- which(!duplicated(key))
    ag <- d[firsts[match(levels(key), as.character(key[firsts]))],
            c("pid", "ci", safe), drop = FALSE]
    ag$z <- as.numeric(agz)
    ag$wcell <- wcell[firsts[match(levels(key),
                                   as.character(key[firsts]))]]

    # Complete panels retain the person-mean between stratum, centred by
    # within cell and class interval. Every person then contributes the same
    # within-cell composition. This marginal adjustment is not valid for
    # unequal panels; the joint model below replaces it for those items.
    cellci <- .factor_cells(data.frame(wcell = ag$wcell, ci = ag$ci),
                            sep = "\r")
    m_cellci <- tapply(ag$z, cellci, mean)
    m_cell <- tapply(ag$z, ag$wcell, mean)
    ctr <- as.numeric(m_cellci[as.character(cellci)])
    miss_ctr <- is.na(ctr)
    if (any(miss_ctr))
      ctr[miss_ctr] <- as.numeric(m_cell[as.character(ag$wcell)[miss_ctr]])
    zc <- ag$z - ctr
    pkey <- factor(ag$pid)
    pz <- tapply(zc, pkey, mean)
    pfirst <- which(!duplicated(pkey))
    pdat <- ag[pfirst[match(levels(pkey), as.character(pkey[pfirst]))],
               c("pid", "ci", bsafe), drop = FALSE]
    pdat$z <- as.numeric(pz)

    # Person means need not have equal precision in unbalanced or incomplete
    # designs. HC3 retains the equal-person estimand while correcting the
    # between-person covariance for heteroskedasticity and leverage.
    robust_terms <- bterms[!vapply(bterms, function(tt)
      "ci" %in% .term_vars(tt), TRUE)]
    ft_b <- .dif_type2(pdat, bterms, variance = "hc3",
                       robust_terms = robust_terms)
    if (mixed && any(table(ag$pid) < nlevels(ag$wcell))) {
      # Unequal panels must be adjusted jointly. Marginal centring of z
      # alone transfers differences in between-factor composition into the
      # occasion adjustment and can manufacture group DIF. Fit the declared
      # model to person-by-cell means, weighting each person's total as one.
      # CR3 retains person clustering for every between term in this branch.
      np <- table(ag$pid)
      ft_b <- .dif_type2(ag, all_terms, variance = "cr3",
        cluster = ag$pid, weights = 1 / as.numeric(np[ag$pid]),
        report_terms = bterms)
      joint_between_items <- c(joint_between_items, colnames(Z)[i])
    }
    ft_w <- NULL
    if (mixed && length(wterms)) {
      # complete within-cell matrix per person; incomplete persons are
      # dropped explicitly (multi-stratum projections on unbalanced data
      # are exactly the murky territory this engine replaces)
      # d retains the nominated factors' levels after item-wise filtering.
      # Do not reduce an item's question when one level has no responses.
      wl <- lapply(wsafe, function(sn) levels(d[[sn]]))
      names(wl) <- wsafe
      wlv <- lapply(wl, length)
      Ywide <- tapply(ag$z, list(factor(ag$pid), ag$wcell), mean)
      # wcell contains only observed combinations. An entirely absent cell
      # makes every panel incomplete relative to the full crossed design.
      # Count it without allocating potentially enormous all-NA columns.
      full_cells <- ncol(Ywide) == prod(unlist(wlv))
      complete <- full_cells & rowSums(is.na(Ywide)) == 0L
      incomplete_note <- incomplete_note + sum(!complete)
      if (sum(complete) >= 6L) {
        Y <- Ywide[complete, , drop = FALSE]
        pd2 <- pdat[match(rownames(Y), as.character(pdat$pid)), ,
                    drop = FALSE]
        ft_w <- .dif_within_tests(Y, pd2, paste(wsafe, collapse = ":"),
                                  wlv, wterms, bterms)
        # The within-test rows already carry their own residual sum of
        # squares for effect-size calculation. Keeping its second generic
        # "Residuals" row would duplicate the between-stratum key and make
        # an otherwise valid mixed-design result impossible to validate or
        # bootstrap.
        if (!is.null(ft_w))
          ft_w <- ft_w[ft_w$term != "Residuals", , drop = FALSE]
      }
    }
    ft <- rbind(ft_b, ft_w)
    if (is.null(ft)) next
    rows[[length(rows) + 1L]] <- data.frame(item = colnames(Z)[i], ft)
  }
  if (!length(rows)) stop("no item yielded an estimable factorial ANOVA")
  terms <- do.call(rbind, rows)
  rownames(terms) <- NULL

  # The requested item-by-term family is fixed before looking at which
  # item-specific tables are estimable. An item absent from one factor level,
  # or a rank-deficient item-specific design, therefore contributes explicit
  # NA rows instead of silently making the remaining Holm adjustment less
  # stringent.
  wanted <- expand.grid(
    item = colnames(Z), term = planned_dif_terms,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  wanted_key <- .factor_keys(wanted)
  observed_key <- .factor_keys(terms[c("item", "term")])
  missing_family <- !wanted_key %in% observed_key
  n_unavailable_terms <- sum(missing_family)
  if (n_unavailable_terms) {
    absent <- wanted[missing_family, , drop = FALSE]
    absent$df <- absent$df_denom <- absent$gg_epsilon <-
      absent$sum_sq <- absent$mean_sq <- absent$F_value <-
      absent$p <- absent$resid_ss <- NA_real_
    terms <- rbind(terms, absent[names(terms)])
  }

  # partial eta-squared per term against its own stratum residual
  terms$eta2_partial <- ifelse(
    terms$term == "Residuals" | is.na(terms$resid_ss), NA_real_,
    terms$sum_sq / (terms$sum_sq + terms$resid_ss))
  terms$resid_ss <- NULL

  # One reported DIF decision can be triggered by any item-by-term test.
  # Adjust them as one family: separate adjustment of the uniform and
  # non-uniform terms gave an approximately 10% all-null probability of at
  # least one flag at nominal 5% in the balanced MFRM validation design.
  terms$p_adj <- NA_real_
  family_member <- !terms$term %in% c("Residuals", "ci")
  sel_test <- family_member & is.finite(terms$p)
  terms$p_adj[sel_test] <- p.adjust(
    terms$p[sel_test], method = p_adjust, n = sum(family_member))
  terms$significant <- !is.na(terms$p_adj) & terms$p_adj < alpha

  # a significant higher-order GROUP interaction supersedes the lower-order
  # group terms built from a subset of its factors (within the same item).
  # Terms crossing the class interval are excluded from the pass: a term's
  # own ci interaction is reported WITH it as non-uniform DIF, so it must
  # not supersede it (an item significant on both would otherwise drop out
  # of the follow-ups entirely).
  terms$superseded <- FALSE
  is_group_t <- !vapply(terms$term, function(t)
    "ci" %in% .term_vars(t), TRUE)
  for (it in unique(terms$item)) {
    all_terms <- which(terms$item == it & is_group_t)
    higher <- all_terms[terms$significant[all_terms]]
    for (lo in all_terms) for (hi in higher) {
      vl <- .term_vars(terms$term[lo]); vh <- .term_vars(terms$term[hi])
      if (length(vl) < length(vh) && all(vl %in% vh))
        terms$superseded[lo] <- TRUE
    }
  }

  # map a term's syntactic stand-ins (f1..fk) back to the nominated factor
  # names by exact whole-token match, so a factor named like a stand-in
  # ("f1") or like the class interval ("ci") cannot be re-substituted or
  # collide. Applied only for display, after all classification is done on
  # the stand-ins.
  relabel <- function(x) vapply(x, function(t) {
    toks <- strsplit(t, ":", fixed = TRUE)[[1]]
    i <- match(toks, safe)
    toks[!is.na(i)] <- vapply(fnames[i[!is.na(i)]], .dif_term_label, "")
    if ("ci" %in% fnames)
      toks[is.na(i) & toks == "ci"] <- "(class interval)"
    paste(toks, collapse = ":")
  }, character(1), USE.NAMES = FALSE)

  # DIF magnitudes in logits for the significant, non-superseded group
  # terms. The ANOVA adjusts every nominated factor jointly, so the reported
  # magnitude must use the same complete design: equal-cell marginal
  # differences for main effects and tensor contrasts for interactions.
  size_tab <- posthoc_tab <- NULL
  posthoc_family_n <- 0
  if (isTRUE(sizes)) {
    ph <- list(); posthoc_fail <- character(0)
    cand <- terms[terms$significant & !terms$superseded &
                  !vapply(terms$term, function(tt)
                    "ci" %in% .term_vars(tt), TRUE), , drop = FALSE]
    for (r in seq_len(nrow(cand))) {
      it <- cand$item[r]; tt <- cand$term[r]
      by_user <- fnames[match(.term_vars(tt), safe)]
      # Declare every opened question before attempting the refit. A failed
      # or unsupported comparison remains part of the multiplicity family,
      # just as it does in btl_dif().
      planned_n <- prod(vapply(by_user, function(fn) {
        choose(nlevels(droplevels(factor(factors[[fn]]))), 2L)
      }, numeric(1)))
      posthoc_family_n <- posthoc_family_n + planned_n
      if (!is.finite(posthoc_family_n) ||
          posthoc_family_n > .Machine$integer.max)
        stop("the planned DIF follow-up family is too large to adjust",
             call. = FALSE)
      dp <- tryCatch(dif_posthoc(
        # `factors` may be an externally supplied, normalized design rather
        # than the columns retained in the fit.  Passing only `fnames` here
        # makes .dif_factors() resolve those names back to fit$factors and
        # silently runs the follow-up on stale group assignments.  Hand the
        # exact normalized data frame used by the omnibus analysis through
        # to the refit so replacement values, external factors, and
        # incomplete/within designs remain aligned.
        fit, it, term = by_user, factors = factors,
        within = fnames[match(wsafe, safe)],
        id = if (fitted_id) NULL else id,
        p_adjust = p_adjust, alpha = alpha),
        error = function(e) e)
      if (inherits(dp, "error"))
        posthoc_fail <- c(posthoc_fail, conditionMessage(dp))
      else
        ph[[length(ph) + 1L]] <- data.frame(
          item = it, term = tt, dp$table, row.names = NULL)
    }
    posthoc_tab <- if (length(ph)) do.call(rbind, ph) else
      data.frame(item = character(), term = character())
    if (nrow(posthoc_tab)) {
      posthoc_tab$p_adj <- .p_adjust_family(
        posthoc_tab$p, method = p_adjust, n = posthoc_family_n)
      posthoc_tab$significant <- !is.na(posthoc_tab$p_adj) &
        posthoc_tab$p_adj < alpha
    }
    size_tab <- posthoc_tab
    size_note <- NULL
    posthoc_note <- if (length(posthoc_fail)) paste0(
      "DIF post-hoc comparisons unavailable for some flagged term(s): ",
      paste(unique(posthoc_fail), collapse = "; ")) else NULL
  } else size_note <- posthoc_note <- NULL

  # compact reading: one row per item and group term, its own effect being
  # uniform DIF and its crossing with the class interval non-uniform DIF
  gterms <- setdiff(unique(terms$term), "Residuals")
  gterms <- gterms[!vapply(gterms, function(tt)
    "ci" %in% .term_vars(tt), TRUE)]
  srows <- list()
  for (it in unique(terms$item)) for (tt in gterms) {
    u <- terms[terms$item == it & terms$term == tt, , drop = FALSE]
    nu <- terms[terms$item == it & terms$term == paste0(tt, ":ci"), ,
                drop = FALSE]
    if (!nrow(u)) next
    # .dif_type2 and .dif_within_tests each emit one row per term; a term
    # reaching here twice would silently break every isTRUE() below, so
    # fail loudly rather than trust that upstream invariant
    if (nrow(u) > 1 || nrow(nu) > 1)
      stop("internal error: term '", tt, "' appears in several error strata")
    srows[[length(srows) + 1L]] <- data.frame(
      item = it, term = tt,
      F_uniform = u$F_value, p_uniform = u$p,
      p_uniform_adj = u$p_adj, eta2_uniform = u$eta2_partial,
      uniform_DIF = isTRUE(u$significant),
      F_nonuniform = if (nrow(nu)) nu$F_value else NA_real_,
      p_nonuniform = if (nrow(nu)) nu$p else NA_real_,
      p_nonuniform_adj = if (nrow(nu)) nu$p_adj else NA_real_,
      eta2_nonuniform = if (nrow(nu)) nu$eta2_partial else NA_real_,
      nonuniform_DIF = nrow(nu) > 0 && isTRUE(nu$significant),
      superseded = isTRUE(u$superseded))
  }
  summary_tab <- do.call(rbind, srows)
  # every term aliased with the class interval leaves nothing to report: a
  # NULL summary would be returned as a malformed object and fail later in
  # the print method and the report
  if (is.null(summary_tab) || !nrow(summary_tab))
    .refuse("no DIF contrast is estimable: the nominated factor(s) are ",
            "aliased with the trait class interval (each interval holds ",
            "one group only), so no within-interval comparison remains")
  rownames(summary_tab) <- NULL
  # Keep the exact factor names aligned with the summary rows. Public labels
  # are deliberately readable strings; downstream refits and the app must
  # not recover structure by splitting those strings on punctuation.
  summary_factors <- lapply(summary_tab$term, function(tt)
    fnames[match(.term_vars(tt), safe)])

  # Retain the collision-free internal term identifiers before replacing
  # them with display labels.  A conditional bootstrap must align the same
  # item-by-term family across refits without trying to parse factor names
  # that may themselves contain punctuation or be called "ci".
  term_ids <- terms$term
  summary_term_ids <- summary_tab$term

  # relabel the stand-ins to the nominated names for display, now that all
  # classification is done
  terms$term <- relabel(terms$term)
  summary_tab$term <- relabel(summary_tab$term)
  if (!is.null(size_tab) && nrow(size_tab))
    size_tab$term <- relabel(size_tab$term)
  if (!is.null(posthoc_tab) && nrow(posthoc_tab))
    posthoc_tab$term <- relabel(posthoc_tab$term)

  notes <- character(0)
  if (!is.null(pooled_note)) notes <- c(notes, pooled_note)
  if (!is.null(drop_frame_note)) notes <- c(notes, drop_frame_note)
  if (!is.null(size_note)) notes <- c(notes, size_note)
  if (!is.null(posthoc_note)) notes <- c(notes, posthoc_note)
  if (incomplete_note > 0L)
    notes <- c(notes, sprintf(
      "%d person-by-item panel(s) missing a within-subject cell were dropped from the within-person tests (their between-person information is retained)",
      incomplete_note))
  if (length(joint_between_items))
    notes <- c(notes, paste(
      "Incomplete panels use joint occasion and person-factor adjustment",
      "with equal total weight per person and person-cluster CR3 covariance",
      "for all between-person terms; F references are approximate. Items:",
      paste(joint_between_items, collapse = ", ")))
  if (n_unavailable_terms > 0L)
    notes <- c(notes, sprintf(
      "%d requested item-term test(s) were not estimable and remain in the adjusted-probability family",
      n_unavailable_terms))
  if (mixed && any(is.na(terms$F_value) & terms$term != "Residuals"))
    notes <- c(notes, paste(
      "term(s) reported NA were not estimable from the retained design",
      "(rank deficiency, insufficient residual information, or incomplete panels)"))
  out <- list(summary = summary_tab, terms = terms,
              summary_factors = summary_factors,
              term_ids = term_ids,
              summary_term_ids = summary_term_ids,
              n_groups = nlevels(as.factor(ci)), within = within,
              algorithm = "joint-between-1",
              factor_names = fnames,
              between_covariance = if (length(joint_between_items))
                "HC3 for uniform factor terms; CR3 for incomplete panels" else
                "HC3 for uniform factor terms",
              effects = effects, alpha = alpha, p_adjust = p_adjust,
              notes = notes,
              fit_signature = .fit_boot_signature(fit),
              bootstrap_design = list(
                factors = factors,
                id = id,
                within = if (length(within)) within else NULL,
                n_groups = nlevels(as.factor(ci)), effects = effects,
                pool_facets = pool_facets, p_adjust = p_adjust,
                alpha = alpha))
  if (isTRUE(sizes)) {
    out$sizes <- size_tab
    out$posthoc <- posthoc_tab
    out$posthoc_family_n <- as.integer(posthoc_family_n)
    out$followup_algorithm <- .dif_followup_algorithm
  }
  out <- .tag_tables(out)
  out$result_signature <- .fit_boot_md5(out)
  class(out) <- "rasch_dif"
  out
}

.validate_primary_dif_tables <- function(dif, unit, excluded, suffix) {
  fail <- function() stop(
    "`dif` is incomplete or internally inconsistent; recompute the DIF analysis",
    call. = FALSE)
  required <- c("summary", "terms", "term_ids", "summary_term_ids",
                "summary_factors", "effects", "alpha", "p_adjust",
                "bootstrap_design", "fit_signature", "result_signature")
  if (!all(required %in% names(dif)) || !is.data.frame(dif$summary) ||
      !is.data.frame(dif$terms) || !is.list(dif$bootstrap_design) ||
      !is.character(dif$result_signature) ||
      length(dif$result_signature) != 1L || is.na(dif$result_signature)) fail()
  unsigned <- unclass(dif)
  unsigned$result_signature <- NULL
  if (!.fit_boot_hash_matches(dif$result_signature, unsigned)) fail()

  if (!is.numeric(dif$alpha) || length(dif$alpha) != 1L ||
      !is.finite(dif$alpha) || dif$alpha <= 0 || dif$alpha >= 1 ||
      !is.character(dif$p_adjust) || length(dif$p_adjust) != 1L ||
      is.na(dif$p_adjust) || !dif$p_adjust %in% stats::p.adjust.methods ||
      !is.character(dif$effects) || length(dif$effects) != 1L ||
      is.na(dif$effects) || !dif$effects %in% c("main", "factorial")) fail()
  design <- dif$bootstrap_design
  if (!all(c("factors", "effects", "p_adjust", "alpha") %in%
           names(design)) ||
      !identical(design$effects, dif$effects) ||
      !identical(design$p_adjust, dif$p_adjust) ||
      !isTRUE(all.equal(design$alpha, dif$alpha, tolerance = 0))) fail()
  if (identical(unit, "item")) {
    if (!is.character(dif$factor_names) || !length(dif$factor_names) ||
        anyNA(dif$factor_names) || any(!nzchar(dif$factor_names)) ||
        anyDuplicated(dif$factor_names) ||
        !is.numeric(dif$n_groups) || length(dif$n_groups) != 1L ||
        !is.finite(dif$n_groups) || dif$n_groups != floor(dif$n_groups) ||
        dif$n_groups < 2L || !identical(design$n_groups, dif$n_groups) ||
        !identical(design$within, if (length(dif$within)) dif$within else NULL))
      fail()
  } else {
    if (!is.character(dif$factors) || !length(dif$factors) ||
        anyNA(dif$factors) || any(!nzchar(dif$factors)) ||
        anyDuplicated(dif$factors) || !is.list(design$factors) ||
        !identical(names(design$factors), dif$factors)) fail()
  }
  term_required <- c(unit, "term", "F_value", "p", "p_adj",
                     "eta2_partial", "significant", "superseded")
  summary_required <- c(
    unit, "term", "F_uniform", "p_uniform", "p_uniform_adj",
    "eta2_uniform", "uniform_DIF", "F_nonuniform", "p_nonuniform",
    "p_nonuniform_adj", "eta2_nonuniform", "nonuniform_DIF",
    "superseded")
  if (!all(term_required %in% names(dif$terms)) ||
      !all(summary_required %in% names(dif$summary)) ||
      !is.character(dif$term_ids) || length(dif$term_ids) != nrow(dif$terms) ||
      anyNA(dif$term_ids) || any(!nzchar(dif$term_ids)) ||
      !is.character(dif$summary_term_ids) ||
      length(dif$summary_term_ids) != nrow(dif$summary) ||
      anyNA(dif$summary_term_ids) || any(!nzchar(dif$summary_term_ids)) ||
      !is.list(dif$summary_factors) ||
      length(dif$summary_factors) != nrow(dif$summary) ||
      any(!vapply(dif$summary_factors, function(x)
        is.character(x) && length(x) && !anyNA(x) && all(nzchar(x)),
        logical(1))) ||
      !is.character(dif$terms[[unit]]) || anyNA(dif$terms[[unit]]) ||
      any(!nzchar(dif$terms[[unit]])) ||
      !is.character(dif$summary[[unit]]) || anyNA(dif$summary[[unit]]) ||
      any(!nzchar(dif$summary[[unit]]))) fail()
  term_key <- .factor_keys(data.frame(
    unit = dif$terms[[unit]], term = dif$term_ids,
    stringsAsFactors = FALSE))
  summary_key <- .factor_keys(data.frame(
    unit = dif$summary[[unit]], term = dif$summary_term_ids,
    stringsAsFactors = FALSE))
  if (anyDuplicated(term_key) || anyDuplicated(summary_key)) fail()

  valid_num <- function(x) is.numeric(x) && !any(is.nan(x))
  valid_p <- function(x) valid_num(x) &&
    all(is.na(x) | (is.finite(x) & x >= 0 & x <= 1))
  same_num <- function(x, y) valid_num(x) && length(x) == length(y) &&
    isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 0,
                     check.attributes = FALSE))
  if (!valid_num(dif$terms$F_value) ||
      any(!is.na(dif$terms$F_value) & dif$terms$F_value < 0) ||
      !valid_p(dif$terms$p) || !valid_p(dif$terms$p_adj) ||
      !valid_p(dif$terms$eta2_partial) ||
      !is.logical(dif$terms$significant) || anyNA(dif$terms$significant) ||
      !is.logical(dif$terms$superseded) || anyNA(dif$terms$superseded)) fail()
  family <- !dif$term_ids %in% excluded
  expected_adj <- rep(NA_real_, nrow(dif$terms))
  usable <- family & is.finite(dif$terms$p)
  expected_adj[usable] <- stats::p.adjust(
    dif$terms$p[usable], method = dif$p_adjust, n = sum(family))
  expected_sig <- !is.na(expected_adj) & expected_adj < dif$alpha
  if (!same_num(dif$terms$p_adj, expected_adj) ||
      !identical(dif$terms$significant, expected_sig)) fail()

  # The compact table is a view of the two corresponding omnibus rows. Do
  # not let reports or restored projects carry a summary whose decisions no
  # longer agree with the term family that was actually adjusted.
  map_num <- list(
    F_uniform = "F_value", p_uniform = "p", p_uniform_adj = "p_adj",
    eta2_uniform = "eta2_partial")
  map_non <- list(
    F_nonuniform = "F_value", p_nonuniform = "p",
    p_nonuniform_adj = "p_adj", eta2_nonuniform = "eta2_partial")
  if (identical(unit, "object")) {
    extra <- c("min_judges", "min_effective_judges",
               "inference_available")
    if (!all(extra %in% names(dif$terms)) ||
        !all(c("min_judges", "min_effective_judges", "uniform_inference",
               "nonuniform_inference") %in% names(dif$summary))) fail()
  }
  for (j in seq_len(nrow(dif$summary))) {
    u <- which(dif$terms[[unit]] == dif$summary[[unit]][j] &
               dif$term_ids == dif$summary_term_ids[j])
    nu <- which(dif$terms[[unit]] == dif$summary[[unit]][j] &
                dif$term_ids == paste0(dif$summary_term_ids[j], suffix))
    if (length(u) != 1L || length(nu) > 1L) fail()
    for (nm in names(map_num))
      if (!same_num(dif$summary[[nm]][j], dif$terms[[map_num[[nm]]]][u]))
        fail()
    if (!identical(dif$summary$uniform_DIF[j],
                   dif$terms$significant[u]) ||
        !identical(dif$summary$superseded[j],
                   dif$terms$superseded[u])) fail()
    if (length(nu)) {
      for (nm in names(map_non))
        if (!same_num(dif$summary[[nm]][j], dif$terms[[map_non[[nm]]]][nu]))
          fail()
      if (!identical(dif$summary$nonuniform_DIF[j],
                     dif$terms$significant[nu])) fail()
    } else {
      for (nm in names(map_non))
        if (!same_num(dif$summary[[nm]][j], NA_real_)) fail()
      if (!identical(dif$summary$nonuniform_DIF[j], FALSE)) fail()
    }
    if (identical(unit, "object")) {
      if (!same_num(dif$summary$min_judges[j],
                    dif$terms$min_judges[u]) ||
          !same_num(dif$summary$min_effective_judges[j],
                    dif$terms$min_effective_judges[u]) ||
          !identical(dif$summary$uniform_inference[j],
                     dif$terms$inference_available[u]) ||
          !identical(dif$summary$nonuniform_inference[j],
            if (length(nu)) dif$terms$inference_available[nu] else NA)) fail()
    }
  }
  invisible(NULL)
}

.validate_dif_result <- function(dif, fit) {
  if (is.null(dif)) return(invisible(NULL))
  if (!inherits(dif, "rasch_dif") || !is.data.frame(dif$summary) ||
      !is.data.frame(dif$terms))
    stop("`dif` must be a current dif_anova() result")
  if (inherits(fit, "rasch_btl"))
    stop("`dif` cannot be exported with a paired-comparison fit")
  if (is.null(dif$fit_signature))
    stop("`dif` predates fitted-model provenance; recompute it from this fit")
  if (!.fit_boot_signature_matches(dif$fit_signature, fit))
    stop("`dif` was computed from a different fitted model")
  .validate_primary_dif_tables(dif, "item", c("Residuals", "ci"), ":ci")
  if (length(dif$within) && !identical(dif$algorithm, "joint-between-1"))
    stop("`dif` predates joint adjustment for incomplete panels; recompute it",
         call. = FALSE)
  if (!.dif_followups_current(dif))
    stop("`dif` follow-ups may use earlier factor values; recompute the DIF analysis",
         call. = FALSE)
  invisible(NULL)
}

#' @export
print.rasch_dif <- function(x, ...) {
  s <- x$summary
  cat(sprintf("DIF by residual analysis of variance (%s; %d class intervals%s)\n",
              if (length(unique(s$term)) > length(unique(s$item)) ||
                  x$effects == "factorial")
                sprintf("%d terms, %s effects", length(unique(s$term)),
                        x$effects) else "one-way",
              x$n_groups[1],
              if (length(x$within))
                sprintf("; within-subject: %s", paste(x$within, collapse = ", "))
              else ""))
  cat("Uniform between-person terms use HC3 covariance; class-interval interactions retain the residual-ANOVA reference for complete panels.\n")
  if (grepl("CR3", x$between_covariance %||% "", fixed = TRUE))
    cat("Incomplete panels use joint adjustment and person-cluster CR3 covariance for all between-person terms.\n")
  show <- s[, c("item", "term", "F_uniform", "p_uniform_adj", "uniform_DIF",
                "F_nonuniform", "p_nonuniform_adj", "nonuniform_DIF")]
  print(.fmt_df(show), row.names = FALSE)
  cat(sprintf("%d uniform, %d non-uniform DIF flag(s) after %s adjustment.\n",
              sum(s$uniform_DIF, na.rm = TRUE),
              sum(s$nonuniform_DIF, na.rm = TRUE), x$p_adjust))
  # a blank row is a test the design could not estimate, not an absence of
  # DIF, and it still counts in the adjusted family: the count line above is
  # unreadable without the notes that say so
  if (length(x$notes)) cat("Notes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}

#' DIF differences between factor levels
#'
#' Resolves an item by one or more person factors and compares the resulting
#' locations. Several factors in \code{by} give pairwise comparisons between
#' their joint cells and can be used to quantify an interaction.
#'
#' @details
#' Let \eqn{\delta_i} contain the resolved locations of item
#' \eqn{i}, and let \eqn{\mathbf{c}_{ab}} place 1 on level \eqn{a}, -1 on
#' level \eqn{b}, and zero elsewhere. The reported difference and its
#' standard error are
#' \deqn{\Delta_{i,ab}=\mathbf{c}_{ab}^{\mathsf T}
#' \delta_i,}
#' \deqn{\operatorname{SE}(\Delta_{i,ab})=
#' \sqrt{\mathbf{c}_{ab}^{\mathsf T}\mathbf{V}_i
#' \mathbf{c}_{ab}},}
#' where \eqn{\mathbf{V}_i} is the full covariance of the resolved
#' locations. Wald probabilities are adjusted over the pairwise family. A
#' comparison with withheld inference remains in that declared family.
#'
#' When person identifiers repeat, the resolved-location covariance uses the
#' person-clustered calibration sandwich and a t reference with the number of
#' independent person clusters minus one degree of freedom. The same reference
#' is used for confidence intervals and the ETS interval-null probability.
#' This permits inference for the logit difference while allowing response
#' rows from the same person to be dependent. For a planned within-person
#' question, \code{\link{dif_contrasts}}
#' remains preferable because it tests the nominated contrast directly from
#' person-level residual contrasts. Inference is withheld if the resolved-
#' location covariance is unavailable or not positive semidefinite.
#'
#' @param fit A fitted object from \code{\link{rasch}} or
#'   \code{\link{rasch_mfrm}}. EFRM fits are excluded because an ordinary
#'   split refit would discard their frame units.
#' @param item Item name or index.
#' @param by One or more person-factor names nominated in the fit (several
#'   names give interaction cells), or a grouping vector/data frame with
#'   one entry per person.
#' @param p_adjust Adjustment over the pairwise comparisons. The default
#'   \code{"holm"} controls familywise error; use \code{"BH"} only for
#'   false-discovery-rate screening. \code{"none"} leaves probabilities
#'   unadjusted.
#' @param alpha Significance level for the adjusted probabilities.
#' @param flag_logits Absolute difference flagged as practically
#'   significant.
#' @param min_n Levels with fewer distinct responders to the item are dropped
#'   (their resolved locations would be too unstable to compare), with a note.
#'   When identifiers repeat, response rows from one person count once within
#'   each level.
#' @return A list of class \code{"rasch_dif_size"}. \code{levels} contains
#'   the resolved location, standard error and sample size for each level.
#'   \code{pairs} contains logit differences, Wald \code{t} statistics,
#'   confidence intervals, raw and adjusted probabilities, and practical
#'   flags. For
#'   dichotomous items it also contains \code{ets}, together with the raw and
#'   adjusted probabilities for exceeding the ETS A boundary; for polytomous
#'   items it contains the descriptive \code{signed_area}.
#'   \code{df} gives the reference degrees of freedom: infinite for independent
#'   response rows and the independent person-cluster count minus one for a
#'   supported repeated-person calibration.
#'   Sampling-uncertainty fields are \code{NA} when the resolved-location
#'   covariance cannot support Wald inference.
#'   The ETS category is also \code{NA} if a probability needed to classify
#'   the contrast is unavailable or its standard error is zero.
#'
#' @section Magnitude conventions:
#' For dichotomous items, \code{ets} applies the ETS A, B and C rules to the
#' itemwise comparison. On the logit scale the magnitude cut-points are
#' \eqn{1/2.35=0.426} and \eqn{1.5/2.35=0.638}. Category A also includes an
#' adjusted test that is not significant. Category C requires a magnitude of
#' at least 0.638 and rejection of the interval null
#' \eqn{|\Delta|\leq 0.426}; B is the remainder. Both probabilities are
#' adjusted by \code{p_adjust} over the requested pairwise family.
#'
#' For a partial credit item with \eqn{m_i} thresholds, the signed area
#' between the two expected-score curves has the closed form
#' \deqn{SA_{ab}=\int\{E_b(X\mid\theta)-E_a(X\mid\theta)\}\,d\theta
#' =\sum_{k=1}^{m_i}(\delta_{iak}-\delta_{ibk})
#' =m_i(\beta_{ia}-\beta_{ib}).}
#' This is returned as \code{signed_area}; a positive value means that level
#' \code{a} has the harder resolved item. It is descriptive and is not
#' given an A/B/C category: score-metric classifications for polytomous DIF
#' are not interchangeable with a PCM logit difference. For pooled MFRM
#' items, the areas use the same precision weight for a facet cell in every
#' group. A comparison is withheld when the groups do not support the same
#' observed response categories.
#' @references
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' Scandinavian Journal of Statistics, 6(2), 65--70.
#'
#' Zieky, M. (1993). Practical questions in the use of DIF statistics in item
#' development. In P. W. Holland and H. Wainer (eds), Differential Item
#' Functioning (pp. 337--364). Erlbaum.
#'
#' Linacre, J. M. and Wright, B. D. (1989). Mantel-Haenszel DIF and PROX are
#' equivalent! Rasch Measurement Transactions, 3(2), 51--53.
#'
#' Cohen, A. S., Kim, S.-H. and Baker, F. B. (1993). Detection of differential
#' item functioning in the graded response model. Applied Psychological
#' Measurement, 17(4), 335--350.
#'
#' Raju, N. S. (1988). The area between two item characteristic curves.
#' Psychometrika, 53(4), 495--502.
#' @seealso \code{\link{dif_anova}} and \code{\link{dif_contrasts}}.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' dif_size(fit, "I3", by = "grp")
#' @export
dif_size <- function(fit, item, by, p_adjust = "holm", alpha = 0.05,
                     flag_logits = 0.5, min_n = 20) {
  if (!is.atomic(item) || !is.null(dim(item)) || length(item) != 1L ||
      is.na(item))
    stop("`item` must name exactly one item")
  .check_dif_args(alpha, p_adjust, flag_logits, min_n)
  if (!inherits(fit, "rasch")) stop("dif_size needs a rasch fit")
  if (inherits(fit, "rasch_efrm"))
    .refuse("resolved DIF magnitudes are not available for EFRM fits; the ",
            "ordinary split refit would discard the fitted frame units")
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; DIF magnitudes are unavailable")
  mfrm_item <- inherits(fit, "rasch_mfrm") && !is.null(fit$virtual_map) &&
    !(item %in% colnames(fit$X)) && item %in% fit$virtual_map$item
  polytomous <- if (mfrm_item) {
    vi <- fit$virtual_map$vkey[fit$virtual_map$item == item]
    any(fit$m[match(vi, colnames(fit$X))] > 1L, na.rm = TRUE)
  } else FALSE
  if (!mfrm_item) {
    i <- .item_idx(fit, item)
    if (is.na(i)) stop("no such item")
    item <- fit$items$item[i]
    polytomous <- fit$m[i] > 1L
  }
  if (.role_columns(by,
                    if (is.null(fit$factors)) character(0) else names(fit$factors),
                    nrow(fit$X))) {
    bad <- if (is.null(fit$factors)) by else setdiff(by, names(fit$factors))
    if (length(bad))
      stop("not a person factor nominated in the fit: ",
           paste(bad, collapse = ", "))
  }
  factors <- .dif_factors(fit, by)
  .check_dif_factor_levels(factors)
  grp <- if (ncol(factors) == 1L) factor(factors[[1]])
         else .factor_cells(factors, sep = ":")
  notes <- character(0)
  # drop levels too thin on this item to resolve
  obs_i <- if (mfrm_item)
    rowSums(!is.na(fit$X[, fit$virtual_map$vkey[
      fit$virtual_map$item == item], drop = FALSE])) > 0L
  else !is.na(fit$X[, i])
  n_lev <- .dif_cell_n(grp, obs_i, fit$person$id)
  thin <- names(n_lev)[n_lev < min_n]
  if (length(thin)) {
    notes <- c(notes, sprintf("level(s) dropped with fewer than %d distinct responders: %s",
                              min_n, paste(thin, collapse = ", ")))
    grp <- factor(ifelse(as.character(grp) %in% thin, NA, as.character(grp)))
  }
  if (nlevels(droplevels(grp)) < 2)
    stop("fewer than two usable levels for item ", item)
  grp <- droplevels(grp)

  inference_df <- NA_real_
  if (mfrm_item) {
    # underlying MFRM item: pooled virtual-level resolution (one joint
    # unstructured refit of the virtual matrix; see .dif_resolve)
    rs <- .dif_resolve(fit, item, grp, min_n)
    if (is.null(rs))
      stop("could not resolve item ", item, " (too little data per level)")
    levs <- rs$levs; loc <- rs$loc; vloc <- rs$vloc; weak_lev <- rs$weak
    category_signature <- rs$category_signature
    area_level <- rs$area
    structure_ok <- !identical(rs$score_compatible, FALSE)
    inference_df <- rs$df %||% NA_real_
    notes <- c(notes, rs$notes)
  } else {
    levs <- levels(grp)
    category_signature <- matrix(vapply(levs, function(lv)
      paste(sort(unique(fit$X[as.character(grp) == lv &
                                !is.na(fit$X[, i]), i])), collapse = ","), ""),
      nrow = 1L, dimnames = list(item, levs))
    expected_signature <- paste(seq.int(0L, fit$m[i]), collapse = ",")
    structure_ok <- length(unique(category_signature[1L, ])) == 1L &&
      identical(unname(category_signature[1L, 1L]), expected_signature)
    if (!structure_ok) {
      # A common logit contrast cannot be obtained without renumbering at
      # least one group's scores. Do not perform that different refit merely
      # to discard its estimate below.
      loc <- rep(NA_real_, length(levs))
      vloc <- matrix(NA_real_, length(levs), length(levs))
      weak_lev <- stats::setNames(rep(FALSE, length(levs)), levs)
      area_level <- rep(NA_real_, length(levs))
    } else {
      refit <- split_items(fit, item, by = grp)
      inference_df <- .dif_refit_df(refit)
      split_names <- paste0(item, " (", levs, ")")
      idx <- match(split_names, refit$items$item)
      if (anyNA(idx))
        stop("resolved item(s) missing after re-analysis (too little data): ",
             paste(split_names[is.na(idx)], collapse = ", "))

      # location covariance from the sandwich: var(mean of a threshold block)
      thr <- refit$thresholds; cv <- refit$est$cov_tau
      block <- lapply(idx, function(k) thr$id[thr$item == k])
      loc <- refit$items$location[idx]
      area_level <- refit$m[idx] * loc
      vloc <- matrix(NA_real_, length(levs), length(levs))
      for (a in seq_along(levs)) for (b in seq_along(levs))
        vloc[a, b] <- mean(cv[block[[a]], block[[b]], drop = FALSE])
      weak_lev <- .dif_weak_levels(refit, as.list(idx)); names(weak_lev) <- levs
      if (any(weak_lev))
        notes <- c(notes, sprintf(
          "location(s) for level(s) %s rest on a near-empty category and are weakly identified; their DIF magnitude and significance are withheld",
          paste(levs[weak_lev], collapse = ", ")))
    }
  }
  vloc_ok <- length(dim(vloc)) == 2L &&
    identical(dim(vloc), c(length(levs), length(levs))) &&
    all(is.finite(vloc)) && .covariance_is_symmetric(vloc) &&
    .covariance_is_psd(vloc)
  reference_ok <- !is.na(inference_df) && inference_df > 0
  covariance_bad <- structure_ok && !vloc_ok
  if (covariance_bad) {
    consequence <- if (any(is.finite(loc))) paste(
      "point differences remain descriptive, but standard errors, confidence",
      "intervals and Wald tests are withheld") else paste(
        "the covariance-weighted point differences and all associated",
        "inference are withheld")
    notes <- c(notes, paste(
      "the resolved-location covariance is unavailable or not positive",
      "semidefinite:", consequence))
  }
  if (structure_ok && vloc_ok && !reference_ok)
    notes <- c(notes, paste(
      "the resolved calibration does not carry enough independent-person",
      "support for covariance-based inference; point differences remain",
      "descriptive, but standard errors, confidence intervals and Wald tests",
      "are withheld"))
  n_item <- as.integer(.dif_cell_n(grp, obs_i, fit$person$id)[levs])
  lev_se <- if (vloc_ok && reference_ok) sqrt(pmax(diag(vloc), 0)) else
    rep(NA_real_, length(levs))
  lev_se[weak_lev] <- NA_real_
  if (covariance_bad) lev_se[] <- NA_real_
  levels_df <- data.frame(level = levs, location = loc,
                          se = lev_se, weak = unname(weak_lev), n = n_item)

  pr <- t(utils::combn(seq_along(levs), 2))
  pair_weak <- weak_lev[pr[, 1]] | weak_lev[pr[, 2]]
  same_categories <- vapply(seq_len(nrow(pr)), function(k)
    all(category_signature[, pr[k, 1]] ==
          category_signature[, pr[k, 2]]), TRUE)
  pair_invalid <- pair_weak | !same_categories | !structure_ok
  if (any(!same_categories) || !structure_ok) {
    bad <- !same_categories | !structure_ok
    bad_pairs <- paste0(levs[pr[bad, 1]], " versus ",
                        levs[pr[bad, 2]])
    notes <- c(notes, paste(
      "DIF magnitude and inference are withheld where resolved groups have",
      "different observed response-category structures or do not retain",
      "the fitted 0:m score structure:",
      paste(bad_pairs, collapse = ", ")))
  }
  pair_var <- if (vloc_ok)
    diag(vloc)[pr[, 1]] + diag(vloc)[pr[, 2]] -
      2 * vloc[cbind(pr[, 1], pr[, 2])] else
    rep(NA_real_, nrow(pr))
  pairs <- data.frame(
    level_a = levs[pr[, 1]], level_b = levs[pr[, 2]],
    difference = loc[pr[, 1]] - loc[pr[, 2]],
    se = if (reference_ok) sqrt(pmax(pair_var, 0)) else
      rep(NA_real_, length(pair_var)))
  # a pair touching a weakly-identified level carries no trustworthy
  # magnitude: withhold its SE and every SE-derived verdict
  pairs$difference[pair_invalid] <- NA_real_
  pairs$se[pair_invalid] <- NA_real_
  if (covariance_bad) pairs$se[] <- NA_real_
  pairs$t <- .wald_ratio(pairs$difference, pairs$se)
  pairs$df <- if (vloc_ok && reference_ok) rep(inference_df, nrow(pairs)) else
    rep(NA_real_, nrow(pairs))
  pairs$p <- 2 * stats::pt(-abs(pairs$t), df = pairs$df)
  pairs$p_adj <- .p_adjust_family(pairs$p, method = p_adjust)
  critical <- stats::qt(0.975, df = pairs$df)
  pairs$lower <- pairs$difference - critical * pairs$se
  pairs$upper <- pairs$difference + critical * pairs$se
  pairs$significant <- ifelse(pair_invalid | covariance_bad, NA,
                              pairs$p_adj < alpha)
  pairs$practical <- ifelse(pair_invalid, NA,
                            abs(pairs$difference) >= flag_logits)
  if (!polytomous) {
    pairs$p_beyond_A <- .ets_p_beyond(pairs$difference, pairs$se, pairs$df)
    pairs$p_beyond_A_adj <- .p_adjust_family(
      pairs$p_beyond_A, method = p_adjust)
    pairs$ets <- .ets_category(pairs$difference, pairs$se, pairs$p_adj,
                               alpha, pairs$p_beyond_A_adj)
  } else {
    pairs$p_beyond_A <- pairs$p_beyond_A_adj <- NA_real_
    pairs$ets <- NA_character_
  }
  pairs$signed_area <- if (polytomous)
    area_level[pr[, 1]] - area_level[pr[, 2]] else NA_real_
  pairs$signed_area[pair_invalid] <- NA_real_

  out <- list(item = item, by = paste(names(factors), collapse = ":"),
              levels = levels_df, pairs = pairs, alpha = alpha,
              p_adjust = p_adjust, flag_logits = flag_logits,
              classification = if (polytomous)
                "PCM signed expected-score area (descriptive)" else "ETS",
              notes = notes)
  out <- .tag_tables(out)
  class(out) <- "rasch_dif_size"
  out
}

#' @export
print.rasch_dif_size <- function(x, ...) {
  cat(sprintf("DIF size for %s by %s (resolved locations, logits)\n",
              x$item, x$by))
  lv <- x$levels; lv[-1] <- lapply(lv[-1], round, 3)
  print(lv, row.names = FALSE)
  pr <- x$pairs
  num <- vapply(pr, is.numeric, TRUE)
  pr[num] <- lapply(pr[num], round, 3)
  pr$significant <- ifelse(pr$significant, "*", "")
  pr$practical <- ifelse(pr$practical, sprintf(">= %.2f", x$flag_logits), "")
  print(pr, row.names = FALSE)
  cat(sprintf("p adjusted by %s over %d pairwise comparison(s); practical criterion %.2f logits\n",
              x$p_adjust, nrow(pr), x$flag_logits))
  if (length(x$notes)) cat("notes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Planned contrasts: the confirmatory alternative to exhaustive post-hoc
# pairwise comparison. The family of questions is derived from the structure
# of the nominated factors (or supplied), estimated on the logit scale from
# resolved item locations, and -- when persons repeat across rows of a
# stacked design -- tested from person-level contrast scores so that
# within-subject dependence is respected.
# ---------------------------------------------------------------------------

# Resolve one item over grouping cells: locations and sandwich covariance.
# A resolved level is weakly identified when its split copy rests on a
# near-empty category: split_items() already flags such thresholds
# (weak = TRUE, se = NA) and leaves the item-location SE NA. A location
# built on such a threshold is a boundary artefact, so its DIF magnitude
# and significance must be withheld rather than recomputed from the ridged
# covariance -- otherwise dif_size()/dif_contrasts() report a fabricated
# finite SE and a spurious 'significant'/'practical' verdict. item_rows is
# a list, one entry per level, of the refit$items row-index(es) whose
# thresholds back that level's location.
.dif_weak_levels <- function(refit, item_rows) {
  wk <- refit$thresholds$weak
  se <- refit$items$se
  vapply(item_rows, function(ks) any(vapply(ks, function(k) {
    (!is.null(wk) && isTRUE(any(wk[refit$thresholds$item == k], na.rm = TRUE))) ||
      (!is.null(se) && k <= length(se) && is.na(se[k]))
  }, logical(1))), logical(1))
}

# Degrees of freedom for a resolved-location covariance. Independent response
# rows retain the asymptotic normal reference (df = Inf). With repeated person
# identifiers, the independent sampling units are the clusters that support
# the accepted sandwich covariance. A malformed or withheld sandwich,
# whether based on independent persons or repeated-person clusters, cannot
# acquire a reference distribution downstream.
.dif_refit_df <- function(refit) {
  support <- refit$est$cluster_support
  # The sandwich records repetition only among clusters that contribute
  # conditional item-pair information. Raw duplicate IDs on all-missing or
  # otherwise uninformative rows must not turn an independent calibration
  # into a finite-cluster analysis.
  if (is.list(support) &&
      (identical(support$repeated, FALSE) || isTRUE(support$repeated))) {
    if (!isTRUE(refit$est$cluster_inference)) return(NA_real_)
    if (identical(support$repeated, FALSE)) return(Inf)
    if (
        !is.numeric(support$n) || length(support$n) != 1L ||
        !is.finite(support$n) || support$n < 2L)
      return(NA_real_)
    return(as.numeric(support$n - 1L))
  }
  # A legacy fit may have raw repeated identifiers but no recorded support.
  # It cannot justify either a cluster df or an asymptotic independent-row
  # reference, so withhold rather than guessing from all response rows.
  repeated_legacy <- isTRUE(refit$repeated_ids) ||
    (!is.null(refit$person$id) && .has_repeated_person_ids(refit$person$id))
  if (repeated_legacy || !is.null(support) ||
      !isTRUE(refit$est$cluster_inference))
    return(NA_real_)
  Inf
}

.dif_resolve <- function(fit, item, grp, min_n) {
  # an UNDERLYING MFRM item resolves at the virtual level: every one of
  # its facet cells is split by the groups in one joint unstructured
  # refit of the virtual matrix (the facet decomposition is not
  # reimposed), and the per-level locations are precision-weighted means
  # over the item's cells with the full covariance carried
  if (inherits(fit, "rasch_mfrm") && !is.null(fit$virtual_map) &&
      !(item %in% colnames(fit$X)) && item %in% fit$virtual_map$item) {
    vm <- fit$virtual_map
    cols <- vm$vkey[vm$item == item]
    notes <- paste0(item, ": resolved at the virtual-item level, pooled ",
                    "over its facet cells (facet structure not reimposed)")
    obs <- rowSums(!is.na(fit$X[, cols, drop = FALSE])) > 0L
    n_lev <- .dif_cell_n(grp, obs, fit$person$id)
    thin <- names(n_lev)[n_lev < min_n]
    if (length(thin)) {
      notes <- c(notes, sprintf(
        "%s: level(s) dropped with fewer than %d distinct responders: %s",
        item, min_n, paste(thin, collapse = ", ")))
      grp <- factor(ifelse(as.character(grp) %in% thin, NA,
                           as.character(grp)))
    }
    grp <- droplevels(grp)
    if (nlevels(grp) < 2) return(NULL)
    levs <- levels(grp)
    category_all <- vapply(levs, function(lv)
      vapply(cols, function(cc)
        paste(sort(unique(fit$X[as.character(grp) == lv &
                                  !is.na(fit$X[, cc]), cc])),
              collapse = ","), ""), character(length(cols)))
    if (is.null(dim(category_all)))
      category_all <- matrix(category_all, nrow = length(cols))
    rownames(category_all) <- cols; colnames(category_all) <- levs
    expected <- vapply(cols, function(cc)
      paste(seq.int(0L, fit$m[match(cc, colnames(fit$X))]), collapse = ","), "")
    compatible <- vapply(seq_along(cols), function(j)
      length(unique(category_all[j, ])) == 1L &&
        identical(unname(category_all[j, 1L]), unname(expected[j])), TRUE)
    if (!any(compatible)) {
      notes <- c(notes, paste0(
        item, ": resolved contrasts withheld because none of its facet ",
        "cells has the same observed response-category structure in every ",
        "group"))
      return(list(
        levs = levs, loc = rep(NA_real_, length(levs)),
        vloc = matrix(NA_real_, length(levs), length(levs)),
        weak = stats::setNames(rep(FALSE, length(levs)), levs),
        m_cell = matrix(NA_real_, nrow(category_all), length(levs),
                        dimnames = dimnames(category_all)),
        category_signature = category_all,
        score_compatible = FALSE,
        area = rep(NA_real_, length(levs)), df = NA_real_, notes = notes))
    }
    if (any(!compatible))
      notes <- c(notes, sprintf(
        "%s: facet cell(s) omitted because groups used different observed response categories: %s",
        item, paste(cols[!compatible], collapse = ", ")))
    cols <- cols[compatible]
    category_all <- category_all[compatible, , drop = FALSE]
    vfit <- fit; class(vfit) <- "rasch"
    vfit$model <- "PCM"   # unstructured virtual thresholds refit as PCM
    refit <- tryCatch(split_items(vfit, cols, by = grp),
                      error = function(e) NULL)
    if (is.null(refit)) return(NULL)
    thr <- refit$thresholds; cv <- refit$est$cov_tau
    if (!.covariance_supports_wald(cv, nrow(thr))) {
      notes <- c(notes, paste0(
        item, ": pooled facet-cell DIF magnitudes are withheld because the ",
        "resolved-threshold covariance is unavailable, asymmetric, or not ",
        "positive semidefinite and therefore cannot define the pooling weights"))
      return(list(
        levs = levs, loc = rep(NA_real_, length(levs)),
        vloc = matrix(NA_real_, length(levs), length(levs)),
        weak = stats::setNames(rep(FALSE, length(levs)), levs),
        m_cell = matrix(NA_real_, nrow(category_all), length(levs),
                        dimnames = dimnames(category_all)),
        category_signature = category_all,
        score_compatible = TRUE, area = rep(NA_real_, length(levs)),
        df = NA_real_, notes = notes))
    }
    # COMMON cells with COMMON weights: every facet cell used must be
    # resolved for EVERY level, and each cell gets one weight shared by
    # all levels, so the cell's facet severity cancels exactly from every
    # level contrast. Group-specific precision weights let severity leak
    # into the DIF magnitude when groups have different facet exposure
    # (a no-DIF design with sex-linked rater allocation read -1.75
    # logits, z = -6.5).
    idx_m <- sapply(levs, function(l)
      match(paste0(cols, " (", l, ")"), refit$items$item))
    if (is.null(dim(idx_m))) idx_m <- matrix(idx_m, nrow = length(cols))
    common <- rowSums(is.na(idx_m)) == 0L
    if (sum(common) < 1L) return(NULL)
    if (any(!common))
      notes <- c(notes, sprintf(
        "%s: facet cell(s) dropped from the magnitude (not resolvable for every level): %s",
        item, paste(cols[!common], collapse = ", ")))
    idx_m <- idx_m[common, , drop = FALSE]
    cols_common <- cols[common]
    category_signature <- category_all[cols_common, , drop = FALSE]
    blocks <- lapply(seq_len(nrow(idx_m)), function(ci)
      lapply(idx_m[ci, ], function(k) thr$id[thr$item == k]))
    # one weight per cell: inverse of the level-averaged location variance
    vr_c <- vapply(blocks, function(bl)
      mean(vapply(bl, function(rws) mean(cv[rws, rws]), 0)), 0)
    # PSD has already been established. Any remaining negative value can only
    # be numerical noise at that matrix's scale.
    w_c <- .inverse_variance_weights(pmax(vr_c, 0))
    w_c <- w_c / sum(w_c)
    cell_loc <- matrix(refit$items$location[idx_m], nrow = nrow(idx_m),
                       ncol = length(levs), dimnames = list(NULL, levs))
    m_cell <- matrix(refit$m[idx_m], nrow = nrow(idx_m),
                     ncol = length(levs), dimnames = list(NULL, levs))
    loc <- vapply(seq_along(levs), function(a)
      sum(w_c * cell_loc[, a]), 0)
    area <- vapply(seq_along(levs), function(a)
      sum(w_c * m_cell[, a] * cell_loc[, a]), 0)
    vloc <- matrix(NA_real_, length(levs), length(levs))
    for (a in seq_along(levs)) for (b in seq_along(levs)) {
      acc <- 0
      for (ca in seq_along(blocks)) for (cb in seq_along(blocks))
        acc <- acc + w_c[ca] * w_c[cb] *
          mean(cv[blocks[[ca]][[a]], blocks[[cb]][[b]], drop = FALSE])
      vloc[a, b] <- acc
    }
    weak_lev <- .dif_weak_levels(refit, lapply(seq_along(levs),
                                               function(a) idx_m[, a]))
    names(weak_lev) <- levs
    if (any(weak_lev))
      notes <- c(notes, sprintf(
        "%s: location(s) for level(s) %s rest on a near-empty category and are weakly identified; their DIF magnitude and significance are withheld",
        item, paste(levs[weak_lev], collapse = ", ")))
    return(list(levs = levs, loc = loc, vloc = vloc, weak = weak_lev,
                m_cell = m_cell, category_signature = category_signature,
                score_compatible = TRUE, area = area,
                df = .dif_refit_df(refit), notes = notes))
  }
  i <- .item_idx(fit, item)
  item <- fit$items$item[i]
  notes <- character(0)
  n_lev <- .dif_cell_n(grp, !is.na(fit$X[, i]), fit$person$id)
  thin <- names(n_lev)[n_lev < min_n]
  if (length(thin)) {
    notes <- c(notes, sprintf(
      "%s: level(s) dropped with fewer than %d distinct responders: %s",
      item, min_n, paste(thin, collapse = ", ")))
    grp <- factor(ifelse(as.character(grp) %in% thin, NA, as.character(grp)))
  }
  grp <- droplevels(grp)
  if (nlevels(grp) < 2) return(NULL)
  levs <- levels(grp)
  category_signature <- matrix(vapply(levs, function(lv)
    paste(sort(unique(fit$X[as.character(grp) == lv &
                              !is.na(fit$X[, i]), i])), collapse = ","), ""),
    nrow = 1L, dimnames = list(item, levs))
  expected_signature <- paste(seq.int(0L, fit$m[i]), collapse = ",")
  structure_ok <- length(unique(category_signature[1L, ])) == 1L &&
    identical(unname(category_signature[1L, 1L]), expected_signature)
  if (!structure_ok) {
    notes <- c(notes, paste(
      item, ": resolved contrasts withheld because groups have different",
      "observed response-category structures"))
    return(list(
      levs = levs, loc = rep(NA_real_, length(levs)),
      vloc = matrix(NA_real_, length(levs), length(levs)),
      weak = stats::setNames(rep(FALSE, length(levs)), levs),
      m_cell = matrix(NA_real_, 1L, length(levs),
                      dimnames = list(item, levs)),
      category_signature = category_signature,
      score_compatible = FALSE,
      area = rep(NA_real_, length(levs)), df = NA_real_, notes = notes))
  }
  # A split refit can fail on this item alone (it cannot preserve the score
  # structure, or it does not converge). Withhold this item's resolution
  # with the reason instead of aborting a whole multi-item analysis; the
  # virtual-item branch above already behaves this way.
  refit <- tryCatch(split_items(fit, item, by = grp),
                    error = function(e) e)
  # The refusal can have any of several causes (an anchored item, a lost
  # score category, no convergence), so carry its own message and mark the
  # refit as never run. The cells' category structures were checked just
  # above and agree; a caller must not report this as a category problem.
  if (inherits(refit, "error"))
    return(list(
      levs = levs, loc = rep(NA_real_, length(levs)),
      vloc = matrix(NA_real_, length(levs), length(levs)),
      weak = stats::setNames(rep(FALSE, length(levs)), levs),
      m_cell = matrix(NA_real_, 1L, length(levs),
                      dimnames = list(item, levs)),
      category_signature = category_signature,
      score_compatible = FALSE, refit_error = conditionMessage(refit),
      area = rep(NA_real_, length(levs)), df = NA_real_,
      notes = c(notes, paste0(item, ": resolved contrasts withheld because ",
                              "the split refit is unavailable: ",
                              conditionMessage(refit)))))
  idx <- match(paste0(item, " (", levs, ")"), refit$items$item)
  if (anyNA(idx)) return(NULL)
  thr <- refit$thresholds; cv <- refit$est$cov_tau
  block <- lapply(idx, function(k) thr$id[thr$item == k])
  loc <- refit$items$location[idx]
  m_cell <- matrix(refit$m[idx], nrow = 1L,
                   dimnames = list(item, levs))
  area <- refit$m[idx] * loc
  vloc <- matrix(NA_real_, length(levs), length(levs))
  for (a in seq_along(levs)) for (b in seq_along(levs))
    vloc[a, b] <- mean(cv[block[[a]], block[[b]], drop = FALSE])
  weak_lev <- .dif_weak_levels(refit, as.list(idx)); names(weak_lev) <- levs
  if (any(weak_lev))
    notes <- c(notes, sprintf(
      "%s: location(s) for level(s) %s rest on a near-empty category and are weakly identified; their DIF magnitude and significance are withheld",
      item, paste(levs[weak_lev], collapse = ", ")))
  list(levs = levs, loc = loc, vloc = vloc, weak = weak_lev,
       m_cell = m_cell, category_signature = category_signature,
       score_compatible = TRUE, area = area,
       df = .dif_refit_df(refit), notes = notes)
}

# contr.poly returns a level that carries no weight in a polynomial
# contrast as a rounding residue (the middle level of an odd-K linear
# trend is about -1e-17, not 0). The support rules read any non-zero
# weight as a required cell, so leave only the weights the contrast
# actually places.
.dif_poly_weights <- function(cp, j, levs) {
  w <- cp[, j]
  w[abs(w) < 1e-10 * max(abs(w))] <- 0
  names(w) <- levs
  w
}

# A factor is treated as ordered when declared ordered or when its levels
# parse as numbers (ages, waves, doses).
.dif_is_ordered <- function(f)
  is.ordered(f) || !any(is.na(suppressWarnings(as.numeric(levels(f)))))

# The leading contrast of a factor: the difference for two levels, the
# linear trend for an ordered factor, none for a nominal many-level factor.
.dif_leading <- function(f) {
  K <- nlevels(f)
  if (K == 2L) {
    w <- c(-1, 1); names(w) <- levels(f)
    list(weights = w,
         label = sprintf("%s - %s", levels(f)[2], levels(f)[1]))
  } else if (.dif_is_ordered(f)) {
    sc <- suppressWarnings(as.numeric(levels(f)))
    cp <- if (!any(is.na(sc))) stats::contr.poly(K, scores = sc)
          else stats::contr.poly(K)
    w <- .dif_poly_weights(cp, 1L, levels(f))
    list(weights = w, label = "linear")
  } else NULL
}

# The planned questions a single factor admits.
.dif_factor_contrasts <- function(f, fname) {
  K <- nlevels(f); out <- list(); labels <- character(0)
  add <- function(label, weights) {
    labels <<- c(labels, label)
    out[[length(out) + 1L]] <<- weights
  }
  if (K == 2L) {
    lead <- .dif_leading(f)
    add(sprintf("%s: %s", fname, lead$label), lead$weights)
  } else if (.dif_is_ordered(f)) {
    sc <- suppressWarnings(as.numeric(levels(f)))
    cp <- if (!any(is.na(sc))) stats::contr.poly(K, scores = sc)
          else stats::contr.poly(K)
    w1 <- .dif_poly_weights(cp, 1L, levels(f))
    add(sprintf("%s: linear", fname), w1)
    w2 <- .dif_poly_weights(cp, 2L, levels(f))
    add(sprintf("%s: quadratic", fname), w2)
  } else if (K <= 4L) {
    pr <- utils::combn(levels(f), 2)
    for (j in seq_len(ncol(pr))) {
      w <- stats::setNames(numeric(K), levels(f))
      w[pr[2, j]] <- 1; w[pr[1, j]] <- -1
      add(sprintf("%s: %s - %s", fname, pr[2, j], pr[1, j]), w)
    }
  } else {
    for (l in levels(f)) {
      w <- stats::setNames(rep(-1 / (K - 1), K), levels(f)); w[l] <- 1
      add(sprintf("%s: %s - others", fname, l), w)
    }
  }
  names(out) <- .dif_unique_labels(labels)
  out
}

# Scale cell weights so the positive and negative parts each sum to one:
# every contrast then reads as a difference between two weighted averages,
# in logits, comparable across contrasts and against the practical flag.
.dif_norm <- function(w) {
  w[is.na(w)] <- 0
  scale <- max(abs(w))
  if (!is.finite(scale) || scale == 0) return(NULL)
  w <- w / scale
  w * 2 / sum(abs(w))
}

# Readable labels can collide when a level itself contains the display
# separator. Keep ordinary labels unchanged, but tag every member of a
# collision so assigning it to a named list cannot replace another question.
.dif_unique_labels <- function(x) {
  x <- as.character(x)
  clash <- duplicated(x) | duplicated(x, fromLast = TRUE)
  if (any(clash))
    x[clash] <- paste0(x[clash], " [comparison ", which(clash), "]")
  if (anyDuplicated(x)) x <- paste0("comparison ", seq_along(x), ": ", x)
  x
}

# Marginal contrasts use only nuisance-factor strata containing every
# non-zero target cell. Averaging the two sides over different observed
# strata changes the estimand and can confound it with a nuisance factor.
.dif_common_support_weights <- function(cellmap, fweights) {
  target <- names(fweights)
  nuisance <- setdiff(names(cellmap), c("cell", target))
  raw <- rep(1, nrow(cellmap)); active <- rep(TRUE, nrow(cellmap))
  for (fn in target) {
    cw <- unname(fweights[[fn]][as.character(cellmap[[fn]])])
    cw[is.na(cw)] <- 0
    raw <- raw * cw
    active <- active & cw != 0
  }
  nkey <- if (length(nuisance))
    as.character(.factor_cells(cellmap[nuisance], sep = "\r"))
  else rep("all", nrow(cellmap))
  required <- prod(vapply(fweights, function(w) sum(w != 0), integer(1)))
  complete_n <- names(which(tapply(active, nkey, sum) == required))
  if (!length(complete_n)) return(NULL)
  raw[!active | !nkey %in% complete_n] <- 0
  .dif_norm(stats::setNames(raw / length(complete_n), cellmap$cell))
}

# Derive the planned family from the factor structure.
.dif_contrast_family <- function(factors, cellmap, within_names) {
  fam <- list(); meta <- list(); labels <- character(0); estimable <- logical(0)
  add <- function(label, weights, metadata) {
    labels <<- c(labels, label)
    estimable <<- c(estimable, !is.null(weights))
    if (!is.null(weights)) {
      fam[[length(fam) + 1L]] <<- weights
      meta[[length(meta) + 1L]] <<- metadata
    }
  }
  for (fname in names(factors)) {
    fc <- .dif_factor_contrasts(factors[[fname]], fname)
    for (nm in names(fc)) {
      w <- .dif_common_support_weights(
        cellmap, stats::setNames(list(fc[[nm]]), fname))
      add(nm, w, list(factors = fname, fweights = fc[nm],
                      within = fname %in% within_names))
    }
  }
  fns <- names(factors)
  if (length(fns) >= 2) for (a in seq_len(length(fns) - 1))
    for (b in seq(a + 1, length(fns))) {
      la <- .dif_leading(factors[[fns[a]]])
      lb <- .dif_leading(factors[[fns[b]]])
      if (is.null(la) || is.null(lb)) next
      w <- .dif_common_support_weights(
        cellmap,
        stats::setNames(list(la$weights, lb$weights), c(fns[a], fns[b])))
      nm <- sprintf("%s(%s) x %s(%s)", fns[a], la$label, fns[b], lb$label)
      add(nm, w, list(factors = c(fns[a], fns[b]),
                      fweights = list(la$weights, lb$weights),
                      within = any(c(fns[a], fns[b]) %in% within_names)))
    }
  planned <- .dif_unique_labels(labels)
  names(fam) <- names(meta) <- planned[estimable]
  list(family = fam, meta = meta, planned = planned,
       planned_n = length(planned))
}

# Test any resolved-cell contrast in a stacked design without treating rows
# from the same person as independent. The supplied cell weights already
# encode the desired marginal comparison. Within each between-person cell we
# first form one weighted residual score per person, then combine the
# independent cell means with a Welch--Satterthwaite reference.
.dif_paired_cell_contrast <- function(z, factors, grp, id, within,
                                      cellmap, weights) {
  id <- .dif_ids(id)
  between <- setdiff(names(factors), within)
  bkey <- if (length(between))
    as.character(.factor_cells(factors[between], sep = "\r"))
  else rep("all", nrow(factors))
  map_bkey <- if (length(between))
    as.character(.factor_cells(cellmap[between], sep = "\r"))
  else rep("all", nrow(cellmap))
  names(weights) <- cellmap$cell
  required_bkeys <- unique(map_bkey[is.finite(weights) & weights != 0])
  if (!length(required_bkeys)) return(NULL)

  people <- split(seq_along(id), id)
  group <- rep(NA_character_, length(people))
  score_num <- rep(NA_real_, length(people))
  for (j in seq_along(people)) {
    r <- people[[j]]
    ok <- is.finite(z[r]) & !is.na(grp[r]) & !is.na(bkey[r])
    if (!any(ok)) next
    r <- r[ok]
    bg <- unique(bkey[r])
    if (length(bg) != 1L) next
    use_cells <- cellmap$cell[map_bkey == bg & weights != 0]
    if (!length(use_cells)) next
    means <- tapply(z[r], as.character(grp[r]), mean)
    if (anyNA(match(use_cells, names(means)))) next
    score_num[j] <- sum(weights[use_cells] * means[use_cells])
    group[j] <- bg
  }
  ok <- is.finite(score_num) & !is.na(group)
  if (sum(ok) < 3L) return(NULL)
  sp <- split(score_num[ok], group[ok])
  # The cell weights define an equal-stratum marginal estimand. Dropping a
  # required between-person stratum merely because its variance cannot be
  # estimated would change that estimand to the surviving strata. Withhold
  # the paired test unless every weighted stratum supplies at least two
  # complete person contrasts.
  if (!all(required_bkeys %in% names(sp))) return(NULL)
  sp <- sp[required_bkeys]
  if (any(lengths(sp) < 2L)) return(NULL)
  means <- vapply(sp, mean, 0)
  vars <- vapply(sp, stats::var, 0)
  ns <- lengths(sp)
  parts <- vars / ns
  vv <- sum(parts)
  if (!is.finite(vv) || vv <= 0) return(NULL)
  df <- vv^2 / sum(parts^2 / (ns - 1))
  stat <- sum(means) / sqrt(vv)
  list(stat = stat, df = df, p = 2 * stats::pt(-abs(stat), df))
}

# Pairwise marginal comparisons for one term. For a main effect these are
# differences between factor levels, averaged equally over complete cells of
# the remaining factors. For an interaction they are tensor products of the
# level differences: difference-in-differences for two factors and the direct
# higher-order analogue beyond two.
.dif_posthoc_family <- function(factors, cellmap, target, within) {
  bad <- setdiff(target, names(factors))
  if (length(bad))
    stop("term factor(s) not found: ", paste(bad, collapse = ", "))
  pairs <- lapply(target, function(fn) {
    lv <- levels(factors[[fn]])
    if (length(lv) < 2L) return(list())
    pr <- utils::combn(lv, 2)
    lapply(seq_len(ncol(pr)), function(j) {
      w <- stats::setNames(numeric(length(lv)), lv)
      w[pr[1, j]] <- -1; w[pr[2, j]] <- 1
      list(weights = w, label = sprintf("%s - %s", pr[2, j], pr[1, j]))
    })
  })
  if (any(!lengths(pairs))) stop("every term factor needs at least two levels")
  grid <- expand.grid(lapply(pairs, seq_along), KEEP.OUT.ATTRS = FALSE)
  nuisance <- setdiff(names(factors), target)
  nkey <- if (length(nuisance))
    as.character(.factor_cells(cellmap[nuisance], sep = "\r"))
  else rep("all", nrow(cellmap))
  family <- meta <- list()
  chosen_grid <- lapply(seq_len(nrow(grid)), function(r)
    lapply(seq_along(target), function(j) pairs[[j]][[grid[r, j]]]))
  planned <- .dif_unique_labels(vapply(chosen_grid, function(chosen)
    paste(vapply(chosen, `[[`, "", "label"), collapse = " x "), ""))
  for (r in seq_len(nrow(grid))) {
    chosen <- chosen_grid[[r]]
    label <- planned[r]
    raw <- rep(1, nrow(cellmap))
    active <- rep(TRUE, nrow(cellmap))
    for (j in seq_along(target)) {
      fw <- chosen[[j]]$weights
      cw <- unname(fw[as.character(cellmap[[target[j]]])])
      cw[is.na(cw)] <- 0
      raw <- raw * cw
      active <- active & cw != 0
    }
    # Marginalise only over nuisance strata containing the complete target
    # contrast. This avoids changing the estimand when an unbalanced design
    # has a structurally absent target cell.
    complete_n <- names(which(tapply(active, nkey, sum) == 2^length(target)))
    raw[!active | !nkey %in% complete_n] <- 0
    if (!length(complete_n)) next
    w <- stats::setNames(raw / length(complete_n), cellmap$cell)
    family[[label]] <- w
    meta[[label]] <- list(
      factors = target,
      fweights = lapply(chosen, `[[`, "weights"),
      within = any(target %in% within))
  }
  # Keep the predeclared family size even when incomplete support makes one
  # of its questions unestimable. Otherwise the remaining Holm probabilities
  # depend on which factorial cells happened to be observed.
  list(family = family, meta = meta, planned = planned,
       planned_n = length(planned))
}

#' Planned DIF contrasts
#'
#' Tests a specified family of one-degree-of-freedom DIF contrasts. By default,
#' contrasts are derived from the factor structure: differences for two-level
#' factors, polynomial trends for ordered factors, and pairwise or
#' level-against-rest comparisons for nominal factors. Leading contrasts are
#' crossed to form two-factor interactions. User-supplied cell weights are
#' also accepted.
#'
#' @details
#' Each logit contrast is calculated from resolved item locations. Weights are
#' scaled so their positive and negative parts each sum to one. With repeated
#' persons, inference uses person-level residual contrast scores with the same
#' cell weights as the resolved estimate. Nuisance-factor cells are averaged
#' equally rather than in proportion to their sample sizes. In an incomplete
#' factorial design, a contrast uses only nuisance strata containing all of
#' its non-zero target cells; an unsupported planned contrast is not estimated
#' but remains in the multiplicity count. Once these weights are defined,
#' every weighted cell must meet \code{min_n} for the item; sparse cells are
#' not dropped and the remaining weights are not renormalised. A level a
#' contrast places no weight on is not required: the middle level of an
#' odd-length linear trend carries weight zero, and so restricts neither the
#' nuisance strata nor the persons the test uses. Independent between-person cells are
#' then combined with a Welch--Satterthwaite reference. If a required between-
#' person cell has fewer than two complete person scores, residual inference is
#' withheld rather than changing the marginal contrast by dropping that cell.
#' The resolved logit
#' estimate is retained, but its
#' calibration-based standard error is withheld because it does not include
#' repeated-person dependence.
#'
#' For independent rows, a contrast with weights \eqn{\mathbf{c}} is
#' \deqn{\Delta_i=\mathbf{c}^{\mathsf T}\delta_i,\qquad
#' \operatorname{SE}(\Delta_i)=
#' \sqrt{\mathbf{c}^{\mathsf T}\mathbf{V}_i\mathbf{c}}.}
#' In a repeated-measures design, a within-person contrast is formed from the
#' standardised residuals,
#' \deqn{s_p=\sum_l c_l z_{pl},}
#' and tested over persons. The complete cell-weight vector is retained for
#' main effects and interactions, so the residual test and resolved estimate
#' address the same marginal contrast. The sign of each residual test is
#' aligned with the resolved logit contrast. Contrasts require a converged
#' calibration. For independent rows, an unavailable or non-positive-
#' semidefinite resolved-location covariance leaves the logit estimate
#' descriptive and causes Wald inference to be withheld.
#' A contrast with withheld inference remains in the adjustment family formed
#' by every requested item and contrast. When the split refit that resolves one
#' item's locations cannot be calibrated, that item's contrasts are withheld
#' with the reason and the other items are unaffected.
#' For an MFRM fit, underlying items are pooled over their facet cells by
#' default. EFRM fits are excluded because the required split refit would
#' discard the frame units.
#'
#' @param fit A fitted object from \code{\link{rasch}} or
#'   \code{\link{rasch_mfrm}}.
#' @param factors A data frame of person factors, a character vector naming
#'   factors nominated in the fit, or a single grouping vector. Defaults to
#'   every factor stored in the fit.
#' @param items Item names or indices to test; all items by default.
#' @param within Names of factors that vary within person (for example
#'   time). Detected automatically when \code{id} is supplied and a factor
#'   varies within an id.
#' @param id Person identifier with one entry per row, or the name of a
#'   nominated factor holding it. By default the identifier stored by the
#'   fitted model is used, so stacked designs retain their pairing.
#' @param contrasts \code{"auto"} (derive the family from the factor
#'   structure) or a named list of numeric cell-weight vectors, each named
#'   by the design-cell labels (factor levels joined by \code{":"}).
#'   Weights are rescaled so the positive and negative parts each sum to
#'   one.
#' @param p_adjust Adjustment across items and contrasts. The default
#'   \code{"holm"} controls familywise error; use \code{"BH"} only for
#'   false-discovery-rate screening. \code{"none"} leaves probabilities
#'   unadjusted.
#' @param alpha Significance level for the adjusted probabilities.
#' @param flag_logits Absolute estimate flagged as practically significant.
#' @param min_n Cells with fewer distinct responders to an item are dropped
#'   from that item's resolution, with a note. When identifiers repeat,
#'   response rows from one person count once within each cell.
#' @return A list of class \code{"rasch_dif_contrasts"}: \code{table} (one row
#'   per item and contrast: estimate in logits, SE, statistic, reference df
#'   (infinite for the normal limit), raw and adjusted p, 95 per cent interval,
#'   \code{significant}, \code{practical}, \code{within}), \code{family}
#'   (the estimable questions with their cell weights), \code{family_n} and
#'   \code{family_n_per_item} (the planned multiplicity counts), the settings,
#'   and any \code{notes}.
#' @references
#' Maxwell, S. E. and Delaney, H. D. (2004). Designing Experiments and
#' Analyzing Data (2nd ed.). Erlbaum.
#'
#' Andrich, D. and Hagquist, C. (2015). Real and artificial differential item
#' functioning in polytomous items. Educational and Psychological Measurement,
#' 75(2), 185--207.
#'
#' Hagquist, C. and Andrich, D. (2017). Recent advances in analysis of
#' differential item functioning in health research using the Rasch model.
#' Health and Quality of Life Outcomes, 15, 181.
#' @seealso \code{\link{dif_anova}} and \code{\link{dif_size}}.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' dif_contrasts(fit, items = c("I3", "I5"))
#' @export
dif_contrasts <- function(fit, factors = NULL, items = NULL, within = NULL,
                          id = NULL, contrasts = "auto", p_adjust = "holm",
                          alpha = 0.05, flag_logits = 0.5, min_n = 20) {
  .check_dif_args(alpha, p_adjust, flag_logits, min_n)
  if (!inherits(fit, "rasch")) stop("dif_contrasts needs a rasch fit")
  if (inherits(fit, "rasch_efrm"))
    .refuse("resolved DIF contrasts are not available for EFRM fits; the ",
            "ordinary split refit would discard the fitted frame units")
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; DIF contrasts are unavailable")
  within <- .check_dif_within(within)
  factors <- .dif_factors(fit, factors)
  .check_dif_factor_levels(factors)
  factors <- as.data.frame(lapply(factors, function(v) {
    f <- droplevels(if (is.ordered(v)) v else factor(v))
    f
  }), check.names = FALSE, stringsAsFactors = FALSE)
  grp <- .factor_cells(factors, sep = ":")
  cellmap <- unique(data.frame(cell = as.character(grp), factors,
                               check.names = FALSE))
  cellmap <- cellmap[match(levels(grp), cellmap$cell), , drop = FALSE]
  usable_id_rows <- .dif_residual_support(fit$residuals)

  # The fitted response-row identifier is the analysis-unit identifier for a
  # repeated-person design. Requiring it again would silently turn an omitted
  # argument into an independent-row analysis.
  if (is.null(id) && !is.null(fit$person$id)) id <- fit$person$id
  if (!is.null(id) && (!is.atomic(id) || !is.null(dim(id))))
    stop("`id` must be an ordinary vector or one fitted factor name",
         call. = FALSE)
  if (is.character(id) && length(id) == 1L && !is.null(fit$factors) &&
      id %in% names(fit$factors)) id <- fit$factors[[id]]
  if (!is.null(id) && length(id) != nrow(factors))
    stop("`id` must have one value per fitted response row")
  if (!is.null(id)) id <- .dif_ids(id)
  paired <- .dif_repeated_support(id, usable_id_rows)
  if (is.null(within) && paired) {
    within <- names(factors)[vapply(names(factors), function(fn)
      .dif_varies_within(factors[[fn]], id, usable_id_rows), TRUE)]
  }
  if (is.null(within)) within <- character(0)
  unknown_w <- setdiff(within, names(factors))
  if (length(unknown_w))
    stop("within-subject factor(s) not among the nominated factors: ",
         paste(unknown_w, collapse = ", "))
  if (length(within) && !paired)
    stop("within-subject factors need repeated person ids (each id ",
         "observed more than once); no id repeats on rows with fitted residuals")
  if (length(within)) {
    varies <- vapply(within, function(fn)
      .dif_varies_within(factors[[fn]], id, usable_id_rows), logical(1))
    if (any(!varies))
      stop("factor(s) declared within-subject never vary within any id: ",
           paste(within[!varies], collapse = ", "))
  }
  if (paired) {
    other <- setdiff(names(factors), within)
    ovaries <- vapply(other, function(fn)
      .dif_varies_within(factors[[fn]], id, usable_id_rows), logical(1))
    if (length(other) && any(ovaries))
      stop("factor(s) vary within persons but are not declared in `within`: ",
           paste(other[ovaries], collapse = ", "),
           "; declare them within-subject so repeated rows are not treated ",
           "as between-person observations")
  }

  family_n_per_item <- NULL
  unavailable_planned <- character(0)
  if (identical(contrasts, "auto")) {
    fam <- .dif_contrast_family(factors, cellmap, within)
    family_n_per_item <- fam$planned_n
    unavailable_planned <- setdiff(fam$planned, names(fam$family))
  } else {
    if (!is.list(contrasts) || is.null(names(contrasts)))
      stop("`contrasts` must be \"auto\" or a named list of cell weights")
    if (anyDuplicated(names(contrasts)))
      stop("duplicate contrast name(s): ",
           paste(unique(names(contrasts)[duplicated(names(contrasts))]),
                 collapse = ", "),
           "; later definitions would be discarded silently")
    if (any(is.na(names(contrasts))) || any(!nzchar(names(contrasts))))
      stop("every contrast needs a non-empty name")
    supplied_meta <- attr(contrasts, "dif_meta", exact = TRUE)
    supplied_family_n <- attr(contrasts, "dif_family_n", exact = TRUE)
    fam <- list(family = list(), meta = list())
    for (nm in names(contrasts)) {
      w <- contrasts[[nm]]
      if (is.null(names(w)) || !all(names(w) %in% cellmap$cell))
        stop("weights of contrast '", nm, "' must be named by design cells: ",
             paste(cellmap$cell, collapse = ", "))
      if (anyDuplicated(names(w)))
        stop("contrast '", nm, "' names cell(s) more than once: ",
             paste(unique(names(w)[duplicated(names(w))]), collapse = ", "),
             "; a repeated cell would silently replace its earlier weight")
      if (!is.numeric(w) || is.complex(w) || !is.null(dim(w)) ||
          !is.null(oldClass(w)) || any(!is.finite(w)))
        stop("weights of contrast '", nm,
             "' must be a plain vector of finite numbers")
      full <- stats::setNames(numeric(nrow(cellmap)), cellmap$cell)
      full[names(w)] <- w
      preserve_scale <- isTRUE(supplied_meta[[nm]]$preserve_scale)
      w <- if (preserve_scale) full else .dif_norm(full)
      if (is.null(w) || !any(w > 0) || !any(w < 0) ||
          abs(sum(w)) > 1e-8)
        stop("contrast '", nm, "' needs positive and negative weights that sum to zero")
      fam$family[[nm]] <- w
      fam$meta[[nm]] <- if (!is.null(supplied_meta[[nm]]))
        supplied_meta[[nm]]
      else list(factors = names(factors), fweights = NULL, within = FALSE)
    }
    family_n_per_item <- length(fam$family)
    if (!is.null(supplied_family_n)) {
      if (!is.numeric(supplied_family_n) || is.complex(supplied_family_n) ||
          !is.null(dim(supplied_family_n)) ||
          !is.null(oldClass(supplied_family_n)) ||
          length(supplied_family_n) != 1L || is.na(supplied_family_n) ||
          !is.finite(supplied_family_n) ||
          supplied_family_n != floor(supplied_family_n) ||
          supplied_family_n < length(fam$family) ||
          supplied_family_n > .Machine$integer.max)
        stop("the internal DIF family size must be one whole number at least as large as the estimable contrast family",
             call. = FALSE)
      family_n_per_item <- as.integer(supplied_family_n)
    }
  }
  if (!length(fam$family)) stop("no contrasts could be formed")

  underlying <- if (inherits(fit, "rasch_mfrm") &&
                    !is.null(fit$virtual_map))
    unique(as.character(fit$virtual_map$item)) else character(0)
  if (is.null(items)) {
    its <- if (length(underlying)) underlying else fit$items$item
  } else {
    if (!is.atomic(items) || !is.null(dim(items)) || !length(items) ||
        anyNA(items))
      stop("`items` must be a non-empty ordinary vector of item names or indices",
           call. = FALSE)
    if (is.factor(items)) items <- as.character(items)
    its <- vapply(items, function(x) {
      if (is.character(x) && length(x) == 1L && x %in% underlying) return(x)
      ii <- .item_idx(fit, x)
      if (length(ii) != 1L || is.na(ii)) return(NA_character_)
      fit$items$item[ii]
    }, character(1))
    if (anyNA(its))
      stop("item(s) not found in the fit: ",
           paste(items[is.na(its)], collapse = ", "))
    if (anyDuplicated(its))
      stop("item(s) named more than once: ",
           paste(unique(its[duplicated(its)]), collapse = ", "),
           "; a repeated item would repeat its hypothesis in the ",
           "multiplicity family")
  }
  Z <- fit$residuals
  notes <- character(0)
  if (length(unavailable_planned))
    notes <- c(notes, paste0(
      "automatic contrast(s) without a complete nuisance-factor stratum were not estimated but remain in the multiplicity family: ",
      paste(unavailable_planned, collapse = ", ")))
  rows <- list()

  for (item in its) {
    # a pooled MFRM item names an underlying item, not a response-cell
    # column; its residual evidence is pooled by .dif_resolve, and the
    # single-column index exists only for ordinary fits
    i <- if (item %in% fit$items$item) .item_idx(fit, item) else NA_integer_
    rs <- .dif_resolve(fit, item, grp, min_n)
    if (!is.null(rs)) notes <- c(notes, rs$notes)
    for (nm in names(fam$family)) {
      w_full <- fam$family[[nm]]
      mt <- fam$meta[[nm]]
      est <- se <- stat <- df <- p <- NA_real_
      if (!is.null(rs)) {
        resolved_cov_ok <- length(dim(rs$vloc)) == 2L &&
          identical(dim(rs$vloc), c(length(rs$levs), length(rs$levs))) &&
          all(is.finite(rs$vloc)) && .covariance_is_symmetric(rs$vloc) &&
          .covariance_is_psd(rs$vloc)
        resolved_reference_ok <- !is.na(rs$df) && rs$df > 0
        w <- w_full[rs$levs]
        w[is.na(w)] <- 0
        # a contrast placing weight on a weakly-identified level rests on a
        # boundary-artefact location: withhold its estimate and SE
        touches_weak <- !is.null(rs$weak) && any(w != 0 & rs$weak)
        preserve_scale <- isTRUE(mt$preserve_scale)
        complete_support <- all(names(w_full)[w_full != 0] %in% rs$levs)
        # A small omitted weight still changes the specified estimand.
        # Require every non-zero cell, for normalised planned contrasts as
        # well as the unscaled post-hoc differences. Never turn a failed
        # five-stratum mean into a four-stratum mean by renormalising it.
        valid_weights <- complete_support &&
          sum(w > 0) > 0 && sum(w < 0) > 0 && abs(sum(w)) < 1e-8
        if (!complete_support)
          notes <- c(notes, paste0(
            item, " [", nm, "]: estimate and inference withheld because ",
            "required contrast cell(s) could not be resolved: ",
            paste(setdiff(names(w_full)[w_full != 0], rs$levs),
                  collapse = ", ")))
        used <- which(w != 0)
        category_ok <- !identical(rs$score_compatible, FALSE) &&
          length(used) >= 2L &&
          all(vapply(seq_len(nrow(rs$category_signature)), function(rr)
            length(unique(rs$category_signature[rr, used])) == 1L, TRUE))
        # When the split refit never ran, the withholding reason is that
        # refusal, already reported once for the whole item; repeating it
        # per contrast as a category mismatch would be false.
        if (valid_weights && !touches_weak && !category_ok &&
            is.null(rs[["refit_error"]])) {
          notes <- c(notes, paste0(
            item, " [", nm, "]: estimate and inference withheld because ",
            "the contrasted cells have different observed response-category ",
            "structures or do not retain the fitted 0:m score structure"))
        }
        if (valid_weights && !touches_weak && category_ok) {
          if (!preserve_scale) w <- .dif_norm(w)
          est <- sum(w * rs$loc)
          if (paired || (resolved_cov_ok && resolved_reference_ok)) {
            if (!paired)
              se <- sqrt(pmax(drop(t(w) %*% rs$vloc %*% w), 0))
          } else {
            reason <- if (resolved_cov_ok && !resolved_reference_ok)
              "does not carry enough independent-person support" else
                "is unavailable or not positive semidefinite"
            notes <- c(notes, paste0(
              item, " [", nm, "]: the resolved-location covariance ",
              reason, "; the estimate is descriptive and Wald inference is withheld"))
          }
        }
      }
      if (!paired) {
        if (is.finite(est)) {
          stat <- .wald_ratio(est, se)
          df <- if (!is.null(rs) && !is.na(rs$df)) rs$df else NA_real_
          p <- 2 * stats::pt(-abs(stat), df = df)
        }
      } else {
        # One person-level calculation covers every repeated design. Using
        # the full resolved-cell weights is essential: a shortcut based on
        # one score per person silently weights nuisance between-person cells
        # by their sample sizes, while the reported resolved estimate averages
        # those cells equally. The test and estimate must address the same
        # marginal contrast.
        attempted_paired <- is.finite(est) && !is.na(i)
        wc <- if (!attempted_paired) NULL else .dif_paired_cell_contrast(
          Z[, i], factors, grp, id, within, cellmap, w_full)
        if (!is.null(wc)) { stat <- -wc$stat; df <- wc$df; p <- wc$p }
        else if (attempted_paired)
          notes <- c(notes, paste0(
            item, " [", nm, "]: person-level residual inference withheld ",
            "because every required between-person stratum must supply at ",
            "least two complete contrasts and their combined sampling ",
            "variance must be positive; ",
            "the resolved logit estimate remains descriptive"))
      }
      if (paired) se <- NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        item = item, contrast = nm, within = isTRUE(mt$within),
        estimate = est, se = se, statistic = stat, df = df, p = p)
    }
  }
  tab <- do.call(rbind, rows)
  family_n <- as.double(family_n_per_item) * length(its)
  if (!is.finite(family_n) || family_n > .Machine$integer.max)
    stop("the planned item-by-contrast family is too large to adjust",
         call. = FALSE)
  tab$p_adj <- .p_adjust_family(tab$p, method = p_adjust, n = family_n)
  critical <- stats::qt(0.975, df = tab$df)
  tab$lower <- tab$estimate - critical * tab$se
  tab$upper <- tab$estimate + critical * tab$se
  tab$significant <- !is.na(tab$p_adj) & tab$p_adj < alpha
  tab$practical <- !is.na(tab$estimate) & abs(tab$estimate) >= flag_logits
  rownames(tab) <- NULL

  fam_df <- data.frame(
    contrast = names(fam$family),
    within = vapply(fam$meta, function(m) isTRUE(m$within), TRUE),
    cells = vapply(fam$family, function(w)
      paste(sprintf("%s %+0.2f", names(w)[w != 0], w[w != 0]),
            collapse = ", "), ""))
  rownames(fam_df) <- NULL

  out <- list(algorithm = "complete-contrast-cells-2",
              table = tab, family = fam_df,
              family_n = as.integer(family_n),
              family_n_per_item = as.integer(family_n_per_item),
              within = within,
              paired = paired, alpha = alpha, p_adjust = p_adjust,
              flag_logits = flag_logits, notes = unique(notes))
  out <- .tag_tables(out)
  class(out) <- "rasch_dif_contrasts"
  out
}

#' Pairwise follow-up comparisons for a DIF term
#'
#' Resolves one item's locations over the complete person-factor design and
#' follows up a selected main effect or interaction. Main effects are pairwise
#' marginal differences. Interactions are differences between those
#' differences, providing a logit-scale magnitude for the interaction itself.
#'
#' @details
#' For levels \eqn{a,b} of one factor, the comparison is
#' \deqn{\Delta_{ba}=\bar\delta_b-\bar\delta_a,}
#' where the bars average equally over complete cells of the other nominated
#' factors. For a two-factor interaction, levels \eqn{a,b} and \eqn{c,d} give
#' \deqn{\Delta_{ba\mathbin{:}dc}=
#' (\delta_{bd}-\delta_{ad})-(\delta_{bc}-\delta_{ac}).}
#' Higher-order interactions use the corresponding tensor-product contrast.
#' Standard errors use the full covariance of the resolved locations.
#'
#' This is the follow-up to a significant DIF term with more than two levels.
#' It reports effects in Rasch logits, adjusts the chosen family of comparisons,
#' and uses person-level scores with the same equal-cell marginal weights in
#' repeated-measures designs. A planned comparison that cannot be estimated
#' because its levels have no common nuisance-factor cell remains in the
#' multiplicity count, although it is omitted from the result table.
#'
#' @param fit A fitted object from \code{\link{rasch}} or
#'   \code{\link{rasch_mfrm}}. EFRM fits are excluded because resolved
#'   comparisons would discard their frame units.
#' @param item Item name or index.
#' @param term A factor name for a main effect, or a character vector of
#'   factor names for an interaction. A single colon-separated string is
#'   also accepted when the factor names themselves contain no colon.
#' @param factors The complete person-factor design, specified as for
#'   \code{\link{dif_contrasts}}. Other factors are retained when calculating
#'   marginal comparisons.
#' @param within Within-person factor names, specified as for
#'   \code{\link{dif_contrasts}}.
#' @param id Person identifiers, specified as for
#'   \code{\link{dif_contrasts}}.
#' @param p_adjust Adjustment over this post-hoc family. The default
#'   \code{"holm"} controls familywise error; use \code{"BH"} only for
#'   false-discovery-rate screening. \code{"none"} leaves probabilities
#'   unadjusted.
#' @param alpha Significance level for adjusted probabilities.
#' @param flag_logits Absolute logit magnitude flagged as practically
#'   important.
#' @param min_n Minimum distinct responders required in a resolved design
#'   cell. When identifiers repeat, response rows from one person count once
#'   within each cell.
#' @return An object of class \code{"rasch_dif_posthoc"}, extending the
#'   \code{\link{dif_contrasts}} result. Its \code{table} contains the pairwise
#'   marginal differences or interaction contrasts, with logit estimates,
#'   standard errors where available, confidence intervals, raw and adjusted
#'   probabilities, and statistical and practical flags.
#' @references Holm, S. (1979). A simple sequentially rejective multiple test
#'   procedure. Scandinavian Journal of Statistics, 6(2), 65--70.
#' @seealso \code{\link{dif_anova}}, \code{\link{dif_size}}, and
#'   \code{\link{dif_contrasts}}.
#' @examples
#' set.seed(1); n <- 800
#' g <- factor(rep(c("A", "B", "C", "D"), each = n / 4))
#' sex <- factor(rep(c("female", "male"), length.out = n))
#' d <- seq(-1.5, 1.5, length.out = 6)
#' sh <- matrix(0, n, 6); sh[g == "D", 2] <- 0.8
#' X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
#' colnames(X) <- paste0("I", 1:6)
#' fit <- rasch(data.frame(X, group = g, sex = sex),
#'              factors = c("group", "sex"))
#' dif_posthoc(fit, "I2", term = "group")
#' @export
dif_posthoc <- function(fit, item, term, factors = NULL, within = NULL,
                        id = NULL, p_adjust = "holm", alpha = 0.05,
                        flag_logits = 0.5, min_n = 20) {
  .check_dif_args(alpha, p_adjust, flag_logits, min_n)
  if (!inherits(fit, "rasch")) stop("dif_posthoc needs a rasch fit")
  if (inherits(fit, "rasch_efrm"))
    .refuse("post-hoc resolved DIF comparisons are not available for EFRM ",
            "fits; the ordinary split refit would discard the fitted frame units")
  within <- .check_dif_within(within)
  if (!is.character(term) || !is.null(dim(term)) || !length(term) ||
      anyNA(term) || any(!nzchar(trimws(term))))
    stop("`term` must contain one or more factor names")
  if (anyDuplicated(term))
    stop("`term` names factor(s) more than once: ",
         paste(unique(term[duplicated(term)]), collapse = ", "))
  if (!is.atomic(item) || !is.null(dim(item)) || length(item) != 1L ||
      is.na(item))
    stop("`item` must name one item; run dif_posthoc() per item so the ",
         "multiplicity adjustment covers one post-hoc family at a time")
  item_name <- is.character(item) && !is.na(item) && nzchar(item)
  item_index <- is.numeric(item) && !is.complex(item) &&
    is.null(dim(item)) && is.null(oldClass(item)) && is.finite(item) &&
    item == floor(item)
  if (!item_name && !item_index)
    stop("`item` must be one non-missing item name or one finite whole-number index")
  item_names <- if (!is.null(fit$items$item)) fit$items$item else colnames(fit$X)
  if (inherits(fit, "rasch_mfrm") && !is.null(fit$virtual_map))
    item_names <- unique(c(as.character(fit$virtual_map$item), item_names))
  ok_item <- if (item_name) item %in% item_names
             else is.finite(item) && item >= 1 &&
               item <= length(fit$items$item %||% colnames(fit$X))
  if (!isTRUE(ok_item))
    stop("item '", item, "' not found in the fit (items: ",
         paste(utils::head(item_names, 8), collapse = ", "),
         if (length(item_names) > 8) ", ..." else "", ")")
  factors <- .dif_factors(fit, factors)
  .check_dif_factor_levels(factors)
  factors <- as.data.frame(lapply(factors, function(v)
    droplevels(if (is.ordered(v)) v else factor(v))),
    check.names = FALSE, stringsAsFactors = FALSE)
  target <- if (all(term %in% names(factors))) term else if (length(term) == 1L)
    .term_vars(term) else term
  fitted_id <- is.null(id)
  if (!is.null(id) && (!is.atomic(id) || !is.null(dim(id))))
    stop("`id` must be an ordinary vector or one fitted factor name",
         call. = FALSE)
  if (is.character(id) && length(id) == 1L && !is.null(fit$factors) &&
      id %in% names(fit$factors)) id <- fit$factors[[id]]
  else if (!is.null(id) && length(id) != nrow(factors))
    # a wrongly sized identifier reaches the within-subject inference and
    # fails there on a base length mismatch, naming neither the argument
    # nor the requirement
    stop("`id` has ", length(id), " entries but the fit has ", nrow(factors),
         " rows; the repeated-measures structure needs one identifier per row")
  if (is.null(id) && !is.null(fit$person$id)) id <- fit$person$id
  if (!is.null(id)) id <- .dif_ids(id)
  usable_id_rows <- .dif_residual_support(fit$residuals)
  paired <- .dif_repeated_support(id, usable_id_rows)
  if (is.null(within) && paired) {
    within <- names(factors)[vapply(names(factors), function(fn)
      .dif_varies_within(factors[[fn]], id, usable_id_rows), TRUE)]
  }
  if (is.null(within)) within <- character(0)
  unknown_within <- setdiff(within, names(factors))
  if (length(unknown_within))
    stop("within-subject factor(s) not found: ",
         paste(unknown_within, collapse = ", "))

  grp <- .factor_cells(factors, sep = ":")
  cellmap <- unique(data.frame(cell = as.character(grp), factors,
                               check.names = FALSE))
  cellmap <- cellmap[match(levels(grp), cellmap$cell), , drop = FALSE]
  fam <- .dif_posthoc_family(factors, cellmap, target, within)
  if (!length(fam$family))
    stop("no complete cells support post-hoc comparisons for term '",
         .dif_term_label(target), "'")
  contrasts <- fam$family
  fam$meta <- lapply(fam$meta, function(x) {
    x$preserve_scale <- TRUE
    x
  })
  attr(contrasts, "dif_meta") <- fam$meta
  attr(contrasts, "dif_family_n") <- fam$planned_n
  out <- dif_contrasts(
    fit, factors = factors, items = item, within = within,
    id = if (fitted_id) NULL else id,
    contrasts = contrasts, p_adjust = p_adjust, alpha = alpha,
    flag_logits = flag_logits, min_n = min_n)
  out$term <- .dif_term_label(target)
  out$type <- if (length(target) > 1L)
    "interaction magnitude" else "pairwise marginal difference"
  unavailable <- setdiff(fam$planned, names(fam$family))
  if (length(unavailable))
    out$notes <- unique(c(out$notes, paste0(
      "post-hoc contrast(s) without a complete nuisance-factor stratum were ",
      "not estimated but remain in the multiplicity family: ",
      paste(unavailable, collapse = ", "))))
  if (nrow(out$table) && all(!is.finite(out$table$estimate)))
    stop("no contrast in the '", out$term, "' family is estimable for item '",
         item, "': a required design cell fell below min_n = ", min_n,
         " responders, the resolved refits were not identified, or the ",
         "contrasted cells used different observed response categories; ",
         "pool sparse levels or check the factor and score coding")
  class(out) <- c("rasch_dif_posthoc", class(out))
  out
}

#' @export
print.rasch_dif_posthoc <- function(x, ...) {
  cat(sprintf("DIF follow-up for %s (%s; %s)\n",
              x$term, x$type, x$p_adjust))
  show <- x$table[, c("item", "contrast", "estimate", "se", "statistic",
                      "p_adj", "significant", "practical")]
  print(.fmt_df(show), row.names = FALSE)
  if (length(x$notes)) cat("\n", paste(x$notes, collapse = "\n"), "\n", sep = "")
  invisible(x)
}

#' @export
print.rasch_dif_contrasts <- function(x, ...) {
  planned <- x$family_n_per_item %||% nrow(x$family)
  question_text <- if (planned == nrow(x$family))
    paste0(planned, " questions") else
      paste0(nrow(x$family), " estimable of ", planned,
             " planned questions")
  cat("Planned DIF contrasts (", question_text, " x ",
      length(unique(x$table$item)), " items; ", x$p_adjust,
      " over the family)\n", sep = "")
  for (r in seq_len(nrow(x$family)))
    cat(sprintf("  %s%s\n", x$family$contrast[r],
                if (x$family$within[r]) "  [within subjects]" else ""))
  if (x$paired)
    cat("Stacked design: tests use person-level residual scores;",
        "logit SEs and intervals are withheld.\n")
  cat("\n")
  tab <- x$table
  show <- tab[, c("item", "contrast", "estimate", "se", "statistic", "p_adj",
                  "significant", "practical")]
  print(.fmt_df(show), row.names = FALSE)
  if (length(x$notes)) cat("\n", paste(x$notes, collapse = "\n"), "\n", sep = "")
  invisible(x)
}
