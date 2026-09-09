# rasch :: testing the invariance a frame model assumes
# ===========================================================================
# rasch_efrm() constrains each item to a single location shared across the
# frames it appears in, scaled by the frame unit. That constraint is the
# model's substantive claim -- the frames differ in unit, not in what the
# items measure -- and it cannot be tested from the fit itself, because the
# fit imposes it. Testing it means calibrating each frame on its own and
# comparing, which is what this function does.
# ===========================================================================

.item_location_covariance <- function(fit) {
  thr <- fit$thresholds
  V <- fit$est$cov_tau
  if (is.null(V) || !nrow(thr) || nrow(V) != nrow(thr)) return(NULL)
  item_names <- fit$items$item
  A <- matrix(0, length(item_names), nrow(thr),
              dimnames = list(item_names, NULL))
  for (i in seq_along(item_names)) {
    rows <- which(thr$item == match(item_names[i], colnames(fit$X)))
    if (length(rows)) A[i, rows] <- 1 / length(rows)
  }
  A %*% V %*% t(A)
}

#' Test item invariance across frames
#'
#' Calibrates each frame separately and compares the locations and
#' discriminations of items administered in more than one frame.
#'
#' @details
#' Let \eqn{\hat\delta_{if}} be the location of item \eqn{i} from a separate
#' calibration of frame \eqn{f}, and let \eqn{\hat\rho_f} be that frame's
#' unit from the fitted EFRM. The common-scale location is
#' \eqn{\hat\delta_{if}^{*}=\hat\delta_{if}/\hat\rho_f}. Because each
#' separate calibration has its own origin, pairwise differences are centred
#' over the common thresholds before testing.
#'
#' The conditional method treats the fitted frame units as fixed.
#' Let \eqn{w_i=m_i/\sum_jm_j}, where \eqn{m_i} is the number of
#' thresholds for item \eqn{i}. If
#' \eqn{C=I-\mathbf{1}\mathbf{w}^{\mathsf T}} centres the common items on the
#' threshold-weighted origin, the covariance of the location differences is
#' \deqn{C\{V_1/\hat\rho_1^2+V_2/\hat\rho_2^2\}C^{\mathsf T}.}
#' This is fast and conditions on the estimated units. The discrimination
#' table gives the difference between the two standardised infit statistics,
#' divided by \eqn{\sqrt{2}}, together with fitted slopes and their ratio.
#' These quantities are descriptive under the conditional method; it does not
#' report discrimination probabilities.
#'
#' With \code{se_method = "bootstrap"}, whole persons are resampled within
#' their observed group, retaining each sampled person's item-set response
#' pattern, and the EFRM and separate frame calibrations are refitted.
#' A replicate is usable only when it retains the observed set of item
#' comparisons, so every centred difference has the same frame origin.
#' Location tests then use the empirical covariance of the centred
#' differences. The discrimination test uses the bootstrap standard error of
#' the log slope ratio. This includes uncertainty in the fitted frame units
#' but is more computationally demanding.
#'
#' Raw and Holm-adjusted probabilities are reported. With conditional
#' uncertainty, Holm adjustment covers the location comparisons. With
#' bootstrap uncertainty, it covers the combined family of location and
#' discrimination comparisons. An unavailable comparison remains in the
#' applicable family. A discrimination probability is unavailable when
#' either separate-frame slope is on its imposed estimation boundary; the
#' ratio remains descriptive and the comparison remains in the Holm family.
#' The summary gives the root mean squared location difference and root mean
#' squared standard error for each set and frame pair. Items from different
#' sets cannot be compared because the sets partition the items.
#' Location differences are relative to the mean difference of the common
#' items. Concentrated DIF can therefore produce non-zero centred contrasts
#' for items that were not themselves shifted. The table identifies the
#' pattern of relative departures; item content or external anchors are needed
#' to determine which items provide the defensible reference.
#' A compared set-by-frame cell must contain at least 50 distinct persons
#' contributing an informative item pair. A pair is informative unless both
#' responses are zero or both are at their item maxima, using the retained
#' items and recoded categories of that frame's separate calibration.
#' Items with weakly determined standard errors in either separate calibration
#' are listed in \code{excluded} rather than tested.
#'
#' A flagged item may be resolved with \code{\link{resolve_frames}} when it
#' remains useful within frames, or removed with \code{\link{drop_items}}
#' when it fits poorly more generally. Either change requires a refit. The
#' invariance tests require a converged frame calibration.
#'
#' @name frame_invariance
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param alpha Significance level used for flags.
#' @param adjust Either \code{"holm"} or \code{"none"}. Both raw and
#'   adjusted probabilities are returned.
#' @param se_method \code{"conditional"} treats the estimated frame units as
#'   fixed; \code{"bootstrap"} refits the complete analysis to whole-person
#'   resamples within group, preserving each person's item-set response
#'   pattern. As in \code{\link{rasch_efrm}}, one response row is required per
#'   person.
#' @param boot_reps Number of bootstrap replicates. At least 30 are required.
#'   At least 90 per cent, and no fewer than 30, must yield the complete set
#'   of comparisons.
#' @param seed Optional bootstrap seed. See \code{\link{rasch_rng}} for
#'   generator support.
#' @return An object of class \code{"rasch_frame_invariance"}. The
#'   \code{locations} and \code{discrimination} tables contain the pairwise
#'   item comparisons; \code{summary} contains set-level RMSD and RMSE
#'   summaries. Under the conditional method, discrimination \code{p},
#'   \code{p_adj}, and \code{flagged} are \code{NA}. \code{excluded} lists
#'   items dropped or rescored by a separate calibration, whose observed
#'   category structures differed between calibrations, or whose
#'   separate-frame estimate was weakly determined.
#'   The remaining components record the multiplicity and uncertainty settings,
#'   including the algorithm identifier, declared comparison-family size
#'   \code{family_n}, and the requested, usable, non-converged and
#'   other-failure bootstrap counts.
#'   \code{bootstrap_stratified} records whether persons were resampled within
#'   group rather than globally.
#' @references
#' Humphry, S. M. (2005). \emph{Maintaining a Common Arbitrary Unit in Social
#' Measurement}. PhD thesis, Murdoch University.
#' @seealso \code{\link{resolve_frames}} to give a flagged item a location
#'   per frame, \code{\link{drop_items}} to remove it altogether, and
#'   \code{\link{rasch_efrm}} for the model whose assumption is tested.
#' @examples
#' d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 1,
#'                    n_groups = 2, group_unit_ratio = 1.4, seed = 2)
#' tr <- attr(d, "truth")
#' fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
#'                   id = "id", boot_reps = 0)
#' frame_invariance(fit)
NULL

.frame_invariance_conditional <- function(fit, strict = TRUE, min_persons = 0L) {
  if (!isTRUE(fit$est$converged)) return(NULL)
  grp <- .frame_group_values(fit)
  glev <- levels(factor(grp))
  sets <- unique(fit$set_of)
  rho <- fit$frames
  vm <- fit$virtual_map
  spec <- fit$refit_spec
  if (is.null(spec)) spec <- list()
  out <- list(); dsc <- list(); excluded <- list(); failures <- character(0)
  failed <- function(set, group, reason) {
    failures <<- c(failures, sprintf("%s/%s (%s)", set, group, reason))
    invisible(NULL)
  }
  for (s in sets) {
    cal <- list()
    observed_groups <- unique(as.character(vm$group[vm$set == s]))
    compared_set <- length(observed_groups) >= 2L
    for (g in glev) {
      vr <- vm$set == s & vm$group == g & vm$vkey %in% colnames(fit$X)
      if (!any(vr)) next
      rows <- !is.na(grp) & as.character(grp) == g
      Xg <- fit$X[rows, vm$vkey[vr], drop = FALSE]
      idg <- fit$person$id[rows]
      colnames(Xg) <- vm$item[vr]
      Xg <- Xg[, colSums(!is.na(Xg)) > 0L, drop = FALSE]
      if (ncol(Xg) < 2L) {
        if (compared_set) failed(s, g, "fewer than two observed items")
        next
      }
      category_signature <- vapply(seq_len(ncol(Xg)), function(j)
        paste(sort(unique(Xg[!is.na(Xg[, j]), j])), collapse = ","), "")
      names(category_signature) <- colnames(Xg)
      fit_error <- NULL
      f <- tryCatch(do.call(rasch, list(
        data = Xg, model = "PCM", id = idg, n_groups = spec$n_groups,
        maxit = spec$maxit %||% 50, tol = spec$tol %||% 1e-7)),
        error = function(e) {
          fit_error <<- conditionMessage(e)
          NULL
        })
      if (is.null(f)) {
        if (compared_set) failed(s, g, paste0("calibration failed: ", fit_error))
        next
      }
      if (!isTRUE(f$est$converged)) {
        if (compared_set) failed(s, g, "calibration did not converge")
        next
      }
      # The separate calibration can drop constant items or recode categories.
      # Its actual score structure, not the pooled EFRM maxima, determines
      # whether a person contributes a non-constant conditional likelihood.
      # Apply the minimum to the observed analysis only: bootstrap replicates
      # must not be selected by whether resampling crosses the support boundary.
      if (compared_set && min_persons > 0L &&
          .frame_person_support(.frame_informative_rows(f$X, f$m),
                                f$person$id) < min_persons) {
        if (!isTRUE(strict)) return(NULL)
        stop("frame-invariance inference needs at least ", min_persons,
             " persons contributing informative item pairs in every compared ",
             "set-by-frame cell; sparse cell(s): ", s, "/", g,
             call. = FALSE)
      }
      # rasch() may legitimately drop a constant item or merge a category
      # that carries no conditional information in this one frame. That item
      # is then unavailable for the frame comparison; it must remain in the
      # declared multiplicity family rather than disappearing because the
      # reduced separate calibration itself happened to converge.
      altered_items <- setdiff(colnames(Xg), colnames(f$X))
      shared_items <- intersect(colnames(Xg), colnames(f$X))
      rescored <- shared_items[!vapply(shared_items, function(nm)
        identical(as.integer(Xg[, nm]), as.integer(f$X[, nm])), logical(1))]
      altered_items <- union(altered_items, rescored)
      r <- rho$rho[rho$set == s & rho$group == g]
      if (length(r) != 1L || !is.finite(r) || r <= 0) {
        if (compared_set) failed(s, g, "fitted frame unit is unavailable")
        next
      }
      V <- .item_location_covariance(f)
      if (!.covariance_supports_wald(V, nrow(f$items))) {
        if (compared_set) failed(s, g, paste(
          "item-location covariance is unavailable, asymmetric, or not",
          "positive semidefinite"))
        next
      }
      good <- is.finite(f$items$location) & is.finite(f$items$se)
      weak_items <- f$items$item[!good]
      if (sum(good) < 2L) {
        if (compared_set) failed(s, g, "fewer than two item locations are estimable")
        next
      }
      V <- V[good, good, drop = FALSE]
      V <- V / r^2
      cal[[g]] <- list(
        table = data.frame(item = f$items$item[good],
          loc = f$items$location[good] / r, se = f$items$se[good] / r,
          n_thresholds = f$m[good],
          category_signature = unname(category_signature[f$items$item[good]]),
          infit = f$items$infit_ms[good], infit_z = f$items$infit_z[good],
          disc = f$items$disc[good], stringsAsFactors = FALSE),
        covariance = V, weak_items = weak_items,
        altered_items = altered_items)
    }
    if (length(failures)) {
      msg <- paste0(
        "frame-invariance inference requires a usable separate calibration ",
        "for every observed frame in a compared item set; unavailable: ",
        paste(unique(failures), collapse = "; "),
        ". Comparisons cannot be selected according to which frame ",
        "calibrations happened to succeed")
      if (isTRUE(strict)) stop(msg, call. = FALSE)
      return(NULL)
    }
    if (length(cal) < 2L) next
    gg <- names(cal)
    for (a in seq_len(length(gg) - 1L)) for (b in (a + 1L):length(gg)) {
      altered_pair <- union(cal[[gg[a]]]$altered_items,
                            cal[[gg[b]]]$altered_items)
      if (length(altered_pair))
        excluded[[length(excluded) + 1L]] <- data.frame(
          set = s, frame_1 = gg[a], frame_2 = gg[b], item = altered_pair,
          reason = paste("category structure changed: dropped or rescored",
                         "in a separate frame calibration"),
          stringsAsFactors = FALSE)
      weak_pair <- setdiff(
        union(cal[[gg[a]]]$weak_items, cal[[gg[b]]]$weak_items),
        altered_pair)
      if (length(weak_pair)) excluded[[length(excluded) + 1L]] <- data.frame(
        set = s, frame_1 = gg[a], frame_2 = gg[b], item = weak_pair,
        reason = "weakly determined in a separate frame calibration",
        stringsAsFactors = FALSE)
      m <- merge(cal[[gg[a]]]$table, cal[[gg[b]]]$table, by = "item",
                 suffixes = c("_1", "_2"))
      m <- m[!m$item %in% altered_pair, , drop = FALSE]
      if (!nrow(m)) next
      comparable <- m$n_thresholds_1 == m$n_thresholds_2 &
        m$category_signature_1 == m$category_signature_2
      if (any(!comparable)) {
        excluded[[length(excluded) + 1L]] <- data.frame(
          set = s, frame_1 = gg[a], frame_2 = gg[b],
          item = m$item[!comparable],
          reason = "different observed category structure",
          stringsAsFactors = FALSE)
        m <- m[comparable, , drop = FALSE]
      }
      if (nrow(m) < 2L) {
        # One location cannot define its own frame-origin contrast. Keep that
        # unavailable comparison in the declared family rather than letting
        # a data-dependent category loss silently reduce multiplicity.
        if (nrow(m)) excluded[[length(excluded) + 1L]] <- data.frame(
          set = s, frame_1 = gg[a], frame_2 = gg[b], item = m$item,
          reason = "fewer than two comparable items to establish the frame origin",
          stringsAsFactors = FALSE)
        next
      }
      raw_d <- m$loc_2 - m$loc_1
      w <- m$n_thresholds_1 / sum(m$n_thresholds_1)
      C <- diag(nrow(m)) - outer(rep(1, nrow(m)), w)
      d <- drop(C %*% raw_d)
      V1 <- cal[[gg[a]]]$covariance
      V2 <- cal[[gg[b]]]$covariance
      V1 <- V1[m$item, m$item, drop = FALSE]
      V2 <- V2[m$item, m$item, drop = FALSE]
      # Var(Cd) = C Sigma C', and C = I - 1w' is not symmetric unless every
      # compared item has the same number of thresholds: without the
      # transpose the polytomous items' standard errors are inflated and the
      # dichotomous items' deflated
      Vd <- C %*% (V1 + V2) %*% t(C)
      se <- sqrt(pmax(diag(Vd), 0))
      loc_wald <- .frame_invariance_wald(d, se)
      z <- loc_wald$statistic
      out[[length(out) + 1L]] <- data.frame(
        set = s, frame_1 = gg[a], frame_2 = gg[b], item = m$item,
        location_1 = m$loc_1, location_2 = m$loc_2,
        difference = d, se = se, statistic = z,
        p = loc_wald$p, stringsAsFactors = FALSE)
      zd <- (m$infit_z_1 - m$infit_z_2) / sqrt(2)
      dsc[[length(dsc) + 1L]] <- data.frame(
        set = s, frame_1 = gg[a], frame_2 = gg[b], item = m$item,
        infit_1 = m$infit_1, infit_2 = m$infit_2,
        infit_z = zd, p = 2 * stats::pnorm(-abs(zd)),
        disc_1 = m$disc_1, disc_2 = m$disc_2,
        disc_ratio = m$disc_2 / m$disc_1,
        disc_boundary = m$disc_1 <= 0.050001 | m$disc_1 >= 4.9999 |
          m$disc_2 <= 0.050001 | m$disc_2 >= 4.9999,
        stringsAsFactors = FALSE)
    }
  }
  excluded <- if (length(excluded)) do.call(rbind, excluded) else
    data.frame(set = character(), frame_1 = character(),
               frame_2 = character(), item = character(),
               reason = character())
  if (!length(out)) {
    if (!nrow(excluded)) return(NULL)
    return(list(
      locations = data.frame(
        set = character(), frame_1 = character(), frame_2 = character(),
        item = character(), location_1 = numeric(), location_2 = numeric(),
        difference = numeric(), se = numeric(), statistic = numeric(),
        p = numeric()),
      discrimination = data.frame(
        set = character(), frame_1 = character(), frame_2 = character(),
        item = character(), infit_1 = numeric(), infit_2 = numeric(),
        infit_z = numeric(), p = numeric(), disc_1 = numeric(),
        disc_2 = numeric(), disc_ratio = numeric(),
        disc_boundary = logical()),
      excluded = excluded))
  }
  list(locations = do.call(rbind, out), discrimination = do.call(rbind, dsc),
       excluded = excluded)
}

# Resample each person-group stratum without relying on sample(x)'s special
# treatment of a length-one numeric vector. If the sole row in a stratum is
# row 17, sample(17, replace = TRUE) samples from 1:17 rather than returning
# row 17; sampling positions within the stratum preserves frame membership.
.frame_stratified_resample <- function(strata) {
  unlist(lapply(strata, function(ii)
    ii[sample.int(length(ii), size = length(ii), replace = TRUE)]),
    use.names = FALSE)
}

# Cluster bootstrap for frame invariance. Repeated rows from one person must
# be sampled together and must keep one bootstrap ID; otherwise both the
# structural EFRM refit and its separate-frame refits treat those rows as
# independent. Where person IDs are nested in groups, stratification retains
# the observed number of persons per group. A person observed in more than
# one group makes that impossible, so sample people globally and retain all
# their group-specific rows.
.frame_cluster_resample <- function(id, group) {
  if (is.null(id) || length(id) != length(group))
    stop("internal frame bootstrap IDs must give one value per response row",
         call. = FALSE)
  cid <- .dif_ids(id)
  grp <- .role_text_values(group)
  clusters <- split(seq_along(cid), cid, drop = TRUE)
  cluster_group <- vapply(clusters, function(ii) {
    z <- unique(grp[ii][!is.na(grp[ii]) & nzchar(grp[ii])])
    if (length(z) == 1L) z else NA_character_
  }, character(1))
  nested <- all(!is.na(cluster_group))
  if (nested) {
    strata <- split(seq_along(clusters), cluster_group, drop = TRUE)
    selected <- unlist(lapply(strata, function(ii)
      ii[sample.int(length(ii), size = length(ii), replace = TRUE)]),
      use.names = FALSE)
  } else {
    selected <- sample.int(length(clusters), size = length(clusters),
                           replace = TRUE)
  }
  blocks <- unname(clusters[selected])
  list(
    rows = unlist(blocks, use.names = FALSE),
    id = rep(sprintf("P%06d", seq_along(blocks)), lengths(blocks)),
    stratified = nested)
}

.frame_ids_cross_groups <- function(id, group) {
  if (is.null(id) || length(id) != length(group)) return(FALSE)
  cid <- .dif_ids(id)
  grp <- .role_text_values(group)
  any(vapply(split(grp, cid, drop = TRUE), function(z) {
    z <- unique(z[!is.na(z) & nzchar(z)])
    length(z) > 1L
  }, logical(1)))
}

# At least one conditional pair has multiple feasible score allocations if
# two scored items are observed and the row is neither all-zero nor all-max.
# Assess this within the current set and frame, not over the person's other
# responses. A one-item response, even in a middle category, supplies no pair.
.frame_informative_rows <- function(X, m) {
  X <- X[, m > 0L, drop = FALSE]
  m <- m[m > 0L]
  rowSums(!is.na(X)) >= 2L &
    rowSums(X, na.rm = TRUE) > 0 &
    rowSums(sweep(X, 2L, m, FUN = "-"), na.rm = TRUE) < 0
}

.frame_person_support <- function(informative, id) {
  if (length(informative) != length(id))
    stop("internal frame support mask has the wrong length", call. = FALSE)
  informative <- !is.na(informative) & informative
  if (!any(informative)) return(0L)
  length(unique(.dif_ids(id)[informative]))
}

# Form a normal-reference statistic only when its estimated uncertainty is
# positive. A constant bootstrap column otherwise turns a non-zero observed
# contrast into Inf and a spurious p = 0.
.frame_invariance_wald <- function(estimate, se) {
  statistic <- .wald_ratio(estimate, se)
  probability <- rep(NA_real_, length(statistic))
  usable <- is.finite(statistic)
  probability[usable] <- 2 * stats::pnorm(-abs(statistic[usable]))
  list(statistic = statistic, p = probability)
}

# A bootstrap contrast must use the same set of items that established the
# observed frame origin. Merely matching the observed rows is not enough: a
# category absent from the observed calibration can reappear in a resample,
# adding an item to the bootstrap centring set and changing the estimand.
.frame_invariance_boot_vector <- function(table, key, value,
                                          transform = identity) {
  if (!is.data.frame(table) || !all(c(
    "set", "frame_1", "frame_2", "item", value) %in% names(table)))
    return(NULL)
  got <- .factor_keys(
    table[, c("set", "frame_1", "frame_2", "item"), drop = FALSE])
  if (length(got) != length(key) || anyDuplicated(got) ||
      !setequal(got, key)) return(NULL)
  out <- transform(table[[value]][match(key, got)])
  if (!is.numeric(out) || length(out) != length(key) ||
      any(!is.finite(out))) return(NULL)
  unname(out)
}

.frame_invariance_probabilities <- function(cmp, dsc, excluded, se_method,
                                            alpha, adjust) {
  if (identical(se_method, "bootstrap")) {
    boundary <- dsc$disc_boundary %in% TRUE
    dsc$statistic[boundary] <- NA_real_
    dsc$p[boundary] <- NA_real_
  }
  allp <- c(cmp$p, dsc$p)
  padj <- rep(NA_real_, length(allp))
  usable <- is.finite(allp)
  n_excluded <- nrow(excluded)
  family_n <- nrow(cmp) + n_excluded + if (se_method == "bootstrap")
    nrow(dsc) + n_excluded else 0L
  padj[usable] <- stats::p.adjust(allp[usable], method = "holm", n = family_n)
  cmp$p_adj <- padj[seq_len(nrow(cmp))]
  dsc$p_adj <- padj[nrow(cmp) + seq_len(nrow(dsc))]
  p_cmp <- if (adjust == "holm") cmp$p_adj else cmp$p
  cmp$flagged <- ifelse(is.finite(p_cmp), p_cmp < alpha, NA)
  p_dsc <- if (adjust == "holm") dsc$p_adj else dsc$p
  dsc$flagged <- ifelse(is.finite(p_dsc), p_dsc < alpha, NA)
  list(locations = cmp, discrimination = dsc, family_n = family_n)
}

#' @rdname frame_invariance
#' @export
frame_invariance <- function(fit, alpha = 0.05, adjust = c("holm", "none"),
                             se_method = c("conditional", "bootstrap"),
                             boot_reps = 200, seed = NULL) {
  if (!inherits(fit, "rasch_efrm"))
    stop("frame_invariance needs a fit from rasch_efrm()")
  if (!is.null(seed)) seed <- .check_whole(seed, "seed", 0)
  if (!isTRUE(fit$est$converged))
    stop("the frame calibration did not converge; invariance tests are unavailable")
  .check_prob(alpha, "alpha")
  adjust <- match.arg(adjust)
  se_method <- match.arg(se_method)
  grp <- .frame_group_values(fit)
  glev <- levels(factor(grp))
  if (length(glev) < 2L)
    stop("item invariance across frames needs at least two person groups: ",
         "with one group each item appears in a single frame, and item sets ",
         "partition the items, so use the item fit statistics within each ",
         "set instead")
  if (identical(se_method, "conditional") &&
      .frame_ids_cross_groups(fit$person$id, grp))
    .refuse("conditional frame-invariance inference is unavailable because ",
            "at least one person appears in more than one frame and the ",
            "cross-frame covariance between separate calibrations is unknown; ",
            "use se_method = \"bootstrap\" to resample whole persons")

  # A separate-frame calibration supplies the covariance used by every item
  # comparison. Small frames produced valid-looking but unstable normal tests
  # in simulation, so require 50 informative persons in every observed
  # set-by-frame cell before reporting invariance probabilities. This pooled
  # screen is preliminary; the separate calibrations below check support
  # again after any item removal or category recoding.
  vm <- fit$virtual_map
  fr <- unique(vm[, c("set", "group")])
  n_frames_by_set <- table(fr$set)
  sparse <- vapply(seq_len(nrow(fr)), function(i) {
    if (n_frames_by_set[fr$set[i]] < 2L) return(FALSE)
    cc <- which(vm$set == fr$set[i] & vm$group == fr$group[i] &
                  vm$vkey %in% colnames(fit$X))
    if (length(cc) < 2L) return(FALSE)
    jj <- match(vm$vkey[cc], colnames(fit$X))
    informative <- .frame_informative_rows(
      fit$X[, jj, drop = FALSE], fit$m[jj])
    .frame_person_support(informative, fit$person$id) < 50L
  }, logical(1))
  if (any(sparse)) stop(
    "frame-invariance inference needs at least 50 persons contributing ",
    "informative item pairs in every compared set-by-frame cell; sparse cell(s): ",
    paste(paste(fr$set[sparse], fr$group[sparse], sep = "/"), collapse = ", "))

  ans <- .frame_invariance_conditional(fit, min_persons = 50L)
  if (is.null(ans))
    stop("no item set is taken by two person groups, so no item appears in ",
         "two frames to be compared")
  cmp <- ans$locations
  dsc <- ans$discrimination
  if (!nrow(cmp))
    .refuse("no frame pair retains at least two items with the same observed ",
            "category structure; frame-invariance comparisons are unavailable")
  reps_used <- 0L
  reps_nonconverged <- 0L
  reps_errors <- 0L
  minimum_usable <- 0L
  if (se_method == "bootstrap") {
    boot_reps <- .check_whole(boot_reps, "boot_reps", 30)
    if (!is.null(seed)) {
      old_seed <- .sim_seed_capture()
      on.exit(.sim_seed_restore(old_seed), add = TRUE)
      set.seed(seed)
    }
    source <- .efrm_source_matrix(fit)
    key4 <- function(x) .factor_keys(
      x[, c("set", "frame_1", "frame_2", "item"), drop = FALSE])
    lk <- key4(cmp)
    dk <- key4(dsc)
    bd <- matrix(NA_real_, boot_reps, nrow(cmp))
    ba <- matrix(NA_real_, boot_reps, nrow(dsc))
    status <- rep("error", boot_reps)
    for (b in seq_len(boot_reps)) {
      sample_b <- .frame_cluster_resample(fit$person$id, grp)
      ii <- sample_b$rows
      fb <- tryCatch(.efrm_refit(
        fit, source[ii, , drop = FALSE], fit$set_of, boot_reps = 0,
        ids = paste0("B", sprintf("%04d", b), sample_b$id),
        factors = fit$factors[ii, , drop = FALSE], se_method = "hybrid"),
        error = function(e) NULL)
      if (is.null(fb)) next
      if (!isTRUE(fb$est$converged) || !.efrm_link_converged(fb)) {
        status[b] <- "nonconverged"
        next
      }
      ib <- .frame_invariance_conditional(fb, strict = FALSE)
      if (is.null(ib)) next
      loc_b <- .frame_invariance_boot_vector(
        ib$locations, lk, "difference")
      disc_b <- .frame_invariance_boot_vector(
        ib$discrimination, dk, "disc_ratio", log)
      if (is.null(loc_b) || is.null(disc_b)) next
      bd[b, ] <- loc_b
      ba[b, ] <- disc_b
      status[b] <- "used"
    }
    good <- rowSums(is.finite(bd)) == ncol(bd) &
      rowSums(is.finite(ba)) == ncol(ba)
    status[status == "used" & !good] <- "error"
    reps_used <- sum(good)
    reps_nonconverged <- sum(status == "nonconverged")
    reps_errors <- sum(status == "error")
    minimum_usable <- .fit_min_boot_success(boot_reps)
    if (reps_used < minimum_usable)
      stop("only ", reps_used, " of ", boot_reps,
           " frame-invariance bootstrap refits succeeded (",
           reps_nonconverged, " did not converge; ", reps_errors,
           " otherwise failed); at least ", minimum_usable,
           " are required")
    cmp$se <- apply(bd[good, , drop = FALSE], 2, stats::sd)
    loc_wald <- .frame_invariance_wald(cmp$difference, cmp$se)
    cmp$statistic <- loc_wald$statistic
    cmp$p <- loc_wald$p
    dsc$log_disc_ratio <- log(dsc$disc_ratio)
    dsc$se_log_disc_ratio <- apply(ba[good, , drop = FALSE], 2, stats::sd)
    disc_wald <- .frame_invariance_wald(
      dsc$log_disc_ratio, dsc$se_log_disc_ratio)
    dsc$statistic <- disc_wald$statistic
    dsc$p <- disc_wald$p
  } else {
    # The standardised-infit comparison was anti-conservative in the
    # validation study (7.1% combined Holm FWER over 2,000 null replicates,
    # with strong item-position dependence). Retain its descriptive value,
    # but do not attach an inferential probability to it.
    dsc$p <- NA_real_
  }
  # Conditional inference covers locations only. The validated bootstrap
  # adds discrimination and controls the two tables as one family.
  inference <- .frame_invariance_probabilities(
    cmp, dsc, ans$excluded, se_method, alpha, adjust)
  cmp <- inference$locations
  dsc <- inference$discrimination
  rownames(cmp) <- NULL
  rownames(dsc) <- NULL

  frame_pairs <- unique(cmp[c("set", "frame_1", "frame_2")])
  smry <- lapply(seq_len(nrow(frame_pairs)), function(i) {
    k <- frame_pairs[i, ]
    same <- function(d) d$set == k$set & d$frame_1 == k$frame_1 &
      d$frame_2 == k$frame_2
    z <- cmp[same(cmp), ]; y <- dsc[same(dsc), ]
    nx <- sum(same(ans$excluded))
    rmsd <- sqrt(mean(z$difference^2))
    rmse <- sqrt(mean(z$se^2))
    ratio <- if (is.finite(rmse) &&
                  rmse > sqrt(.Machine$double.eps) * pmax(1, rmsd))
      rmsd / rmse else NA_real_
    data.frame(set = z$set[1], frame_1 = z$frame_1[1], frame_2 = z$frame_2[1],
      n_items = nrow(z), n_excluded = nx,
      rmsd = rmsd, rmse = rmse, ratio = ratio,
      n_location = sum(z$flagged %in% TRUE),
      n_discrimination = if (identical(se_method, "bootstrap"))
        sum(y$flagged %in% TRUE) else NA_integer_,
      stringsAsFactors = FALSE)
  })
  smry <- if (length(smry)) do.call(rbind, smry) else data.frame(
    set = character(), frame_1 = character(), frame_2 = character(),
    n_items = integer(), n_excluded = integer(), rmsd = numeric(),
    rmse = numeric(), ratio = numeric(), n_location = integer(),
    n_discrimination = integer())
  rownames(smry) <- NULL
  out <- .tag_tables(list(locations = cmp, discrimination = dsc,
                          summary = smry, excluded = ans$excluded,
                          algorithm = "frame-invariance-complete-family-1",
                          alpha = alpha, adjust = adjust,
                          se_method = se_method,
                          family_n = inference$family_n,
                          boot_reps = if (se_method == "bootstrap")
                            as.integer(boot_reps) else 0L,
                          boot_reps_used = reps_used,
                          boot_reps_nonconverged = reps_nonconverged,
                          boot_reps_errors = reps_errors,
                          boot_minimum_usable = minimum_usable,
                          bootstrap_stratified = if (se_method == "bootstrap")
                            !.frame_ids_cross_groups(fit$person$id, grp) else
                            NA,
                          seed = seed,
                          fit_signature = .fit_boot_signature(fit)))
  out$result_signature <- .fit_boot_md5(out)
  class(out) <- c("rasch_frame_invariance", "list")
  out
}

.validate_frame_invariance <- function(invariance, fit) {
  if (is.null(invariance)) return(invisible(NULL))
  if (!inherits(fit, "rasch_efrm"))
    stop("`invariance` is only available for an EFRM fit")
  signature <- invariance$result_signature
  unsigned <- unclass(invariance)
  unsigned$result_signature <- NULL
  if (!inherits(invariance, "rasch_frame_invariance") ||
      !identical(invariance$algorithm,
                 "frame-invariance-complete-family-1") ||
      !is.data.frame(invariance$summary) ||
      !is.data.frame(invariance$locations) ||
      !is.data.frame(invariance$discrimination) ||
      !is.character(signature) || length(signature) != 1L || is.na(signature) ||
      !.fit_boot_hash_matches(signature, unsigned) ||
      !.fit_boot_signature_matches(invariance$fit_signature, fit))
    stop("`invariance` must be a frame_invariance() result from this fitted model")
  counts <- c("boot_reps", "boot_reps_used", "boot_reps_nonconverged",
              "boot_reps_errors", "boot_minimum_usable")
  whole <- function(x)
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0 &&
      x == floor(x)
  if (!all(counts %in% names(invariance)) ||
      any(!vapply(invariance[counts], whole, logical(1))) ||
      invariance$boot_reps_used + invariance$boot_reps_nonconverged +
        invariance$boot_reps_errors != invariance$boot_reps ||
      (identical(invariance$se_method, "bootstrap") &&
       (invariance$boot_minimum_usable !=
          .fit_min_boot_success(invariance$boot_reps) ||
        invariance$boot_reps_used < invariance$boot_minimum_usable)) ||
      (!identical(invariance$se_method, "bootstrap") &&
       any(unlist(invariance[counts]) != 0)))
    stop("`invariance` has inconsistent bootstrap accounting; recompute it ",
         "with frame_invariance()")
  stratified_ok <- "bootstrap_stratified" %in% names(invariance) &&
    is.logical(invariance$bootstrap_stratified) &&
    length(invariance$bootstrap_stratified) == 1L &&
    if (identical(invariance$se_method, "bootstrap"))
      !is.na(invariance$bootstrap_stratified) else
      is.na(invariance$bootstrap_stratified)
  if (!stratified_ok)
    stop("`invariance` has inconsistent bootstrap resampling provenance; ",
         "recompute it with frame_invariance()")
  expected_family <- nrow(invariance$locations) + nrow(invariance$excluded) +
    if (identical(invariance$se_method, "bootstrap"))
      nrow(invariance$discrimination) + nrow(invariance$excluded) else 0L
  if (!is.numeric(invariance$family_n) || length(invariance$family_n) != 1L ||
      !is.finite(invariance$family_n) || invariance$family_n != expected_family)
    stop("`invariance` has inconsistent multiplicity accounting; recompute it ",
         "with frame_invariance()")
  invisible(invariance)
}

#' @export
print.rasch_frame_invariance <- function(x, ...) {
  adj <- if (is.null(x$adjust)) "holm" else x$adjust
  pcol <- if (adj == "holm") "p_adj" else "p"
  rule <- if (adj == "holm") "Holm-adjusted" else "unadjusted, screening"
  cat("Item invariance across frames (each frame calibrated separately)\n\n")
  cat("Uncertainty:", if (identical(x$se_method, "bootstrap"))
    sprintf(paste0(if (isTRUE(x$bootstrap_stratified))
                     "whole-person bootstrap within person group " else
                     "whole-person bootstrap across persons ",
                   "(%d/%d usable; ",
                   "%d non-converged; %d other failures)"),
            x$boot_reps_used, x$boot_reps, x$boot_reps_nonconverged,
            x$boot_reps_errors) else "conditional on the fitted frame units",
    "\n\n")
  print(.fmt_df(x$summary), row.names = FALSE)
  if (!is.null(x$excluded) && nrow(x$excluded)) {
    why <- table(x$excluded$reason)
    cat(sprintf("\n%d item comparison(s) excluded:\n", nrow(x$excluded)))
    for (i in seq_along(why))
      cat(sprintf("  %d %s\n", unname(why[i]), names(why)[i]))
  }
  cat("\nrmsd/rmse above 1 indicates item behaviour the frame units do not account for\n")
  fl <- x$locations[x$locations$flagged %in% TRUE, , drop = FALSE]
  if (nrow(fl)) {
    cat(sprintf(paste0("\nLocation differs across frames in %d item ",
                       "comparison(s) (%d item(s)) at alpha = %.2f (%s):\n"),
                nrow(fl), length(unique(fl$item)), x$alpha, rule))
    print(.fmt_df(fl[, c("set", "frame_1", "frame_2", "item", "difference",
                         "se", "statistic", pcol)]), row.names = FALSE)
  } else {
    cat(sprintf("\nNo available item-location comparison differs across frames at alpha = %.2f (%s).\n",
                x$alpha, rule))
  }
  if (!identical(x$se_method, "bootstrap")) {
    cat("\nThe discrimination comparisons are descriptive:\n")
    cols <- c("set", "frame_1", "frame_2", "item", "infit_1", "infit_2",
              "infit_z", "disc_1", "disc_2", "disc_ratio", "disc_boundary")
    print(.fmt_df(x$discrimination[, intersect(cols,
                                                names(x$discrimination)),
                                           drop = FALSE]), row.names = FALSE)
    cat("Use se_method = \"bootstrap\" for discrimination probabilities.\n")
    return(invisible(x))
  }
  nb <- sum(x$discrimination$disc_boundary %in% TRUE)
  if (nb)
    cat(sprintf("\n%d discrimination probability/probabilities withheld because a fitted slope is on its boundary.\n",
                nb))
  fd <- x$discrimination[x$discrimination$flagged %in% TRUE, , drop = FALSE]
  if (nrow(fd)) {
    cat(sprintf(paste0("\nDiscrimination differs across frames in %d item ",
                       "comparison(s) (%d item(s)):\n"),
                nrow(fd), length(unique(fd$item))))
    cols <- c("set", "frame_1", "frame_2", "item", "log_disc_ratio",
              "se_log_disc_ratio", "statistic", pcol, "disc_1", "disc_2",
              "disc_ratio", "disc_boundary")
    print(.fmt_df(fd[, intersect(cols, names(fd)), drop = FALSE]),
          row.names = FALSE)
  } else {
    cat("\nNo available discrimination comparison differs across frames.\n")
  }
  invisible(x)
}
