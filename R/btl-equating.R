# rasch :: common-object equating for paired comparisons
# ===========================================================================
# The paired-comparison analogue of equate_tests(). Two Bradley-Terry-Luce
# calibrations that share a set of common objects -- the same scripts,
# performances, or products judged by different panels, or by the same panel
# in different years -- each fix their own origin by the sum-zero constraint.
# Because the two constraints are imposed over DIFFERENT object sets, the two
# origins do not coincide even when the objects are unchanged: each scale is
# centred on the mean of a different collection. Equating therefore normally
# estimates the scale shift between the origins (the precision-weighted mean
# difference over the common objects) and then tests each common object against
# the shifted identity line. Its variance includes both the covariance induced
# by each sum-zero calibration and estimation of the shift from those same
# objects; it is not the naive sum of two marginal variances. When external
# anchors have already fixed a shared origin, shift = "none" instead tests the
# raw differences and uses their marginal variances.
#
# Objects that survive define the equating link and carry the second panel's
# whole scale onto the first; objects that fail show drift -- a script the two
# panels valued differently, or a standard that moved between years -- and
# should be reviewed before the link is trusted. This is the standards-
# maintenance use of comparative judgement (Bramley 2007): a common set of
# anchor scripts lets panels judged apart be placed on one scale.
# ===========================================================================

# Coerce the second calibration to (object, location, se). It may be another
# btl fit or a "bank" -- a data frame of previously banked object locations,
# the paired-comparison counterpart of an item bank.
.btl_equate_ref <- function(reference) {
  if (inherits(reference, "rasch_btl")) {
    tab <- reference$objects
    # an extrapolated boundary location is a reporting value, not a
    # calibrated estimate; it takes no part in equating
    if ("extreme" %in% names(tab)) tab <- tab[!tab$extreme, ]
    return(data.frame(object = as.character(tab$object),
                      location = tab$location,
                      se = tab$se,
                      stringsAsFactors = FALSE))
  }
  if (is.data.frame(reference) || (is.list(reference) && !is.null(names(reference)))) {
    .check_column_names(reference)
    reference <- as.data.frame(reference, stringsAsFactors = FALSE)
    .check_column_names(reference)
    if (!all(c("object", "location") %in% names(reference)))
      stop("a bank needs columns 'object' and 'location' (and ideally 'se')")
    if (!"se" %in% names(reference)) reference$se <- NA_real_
    out <- data.frame(object = .role_text_values(reference$object),
                      location = .bank_numeric(reference$location, "location"),
                      se = .bank_numeric(reference$se, "se"),
                      stringsAsFactors = FALSE)
    if (anyNA(out$object) || any(!nzchar(out$object)))
      stop("bank object names must be non-missing and non-empty")
    if (anyDuplicated(out$object))
      stop("bank object names must be unique: ",
           paste(unique(out$object[duplicated(out$object)]), collapse = ", "))
    if (any(!is.finite(out$location)))
      stop("bank object locations must be finite")
    if (any(!is.na(out$se) & (!is.finite(out$se) | out$se < 0)))
      stop("bank standard errors must be non-negative finite values or NA")
    return(out)
  }
  stop("`fit2` must be a btl fit or a bank data frame (object, location, se)")
}

.btl_equate_bank_cov <- function(reference, ids) {
  C <- attr(reference, "cov_location", exact = TRUE)
  if (is.null(C)) return(NULL)
  if (!is.matrix(C) || !is.numeric(C) || any(!is.finite(C)) ||
      !identical(dim(C), c(length(ids), length(ids))))
    stop("attr(fit2, 'cov_location') must be a finite numeric square matrix ",
         "with one row and column per bank object")
  if (!.covariance_is_symmetric(C))
    stop("attr(fit2, 'cov_location') must be symmetric")
  if (!is.null(rownames(C)) || !is.null(colnames(C))) {
    if (is.null(rownames(C)) || is.null(colnames(C)) ||
        anyNA(match(ids, rownames(C))) || anyNA(match(ids, colnames(C))))
      stop("named bank covariance rows and columns must match every bank object")
    C <- C[ids, ids, drop = FALSE]
  }
  if (!.covariance_is_psd(C))
    stop("attr(fit2, 'cov_location') must be positive semidefinite")
  C
}

# Residual degrees of freedom carried by a calibration covariance. A clustered
# BTL fit has judges as its independent sampling units. An external bank may
# supply the corresponding value explicitly; otherwise its covariance is
# treated as asymptotically normal.
.btl_equate_cov_df <- function(x, objects = NULL) {
  # These linking errors condition on estimated stage-one quantities. A
  # fixed origin removes shift uncertainty, not the omitted frame-parameter
  # uncertainty. They therefore cannot become full sampling errors merely
  # because shift = "none" needs only marginal variances.
  if (inherits(x, "rasch_btl_efrm") &&
      identical(x$se_method, "conditional")) return(NA_real_)
  if (inherits(x, "rasch_btl") && !inherits(x, "rasch_btl_efrm")) {
    requested <- as.character(objects)
    fixed <- names(x$anchors)
    # A named external anchor is a constant, not an estimate. Its variance
    # and covariance are exactly zero even when the free sandwich is
    # unsupported, so an all-anchor contrast has the limiting (exact) df.
    if (length(requested) && length(fixed) &&
        all(!is.na(requested) & requested %in% fixed))
      return(Inf)
    cl <- x[["cl"]]
    if (!is.null(cl) && identical(cl$inference_available, FALSE))
      return(NA_real_)
  }
  # Ordinary BTL covariance is judge-clustered when a judge role was fitted.
  # A BTL-EFRM fit always records judges, but only the judge bootstrap treats
  # them as the sampling units. Its parametric bootstrap is comparison-level,
  # so applying J - 1 degrees of freedom to it would mix two uncertainty
  # schemes and make the equating tests needlessly conservative.
  judge_based <- if (inherits(x, "rasch_btl_efrm"))
    identical(x$se_method, "judge_bootstrap") else
      inherits(x, "rasch_btl") && isTRUE(x$clustered)
  if (judge_based && !is.null(x$comparisons$judge) &&
      any(!is.na(x$comparisons$judge))) {
    df <- max(length(unique(
      x$comparisons$judge[!is.na(x$comparisons$judge)]
    )) - 1L, 1L)
    if (inherits(x, "rasch_btl_efrm") &&
        identical(x$se_method, "judge_bootstrap")) {
      # Common-scale locations outside the reference set inherit alpha and
      # kappa uncertainty along that set's path to the reference. They cannot
      # use more denominator information in equating than the fit allowed for
      # those unit parameters. Panel support also enters every stage-one fit.
      if (!is.null(x$phi_table$df)) {
        panel_df <- unique(x$phi_table$df)
        if (anyNA(panel_df)) return(NA_real_)
        df <- min(df, panel_df)
      }
      if (!is.null(objects) && length(objects) &&
          all(c("object", "set") %in% names(x$objects)) &&
          all(c("set", "df") %in% names(x$alpha_table))) {
        object_set <- x$objects$set[match(objects, x$objects$object)]
        reference <- x$alpha_table$set[1L]
        linked <- unique(object_set[!is.na(object_set) &
                                      object_set != reference])
        if (length(linked)) {
          set_df <- x$alpha_table$df[match(linked, x$alpha_table$set)]
          if (anyNA(set_df)) return(NA_real_)
          df <- min(df, set_df)
        }
      }
    }
    return(as.numeric(df))
  }
  z <- attr(x, "df_location", exact = TRUE)
  if (is.null(z)) return(Inf)
  if (!is.numeric(z) || is.complex(z) || length(z) != 1L ||
      !is.null(dim(z)) || !is.null(oldClass(z)) || is.na(z) ||
      !is.finite(z) || z <= 0)
    stop("attr(fit2, 'df_location') must be one positive numeric degree ",
         "of freedom", call. = FALSE)
  as.numeric(z)
}

#' Equate two paired-comparison calibrations through their common objects
#'
#' Places two Bradley--Terry--Luce calibrations on a common origin using their
#' shared objects, then tests the shared objects for drift. The second
#' calibration may be a fitted model or an object bank.
#'
#' @details
#' Let \eqn{d_j} be the location difference for common object \eqn{j} and
#' \eqn{v_j} its marginal variance. With \code{shift = "mean"}, the origin
#' shift is the precision-weighted mean
#' \deqn{\hat s=\frac{\sum_j d_j/v_j}{\sum_j 1/v_j}.}
#' If fewer than two common objects have usable variances but at least two have
#' finite locations, their unweighted mean difference is returned as a
#' descriptive fallback and recorded in \code{shift_method}.
#' Every error built from those precision weights conditions on them as if
#' the two calibrations' standard errors were known: \code{shift_se}, each
#' drift contrast's \code{se_diff} (and so its probability and
#' \code{drifting} flag), and the equated location errors that carry the
#' shift. Weight uncertainty adds a positive term all of them omit, so they
#' understate uncertainty -- intervals under-cover, drift probabilities run
#' small -- when the calibrations rest on few judges; a judge resample of
#' both calibrations is the weight-aware alternative.
#' An exact common anchor determines the shift even when it is the only
#' common object with usable uncertainty.
#' Each object is tested using its shifted difference \eqn{d_j-\hat s}. The
#' covariance calculation retains the dependence induced by the sum-zero
#' constraints. Drift tests then require independent calibrations and at least
#' three common objects with usable, positive-semidefinite joint covariance
#' information. Two common objects identify a descriptive origin shift, but
#' do not support an object-drift test. With \code{shift = "none"}, the origin
#' is fixed before the comparison and each object's variance is the sum of its
#' two marginal variances; joint covariance information and a three-object
#' link are unnecessary. One common object is sufficient for that fixed-origin
#' comparison; estimating a shift still requires at least two.
#' A judge-clustered ordinary BTL covariance, or a BTL--EFRM covariance from
#' the judge bootstrap, uses finite judge-cluster degrees of freedom. A
#' contrast involving only fixed external anchors has exact zero covariance
#' from that calibration and therefore uses infinite degrees of freedom even
#' when inference for its estimated objects is unavailable. A
#' BTL--EFRM location outside the reference set is also limited by the
#' weakest edge on its strongest supported path to that reference. A
#' comparison-level parametric-bootstrap BTL--EFRM covariance uses the
#' asymptotic normal reference instead. Conditional frame errors are
#' preliminary and do not support drift inference, including comparisons on
#' a fixed origin.
#' Binary fits have no threshold parameters, so their recorded threshold
#' structure does not affect compatibility. Polytomous fits must use the
#' same category scale and threshold structure.
#'
#' The \code{equated} table includes uncertainty in the estimated shift.
#' For independent calibrations, with \eqn{y_j=b_j+\hat s},
#' \deqn{\operatorname{Var}(y_j)=\operatorname{Var}(b_j)+
#' \operatorname{Var}(\hat s)+2\operatorname{Cov}(b_j,\hat s).}
#' These SEs are withheld if joint uncertainty is unavailable. A fixed shift
#' (\code{shift = "none"} or an exact common anchor) leaves supported original
#' SEs unchanged. SEs from a conditional frame reference are withheld in the
#' equated bank. When available, the table carries its full covariance in
#' \code{attr(equated, "cov_location")} and conservative finite sampling-unit
#' degrees of freedom in \code{attr(equated, "df_location")}. An equated bank
#' is not independent of either calibration used to construct it.
#'
#' The common-object set should contain a stable majority. If most common
#' objects move in the same direction, the estimated shift follows them and
#' stable objects can appear to drift. In that case, repeat the equating with a
#' substantively justified anchor set.
#'
#' @param fit1 A fitted object from \code{\link{btl}}: the calibration whose
#'   scale (origin) the equating targets.
#' @param fit2 A second \code{\link{btl}} fit, or a bank: a data frame with
#'   columns \code{object}, \code{location}, and optionally \code{se}; object
#'   names and column names must be unique. Numeric fields may be numeric
#'   columns, numeric text, or factors with numeric labels; other column
#'   classes are refused. Locations must be finite. Bank-based drift inference
#'   with an estimated mean shift requires the joint location covariance as a square matrix in
#'   \code{attr(fit2, "cov_location")}, ordered like the bank rows (or named by
#'   object), unless the bank is treated as fixed with zero SEs. Marginal
#'   standard errors are sufficient with \code{shift = "none"}. A bank whose
#'   covariance was estimated from a finite number of independent sampling
#'   units may carry their residual degrees of freedom in
#'   \code{attr(fit2, "df_location")} as one positive numeric value. For a polytomous
#'   fit the bank must carry
#'   \code{attr(bank, "m")} matching the number of fitted score steps.
#' @param alpha Significance level for the (multiplicity-adjusted) drift tests.
#' @param p_adjust Adjustment for the common-object tests, passed to
#'   \code{stats::p.adjust}. The default is \code{"holm"}. A common object
#'   remains in the family when its drift probability is unavailable.
#' @param independent Whether the calibrations have independent judges and
#'   comparisons. For two fitted objects the default \code{NULL} withholds
#'   drift tests until independence is stated explicitly. Bank tables are
#'   treated as independent unless \code{FALSE} is supplied. Dependent
#'   calibrations require a joint or paired bootstrap for inference.
#' @param shift \code{"mean"} (default) estimates the origin shift from the
#'   common objects; \code{"none"} compares raw locations when both
#'   calibrations have already been placed on the same externally anchored
#'   scale.
#' @return A list of class \code{"rasch_btl_equate"}: the comparison
#'   \code{table} (per common object: object, both locations and standard
#'   errors, their \code{difference}, the \code{shifted_difference} against the
#'   estimated origin, the pooled \code{se_diff}, \code{t}, raw and adjusted
#'   \code{p}, and the \code{drifting} flag); the estimated \code{shift}, its
#'   \code{shift_method} and \code{shift_se}; \code{equated}, the second
#'   calibration's full object table re-expressed on \code{fit1}'s scale; the
#'   number of common objects
#'   \code{n_common}; the number usable for inference \code{n_inference};
#'   whether inference was available \code{inferential}; \code{alpha};
#'   \code{p_adjust}; the requested \code{shift_setting}; and \code{notes}.
#'   A drift probability is withheld when its contrast has zero estimated
#'   uncertainty. Such an object remains in the multiplicity family.
#' @references Bramley, T. (2007). Paired comparison methods. In P. Newton,
#'   J. Baird, H. Goldstein, H. Patrick, & P. Tymms (Eds.), \emph{Techniques
#'   for monitoring the comparability of examination standards} (pp. 246-294).
#'   London: Qualifications and Curriculum Authority.
#' @examples
#' set.seed(1)
#' beta <- setNames(seq(-2, 2, length.out = 8), paste0("O", 1:8))
#' sim <- function(objs) {
#'   pr <- t(utils::combn(objs, 2))
#'   d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
#'   d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#'   btl(d, "a", "b", "win")
#' }
#' eq <- btl_equate(sim(paste0("O", 1:7)), sim(paste0("O", 2:8)),
#'                   independent = TRUE)
#' eq$table
#' @export
btl_equate <- function(fit1, fit2, alpha = 0.05, p_adjust = "holm",
                       independent = NULL, shift = c("mean", "none")) {
  shift <- match.arg(shift)
  .check_prob(alpha, "alpha")
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      !is.null(dim(p_adjust)) || !is.null(oldClass(p_adjust)) ||
      !p_adjust %in% stats::p.adjust.methods)
    stop("p_adjust must name a method in stats::p.adjust.methods")
  if (inherits(fit1, "rasch_btl_explanatory") ||
      inherits(fit2, "rasch_btl_explanatory"))
    stop("object drift is not defined for an explanatory comparison fit: ",
         "the object locations are functions of their predictors. Equate ",
         "the unrestricted calibrations")
  if (!is.null(independent) && (length(independent) != 1L ||
      !is.logical(independent) || !is.null(dim(independent)) ||
      !is.null(oldClass(independent)) || is.na(independent)))
    stop("independent must be NULL, TRUE, or FALSE")
  if (!inherits(fit1, "rasch_btl"))
    stop("`fit1` must be a paired-comparison (btl) fit")
  if (inherits(fit2, "rasch") && !inherits(fit2, "rasch_btl"))
    stop("`fit2` must be a paired-comparison fit or object bank; a ",
         "person-by-item Rasch calibration is on a different response scale")
  # equating a non-converged calibration carries its unidentified locations
  # and understated covariance straight into the drift table: refuse it
  for (nm in c("fit1", "fit2")) {
    f <- get(nm)
    if (inherits(f, "rasch_btl") && !isTRUE(f$converged))
      stop("`", nm, "` did not converge (its comparison design does not ",
           "identify some object locations); resolve that before equating -- ",
           "the drift statistics would inherit boundary estimates and ",
           "understated standard errors", call. = FALSE)
  }
  if (inherits(fit2, "rasch_btl")) {
    # Binary fits have no threshold parameters. Ordinary BTL records "free",
    # whereas a frame fit records
    # "dichotomous"; those labels do not change the binary response model.
    if (!identical(fit1$m, fit2$m) ||
        !identical(as.character(fit1$categories),
                   as.character(fit2$categories)) ||
        (fit1$m > 1L && !identical(fit1$thr_structure, fit2$thr_structure)))
      stop("the paired-comparison calibrations use incompatible response ",
           "scales or threshold structures; shift-only equating requires ",
           "the same model and category scale")
  } else if (fit1$m > 1L) {
    fm <- attr(fit2, "m", exact = TRUE)
    if (!is.numeric(fm) || is.complex(fm) || length(fm) != 1L ||
        !is.null(dim(fm)) || !is.null(oldClass(fm)) || is.na(fm) ||
        !is.finite(fm) || fm != floor(fm) || fm > .Machine$integer.max ||
        !identical(as.integer(fm), as.integer(fit1$m)))
      stop("a polytomous-comparison bank must carry attr(bank, 'm') matching ",
           "the fitted number of score steps; otherwise scale compatibility ",
           "cannot be established")
  }
  cur_tab <- fit1$objects
  if ("extreme" %in% names(cur_tab)) cur_tab <- cur_tab[!cur_tab$extreme, ]
  cur <- data.frame(object = as.character(cur_tab$object),
                    location = cur_tab$location,
                    se = cur_tab$se, stringsAsFactors = FALSE)
  ref <- .btl_equate_ref(fit2)
  # If bank metadata are supplied, validate them even when another design
  # condition later withholds inference. Silently ignoring a malformed df
  # would let the same bank switch reference distribution across analyses.
  if (!inherits(fit2, "rasch_btl")) invisible(.btl_equate_cov_df(fit2))
  common <- intersect(cur$object, ref$object)
  min_common <- if (identical(shift, "mean")) 2L else 1L
  if (length(common) < min_common)
    stop(if (identical(shift, "mean"))
      paste("need at least two common objects to estimate an origin shift",
            "between paired-comparison scales") else
      paste("need at least one common object for a fixed-origin",
            "paired-comparison comparison"))
  a <- cur[match(common, cur$object), ]
  b <- ref[match(common, ref$object), ]
  bank_cov <- if (inherits(fit2, "rasch_btl")) NULL else
    .btl_equate_bank_cov(fit2, ref$object)
  if (!is.null(bank_cov)) {
    cov_se <- sqrt(pmax(diag(bank_cov), 0))
    stated <- is.finite(ref$se)
    if (any(stated & !.se_covariance_agree(ref$se, cov_se)))
      stop("the bank standard errors must agree with the diagonal of ",
           "attr(fit2, 'cov_location')")
    ref$se[!stated] <- cov_se[!stated]
    b <- ref[match(common, ref$object), ]
  }
  d <- a$location - b$location
  v <- a$se^2 + b$se^2
  finite_loc <- is.finite(d)
  exact1 <- !inherits(fit1, "rasch_btl_efrm") &
    common %in% names(fit1$anchors)
  exact2 <- if (inherits(fit2, "rasch_btl"))
    !inherits(fit2, "rasch_btl_efrm") & common %in% names(fit2$anchors)
  else is.finite(b$se) & b$se == 0
  exact_pair <- exact1 & exact2
  # A zero sandwich diagonal is not an exact calibration. It can enter the
  # link with zero pooled variance only when both object locations were
  # explicitly fixed; otherwise it is excluded from weighting and testing.
  usable <- is.finite(d) & is.finite(v) & (v > 0 | exact_pair)
  zero_not_exact <- is.finite(d) & is.finite(v) & v == 0 & !exact_pair
  .check_exact_equating_origin(d, exact_pair,
                               estimates_shift = identical(shift, "mean"),
                               label = "object")
  independent_ok <- if (is.null(independent)) !inherits(fit2, "rasch_btl")
                    else isTRUE(independent)
  estimates_shift <- identical(shift, "mean")
  joint_cov_1 <- !estimates_shift || !is.null(fit1$cov_beta) ||
    all(a$se[usable] == 0)
  joint_cov_2 <- if (!estimates_shift) TRUE else if (inherits(fit2, "rasch_btl"))
    !is.null(fit2$cov_beta) || all(b$se[usable] == 0)
  else
    !is.null(bank_cov) || all(b$se[usable] == 0)
  joint_cov_ok <- joint_cov_1 && joint_cov_2
  min_inference <- if (estimates_shift) 3L else 1L
  df1 <- .btl_equate_cov_df(fit1, common[usable])
  df2 <- .btl_equate_cov_df(fit2, common[usable])
  conditional1 <- inherits(fit1, "rasch_btl_efrm") &&
    identical(fit1$se_method, "conditional")
  conditional2 <- inherits(fit2, "rasch_btl_efrm") &&
    identical(fit2$se_method, "conditional")
  df_available <- !is.na(df1) && !is.na(df2)
  inferential <- independent_ok && sum(usable) >= min_inference &&
    joint_cov_ok && df_available
  # Two stable objects identify an origin shift; three are required only to
  # distinguish individual drift from that estimated shift. Do not let the
  # inferential threshold replace the documented link estimator.
  if (!estimates_shift) {
    w <- numeric(0)
    c0 <- 0
    shift_method <- "none"
  } else if (sum(usable) >= 2L || any(exact_pair & usable)) {
    w <- .inverse_variance_weights(v[usable])
    # precision-weighted mean difference: the shift between the two sum-zero
    # origins, best estimated where both calibrations are most certain
    c0 <- sum(w * d[usable]) / sum(w)
    shift_method <- "precision-weighted"
  } else {
    if (sum(finite_loc) < 2L)
      stop("need at least two common objects with finite locations to estimate ",
           "an origin shift")
    w <- numeric(0)
    c0 <- mean(d[finite_loc])
    shift_method <- "unweighted"
  }
  # the common objects' location estimates are CORRELATED within each
  # sum-zero calibration, so Var(c0) = u' (Sigma1 + Sigma2) u with
  # u = w / sum(w), taken from the stored sandwich covariances. A bank
  # contributes an explicitly attached joint covariance; without it,
  # inference is withheld unless the bank is fixed (all SEs zero).
  covsub <- function(fit, objs_c) {
    if (inherits(fit, "rasch_btl") && !is.null(fit$cov_beta)) {
      if (!is.null(rownames(fit$cov_beta)))
        return(fit$cov_beta[objs_c, objs_c, drop = FALSE])
      i <- match(objs_c, as.character(fit$objects$object))
      fit$cov_beta[i, i, drop = FALSE]
    } else {
      C <- .btl_equate_bank_cov(fit, ref$object)
      if (is.null(C)) NULL else {
        ii <- match(objs_c, ref$object)
        C[ii, ii, drop = FALSE]
      }
    }
  }
  exact_link <- estimates_shift && any(exact_pair & usable)
  shift_se <- if (!estimates_shift || exact_link) 0 else NA_real_
  se_diff <- t <- df <- p <- p_adj <-
    rep(NA_real_, length(common)); drifting <- rep(NA, length(common))
  testable <- rep(FALSE, length(common))
  covariance_invalid <- FALSE
  if (inferential && !estimates_shift) {
    se_diff[usable] <- sqrt(pmax(v[usable], 0))
    t[usable] <- .wald_ratio(d[usable], se_diff[usable])
    v1 <- a$se[usable]^2
    v2 <- b$se[usable]^2
    den <- if (is.finite(df1)) v1^2 / df1 else rep(0, length(v1))
    if (is.finite(df2)) den <- den + v2^2 / df2
    dfs <- ifelse(den > 0, (v1 + v2)^2 / den, Inf)
    df[usable] <- dfs
    p[usable] <- 2 * stats::pt(-abs(t[usable]), df = dfs)
    testable <- usable & is.finite(p)
    if (any(testable)) p_adj[testable] <- p.adjust(
      p[testable], method = p_adjust, n = length(common))
    drifting[testable] <- p_adj[testable] < alpha
  } else if (inferential) {
    u <- w / sum(w)
    S1all <- tryCatch(covsub(fit1, common), error = function(e) NULL)
    S1 <- if (is.null(S1all)) diag(a$se[usable]^2, sum(usable))
          else S1all[usable, usable, drop = FALSE]
    S2all <- tryCatch(covsub(fit2, common), error = function(e) NULL)
    S2 <- if (is.null(S2all)) diag(b$se[usable]^2, sum(usable))
          else S2all[usable, usable, drop = FALSE]
    covariance_extract_failed <-
      (!is.null(fit1$cov_beta) && is.null(S1all)) ||
      ((inherits(fit2, "rasch_btl") && !is.null(fit2$cov_beta) ||
        !inherits(fit2, "rasch_btl") && !is.null(bank_cov)) &&
       is.null(S2all))
    if (covariance_extract_failed ||
        !.covariance_supports_wald(S1, sum(usable)) ||
        !.covariance_supports_wald(S2, sum(usable))) {
      inferential <- FALSE
      covariance_invalid <- TRUE
    } else {
      shift_se <- sqrt(pmax(drop(t(u) %*% (S1 + S2) %*% u), 0))
      # each drift test compares d_i - c0, and c0 is estimated from the SAME
      # common objects: Var(d_i - c0) = [(I - 1u') Sigma (I - u 1')]_ii
      #   = Sigma_ii - 2 (Sigma u)_i + u' Sigma u,
      # not the naive Sigma_ii -- ignoring the estimated shift (and the
      # within-calibration covariance) mis-states every drift p-value
      Sg <- S1 + S2
      Su <- drop(Sg %*% u)
      var_d <- pmax(diag(Sg) - 2 * Su + drop(t(u) %*% Su), 0)
      se_diff[usable] <- sqrt(var_d)
      t[usable] <- .wald_ratio(d[usable] - c0, se_diff[usable])
      # Welch-Satterthwaite reference for independent fitted calibrations.
      # Each shifted contrast h = e_i - u receives the finite-judge
      # contribution from each panel separately; a non-clustered fit or
      # fixed/external bank has infinite df and contributes no denominator.
      Hc <- diag(length(u)) - matrix(u, nrow = length(u), ncol = length(u),
                                    byrow = TRUE)
      v1 <- pmax(diag(Hc %*% S1 %*% t(Hc)), 0)
      v2 <- pmax(diag(Hc %*% S2 %*% t(Hc)), 0)
      den <- if (is.finite(df1)) v1^2 / df1 else rep(0, length(v1))
      if (is.finite(df2)) den <- den + v2^2 / df2
      dfs <- ifelse(den > 0, (v1 + v2)^2 / den, Inf)
      df[usable] <- dfs
      p[usable] <- 2 * stats::pt(-abs(t[usable]), df = dfs)
      testable <- usable & is.finite(p)
      if (any(testable)) p_adj[testable] <- p.adjust(
        p[testable], method = p_adjust, n = length(common))
      drifting[testable] <- p_adj[testable] < alpha
    }
  }
  zero_uncertainty <- inferential & usable & !is.finite(p)
  tab <- data.frame(object = common,
                    location_1 = a$location, se_1 = a$se,
                    location_2 = b$location, se_2 = b$se,
                    difference = d, shifted_difference = d - c0,
                    se_diff = se_diff, t = t, df = df,
                    p = p, p_adj = p_adj,
                    drifting = drifting, stringsAsFactors = FALSE)
  rownames(tab) <- NULL
  # Re-expression changes uncertainty as well as location. Under independent
  # calibrations, Cov(b_j, shift) = -Cov(b_j, b_common) u. A marginal-SE-only
  # bank cannot supply this cross term, even for its non-common objects.
  known_shift <- !estimates_shift || exact_link
  reference_df <- .btl_equate_cov_df(fit2, ref$object)
  reference_cov <- tryCatch(covsub(fit2, ref$object), error = function(e) NULL)
  if (is.null(reference_cov) && all(is.finite(ref$se) & ref$se == 0) &&
      (!inherits(fit2, "rasch_btl") || all(ref$object %in% names(fit2$anchors))))
    reference_cov <- matrix(0, nrow(ref), nrow(ref))
  full_cov_ok <- !is.na(reference_df) && all(is.finite(ref$se)) &&
    .covariance_supports_wald(reference_cov, nrow(ref))
  equated_cov <- NULL
  equated_se <- if (known_shift && !conditional2) ref$se else
    rep(NA_real_, nrow(ref))
  if (known_shift && full_cov_ok) {
    equated_cov <- reference_cov
  } else if (!known_shift && independent_ok && is.finite(shift_se) &&
             full_cov_ok && length(w)) {
    u <- w / sum(w)
    k <- drop(reference_cov[, match(common[usable], ref$object), drop = FALSE] %*% u)
    candidate <- reference_cov + shift_se^2 - outer(k, rep(1, length(k))) -
      outer(rep(1, length(k)), k)
    candidate <- (candidate + t(candidate)) / 2
    if (.covariance_supports_wald(candidate, nrow(ref))) {
      equated_cov <- candidate
      equated_se <- sqrt(pmax(diag(candidate), 0))
    }
  }
  equated <- data.frame(object = ref$object,
                        location = ref$location + c0,
                        se = equated_se, stringsAsFactors = FALSE)
  rownames(equated) <- NULL
  if (!is.null(equated_cov)) {
    dimnames(equated_cov) <- list(ref$object, ref$object)
    attr(equated, "cov_location") <- equated_cov
  }
  equated_df <- if (known_shift) reference_df else min(df1, reference_df)
  if (is.finite(equated_df) && any(is.finite(equated_se)))
    attr(equated, "df_location") <- equated_df
  if (fit1$m > 1L) attr(equated, "m") <- fit1$m
  notes <- if (estimates_shift) sprintf(paste0(
    "Origins differ because each calibration is sum-zero over its own object ",
    "set; a shift of %.3f logits aligns fit2 to fit1."), c0) else
      paste("The calibrations were compared on their fixed common origin;",
            "no common-object shift was estimated.")
  if (exact_link)
    notes <- c(notes, "Exact common anchors determine the origin shift.")
  # u is built from the calibrations' own estimated standard errors, so every
  # quadratic form in u conditions on weights that are themselves estimated:
  # the shift error, each drift contrast error, and the equated errors that
  # carry the shift. The omitted term is positive, and it grows as the
  # standard errors become noisier, so name every affected quantity that this
  # call actually reports and state the direction, rather than manufacture a
  # correction.
  if (shift_method == "precision-weighted" && !exact_link &&
      is.finite(shift_se)) {
    weighted <- c("the shift standard error",
                  if (any(testable)) "each drift contrast standard error",
                  if (any(is.finite(equated_se)))
                    "the equated location standard errors")
    notes <- c(notes, paste0(
      "Every error built from the estimated precision weights (",
      paste(weighted, collapse = "; "),
      ") treats the two calibrations' standard errors as known. The omitted ",
      "weight-uncertainty term is positive, and it grows as those standard ",
      "errors become noisier: with few judges per panel these errors ",
      "understate the uncertainty and their intervals under-cover",
      if (any(testable)) ", and the drift probabilities run small" else "",
      ". Resample the judges of both calibrations for weight-aware errors."))
  }
  if (conditional1 || conditional2)
    notes <- c(notes, paste(
      "Drift tests are withheld because conditional frame errors do not",
      "propagate all estimation uncertainty; use bootstrap errors for inference."))
  if (conditional2)
    notes <- c(notes, paste(
      "Equated location standard errors are withheld because the reference",
      "has preliminary conditional frame errors. Original conditional SEs",
      "remain in the input calibration."))
  if (!conditional2 && !known_shift && anyNA(equated_se))
    notes <- c(notes, paste(
      "Equated location standard errors are withheld because the estimated",
      "shift and reference locations lack supported joint uncertainty.",
      "Original-scale standard errors remain in the input calibration."))
  if (any(zero_not_exact))
    notes <- c(notes, paste0(
      "Objects with zero pooled estimated variance but without exact anchors ",
      "were excluded from weighting and drift tests: ",
      paste(common[zero_not_exact], collapse = ", "), "."))
  if (covariance_invalid)
    notes <- c(notes, paste(
      "Drift tests are withheld because a fitted object-location covariance",
      "is unavailable or not positive semidefinite."))
  if (!df_available && !(conditional1 || conditional2))
    notes <- c(notes, paste(
      "Drift tests are withheld because a contributing panel or set link",
      "does not have enough independent-judge support for inference."))
  if (shift_method == "unweighted") {
    notes <- c(notes, paste(
      "Fewer than two common objects had usable variances; the reported",
      "shift is the unweighted mean of the finite location differences and is descriptive."))
    no_se <- finite_loc & !usable
    if (any(no_se)) notes <- c(notes, paste0(
        "Objects included in the descriptive shift but excluded from drift tests ",
        "because their uncertainty was unavailable or was a non-exact zero: ",
        paste(common[no_se], collapse = ", "), "."))
  } else if (any(!usable)) {
    notes <- c(notes, paste0(if (estimates_shift)
      "Objects excluded from the precision-weighted shift and drift tests " else
      "Drift tests unavailable for common objects ",
      "because their locations or uncertainty were unavailable or unusable: ",
      paste(common[!usable], collapse = ", "), "."))
  }
  if (any(!finite_loc))
    notes <- c(notes, paste0(
      "Objects excluded from the shift and drift tests because their locations ",
      "were unavailable: ", paste(common[!finite_loc], collapse = ", "), "."))
  if (any(zero_uncertainty))
    notes <- c(notes, paste0(
      "Drift probability/probabilities withheld for zero contrast ",
      "uncertainty: ",
      paste(common[zero_uncertainty], collapse = ", "), "."))
  if (any(testable & is.finite(df)))
    notes <- c(notes, paste(
      "Drift probabilities use contrast-specific Welch-Satterthwaite",
      "degrees of freedom for the finite judge-cluster covariances."))
  if (is.null(independent) && inherits(fit2, "rasch_btl"))
    notes <- c(notes, paste(
      "Drift tests withheld because independence between fitted",
      "calibrations was not stated; set independent = TRUE only for",
      "independent judges and comparisons."))
  if (identical(independent, FALSE))
    notes <- c(notes, paste(
      "Drift tests withheld for dependent calibrations because cross-fit",
      "covariance is unavailable; use a joint or paired bootstrap."))
  if (independent_ok && sum(usable) >= min_inference && !joint_cov_ok)
    notes <- c(notes, paste(
      "Drift tests withheld because every calibration with non-zero",
      "marginal SEs needs its joint object-location covariance; marginal",
      "SEs do not carry the calibration-origin covariance. For a bank,",
      "supply this in attr(fit2, 'cov_location'). Frame-dependent fits",
      "supply it when bootstrap standard errors are used."))
  if (independent_ok && sum(usable) < min_inference)
    notes <- c(notes, if (estimates_shift) paste(
      "Drift tests withheld because at least three common objects with",
      "standard errors are required.") else paste(
        "Drift tests withheld because no common object has locations and",
        "standard errors in both calibrations."))
  if (any(drifting %in% TRUE))
    notes <- c(notes, sprintf(
      "%d common object(s) drift beyond the shifted link: %s",
      sum(drifting %in% TRUE), paste(common[drifting %in% TRUE], collapse = ", ")))
  inferential <- inferential && any(testable)
  structure(class = "rasch_btl_equate",
            list(table = tab, shift = c0, shift_method = shift_method,
                 shift_se = shift_se,
                 equated = equated, n_common = length(common),
                 n_inference = sum(testable), inferential = inferential,
                 alpha = alpha, p_adjust = p_adjust, shift_setting = shift,
                 notes = notes))
}

#' Plot a paired-comparison equating comparison
#'
#' Scatter of the two calibrations' common-object locations with the shifted
#' identity line, per-object contrast intervals at the requested confidence
#' level, and a dotted guide band at their average half-width; objects that drift (after
#' the multiplicity adjustment) are highlighted and labelled. The counterpart
#' of \code{\link{plot_equate}} for Bradley-Terry-Luce scales.
#'
#' @param fit1 A fitted object from \code{\link{btl}}.
#' @param fit2 A second \code{\link{btl}} fit, or a bank data frame with columns
#'   \code{object}, \code{location}, and optionally \code{se}.
#' @param ... Passed to \code{\link{btl_equate}} (e.g. \code{alpha},
#'   \code{p_adjust}).
#' @return Called for its plotting side effect; invisibly the
#'   \code{\link{btl_equate}} result.
#' @examples
#' set.seed(1)
#' beta <- setNames(seq(-2, 2, length.out = 8), paste0("O", 1:8))
#' sim <- function(objs) {
#'   pr <- t(utils::combn(objs, 2))
#'   d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
#'   d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#'   btl(d, "a", "b", "win")
#' }
#' plot_btl_equate(sim(paste0("O", 1:7)), sim(paste0("O", 2:8)),
#'                  independent = TRUE)
#' @export
plot_btl_equate <- function(fit1, fit2, ...) {
  eq <- btl_equate(fit1, fit2, ...)
  tab <- eq$table
  paired <- is.finite(tab$location_1) & is.finite(tab$location_2)
  if (!any(paired))
    .refuse("no common object has finite locations in both calibrations; ",
            "there is nothing to display")
  cor_paired <- if (sum(paired) >= 2L &&
                    stats::sd(tab$location_1[paired]) > 0 &&
                    stats::sd(tab$location_2[paired]) > 0)
    stats::cor(tab$location_1[paired], tab$location_2[paired]) else NA_real_
  rng <- range(c(tab$location_1[paired], tab$location_2[paired])) +
    c(-0.4, 0.4)
  op <- .rr_canvas(rng, rng, "Calibration 2 location (logits)",
                   "Calibration 1 location (logits)",
                   sprintf("%d common objects, shift %.3f, r = %.3f",
                           eq$n_common, eq$shift, cor_paired),
                   grid_x = TRUE)
  on.exit(par(op))
  abline(eq$shift, 1, col = .rr$ink, lwd = 2)
  half <- .equate_interval_halfwidth(tab, alpha = eq$alpha)
  band_rows <- paired & is.finite(half)
  band <- if (any(band_rows)) mean(half[band_rows]) else NA_real_
  if (is.finite(band)) {
    abline(eq$shift + band, 1, lty = 3, col = .rr$soft)
    abline(eq$shift - band, 1, lty = 3, col = .rr$soft)
  }
  hs <- paired & is.finite(half)
  segments(tab$location_2[hs], tab$location_1[hs] - half[hs],
           tab$location_2[hs], tab$location_1[hs] + half[hs],
           col = paste0(.rr$soft, "88"))
  points(tab$location_2[paired], tab$location_1[paired], pch = 21, cex = 1.6,
         bg = ifelse(tab$drifting[paired] %in% TRUE, .rr$red, .rr$blue),
         col = "white", lwd = 1.2)
  dr <- paired & tab$drifting %in% TRUE
  if (any(dr))
    text(tab$location_2[dr], tab$location_1[dr],
         tab$object[dr], pos = 3, offset = 0.5, cex = 0.75,
         col = .rr$red)
  invisible(eq)
}

#' @export
print.rasch_btl_equate <- function(x, ...) {
  tab <- x$table
  paired <- is.finite(tab$location_1) & is.finite(tab$location_2)
  cor_paired <- if (sum(paired) >= 2L &&
                    stats::sd(tab$location_1[paired]) > 0 &&
                    stats::sd(tab$location_2[paired]) > 0)
    stats::cor(tab$location_1[paired], tab$location_2[paired]) else NA_real_
  method <- if (is.null(x$shift_method)) "method unavailable" else x$shift_method
  cat(sprintf(paste0("Common-object equating over %d object(s): shift %.3f ",
                     "(%s; SE %s), correlation %s, RMSD %s\n"),
              x$n_common, x$shift, method,
              if (is.finite(x$shift_se)) sprintf("%.3f", x$shift_se) else "withheld",
              if (is.finite(cor_paired)) sprintf("%.3f", cor_paired)
              else "unavailable",
              if (any(is.finite(tab$shifted_difference))) sprintf("%.3f",
                sqrt(mean(tab$shifted_difference^2, na.rm = TRUE)))
              else "unavailable"))
  core <- c("object", "location_1", "location_2", "shifted_difference", "t",
            "p_adj", "drifting")
  print(.fmt_df(tab[, intersect(core, names(tab))]), row.names = FALSE)
  if (isTRUE(x$inferential))
    cat(sprintf("%d object(s) drift beyond the %s-adjusted %.0f%% level.\n",
                sum(tab$drifting %in% TRUE), x$p_adjust, 100 * (1 - x$alpha)))
  else cat("Drift inference withheld; see $notes.\n")
  cat("(standard errors and unadjusted columns on $table; fit2 on fit1's scale in $equated)\n")
  invisible(x)
}
