# rasch :: test equating
# ===========================================================================
# Comparison of common-item locations across two separately analysed
# datasets. Because each analysis fixes its own origin, the comparison
# allows for a scale shift (estimated by precision-weighted mean difference)
# and, when the calibrations are independent, tests each common item against
# the shifted identity line. The variance includes estimation of the shift
# from the same correlated item locations; it is not the naive sum of two
# marginal variances.
#
# Items that survive define the equating link; items that fail show item
# drift and should be dropped from the link (or anchored individually).
# ===========================================================================

# Item-location covariance over a set of items: block means of the fit's
# threshold covariance (recentred parameterisation), which carries the
# negative correlations induced by the identification constraint. A bank
# table contributes its attached joint covariance. Marginal standard errors
# alone cannot recover the correlations needed when the origin shift is
# estimated from the common items. They are sufficient for shift = "none",
# where every item is tested on an origin fixed before this comparison.
.equate_loc_cov <- function(obj, items) {
  if (inherits(obj, "rasch")) {
    # A fitted calibration needs its joint threshold covariance here. Do not
    # silently reconstruct a diagonal matrix from marginal SEs: that loses
    # the centring covariance induced by the fitted origin and understates or
    # overstates the variance after estimating the common-item shift.
    if (is.null(obj$est$cov_tau)) return(NULL)
    thr <- obj$est$thr
    # thr$item holds integer item positions in column order, which is the
    # order of obj$items: the match below is by position, not name
    idx <- match(items, obj$items$item)
    rows <- lapply(idx, function(i) thr$id[thr$item == i])
    S <- matrix(0, length(items), length(items))
    for (i in seq_along(rows)) for (j in seq_along(rows))
      S[i, j] <- mean(obj$est$cov_tau[rows[[i]], rows[[j]]])
    # A location is an exact constant only when every threshold contributing
    # to it is anchored. This includes a location anchor, whose individual
    # threshold spread remains estimated. Its row and column are nevertheless
    # exactly zero on the item-location scale, even when an unsupported
    # clustered covariance makes those individual threshold entries NA.
    exact <- vapply(idx, function(i) {
      z <- thr$anchored[thr$item == i]
      length(z) && all(z)
    }, logical(1))
    if (any(exact)) {
      S[exact, ] <- 0
      S[, exact] <- 0
    }
    return(S)
  }
  ref <- .equate_ref(obj)
  supplied <- .equate_bank_cov(obj, ref$item)
  if (!is.null(supplied)) {
    ii <- match(items, ref$item)
    return(supplied[ii, ii, drop = FALSE])
  }
  se <- ref$se[match(items, ref$item)]
  diag(se^2, length(items))
}

# Residual degrees of freedom carried by an item-location covariance. A
# repeated-person calibration has person clusters as its independent sampling
# units. Its cluster covariance is usable only when the calibration guard
# accepted it. An external bank can declare finite sampling-unit degrees of
# freedom; without that metadata its covariance is treated as asymptotically
# normal.
.equate_cov_df <- function(x) {
  if (inherits(x, "rasch")) {
    support <- x$est$cluster_support
    # This support object is derived from informative conditional item-pair
    # contributions. Raw IDs are not authoritative: repeated all-missing rows
    # must not turn an otherwise row-independent covariance into a finite-df
    # cluster covariance.
    if (is.list(support)) {
      if (!isTRUE(x$est$cluster_inference)) return(NA_real_)
      if (!(identical(support$repeated, FALSE) ||
            identical(support$repeated, TRUE)))
        return(NA_real_)
      if (identical(support$repeated, FALSE)) return(Inf)
    } else {
      repeated <- isTRUE(x$repeated_ids) ||
        (!is.null(x$person$id) && .has_repeated_person_ids(x$person$id))
      return(if (repeated) NA_real_ else Inf)
    }
    if (!is.numeric(support$n) || length(support$n) != 1L ||
        !is.finite(support$n) || support$n < 2)
      return(NA_real_)
    return(as.numeric(support$n - 1))
  }
  z <- attr(x, "df_location", exact = TRUE)
  if (is.null(z)) return(Inf)
  if (!is.numeric(z) || is.complex(z) || length(z) != 1L ||
      !is.null(dim(z)) || !is.null(oldClass(z)) || is.na(z) ||
      !is.finite(z) || z <= 0)
    stop("attr(reference, 'df_location') must be one positive numeric degree ",
         "of freedom", call. = FALSE)
  as.numeric(z)
}

.equate_exact_locations <- function(x, items) {
  if (!inherits(x, "rasch")) return(rep(FALSE, length(items)))
  idx <- match(items, x$items$item)
  vapply(idx, function(i) {
    if (is.na(i)) return(FALSE)
    z <- x$est$thr$anchored[x$est$thr$item == i]
    length(z) && all(z)
  }, logical(1))
}

.equate_interval_halfwidth <- function(tab, alpha = 0.05) {
  out <- rep(NA_real_, nrow(tab))
  if (!all(c("se_diff", "df") %in% names(tab))) return(out)
  ok <- is.finite(tab$se_diff) & tab$se_diff >= 0 &
    !is.na(tab$df) & tab$df > 0
  if (any(ok))
    out[ok] <- stats::qt(1 - alpha / 2, df = tab$df[ok]) * tab$se_diff[ok]
  out
}

.equate_bank_cov <- function(reference, ids) {
  C <- attr(reference, "cov_location", exact = TRUE)
  if (is.null(C)) return(NULL)
  if (!is.matrix(C) || !is.numeric(C) || any(!is.finite(C)) ||
      !identical(dim(C), c(length(ids), length(ids))))
    stop("attr(reference, 'cov_location') must be a finite numeric square ",
         "matrix with one row and column per bank item")
  if (!.covariance_is_symmetric(C))
    stop("attr(reference, 'cov_location') must be symmetric")
  if (!is.null(rownames(C)) || !is.null(colnames(C))) {
    if (is.null(rownames(C)) || is.null(colnames(C)) ||
        anyNA(match(ids, rownames(C))) || anyNA(match(ids, colnames(C))))
      stop("named bank covariance rows and columns must match every bank item")
    C <- C[ids, ids, drop = FALSE]
  }
  if (!.covariance_is_psd(C))
    stop("attr(reference, 'cov_location') must be positive semidefinite")
  C
}

.equate_ref <- function(reference) {
  if (inherits(reference, "rasch"))
    return(data.frame(item = reference$items$item,
                      location = reference$items$location,
                      se = reference$items$se,
                      max = reference$items$max))
  .check_column_names(reference)
  reference <- as.data.frame(reference)
  .check_column_names(reference)
  if (!all(c("item", "location") %in% names(reference)))
    stop("reference needs columns item, location (and ideally se)")
  # Each default below is a length-one value, which `$<-.data.frame` refuses
  # on a zero-row frame with an implementation message. A reference with no
  # items cannot link anything: say so before filling its optional columns.
  if (nrow(reference) == 0L)
    stop("reference has no items; it needs at least one item to equate")
  if (!"se" %in% names(reference)) reference$se <- NA_real_
  if (!"max" %in% names(reference)) reference$max <- NA_integer_
  out <- reference[, c("item", "location", "se", "max")]
  out$item <- .role_text_values(out$item)
  out$location <- .bank_numeric(out$location, "location")
  out$se <- .bank_numeric(out$se, "se")
  out$max <- .bank_numeric(out$max, "max")
  if (anyNA(out$item) || any(!nzchar(out$item)))
    stop("reference item names must be non-missing and non-empty")
  if (anyDuplicated(out$item))
    stop("reference item names must be unique: ",
         paste(unique(out$item[duplicated(out$item)]), collapse = ", "))
  if (any(!is.finite(out$location)))
    stop("reference locations must be finite")
  if (any(!is.na(out$se) & (!is.finite(out$se) | out$se < 0)))
    stop("reference standard errors must be non-negative finite values or NA")
  if (any(!is.na(out$max) & (!is.finite(out$max) | out$max < 1 |
                             out$max != floor(out$max))))
    stop("reference max values must be positive whole numbers or NA")
  out
}

# Common locations with zero pooled variance are fixed constants, not noisy
# estimates.  They therefore have to define one coherent origin.  Averaging
# incompatible exact anchors would move every estimated item away from its
# proper scale and then make the free items look as though they drifted, while
# the contradictory anchors themselves receive NA Wald probabilities because
# their contrast variance is zero.  Refuse that premise before estimating the
# shift.  On a declared fixed origin, even one non-zero exact difference is a
# contradiction.
.check_exact_equating_origin <- function(d, exact, estimates_shift,
                                         label = "item") {
  if (!is.logical(exact) || length(exact) != length(d) || anyNA(exact))
    stop("internal exact-anchor mask is invalid", call. = FALSE)
  fixed <- exact & is.finite(d)
  if (!any(fixed)) return(invisible(NULL))
  exact <- d[fixed]
  tol <- sqrt(.Machine$double.eps) * max(1, abs(exact))
  inconsistent <- if (estimates_shift)
    length(exact) > 1L && diff(range(exact)) > tol
  else any(abs(exact) > tol)
  if (inconsistent) {
    noun <- if (sum(fixed) == 1L) label else paste0(label, "s")
    stop("common ", noun, " with zero pooled uncertainty imply incompatible ",
         if (estimates_shift) "origin shifts" else "values on the declared fixed origin",
         "; reconcile the external anchors before equating", call. = FALSE)
  }
  invisible(NULL)
}

#' Equate two test calibrations through their common items
#'
#' Places two calibrations on a common origin using their shared items, then
#' tests the shared items for drift. The reference may be a fitted model or an
#' item bank.
#'
#' @details
#' Let \eqn{d_j} be the location difference for common item \eqn{j} and
#' \eqn{v_j} its marginal variance. With \code{shift = "mean"} the scale
#' shift is the precision-weighted mean
#' \deqn{\hat s=\frac{\sum_j d_j/v_j}{\sum_j 1/v_j},}
#' and each item is tested using \eqn{d_j-\hat s} with a variance that
#' accounts for the estimated shift through the items' joint covariance.
#' If fewer than two common items have usable variances but at least two have
#' finite locations, the function returns their unweighted mean difference as
#' a descriptive fallback and records \code{shift_method = "unweighted"}.
#' An exact common anchor determines the shift even when it is the only
#' common item with usable uncertainty; unavailable SEs do not override it.
#' When the shift is estimated, drift inference requires independent
#' calibrations and at least three common items with usable,
#' positive-semidefinite joint covariance information. With
#' \code{shift = "none"}, the origin is fixed before the comparison and each
#' item's variance is the sum of its two marginal variances; joint covariance
#' information and a three-item link are then unnecessary. One common item is
#' sufficient for that fixed-origin comparison; estimating a shift still
#' requires at least two. Otherwise the function returns a descriptive link.
#' A fitted calibration's empirical covariance must also pass the
#' informative-person count, effective-support and projected-rank checks used
#' by \code{\link{rasch}}. Supported independent rows use the limiting normal
#' reference.
#' When either covariance comes from repeated person clusters, drift
#' probabilities use contrast-specific Welch--Satterthwaite degrees of
#' freedom. The corresponding residual degrees of freedom are the number of
#' independent person clusters minus one. A fixed anchor contributes zero
#' variance and does not consume cluster degrees of freedom.
#' Fitted calibrations must have converged.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param reference A second \code{\link{rasch}} fit, or a data frame with
#'   columns \code{item}, \code{location}, and optionally \code{se}. Item
#'   names and column names must be unique. Numeric fields may be numeric
#'   columns, numeric text, or factors with numeric labels; other column
#'   classes are refused. Locations must be finite. For bank-based drift inference
#'   with an estimated mean shift, attach the bank's joint item-location
#'   covariance as a square matrix in \code{attr(reference, "cov_location")},
#'   ordered like the bank rows (or named by item); marginal SEs alone do not
#'   carry the centring covariance. They are sufficient with
#'   \code{shift = "none"}.
#'   A bank treated as fixed may instead have zero SEs. A polytomous bank must
#'   also include \code{max}, the maximum item score. A bank covariance
#'   estimated from a finite number of independent sampling units may carry
#'   its positive residual degrees of freedom in
#'   \code{attr(reference, "df_location")}.
#' @param shift \code{"mean"} (default) allows a scale shift between the two
#'   analyses; \code{"none"} compares raw locations, appropriate when both
#'   analyses are already on a shared (anchored) scale.
#' @param independent Whether the two calibrations use independent sampling
#'   units. For two fitted objects this must be stated explicitly: the default
#'   \code{NULL} withholds drift tests because cross-fit covariance is
#'   otherwise unknown. A bank table is treated as independent unless
#'   \code{FALSE} is supplied. When \code{FALSE}, descriptive equating is
#'   returned but inferential drift columns are withheld.
#' @return A list with the comparison \code{table} (locations, standard
#'   errors, difference, its standard error, t statistic, reference degrees
#'   of freedom, raw and Holm-adjusted p, drift flag), the
#'   estimated \code{shift} and its \code{shift_method},
#'   the location \code{correlation}, the root mean square difference after
#'   shifting (\code{rmsd}), the number of common items \code{n_common}, the
#'   number with usable standard errors \code{n}, and whether drift inference
#'   was available (\code{inferential}). The \code{note} component records
#'   exclusions and the reason inference was withheld, where applicable.
#'   An individual drift probability is also withheld when its contrast has
#'   zero estimated uncertainty.
#'   Common items with unavailable drift probabilities remain in the
#'   multiplicity family.
#' @examples
#' set.seed(1); d <- seq(-1.5, 1.5, length.out = 8)
#' mk <- function() {
#'   X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
#'   colnames(X) <- paste0("I", 1:8); rasch(X)
#' }
#' eq <- equate_tests(mk(), mk(), independent = TRUE)
#' eq$table
#' @export
equate_tests <- function(fit, reference, shift = c("mean", "none"),
                         independent = NULL) {
  shift <- match.arg(shift)
  if (!is.null(independent) && (length(independent) != 1L ||
      !is.logical(independent) || !is.null(dim(independent)) ||
      !is.null(oldClass(independent)) || is.na(independent)))
    stop("independent must be NULL, TRUE, or FALSE")
  if (!inherits(fit, "rasch") ||
      inherits(fit, "rasch_efrm") || inherits(fit, "rasch_mfrm"))
    stop("fit must be an ordinary person-by-item Rasch calibration")
  if (inherits(reference, "rasch_btl"))
    stop("reference must be an ordinary Rasch calibration or item bank; ",
         "a paired-comparison calibration is on a different response scale")
  if (inherits(reference, "rasch_efrm") || inherits(reference, "rasch_mfrm"))
    stop("reference must be an ordinary Rasch calibration or item bank")
  # under an explanatory design an item location is not an item parameter:
  # items sharing a design cell are forced to one location and carry the
  # design coefficients' uncertainty, so drift can be neither localised nor
  # tested at its stated size
  if (inherits(fit, "rasch_explanatory") ||
      inherits(reference, "rasch_explanatory"))
    stop("item drift is not defined for an explanatory calibration: the ",
         "item locations are functions of their predictors, so a drifted ",
         "item is smeared over every item sharing its design cell and the ",
         "standard errors are the design coefficients'. Equate the ",
         "unrestricted calibrations, or test the design coefficients")
  if (!isTRUE(fit$est$converged))
    stop("the current calibration did not converge; equating is unavailable")
  if (inherits(reference, "rasch") && !isTRUE(reference$est$converged))
    stop("the reference calibration did not converge; equating is unavailable")
  ref <- .equate_ref(reference)
  cur <- data.frame(item = fit$items$item, location = fit$items$location,
                    se = fit$items$se, max = fit$items$max)
  common <- intersect(cur$item, ref$item)
  min_common <- if (identical(shift, "mean")) 2L else 1L
  if (length(common) < min_common)
    stop(if (identical(shift, "mean"))
      "need at least two common items to estimate an origin shift" else
      "need at least one common item for a fixed-origin comparison")
  a <- cur[match(common, cur$item), ]
  b <- ref[match(common, ref$item), ]
  bank_cov <- if (inherits(reference, "rasch")) NULL else
    .equate_bank_cov(reference, ref$item)
  # Validate supplied finite-sample metadata even if another design condition
  # later withholds inference; it must not be silently ignored in one branch.
  df1 <- .equate_cov_df(fit)
  df2 <- .equate_cov_df(reference)
  if (!is.null(bank_cov)) {
    cov_se <- sqrt(pmax(diag(bank_cov), 0))
    stated <- is.finite(ref$se)
    if (any(stated & !.se_covariance_agree(ref$se, cov_se)))
      stop("the bank standard errors must agree with the diagonal of ",
           "attr(reference, 'cov_location')")
    # The joint covariance supplies every marginal variance. Retain a stated
    # SE after checking it against the diagonal, and complete any genuinely
    # missing entries from that same diagonal.
    ref$se[!stated] <- cov_se[!stated]
    b <- ref[match(common, ref$item), ]
  }
  known_max <- is.finite(b$max)
  if (any(known_max & a$max != b$max))
    stop("common items use different maximum scores: ",
         paste(common[known_max & a$max != b$max], collapse = ", "))
  if (any(a$max > 1L) && !all(known_max))
    stop("a polytomous item bank must include a 'max' column so scoring-scale ",
         "compatibility can be checked")
  d <- a$location - b$location
  v <- a$se^2 + b$se^2
  finite_loc <- is.finite(d)
  exact1 <- .equate_exact_locations(fit, common)
  exact2 <- .equate_exact_locations(reference, common)
  exact_reference <- if (inherits(reference, "rasch")) exact2 else
    is.finite(b$se) & b$se == 0
  exact_pair <- exact1 & exact_reference
  # A numerical zero from an estimated sandwich is not infinite precision.
  # Admit zero pooled variance only when both locations are declared exact;
  # otherwise it must neither dominate the origin shift nor enter a Wald
  # contrast. Unavailable rows remain in the output with NA test columns.
  usable <- is.finite(d) & is.finite(v) & (v > 0 | exact_pair)
  zero_not_exact <- is.finite(d) & is.finite(v) & v == 0 & !exact_pair
  notes <- character(0)
  w <- numeric(0)
  .check_exact_equating_origin(d, exact_pair,
                               estimates_shift = identical(shift, "mean"),
                               label = "item")
  if (any(zero_not_exact))
    notes <- c(notes, sprintf(paste(
      "common item(s) with zero pooled estimated variance but without exact",
      "anchors were excluded from weighting and drift tests: %s"),
      paste(common[zero_not_exact], collapse = ", ")))
  if (shift == "mean" && (sum(usable) >= 2L || any(exact_pair & usable))) {
    w <- .inverse_variance_weights(v[usable])
    c0 <- sum(w * d[usable]) / sum(w)
    shift_method <- "precision-weighted"
  } else if (shift == "mean") {
    if (sum(finite_loc) < 2L)
      stop("need at least two common items with finite locations to estimate ",
           "an origin shift")
    c0 <- mean(d[finite_loc])
    shift_method <- "unweighted"
    notes <- c(notes, paste(
      "fewer than two common items had usable variances; the reported shift",
      "is the unweighted mean of the finite location differences and is descriptive"))
  } else {
    c0 <- 0
    shift_method <- "none"
  }
  if (shift_method == "precision-weighted" && any(!usable)) {
    notes <- c(notes, sprintf(
      "common item(s) excluded from the precision-weighted shift and drift tests (location or uncertainty unavailable or unusable): %s",
      paste(common[!usable], collapse = ", ")))
  } else if (shift_method == "unweighted") {
    no_se <- finite_loc & !usable
    if (any(no_se)) notes <- c(notes, sprintf(
      "common item(s) included in the descriptive shift but excluded from drift tests because their uncertainty was unavailable or was a non-exact zero: %s",
      paste(common[no_se], collapse = ", ")))
    if (any(!finite_loc)) notes <- c(notes, sprintf(
      "common item(s) excluded from the shift and drift tests because their locations were unavailable: %s",
      paste(common[!finite_loc], collapse = ", ")))
  } else if (any(!usable)) {
    notes <- c(notes, sprintf(
      "drift tests unavailable for common item(s) with unavailable or unusable locations or uncertainty: %s",
      paste(common[!usable], collapse = ", ")))
  }
  independent_ok <- if (is.null(independent)) !inherits(reference, "rasch")
                    else isTRUE(independent)
  estimates_shift <- identical(shift, "mean")
  joint_cov_ok <- !estimates_shift || inherits(reference, "rasch") ||
    !is.null(bank_cov) || all(b$se[usable] == 0)
  min_inference <- if (estimates_shift) 3L else 1L
  # An unavailable cluster df can be bypassed only by a genuinely exact
  # anchored location, never merely by a numerical zero SE. This prevents a
  # legacy one-cluster sandwich from masquerading as a fixed calibration.
  df_ok_1 <- !is.na(df1) | exact1
  df_ok_2 <- !is.na(df2) | exact2
  df_available <- if (estimates_shift)
    (all(df_ok_1[usable]) && all(df_ok_2[usable])) else
      any(usable & df_ok_1 & df_ok_2)
  if (is.null(independent) && inherits(reference, "rasch"))
    notes <- c(notes, paste(
      "drift tests withheld because independence between the two fitted",
      "calibrations was not stated; set independent = TRUE only for",
      "independent sampling units"))
  if (identical(independent, FALSE))
    notes <- c(notes, paste(
      "drift tests withheld for dependent calibrations because cross-fit",
      "covariance is not available; use a joint or paired bootstrap"))
  if (independent_ok && sum(usable) >= min_inference && !joint_cov_ok)
    notes <- c(notes, paste(
      "drift tests withheld because a bank with non-zero marginal SEs needs",
      "its joint item-location covariance in attr(reference, 'cov_location');",
      "marginal SEs do not carry the calibration-origin covariance"))
  if (independent_ok && sum(usable) < min_inference) {
    notes <- c(notes, if (estimates_shift) paste(
      "drift tests withheld because at least three common items with",
      "standard errors are needed to distinguish item drift from the link")
    else paste(
      "drift tests withheld because no common item has locations and",
      "standard errors in both calibrations"))
  }
  if (independent_ok && sum(usable) >= min_inference && !df_available)
    notes <- c(notes, paste(
      "drift tests are withheld because a contributing calibration does not",
      "have enough independent-person support for inference"))
  unsupported_df <- usable & !(df_ok_1 & df_ok_2)
  if (!estimates_shift && any(unsupported_df) && df_available)
    notes <- c(notes, sprintf(paste(
      "drift tests are unavailable for common item(s) whose contributing",
      "calibration lacks independent-person support: %s"),
      paste(common[unsupported_df], collapse = ", ")))
  inferential <- independent_ok && sum(usable) >= min_inference &&
    joint_cov_ok && df_available
  # the shift c0 is estimated from the same common items the drift tests
  # then examine, and each calibration's locations are correlated through
  # its identification constraint: the drift denominators use
  # Var(d_i - c0) = [(I - 1u') Sigma (I - u 1')]_ii over the usable items,
  # with Sigma the sum of the two calibrations' item-location covariances
  var_d <- rep(NA_real_, length(common))
  var_1 <- var_2 <- rep(NA_real_, length(common))
  if (inferential && !estimates_shift) {
    # The origin was fixed outside this comparison. No common-item shift is
    # estimated, so cross-item covariance cannot enter an individual drift
    # contrast; its variance is simply Var(delta_1) + Var(delta_2).
    test_support <- usable & df_ok_1 & df_ok_2
    var_1[test_support] <- a$se[test_support]^2
    var_2[test_support] <- b$se[test_support]^2
    var_d[test_support] <- pmax(var_1[test_support] + var_2[test_support], 0)
  } else if (inferential) {
    S1all <- tryCatch(.equate_loc_cov(fit, common), error = function(e) NULL)
    S2all <- tryCatch(.equate_loc_cov(reference, common),
                      error = function(e) NULL)
    S1 <- if (is.null(S1all)) NULL else S1all[usable, usable, drop = FALSE]
    S2 <- if (is.null(S2all)) NULL else S2all[usable, usable, drop = FALSE]
    covariance_ok <- .covariance_supports_wald(S1, sum(usable)) &&
      .covariance_supports_wald(S2, sum(usable))
    if (!covariance_ok) {
      inferential <- FALSE
      notes <- c(notes, paste(
        "drift tests withheld because a fitted item-location covariance is",
        "unavailable or not positive semidefinite"))
    } else {
      u <- w / sum(w)
      Hc <- diag(length(u)) - matrix(u, nrow = length(u), ncol = length(u),
                                    byrow = TRUE)
      V1 <- Hc %*% S1 %*% t(Hc)
      V2 <- Hc %*% S2 %*% t(Hc)
      var_1[usable] <- pmax(diag(V1), 0)
      var_2[usable] <- pmax(diag(V2), 0)
      var_d[usable] <- pmax(var_1[usable] + var_2[usable], 0)
    }
  }
  se_diff <- sqrt(var_d)
  t <- if (inferential) .wald_ratio(d - c0, se_diff) else
    rep(NA_real_, length(common))
  df <- rep(NA_real_, length(common))
  if (inferential) {
    den <- if (is.finite(df1)) var_1^2 / df1 else
      rep(0, length(common))
    if (is.finite(df2)) den <- den + var_2^2 / df2
    supported <- is.finite(var_1) & is.finite(var_2)
    df[supported] <- ifelse(den[supported] > 0,
      (var_1[supported] + var_2[supported])^2 / den[supported], Inf)
  }
  p <- 2 * stats::pt(-abs(t), df = df)
  n <- sum(usable)
  p_adj <- rep(NA_real_, length(p))
  testable <- inferential & usable & is.finite(p)
  if (any(testable)) p_adj[testable] <- p.adjust(
    p[testable], method = "holm", n = length(common))
  zero_uncertainty <- inferential & usable & !is.finite(p)
  if (any(zero_uncertainty))
    notes <- c(notes, paste0(
      "drift probability/probabilities withheld for zero contrast ",
      "uncertainty: ",
      paste(common[zero_uncertainty], collapse = ", ")))
  inferential <- inferential && any(testable)
  tab <- data.frame(item = common,
                    location_1 = a$location, se_1 = a$se,
                    location_2 = b$location, se_2 = b$se,
                    difference = d, adj_difference = d - c0,
                    se_diff = se_diff, t = t, df = df,
                    p = p, p_adj = p_adj,
                    drift = ifelse(is.na(p_adj), NA, p_adj < 0.05))
  rownames(tab) <- NULL
  cor_link <- if (sum(finite_loc) >= 2L &&
                   stats::sd(a$location[finite_loc]) > 0 &&
                   stats::sd(b$location[finite_loc]) > 0)
    stats::cor(a$location[finite_loc], b$location[finite_loc]) else NA_real_
  if (any(testable & is.finite(df)))
    notes <- c(notes, paste(
      "drift probabilities use contrast-specific Welch-Satterthwaite",
      "degrees of freedom for finite person-cluster covariances"))
  structure(class = "rasch_equate", list(table = tab, shift = c0,
       shift_method = shift_method,
       correlation = cor_link,
       rmsd = sqrt(mean((d[finite_loc] - c0)^2)), n = n,
       n_common = length(common), inferential = inferential,
       note = if (length(notes)) paste(notes, collapse = "; ") else NULL))
}

#' Plot a test-equating comparison
#'
#' Scatter of the two calibrations' common-item locations with the shifted
#' identity line, per-item 95 per cent contrast intervals, and a dotted guide
#' band at their average half-width; drifting items
#' (Holm-adjusted) are highlighted and labelled.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param reference A second \code{\link{rasch}} fit, or a data frame with
#'   columns \code{item}, \code{location}, and optionally \code{se}; a
#'   polytomous bank also needs \code{max}.
#' @param shift Passed to \code{\link{equate_tests}}.
#' @param independent Passed to \code{\link{equate_tests}}.
#' @return Called for its plotting side effect; invisibly the
#'   \code{\link{equate_tests}} result.
#' @examples
#' set.seed(1); d <- seq(-1.5, 1.5, length.out = 8)
#' mk <- function() {
#'   X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
#'   colnames(X) <- paste0("I", 1:8); rasch(X)
#' }
#' plot_equate(mk(), mk(), independent = TRUE)
#' @export
plot_equate <- function(fit, reference, shift = c("mean", "none"),
                        independent = NULL) {
  eq <- equate_tests(fit, reference, shift, independent = independent)
  tab <- eq$table
  paired <- is.finite(tab$location_1) & is.finite(tab$location_2)
  if (!any(paired))
    .refuse("no common item has finite locations in both calibrations; ",
            "there is nothing to display")
  rng <- range(c(tab$location_1[paired], tab$location_2[paired])) +
    c(-0.4, 0.4)
  op <- .rr_canvas(rng, rng, "Reference location (logits)",
                   "Current location (logits)",
                   sprintf("%d common items, shift %.3f, r = %.3f",
                           eq$n_common, eq$shift, eq$correlation),
                   grid_x = TRUE)
  on.exit(par(op))
  abline(eq$shift, 1, col = .rr$ink, lwd = 2)
  # Excluded items still appear as points. Inferential intervals use each
  # contrast's own Welch--Satterthwaite critical value (or the normal limit
  # when df is infinite), and disappear when inference was withheld.
  half <- .equate_interval_halfwidth(tab)
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
  dr <- tab$drift %in% TRUE
  points(tab$location_2, tab$location_1, pch = 21, cex = 1.6,
         bg = ifelse(dr, .rr$red, .rr$blue), col = "white", lwd = 1.2)
  if (any(dr))
    text(tab$location_2[dr], tab$location_1[dr],
         tab$item[dr], pos = 3, offset = 0.5, cex = 0.75, col = .rr$red)
  invisible(eq)
}


#' @export
print.rasch_equate <- function(x, ...) {
  method <- if (is.null(x$shift_method)) "method unavailable" else x$shift_method
  cat(sprintf("Common-item equating over %d item(s): shift %.3f (%s), correlation %.3f, RMSD %.3f\n",
              x$n_common, x$shift, method, x$correlation, x$rmsd))
  core <- c("item", "location_1", "location_2", "adj_difference", "t",
            "p_adj", "drift")
  print(.fmt_df(x$table[, intersect(core, names(x$table))]), row.names = FALSE)
  if (!isTRUE(x$inferential)) cat("Drift inference withheld; see note below.\n")
  cat("(standard errors and unadjusted columns on $table)\n")
  if (!is.null(x$note)) cat("Note:", x$note, "\n")
  invisible(x)
}
