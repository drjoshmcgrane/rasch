# rasch :: top-level analysis
# ===========================================================================
# rasch() ties the engine together: data preparation (with category
# collapsing and constant-item removal, recorded as notes), pairwise
# conditional ML item estimation, Warm WLE person estimation per
# missing-data pattern, the full test-of-fit suite, and the score-to-measure
# table. An ID variable and any number of person factors carry through to
# the person estimates.
# ===========================================================================

# Map missing-data codes (and any negative score) to NA. Valid item scores
# are non-negative integers from zero; by long-standing convention -1 marks a
# missing response, so any value below zero is read as missing.
.check_na_codes <- function(na_codes) {
  if (!(is.numeric(na_codes) || is.character(na_codes)) ||
      is.complex(na_codes) || !is.null(dim(na_codes)) ||
      !is.null(oldClass(na_codes)) || anyNA(na_codes))
    stop("`na_codes` must be a plain numeric or character vector without missing values",
         call. = FALSE)
  if (is.numeric(na_codes) &&
      (any(!is.finite(na_codes)) || any(na_codes != floor(na_codes)) ||
       any(abs(na_codes) > .Machine$integer.max)))
    stop("numeric `na_codes` must be finite integer score values",
         call. = FALSE)
  if (is.character(na_codes) &&
      any(!nzchar(trimws(na_codes))))
    stop("character `na_codes` must be non-empty values", call. = FALSE)
  invisible(na_codes)
}

# Match declared codes before numeric score validation. This preserves text
# codes such as "." and "09", while also treating numerically equivalent
# representations (9, "9", and "09") alike. Undeclared negative values are
# handled separately because valid Rasch scores begin at zero.
.missing_code_mask <- function(v, na_codes) {
  .check_na_codes(na_codes)
  txt <- as.character(v)
  codes <- as.character(na_codes)
  exact <- !is.na(v) & txt %in% codes
  vn <- suppressWarnings(as.numeric(txt))
  cn <- suppressWarnings(as.numeric(codes))
  cn <- cn[is.finite(cn)]
  numeric_match <- !is.na(vn) & is.finite(vn) & length(cn) > 0L & vn %in% cn
  exact | numeric_match
}

.apply_na_codes <- function(v, na_codes) {
  v[v %in% na_codes | (!is.na(v) & v < 0)] <- NA
  v
}

# Prepare the item matrix: integer scores from 0, consecutive observed
# categories, no constant items. Returns the matrix plus human-readable notes.
# `anchors`, with item names already resolved, exempts the items whose every
# threshold is fixed from the PCM category merge: nothing is estimated for
# them, so they need no conditional information. `scored_items` have already
# had their raw missing codes removed by keyed scoring; a generated score
# must not be mistaken for a raw missing code with the same numeric value.
.prepare_X <- function(X, na_codes = -1, model = "PCM", anchors = NULL,
                        scored_items = character(0)) {
  .check_na_codes(na_codes)
  notes <- character(0)
  X <- as.matrix(X)
  if (is.complex(X))
    stop("complex response scores are not supported; scores must be real integer counts",
         call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- sprintf("I%02d", seq_len(ncol(X)))
  if (anyNA(colnames(X)) || any(!nzchar(trimws(colnames(X)))))
    stop("item column names must be non-missing and non-empty (not whitespace-only)")
  if (anyDuplicated(colnames(X)))
    stop("item column names must be unique: ",
         paste(unique(colnames(X)[duplicated(colnames(X))]), collapse = ", "))
  raw_code <- matrix(.missing_code_mask(as.vector(X), na_codes),
                     nrow(X), ncol(X), dimnames = dimnames(X))
  raw_code[, colnames(X) %in% scored_items] <- FALSE
  Xn <- suppressWarnings(apply(X, 2, function(col) as.numeric(as.character(col))))
  dim(Xn) <- dim(X); dimnames(Xn) <- dimnames(X)
  too_large <- !raw_code & !is.na(Xn) & is.finite(Xn) &
    (Xn > .Machine$integer.max | Xn < -.Machine$integer.max)
  if (any(too_large))
    stop("score(s) outside the supported integer range in: ",
         paste(colnames(X)[colSums(too_large) > 0], collapse = ", "),
         "; rescore the response categories before analysis", call. = FALSE)
  Xi <- suppressWarnings(apply(X, 2, function(col) as.integer(as.character(col))))
  dim(Xi) <- dim(X); dimnames(Xi) <- dimnames(X)
  # An infinite numeric value is corrupt response data, not a missing-data
  # marker. as.integer(Inf) returns NA, which previously sent it down the
  # ordinary "non-numeric entries set to missing" path and silently removed it.
  # `NaN` is also `NA` in R, so test the converted values directly rather
  # than guarding the check with !is.na(X). Genuine NA is not NaN.
  nonfinite <- !raw_code & (is.infinite(Xn) | is.nan(Xn))
  if (any(nonfinite))
    stop("non-finite score(s) in: ",
         paste(colnames(X)[colSums(nonfinite) > 0], collapse = ", "),
         "; use NA or a declared missing-data code for missing responses",
         call. = FALSE)
  # as.integer() TRUNCATES fractional values (1.9 -> 1) without a warning:
  # that silently alters response data, so it must be an error, not a note
  fractional <- !raw_code & !is.na(Xn) & !is.na(Xi) & Xn != Xi
  frac <- colSums(fractional) > 0
  if (any(frac))
    stop("non-integer score(s) in: ",
         paste(colnames(X)[frac], collapse = ", "),
         " (e.g. ", format(Xn[fractional][1]),
         "); Rasch categories are integer counts -- round or rescore ",
         "explicitly before analysis")
  bad_num <- colSums(!is.na(X) & !raw_code & is.na(Xi)) > 0
  if (any(bad_num))
    notes <- c(notes, paste0("non-numeric entries set to missing in: ",
                             paste(colnames(X)[bad_num], collapse = ", ")))
  n_na <- sum(raw_code | (!is.na(Xi) & Xi < 0))
  Xi[raw_code] <- NA_integer_
  # Declared codes were matched to the raw values above. At this point only
  # the negative-score convention remains; reapplying positive codes would
  # erase valid scores produced by a key.
  Xi[] <- .apply_na_codes(Xi, integer(0))
  if (n_na > 0) {
    codes <- paste(unique(c(na_codes, "negative")), collapse = ", ")
    notes <- c(notes, sprintf("%d cell(s) with a missing-data code (%s) set to missing",
                              n_na, codes))
  }
  X <- Xi
  const <- apply(X, 2, function(col) length(unique(col[!is.na(col)])) < 2)
  if (any(const)) {
    notes <- c(notes, paste0("dropped constant item(s): ",
                             paste(colnames(X)[const], collapse = ", ")))
    X <- X[, !const, drop = FALSE]
  }
  if (ncol(X) < 2) stop("need at least two non-constant items")
  for (i in seq_len(ncol(X))) {
    v <- X[, i]; obs <- sort(unique(v[!is.na(v)]))
    full <- seq(0L, max(obs))
    if (!identical(obs, full)) {
      X[, i] <- match(v, obs) - 1L
      notes <- c(notes, sprintf("item %s rescored: observed categories [%s] mapped to 0:%d",
                                colnames(X)[i], paste(obs, collapse = ","),
                                length(obs) - 1L))
    }
  }
  if (model == "PCM") {
    merged <- .merge_uninformative_categories(X, keep = .fully_anchored(X, anchors))
    X <- merged$X; notes <- c(notes, merged$notes)
    if (ncol(X) < 2) stop("need at least two non-constant items")
  }
  list(X = X, notes = notes)
}

# The items whose every threshold `anchors` fixes: a numeric k on each of
# 1..m_i, or a location anchor (k = NA) on a dichotomous item, which pcml()
# converts to its single threshold. Average anchoring fixes nothing, and a
# location anchor on a polytomous item leaves its thresholds free, so
# neither exempts an item. Anchor items are matched by name.
.fully_anchored <- function(X, anchors) {
  if (is.null(anchors) || !is.data.frame(anchors) ||
      !all(c("item", "k") %in% names(anchors)) ||
      ("average" %in% names(anchors) && any(anchors$average %in% TRUE)))
    return(character(0))
  a_item <- as.character(anchors$item)
  present <- intersect(unique(a_item), colnames(X))
  if (!length(present)) return(character(0))
  fixed <- vapply(present, function(nm) {
    mi <- max(X[, nm], na.rm = TRUE)
    k <- anchors$k[a_item == nm]
    if (all(is.na(k))) mi == 1L
    else !anyNA(k) && setequal(k, seq_len(mi))
  }, TRUE)
  present[fixed]
}

# The categories of each item that at least one informative response pattern
# observes: a person with two or more observed items whose raw score lies
# strictly between zero and the maximum. Only those patterns enter the
# pairwise conditional likelihood, so only those categories carry
# conditional information about the item's thresholds.
.conditional_categories <- function(X) {
  obs_mask <- !is.na(X)
  mx <- apply(X, 2, max, na.rm = TRUE)
  raw <- rowSums(X, na.rm = TRUE)
  max_raw <- as.numeric(obs_mask %*% mx)
  informative <- rowSums(obs_mask) >= 2L & raw > 0 & raw < max_raw
  lapply(seq_len(ncol(X)), function(i) {
    v <- X[, i]
    sort(unique(v[!is.na(v) & informative]))
  })
}

# Categories that occupy a non-zero interval on the upper envelope of the
# partial-credit category logits. Evaluating a regular theta grid can miss a
# genuinely modal category when two adjacent thresholds are close together.
# Category k beats every lower category above the largest mean of the
# thresholds leading up to it, and every higher category below the smallest
# mean of the thresholds leading on from it; it is modal when that interval
# has positive width. The tolerance folds floating-point noise: thresholds
# that are equal by construction give crossing points that differ by an ulp,
# and probing between them would report a category modal at a single point.
.modal_score_categories <- function(tau, tol = 1e-9) {
  m <- length(tau)
  if (!m) return(0L)
  cs <- c(0, cumsum(tau))
  modal <- c(TRUE, vapply(seq_len(m - 1L), function(k) {
    below <- seq_len(k) - 1L
    above <- (k + 1L):m
    lower <- max((cs[k + 1L] - cs[below + 1L]) / (k - below))
    upper <- min((cs[above + 1L] - cs[k + 1L]) / (above - k))
    is.finite(lower) && is.finite(upper) &&
      upper - lower > tol * max(1, abs(lower), abs(upper))
  }, logical(1)), TRUE)
  which(modal) - 1L
}

# A category observed only in extreme response patterns -- every observed
# item at its minimum, or every one at its maximum -- or only by persons who
# answered a single item carries no pairwise conditional information, exactly
# as an unobserved category does: its PCM threshold diverges and the projected
# information matrix goes singular. Merge such a category into its neighbour
# and say so. A merge can lower a maximum score, and with it change which
# persons are extreme, so repeat until the coding is stable. Under the RSM the
# thresholds are shared across items and the category stays identified, so
# the caller skips this step. Items named in `keep` -- every threshold fixed
# by an anchor -- are left as coded: their thresholds are not estimated, so
# an uninformative category costs them nothing, and pcml() exempts them
# from its own uninformative-category check for the same reason.
.merge_uninformative_categories <- function(X, keep = character(0)) {
  notes <- character(0)
  repeat {
    changed <- FALSE
    mx <- apply(X, 2, max, na.rm = TRUE)
    cc <- .conditional_categories(X)
    drop <- logical(ncol(X))
    for (i in seq_len(ncol(X))) {
      if (colnames(X)[i] %in% keep) next
      v <- X[, i]
      cond <- cc[[i]]
      full <- seq(0L, mx[i])
      if (identical(cond, full)) next
      changed <- TRUE
      if (length(cond) < 2L) { drop[i] <- TRUE; next }
      # each uninformative category joins the nearest informative category
      # below it (above it at the bottom of the scale)
      target <- vapply(full, function(k) {
        lower <- cond[cond <= k]
        if (length(lower)) max(lower) else min(cond)
      }, 0L)
      X[, i] <- (match(target, cond) - 1L)[v + 1L]
      lost <- setdiff(full, cond)
      notes <- c(notes, sprintf(paste0(
        "item %s rescored: categor%s %s observed only in extreme response ",
        "patterns (no conditional information) merged with the adjacent ",
        "categor%s; categories mapped to 0:%d"),
        colnames(X)[i], if (length(lost) > 1L) "ies" else "y",
        paste(lost, collapse = ","), if (length(lost) > 1L) "ies" else "y",
        length(cond) - 1L))
    }
    if (any(drop)) {
      notes <- c(notes, paste0(
        "dropped item(s) with no conditional information (every response ",
        "outside one category comes from an extreme response pattern): ",
        paste(colnames(X)[drop], collapse = ", ")))
      X <- X[, !drop, drop = FALSE]
      if (ncol(X) < 2) break
    }
    if (!changed) break
  }
  list(X = X, notes = notes)
}

#' Fit a Rasch model
#'
#' Fits the partial credit model (PCM) or rating scale model (RSM) by pairwise
#' conditional maximum likelihood. Person locations are Warm weighted
#' likelihood estimates. The fitted object contains item and person fit,
#' targeting, reliability, threshold diagnostics, residuals, and a
#' score-to-measure table.
#'
#' @details
#' For scores \eqn{x=0,\ldots,m_i}, the PCM is
#' \deqn{P(X_{ni}=x)=\frac{\exp\{x\theta_n-\sum_{k=1}^{x}\delta_{ik}\}}
#' {\sum_{y=0}^{m_i}\exp\{y\theta_n-\sum_{k=1}^{y}\delta_{ik}\}}.}
#' The RSM constrains \eqn{\delta_{ik}=\beta_i+\tau_k}, where \eqn{\beta_i}
#' is the item location and \eqn{\tau_k} is common across items. Dichotomous
#' items are the one-threshold case of the PCM.
#'
#' Pairwise conditioning removes \eqn{\theta_n} from the item likelihood.
#' Missing responses are omitted from pairwise contributions, and person
#' measures are estimated within each observed item pattern. The observed
#' item-pair graph must identify a common scale. This covers planned linked
#' designs and ignorable missingness; informative missingness can still bias
#' the estimates.
#'
#' Item-parameter uncertainty uses the empirical Godambe sandwich over
#' independent persons, or over person clusters when IDs repeat. It is
#' withheld unless at least 10 contributing units, at least 8 effective units,
#' more effective units than fitted parameters, and a full-rank score
#' covariance support the fitted directions. Effective support reflects the
#' number of informative conditional item pairs contributed by each unit;
#' rows without one do not count. Point estimates and exact anchors remain
#' available when uncertainty is withheld.
#'
#' The fit residual is the log-of-mean-square statistic described by Andrich
#' and Marais (2019, ch. 23). Positive values indicate under-discrimination and
#' negative values indicate over-discrimination. Its standard-normal reading,
#' the item-trait chi-square and the class-interval F test are asymptotic
#' approximations. For ordinary Rasch, PCM and RSM fits,
#' \code{\link{fit_bootstrap}} supplies calibrated probabilities.
#'
#' Multiple-choice responses may be scored from a named item-to-key vector,
#' an item/key table, or an item/option/score table. A slash separates
#' alternative correct options. The third form assigns integer category
#' scores to nominated options and fits the resulting item as polytomous;
#' unlisted options score zero. Raw responses are retained in \code{fit$mc}
#' for distractor analysis.
#'
#' @param data Persons-by-items integer score matrix (categories from 0), or a
#'   data frame also containing ID and person-factor columns. Missing values
#'   are allowed subject to the identification and ignorability conditions
#'   described above.
#' @param model Either \code{"PCM"} (partial credit) or \code{"RSM"} (rating
#'   scale).
#' @param id Optional name of an ID column in \code{data}, or a vector of IDs.
#'   Repeated values cluster the item-parameter sandwich covariance and define
#'   the person unit in repeated-measures DIF. The ordinary item-fit reference
#'   distributions are row-based, so their probabilities are withheld when an
#'   ID occurs on more than one response row. The support conditions described
#'   above then apply to person clusters rather than response rows.
#' @param factors Optional character vector of person-factor column names in
#'   \code{data} (for DIF analysis), a data frame of factors, or one grouping
#'   vector with one entry per data row.
#' @param items Optional item column names or numeric column indices. By
#'   default, columns named in \code{id} or \code{factors} are excluded.
#'   A separate factor data frame may share item names when \code{items} is
#'   explicit. Without it, matching names exclude columns whose values agree;
#'   conflicting values are refused.
#' @param n_groups Number of class intervals for the item-trait chi-square
#'   and ANOVA item fit. The default \code{NULL} applies the rule of Andrich
#'   and Marais (2019, ch. 15): as
#'   many intervals of at least 50 non-extreme persons as the sample allows,
#'   at most 10, at least 2. The resolved value is stored in
#'   \code{fit$n_groups}.
#' @param anchors Optional anchor table for equating: a data frame with
#'   columns \code{item}, \code{k}, and \code{tau}, and optionally
#'   \code{average = TRUE} for average item anchoring; see
#'   \code{\link{pcml}}. Column names must be unique. Anchors determine the
#'   scale origin. Anchor values are treated as fixed, so their uncertainty
#'   is not included in the fitted standard errors.
#' @param na_codes Numeric or character values to read as missing. They are
#'   matched before scores are converted to numbers, including numerically
#'   equivalent labels (for example, \code{"09"} matches a score of 9).
#'   Defaults to \code{-1}, the conventional missing-response code; any
#'   negative score is also treated as missing, since valid category scores
#'   start at zero. For keyed items, codes apply to raw answer options, not
#'   to the scores assigned by the key. Structural refits use the prepared
#'   scores and do not apply the original raw codes again.
#' @param maxit,tol Newton-Raphson iteration cap and convergence
#'   tolerance of the pairwise conditional estimation.
#' @param key Optional multiple-choice key: a named item-to-option vector, an
#'   item/key table, or an item/option/score table. Table column names must be
#'   unique. See Details.
#' @param pc_components \code{NULL} (the default) estimates all PCM thresholds
#'   freely. Values from 1 to 4 use the principal-components form in
#'   \code{\link{pcml_pc}}: location, then spread, skewness, and kurtosis.
#'   This can stabilise sparse categories. Component estimates are stored in
#'   the estimation details. Available for PCM fits without anchors.
#' @section Estimated item discrimination:
#' The item summary includes a post-estimation slope \code{disc}. For item
#' \eqn{i}, it maximises that item's response likelihood over \eqn{a_i} while
#' holding the fitted person locations and thresholds fixed:
#' \deqn{\hat a_i=\arg\max_{a_i}
#'   \sum_n\log P(X_{ni}=x_{ni}\mid\hat\theta_n,\hat\delta_i,a_i).}
#' The same slope multiplies every threshold of a polytomous item. It is a
#' descriptive index, not a freely estimated parameter of the Rasch model,
#' and no sampling standard error or hypothesis test is attached to it.
#' @section Item-fit probabilities:
#' The item-trait chi-square assesses invariance over class intervals, but its
#' asymptotic reference treats the estimated person locations used to form
#' those intervals as known. Its calibration therefore changes with sample
#' size and test length. The class-interval ANOVA and standardised residual
#' readings are approximate for the same reason. The item table retains their
#' raw and Holm-adjusted probabilities as descriptive diagnostics. Each
#' adjustment retains the full item family when one probability is
#' unavailable.
#' \code{\link{fit_bootstrap}} re-estimates every replicate and should be used
#' for item-level inference where it is available. With repeated IDs, the
#' ordinary asymptotic probabilities are withheld and \code{fit_bootstrap()}
#' is unavailable because neither reference models within-person dependence;
#' the residuals and fit statistics remain descriptive.
#' @return An object of class \code{"rasch"}. Its principal components are
#'   the item summary, threshold table, person table, score table, residuals,
#'   reliability, targeting, item-trait statistics, threshold diagnostics,
#'   and estimation details. The component \code{summary_stats} contains the
#'   distribution summaries, fit-location correlations, and the cell
#'   degrees-of-freedom factor. The item summary carries a \code{disc}
#'   column described below. \code{repeated_ids} records whether a person
#'   contributes more than one informative calibration row;
#'   \code{repeated_residual_ids} records repetition among rows contributing
#'   fitted residuals, which governs the row-based fit references. If estimation
#'   does not converge, locations and
#'   residual patterns are retained for diagnosis, but standard errors,
#'   separation indices and inferential probabilities are \code{NA}.
#' @references
#' Rasch, G. (1960). Probabilistic Models for Some Intelligence and
#' Attainment Tests. Copenhagen: Danish Institute for Educational Research.
#' (Expanded edition, 1980, Chicago: University of Chicago Press.)
#'
#' Rasch, G. (1961). On general laws and the meaning of measurement in
#' psychology. In Proceedings of the Fourth Berkeley Symposium on
#' Mathematical Statistics and Probability (Vol. 4, pp. 321--333).
#' Berkeley: University of California Press.
#'
#' Andrich, D. and Luo, G. (2003). Conditional pairwise estimation in the
#' Rasch model for ordered response categories using principal components.
#' Journal of Applied Measurement, 4(3), 205--221.
#'
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. Psychometrika, 54(3), 427--450.
#' @seealso \code{\link{rasch_mfrm}}, \code{\link{rasch_efrm}},
#'   \code{\link{btl}}, \code{\link{dif_anova}},
#'   \code{\link{test_information}}, and \code{\link{run_app}}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(X, model = "PCM")
#' fit$items
#' fit$psi$PSI
#' @export
rasch <- function(data, model = c("PCM", "RSM"), id = NULL, factors = NULL,
                  items = NULL, n_groups = NULL, anchors = NULL,
                  na_codes = -1, key = NULL, pc_components = NULL,
                  maxit = 60, tol = 1e-8) {
  .check_column_names(data)
  n_groups_requested <- n_groups
  # name for a factors= vector passed by value (not by column name)
  .factors_sym <- substitute(factors)
  .factors_label <- if (is.name(.factors_sym))
    as.character(.factors_sym) else "factor"
  model <- match.arg(model)
  if (!is.null(id) && (!is.atomic(id) || !is.null(dim(id))))
    stop("`id` must name one data column or be a plain vector with one value per row",
         call. = FALSE)
  if (!is.null(items) &&
      (!(is.character(items) || is.numeric(items)) || is.complex(items) ||
       !is.null(dim(items)) || !is.null(oldClass(items)) || !length(items) ||
       anyNA(items)))
    stop("`items` must be a non-empty plain vector of item names or indices",
         call. = FALSE)
  if (!is.null(factors) && !is.data.frame(factors) &&
      (!is.atomic(factors) || !is.null(dim(factors))))
    stop("`factors` must be a data frame, column names, or a plain vector with one value per row",
         call. = FALSE)
  if (!is.null(n_groups))
    n_groups <- .check_whole(n_groups, "n_groups", 2)
  .check_controls(maxit, tol)
  if (!is.null(pc_components)) {
    if (model != "PCM")
      stop("pc_components applies to the PCM only")
    if (!is.null(anchors))
      stop("pc_components cannot be combined with anchors")
  }
  if (!is.null(anchors)) {
    if (!is.data.frame(anchors))
      stop("anchors must be a data frame with columns item, k, tau")
    .check_column_names(anchors)
    missing_anchor_columns <- setdiff(c("item", "k", "tau"), names(anchors))
    if (length(missing_anchor_columns))
      stop("anchors must contain columns item, k, tau; missing: ",
           paste(missing_anchor_columns, collapse = ", "))
    if (!nrow(anchors))
      stop("anchors must contain at least one threshold or item-location row")
    anchors <- .pcml_anchor_columns(anchors)
  }

  # --- split data frame into ID, factors, and item columns ---------------
  id_vec <- NULL; fac_df <- NULL
  # a simulated dataset carries its own person identifier: use it, so the
  # documented bare call rasch(simulate_rasch(...)) keeps person ids
  if (inherits(data, "rasch_sim") && is.null(id) && "id" %in% names(data))
    id <- "id"
  if (is.data.frame(data)) {
    nm <- names(data)
    id_is_col <- is.character(id) && length(id) == 1L
    # a misspelled column name must be an error, never a silent fallback:
    # dropping it quietly produces a valid-looking analysis of the wrong data
    if (id_is_col) {
      if (!id %in% nm)
        stop("id column '", id, "' not found in the data")
      id_vec <- data[[id]]
    } else if (!is.null(id)) {
      # a supplied id vector must line up with the data; a length mismatch
      # (e.g. a stale upstream vector) must error, not be silently dropped
      if (length(id) != nrow(data))
        stop("`id` has ", length(id), " entries but the data has ",
             nrow(data), " rows")
      id_vec <- id
    }
    # Character input is ambiguous: all existing column names means the
    # documented column-name form; a row-length character vector is a
    # grouping vector passed by value. A short non-matching character input
    # remains a misspelled-column error rather than silently changing modes.
    factors_are_cols <- .role_columns(factors, nm, nrow(data))
    factors_by_value <- !is.null(factors) && is.atomic(factors) &&
      !factors_are_cols
    if (factors_are_cols) {
      miss <- setdiff(factors, nm)
      if (length(miss))
        stop("factor column(s) not found in the data: ",
             paste(miss, collapse = ", "))
      if (anyDuplicated(factors))
        stop("factor column(s) named more than once: ",
             paste(unique(factors[duplicated(factors)]), collapse = ", "))
      fac_df <- data[, factors, drop = FALSE]
    } else if (is.data.frame(factors)) {
      if (nrow(factors) != nrow(data))
        stop("`factors` data frame has ", nrow(factors), " rows but the data ",
             "has ", nrow(data), " rows")
      if (anyDuplicated(names(factors)))
        stop("duplicate factor column name(s): ",
             paste(unique(names(factors)[duplicated(names(factors))]),
                   collapse = ", "))
      clash <- intersect(names(factors), nm)
      different <- clash[!vapply(clash, function(cn)
        .same_role_values(factors[[cn]], data[[cn]]), logical(1))]
      if (length(different) && is.null(items))
        stop("external factor column(s) share item-data names but contain ",
             "different values: ", paste(different, collapse = ", "),
             ". Rename the external factor column(s), or name the item ",
             "columns explicitly with items=", call. = FALSE)
      fac_df <- factors
    } else if (factors_by_value) {
      # a factors= grouping vector passed by VALUE (not by column name):
      # accept it when it lines up, rather than silently ignoring it and
      # leaving any same-named data column to be treated as an item
      if (!is.atomic(factors) || length(factors) != nrow(data))
        stop("`factors` must be column name(s) in the data, a data frame ",
             "with one row per data row, or a vector with one entry per row")
      fac_df <- stats::setNames(data.frame(factors, stringsAsFactors = FALSE),
                                .factors_label)
    } else if (!is.null(factors)) {
      stop("`factors` must be column name(s) in the data, a data frame ",
           "with one row per data row, or a vector with one entry per row")
    }
    # a data column whose values are identical to a by-value factors vector
    # is almost certainly that same variable: exclude it so it is not also
    # scored as a numeric item
    val_factor_cols <- if (factors_by_value)
      nm[vapply(data, function(col)
        length(col) == length(factors) &&
          .same_role_values(col, factors), logical(1))] else NULL
    val_id_cols <- if (!is.null(id) && !id_is_col)
      nm[vapply(data, function(col)
        length(col) == length(id) &&
          .same_role_values(col, id), logical(1))] else NULL
    # a data column identical to a by-value role vector may be that same
    # variable, or a genuine item whose responses happen to agree. Deciding
    # silently risks the wrong analysis either way -- a real item would
    # vanish from the fit -- so an ambiguous match is refused unless items=
    # states which columns are items (the rule rasch_efrm() applies)
    val_matched <- unique(c(val_id_cols, val_factor_cols))
    if (is.null(items) && length(val_matched))
      stop("data column(s) identical to a supplied role vector: ",
           paste(val_matched, collapse = ", "),
           ". If they are the same variable, drop them from the data or ",
           "name the item columns with items=; a genuine item identical ",
           "to a role must be listed in items=")
    drop_cols <- c(if (id_is_col) id else NULL,
                   if (factors_are_cols) factors else NULL,
                   # Without an explicit item selector, an externally
                   # supplied factor data frame whose column names also
                   # appear in `data` almost certainly refers to those
                   # columns.  With `items=`, the separate factor frame is
                   # unambiguous and the selected data columns are genuine
                   # items; selected in-data roles remain guarded above.
                   if (is.data.frame(factors) && is.null(items))
                     intersect(names(factors), nm) else NULL)
    # identifier-named columns must never be silently SCORED as items: the
    # stacked/racked reshapes emit id/row_id/time columns, and calling
    # rasch(stacked) without id = "id" would otherwise rescore a numeric
    # person identifier as a many-category item with a valid-looking
    # report. Only numeric-convertible columns can be scored, so only they
    # are refused; character identifiers keep the old dropped-with-a-note
    # path (they can never silently enter the item matrix).
    ident_like <- intersect(c("id", "row_id", "time", "person"), nm)
    ident_like <- setdiff(ident_like, c(drop_cols,
                                        if (is.character(items)) items))
    ident_like <- ident_like[vapply(ident_like, function(cn) {
      v <- data[[cn]]
      vn <- suppressWarnings(as.numeric(as.character(v)))
      any(!is.na(vn))
    }, logical(1))]
    if (is.null(items) && length(ident_like))
      stop("the data contain identifier-like column(s) not assigned a role: ",
           paste(ident_like, collapse = ", "),
           " -- pass them via id=/factors= and name the item columns with ",
           "items= (for stack_data output: rasch(stacked, id = \"id\", ",
           "factors = \"time\", items = <the item columns>)), or drop them")
    item_cols <- if (is.null(items)) setdiff(nm, drop_cols)
    else if (is.character(items)) {
      miss <- setdiff(items, nm)
      if (length(miss))
        stop("item column(s) not found in the data: ",
             paste(miss, collapse = ", "))
      items
    } else {
      if (!is.numeric(items) || any(!is.finite(items)) ||
          any(items != floor(items)) || any(items < 1) ||
          any(items > length(nm)))
        stop("numeric `items` indices must be whole numbers between 1 and ",
             length(nm))
      nm[as.integer(items)]
    }
    # an explicit items= must not silently pull in an id/factor column (a
    # positional items = 1:k over an id-first layout would score the id as
    # an item and drop a real one), nor name the same column twice
    if (!is.null(items)) {
      clash <- intersect(item_cols, drop_cols)
      if (length(clash))
        stop("items= includes id/factor column(s): ",
             paste(clash, collapse = ", "),
             " -- name only item columns, or drop them from id=/factors=")
    }
    dup <- item_cols[duplicated(item_cols)]
    if (length(dup))
      stop("item column(s) named more than once: ",
           paste(unique(dup), collapse = ", "))
    X <- as.matrix(data[, item_cols, drop = FALSE])
  } else {
    X <- as.matrix(data)
    # matrix input: items= selects columns exactly as it does for a data
    # frame; ignoring it silently would score every column while the caller
    # believes a subset was fitted
    if (!is.null(items)) {
      if (is.character(items)) {
        miss <- setdiff(items, colnames(X))
        if (length(miss))
          stop("item column(s) not found in the data: ",
               paste(miss, collapse = ", "))
        X <- X[, items, drop = FALSE]
      } else {
        if (!is.numeric(items) || any(!is.finite(items)) ||
            any(items != floor(items)) || any(items < 1) ||
            any(items > ncol(X)))
          stop("numeric `items` indices must be whole numbers between 1 and ",
               ncol(X))
        X <- X[, as.integer(items), drop = FALSE]
      }
    }
    if (!is.null(id)) {
      if (length(id) != nrow(X))
        stop("`id` has ", length(id), " entries but the data has ",
             nrow(X), " rows")
      id_vec <- id
    }
    if (is.data.frame(factors)) {
      if (nrow(factors) != nrow(X))
        stop("`factors` data frame has ", nrow(factors), " rows but the data ",
             "has ", nrow(X), " rows")
      fac_df <- factors
    } else if (!is.null(factors)) {
      if (!is.atomic(factors) || length(factors) != nrow(X))
        stop("`factors` must be a data frame with one row per data row, or ",
             "a vector with one entry per row")
      fac_df <- stats::setNames(data.frame(factors, stringsAsFactors = FALSE),
                                .factors_label)
    }
  }
  if (is.null(id_vec)) id_vec <- seq_len(nrow(X))
  id_vec <- .canonical_role_column(id_vec)
  if (!is.null(fac_df)) fac_df[] <- lapply(fac_df, .canonical_role_column)

  # score multiple-choice items against the key, keeping the raw responses
  mc <- NULL
  if (!is.null(key)) {
    key <- .resolve_key(key)
    sc <- .score_mc(X, key, na_codes = na_codes)
    X[, colnames(sc$scored)] <- sc$scored
    mc <- list(key = sc$key, map = sc$map, raw = sc$raw)
  }

  items_before_prep <- colnames(X)
  if (!is.null(anchors)) {
    # a numeric anchor index means the caller's column, so it must resolve
    # against the data as supplied: preparation can drop a constant item,
    # and resolving against the surviving columns would silently anchor a
    # different item. The names are resolved before preparation so that
    # the category merge can leave a fully anchored item alone
    a_names <- if (is.character(anchors$item) || is.factor(anchors$item))
      as.character(anchors$item)
    else {
      ai <- anchors$item
      if (!is.numeric(ai) || any(!is.finite(ai)) || any(ai != floor(ai)) ||
          any(ai < 1) || any(ai > length(items_before_prep)))
        stop("numeric anchor item indices must be whole numbers between 1 and ",
             length(items_before_prep))
      items_before_prep[as.integer(ai)]
    }
    if (anyDuplicated(paste(a_names, anchors$k)))
      stop("duplicate anchor for the same item threshold: ",
           paste(unique(a_names[duplicated(paste(a_names, anchors$k))]),
                 collapse = ", "))
    anchors$item <- a_names
  }
  prep <- .prepare_X(X, na_codes = na_codes, model = model, anchors = anchors,
                     scored_items = if (is.null(mc)) character(0) else
                       colnames(mc$raw))
  X <- prep$X
  if (!is.null(mc)) {
    binary_key <- vapply(mc$map, function(z) max(z) == 1L, logical(1))
    if (any(binary_key))
      prep$notes <- c(prep$notes, sprintf(
        "%d item(s) scored 0/1 against the key", sum(binary_key)))
    if (any(!binary_key))
      prep$notes <- c(prep$notes, sprintf(
        "%d item(s) scored with polytomous option-score maps",
        sum(!binary_key)))
    gone <- setdiff(colnames(mc$raw), colnames(X))
    if (length(gone)) mc$raw <- mc$raw[, setdiff(colnames(mc$raw), gone), drop = FALSE]
  }

  if (!is.null(anchors)) {
    gone <- setdiff(a_names, colnames(X))
    if (length(gone))
      stop("anchored item(s) not present after data preparation: ",
           paste(gone, collapse = ", "))
    # match the note's own "item <name> rescored" prefix: a bare name search
    # would let anchor I1 trip over a note about I10
    resc <- vapply(prep$notes, function(n)
      any(vapply(paste0("item ", a_names, " rescored"), grepl, TRUE,
                 x = n, fixed = TRUE)), TRUE)
    if (any(resc))
      stop("anchored item(s) were rescored during data preparation; ",
           "anchor values would no longer match the threshold numbering")
    prep$notes <- c(prep$notes,
                    if ("average" %in% names(anchors) && all(anchors$average %in% TRUE))
                      sprintf("scale origin from the average location of %d anchor item(s)",
                              nrow(anchors))
                    else sprintf("%d anchor constraint(s); scale origin from anchors",
                                 nrow(anchors)))
  }

  # --- item estimation ----------------------------------------------------
  est <- if (is.null(pc_components))
    .pcml_fit(X, model = model, anchors = anchors, maxit = maxit, tol = tol,
              cluster = id_vec)
  else .pcml_pc_fit(X, n_components = pc_components, maxit = maxit,
                    tol = tol, cluster = id_vec)
  repeated_ids <- .est_has_repeated_calibration_units(est, id_vec)
  if (repeated_ids)
    prep$notes <- c(prep$notes, paste(
      "item-parameter sandwich covariance clustered by repeated person id;",
      "point estimates retain every response row"))
  if (!is.null(pc_components))
    prep$notes <- c(prep$notes,
                    sprintf("thresholds estimated through %d principal component(s); see est$components",
                            pc_components))
  if (!isTRUE(est$converged))
    warning("estimation did NOT converge in ", est$iterations,
            " iterations: estimates, standard errors, fit statistics, and ",
            "p-values are unreliable -- increase maxit or check the data ",
            "for unanswerable structure", call. = FALSE)
  .check_factor_frame(fac_df)
  fit <- .assemble_fit(model, X, est, id_vec, fac_df, n_groups,
                       c(prep$notes, est$notes))
  fit$mc <- mc
  # Keep the arguments that define the fitted model. Post-fit operations such
  # as drop_items() must not silently change the identification, threshold
  # parameterisation, fit grouping, or optimiser controls when they refit.
  # Anchors are stored by item name so their meaning survives column removal.
  anchors_named <- anchors
  if (!is.null(anchors_named) &&
      !(is.character(anchors_named$item) || is.factor(anchors_named$item)))
    anchors_named$item <- colnames(X)[as.integer(anchors_named$item)]
  key_spec <- NULL
  if (!is.null(mc)) {
    key_spec <- do.call(rbind, lapply(names(mc$map), function(it) {
      data.frame(item = it, option = names(mc$map[[it]]),
                 score = unname(mc$map[[it]]), stringsAsFactors = FALSE)
    }))
    rownames(key_spec) <- NULL
  }
  fit$refit_spec <- list(
    model = model, n_groups = n_groups_requested,
    anchors = anchors_named, na_codes = na_codes, key = key_spec,
    pc_components = pc_components, maxit = maxit, tol = tol)
  fit
}

# Post-estimation pipeline shared by rasch(), rasch_mfrm(), and rasch_efrm():
# person estimation, residuals, the full fit suite, and the assembled tables.
# disc is an optional per-column discrimination (frame unit) vector; with
# unequal discriminations the raw score is no longer sufficient, so person
# estimation switches to the weighted-score routine and the score table is
# replaced by per-unit score curves.
.has_repeated_person_ids <- function(id) {
  if (is.null(id)) return(FALSE)
  z <- .role_text_values(id)
  present <- !is.na(z) & nzchar(z)
  anyDuplicated(z[present]) > 0L
}

# Row-based diagnostics have a different sampling-unit requirement from the
# calibration sandwich. A row can carry a fitted residual even when it has no
# informative conditional item pair, so it must still count when deciding
# whether a row-independent fit reference or null generator is defensible.
.has_repeated_residual_units <- function(fit) {
  if (!is.list(fit) || is.null(fit$person$id) ||
      !is.matrix(fit$residuals) ||
      nrow(fit$residuals) != length(fit$person$id))
    return(.has_repeated_person_ids(
      if (is.list(fit) && !is.null(fit$person)) fit$person$id else NULL))
  contributes <- rowSums(is.finite(fit$residuals)) > 0L
  z <- .role_text_values(fit$person$id)
  present <- contributes & !is.na(z) & nzchar(z)
  anyDuplicated(z[present]) > 0L
}

# Whether more than one informative calibration row belongs to a person.
# Current PCML fits record this after excluding rows with no informative
# conditional item pair. Raw IDs are used only for legacy or structural fits
# that do not carry this support record.
.est_has_repeated_calibration_units <- function(est, id = NULL) {
  support <- est$cluster_support
  if (is.list(support) && is.logical(support$repeated) &&
      length(support$repeated) == 1L && !is.na(support$repeated))
    return(isTRUE(support$repeated))
  .has_repeated_person_ids(id)
}

.has_repeated_calibration_units <- function(fit) {
  if (!is.list(fit)) return(FALSE)
  id <- if (!is.null(fit$person)) fit$person$id else NULL
  .est_has_repeated_calibration_units(fit$est %||% list(), id)
}

.assemble_fit <- function(model, X, est, id_vec, fac_df, n_groups,
                          notes, disc = NULL) {
  m <- est$m; L <- ncol(X)
  repeated_ids <- .est_has_repeated_calibration_units(est, id_vec)
  thr <- est$thr
  tau_list <- lapply(seq_len(L), function(i) thr$tau[thr$item == i])
  names(tau_list) <- colnames(X)
  equal_disc <- is.null(disc) || length(unique(disc)) == 1L
  disc_v <- if (is.null(disc)) rep(1, L) else disc

  # --- person estimation and residuals -------------------------------------
  person <- if (equal_disc) .person_estimates(X, tau_list, disc = disc_v[1])
            else .efrm_person_estimates(X, tau_list, disc_v)
  mo <- .moment_arrays(person$theta, tau_list, disc = disc_v)
  Z <- (X - mo$E) / sqrt(mo$V)
  colnames(Z) <- colnames(X)
  repeated_residual_ids <- {
    z <- .role_text_values(id_vec)
    present <- rowSums(is.finite(Z)) > 0L & !is.na(z) & nzchar(z)
    anyDuplicated(z[present]) > 0L
  }

  # --- fit statistics ------------------------------------------------------
  # an item scored at its floor or ceiling by every non-extreme person has
  # no finite location; exclude it from person fit as extreme persons are
  # excluded from item fit
  m_i <- if (!is.null(est$m)) est$m else apply(X, 2, max, na.rm = TRUE)
  item_extreme <- vapply(seq_len(ncol(X)), function(j) {
    col <- X[!person$extreme, j]
    tot <- sum(col, na.rm = TRUE); nn <- sum(!is.na(col))
    nn == 0 || tot == 0 || tot == nn * m_i[j]
  }, logical(1))
  ifit <- .item_fit(X, Z, mo, disc = if (is.null(disc)) NULL else disc_v,
                    extreme = person$extreme)
  pfit <- .person_fit(X, Z, mo, disc = if (is.null(disc)) NULL else disc_v,
                      item_extreme = item_extreme)
  n_par <- if (is.null(est$n_parameters)) nrow(est$thr) - 1L else est$n_parameters
  rf <- .fitres(Z, mo, person$extreme, n_par)
  ng_req <- n_groups                       # NULL = the automatic rule
  ci <- .class_intervals(person$theta, person$extreme, n_groups)
  n_groups <- attr(ci, "n_groups")
  # class intervals compiled per item when data are missing (the automatic
  # per-item basis; Andrich & Marais 2019, ch. 15), so every item is tested
  # over intervals of
  # its own responders (with the group-count rule applied per item when
  # automatic)
  ci_list <- if (anyNA(X))
    .class_intervals_by_item(X, person$theta, person$extreme, ng_req)
  else NULL
  it <- .item_trait(X, mo, ci, ci_list = ci_list)
  ia <- .item_anova(Z, ci, person$extreme, ci_list = ci_list)
  if (repeated_residual_ids) {
    for (nm in intersect(c("p", "p_adj", "p_bonf"), names(it)))
      it[[nm]][] <- NA_real_
    for (nm in intersect(c("p", "p_adj", "p_bonf"), names(ia)))
      ia[[nm]][] <- NA_real_
    notes <- c(notes, paste(
      "item-trait and class-interval ANOVA probabilities withheld because",
      "their row-based references do not model within-person dependence;",
      "fit statistics remain descriptive"))
  }
  psi <- .psi(person$theta, person$se)
  psi_noext <- .psi(person$theta, person$se, keep = !person$extreme)
  alpha <- .alpha(X)

  # --- assembled person table ----------------------------------------------
  parts <- list(data.frame(id = id_vec), fac_df, person,
                data.frame(infit_ms = pfit$infit_ms, outfit_ms = pfit$outfit_ms,
                           infit_z = pfit$infit_z, outfit_z = pfit$outfit_z,
                           fit_resid = rf$persons$fit_resid,
                           natural_resid = rf$persons$natural,
                           df_fit = rf$persons$df, class_interval = ci))
  person <- do.call(cbind, parts[!vapply(parts, is.null, TRUE)])
  rownames(person) <- NULL

  # --- item table ----------------------------------------------------------
  loc <- vapply(tau_list, mean, 0)
  weak_thr <- if (is.null(thr$weak)) rep(FALSE, nrow(thr)) else thr$weak
  se_loc <- vapply(seq_len(L), function(i) {
    rows <- thr$id[thr$item == i]
    anchored_i <- thr$anchored[thr$item == i]
    # A location anchor fixes the mean of an item's thresholds exactly even
    # when a sparse category makes their individual spread estimates weak.
    # Preserve the exact zero variance of that mean; only the free threshold
    # SEs are unavailable.  Checking weakness first would turn a known anchor
    # into NA and remove it from fixed-origin equating downstream.
    if (length(anchored_i) && all(anchored_i)) {
      vloc <- mean(est$cov_tau[rows, rows])
      if (is.finite(vloc)) return(sqrt(max(vloc, 0)))
      # When repeated-person support is insufficient, the estimated
      # covariance is deliberately unavailable. A location anchor nevertheless
      # fixes this particular mean exactly; do not make the exact constraint
      # look uncertain merely because its free threshold spread has no SE.
      if (isTRUE(est$converged) && identical(est$cluster_inference, FALSE))
        return(0)
      return(NA_real_)
    }
    # a weakly determined threshold (sparse adjacent category) makes the
    # ridged covariance block spuriously small: report NA, not a number
    if (any(weak_thr[thr$item == i])) return(NA_real_)
    # anchored items have a structurally zero variance that floating-point
    # noise can render as a tiny negative number on some BLAS builds
    sqrt(max(mean(est$cov_tau[rows, rows]), 0))
  }, 0)
  items_df <- data.frame(item = colnames(X), max = m, location = loc,
                         se = se_loc,
                         disc = .item_discrim(person$theta, X, tau_list,
                                              person$extreme),
                         fit_resid = rf$items$fit_resid, df_fit = rf$items$df,
                         natural_resid = rf$items$natural,
                         infit_ms = ifit$infit_ms, outfit_ms = ifit$outfit_ms,
                         infit_z = ifit$infit_z, outfit_z = ifit$outfit_z,
                         chisq = it$chisq, df = it$df, p = it$p,
                         p_adj = it$p_adj, p_bonf = it$p_bonf,
                         F_anova = ia$F_anova, p_anova = ia$p,
                         p_anova_adj = ia$p_adj,
                         p_anova_bonf = ia$p_bonf)
  rownames(items_df) <- NULL

  # --- score table (complete responders; raw score is only sufficient when
  # --- discriminations are equal) ---------------------------------------------
  sc <- if (equal_disc) {
    pe_full <- person_wle(tau_list, disc = disc_v[1])
    data.frame(score = 0:sum(m), theta = unname(pe_full$theta),
               se = unname(pe_full$se))
  } else NULL
  if (is.null(sc))
    notes <- c(notes, "person measures use the weighted score; per-group score curves replace the raw-score table (see score_curves)")

  # --- threshold diagnostics --------------------------------------------------
  td <- lapply(seq_len(L), function(i) {
    tau_i <- tau_list[[i]]
    modal <- .modal_score_categories(tau_i)
    list(item = colnames(X)[i], thresholds = tau_i,
         ordered = all(diff(tau_i) > 0) || length(tau_i) == 1L,
         reversed_at = which(diff(tau_i) <= 0) + 1L,
         never_modal_categories = setdiff(0:length(tau_i), modal),
         category_counts = as.integer(table(factor(X[, i], levels = 0:length(tau_i)))))
  })
  names(td) <- colnames(X)

  separation_quality <- .separation_quality(psi$PSI)
  total_ok <- is.finite(it$chisq) & is.finite(it$df) & it$df > 0
  total_chisq <- if (any(total_ok)) sum(it$chisq[total_ok]) else NA_real_
  total_df <- if (any(total_ok)) sum(it$df[total_ok]) else NA_integer_
  if (!isTRUE(est$converged)) {
    # Keep locations and residual patterns so a stalled optimisation can be
    # diagnosed, but do not attach uncertainty or hypothesis tests to its
    # last numerical iterate.
    thr$se[] <- NA_real_
    person$se[] <- NA_real_
    if (!is.null(sc)) sc$se[] <- NA_real_
    items_df$se[] <- NA_real_
    for (nm in intersect(c("p", "p_adj", "p_bonf", "p_anova",
                           "p_anova_adj", "p_anova_bonf"), names(items_df)))
      items_df[[nm]][] <- NA_real_
    for (nm in intersect(c("p", "p_adj", "p_bonf"), names(it)))
      it[[nm]][] <- NA_real_
    for (nm in intersect(c("p", "p_adj", "p_bonf"), names(ia)))
      ia[[nm]][] <- NA_real_
    psi <- .psi(person$theta, person$se)
    psi_noext <- .psi(person$theta, person$se, keep = !person$extreme)
    separation_quality <- "unknown"
  }
  estimated_item <- vapply(seq_len(L), function(i) {
    ai <- thr$anchored[thr$item == i]
    !length(ai) || any(!ai)
  }, logical(1))
  out <- list(model = model, X = X, m = m, items = items_df, thresholds = thr,
              tau_list = tau_list, person = person, score_table = sc,
              person_scoring_algorithm = "max-wle-1",
              residuals = Z, moments = mo, n_groups = n_groups,
              ci_item = ci_list,
              item_trait = it, item_anova = ia,
              psi = psi, psi_noext = psi_noext,
              isi = .psi(items_df$location, items_df$se,
                         keep = estimated_item),
              alpha = alpha,
              targeting = .targeting(person, thr),
              repeated_ids = repeated_ids,
              repeated_residual_ids = repeated_residual_ids,
              separation_quality = separation_quality,
              power_of_fit = separation_quality,
              total_chisq = total_chisq, total_df = total_df,
              total_chisq_p = if (isTRUE(est$converged) &&
                                    !repeated_residual_ids &&
                                    is.finite(total_chisq))
                pchisq(total_chisq, total_df, lower.tail = FALSE) else NA_real_,
              item_fit_summary = .dist_stats(rf$items$fit_resid),
              person_fit_summary = .dist_stats(rf$persons$fit_resid),
              summary_stats = list(
                item_location = .dist_stats(items_df$location),
                person_location = .dist_stats(person$theta),
                person_location_noext = .dist_stats(person$theta[!person$extreme]),
                cor_item_fit_location = .safe_cor(items_df$location,
                                                  rf$items$fit_resid),
                cor_person_fit_location = .safe_cor(person$theta,
                                                    rf$persons$fit_resid),
                df_factor = rf$f_cell),
              thresholds_diag = td, est = est, notes = notes,
              factors = fac_df, disc = disc)
  out <- .tag_tables(out)
  class(out) <- "rasch"
  out
}

#' @export
print.rasch <- function(x, ...) {
  separation_quality <- x$separation_quality %||% x$power_of_fit %||%
    .separation_quality(x$psi$PSI)
  cat(sprintf("rasch %s analysis: %d items, %d persons\n",
              x$model, ncol(x$X), nrow(x$X)))
  cat(sprintf("Pairwise conditional ML (%s): %s in %d iterations\n",
              if (is.null(x$est$components)) "Zwinderman"
              else "Andrich & Luo principal components",
              if (x$est$converged) "converged" else "NOT converged",
              x$est$iterations))
  cat(sprintf("PSI %.3f (no extremes %.3f), item SI %.3f, alpha %.3f%s, separation quality: %s\n",
              x$psi$PSI, x$psi_noext$PSI, x$isi$PSI, x$alpha$alpha,
              if (isFALSE(x$alpha$applicable))
                sprintf(" [complete cases only, n = %d]", x$alpha$n) else "",
              separation_quality))
  total_p <- if (.has_repeated_residual_units(x)) NA_real_
    else x$total_chisq_p
  fit_probability <- if (is.finite(total_p))
    paste0("p = ", .fmt_p(total_p)) else "probability unavailable"
  cat(sprintf(paste0("Approximate asymptotic total item-trait chi-square ",
                     "%.3f on %d df, %s\n"),
              x$total_chisq, x$total_df, fit_probability))
  if (length(x$notes)) cat(sprintf("Notes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

#' @export
summary.rasch <- function(object, ...) {
  x <- object
  structural <- inherits(x, c("rasch_mfrm", "rasch_efrm"))
  unit <- if (structural) "Response-cell" else "Item"
  units <- if (structural) "response cells" else "items"
  print(x)
  cat(sprintf("\nTargeting: person mean %.3f (SD %.3f); %sthresholds span %.3f to %.3f\n",
              x$targeting$person_mean, x$targeting$person_sd,
              if (structural) "calibration " else "",
              x$targeting$threshold_range[1], x$targeting$threshold_range[2]))
  cat(sprintf("%s fit residual mean %.3f SD %.3f (skew %.2f, kurt %.2f); person fit residual mean %.3f SD %.3f (skew %.2f, kurt %.2f)\n",
              unit,
              x$item_fit_summary$mean, x$item_fit_summary$sd,
              x$item_fit_summary$skewness, x$item_fit_summary$kurtosis,
              x$person_fit_summary$mean, x$person_fit_summary$sd,
              x$person_fit_summary$skewness, x$person_fit_summary$kurtosis))
  cat(sprintf("Fit residual-location correlation: %s %.3f, persons %.3f; cell df factor %.3f\n",
              units,
              x$summary_stats$cor_item_fit_location,
              x$summary_stats$cor_person_fit_location,
              x$summary_stats$df_factor))
  item_p <- if (.has_repeated_residual_units(x))
    rep(NA_real_, nrow(x$items)) else x$items$p_adj
  inference <- .inference_count(item_p)
  cat(sprintf("%s with approximate asymptotic Holm p < 0.05: %s\n\n",
              if (structural) "Response cells" else "Items",
              inference$text))
  core <- c("item", "max", "location", "se", "fit_resid", "infit_ms",
            "outfit_ms", "chisq", "df", "p_adj")
  print(.fmt_df(x$items[, intersect(core, names(x$items))]), row.names = FALSE)
  cat("(further columns on fit$items: natural and standardised forms,\n",
      " ANOVA fit, Bonferroni probabilities)\n", sep = "")
  dis <- vapply(x$thresholds_diag, function(d) !d$ordered, TRUE) &
    vapply(x$thresholds_diag, function(d) length(d$thresholds) > 1L, TRUE)
  if (any(dis)) cat(sprintf("\nDisordered %sthresholds: %s\n",
                            if (structural) "response-cell " else "",
                            paste(names(dis)[dis], collapse = ", ")))
  invisible(x)
}
