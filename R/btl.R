# rasch :: Bradley-Terry-Luce paired comparisons
# ===========================================================================
# The Bradley-Terry-Luce model (Bradley & Terry 1952; Luce 1959) for paired
# comparisons: Pr{a beats b} = exp(beta_a - beta_b) / (1 + exp(...)). This
# is the conditional form of the dichotomous Rasch model (Rasch 1961;
# Andrich 1978): given that one of two items is answered correctly, the
# probability the harder one was is exactly of BTL form, and the package's
# pairwise conditional estimation maximises precisely such a likelihood on
# the pair-conditional counts. BTL therefore sits inside the same family
# with the person parameter replaced by exchangeable comparison
# replications, and it is estimated here by the same conventions as
# everything else: Newton-Raphson on the conditional likelihood, sum-zero
# identification, and Godambe sandwich standard errors, clustered by judge
# when judges are identified. Fit follows the same residual logic as the
# rest of the package (Andrich & Marais 2019, ch. 23): per
# comparison z = (y - P)/sqrt(PQ); objects and judges carry the
# log-of-mean-square fit residual over their comparisons with apportioned
# degrees of freedom, and the classical pairwise chi-square compares the
# observed and expected win proportions of every pair.
# ===========================================================================

# Exposure and carry-over covariates from each judge's own history: for
# comparison r, exposure is 1(judge saw object_a before) - 1(saw object_b
# before); carry-over differences the judge's mean prior verdicts on the
# two objects (oriented to each object, scaled to [-1, 1], zero when
# unseen). Both enter the exponent like the location difference does, so
# the dependence is measured in logits (Davidson & Beaver 1977 order-effect
# device; response-dependence logic of Marais & Andrich 2008).
.btl_exposure <- function(a, b, x, m, jd, ord, w = rep(1, length(a))) {
  # the history is built by walking each judge's comparisons in sequence,
  # so the sequence has to order them. Tied values are broken by row order,
  # and the same data read in a different order would then carry different
  # exposure and carry-over covariates.
  dup <- duplicated(data.frame(judge = jd, ord = ord,
                               stringsAsFactors = FALSE))
  if (any(dup))
    stop("the judging sequence repeats within judge(s): ",
         paste(utils::head(unique(jd[dup]), 5), collapse = ", "),
         if (length(unique(jd[dup])) > 5) ", ..." else "",
         "; a repeated value leaves the order of those comparisons to the ",
         "row order of the data", call. = FALSE)
  R <- length(a)
  Fa <- Fb <- Wa <- Wb <- numeric(R)
  cnt <- new.env(hash = TRUE, parent = emptyenv())
  tot <- new.env(hash = TRUE, parent = emptyenv())
  gets <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
  key_a <- .factor_keys(data.frame(judge = jd, object = a,
                                   stringsAsFactors = FALSE))
  key_b <- .factor_keys(data.frame(judge = jd, object = b,
                                   stringsAsFactors = FALSE))
  for (r in order(jd, ord)) {
    ka <- key_a[r]; kb <- key_b[r]
    na_ <- gets(cnt, ka); nb_ <- gets(cnt, kb)
    Fa[r] <- as.numeric(na_ > 0); Fb[r] <- as.numeric(nb_ > 0)
    if (na_ > 0) Wa[r] <- gets(tot, ka) / na_
    if (nb_ > 0) Wb[r] <- gets(tot, kb) / nb_
    # count-weighted rows stand for w identical comparisons, so they enter
    # the history with weight w, consistent with the weighted likelihood
    cnt[[ka]] <- na_ + w[r]; cnt[[kb]] <- nb_ + w[r]
    tot[[ka]] <- gets(tot, ka) + w[r] * (2 * x[r] / m - 1)
    tot[[kb]] <- gets(tot, kb) + w[r] * (2 * (m - x[r]) / m - 1)
  }
  cbind(exposure = Fa - Fb, carry_over = Wa - Wb)
}

# union-find connectivity over the comparison graph
.btl_components <- function(K, ia, ib) {
  parent <- seq_len(K)
  find <- function(x) { while (parent[x] != x) x <- parent[x]; x }
  for (r in seq_along(ia)) {
    ra <- find(ia[r]); rb <- find(ib[r])
    if (ra != rb) parent[ra] <- rb
  }
  vapply(seq_len(K), find, 1L)
}

#' Fit comparative judgement models to paired comparisons
#'
#' Fits the Bradley--Terry--Luce model to dichotomous comparisons or its
#' ordered-response extension (Tutz 1986). Object, judge, and pair fit are
#' reported with an object separation index and design diagnostics.
#'
#' @details
#' For objects \eqn{a} and \eqn{b}, the dichotomous model is
#' \deqn{P(a\succ b)=\frac{\exp(\beta_a)}
#' {\exp(\beta_a)+\exp(\beta_b)}.}
#' This is the conditional form of the dichotomous Rasch model (Andrich 1978).
#' For an ordered response \eqn{Y=0,\ldots,m},
#' \deqn{\log\{P(Y=r)/P(Y=r-1)\}=\beta_a-\beta_b-\tau_r,}
#' with thresholds constrained to be symmetric under reversal of presentation
#' order. Two categories reproduce the dichotomous model.
#'
#' Locations are identified by a sum-zero constraint unless anchors are
#' supplied. The comparison graph must be connected, and the directed win
#' graph must be strongly connected for all free locations to be finite
#' (Ford 1957). Boundary objects are removed when this leaves an identified
#' model; otherwise fitting stops.
#'
#' Standard errors use the Godambe sandwich covariance. Inference is withheld
#' if the empirical sampling-unit score covariance does not identify every
#' fitted direction. When \code{judge} is supplied, the covariance is
#' clustered by judge and inference is also withheld when there are fewer
#' than ten judges, fewer than eight effective judges, or no residual cluster
#' degrees of freedom. A caution is attached
#' when the effective count is below 9.5 or one judge supplies more than 20
#' per cent of the comparisons.
#' The pairwise chi-square remains descriptive when judges are identified
#' because its row-based chi-square reference does not model within-judge
#' dependence. \code{\link{fit_bootstrap}} supplies a parametric calibration
#' under the fitted response model.
#'
#' Dichotomous data may be supplied as a winner, with ties dropped or divided
#' equally between the two outcomes. Ordered data may instead be supplied
#' directly as scores from 0 to \eqn{m}, or assembled from the winner and an
#' ordered margin of victory. Plain factors are refused because alphabetical
#' ordering can reverse the response scale. The \code{"pc"} threshold option
#' retains the symmetric spread component, which can stabilise thin categories.
#'
#' If comparison order is supplied, exposure and carry-over effects are
#' estimated from each judge's preceding comparisons. The \code{position}
#' term estimates a first-presentation effect. These coefficients enter the
#' model jointly with the object locations and are reported in logits. They
#' are refused when the comparison design confounds them exactly with the
#' object-location contrasts. The carry-over estimate and clustered SE
#' remain descriptive below 30 judges;
#' its probability is withheld because null calibration is mildly
#' anti-conservative at smaller judge counts. Raw probabilities are retained,
#' but simultaneous decisions across the fitted dependence effects use Holm's
#' familywise adjustment in \code{dependence$p_adj}.
#' Anchors fix nominated object locations and replace the sum-zero origin.
#'
#' @param data A data frame with one comparison per row.
#' @param object_a,object_b Names of the columns holding the two objects
#'   compared. Columns used for the comparison roles must be distinct.
#' @param winner Name of the column holding the winner of each row: its
#'   value must equal one of the two objects. \code{"tie"} and \code{"draw"}
#'   mark ties. Do not supply both \code{winner} and \code{response}.
#' @param margin Optional ordered margin-of-victory column, combined with
#'   \code{winner} to construct an orientation-invariant response. Use an
#'   ordered factor (levels from the smallest to largest margin) or a positive
#'   numeric magnitude. Margins on ties and rows excluded from the analysis
#'   are ignored.
#' @param thresholds \code{"free"} (default) estimates every symmetric
#'   threshold; \code{"pc"} retains only the symmetric spread component.
#' @param response Optional ordered response favouring \code{object_a} over
#'   \code{object_b}: an ordered factor from least to greatest preference for
#'   \code{object_a}, or integer scores \code{0..m}.
#' @param judge Optional name of a judge column; enables the judge fit
#'   table and clusters the sandwich standard errors by judge.
#' @param order Optional column giving each judge's comparison sequence;
#'   requires \code{judge}. See Details. Incompatible with
#'   \code{ties = "half"}.
#' @param position If \code{TRUE}, estimate a first-presentation advantage,
#'   treating \code{object_a} as the first object in each comparison.
#' @param anchors Optional named numeric vector of fixed object locations.
#'   Anchored objects have standard error zero and must not be boundary objects.
#'   The values are treated as fixed: uncertainty from an earlier calibration
#'   is not included in the returned covariance or standard errors. Several
#'   anchors impose their stated relative spacing as well as the scale origin.
#' @param count Optional name of a column of replication counts (a row
#'   standing for several identical comparisons). Counts greater than one
#'   cannot be combined with \code{order}, because a compressed row does not
#'   retain the sequence of the comparisons it represents.
#' @param ties How to treat ties in the dichotomous analysis:
#'   \code{"drop"} (default, removed with a note), \code{"half"} (half a
#'   win each way, a common pragmatic device; the two halves remain one
#'   sampling unit in the sandwich because they are not independent
#'   Bernoulli trials), or
#'   \code{"error"}. With polytomous responses, code ties as a middle
#'   category instead.
#' @param maxit,tol Newton-Raphson iteration cap and convergence tolerance.
#' @param .object_design Internal object-location design used by
#'   \code{\link{btl_explanatory}}.
#' @return A \code{"rasch_btl"} object. Principal components are
#'   \code{objects}, \code{pairs}, \code{judges}, the total pair-fit test,
#'   \code{osi}, \code{loglik}, composite-likelihood information \code{cl},
#'   convergence details, and \code{notes}. Ordered-response fits also contain
#'   \code{thresholds}, \code{m}, and \code{categories}. Fits using
#'   \code{order} contain \code{dependence} and \code{dependence_data}; the
#'   former reports raw \code{p} and Holm-adjusted \code{p_adj}.
#'   A non-converged fit retains estimates and residual patterns for diagnosis
#'   but withholds standard errors, separation indices and probabilities.
#'   \code{comparisons} contains the rows used by the fitted likelihood;
#'   \code{observed_comparisons} retains all otherwise usable rows before
#'   free boundary objects are set aside, for observed-data descriptions such
#'   as transitivity.
#'   An undefeated or winless object is set aside from estimation, as an
#'   extreme person is in a Rasch calibration, and reported in
#'   \code{objects} with \code{extreme = TRUE} at an extrapolated location:
#'   the profile solution with its score moved half a point inside the
#'   boundary against the calibrated scale. Its standard error and fit are
#'   withheld and the row takes no part in inference or equating.
#' @references
#' Bradley, R. A. and Terry, M. E. (1952). Rank analysis of incomplete block
#' designs: I. The method of paired comparisons. Biometrika, 39, 324--345.
#'
#' Luce, R. D. (1959). Individual Choice Behavior. Wiley.
#'
#' Andrich, D. (1978). Relationships between the Thurstone and Rasch
#' approaches to item scaling. Applied Psychological Measurement, 2,
#' 451--462.
#'
#' Tutz, G. (1986). Bradley-Terry-Luce models with an ordered response.
#' Journal of Mathematical Psychology, 30(3), 306--316.
#'
#' Agresti, A. (1992). Analysis of ordinal paired comparison data. Journal of
#' the Royal Statistical Society C, 41(2), 287--297.
#'
#' Davidson, R. R. (1970). On extending the Bradley-Terry model to accommodate
#' ties in paired comparison experiments. Journal of the American Statistical
#' Association, 65(329), 317--328.
#'
#' Ford, L. R. (1957). Solution of a ranking problem from binary comparisons.
#' American Mathematical Monthly, 64(8), 28--33.
#'
#' Davidson, R. R. and Beaver, R. J. (1977). On extending the Bradley-Terry
#' model to incorporate within-pair order effects. Biometrics, 33(4),
#' 693--702.
#' @seealso \code{\link{btl_dif}}, \code{\link{btl_efrm}},
#'   \code{\link{btl_information}}, \code{\link{btl_transitivity}}, and
#'   \code{\link{simulate_btl}}.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pairs <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pairs[, 1], each = 30),
#'                 b = rep(pairs[, 2], each = 30))
#' p <- plogis(beta[d$a] - beta[d$b])
#' d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' btl(d, object_a = "a", object_b = "b", winner = "win")
#' @export
btl <- function(data, object_a, object_b, winner = NULL, response = NULL,
                margin = NULL, judge = NULL, count = NULL, order = NULL,
                position = FALSE, anchors = NULL,
                ties = c("drop", "half", "error"),
                thresholds = c("free", "pc"), maxit = 60, tol = 1e-8,
                .object_design = NULL) {
  .check_column_names(data)
  .check_controls(maxit, tol)
  if (!isTRUE(position) && !isFALSE(position))
    stop("`position` must be TRUE or FALSE")
  ties <- match.arg(ties)
  thresholds <- match.arg(thresholds)
  data <- as.data.frame(data)
  if (is.null(winner) && is.null(response))
    stop("give either `winner` (dichotomous) or `response` (polytomous)")
  # the polytomous branch reads `response` and never looks at `winner`, so
  # accepting both would fit one outcome while the caller states another
  if (!is.null(winner) && !is.null(response))
    stop("give either `winner` (dichotomous) or `response` (polytomous), ",
         "not both: the response would be fitted and the winner ignored")
  if (!is.null(margin) && is.null(winner))
    stop("`margin` requires `winner`")
  if (!is.null(order) && is.null(judge))
    stop("`order` requires `judge`: exposure is a within-judge history")
  for (nm in c("object_a", "object_b", "winner", "response", "margin",
               "judge", "count", "order")) {
    v <- get(nm, inherits = FALSE)
    if (!is.null(v)) .check_reshape_column(data, v, nm)
  }
  role_columns <- unlist(Filter(Negate(is.null),
                                list(object_a = object_a,
                                     object_b = object_b,
                                     winner = winner,
                                     response = response,
                                     margin = margin,
                                     judge = judge,
                                     count = count,
                                     order = order)), use.names = TRUE)
  repeated_roles <- unique(role_columns[duplicated(role_columns)])
  if (length(repeated_roles))
    stop("comparison role columns must be distinct; repeated: ",
         paste(repeated_roles, collapse = ", "))
  a <- .role_text_values(data[[object_a]])
  b <- .role_text_values(data[[object_b]])
  if (any(!is.na(a) & !nzchar(a)) || any(!is.na(b) & !nzchar(b)))
    stop("blank object identifier(s) in ", object_a, "/", object_b,
         "; a whitespace-only name is not an object")
  jd <- if (is.null(judge)) NULL else .role_text_values(data[[judge]])
  if (!is.null(jd) && any(!is.na(jd) & !nzchar(jd)))
    stop("blank judge identifier(s) in ", judge,
         "; a whitespace-only name is not a judge")
  # count/order must be read as their labelled values: as.numeric() on a
  # factor returns level codes (1..k in label order), silently replacing
  # the real counts/positions -- coerce through the character labels and
  # refuse anything non-numeric
  .num_col <- function(x, what) {
    # as.numeric() on a complex vector keeps the real part and discards the
    # imaginary one without a word, so a complex column is refused outright
    if (is.complex(x))
      stop("`", what, "` is complex; it must hold real numeric values")
    v <- suppressWarnings(as.numeric(if (is.factor(x)) as.character(x) else x))
    if (any(is.na(v) & !is.na(x)))
      stop("`", what, "` has non-numeric value(s); it must be numeric")
    v
  }
  w <- if (is.null(count)) rep(1, nrow(data)) else .num_col(data[[count]], count)
  if (!is.null(count) &&
      (any(!is.finite(w[!is.na(w)])) || any(w[!is.na(w)] < 0) ||
       any(w[!is.na(w)] != floor(w[!is.na(w)]))))
    stop("`", count, "` must hold whole non-negative replication counts; ",
         "zero-count rows are dropped")
  ord <- if (is.null(order)) NULL else .num_col(data[[order]], order)
  if (!is.null(ord) && any(!is.finite(ord[!is.na(ord)])))
    stop("`", order, "` must hold finite sequence values")
  check_order_counts <- function(ord, w) {
    if (!is.null(ord) && any(w > 1, na.rm = TRUE))
      stop("count-compressed rows cannot be combined with `order`: each ",
           "ordered comparison must occupy its own row with a unique sequence ",
           "value so exposure and carry-over histories are defined", call. = FALSE)
  }
  notes <- character(0)
  if (!is.null(anchors)) {
    if (!is.numeric(anchors) || is.complex(anchors) ||
        !is.null(dim(anchors)) || !is.null(oldClass(anchors)) ||
        !length(anchors) || is.null(names(anchors)) ||
        anyNA(names(anchors)) || any(!nzchar(trimws(names(anchors)))))
      stop("`anchors` must be a non-empty named numeric vector ",
           "(names = object names)")
    names(anchors) <- .role_text_values(names(anchors))
    if (anyDuplicated(names(anchors)))
      stop("duplicate anchor(s) after trimming for: ",
           paste(unique(names(anchors)[duplicated(names(anchors))]),
                 collapse = ", "),
           "; each object may be anchored once")
    if (any(!is.finite(anchors)))
      stop("anchor value(s) must be finite: ",
           paste(names(anchors)[!is.finite(anchors)], collapse = ", "))
    # every anchor name must match an object: a misspelled name silently
    # dropped would leave the intended object free (se != 0) while the user
    # believes it is fixed -- refuse and name the offenders
    bad_anch <- setdiff(names(anchors), unique(c(a, b)))
    if (length(bad_anch))
      stop("`anchors` name(s) do not match any object in the data: ",
           paste(bad_anch, collapse = ", "))
  }
  # a constant, object_a-oriented covariate for the first-position advantage;
  # appended to the dependence design (alone, or beside exposure/carry-over)
  add_pos <- function(Z, n)
    if (isTRUE(position)) cbind(Z, position = rep(1, n)) else Z

  if (!is.null(response)) {
    xr <- data[[response]]
    # a row with no weight is not data. An ordered factor STATES its scale
    # in its levels, so those stand whatever the weights are (an empty
    # declared extreme is an identifiability question, handled downstream);
    # a numeric response has no declared range, so its top category must
    # come from the rows the analysis keeps, or one zero-count row would
    # change m in an otherwise identical model.
    usable <- !is.na(a) & !is.na(b) & a != b & !is.na(w) & w > 0
    if (!is.null(jd)) usable <- usable & !is.na(jd)
    if (!is.null(ord)) usable <- usable & !is.na(ord)
    if (is.factor(xr)) {
      # a plain factor's alphabetical level order would silently define the
      # response scale (and can reverse it); the order must be explicit
      if (!is.ordered(xr))
        stop("a polytomous response factor must be ORDERED ",
             "(factor(..., ordered = TRUE) with levels from worst to ",
             "best), or supply integer scores 0..m: an unordered ",
             "factor's alphabetical levels would silently define -- and ",
             "can reverse -- the response scale")
      cats <- levels(xr); x <- as.integer(xr) - 1L
    } else {
      .check_integer_scores(xr[usable], "the polytomous response")
      xn <- suppressWarnings(as.numeric(as.character(xr)))
      x <- as.integer(xn)
      if (any(x[usable] < 0, na.rm = TRUE))
        stop("polytomous responses must be non-negative integers 0..m")
      if (!any(usable & !is.na(x)))
        stop("no usable comparisons: every response is missing")
      cats <- as.character(0:max(x[usable], na.rm = TRUE))
      x[!usable] <- NA_integer_          # a dropped score is not a category
    }
    keep <- usable & !is.na(x)
    if (any(!keep)) {
      notes <- c(notes, sprintf(
        "%d row(s) dropped (missing, zero-count, or self-comparison)",
        sum(!keep)))
      a <- a[keep]; b <- b[keep]; x <- x[keep]; w <- w[keep]
      if (!is.null(jd)) jd <- jd[keep]
      if (!is.null(ord)) ord <- ord[keep]
    }
    if (!length(a)) stop("no usable comparisons")
    check_order_counts(ord, w)
    Z <- if (is.null(ord)) NULL else
      .btl_exposure(a, b, x, length(cats) - 1L, jd, ord, w)
    Z <- add_pos(Z, length(a))
    return(.btl_graded(a, b, x, jd, w, cats, maxit, tol, notes,
                       thr = thresholds, Z = Z, ord = ord, anchors = anchors,
                       object_design = .object_design))
  }

  if (!is.null(margin)) {
    # winner + margin entry: orientation-free by construction. The polytomous
    # response is assembled from "who won" and "by how much"; a winner value
    # matching neither object is a tie and becomes the middle category.
    mg <- data[[margin]]
    # the margin's level order defines the polytomous scale, so it must be
    # explicit: an ordered factor (smallest to largest margin) or a numeric
    # magnitude. A plain factor's -- or a character column's -- alphabetical
    # order can silently reverse which margin counts as the big win.
    if (is.factor(mg) && !is.ordered(mg))
      stop("`margin` must be an ORDERED factor (factor(..., ordered = ",
           "TRUE), levels smallest to largest margin) or a numeric ",
           "magnitude: a plain factor's alphabetical level order would ",
           "silently define -- and can reverse -- the margin scale")
    if (is.character(mg))
      stop("`margin` is a character column; supply an ORDERED factor ",
           "(levels smallest to largest margin) or a numeric magnitude, ",
           "so the margin order is explicit rather than alphabetical")
    # anything else that merely sorts -- a logical, a complex, a Date --
    # would produce plausible-looking categories from a scale that is not a
    # margin
    if (!is.factor(mg)) {
      if (!is.numeric(mg) || is.complex(mg) || !is.null(dim(mg)) ||
          !is.null(oldClass(mg)))
        stop("`margin` must be an ORDERED factor (levels smallest to ",
             "largest margin) or a plain numeric magnitude; a ",
             paste(class(mg), collapse = "/"), " column is not a margin")
    }
    wn <- .role_text_values(data[[winner]])
    is_a <- !is.na(wn) & wn == a
    is_b <- !is.na(wn) & wn == b
    tie <- !is.na(wn) & !is_a & !is_b & tolower(wn) %in% c("tie", "draw")
    miss_wn <- !is.na(wn) & !is_a & !is_b & !tie
    # the scale is read from the rows the analysis keeps: a row with no
    # weight, no winner, or a winner naming neither object is not data, and
    # letting one contribute a margin level would add a category nothing was
    # judged in -- or fail the fit outright
    usable <- !is.na(a) & !is.na(b) & a != b & !is.na(w) & w > 0 &
      !is.na(wn) & !miss_wn
    if (!is.null(jd)) usable <- usable & !is.na(jd)
    if (!is.null(ord)) usable <- usable & !is.na(ord)
    # a tie carries no margin -- it becomes the middle category -- so its
    # margin value must not open win and loss categories nothing was
    # judged in, which would leave both extremes empty and stop the fit
    scale_rows <- usable & !tie
    # Only retained non-ties define the margin scale. Dropped rows and ties
    # carry no analysed margin; for an analysed win, zero is itself a tie and
    # cannot be labelled as a positive margin of victory.
    if (!is.factor(mg)) {
      kept_margin <- mg[scale_rows & !is.na(mg)]
      if (any(!is.finite(kept_margin)) || any(kept_margin <= 0))
        stop("retained non-tie `margin` magnitudes must be finite and ",
             "positive: zero is a tie, not a margin of victory")
    }
    lv <- if (is.factor(mg)) levels(droplevels(mg[scale_rows])) else
      as.character(sort(unique(mg[scale_rows & !is.na(mg)])))
    q <- length(lv)
    if (q < 1L) stop("`margin` has no usable levels")
    mgi <- match(as.character(mg), lv)
    if (any(miss_wn))
      notes <- c(notes, sprintf(
        "%d row(s) with winner matching neither object treated as missing",
        sum(miss_wn)))
    # a tie only in a dropped row is not a tie in the data either
    ties_present <- any(tie & usable)
    keep <- usable & (tie | !is.na(mgi))
    if (!is.null(jd)) keep <- keep & !is.na(jd)
    if (!is.null(ord)) keep <- keep & !is.na(ord)
    if (any(!keep)) {
      notes <- c(notes, sprintf(
        "%d row(s) dropped (missing winner or margin, zero-count, or self-comparison)",
        sum(!keep)))
      a <- a[keep]; b <- b[keep]; w <- w[keep]
      mgi <- mgi[keep]; is_a <- is_a[keep]; is_b <- is_b[keep]
      tie <- tie[keep]
      if (!is.null(jd)) jd <- jd[keep]
      if (!is.null(ord)) ord <- ord[keep]
    }
    if (!length(a)) stop("no usable comparisons")
    check_order_counts(ord, w)
    base <- q - 1L + as.integer(ties_present)
    x <- ifelse(is_a, base + mgi, ifelse(is_b, q - mgi, q))
    cats <- c(paste0("worse by ", rev(lv)), if (ties_present) "tie",
              paste0("better by ", lv))
    if (ties_present)
      notes <- c(notes, sprintf("%d tie(s) placed in the middle category",
                                sum(tie)))
    Z <- if (is.null(ord)) NULL else
      .btl_exposure(a, b, as.integer(x), length(cats) - 1L, jd, ord, w)
    Z <- add_pos(Z, length(a))
    return(.btl_graded(a, b, as.integer(x), jd, w, cats, maxit, tol, notes,
                       thr = thresholds, Z = Z, ord = ord, anchors = anchors,
                       object_design = .object_design))
  }

  wn <- .role_text_values(data[[winner]])
  keep <- !is.na(a) & !is.na(b) & !is.na(wn) & a != b & !is.na(w) & w > 0
  if (!is.null(jd)) keep <- keep & !is.na(jd)
  if (!is.null(ord)) keep <- keep & !is.na(ord)
  if (any(!keep)) {
    notes <- c(notes, sprintf("%d row(s) dropped (missing, zero-count, or self-comparison)",
                              sum(!keep)))
    a <- a[keep]; b <- b[keep]; wn <- wn[keep]; w <- w[keep]
    if (!is.null(jd)) jd <- jd[keep]
    if (!is.null(ord)) ord <- ord[keep]
  }
  if (!length(a)) stop("no usable comparisons")
  if (!is.null(ord) && ties == "half")
    stop("exposure analysis is incompatible with ties = 'half';",
         " drop ties or code them as a polytomous middle category")

  # outcome: 1 = a wins, 0 = b wins; an explicit "tie"/"draw" entry is a
  # tie; anything else matching neither object is missing, not a tie
  y <- ifelse(wn == a, 1, ifelse(wn == b, 0, NA))
  is_tie <- is.na(y) & tolower(wn) %in% c("tie", "draw")
  miss <- is.na(y) & !is_tie
  if (any(miss)) {
    notes <- c(notes, sprintf(
      "%d row(s) with winner matching neither object treated as missing",
      sum(miss)))
    sel <- !miss
    a <- a[sel]; b <- b[sel]; y <- y[sel]; w <- w[sel]
    if (!is.null(jd)) jd <- jd[sel]
    if (!is.null(ord)) ord <- ord[sel]
    if (!length(a)) stop("no usable comparisons")
  }
  row_cluster <- row_replicates <- NULL
  if (anyNA(y)) {
    n_tie <- sum(is.na(y))
    if (ties == "error") stop(n_tie, " tie(s) present; set ties = 'drop' or 'half'")
    if (ties == "drop") {
      notes <- c(notes, sprintf("%d tie(s) dropped", n_tie))
      sel <- !is.na(y)
      a <- a[sel]; b <- b[sel]; y <- y[sel]; w <- w[sel]
      if (!is.null(jd)) jd <- jd[sel]
      if (!is.null(ord)) ord <- ord[sel]
    } else {
      notes <- c(notes, sprintf("%d tie(s) scored half a win each way (halves are not independent trials)",
                                n_tie))
      t_i <- which(is.na(y))
      # Keep the two half rows tied to their original comparison for the
      # unclustered sandwich. Their scores are two parts of one pragmatic
      # tie contribution, not two independent Bernoulli observations.
      row_cluster <- seq_along(y)
      row_replicates <- w
      a <- c(a, a[t_i]); b <- c(b, b[t_i])
      y[t_i] <- 1; y <- c(y, rep(0, length(t_i)))
      w[t_i] <- w[t_i] / 2; w <- c(w, w[t_i])
      row_cluster <- c(row_cluster, row_cluster[t_i])
      row_replicates <- c(row_replicates, row_replicates[t_i])
      if (!is.null(jd)) jd <- c(jd, jd[t_i])
    }
  }

  if (!is.null(ord) || isTRUE(position)) {
    check_order_counts(ord, w)
    # exposure/position covariates route through the polytomous engine, whose two-
    # category case reproduces the dichotomous analysis exactly
    Z <- if (is.null(ord)) NULL else
      .btl_exposure(a, b, as.integer(y), 1L, jd, ord, w)
    Z <- add_pos(Z, length(a))
    return(.btl_graded(a, b, as.integer(y), jd, w, c("0", "1"), maxit, tol,
                       notes, thr = "free", Z = Z, ord = ord,
                       anchors = anchors, object_design = .object_design,
                       row_cluster = row_cluster,
                       row_replicates = row_replicates))
  }

  # the two-category polytomous engine IS the dichotomous conditional model
  # (their equivalence is tested to machine precision), so one estimator
  # serves both routes; m == 1 results are presented as wins / win
  # proportions inside .btl_graded
  .btl_graded(a, b, as.integer(y), jd, w, c("0", "1"), maxit, tol,
              notes, thr = "free", anchors = anchors,
              object_design = .object_design,
              row_cluster = row_cluster %||% NULL,
              row_replicates = row_replicates %||% NULL)
}

#' @export
print.rasch_btl <- function(x, ...) {
  cat(sprintf("Bradley-Terry-Luce analysis: %d objects, %.0f comparisons%s\n",
              nrow(x$objects), x$n_comparisons,
              if (!is.null(x$judges)) sprintf(", %d judges", nrow(x$judges)) else ""))
  cat(sprintf("Conditional ML: %s in %d iterations; sandwich SEs%s\n",
              if (x$converged) "converged" else "NOT converged", x$iterations,
              if (x$clustered) " clustered by judge" else ""))
  fit_probability <- if (is.finite(x$total_p))
    paste0("p = ", .fmt_p(x$total_p)) else "probability unavailable"
  cat(sprintf(paste0("Object separation index %.3f; pairwise chi-square ",
                     "%.2f on %d df, %s\n"),
              x$osi$PSI, x$total_chisq, x$total_df, fit_probability))
  if (!is.null(x$anchors))
    cat(sprintf("Anchored at %d object(s) (se = 0): %s\n",
                length(x$anchors), paste(names(x$anchors), collapse = ", ")))
  if (!is.null(x$dependence)) {
    for (r in seq_len(nrow(x$dependence))) {
      # position is a static first-presented advantage, not a within-judge
      # history effect, so it is not labelled "Within-judge"
      lab <- if (x$dependence$effect[r] == "position")
        "First-position advantage" else
        paste("Within-judge", gsub("_", "-", x$dependence$effect[r]))
      # three saved-fit schemas: `t` (current, t reference); `z` with `df`
      # (1.11.4 transitional: t-based inference under the old name); `z`
      # without `df` (older fits whose p-values used the NORMAL reference:
      # keep their own label rather than misrepresent their inference)
      if (!is.null(x$dependence$t)) {
        st_lab <- "t"; st_r <- x$dependence$t[r]
      } else if (!is.null(x$dependence$df)) {
        st_lab <- "t"; st_r <- x$dependence$z[r]
      } else {
        st_lab <- "z"; st_r <- x$dependence$z[r]
      }
      p_label <- if (!is.null(x$dependence$p_adj) &&
                       length(x$dependence$p_adj) >= r &&
                       is.finite(x$dependence$p_adj[r]))
        paste0("Holm p = ", .fmt_p(x$dependence$p_adj[r])) else
          "Holm p unavailable"
      cat(sprintf("%s: %.3f logits (SE %.3f, %s = %.2f, %s)\n",
                  lab, x$dependence$estimate[r], x$dependence$se[r],
                  st_lab, st_r, p_label))
    }
  }
  if (!is.null(x$thresholds)) {
    cat(sprintf("Polytomous comparisons in %d categories%s; symmetric thresholds: %s\n",
                x$m + 1L,
                if (!is.null(x$categories) &&
                    !all(x$categories == as.character(0:x$m)))
                  paste0(" (", paste(x$categories, collapse = " < "), ")")
                else "",
                paste(sprintf("%.3f", x$thresholds$tau), collapse = ", ")))
  }
  print(.fmt_df(x$objects[, intersect(c("object", "location", "se",
                                        "comparisons", "wins", "score",
                                        "fit_resid", "extreme"),
                                      names(x$objects))]), row.names = FALSE)
  if (!is.null(x$judges)) {
    mis <- x$judges[!is.na(x$judges$fit_resid) & abs(x$judges$fit_resid) > 2.5, ]
    cat(sprintf("Judges beyond |fit residual| 2.5: %d%s\n", nrow(mis),
                if (nrow(mis)) paste0(" (", paste(mis$judge, collapse = ", "), ")")
                else ""))
  }
  if (length(x$notes)) cat(sprintf("Notes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

#' Plot Bradley-Terry-Luce object locations
#'
#' Caterpillar plot of object locations with 95 per cent error bars. The
#' intervals use the reference degrees of freedom carried by the fitted
#' covariance: a t reference for judge-based uncertainty and the limiting
#' normal reference otherwise. An interval is omitted when its covariance
#' does not support inference. Objects beyond the specified fit-residual band
#' are marked.
#'
#' @param fit An object from \code{\link{btl}}.
#' @param band Absolute fit-residual value beyond which an object is
#'   highlighted.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pairs <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pairs[, 1], each = 30),
#'                 b = rep(pairs[, 2], each = 30))
#' p <- plogis(beta[d$a] - beta[d$b])
#' d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' plot_btl(btl(d, "a", "b", "win"))
#' @export
plot_btl <- function(fit, band = 2.5) {
  .check_btl_display_fit(fit)
  .check_band(band)
  d <- fit$objects[order(fit$objects$location), ]
  k <- nrow(d)
  ref_df <- vapply(d$object, function(object)
    .btl_equate_cov_df(fit, object), numeric(1))
  critical <- stats::qt(0.975, df = ref_df)
  has_se <- is.finite(d$se) & is.finite(critical)
  lower <- upper <- rep(NA_real_, k)
  lower[has_se] <- d$location[has_se] - critical[has_se] * d$se[has_se]
  upper[has_se] <- d$location[has_se] + critical[has_se] * d$se[has_se]
  xerr <- c(lower[has_se], upper[has_se], d$location)
  xlim <- if (any(is.finite(xerr))) range(xerr, na.rm = TRUE)
          else range(d$location)
  op <- .rr_canvas(xlim + c(-0.15, 0.15) * diff(xlim), c(0.5, k + 0.5),
                   "Location (logits)", "", grid_y = FALSE, grid_x = TRUE,
                   yaxis = FALSE)
  on.exit(par(op))
  mis <- !is.na(d$fit_resid) & abs(d$fit_resid) > band
  segments(lower[has_se], which(has_se), upper[has_se], which(has_se),
           col = ifelse(mis[has_se], .rr$red, .rr$soft), lwd = 2.2)
  points(d$location, seq_len(k), pch = 21, cex = 1.6, lwd = 1.2,
         bg = ifelse(mis, .rr$red, .rr$blue), col = "white")
  text(d$location, seq_len(k), d$object, pos = 3, offset = 0.55, cex = 0.8,
       col = .rr$ink)
  if (any(mis))
    .rr_legend("bottomright", sprintf("|fit residual| > %.1f", band),
               pch = 21, pt.bg = .rr$red, col = "white", pt.cex = 1.4)
  invisible(NULL)
}

# A failed optimiser can leave a complete-looking final iterate in the fit.
# It is useful for diagnosis but is not a calibrated scale and must not reach
# a public plot as though it were an estimate.
.check_btl_display_fit <- function(fit) {
  if (!inherits(fit, "rasch_btl"))
    stop("not a paired-comparison (btl) fit", call. = FALSE)
  if (!isTRUE(fit$converged))
    stop("the paired-comparison calibration did not converge; fitted displays are unavailable",
         call. = FALSE)
  invisible(fit)
}

# ---------------------------------------------------------------------------
# Polytomous paired comparisons: the adjacent-categories (Rasch-type) ordinal
# extension of BTL (Tutz 1986; Agresti 1992). The response is one of m+1
# ordered categories from object_a's perspective ("much worse" ... "much
# better"); category probabilities follow a partial-credit structure on the
# difference beta_a - beta_b with thresholds constrained symmetric,
# tau_k = -tau_{m+1-k}, so the model is invariant to presentation order and
# judge tendencies cancel. m = 1 is exactly BTL; m = 2 has the Davidson
# (1970) ties structure after the endpoint log-strength is mapped to 2*beta.
# Estimation, identification, sandwich errors, and fit
# follow the package conventions established in btl().
# ---------------------------------------------------------------------------
.btl_graded <- function(a, b, x, jd, w, cats, maxit, tol, notes,
                        thr = "free", Z = NULL, ord = NULL, anchors = NULL,
                        object_design = NULL, row_cluster = NULL,
                        row_replicates = NULL, withheld_notes = TRUE) {
  # withheld_notes: emit the notes that explain a statistic this engine
  # withholds from a table the caller may never report. A refit whose
  # pairwise fit and dependence table are discarded must not carry them,
  # or engine boilerplate is prefixed as a finding about that caller's
  # own analysis.
  m <- length(cats) - 1L
  if (m < 1L) stop("polytomous responses need at least two categories")
  if (is.null(row_cluster)) row_cluster <- seq_along(a)
  if (is.null(row_replicates)) row_replicates <- w
  if (length(row_cluster) != length(a) ||
      length(row_replicates) != length(a) ||
      any(!is.finite(row_replicates)) || any(row_replicates <= 0))
    stop("internal independent-comparison bookkeeping is inconsistent",
         call. = FALSE)
  # identifiability: empty EXTREME categories leave no finite spread (the
  # data are evidence of infinite spread, as a zero raw score is of an
  # infinite person location); empty interior categories are unidentified
  # under free thresholds but pooled over by the principal-component
  # structure
  check_cats <- function(x, note_interior) {
    xs <- c(x, m - x)
    emp <- which(tabulate(xs + 1L, m + 1L) == 0) - 1L
    if (any(emp %in% c(0L, m)))
      stop("extreme category never used (in either orientation): ",
           paste(cats[intersect(emp, c(0L, m)) + 1L], collapse = ", "),
           "; no finite threshold estimate exists - collapse categories")
    if (length(emp) && thr == "free")
      stop("interior category never used (in either orientation): ",
           paste(cats[emp + 1L], collapse = ", "),
           "; use thresholds = 'pc' (pooled principal-component structure)",
           " or collapse categories")
    if (length(emp) && note_interior)
      notes <<- c(notes, sprintf(
        "interior category unused (%s); thresholds pooled by the principal-component structure",
        paste(cats[emp + 1L], collapse = ", ")))
    invisible(NULL)
  }
  check_cats(x, note_interior = TRUE)

  # objects whose every response sits at the boundary have no finite
  # estimate, as extreme persons are set aside in a Rasch calibration.
  # Their comparisons are kept so a reporting location can be
  # extrapolated against the calibrated scale after estimation.
  a0 <- a; b0 <- b; x0 <- x; w0 <- w; Z0 <- Z
  jd0 <- jd; ord0 <- ord
  removed_ext <- character(0)
  removed_any <- FALSE
  repeat {
    objs <- sort(unique(c(a, b)))
    T_of <- setNames(numeric(length(objs)), objs); N_of <- T_of
    for (r in seq_along(a)) {
      T_of[a[r]] <- T_of[a[r]] + w[r] * x[r]
      T_of[b[r]] <- T_of[b[r]] + w[r] * (m - x[r])
      N_of[a[r]] <- N_of[a[r]] + w[r] * m
      N_of[b[r]] <- N_of[b[r]] + w[r] * m
    }
    ext <- names(T_of)[T_of == 0 | T_of == N_of]
    if (!length(ext)) break
    # an anchored object at a boundary is held at a known location, so it is
    # not silently removed the way a free extreme object is; equating cannot
    # proceed on a scale whose anchor has no comparisons to place it against
    if (!is.null(anchors) && any(ext %in% names(anchors)))
      stop("anchored object(s) at a response boundary (undefeated or winless): ",
           paste(intersect(ext, names(anchors)), collapse = ", "),
           "; an anchored object cannot be removed - drop it from `anchors`",
           " or supply comparisons that place it")
    notes <- c(notes, sprintf(
      "object(s) at a response boundary set aside from estimation: %s; reported at an extrapolated location (score half a point inside the boundary), standard error withheld",
      paste(ext, collapse = ", ")))
    removed_ext <- c(removed_ext, ext)
    sel <- !(a %in% ext) & !(b %in% ext)
    a <- a[sel]; b <- b[sel]; x <- x[sel]; w <- w[sel]
    row_cluster <- row_cluster[sel]
    row_replicates <- row_replicates[sel]
    if (!is.null(jd)) jd <- jd[sel]
    if (!is.null(ord)) ord <- ord[sel]
    if (!is.null(Z)) Z <- Z[sel, , drop = FALSE]
    removed_any <- TRUE
    if (!length(a)) stop("no comparisons remain after removing extreme objects")
  }
  # removing boundary objects can itself empty a category; re-check so the
  # user gets the collapse-categories error, not a singular Newton step
  if (removed_any) check_cats(x, note_interior = FALSE)
  objs <- sort(unique(c(a, b)))
  K <- length(objs)
  if (K < 3) stop("need at least three comparable objects")
  ia <- match(a, objs); ib <- match(b, objs)
  comp <- .btl_components(K, ia, ib)
  if (length(unique(comp)) > 1) {
    parts <- split(objs, comp)
    stop("the comparison graph is disconnected; components: ",
         paste(vapply(parts, paste, "", collapse = ","), collapse = " | "))
  }
  # Ford's (1957) existence condition: finite maximum-likelihood locations
  # exist only when the directed win graph is strongly connected. A single
  # undefeated (or winless) object was removed above; a CLUSTER that never
  # concedes a point to the rest is the collective form of the same
  # boundary, and the likelihood pushes the divide to infinity while the
  # optimiser stops at enormous finite values that look converged.
  Wp <- matrix(0, K, K)
  pts_a <- w * x; pts_b <- w * (m - x)
  for (r in seq_along(ia)) {
    Wp[ia[r], ib[r]] <- Wp[ia[r], ib[r]] + pts_a[r]
    Wp[ib[r], ia[r]] <- Wp[ib[r], ia[r]] + pts_b[r]
  }
  adj <- Wp > 0
  reach <- function(Amat, start) {
    seen <- rep(FALSE, K); seen[start] <- TRUE; front <- start
    while (length(front)) {
      nxt <- which(colSums(Amat[front, , drop = FALSE]) > 0 & !seen)
      seen[nxt] <- TRUE; front <- nxt
    }
    seen
  }
  anch_idx <- if (!is.null(anchors)) which(objs %in% names(anchors))
              else integer(0)
  if (!length(anch_idx)) {
    if (!(all(reach(adj, 1L)) && all(reach(t(adj), 1L)))) {
      fwd <- reach(adj, 1L); bwd <- reach(t(adj), 1L)
      sep <- objs[!(fwd & bwd)]
      stop("the win graph is not strongly connected: object(s) ",
           paste(sep, collapse = ", "), " never concede points to (or never ",
           "take points from) the rest, so finite maximum-likelihood ",
           "locations do not exist (Ford 1957) -- remove the separated ",
           "object(s), or collect comparisons (or ties/margins) that cross ",
           "the divide in the missing direction", call. = FALSE)
    }
  } else {
    # anchored fits relax Ford: an anchor pins its object, so a divergent
    # (recession) direction exists only for objects not tied to an anchor
    # in BOTH constraint directions. With edge j -> i whenever i took
    # points off j (the constraint d_i >= d_j on any recession direction
    # d, with d = 0 at anchors), an object is safe iff an anchor is
    # reachable from it (bounded above) and it is reachable from an
    # anchor (bounded below); otherwise its up- or down-set contains no
    # anchor and the likelihood recedes along that set
    up <- reach(adj, anch_idx); down <- reach(t(adj), anch_idx)
    bad <- !(up & down)
    if (any(bad))
      stop("object(s) ", paste(objs[bad], collapse = ", "), " are not ",
           "tied to an anchor in both win directions: even with the ",
           "anchors fixed, their likelihood recedes to a boundary and no ",
           "finite maximum exists -- anchor an object in the separated ",
           "cluster or collect comparisons that cross the divide",
           call. = FALSE)
  }

  # symmetric-threshold map: tau = Cmat %*% tfree, tau_k = -tau_{m+1-k}.
  # Under thr = "pc" the free symmetric parameters are further pooled to
  # the spread (linear) component alone - Andrich's principal-component
  # structure with the even (skewness) components structurally zero under
  # symmetry - so sparse categories borrow strength from every response.
  if (thr == "pc" && m >= 2L) {
    v1 <- seq_len(m) - (m + 1) / 2
    Cmat <- cbind(v1 / sqrt(sum(v1^2)))
    q <- 1L
  } else {
    q <- m %/% 2L
    Cmat <- matrix(0, m, q)
    for (k in seq_len(q)) { Cmat[k, k] <- 1; Cmat[m + 1L - k, k] <- -1 }
  }
  # location design B (K x nb) and fixed offset beta0 (length K), so that
  # beta = Bmat %*% bfree + beta0. Without anchors this is the sum-zero map
  # (bit-identical to the pre-anchor path); with anchors it is a selection
  # matrix -- identity rows for the free objects, zero rows for the anchored
  # ones -- plus the anchored values as the offset, so the anchored locations
  # are held fixed and the rest float with no sum-zero constraint, the origin
  # and scale coming from the anchors (as in an anchored rasch() calibration).
  anch <- NULL
  object_parameter_names <- NULL
  if (!is.null(object_design)) {
    if (!is.null(anchors))
      stop("explanatory object restrictions cannot be combined with anchors")
    if (!is.list(object_design) || !is.matrix(object_design$B) ||
        is.null(rownames(object_design$B)))
      stop("the internal object design must contain a named matrix `B`")
    miss <- setdiff(objs, rownames(object_design$B))
    if (length(miss))
      stop("object predictor metadata are missing after data preparation: ",
           paste(miss, collapse = ", "))
    Bmat <- object_design$B[objs, , drop = FALSE]
    beta0 <- object_design$offset %||% setNames(numeric(nrow(object_design$B)),
                                                rownames(object_design$B))
    if (is.null(names(beta0)))
      stop("the internal explanatory object offset must be named")
    beta0 <- as.numeric(beta0[objs])
    if (anyNA(beta0)) stop("the explanatory object offset is incomplete")
    Bcheck <- if (ncol(Bmat))
      sweep(Bmat, 2L, .design_column_scale(Bmat), `/`) else Bmat
    if (!ncol(Bmat) || qr(Bcheck, tol = 1e-10)$rank < ncol(Bmat) ||
        qr(cbind(1, Bcheck), tol = 1e-10)$rank < ncol(Bmat) + 1L)
      stop("the explanatory object design is not identified after data preparation")
    if (ncol(Bmat) > K - 1L)
      stop("the explanatory object design has more parameters than the free calibration")
    object_parameter_names <- colnames(Bmat) %||%
      paste0("object_effect", seq_len(ncol(Bmat)))
    notes <- c(notes, sprintf("object locations constrained by %d explanatory effect(s)",
                              ncol(Bmat)))
  } else if (is.null(anchors)) {
    Bmat <- rbind(diag(K - 1L), rep(-1, K - 1L))
    beta0 <- numeric(K)
  } else {
    unavailable_anchors <- setdiff(names(anchors), objs)
    if (length(unavailable_anchors))
      stop("anchored object(s) have no usable comparisons after data ",
           "preparation: ", paste(unavailable_anchors, collapse = ", "),
           "; remove them from `anchors` or supply informative comparisons",
           call. = FALSE)
    anch <- anchors
    apos <- match(names(anch), objs)
    free <- setdiff(seq_len(K), apos)
    if (!length(free)) stop("every object is anchored; nothing to estimate")
    Bmat <- matrix(0, K, length(free))
    Bmat[cbind(free, seq_along(free))] <- 1
    beta0 <- numeric(K); beta0[apos] <- as.numeric(anch)
    # start the free objects AT the anchored scale's centre, as pcml's
    # anchored path shifts its start values: beginning at 0 when the anchors
    # sit several logits away sends the undamped Newton step into the
    # logistic tails and diverges, though the maximum is finite
    beta0[free] <- mean(as.numeric(anch))
    notes <- c(notes, sprintf(
      "%d object(s) anchored; scale origin from anchors", length(anch)))
  }
  # Keep the public design in its supplied units, but fit with bounded
  # columns so damping, convergence and rank checks do not depend on them.
  Bmat_original <- Bmat
  object_scale <- .design_column_scale(Bmat)
  Bmat <- sweep(Bmat, 2L, object_scale, `/`)
  nb <- ncol(Bmat)
  pz <- if (is.null(Z)) 0L else ncol(Z)
  Zfull <- Z                         # all effect columns, for the diagnostic table
  if (pz) {
    keepz <- colSums(abs(Z)) > 0
    if (!all(keepz)) {
      notes <- c(notes, sprintf(
        "dependence effect(s) with no informative comparisons dropped: %s",
        paste(colnames(Z)[!keepz], collapse = ", ")))
      Z <- Z[, keepz, drop = FALSE]; pz <- ncol(Z)
    }
    if (pz) {
      # An order or dependence covariate can be exactly reproduced by the
      # object contrasts in a thin comparison graph. For example, a constant
      # first-position column is an object-potential difference on a tree.
      # The requested effect is then not separately identified; a numerical
      # inverse or ridge would only choose an arbitrary decomposition.
      Dloc <- Bmat[ia, , drop = FALSE] - Bmat[ib, , drop = FALSE]
      augmented <- cbind(Dloc, Z)
      if (qr(augmented, tol = 1e-10)$rank < ncol(augmented))
        stop("the position/dependence effects are not separately identified ",
             "from the object locations in this comparison design; randomise ",
             "presentation order or add comparison cycles before estimating ",
             "these effects", call. = FALSE)
    }
  }
  np <- nb + q + pz
  dep <- numeric(pz)
  sc <- 0:m

  # per-row moments for current parameters: probabilities, E, V, mu4,
  # survivor S_k = P(X >= k) and EXc_k = E[X 1(X >= k)], k = 1..m
  moments <- function(beta, tfree, dep) {
    tau <- if (q) drop(Cmat %*% tfree) else numeric(m)
    d <- beta[ia] - beta[ib]
    if (pz) d <- d + drop(Z %*% dep)
    eta <- outer(d, sc) - matrix(rep(c(0, cumsum(tau)), each = length(d)),
                                 length(d), m + 1L)
    eta <- eta - apply(eta, 1, max)
    P <- exp(eta); P <- P / rowSums(P)
    E <- drop(P %*% sc)
    V <- drop(P %*% sc^2) - E^2
    mu4 <- drop(P %*% sc^4) - 4 * E * drop(P %*% sc^3) +
      6 * E^2 * drop(P %*% sc^2) - 3 * E^4
    S <- t(apply(P, 1, function(p) rev(cumsum(rev(p)))))[, -1L, drop = FALSE]
    EXc <- t(apply(P * rep(sc, each = length(d)), 1,
                   function(p) rev(cumsum(rev(p)))))[, -1L, drop = FALSE]
    list(P = P, E = E, V = pmax(V, 1e-12), mu4 = mu4, S = S, EXc = EXc,
         tau = tau)
  }
  cumInd <- outer(x, seq_len(m), ">=") * 1

  # accumulate per-comparison rows (vector or matrix) into indexed slots:
  # rowsum() replaces the interpreted per-row loops that dominated large fits
  acc <- function(v, idx, nrows) {
    out <- matrix(0, nrows, if (is.matrix(v)) ncol(v) else 1L)
    rs <- rowsum(v, idx)
    out[as.integer(rownames(rs)), ] <- rs
    out
  }

  # generic gradient/Hessian over theta = (beta_red, tfree, dep): the
  # covariates enter the exponent multiplied by the score, exactly as the
  # location difference does, so every block follows the same moments
  gH <- function(mo) {
    resE <- w * (x - mo$E)
    g_beta_full <- drop(acc(resE, ia, K) - acc(resE, ib, K))
    g <- drop(crossprod(Bmat, g_beta_full))
    if (q) g <- c(g, drop(crossprod(Cmat, colSums(w * (mo$S - cumInd)))))
    if (pz) g <- c(g, drop(crossprod(Z, resE)))
    H <- matrix(0, np, np)
    hv <- w * mo$V
    Hbb <- matrix(0, K, K)
    diag(Hbb) <- drop(acc(hv, ia, K) + acc(hv, ib, K))
    # off-diagonal cells accumulated over the unordered pair, so both
    # presentation orders of the same pair land in the same cell
    lo <- pmin(ia, ib); hi <- pmax(ia, ib)
    hp <- rowsum(hv, (lo - 1L) * K + hi)
    kk <- as.integer(rownames(hp))
    i0 <- (kk - 1L) %/% K + 1L; j0 <- (kk - 1L) %% K + 1L
    Hbb[cbind(i0, j0)] <- Hbb[cbind(i0, j0)] - hp
    Hbb[cbind(j0, i0)] <- Hbb[cbind(j0, i0)] - hp
    H[1:nb, 1:nb] <- crossprod(Bmat, Hbb %*% Bmat)
    if (q) {
      CovXc <- mo$EXc - mo$E * mo$S              # rows: Cov(X, 1(X>=k))
      wc <- (w * CovXc) %*% Cmat
      Hbt_full <- acc(wc, ib, K) - acc(wc, ia, K)
      ti <- (nb + 1L):(nb + q)
      H[1:nb, ti] <- crossprod(Bmat, Hbt_full)
      H[ti, 1:nb] <- t(H[1:nb, ti])
      # sum_r w_r Cov(1(X>=i), 1(X>=j)) = ws[max(i,j)] - crossprod term,
      # since S_r[max(i, j)] depends only on the larger index
      ws <- colSums(w * mo$S)
      Mcc <- outer(seq_len(m), seq_len(m), function(i, j) ws[pmax(i, j)]) -
        crossprod(mo$S, w * mo$S)
      H[ti, ti] <- crossprod(Cmat, Mcc %*% Cmat)
    }
    if (pz) {
      zi <- (nb + q + 1L):np
      wv <- w * mo$V
      H[zi, zi] <- crossprod(Z, Z * wv)
      Hbz_full <- acc(Z * wv, ia, K) - acc(Z * wv, ib, K)
      H[1:nb, zi] <- crossprod(Bmat, Hbz_full)
      H[zi, 1:nb] <- t(H[1:nb, zi])
      if (q) {
        CovXc <- mo$EXc - mo$E * mo$S
        wc <- (w * CovXc) %*% Cmat
        Htz <- -crossprod(wc, Z)                     # q x pz
        H[(nb + 1L):(nb + q), zi] <- Htz
        H[zi, (nb + 1L):(nb + q)] <- t(Htz)
      }
    }
    list(g = g, H = H, resE = resE)
  }

  beta <- beta0; tfree <- numeric(q)
  repeat {
    for (it in seq_len(maxit)) {
      mo <- moments(beta, tfree, dep)
      gh <- gH(mo)
      step <- solve(gh$H, gh$g)
      # trust-region damp: an undamped Newton step of many logits overshoots
      # into the flat logistic tails (seen with distant anchors); inert for
      # ordinary fits, whose steps are far smaller
      location_step <- drop(Bmat %*% step[seq_len(nb)])
      ms <- max(abs(c(location_step, step[-seq_len(nb)])))
      if (is.finite(ms) && ms > 5) step <- step * (5 / ms)
      beta <- beta + drop(Bmat %*% step[1:nb])
      if (q) tfree <- tfree + step[(nb + 1L):(nb + q)]
      if (pz) dep <- dep + step[(nb + q + 1L):np]
      if (max(abs(c(drop(Bmat %*% step[seq_len(nb)]),
                    step[-seq_len(nb)]))) < tol) break
    }
    # a dependence effect running to a boundary is separation: its informative
    # comparisons all point one way, so the data are evidence of an infinite
    # effect (as an all-wins object is of an infinite location). The gradient
    # plateaus there, so the usual test would report convergence with an
    # arbitrary, maxit-dependent estimate; instead the column is set aside
    # with a note and the model refitted without it.
    runaway <- if (pz) which(abs(dep) > 10) else integer(0)
    if (!length(runaway)) break
    notes <- c(notes, sprintf(
      "dependence effect(s) separated (one-sided informative comparisons) and dropped: %s",
      paste(colnames(Z)[runaway], collapse = ", ")))
    Z <- Z[, -runaway, drop = FALSE]; pz <- ncol(Z)
    np <- nb + q + pz
    dep <- numeric(pz)
    beta <- beta0; tfree <- numeric(q)
  }
  mo <- moments(beta, tfree, dep)
  tau <- mo$tau
  loglik <- sum(w * log(pmax(mo$P[cbind(seq_along(x), x + 1L)], 1e-300)))
  gh <- gH(mo)
  resE <- gh$resE
  # The gradient uses bounded design columns and is per count-weighted
  # comparison: neither predictor units nor replication sets this tolerance.
  converged <- isTRUE(max(abs(gh$g)) < 1e-6 * sum(w))
  if (!converged)
    warning("btl estimation did NOT converge in ", it, " iterations: ",
            "estimates, standard errors, and fit statistics are ",
            "unreliable -- increase maxit or check the comparison design",
            call. = FALSE)

  # Godambe sandwich over the full parameter, clustered by judge (each
  # cluster's score contributions accumulated by rowsum, not per-row loops)
  cl <- if (is.null(jd)) as.character(row_cluster) else jd
  ucl <- unique(cl)
  nc <- length(ucl); cidx <- match(cl, ucl)
  # the clustered sandwich estimates the meat from between-judge variation:
  # with one judge it is identically ~zero (SEs collapse to ~1e-16), and
  # with very few judges it understates -- refuse the former, note the latter
  if (!is.null(jd)) {
    if (nc < 2L)
      stop("judge-clustered standard errors need at least 2 judges (got ",
           nc, "); with a single judge drop judge= so comparisons are ",
           "treated as independent, or supply more judges")
  }
  Gm <- matrix(0, nc, np, dimnames = list(ucl, NULL))
  # beta block: per-cluster sums of resE into the winner / loser slots,
  # laid out cluster-major so one rowsum fills the (cluster x object) grid
  # Without judge clusters, a count-weighted row represents that many
  # independent repeats. Divide its aggregate score by sqrt(count) before
  # forming the meat. For a half-scored tie, both expanded rows share the
  # original comparison cluster and original count, so their scaled scores
  # are summed before squaring.
  meat_scale <- if (is.null(jd)) 1 / sqrt(row_replicates) else
    rep(1, length(resE))
  resE_meat <- resE * meat_scale
  gA <- drop(acc(resE_meat, (cidx - 1L) * K + ia, nc * K))
  gB <- drop(acc(resE_meat, (cidx - 1L) * K + ib, nc * K))
  Gm[, 1:nb] <- t(matrix(gA - gB, nrow = K)) %*% Bmat
  if (q) {
    st_tau <- (w * meat_scale * (mo$S - cumInd)) %*% Cmat
    Gm[, (nb + 1L):(nb + q)] <- acc(st_tau, cidx, nc)
  }
  if (pz) Gm[, (nb + q + 1L):np] <-
    acc(Z * resE_meat, cidx, nc)
  H <- gh$H
  # Identification of the free-parameter information. Two failure modes:
  #  (i) genuine singularity (duplicate/degenerate objects, an exactly
  #      zero-information direction): rcond collapses to ~0 -- refuse.
  #  (ii) (quasi-)separation of a SUBSET of objects: a cluster linked to
  #      the rest by too few informative comparisons (e.g. cross-divide
  #      comparisons all at the ceiling category, with one near-ceiling
  #      concession that makes the win graph strongly connected so Ford
  #      passes). The between-cluster contrast is then near-flat and its
  #      location runs to the trust-region boundary while the ridged
  #      inverse reports plausible-looking SEs. Unlike (i) the information
  #      is only mildly ill-conditioned (the ceiling observations inject a
  #      little curvature), so no single rcond threshold separates it from
  #      a legitimately sparse design without risking false refusals. The
  #      robust signature is the CONJUNCTION of an ill-conditioned
  #      direction and a location driven to the boundary along it; report
  #      it as non-convergence with the affected SEs withheld, per the
  #      standard remedy for a boundary estimate.
  rc <- tryCatch(rcond(H), error = function(e) 0)
  h_values <- tryCatch(eigen((H + t(H)) / 2, symmetric = TRUE,
                             only.values = TRUE)$values,
                       error = function(e) NA_real_)
  h_scale <- if (all(is.finite(h_values)) && length(h_values))
    max(abs(h_values)) else NA_real_
  h_positive <- is.finite(h_scale) && h_scale > 0 &&
    min(h_values) > sqrt(.Machine$double.eps) * h_scale
  if (!(is.finite(rc) && rc > 1e-8) || !h_positive)
    stop("the information matrix is singular or not positive definite (reciprocal condition number ",
         format(rc, digits = 3), "): an object location is not identified ",
         "-- typically duplicate objects or a comparison design with a ",
         "zero-information direction. Add comparisons that place the ",
         "affected object(s), or anchor them", call. = FALSE)
  Hi <- solve(H)
  # CR1 removes the leading finite-cluster scale bias, but it cannot make a
  # rank-deficient meat estimable and is not a substitute for a small-cluster
  # variance correction. With fewer than ten judges, or no more judges than
  # fitted parameters, retain the point estimates and descriptive fit but
  # withhold covariance-based inference.
  cr1 <- if (!is.null(jd) && nc > 1L) nc / (nc - 1) else 1
  # The number of sampling units is only a necessary rank condition. Judges,
  # or independent rows, can supply repeated score vectors, in which case
  # the empirical meat still has no variance in some fitted directions.
  # Assess the actual score matrix in the free-parameter space. The relative
  # sqrt(eps) tolerance matches the positive-definiteness guard above and
  # treats numerically vanishing score directions as unavailable inference.
  sval <- tryCatch(svd(Gm, nu = 0L, nv = 0L)$d,
                    error = function(e) NA_real_)
  score_rank <- if (!length(sval) || any(!is.finite(sval)) ||
                    max(sval) <= 0) 0L
                else sum(sval > max(sval) * sqrt(.Machine$double.eps))
  rank_deficient <- score_rank < np
  # effective number of clusters under unequal allocation: the inverse
  # Simpson index of each judge's share of the comparisons. The clustered
  # sandwich is calibrated when the work spreads across judges but not when
  # it concentrates: in simulation (t reference, df = clusters - 1;
  # tools/simval/studies/followups/btl_share_sweep.R) the null rejection
  # is nominal for balanced designs at 10-50 judges (~4.5-5.6% at
  # effective counts 9.7+), ~6% at 8 effective, 7.2-7.5% at 6-7, and ~9%
  # at ~4 -- regardless of the raw judge count. Withhold below 8 effective
  # clusters (where the documented inflation reaches ~7%), annotate in
  # [8, 9.5) (mild, ~6%); balanced designs sit within multinomial noise of
  # their judge count (~9.7-9.9 at J=10, calibrated) and stay silent. The
  # parameter-count condition also uses the effective count: a covariance
  # concentrated on few effective clusters cannot support more parameters
  # than a balanced one could.
  nc_eff <- if (!is.null(jd)) {
    shr <- tapply(w, jd, sum); shr <- shr / sum(shr)
    1 / sum(shr^2)
  } else Inf
  # simulation-only escape hatch: the calibration sweep that SET the
  # concentration thresholds (tools/simval/studies/followups/
  # btl_share_sweep.R) must be able to measure the withheld
  # 4-7-effective-judge region, or its evidence could never be reproduced
  # against the guarded package. The override lifts ONLY the
  # concentration conditions (effective count and its parameter
  # comparison); the nominal cluster-count and rank conditions can never
  # be bypassed -- with too few clusters or a rank-deficient meat there
  # is no covariance worth measuring, only fabrication.
  guard_off <- isTRUE(getOption("rasch.btl_guard_override", FALSE))
  eff_underparam <- !is.null(jd) && nc_eff <= np
  cluster_inference <- !rank_deficient && (is.null(jd) ||
    (nc >= 10L && nc > np &&
     (guard_off || (nc_eff >= 8 && !eff_underparam))))
  if (!cluster_inference) {
    if (is.null(jd)) {
      notes <- c(notes, sprintf(
        paste0("the empirical independent-comparison score covariance has ",
               "rank %d for %d fitted parameters; sandwich uncertainty and ",
               "covariance-based inference are withheld, while point ",
               "estimates and fit summaries remain descriptive"),
        score_rank, np))
    } else {
      notes <- c(notes, sprintf(
      paste0("%d judge clusters (%.1f effective) for %d parameters: ",
             "cluster-robust inference is withheld%s; point estimates and ",
             "fit summaries remain descriptive -- use at least 10 judges ",
             "(and at least 8 effective: spread the comparisons rather ",
             "than concentrating them on few judges) and more effective ",
             "judges than fitted parameters, or a design-level bootstrap"),
      nc, nc_eff, np,
      if (rank_deficient) " because the empirical covariance is rank-deficient"
      else if (nc < 10L || nc <= np) " because the cluster count is too small"
      else if (eff_underparam) paste0(" because the effective cluster count ",
        "does not exceed the parameter count (a concentration heuristic, ",
        "not a statement of mathematical rank deficiency)")
      else " because the comparison allocation concentrates on too few judges"))
    }
  }
  max_share <- if (!is.null(jd)) {
    shr <- tapply(w, jd, sum); max(shr) / sum(shr)
  } else 0
  if (cluster_inference && !is.null(jd) &&
      (nc_eff < 9.5 || max_share > 0.2))
    notes <- c(notes, sprintf(
      paste0("comparison allocation is uneven across judges (%.1f effective ",
             "clusters from %d; largest single-judge share %.0f%%): ",
             "clustered standard errors may be mildly anti-conservative"),
      nc_eff, nc, 100 * max_share))
  covth <- Hi %*% (crossprod(Gm) * cr1) %*% Hi
  # composite-likelihood information ingredients: tr(H^-1 J) = tr(covth H)
  # is the effective parameter count of the Godambe penalty (Varin & Vidoni
  # 2005). H is the positive sensitivity checked above, so a negative trace
  # is invalid rather than a sign convention to hide with abs(). Independent
  # units are judges when clustered, else the count-weighted comparisons.
  # when clustered inference is withheld (few clusters), the Godambe meat
  # is rank-deficient by construction, so the effective parameter count --
  # and any information criterion built on it -- is withheld with it,
  # exactly as the documentation states
  eff_params <- if (cluster_inference) sum(diag(covth %*% H)) else NA_real_
  if (!is.finite(eff_params) || eff_params < 0) eff_params <- NA_real_
  cl_info <- list(eff_params = eff_params,
                  n_units = if (is.null(jd)) sum(w) else length(ucl),
                  n_units_effective = if (is.null(jd)) NA_real_ else nc_eff,
                  n_parameters = np,
                  score_rank = score_rank,
                  inference_available = cluster_inference)
  # anchored objects have a zero row in Bmat, so their location variance is
  # structurally zero (se == 0): the location is a fixed constant, not an estimate
  cov_beta <- Bmat %*% covth[1:nb, 1:nb, drop = FALSE] %*% t(Bmat)
  # rows carry the calibrated objects' names so downstream consumers align
  # by name; the reporting table may hold extra extrapolated boundary rows
  dimnames(cov_beta) <- list(objs, objs)
  se <- sqrt(pmax(diag(cov_beta), 0))
  if (!cluster_inference) {
    se[] <- NA_real_
    if (length(anch_idx)) se[anch_idx] <- 0
  }
  # (quasi-)separation of an object subset: inspect the weak eigendirections
  # of the information and ask whether the fitted locations have run at
  # least 3 logits ALONG one of those directions. Absolute locations cannot
  # be used here: on an anchored scale a pure translation (anchors 0 -> 100)
  # changes every beta but changes neither the likelihood nor identification.
  # Mapping the weak parameter directions through Bmat also prevents a weak
  # threshold/dependence direction from spuriously condemning a large but
  # well-placed object. Anchored rows have zero loading and can never be
  # labelled as estimates that ran away.
  sep_run <- rep(FALSE, K)
  if (rc < 1e-2) {
    eh <- eigen((H + t(H)) / 2, symmetric = TRUE)
    ev_max <- max(abs(eh$values))
    weak_dir <- if (is.finite(ev_max) && ev_max > 0)
      which(abs(eh$values) / ev_max < 1e-2) else integer(0)
    beta_rel <- beta - if (length(anch_idx)) mean(beta[anch_idx]) else mean(beta)
    for (jj in weak_dir) {
      loc_dir <- drop(Bmat %*% eh$vectors[seq_len(nb), jj])
      ld2 <- sum(loc_dir^2)
      if (!is.finite(ld2) || ld2 < 1e-10) next
      run_coef <- sum(loc_dir * beta_rel) / ld2
      run_loc <- run_coef * loc_dir
      # Require an actual three-logit displacement in object-location
      # space, not merely a large coefficient caused by a direction whose
      # location loading is numerically tiny.
      if (all(is.finite(run_loc)) && max(abs(run_loc)) >= 3) {
        affected <- abs(run_loc) >= 0.25 * max(abs(run_loc))
        if (length(anch_idx)) affected[anch_idx] <- FALSE
        sep_run <- sep_run | affected
      }
    }
  }
  if (any(sep_run)) {
    converged <- FALSE
    # The weak direction belongs to the joint location/threshold/dependence
    # solution. Once that solution is rejected, none of its sandwich
    # uncertainty is inferentially usable: retaining the unaffected-looking
    # marginal SEs, or the composite effective parameter count, presents the
    # covariance of a fit we have just declared invalid.
    cluster_inference <- FALSE
    cl_info$eff_params <- NA_real_
    cl_info$inference_available <- FALSE
    se[] <- NA_real_
    if (length(anch_idx)) se[anch_idx] <- 0
    notes <- c(notes, sprintf(
      "object(s) %s have run to the location boundary and the design does not identify them (a cluster linked to the rest by too few informative comparisons -- e.g. cross-divide comparisons all at an extreme category); their standard errors are withheld and the fit is marked not converged. Add comparisons that place them, or anchor them",
      paste(objs[sep_run], collapse = ", ")))
  }
  if (!isTRUE(converged) && isTRUE(cluster_inference)) {
    # A solver that stops short of convergence has not supplied an estimated
    # covariance, even when its last Hessian happens to be invertible.
    cluster_inference <- FALSE
    cl_info$eff_params <- NA_real_
    cl_info$inference_available <- FALSE
    se[] <- NA_real_
    if (length(anch_idx)) se[anch_idx] <- 0
  }
  dependence <- NULL
  if (pz) {
    zi <- (nb + q + 1L):np
    dse <- sqrt(pmax(diag(covth)[zi], 0))
    if (!cluster_inference) dse[] <- NA_real_
    # Clustered Wald statistics get a t reference with G - 1 degrees of
    # freedom (the standard few-cluster correction) rather than normal
    # theory, so five judges give honestly wide p-values
    t_df <- if (!cluster_inference) NA_real_ else if (!is.null(jd))
      max(nc - 1L, 1L) else Inf
    dep_stat <- .wald_ratio(dep, dse)
    dependence <- data.frame(
      effect = colnames(Z), estimate = dep, se = dse,
      t = dep_stat, df = t_df,
      p = 2 * pt(-abs(dep_stat), df = t_df),
      # count-weighted: the number of comparisons (not rows) that carry
      # information about each effect
      n_informative = vapply(seq_len(ncol(Z)), function(j)
        sum(w[Z[, j] != 0]), 0))
    # The sequential carry-over statistic is more persistent within judge than
    # the position and exposure covariates. Two independent null batteries at
    # 14 judges rejected 7.5% and 8.25% of the time; the rate returned to 5.25%
    # at 30 judges. Retain its estimate and clustered SE for description, but
    # do not attach an uncalibrated probability below the validated boundary.
    carry_small <- !is.null(jd) && nc < 30L &
      dependence$effect == "carry_over"
    if (any(carry_small)) {
      dependence$df[carry_small] <- NA_real_
      dependence$p[carry_small] <- NA_real_
      # Same rule as the pairwise note: it explains a withheld entry of this
      # dependence table, so only a caller that reports the table carries it.
      if (withheld_notes) notes <- c(notes, paste0(
        "carry-over probability withheld with fewer than 30 judges: null ",
        "simulation found mild anti-conservatism at 14 judges; the estimate ",
        "and clustered standard error remain descriptive"))
    }
    dependence$p_adj <- NA_real_
    usable_p <- is.finite(dependence$p)
    dependence$p_adj[usable_p] <- stats::p.adjust(
      dependence$p[usable_p], method = "holm", n = nrow(dependence))
    dependence$significant <- ifelse(
      is.finite(dependence$p_adj), dependence$p_adj < 0.05, NA)
    rownames(dependence) <- NULL
  }
  thresholds <- NULL; components <- NULL
  if (q) {
    ti <- (nb + 1L):(nb + q)
    cov_tau <- Cmat %*% covth[ti, ti, drop = FALSE] %*% t(Cmat)
    thresholds <- data.frame(threshold = seq_len(m), tau = tau,
                             se = sqrt(pmax(diag(cov_tau), 0)))
    if (!cluster_inference) thresholds$se[] <- NA_real_
    # principal-component decomposition of the threshold structure: the
    # odd components (spread; kurtosis from five thresholds up) carry the
    # symmetric structure, the even skewness component is structurally
    # zero under presentation-order symmetry
    v1 <- seq_len(m) - (m + 1) / 2
    v1 <- v1 / sqrt(sum(v1^2))
    comp_rows <- list(data.frame(
      component = "spread", estimate = sum(v1 * tau),
      se = sqrt(pmax(drop(t(v1) %*% cov_tau %*% v1), 0))))
    if (m >= 4L) {
      v3 <- (seq_len(m) - (m + 1) / 2)^3
      v3 <- v3 - sum(v3 * v1) * v1
      v3 <- v3 / sqrt(sum(v3^2))
      comp_rows[[2]] <- data.frame(
        component = "kurtosis", estimate = sum(v3 * tau),
        se = if (thr == "pc") 0 else
          sqrt(pmax(drop(t(v3) %*% cov_tau %*% v3), 0)))
    }
    components <- do.call(rbind, comp_rows)
    rownames(components) <- NULL
  } else if (m > 1L) {
    thresholds <- data.frame(threshold = seq_len(m), tau = tau, se = 0)
  }

  # fit: per-comparison z; objects and judges pool their cells
  z <- (x - mo$E) / sqrt(mo$V)
  c4v <- mo$mu4 / mo$V^2 - 1
  n_rows <- sum(w)
  f_cell <- (n_rows - np) / n_rows
  pool <- function(sel) {
    if (sum(w[sel]) < 3)
      return(list(infit_ms = NA_real_, outfit_ms = NA_real_,
                  fit_resid = NA_real_, df = NA_real_, n = sum(w[sel])))
    y2 <- sum(w[sel] * z[sel]^2); f <- f_cell * sum(w[sel])
    # information-weighted infit, as in the dichotomous path but over the
    # polytomous response variance
    wv <- sum(w[sel] * mo$V[sel])
    infit <- if (wv > 1e-12)
      sum(w[sel] * z[sel]^2 * mo$V[sel]) / (f_cell * wv) else NA_real_
    v <- sum(w[sel] * c4v[sel])
    fr <- if (v > 1e-8 && y2 > 0) f * (log(y2) - log(f)) / sqrt(v) else NA_real_
    list(infit_ms = infit, outfit_ms = y2 / f, fit_resid = fr, df = f,
         n = sum(w[sel]))
  }
  ofit <- lapply(seq_len(K), function(k) pool(ia == k | ib == k))
  score_of <- vapply(seq_len(K), function(k)
    sum(w[ia == k] * x[ia == k]) + sum(w[ib == k] * (m - x[ib == k])), 0)
  objects <- data.frame(object = objs, location = beta, se = se,
                        comparisons = vapply(ofit, `[[`, 0, "n"),
                        score = score_of,
                        infit_ms = vapply(ofit, `[[`, 0, "infit_ms"),
                        outfit_ms = vapply(ofit, `[[`, 0, "outfit_ms"),
                        fit_resid = vapply(ofit, `[[`, 0, "fit_resid"),
                        df_fit = vapply(ofit, `[[`, 0, "df"))
  rownames(objects) <- NULL

  # a set-aside boundary object is reported at an extrapolated location:
  # the profile solution with its score moved half a point inside the
  # boundary, against the calibrated objects, thresholds, and dependence
  # effects held fixed. As with an extreme person measure, the location is
  # reported for completeness; the standard error and fit are withheld and
  # the row takes no part in estimation or inference.
  objects$extreme <- FALSE
  if (length(removed_ext)) {
    bl <- setNames(objects$location, objects$object)
    sc0 <- 0:m
    ctau <- c(0, cumsum(tau))
    exp_score <- function(theta, r, as_a) {
      d <- if (as_a) theta - bl[[b0[r]]] else bl[[a0[r]]] - theta
      # dependence effects can be dropped after the boundary rows were set
      # aside (no informative comparisons, or separation); the extrapolation
      # uses the retained effects only, aligned by name
      if (!is.null(Z0) && !is.null(Z) && length(dep) && ncol(Z))
        d <- d + drop(Z0[r, colnames(Z), drop = FALSE] %*% dep)
      eta <- d * sc0 - ctau; eta <- eta - max(eta)
      p <- exp(eta) / sum(exp(eta))
      if (as_a) sum(p * sc0) else m - sum(p * sc0)
    }
    for (e in removed_ext) {
      ra <- which(a0 == e & b0 %in% names(bl))
      rb <- which(b0 == e & a0 %in% names(bl))
      if (!length(ra) && !length(rb)) next  # linked only to other extremes
      Tn <- sum(w0[ra] * x0[ra]) + sum(w0[rb] * (m - x0[rb]))
      Nn <- (sum(w0[ra]) + sum(w0[rb])) * m
      Tstar <- if (Tn <= 0) 0.5 else Nn - 0.5
      g <- function(theta)
        sum(vapply(ra, function(r) w0[r] * exp_score(theta, r, TRUE), 0)) +
        sum(vapply(rb, function(r) w0[r] * exp_score(theta, r, FALSE), 0)) -
        Tstar
      lim <- range(bl) + c(-12, 12)
      loc <- tryCatch(stats::uniroot(g, interval = lim, tol = 1e-8)$root,
                      error = function(err) NA_real_)
      if (!is.finite(loc)) {
        notes <- c(notes, sprintf(
          "no stable extrapolated location for %s; the object is omitted from the object table",
          e))
        next
      }
      row <- objects[1, ]
      row[1, seq_along(row)] <- NA
      row$object <- e
      row$location <- loc
      row$comparisons <- sum(w0[ra]) + sum(w0[rb])
      row$score <- Tn
      row$extreme <- TRUE
      objects <- rbind(objects, row)
    }
    objects <- objects[order(objects$object), ]
    rownames(objects) <- NULL
  }

  object_coefficients <- NULL
  if (!is.null(object_design)) {
    bhat <- drop(solve(crossprod(Bmat), crossprod(Bmat, beta - beta0))) /
      object_scale
    bse <- sqrt(pmax(diag(covth)[seq_len(nb)], 0)) / object_scale
    if (!cluster_inference) bse[] <- NA_real_
    stat <- .wald_ratio(bhat, bse)
    ref_df <- if (!cluster_inference) NA_real_ else if (!is.null(jd))
      max(nc - 1L, 1L) else Inf
    prob <- 2 * stats::pt(-abs(stat), df = ref_df)
    object_coefficients <- data.frame(
      term = object_parameter_names, estimate = bhat, se = bse,
      t = stat, df = ref_df, p = prob, stringsAsFactors = FALSE)
    object_coefficients$p_adj <- .p_adjust_family(
      object_coefficients$p, method = "holm")
    rownames(object_coefficients) <- object_coefficients$term
  }

  judges <- NULL
  if (!is.null(jd)) {
    ju <- sort(unique(jd))
    jfit <- lapply(ju, function(j) pool(jd == j))
    judges <- data.frame(judge = ju,
                         n = vapply(jfit, `[[`, 0, "n"),
                         infit_ms = vapply(jfit, `[[`, 0, "infit_ms"),
                         outfit_ms = vapply(jfit, `[[`, 0, "outfit_ms"),
                         fit_resid = vapply(jfit, `[[`, 0, "fit_resid"),
                         df_fit = vapply(jfit, `[[`, 0, "df"))
    rownames(judges) <- NULL
  }

  # pairwise goodness of fit on the oriented mean response
  key <- ifelse(ia < ib, paste(ia, ib), paste(ib, ia))
  x_lo <- ifelse(ia < ib, x, m - x)
  E_lo <- ifelse(ia < ib, mo$E, m - mo$E)
  n_pair <- tapply(w, key, sum)
  obs_m <- tapply(w * x_lo, key, sum) / n_pair
  exp_m <- tapply(w * E_lo, key, sum) / n_pair
  v_pair <- tapply(w * mo$V, key, sum)
  zp <- tapply(w * (x_lo - E_lo), key, sum) / sqrt(pmax(v_pair, 1e-12))
  idx <- do.call(rbind, strsplit(names(n_pair), " "))
  pairs <- data.frame(object_a = objs[as.integer(idx[, 1])],
                      object_b = objs[as.integer(idx[, 2])],
                      n = as.numeric(n_pair),
                      obs_mean = as.numeric(obs_m),
                      exp_mean = as.numeric(exp_m),
                      residual = as.numeric(zp),
                      chisq = as.numeric(zp)^2)
  rownames(pairs) <- NULL
  used <- pairs$n >= 2
  # df: informative pairs minus ALL estimated parameters (locations,
  # thresholds, and any position/dependence covariates); when nothing
  # testable remains the total is NA, not a manufactured df = 1
  total_chisq <- sum(pairs$chisq[used])
  total_df <- sum(used) - np
  if (total_df < 1L) { total_chisq <- NA_real_; total_df <- NA_integer_ }
  # The note explains a withheld total_p, so only a caller that reports that
  # probability should carry it.
  if (!is.null(jd) && withheld_notes)
    notes <- c(notes, paste(
      "the pairwise chi-square probability is withheld because its row-based",
      "reference does not model within-judge dependence; the statistic remains",
      "descriptive, and fit_bootstrap() provides calibration under the fitted",
      "response model"))
  osi <- if (!cluster_inference)
    list(PSI = NA_real_, separation = NA_real_, strata = NA_real_,
         var_theta = NA_real_, mean_error_var = NA_real_, n = 0L)
  else .psi(objects$location, objects$se)

  # two categories ARE the dichotomous conditional model, so an m == 1 fit is
  # presented in dichotomous terms: the score is the win count and the mean
  # polytomous responses are win proportions (one estimator serves both routes)
  if (m == 1L) {
    names(objects)[names(objects) == "score"] <- "wins"
    names(pairs)[names(pairs) == "obs_mean"] <- "obs_prop"
    names(pairs)[names(pairs) == "exp_mean"] <- "exp_prop"
  }

  # the per-comparison history covariates, so the dependence effects can be
  # interrogated: a comparison is informative for an effect when its covariate
  # is non-zero (the two objects' histories differ)
  dependence_data <- if (is.null(Zfull)) NULL else {
    dd <- data.frame(judge = if (is.null(jd)) NA_character_ else jd,
                     order = if (is.null(ord)) NA_real_ else ord,
                     object_a = a, object_b = b, response = x, weight = w,
                     stringsAsFactors = FALSE)
    # whichever covariate columns the design carried (exposure, carry_over,
    # position), added by name so any subset works
    for (cn in colnames(Zfull)) dd[[cn]] <- Zfull[, cn]
    dd <- dd[order(dd$judge, dd$order), ]; rownames(dd) <- NULL; dd
  }
  if (!cluster_inference) {
    if (!is.null(components)) components$se[] <- NA_real_
    cov_beta[,] <- NA_real_
    if (length(anch_idx)) {
      # Fixed anchors have zero covariance with themselves and every free
      # estimate. Only the unsupported free-by-free block is unavailable.
      cov_beta[anch_idx, ] <- 0
      cov_beta[, anch_idx] <- 0
    }
  }
  observed_comparisons <- data.frame(
    object_a = a0, object_b = b0, response = x0, weight = w0,
    judge = if (is.null(jd0)) NA_character_ else jd0,
    stringsAsFactors = FALSE)
  if (!is.null(ord0)) observed_comparisons$order <- ord0

  out <- list(objects = objects, thresholds = thresholds,
              components = components, thr_structure = thr,
              dependence = dependence, dependence_data = dependence_data,
              pairs = pairs,
              judges = judges, m = m, categories = cats,
              total_chisq = total_chisq, total_df = total_df,
              total_p = if (isTRUE(converged) && isTRUE(cluster_inference) &&
                             is.null(jd) && is.finite(total_df))
                pchisq(total_chisq, total_df, lower.tail = FALSE) else NA_real_,
              osi = osi, loglik = loglik, iterations = it,
              converged = converged, n_comparisons = n_rows,
              clustered = !is.null(jd), cov_beta = cov_beta, cl = cl_info,
              location_design = Bmat_original, location_offset = beta0,
              object_design = if (is.null(object_design)) NULL else Bmat_original,
              object_offset = if (is.null(object_design)) NULL else beta0,
              object_coefficients = object_coefficients,
              sensitivity = H * outer(c(object_scale, rep(1, q + pz)),
                                      c(object_scale, rep(1, q + pz))),
              cov_parameters = if (isTRUE(converged) &&
                                    isTRUE(cluster_inference))
                covth / outer(c(object_scale, rep(1, q + pz)),
                              c(object_scale, rep(1, q + pz))) else {
                z <- covth; z[,] <- NA_real_; z
              },
              fitted_prob = mo$P,
              refit_spec = list(
                thresholds = thr, anchors = anch,
                maxit = maxit, tol = tol,
                has_order = !is.null(ord),
                position = !is.null(Zfull) && "position" %in% colnames(Zfull),
                object_design = object_design),
              comparisons = {
                cmp <- data.frame(object_a = a, object_b = b,
                                  response = x, weight = w,
                                  judge = if (is.null(jd))
                                    NA_character_ else jd)
                if (!is.null(ord)) cmp$order <- ord
                # row-aligned history covariates, so downstream analyses
                # (btl_dif) can hold the fitted dependence effects fixed
                if (!is.null(Zfull))
                  for (cn in colnames(Zfull)) cmp[[cn]] <- Zfull[, cn]
                cmp
              },
              observed_comparisons = observed_comparisons,
              # Row-aligned bookkeeping for the independent comparison
              # units.  Ordinarily each stored row is one unit repeated
              # `weight` times.  A half-scored tie instead has two stored
              # rows in the same unit; retaining that allocation lets model
              # comparison distinguish one tie from two independent,
              # opposing judgements while recognising count compression as
              # a change of representation only.
              independent_allocation = if (is.null(jd)) data.frame(
                unit = as.character(row_cluster),
                replicates = as.numeric(row_replicates),
                stringsAsFactors = FALSE) else NULL,
              anchors = anch,
              notes = notes)
  out <- .tag_tables(out)
  class(out) <- "rasch_btl"
  out
}

#' Plot polytomous-comparison category curves
#'
#' For a polytomous paired-comparison fit, the probability of each response
#' category as a function of the location difference
#' \code{beta_a - beta_b}, with the symmetric threshold structure marked.
#' The display is the paired-comparison counterpart of the category
#' probability curves of a polytomous item.
#' A fitted position effect is included for object a presented first. For a
#' history-dependent fit, the horizontal axis is the full linear predictor,
#' including position and history terms, rather than the object difference alone.
#'
#' @param fit A polytomous fit from \code{\link{btl}} (with \code{response}).
#' @param grid Difference grid, in logits; the full linear-predictor grid for
#'   history-dependent fits.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
#' P <- vapply(seq_len(nrow(d)), function(r)
#'   item_moments(beta[d$a[r]] - beta[d$b[r]], c(-1, 0, 1))$P, numeric(4))
#' d$grade <- apply(P, 2, function(p) sample(0:3, 1, prob = p))
#' plot_btl_categories(btl(d, "a", "b", response = "grade"))
#' @export
plot_btl_categories <- function(fit, grid = seq(-4, 4, 0.05)) {
  .check_btl_display_fit(fit)
  .check_grid(grid)
  if (is.null(fit$m) || fit$m < 2L)
    stop("category curves need a polytomous fit (three or more categories)")
  tau <- fit$thresholds$tau
  history <- .btl_information_history(fit)
  position <- if (history) 0 else .btl_information_position(fit)
  P <- vapply(grid, function(d) item_moments(d + position, tau)$P,
              numeric(fit$m + 1L))
  op <- .rr_canvas(range(grid), c(0, 1),
                   if (history) "Full linear predictor (logits)" else
                     "Location difference (logits)", "Category probability")
  on.exit(par(op))
  abline(v = tau - position, lty = 3, col = .rr$soft)
  labs <- if (!is.null(fit$categories)) fit$categories else
    as.character(0:fit$m)
  for (cat in 0:fit$m)
    lines(grid, P[cat + 1L, ], lwd = 2.6,
          col = .rr$pal[cat %% length(.rr$pal) + 1L])
  .rr_legend("right", labs, lwd = 2.6,
             col = .rr$pal[(0:fit$m) %% length(.rr$pal) + 1L])
  invisible(NULL)
}

#' Plot an object characteristic curve
#'
#' Plots the expected response for one object against opponent location, with
#' the observed mean response against each sufficiently observed opponent.
#' For dichotomous fits the curve is the win probability; for ordered fits it
#' is the expected response.
#' With fitted position or history effects, model values are weighted means
#' of the fitted expectations for each observed opponent, joined in location
#' order. Judge-group overlays use each group's own comparison context.
#'
#' @param fit An object from \code{\link{btl}}.
#' @param object Object name.
#' @param group Optional judge grouping for a DIF overlay: either one value
#'   per comparison row of \code{fit$comparisons} or a vector named by
#'   every judge in the fit. Observed means are then drawn separately per group, as
#'   \code{\link{plot_icc}} draws person groups.
#' @param grid Opponent-location grid, in logits.
#' @param min_n An opponent's observed point is drawn only when the object
#'   (or, in the grouped display, that judge group) met it at least this many
#'   times; sparser pairs from incomplete or unbalanced designs are omitted.
#' @return Called for its plotting side effect; invisibly the names of the
#'   opponents drawn (the ungrouped display), or \code{NULL} for the grouped
#'   display.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
#' p <- plogis(beta[d$a] - beta[d$b])
#' d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' plot_btl_icc(btl(d, "a", "b", winner = "win"), "C")
#' @export
plot_btl_icc <- function(fit, object, group = NULL, grid = NULL,
                         min_n = 10) {
  .check_btl_display_fit(fit)
  if (missing(object) || length(object) != 1L || is.na(object) ||
      !is.atomic(object) || !is.null(dim(object)))
    stop("`object` must name exactly one object")
  object <- .role_text_values(object)
  if (!nzchar(object)) stop("`object` must name exactly one object")
  if (!is.null(grid)) .check_grid(grid)
  min_n <- .check_whole(min_n, "min_n", 1)
  if (inherits(fit, "rasch_btl_efrm")) {
    if (!is.null(group))
      stop("judge-group DIF curves are not defined after a frame adjustment; ",
           "inspect the fitted panel-by-set frame curves, or use the ",
           "equal-unit fit for judge-group DIF")
    return(.plot_btl_efrm_icc(fit, object, grid = grid, min_n = min_n))
  }
  ob_all <- fit$objects
  if (!object %in% ob_all$object) stop("no such object: ", object)
  if (isTRUE(ob_all$extreme[ob_all$object == object]))
    .refuse(object, " was set aside at a response boundary; it has no ",
            "fitted curve. Its reported location is an extrapolation ",
            "for display only")
  ob <- if ("extreme" %in% names(ob_all))
    ob_all[!(ob_all$extreme %in% TRUE), , drop = FALSE] else ob_all
  m <- if (is.null(fit$m)) 1L else fit$m
  tau <- if (!is.null(fit$thresholds)) fit$thresholds$tau else numeric(1)
  b_o <- ob$location[ob$object == object]
  if (is.null(grid)) {
    rng <- range(ob$location) + c(-1, 1)
    grid <- seq(rng[1], rng[2], length.out = 201)
  }
  Ecurve <- vapply(grid, function(t) item_moments(b_o - t, tau)$E, 0)
  cm <- fit$comparisons
  conditional_curve <- is.data.frame(fit$dependence) && nrow(fit$dependence) > 0L
  expected_row <- if (conditional_curve) .btl_fitted_score_moments(fit)$E else NULL
  gv <- NULL
  if (!is.null(group)) {
    gv <- if (!is.null(names(group))) {
      if (anyNA(names(group)) || any(!nzchar(trimws(names(group)))))
        stop("the group map must use non-missing judge names")
      names(group) <- .role_text_values(names(group))
      if (anyDuplicated(names(group)))
        stop("duplicate judge(s) in the group map after trimming: ",
             paste(unique(names(group)[duplicated(names(group))]),
                   collapse = ", "),
             "; each judge may carry one value")
      observed <- unique(cm$judge[!is.na(cm$judge)])
      absent <- setdiff(observed, names(group))
      if (length(absent))
        stop("judge(s) missing from the group map: ",
             paste(utils::head(absent, 5L), collapse = ", "))
      # a map built from the source data legitimately names judges the fit
      # set aside; they take no part in the curve, so they are ignored
      unname(.role_text_values(group)[match(cm$judge, names(group))])
    } else {
      if (length(group) != nrow(cm))
        stop("`group` must have one entry per comparison or be named by judge")
      .role_text_values(group)
    }
    gv[!is.na(gv) & !nzchar(gv)] <- NA_character_
    if (!any(!is.na(gv)))
      stop("`group` has no observed judge-group values in the fitted comparisons")
  }
  sel_a <- cm$object_a == object
  sel_b <- cm$object_b == object
  opp <- c(cm$object_b[sel_a], cm$object_a[sel_b])
  resp <- c(cm$response[sel_a], m - cm$response[sel_b])
  expected <- if (conditional_curve)
    c(expected_row[sel_a], m - expected_row[sel_b]) else NULL
  wt <- c(cm$weight[sel_a], cm$weight[sel_b])
  gg <- if (is.null(gv)) NULL else c(gv[sel_a], gv[sel_b])
  keep <- opp %in% ob$object
  if (!is.null(gg)) keep <- keep & !is.na(gg)
  opp <- opp[keep]; resp <- resp[keep]; wt <- wt[keep]
  if (conditional_curve) expected <- expected[keep]
  if (!is.null(gg)) {
    gg <- gg[keep]
    if (!length(gg))
      stop("`group` has no observed judge-group values for object '", object,
           "' against a calibrated opponent")
  }
  obs <- data.frame(
    opponent = tapply(opp, opp, `[`, 1),
    loc = ob$location[match(names(tapply(wt, opp, sum)), ob$object)],
    mean = as.numeric(tapply(wt * resp, opp, sum) / tapply(wt, opp, sum)),
    n = as.numeric(tapply(wt, opp, sum)))
  if (conditional_curve)
    obs$model <- as.numeric(tapply(wt * expected, opp, sum) / tapply(wt, opp, sum))
  op <- .rr_canvas(range(grid), c(0, m), "Opponent location (logits)",
                   if (m == 1L) "Probability preferred" else
                     "Expected polytomous response",
                   sprintf("%s  (location %.3f)", object, b_o))
  on.exit(par(op))
  if (!conditional_curve) lines(grid, Ecurve, lwd = 3, col = .rr$ink)
  abline(v = b_o, lty = 3, col = .rr$soft)
  if (is.null(gg)) {
    # a comparator is shown only when the object met it enough times for the
    # observed proportion to be informative; sparser pairs (incomplete or
    # unbalanced designs) are omitted rather than plotted as noise
    shown <- obs[obs$n >= min_n, , drop = FALSE]
    n_omit <- nrow(obs) - nrow(shown)
    if (conditional_curve && nrow(shown)) {
      oo <- order(shown$loc)
      lines(shown$loc[oo], shown$model[oo], lwd = 3, col = .rr$ink)
      points(shown$loc, shown$model, pch = 3, col = .rr$ink)
    }
    if (nrow(shown)) {
      points(shown$loc, shown$mean, pch = 21, bg = .rr$blue,
             col = "white", cex = 1.5, lwd = 1.2)
      text(shown$loc, shown$mean, shown$opponent, pos = 3, offset = 0.45,
           cex = 0.72, col = .rr$soft)
    }
    .rr_legend("topright",
               c("Model", "Observed",
                 if (n_omit)
                   sprintf("%d omitted (< %d comparisons)", n_omit, min_n)),
               lwd = c(3, NA, if (n_omit) NA),
               pch = c(NA, 21, if (n_omit) NA),
               pt.bg = c(NA, .rr$blue, if (n_omit) NA),
               col = c(.rr$ink, "white", if (n_omit) .rr$soft),
               pt.cex = 1.3)
  } else {
    # the graphical DIF display: per-opponent means drawn separately for
    # each judge group, as plot_icc draws person groups. A group's point for
    # an opponent is shown only where that group met it enough times
    levs <- sort(unique(gg))
    for (li in seq_along(levs)) {
      sel <- gg == levs[li]
      nn <- tapply(wt[sel], opp[sel], sum)
      om <- tapply(wt[sel] * resp[sel], opp[sel], sum) / nn
      om <- om[nn >= min_n]
      ol <- ob$location[match(names(om), ob$object)]
      colr <- .rr$pal[(li - 1L) %% length(.rr$pal) + 1L]
      oo <- order(ol)
      if (conditional_curve) {
        em <- tapply(wt[sel] * expected[sel], opp[sel], sum) / nn
        em <- em[names(om)]
        lines(ol[oo], em[oo], col = colr, lwd = 2.2)
        points(ol, em, pch = 3, col = colr)
      }
      lines(ol[oo], om[oo], col = colr, lwd = 1.4, lty = 3)
      points(ol, om, pch = 21, bg = colr, col = "white", cex = 1.4,
             lwd = 1.1)
    }
    if (conditional_curve) {
      .rr_legend("topright", levs, lwd = rep(2.2, length(levs)),
                 col = .rr$pal[(seq_along(levs) - 1L) %% length(.rr$pal) + 1L])
      .rr_legend("bottomleft", c("Model", "Observed"),
                 lwd = c(2.2, 1.4), lty = c(1, 3), pch = c(3, 21),
                 col = c(.rr$ink, .rr$ink), pt.bg = c(NA, .rr$blue))
    } else .rr_legend("topright", c("Model", levs), lwd = c(3, rep(1.4, length(levs))),
               lty = c(1, rep(3, length(levs))),
               pch = c(NA, rep(21, length(levs))),
               pt.bg = c(NA, .rr$pal[(seq_along(levs) - 1L) %%
                                     length(.rr$pal) + 1L]),
               col = c(.rr$ink, .rr$pal[(seq_along(levs) - 1L) %%
                                       length(.rr$pal) + 1L]), pt.cex = 1.2)
  }
  # the opponents actually drawn (ungrouped display), for inspection and tests
  invisible(if (is.null(gg)) obs[obs$n >= min_n, "opponent"] else NULL)
}

#' Plot a within-judge dependence effect
#'
#' The graphical display of a paired-comparison dependence effect, the
#' counterpart of the DIF characteristic curve. For every comparison the
#' departure of the observed response from what the object locations alone
#' predict is taken, and the contribution of the \emph{other} dependence
#' effect is removed (a partial-residual display); these departures are then
#' averaged in bins of the effect's own history covariate and plotted against
#' it, with the model's fitted contribution overlaid. Observed points that
#' rise with the covariate along the fitted line are the effect the
#' coefficient summarises; a flat, scattered cloud means the estimate rests on
#' little. Only the informative comparisons (a non-zero covariate: the two
#' objects' histories differ) carry the effect, and the count in each bin is
#' printed so a thin exposure tail is visible.
#'
#' @param fit An object from \code{\link{btl}} fitted with an \code{order}
#'   column, so \code{fit$dependence_data} is present.
#' @param effect Which effect to display: \code{"exposure"} (the seen-before
#'   advantage, the default) or \code{"carry_over"} (response dependence).
#' @param bins Number of covariate bins for the continuous carry-over display;
#'   exposure takes its three natural levels (-1, 0, +1).
#' @return Called for its plotting side effect; invisibly a data frame of the
#'   binned covariate value, observed and fitted departure, and bin count.
#' @references Davidson, R. R., & Beaver, R. J. (1977). On extending the
#'   Bradley-Terry model to incorporate within-pair order effects.
#'   \emph{Biometrics}, 33(4), 693-702.
#' @examples
#' set.seed(1)
#' beta <- c(A = -0.8, B = -0.2, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
#' d$judge <- sample(sprintf("J%02d", 1:8), nrow(d), TRUE)
#' d <- d[order(d$judge), ]; d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' f <- btl(d, "a", "b", winner = "win", judge = "judge", order = "t")
#' plot_btl_dependence(f, "carry_over")
#' @export
plot_btl_dependence <- function(fit, effect = c("exposure", "carry_over"),
                                bins = 6) {
  .check_btl_display_fit(fit)
  if (missing(effect)) effect <- "exposure"
  if (!is.character(effect) || length(effect) != 1L || is.na(effect))
    stop("`effect` must be one of \"exposure\" or \"carry_over\"")
  effect <- match.arg(effect)
  bins <- .check_whole(bins, "bins", 2)
  dd <- fit$dependence_data
  if (is.null(dd))
    stop("no dependence data: fit btl() with an `order` (and `judge`) column")
  if (is.null(fit$dependence))
    stop("no dependence effect was estimable (see fit$notes): ",
         paste(fit$notes, collapse = "; "))
  eff <- fit$dependence[fit$dependence$effect == effect, ]
  if (!nrow(eff))
    stop("effect not estimated (see fit$notes): ", effect)
  bl <- setNames(fit$objects$location, fit$objects$object)
  m <- if (is.null(fit$m)) 1L else fit$m
  tau <- if (!is.null(fit$thresholds)) fit$thresholds$tau else numeric(1)
  dep <- setNames(fit$dependence$estimate, fit$dependence$effect)
  phi <- if ("exposure" %in% names(dep)) dep[["exposure"]] else 0
  psi <- if ("carry_over" %in% names(dep)) dep[["carry_over"]] else 0
  # the first-position advantage, when fitted, shifts every comparison's
  # linear predictor; the expected-score map is nonlinear, so leaving it out
  # of the baseline would displace both the observed and the fitted departure
  pos <- if ("position" %in% names(dep)) dep[["position"]] else 0
  pos_cov <- if (is.null(dd$position)) 0 else dd$position
  Emom <- function(v) vapply(v, function(t) item_moments(t, tau)$E, 0)

  # partial-residual display: hold everything but this effect at its fitted
  # value, so the plotted departure isolates this covariate's contribution
  base <- unname(bl[dd$object_a] - bl[dd$object_b])
  lin_full <- base + phi * dd$exposure + psi * dd$carry_over + pos * pos_cov
  cov <- dd[[effect]]
  coef_e <- if (effect == "exposure") phi else psi
  E_other <- Emom(lin_full - coef_e * cov)
  obs <- dd$response - E_other                 # observed departure
  fit_lift <- Emom(lin_full) - E_other         # model-fitted departure

  wt <- if (is.null(dd$weight)) rep(1, nrow(dd)) else dd$weight
  if (effect == "exposure") {
    g <- factor(cov, levels = sort(unique(cov))); xb <- as.numeric(levels(g))
  } else {
    bins <- max(2L, as.integer(bins))
    br <- unique(stats::quantile(cov, seq(0, 1, length.out = bins + 1L),
                                 na.rm = TRUE))
    # a heavy mass point (many unseen pairs at 0) can collapse the quantile
    # breaks; fall back to equal-width bins rather than per-value singletons
    if (length(br) <= 2L)
      br <- unique(pretty(range(cov, na.rm = TRUE), bins))
    g <- if (length(br) > 2L) cut(cov, br, include.lowest = TRUE) else
      factor(cov)
    xw <- tapply(wt * cov, g, sum); xb <- as.numeric(xw / tapply(wt, g, sum))
  }
  # count-weighted rows stand for several comparisons: weighted bin means,
  # and the printed n is the number of comparisons, not rows
  ob <- as.numeric(tapply(wt * obs, g, sum) / tapply(wt, g, sum))
  fb <- as.numeric(tapply(wt * fit_lift, g, sum) / tapply(wt, g, sum))
  nb <- as.numeric(tapply(wt, g, sum))
  keep <- !is.na(xb) & !is.na(ob) & nb > 0
  xb <- xb[keep]; ob <- ob[keep]; fb <- fb[keep]; nb <- nb[keep]
  oo <- order(xb); xb <- xb[oo]; ob <- ob[oo]; fb <- fb[oo]; nb <- nb[oo]

  lab <- gsub("_", "-", effect)
  yl <- range(c(ob, fb, 0), na.rm = TRUE); yl <- yl + c(-1, 1) * 0.08 * (diff(yl) + 1e-6)
  xr <- range(xb); xl <- xr + c(-1, 1) * (0.12 * diff(xr) + 0.05)
  p_label <- if ("p_adj" %in% names(eff) && is.finite(eff$p_adj))
    paste0("Holm p = ", .fmt_p(eff$p_adj)) else "Holm p unavailable"
  op <- .rr_canvas(xl, yl, sprintf("%s covariate", lab),
                   if (m == 1L) "Observed - expected win probability"
                   else "Observed - expected response",
                   sprintf("%s dependence: %.3f logits (SE %.3f, %s)",
                           lab, eff$estimate, eff$se, p_label),
                   grid_x = TRUE)
  on.exit(par(op))
  abline(h = 0, lty = 3, col = .rr$soft)
  lines(xb, fb, lwd = 2.6, col = .rr$ink)
  points(xb, ob, pch = 21, bg = .rr$blue, col = "white", cex = 1.7, lwd = 1.2)
  text(xb, ob, nb, pos = 3, offset = 0.6, cex = 0.65, col = .rr$soft)
  .rr_legend("topleft", c("Model", "Observed"),
             lwd = c(2.6, NA), pch = c(NA, 21), pt.bg = c(NA, .rr$blue),
             col = c(.rr$ink, "white"), pt.cex = 1.3)
  invisible(data.frame(covariate = xb, observed = ob, fitted = fb, n = nb))
}

# ---------------------------------------------------------------------------
# DIF for paired comparisons: object-by-judge-group interaction. Judge
# severity cancels within a comparison, so group membership can only reach
# the measurement through object-specific preference - which is DIF, tested
# here by the package's two standard routes: a residual analysis of
# variance (group crossed with opponent-strength bands, the class-interval
# analogue) and resolved locations in logits (the object split into one
# copy per judge group inside a joint fit; Dittrich, Hatzinger &
# Katzenbeisser 1998 model these judge-covariate-by-object terms in the
# log-linear frame).
# ---------------------------------------------------------------------------

# Denominator reference for a resolved BTL-DIF contrast. The ordinary
# two-cell comparison retains its simulation-calibrated Welch reference.
# With more cells, the individual cluster contributions needed for a general
# Satterthwaite calculation are not recoverable from cell effective counts;
# anchor the reference to the weakest contributing cell instead.
.btl_dif_contrast_df <- function(effective_judges) {
  if (length(effective_judges) == 2L)
    return(max(sum(effective_judges) - 2, 1))
  max(min(effective_judges) - 1, 1)
}

#' DIF analysis for paired comparisons
#'
#' Tests whether object locations differ across groups of judges. Several
#' judge factors can be fitted jointly, with optional factor-by-factor
#' interactions. Uniform DIF is a judge-factor effect; non-uniform DIF is its
#' interaction with opponent-strength band.
#'
#' @details
#' Judges are the independent units. For each object, oriented residuals are
#' aggregated to one weighted mean per judge and opponent band. A split-plot
#' analysis then tests judge factors between judges and band effects within
#' judges. Each factor level requires at least two judges. Confirmatory Wald
#' tests are available only when the base fit supplies a valid judge-clustered
#' covariance. The base paired-comparison calibration must have converged.
#' BTL-EFRM fits are not accepted: the ordinary residual and resolution
#' models do not contain the fitted panel and set units.
#'
#' A significant uniform term is followed by a joint refit in which the object
#' has one location per cell of the complete judge-factor design. Main-effect
#' magnitudes average these cells equally over the other factors. Interaction
#' magnitudes are differences between differences, with the corresponding
#' higher-order tensor contrast beyond two factors. A contributing cell needs
#' at least eight effective judges for inference; otherwise its location and
#' contrasts remain descriptive. Higher-order terms supersede their component
#' terms. Two-cell contrasts retain the Welch reference used by the ordinary
#' pairwise comparison. Contrasts spanning more than two fitted cells use the
#' effective-judge count in their least-supported cell as a conservative
#' denominator reference.
#' Models fitted with \code{order} retain the exposure and carry-over effects
#' in both the residual analysis and refit.
#' Between-judge tests use HC3 covariance so unequal comparison workloads do
#' not impose equal precision on judge means. Omnibus probabilities require
#' at least eight judges and eight effective judges in every factor cell.
#' Holm adjustment is the default; \code{"BH"} remains available for
#' false-discovery-rate screening. A reported object-by-term test remains in
#' the adjustment family when its probability is unavailable.
#'
#' Objects are resolved one at a time against the common locations of the
#' remaining objects. With DIF in several objects, this can induce compensating
#' apparent DIF in invariant objects (Andrich and Hagquist 2012, 2015).
#' An externally anchored object is not resolved: fixing each of its copies at
#' the same anchor would define their difference as zero. Anchors on the other
#' objects are retained in the joint refit. If a resolved-location covariance
#' is unavailable or not positive semidefinite, the locations and differences
#' remain descriptive but their uncertainty and tests are withheld.
#' A note raised by a resolution refit is reported against that object.
#' Engine notes explaining a statistic this analysis never reports -- the
#' refit's pairwise chi-square probability and its dependence table's
#' carry-over probability -- are not carried; the base fit retains them.
#'
#' @param fit An ordinary paired-comparison fit from \code{\link{btl}}.
#' @param factors One judge factor, or a named list containing several. Each
#'   factor may have one value per comparison row or be a vector named by
#'   every judge in the fit.
#' @param objects Objects to test; all by default.
#' @param effects \code{"main"} (default) models several factors additively
#'   (each factor's main effect and its band interaction); \code{"factorial"}
#'   also crosses the factors with one another.
#' @param p_adjust Multiplicity adjustment over all object-by-term tests;
#'   the resolved-contrast probabilities are adjusted separately in one pool
#'   over all objects, terms, and contrasts.
#' @param alpha Significance level for adjusted probabilities.
#' @param flag_logits Absolute resolved difference flagged as practically
#'   significant.
#' @param min_n Term cells with fewer comparisons involving the object are
#'   dropped from its resolution, with a note.
#' @param maxit,tol Newton controls for the resolution refits.
#' @return A list of class \code{"rasch_btl_dif"}: \code{summary} (one row per
#'   object and group term with the uniform F, adjusted p and partial
#'   eta-squared -- the term itself -- the non-uniform ones -- the term
#'   crossed with the opponent band -- plus \code{uniform_DIF},
#'   \code{nonuniform_DIF} and \code{superseded} flags); \code{terms} (the
#'   full per-object analysis-of-variance table, including its raw and
#'   effective judge support); \code{levels} (resolved
#'   location, SE, comparison count, judge count and effective judge count per
#'   object, term and complete-design cell); \code{sizes} (per object, term
#'   and marginal or interaction contrast: difference in logits, judge
#'   support for both sides, SE, t, degrees of
#'   freedom, adjusted p, significance and practical flags); \code{effects},
#'   \code{factors}, \code{alpha}, \code{p_adjust}, \code{flag_logits}, and
#'   \code{notes}. \code{size_family_n} records the complete planned
#'   resolved-contrast family, including unavailable comparisons.
#'   \code{summary_factors} retains the factor membership of each displayed
#'   term.
#' @references Andrich, D., & Hagquist, C. (2012). Real and artificial
#'   differential item functioning. \emph{Journal of Educational and
#'   Behavioral Statistics}, 37(3), 387-416.
#'
#'   Dittrich, R., Hatzinger, R., & Katzenbeisser, W. (1998).
#'   Modelling the effect of subject-specific covariates in paired
#'   comparison studies with an application to university rankings.
#'   \emph{Journal of the Royal Statistical Society C}, 47(4), 511-525.
#'
#'   MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent
#'   covariance matrix estimators with improved finite sample properties.
#'   \emph{Journal of Econometrics}, 29(3), 305--325.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 100), b = rep(pr[, 2], each = 100),
#'                 judge = sample(sprintf("J%02d", 1:20), 600, TRUE))
#' shift <- ifelse(d$judge %in% sprintf("J%02d", 1:10) & d$a == "C", 0.9,
#'          ifelse(d$judge %in% sprintf("J%02d", 1:10) & d$b == "C", -0.9, 0))
#' p <- plogis(beta[d$a] - beta[d$b] + shift)
#' d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' f <- btl(d, "a", "b", winner = "win", judge = "judge")
#' grp <- setNames(rep(c("g1", "g2"), each = 10), sprintf("J%02d", 1:20))
#' btl_dif(f, grp, objects = "C")
#' @export
btl_dif <- function(fit, factors, objects = NULL,
                    effects = c("main", "factorial"),
                    p_adjust = "holm", alpha = 0.05, flag_logits = 0.5,
                    min_n = 20, maxit = 60, tol = 1e-8) {
  .check_dif_args(alpha, p_adjust, flag_logits, min_n)
  .check_controls(maxit, tol)
  if (!inherits(fit, "rasch_btl"))
    stop("btl_dif needs a paired-comparison fit from btl()")
  if (inherits(fit, "rasch_btl_efrm"))
    stop("judge-group DIF is not defined after a BTL-EFRM frame adjustment; ",
         "the ordinary residual and resolution models do not carry the ",
         "fitted panel and set units. Examine frame-specific fit, or fit the ",
         "equal-frame BTL model for an explicitly conditional DIF analysis")
  if (inherits(fit, "rasch_btl_explanatory"))
    stop("judge-group DIF is not defined for an explanatory comparison fit: ",
         "the object locations are functions of their covariates, so the ",
         "resolved copies would either be forced equal by the design -- ",
         "defining their DIF as zero -- or estimated without it. Fit the ",
         "unrestricted btl() model for a DIF analysis")
  if (!isTRUE(fit$converged))
    stop("the paired-comparison calibration did not converge; DIF inference is unavailable")
  effects <- match.arg(effects)
  cm <- fit$comparisons
  if (is.null(cm)) stop("the fit carries no comparisons")
  if (all(is.na(cm$judge)))
    stop("btl_dif needs judge identifiers: judges are the independent units")
  if (!isTRUE(fit$cl$inference_available))
    stop("the base fit does not support cluster-robust inference; btl_dif ",
         "requires at least 10 judges, at least 8 effective judges, and ",
         "more effective judges than fitted parameters, with full empirical ",
         "score-covariance rank; spread comparisons across judges or simplify ",
         "the fitted model, and inspect fit$notes for the failed condition")
  # a single grouping is promoted to a one-factor list; several judge factors
  # are modelled jointly (main effects by default, interactions if asked)
  if (!is.list(factors)) factors <- list(group = factors)
  if (!length(factors))
    stop("`factors` must name at least one judge factor")
  if (is.null(names(factors)))
    names(factors) <- paste0("factor", seq_along(factors))
  else if (anyNA(names(factors)) || any(!nzchar(trimws(names(factors)))))
    stop("every judge factor needs a non-empty name")
  if (anyDuplicated(names(factors)))
    stop("duplicate factor name(s): ",
         paste(unique(names(factors)[duplicated(names(factors))]),
               collapse = ", "))
  fnames <- names(factors)
  unfitted <- character(0)
  gvs <- lapply(factors, function(g) {
    if (!is.atomic(g) || !length(g) || !is.null(dim(g)))
      stop("each judge factor must be a plain vector, optionally named by judge",
           call. = FALSE)
    if (!is.null(names(g))) {
      if (anyNA(names(g)) || any(!nzchar(trimws(names(g)))))
        stop("named judge factors must use non-missing judge names")
      names(g) <- .role_text_values(names(g))
      if (anyDuplicated(names(g)))
        stop("duplicate judge(s) in a named factor after trimming: ",
             paste(unique(names(g)[duplicated(names(g))]), collapse = ", "),
             "; each judge may carry one value")
      observed <- unique(cm$judge[!is.na(cm$judge)])
      absent <- setdiff(observed, names(g))
      extra <- setdiff(names(g), observed)
      if (length(absent))
        stop("judge(s) missing from a named factor: ",
             paste(utils::head(absent, 5L), collapse = ", "),
             "; every judge needs an explicit factor entry")
      # a map built from the source data legitimately names judges the fit
      # set aside -- one who only ever tied, or whose rows were dropped --
      # and they take no part in the analysis, so they are ignored and
      # reported rather than refused
      if (length(extra)) unfitted <<- union(unfitted, extra)
      unname(.role_text_values(g)[match(cm$judge, names(g))])
    } else {
      if (length(g) != nrow(cm))
        stop("each factor needs one value per comparison or names by judge")
      .role_text_values(g)
    }
  })
  gvs <- lapply(gvs, function(g) {
    g[!is.na(g) & !nzchar(g)] <- NA_character_
    g
  })
  # judge-group DIF tests judge ATTRIBUTES: a row-wise factor that varies
  # within a judge has no judge-level value, and the judge-level analysis
  # would silently take whichever row came first
  for (j in seq_along(gvs)) {
    nvar <- tapply(gvs[[j]], cm$judge, function(v)
      length(unique(v[!is.na(v)])))
    if (any(nvar > 1L, na.rm = TRUE))
      stop("factor '", fnames[j], "' varies within judge(s) ",
           paste(names(nvar)[which(nvar > 1L)], collapse = ", "),
           ": judge-group DIF needs judge-constant factors")
  }
  ok <- Reduce(`&`, lapply(gvs, function(g) !is.na(g)))
  safe <- paste0("f", seq_along(fnames))            # syntactic stand-ins
  op <- if (effects == "factorial") " * " else " + "
  tvars <- function(t) strsplit(t, ":", fixed = TRUE)[[1]]
  # Map whole stand-in tokens back to the nominated factor names for messages
  # and final tables. A user factor named "band" remains distinct from the
  # internal opponent-strength band.
  relab <- function(x) vapply(x, function(t) {
    toks <- strsplit(t, ":", fixed = TRUE)[[1]]
    i <- match(toks, safe)
    toks[!is.na(i)] <- vapply(fnames[i[!is.na(i)]], .dif_term_label, "")
    if ("band" %in% fnames)
      toks[is.na(i) & toks == "band"] <- "(opponent band)"
    paste(toks, collapse = ":")
  }, character(1), USE.NAMES = FALSE)
  planned_dif_terms <- setdiff(attr(stats::terms(stats::as.formula(
    paste("z ~ (", paste(safe, collapse = op), ") * band"))),
    "term.labels"), "band")

  m <- if (is.null(fit$m)) 1L else fit$m
  cats <- if (!is.null(fit$categories)) fit$categories else c("0", "1")
  thr <- if (!is.null(fit$thr_structure)) fit$thr_structure else "free"
  tau <- if (!is.null(fit$thresholds)) fit$thresholds$tau else numeric(1)
  excluded_extreme <- if (is.null(fit$objects$extreme)) character(0) else
    as.character(fit$objects$object[fit$objects$extreme %in% TRUE])
  calibrated_objects <- setdiff(as.character(fit$objects$object),
                                excluded_extreme)
  its <- if (is.null(objects)) calibrated_objects else {
    if (!is.atomic(objects) || !is.null(dim(objects)) || !length(objects) ||
        anyNA(objects))
      stop("`objects` must contain one or more non-missing object names")
    objects <- .role_text_values(objects)
    if (any(!nzchar(objects)))
      stop("`objects` must contain one or more non-missing object names")
    unknown <- setdiff(objects, fit$objects$object)
    if (length(unknown))
      stop("object(s) not in the fit: ", paste(unknown, collapse = ", "))
    if (anyDuplicated(objects))
      stop("object(s) named more than once: ",
           paste(unique(objects[duplicated(objects)]), collapse = ", "),
           "; a repeated object would repeat its tests in the ",
           "multiplicity family")
    boundary <- intersect(objects, excluded_extreme)
    if (length(boundary))
      .refuse("object(s) set aside at a response boundary have no fitted DIF ",
              "model: ", paste(boundary, collapse = ", "))
    objects
  }
  bl <- setNames(fit$objects$location, fit$objects$object)
  jd_all <- if (all(is.na(cm$judge))) NULL else cm$judge

  # base-fit moments per comparison, including any fitted within-judge
  # dependence effects: leaving them out would push the dependence structure
  # into the residuals, where a judge-level factor would absorb it as
  # spurious DIF (and the resolved locations would absorb it as spurious
  # magnitude)
  d0 <- bl[cm$object_a] - bl[cm$object_b]
  Zc <- NULL
  if (!is.null(fit$dependence) &&
      all(fit$dependence$effect %in% names(cm))) {
    Zc <- as.matrix(cm[, fit$dependence$effect, drop = FALSE])
    d0 <- d0 + drop(Zc %*% fit$dependence$estimate)
  }
  sc <- 0:m
  eta <- outer(unname(d0), sc) -
    matrix(rep(c(0, cumsum(tau)), each = nrow(cm)), nrow(cm), m + 1L)
  eta <- eta - apply(eta, 1, max)
  P <- exp(eta); P <- P / rowSums(P)
  E <- drop(P %*% sc); V <- pmax(drop(P %*% sc^2) - E^2, 1e-12)

  # per object: the residual ANOVA z ~ (f1 [+/*] fk) * band, one row per term
  notes <- character(0); term_rows <- list(); caution_count <- 0L
  if (length(excluded_extreme))
    notes <- c(notes, paste0(
      "object(s) set aside at a response boundary were excluded because ",
      "their displayed locations are extrapolations, not fitted DIF ",
      "parameters: ", paste(excluded_extreme, collapse = ", ")))
  if (length(unfitted))
    notes <- c(notes, paste0(
      "judge(s) named in a factor but not in the fitted comparisons, and ",
      "ignored: ", paste(utils::head(sort(unfitted), 5L), collapse = ", "),
      if (length(unfitted) > 5L) ", ..." else ""))
  for (o in its) {
    sel_a <- cm$object_a == o & ok
    sel_b <- cm$object_b == o & ok
    zo <- c((cm$response[sel_a] - E[sel_a]) / sqrt(V[sel_a]),
            -(cm$response[sel_b] - E[sel_b]) / sqrt(V[sel_b]))
    opp <- c(cm$object_b[sel_a], cm$object_a[sel_b])
    wo <- c(cm$weight[sel_a], cm$weight[sel_b])
    gcols <- lapply(gvs, function(g) c(g[sel_a], g[sel_b]))
    oloc <- bl[opp]
    keep <- !is.na(oloc)
    zo <- zo[keep]; oloc <- oloc[keep]; wo <- wo[keep]
    gcols <- lapply(gcols, function(g) g[keep])
    # a count-weighted row stands for `weight` identical comparisons: the
    # residual z then has variance 1/weight, so weighted least squares is
    # exactly the expanded-rows analysis
    n_o <- sum(wo)
    jo <- c(cm$judge[sel_a], cm$judge[sel_b])[keep]
    d <- data.frame(z = zo, w = wo, judge_unit = factor(jo))
    for (j in seq_along(fnames)) d[[safe[j]]] <- factor(gcols[[j]])
    if (n_o < 10 || any(vapply(safe, function(s)
      nlevels(droplevels(d[[s]])) < 2, TRUE))) next
    nb <- if (n_o >= 90) 3L else if (n_o >= 40) 2L else 1L
    if (nb > 1L) {
      # band breaks at count-weighted quantiles of the opponent location:
      # invariant to whether comparisons arrive expanded or count-weighted
      # (row-rank cuts were not), so aggregated and expanded data land in
      # identical bands
      ord <- order(oloc); cw <- cumsum(wo[ord]) / sum(wo)
      br <- unique(vapply(seq_len(nb - 1L) / nb, function(q)
        oloc[ord][which(cw >= q - 1e-12)[1]], 0))
      d$band <- factor(findInterval(oloc, br, left.open = TRUE) + 1L)
      if (nlevels(droplevels(d$band)) < 2L) d$band <- NULL
    }
    if (is.null(d$band)) nb <- 1L
    # JUDGES, not comparisons, are the independent units here: a judge's
    # comparisons share that judge's idiosyncratic preferences, and testing
    # judge-level factors against comparison-level residuals
    # pseudo-replicates -- a null simulation with judge heterogeneity and
    # arbitrary groups falsely flagged uniform DIF in 6 of 10 datasets.
    # Aggregate to one weighted-mean residual per judge (per opponent
    # band), then a split-plot aov with the judge as the error unit: group
    # terms are tested between judges, band terms and their interactions
    # within judges -- the same design logic as the mixed-design person
    # DIF ANOVA.
    cellkey <- if (nb > 1L)
      .factor_cells(data.frame(judge = d$judge_unit, band = d$band), sep = "\r")
               else droplevels(d$judge_unit)
    zbar <- tapply(d$z * d$w, cellkey, sum) / tapply(d$w, cellkey, sum)
    firsts <- which(!duplicated(cellkey))
    ag <- d[firsts[match(levels(cellkey), as.character(cellkey[firsts]))],
            c("judge_unit", safe, if (nb > 1L) "band"), drop = FALSE]
    ag$z <- as.numeric(zbar)
    # between-judge tests need at least two judges per level and enough
    # judges overall to leave residual degrees of freedom
    nj_ok <- all(vapply(safe, function(sn)
      min(tapply(as.character(ag$judge_unit), ag[[sn]],
                 function(x) length(unique(x)))) >= 2L, 0) >= 1)
    if (!nj_ok || length(unique(ag$judge_unit)) < 4L) next
    # the same order-invariant machinery as the person DIF ANOVA:
    # between-judge terms by Type II sums of squares on band-centred judge
    # margins (sequential aov let entry order decide which correlated
    # judge factor flagged), band-crossing terms on the judge-by-band mean
    # matrix through orthonormal contrasts with the Greenhouse-Geisser
    # correction
    rhs_terms <- attr(stats::terms(stats::as.formula(
      paste("z ~ (", paste(safe, collapse = op), ")",
            if (nb > 1L) "* band" else ""))), "term.labels")
    bterms_o <- rhs_terms[!vapply(rhs_terms, function(tt)
      "band" %in% .term_vars(tt), TRUE)]
    wterms_o <- setdiff(rhs_terms, bterms_o)
    if (nb > 1L) {
      bmn <- tapply(ag$z, ag$band, mean)
      zc <- ag$z - as.numeric(bmn[as.character(ag$band)])
    } else zc <- ag$z
    jk <- factor(ag$judge_unit)
    pzj <- tapply(zc, jk, mean)
    jfirst <- which(!duplicated(jk))
    pdat_o <- ag[jfirst[match(levels(jk), as.character(jk[jfirst]))],
                 c("judge_unit", safe), drop = FALSE]
    pdat_o$z <- as.numeric(pzj)
    # Group comparisons are made between judges. HC3 protects the Type-II
    # tests against the unequal precision of judge means when comparison
    # workloads differ.
    ft_b <- .dif_type2(pdat_o, bterms_o, variance = "hc3")
    ft_w <- NULL
    if (nb > 1L && length(wterms_o)) {
      Yw <- tapply(ag$z, list(jk, ag$band), mean)
      compl <- rowSums(is.na(Yw)) == 0L
      if (sum(compl) >= 6L) {
        Yb <- Yw[compl, , drop = FALSE]
        pd2 <- pdat_o[match(rownames(Yb),
                            as.character(pdat_o$judge_unit)), ,
                      drop = FALSE]
        ft_w <- .dif_within_tests(Yb, pd2, "band",
                                  list(band = ncol(Yw)), wterms_o,
                                  bterms_o)
      }
    }
    ft <- rbind(ft_b, ft_w)
    if (is.null(ft)) next
    ft <- ft[ft$term != "Residuals", , drop = FALSE]
    # The global cluster guard does not ensure adequate support in every
    # level of a DIF factor. For each tested term, count the judges and their
    # Kish effective number from this object's comparison workloads. The
    # estimate remains visible when a cell is sparse, but its probability is
    # withheld below the calibrated eight-effective-judge boundary.
    jw <- tapply(d$w, d$judge_unit, sum)
    support <- lapply(ft$term, function(tt) {
      vv <- setdiff(tvars(tt), "band")
      if (!length(vv)) return(c(raw = Inf, effective = Inf))
      cells <- .factor_cells(pdat_o[, vv, drop = FALSE], sep = "\r")
      vals <- lapply(levels(cells), function(cc) {
        ids <- as.character(pdat_o$judge_unit[cells == cc])
        ww <- unname(jw[ids]); ww <- ww[is.finite(ww) & ww > 0]
        c(raw = length(ww), effective = if (length(ww))
          sum(ww)^2 / sum(ww^2) else 0)
      })
      Reduce(pmin, vals)
    })
    ft$min_judges <- vapply(support, `[[`, 0, "raw")
    ft$min_effective_judges <- vapply(support, `[[`, 0, "effective")
    ft$inference_available <- is.infinite(ft$min_effective_judges) |
      (ft$min_judges >= 8L &
         ft$min_effective_judges >= 8 - sqrt(.Machine$double.eps))
    withheld <- !ft$inference_available & ft$term != "band"
    if (any(withheld)) {
      ft$F_value[withheld] <- NA_real_
      ft$p[withheld] <- NA_real_
      notes <- c(notes, sprintf(
        "%s: DIF inference withheld for term(s) below eight judges or eight effective judges in a factor cell: %s",
        o, paste(relab(ft$term[withheld]), collapse = ", ")))
    }
    caution <- ft$inference_available & is.finite(ft$min_effective_judges) &
      ft$min_effective_judges < 9.5
    caution_count <- caution_count + sum(caution)
    for (k in seq_len(nrow(ft)))
      term_rows[[length(term_rows) + 1L]] <- data.frame(
        object = o, term = ft$term[k], df = ft$df[k], sum_sq = ft$sum_sq[k],
        df_denom = ft$df_denom[k], F_value = ft$F_value[k], p = ft$p[k],
        min_judges = ft$min_judges[k],
        min_effective_judges = ft$min_effective_judges[k],
        inference_available = ft$inference_available[k],
        resid_ss = ft$resid_ss[k])
  }
  if (caution_count > 0L)
    notes <- c(notes, sprintf(
      "%d object-term test(s) have 8.0--9.4 effective judges in their smallest cell; see min_effective_judges and interpret these results cautiously",
      caution_count))
  if (!length(term_rows)) stop("no object yielded an estimable DIF ANOVA")
  terms <- do.call(rbind, term_rows); rownames(terms) <- NULL
  # Fix the object-by-term family before item-specific support and opponent
  # banding are inspected. Sparse or aliased questions remain as explicit NA
  # rows, rather than disappearing and reducing the Holm denominator for all
  # other objects.
  wanted <- expand.grid(
    object = its, term = planned_dif_terms,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  wanted_key <- .factor_keys(wanted)
  observed_key <- .factor_keys(terms[c("object", "term")])
  missing_family <- !wanted_key %in% observed_key
  n_unavailable_terms <- sum(missing_family)
  if (n_unavailable_terms) {
    absent <- wanted[missing_family, , drop = FALSE]
    absent$df <- absent$sum_sq <- absent$df_denom <- absent$F_value <-
      absent$p <- absent$min_judges <- absent$min_effective_judges <-
      absent$resid_ss <- NA_real_
    absent$inference_available <- FALSE
    terms <- rbind(terms, absent[names(terms)])
    notes <- c(notes, sprintf(
      "%d requested object-term test(s) were not estimable and remain in the adjusted-probability family",
      n_unavailable_terms))
  }
  terms$eta2_partial <- terms$sum_sq / (terms$sum_sq + terms$resid_ss)
  terms$resid_ss <- NULL
  # Uniform and non-uniform flags feed one reported DIF decision, so the
  # object-by-term tests form one multiplicity family.
  terms$p_adj <- NA_real_
  family_member <- terms$term != "band"
  sel_test <- family_member & is.finite(terms$p)
  terms$p_adj[sel_test] <- p.adjust(
    terms$p[sel_test], method = p_adjust, n = sum(family_member))
  terms$significant <- !is.na(terms$p_adj) & terms$p_adj < alpha
  # a significant higher-order GROUP term supersedes lower-order group terms
  # built from a subset of its factors, within the same object. Band-crossing
  # terms are excluded from the pass: a term's own band interaction is
  # reported WITH it (as non-uniform DIF), so it must not supersede it.
  terms$superseded <- FALSE
  is_group <- !vapply(terms$term, function(t) "band" %in% tvars(t), TRUE)
  for (ob in unique(terms$object)) {
    all_terms <- which(terms$object == ob & is_group)
    higher <- all_terms[terms$significant[all_terms]]
    for (lo in all_terms) for (hi in higher) {
      vl <- tvars(terms$term[lo]); vh <- tvars(terms$term[hi])
      if (length(vl) < length(vh) && all(vl %in% vh))
        terms$superseded[lo] <- TRUE
    }
  }
  # compact reading: one row per object and group term (a term not crossing the
  # opponent band), its own effect uniform DIF and its band crossing
  # non-uniform DIF. Classified on the stand-in tokens, so a factor named
  # "band" (held as f_j) is never confused with the band variable.
  gterms <- unique(terms$term)   # Residuals never reaches term_rows
  gterms <- gterms[!vapply(gterms, function(t) "band" %in% tvars(t), TRUE)]
  srows <- list()
  for (ob in unique(terms$object)) for (tt in gterms) {
    u <- terms[terms$object == ob & terms$term == tt, , drop = FALSE]
    if (!nrow(u)) next
    nu <- terms[terms$object == ob & terms$term == paste0(tt, ":band"), ,
                drop = FALSE]
    srows[[length(srows) + 1L]] <- data.frame(
      object = ob, term = tt,
      F_uniform = u$F_value, p_uniform = u$p, p_uniform_adj = u$p_adj,
      eta2_uniform = u$eta2_partial, uniform_DIF = isTRUE(u$significant),
      min_judges = u$min_judges,
      min_effective_judges = u$min_effective_judges,
      uniform_inference = u$inference_available,
      F_nonuniform = if (nrow(nu)) nu$F_value else NA_real_,
      p_nonuniform = if (nrow(nu)) nu$p else NA_real_,
      p_nonuniform_adj = if (nrow(nu)) nu$p_adj else NA_real_,
      eta2_nonuniform = if (nrow(nu)) nu$eta2_partial else NA_real_,
      nonuniform_DIF = nrow(nu) > 0 && isTRUE(nu$significant),
      nonuniform_inference = if (nrow(nu)) nu$inference_available else NA,
      superseded = isTRUE(u$superseded))
  }
  summary_tab <- if (length(srows)) do.call(rbind, srows) else NULL
  summary_factors <- if (is.null(summary_tab)) list() else
    lapply(summary_tab$term, function(tt) fnames[match(tvars(tt), safe)])
  # Preserve collision-free term identifiers before replacing the internal
  # stand-ins with display labels.  Judge maps, rather than row-wise factor
  # vectors, remain valid when a fitted-null replicate expands a comparison
  # row into several observed response categories.
  term_ids <- terms$term
  summary_term_ids <- if (is.null(summary_tab)) character(0) else
    summary_tab$term
  judges_used <- unique(cm$judge[!is.na(cm$judge)])
  factor_maps <- lapply(seq_along(gvs), function(j) {
    first <- match(judges_used, cm$judge)
    stats::setNames(gvs[[j]][first], judges_used)
  })
  names(factor_maps) <- fnames

  # resolution: for each flagged, non-superseded group term, resolve the object
  # over the COMPLETE judge-factor design. The omnibus terms were fitted
  # jointly, so resolving only the selected term would pool nuisance factors
  # in their observed proportions and change the estimand in an unbalanced
  # design. Main effects are equal-cell marginal differences; interactions
  # are the corresponding tensor contrasts.
  full_factors <- as.data.frame(lapply(gvs, function(x) factor(x)),
                                check.names = FALSE)
  names(full_factors) <- fnames
  full_cell <- .factor_cells(full_factors, sep = ":")
  full_map <- unique(data.frame(cell = as.character(full_cell), full_factors,
                                check.names = FALSE))
  full_map <- full_map[match(levels(full_cell), full_map$cell), , drop = FALSE]
  lev_rows <- list(); sz_rows <- list(); size_family_n <- 0L
  flagged <- if (is.null(summary_tab)) integer(0) else
    which(summary_tab$uniform_DIF & !summary_tab$superseded)
  for (r in flagged) {
    ob <- summary_tab$object[r]; tt <- summary_tab$term[r]; ttd <- relab(tt)
    jf <- match(tvars(tt), safe)
    target <- fnames[jf]
    # Fix every follow-up question opened by this significant term before
    # inspecting the object-specific cell counts. If no cell survives min_n,
    # all questions are unavailable but still belong to the family that also
    # adjusts the other resolved terms.
    planned_n <- prod(vapply(target, function(fn)
      choose(nlevels(full_factors[[fn]]), 2L), numeric(1)))
    new_family_n <- size_family_n + planned_n
    if (!is.finite(new_family_n) || new_family_n > .Machine$integer.max)
      stop("the planned paired-comparison DIF follow-up family is too large to adjust",
           call. = FALSE)
    size_family_n <- as.integer(new_family_n)
    cell <- as.character(full_cell)
    inv <- ok & (cm$object_a == ob | cm$object_b == ob)
    # cell sizes in comparisons (count-weighted), not rows
    lev_n <- tapply(cm$weight[inv], cell[inv], sum)
    use_lev <- names(lev_n)[lev_n >= min_n]
    if (length(use_lev) < 2) {
      notes <- c(notes, sprintf(
        paste0("%s [%s]: fewer than two cells with %d+ comparisons; not ",
               "resolved, and its %d unavailable follow-up question(s) ",
               "remain in the multiplicity family"),
        ob, ttd, min_n, planned_n))
      next
    }
    if (length(use_lev) < length(lev_n))
      notes <- c(notes, sprintf(
        "%s [%s]: cell(s) dropped with fewer than %d comparisons: %s",
        ob, ttd, min_n, paste(setdiff(names(lev_n), use_lev), collapse = ", ")))
    # Declare this term's contrast family before attempting its resolution.
    # An anchor, name collision or failed refit makes the corresponding
    # probabilities unavailable; it does not remove those planned questions
    # from the pooled multiplicity family promised by the public interface.
    map_use <- full_map[match(use_lev, full_map$cell), , drop = FALSE]
    fam <- tryCatch(.dif_posthoc_family(full_factors, map_use, target,
                                        within = character(0)),
                    error = function(e) e)
    if (inherits(fam, "error")) {
      notes <- c(notes, sprintf("%s [%s]: %s", ob, ttd,
                                conditionMessage(fam)))
      next
    }
    if (!identical(as.integer(fam$planned_n), as.integer(planned_n)))
      stop("internal paired-comparison DIF family size changed during resolution",
           call. = FALSE)
    if (!length(fam$family)) {
      notes <- c(notes, sprintf(
        "%s [%s]: no complete nuisance-factor cells support a post-hoc contrast; the unavailable questions remain in the multiplicity family",
        ob, ttd))
      next
    }
    unavailable <- setdiff(fam$planned, names(fam$family))
    if (length(unavailable))
      notes <- c(notes, sprintf(
        "%s [%s]: post-hoc contrast(s) without a complete nuisance-factor stratum were not estimated but remain in the multiplicity family: %s",
        ob, ttd, paste(unavailable, collapse = ", ")))
    if (!is.null(fit$anchors) && ob %in% names(fit$anchors)) {
      notes <- c(notes, sprintf(
        "%s [%s]: the object is externally anchored and cannot be resolved; fixing every copy at the anchor would define its DIF as zero",
        ob, ttd))
      next
    }
    rsel <- ok & (!(cm$object_a == ob | cm$object_b == ob) | cell %in% use_lev)
    a2 <- cm$object_a[rsel]; b2 <- cm$object_b[rsel]; c2 <- cell[rsel]
    made <- paste0(ob, " (", unique(c2[a2 == ob | b2 == ob]), ")")
    clash <- intersect(made, as.character(fit$objects$object))
    if (length(clash)) {
      notes <- c(notes, sprintf(
        "%s [%s]: the resolved name(s) %s already name object(s) in the fit; magnitude withheld rather than merging them",
        ob, ttd, paste(clash, collapse = ", ")))
      next
    }
    a2 <- ifelse(a2 == ob, paste0(ob, " (", c2, ")"), a2)
    b2 <- ifelse(b2 == ob, paste0(ob, " (", c2, ")"), b2)
    # the refit keeps the fitted dependence structure: the history covariates
    # are fixed by the original judgment sequence, so they pass through as-is
    rf <- tryCatch(.btl_graded(
      a2, b2, cm$response[rsel], if (is.null(jd_all)) NULL else jd_all[rsel],
      cm$weight[rsel], cats, maxit, tol, character(0), thr = thr,
      Z = if (is.null(Zc)) NULL else Zc[rsel, , drop = FALSE],
      anchors = fit$anchors, withheld_notes = FALSE),
      error = function(e) NULL)
    if (is.null(rf) || !isTRUE(rf$converged)) {
      notes <- c(notes, sprintf(
        "%s [%s]: the resolved calibration did not converge; magnitude withheld",
        ob, ttd))
      next
    }
    if (length(rf$notes))
      notes <- c(notes, sprintf("%s [%s]: %s", ob, ttd, rf$notes))
    copies <- paste0(ob, " (", use_lev, ")")
    idx <- match(copies, rf$objects$object)
    if (anyNA(idx) || any(rf$objects$extreme[idx] %in% TRUE)) {
      notes <- c(notes, sprintf("%s [%s]: resolved copies missing", ob, ttd))
      next
    }
    loc <- rf$objects$location[idx]
    vv <- tryCatch(rf$cov_beta[copies, copies, drop = FALSE],
                   error = function(e)
                     matrix(NA_real_, length(copies), length(copies)))
    vv_ok <- identical(dim(vv), c(length(copies), length(copies))) &&
      all(is.finite(vv)) && .covariance_is_symmetric(vv) &&
      .covariance_is_psd(vv)
    # A resolved cell is a judge-level group estimate. The base fit's global
    # cluster guard is not enough when a many-level factor leaves only a few
    # judges in each cell: simulations with four to six judges per level gave
    # materially anti-conservative pairwise tests despite 12--20 judges
    # overall. Apply the same calibrated eight-effective-cluster boundary to
    # EACH resolved level, using only that object's comparisons in the level.
    # The locations remain useful descriptively when inference is withheld.
    lev_j <- lev_eff <- setNames(numeric(length(use_lev)), use_lev)
    for (lv in use_lev) {
      rr <- inv & cell == lv
      jw <- tapply(cm$weight[rr], cm$judge[rr], sum)
      jw <- jw[is.finite(jw) & jw > 0]
      lev_j[lv] <- length(jw)
      sh <- jw / sum(jw)
      lev_eff[lv] <- if (length(sh)) 1 / sum(sh^2) else 0
    }
    lev_ok <- lev_j >= 8 & lev_eff >= 8 - sqrt(.Machine$double.eps)
    if (!isTRUE(rf$cl$inference_available)) {
      lev_ok[] <- FALSE
      notes <- c(notes, sprintf(
        "%s [%s]: the resolved calibration does not support cluster-robust inference; locations are descriptive",
        ob, ttd))
    }
    if (!vv_ok) {
      lev_ok[] <- FALSE
      notes <- c(notes, sprintf(
        "%s [%s]: the resolved-location covariance is unavailable or not positive semidefinite; locations and differences are descriptive",
        ob, ttd))
    }
    if (any(!lev_ok))
      notes <- c(notes, sprintf(
        "%s [%s]: pairwise inference is withheld for level(s) below eight effective judges: %s; resolved locations and differences remain descriptive",
        ob, ttd, paste(sprintf("%s (%.1f)", use_lev[!lev_ok],
                               lev_eff[!lev_ok]), collapse = ", ")))
    caution <- lev_ok & lev_eff < 9.5
    if (any(caution))
      notes <- c(notes, sprintf(
        "%s [%s]: level(s) %s have 8.0--9.4 effective judges; interpret pairwise inference cautiously",
        ob, ttd, paste(use_lev[caution], collapse = ", ")))
    lev_se <- if (vv_ok) sqrt(pmax(diag(vv), 0)) else
      rep(NA_real_, length(use_lev))
    lev_se[!lev_ok] <- NA_real_
    lev_rows[[length(lev_rows) + 1L]] <- data.frame(
      object = ob, term = tt, level = use_lev, location = loc,
      se = lev_se, n = as.numeric(lev_n[use_lev]),
      n_judges = unname(lev_j), effective_judges = unname(lev_eff))
    for (nm in names(fam$family)) {
      w <- fam$family[[nm]][use_lev]
      w[is.na(w)] <- 0
      pos <- which(w > 0); neg <- which(w < 0); used <- c(pos, neg)
      if (!length(pos) || !length(neg)) next
      contrast_ok <- all(lev_ok[used])
      contrast_se <- if (contrast_ok)
        sqrt(max(drop(t(w) %*% vv %*% w), 0)) else NA_real_
      # The established two-cell comparison uses its validated Welch
      # reference. There is no corresponding Welch formula based only on
      # cell effective counts for a contrast spanning three or more fitted
      # cells, so use the weakest contributing cell as a conservative
      # small-cluster reference instead of treating pooled cell counts as
      # independent degrees of freedom.
      contrast_df <- if (!contrast_ok) {
        NA_real_
      } else {
        .btl_dif_contrast_df(lev_eff[used])
      }
      sz_rows[[length(sz_rows) + 1L]] <- data.frame(
        object = ob, term = tt, contrast = nm,
        level_a = paste(use_lev[pos], collapse = ", "),
        level_b = paste(use_lev[neg], collapse = ", "),
        difference = sum(w * loc),
        n_judges_a = sum(lev_j[pos]), n_judges_b = sum(lev_j[neg]),
        effective_judges_a = sum(lev_eff[pos]),
        effective_judges_b = sum(lev_eff[neg]),
        se = contrast_se, df = contrast_df,
        stringsAsFactors = FALSE)
    }
  }
  # a summary row can be flagged yet carry no magnitude row; say so rather
  # than leave the omission silent
  if (!is.null(summary_tab)) {
    for (r in which(summary_tab$uniform_DIF & summary_tab$superseded))
      notes <- c(notes, sprintf(
        "%s [%s]: uniform DIF superseded by a higher-order term; not resolved separately",
        summary_tab$object[r], relab(summary_tab$term[r])))
    for (r in which(summary_tab$nonuniform_DIF & !summary_tab$uniform_DIF))
      notes <- c(notes, sprintf(
        "%s [%s]: non-uniform DIF only; no single location difference summarises it, so no magnitude row is reported",
        summary_tab$object[r], relab(summary_tab$term[r])))
  }
  levels_df <- if (length(lev_rows)) do.call(rbind, lev_rows) else NULL
  sizes <- if (length(sz_rows)) do.call(rbind, sz_rows) else NULL
  if (!is.null(sizes)) {
    sizes$t <- .wald_ratio(sizes$difference, sizes$se)
    sizes$p <- 2 * stats::pt(-abs(sizes$t), df = sizes$df)
    sizes$p_adj <- .p_adjust_family(sizes$p, method = p_adjust,
                                    n = size_family_n)
    sizes$significant <- ifelse(is.finite(sizes$p_adj),
                                sizes$p_adj < alpha, NA)
    sizes$practical <- abs(sizes$difference) >= flag_logits
    rownames(sizes) <- NULL
  }
  if (!is.null(levels_df)) rownames(levels_df) <- NULL
  # relabel stand-ins to the factor names for display, now that all term
  # classification and resolution are done
  terms$term <- relab(terms$term)
  if (!is.null(summary_tab)) {
    summary_tab$term <- relab(summary_tab$term); rownames(summary_tab) <- NULL
  }
  if (!is.null(sizes)) sizes$term <- relab(sizes$term)
  if (!is.null(levels_df)) levels_df$term <- relab(levels_df$term)
  out <- list(summary = summary_tab, terms = terms, levels = levels_df,
              sizes = sizes, summary_factors = summary_factors,
              term_ids = term_ids,
              summary_term_ids = summary_term_ids,
              size_family_n = size_family_n,
              effects = effects, factors = fnames,
              alpha = alpha, p_adjust = p_adjust, flag_logits = flag_logits,
              notes = unique(notes),
              fit_signature = .fit_boot_signature(fit),
              bootstrap_design = list(
                factors = factor_maps, objects = its, effects = effects,
                p_adjust = p_adjust, alpha = alpha,
                flag_logits = flag_logits, min_n = min_n,
                maxit = maxit, tol = tol))
  out <- .tag_tables(out)
  out$result_signature <- .fit_boot_md5(out)
  class(out) <- "rasch_btl_dif"
  out
}

#' @export
print.rasch_btl_dif <- function(x, ...) {
  nf <- length(x$factors)
  cat(sprintf("DIF for paired comparisons: %d factor(s) [%s], %s effects\n",
              nf, paste(x$factors, collapse = ", "), x$effects))
  cat("Residual ANOVA per object and term (uniform = term; non-uniform = term x opponent band)\n")
  print(.fmt_df(x$summary[, c("object", "term", "F_uniform", "p_uniform_adj",
                              "uniform_DIF", "F_nonuniform",
                              "p_nonuniform_adj", "nonuniform_DIF")]),
        row.names = FALSE)
  if (!is.null(x$sizes)) {
    cat(sprintf("\nResolved locations (logits; %s over %d planned comparison(s); practical %.2f)\n",
                x$p_adjust, x$size_family_n %||% nrow(x$sizes),
                x$flag_logits))
    cols <- c("object", "term", "level_a", "level_b", "difference",
              "n_judges_a", "n_judges_b", "effective_judges_a",
              "effective_judges_b", "se", "t", "df", "p_adj",
              "significant", "practical")
    print(.fmt_df(x$sizes[, intersect(cols, names(x$sizes)), drop = FALSE]),
          row.names = FALSE)
  }
  if (length(x$notes)) cat("Notes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}
