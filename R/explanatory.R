# rasch :: explanatory item and threshold models
# =============================================================================
# The LLTM and LPCM retain the Rasch response function and replace freely
# estimated thresholds by linear functions of observed item or threshold
# characteristics. Estimation uses the same pairwise conditional likelihood as
# pcml(); only its threshold design matrix changes. Fixed departures are added
# as nominated design columns, so relaxation remains a conditional Rasch model.
# =============================================================================

.explanatory_projector <- function(m, thr) {
  L <- length(m); M <- nrow(thr)
  a <- 1 / (L * m[thr$item])
  diag(M) - matrix(1, M, 1L) %*% t(a)
}

# Apply the origin constraint after constructing the model matrix. Subtract a
# reference entry first, so an arbitrary large column offset is removed before
# averaging. Do not centre the predictors before model.matrix(): that could
# change a formula containing interactions without their main effects.
.explanatory_centre <- function(mm, weights = rep(1 / nrow(mm), nrow(mm))) {
  relative <- sweep(mm, 2L, mm[1L, ], `-`)
  out <- sweep(relative, 2L, drop(crossprod(weights, relative)), `-`)
  if (any(!is.finite(out)))
    stop("the centred explanatory design is outside the representable range; ",
         "rescale the predictors", call. = FALSE)
  out
}

# model.matrix() can give different design columns the same label: for
# example factor x's level B and a numeric predictor called xB. Retain the
# design, but give every coefficient an unambiguous stored name. New fixed
# departures also avoid names already present in the active model.
.explanatory_unique_columns <- function(B, previous = character()) {
  if (ncol(B))
    colnames(B) <- utils::tail(make.unique(c(previous, colnames(B))), ncol(B))
  B
}

.explanatory_ordinal_contrasts <- function(x) {
  k <- nlevels(x)
  if (k < 2L) stop("an ordinal predictor needs at least two observed levels")
  out <- outer(seq_len(k), seq_len(k - 1L), `>`) * 1
  # make.names() is not one-to-one: levels such as "a b", "a-b" and "a.b"
  # all collapse to the same token. The adjacent-contrast number is part of
  # the parameter name so distinct ordered transitions remain distinct even
  # when their readable labels have the same syntactic form.
  colnames(out) <- paste0(
    "adjacent_", seq_len(k - 1L), "_",
    make.names(levels(x)[-1L]), "_vs_", make.names(levels(x)[-k]))
  out
}

# Item names are response-column selectors, not labels to canonicalise.
.explanatory_match_items <- function(x, items) {
  missing <- is.na(x)
  x <- as.character(x)
  x[missing | is.na(x)] <- NA_character_
  vapply(x, function(nm) {
    if (is.na(nm) || nm %in% items) return(nm)
    hit <- items[trimws(items) == trimws(nm)]
    if (length(hit) > 1L)
      stop("ambiguous item name: ", nm, "; use the exact response-column name",
           call. = FALSE)
    if (length(hit) == 1L) hit else nm
  }, character(1), USE.NAMES = FALSE)
}

.explanatory_metadata <- function(predictors, formula, X,
                                  level = c("item", "threshold")) {
  level <- match.arg(level)
  if (!is.data.frame(predictors))
    stop("`predictors` must be a data frame")
  .check_column_names(predictors)
  if (!"item" %in% names(predictors))
    stop("`predictors` needs an `item` column")
  if (!inherits(formula, "formula") || length(formula) != 2L)
    stop("`formula` must be one-sided, for example ~ format + threshold")

  items <- colnames(X); m <- apply(X, 2L, max, na.rm = TRUE)
  thr <- threshold_index(m)
  index <- data.frame(item = items[thr$item],
                      threshold_number = thr$k,
                      stringsAsFactors = FALSE)
  predictors$item <- .explanatory_match_items(predictors$item, items)
  if (anyNA(predictors$item) || any(!nzchar(trimws(predictors$item))))
    stop("predictor item names must be non-missing and non-empty")
  unknown <- setdiff(unique(predictors$item), items)
  if (length(unknown))
    stop("predictor item(s) are not present in the fitted response data: ",
         paste(unknown, collapse = ", "))

  if (level == "item") {
    if (anyDuplicated(predictors$item))
      stop("item-level predictors need exactly one row per item")
    missing <- setdiff(items, predictors$item)
    if (length(missing))
      stop("item-level predictors are missing: ", paste(missing, collapse = ", "))
    meta <- predictors[match(index$item, predictors$item), , drop = FALSE]
    meta$threshold_number <- index$threshold_number
  } else {
    if (!"threshold" %in% names(predictors) &&
        !"threshold_number" %in% names(predictors))
      stop("threshold-level predictors need a `threshold` column")
    kn <- if ("threshold" %in% names(predictors))
      predictors$threshold else predictors$threshold_number
    kn_num <- suppressWarnings(as.numeric(as.character(kn)))
    if (anyNA(kn_num) || any(!is.finite(kn_num)) ||
        any(kn_num != floor(kn_num)) || any(kn_num < 1))
      stop("threshold numbers must be positive integers")
    kn <- as.integer(kn_num)
    key <- paste(predictors$item, kn, sep = "\r")
    if (anyDuplicated(key))
      stop("threshold-level predictors need one row per item and threshold")
    wanted <- paste(index$item, index$threshold_number, sep = "\r")
    missing <- which(!wanted %in% key)
    extra <- setdiff(key, wanted)
    if (length(missing))
      stop("threshold-level predictors are missing: ",
           paste(paste0(index$item[missing], " threshold ",
                        index$threshold_number[missing]), collapse = ", "))
    if (length(extra))
      stop("predictor rows do not correspond to an observed item threshold: ",
           paste(gsub("\r", " threshold ", extra, fixed = TRUE),
                 collapse = ", "))
    meta <- predictors[match(wanted, key), , drop = FALSE]
    meta$threshold_number <- index$threshold_number
  }
  meta$item <- index$item
  meta$threshold <- factor(index$threshold_number,
                           levels = sort(unique(index$threshold_number)))
  rownames(meta) <- NULL

  reserved <- c("item", "threshold", "threshold_number")
  for (nm in setdiff(names(meta), reserved)) {
    if (is.character(meta[[nm]]) || is.logical(meta[[nm]]))
      meta[[nm]] <- factor(meta[[nm]])
    # a level no item carries -- after a subset, a drop, or simply a level
    # the predictor table declares and never uses -- contributes an
    # all-zero column and makes an identified design look rank-deficient
    else if (is.factor(meta[[nm]])) meta[[nm]] <- droplevels(meta[[nm]])
  }
  mf <- tryCatch(stats::model.frame(formula, data = meta,
                                    na.action = stats::na.fail),
                 error = function(e) stop("cannot construct the explanatory ",
                   "model: ", conditionMessage(e), call. = FALSE))
  # model.matrix() omits offsets; no fixed contribution is carried through
  # explanatory estimation or its refits, so accepting one would change the
  # requested model silently. Inspect terms rather than predictor names.
  if (length(attr(attr(mf, "terms"), "offset")))
    stop("formula offsets are not supported by explanatory models; remove offset() terms",
         call. = FALSE)
  # Unused metadata must not constrain the requested model.
  for (nm in setdiff(names(mf), reserved))
    if (is.ordered(mf[[nm]]))
      contrasts(mf[[nm]]) <- .explanatory_ordinal_contrasts(mf[[nm]])
  mm <- tryCatch(stats::model.matrix(formula, data = mf),
                 error = function(e) stop("cannot construct the explanatory ",
                   "model matrix: ", conditionMessage(e), call. = FALSE))
  if (!ncol(mm)) stop("the explanatory formula produced no predictors")
  if (any(!is.finite(mm)))
    stop("the explanatory predictors produce non-finite model-matrix values")
  mm <- .explanatory_unique_columns(mm)

  A <- .explanatory_projector(m, thr)
  B0 <- .explanatory_centre(mm, 1 / (length(m) * m[thr$item]))
  zero <- colSums(B0 != 0) == 0L
  if (any(zero & colnames(mm) != "(Intercept)"))
    stop("predictor(s) have no estimable variation after fixing the scale ",
         "origin: ", paste(colnames(mm)[zero &
           colnames(mm) != "(Intercept)"], collapse = ", "))
  B <- B0[, !zero, drop = FALSE]
  mm_keep <- mm[, !zero, drop = FALSE]
  if (!ncol(B)) stop("the explanatory formula contains only an intercept")
  qb <- qr(sweep(B, 2L, .design_column_scale(B), `/`), tol = 1e-10)
  if (qb$rank < ncol(B)) {
    aliased <- colnames(B)[qb$pivot[seq.int(qb$rank + 1L, ncol(B))]]
    stop("the explanatory design is not identified; aliased term(s): ",
         paste(aliased, collapse = ", "))
  }
  if (ncol(B) > nrow(B) - 1L)
    stop("the explanatory design has more parameters than the free calibration")
  colnames(B) <- colnames(mm_keep)
  list(B = B, matrix = mm_keep, metadata = meta,
       source_predictors = predictors, threshold_index = thr,
       m = m, level = level, formula = formula, projector = A)
}

.pcml_design <- function(X, B, parameter_names = colnames(B), maxit = 60,
                         tol = 1e-8, cluster = NULL) {
  X <- as.matrix(X); .check_integer_scores(X, "the score matrix")
  storage.mode(X) <- "integer"
  m <- apply(X, 2L, max, na.rm = TRUE); L <- ncol(X)
  thr <- threshold_index(m); M <- nrow(thr)
  inames <- colnames(X) %||% paste0("V", seq_len(L))
  if (!is.matrix(B) || nrow(B) != M)
    stop("the explanatory design must have one row per fitted threshold")
  if (!ncol(B))
    stop("the explanatory design matrix is not full column rank")
  bs <- .design_column_scale(B)
  B_work <- sweep(B, 2L, bs, `/`)
  if (qr(B_work, tol = 1e-10)$rank < ncol(B))
    stop("the explanatory design matrix is not full column rank")
  pairs <- .pair_counts(X, m)
  .pcml_check_connected(pairs, L, inames)
  weak <- .pcml_weak_thresholds(X, m, thr, inames)
  st <- .start_tau(X, thr)
  beta0 <- tryCatch(qr.solve(B_work, st, tol = 1e-10) / bs,
                    error = function(e) rep(0, ncol(B)))
  beta0[!is.finite(beta0)] <- 0
  sol <- .pcml_solve(X, thr, m, B, beta0, maxit = maxit, tol = tol,
                     pairs = pairs, cluster = cluster)
  repeated <- !is.null(sol$cluster_support) &&
    isTRUE(sol$cluster_support$repeated)
  if (isTRUE(sol$converged) && isTRUE(sol$cluster_inference) && repeated) {
    cov_small <- .pcml_linearised_cluster_cov(
      X, thr, m, sol$tau, pairs, B, sol$H_beta, cluster)
    if (is.null(cov_small)) {
      sol$cluster_inference <- FALSE
      sol$cov_beta[,] <- NA_real_
      sol$cov_tau[,] <- NA_real_
      sol$se_tau[] <- NA_real_
      sol$cluster_note <- paste(
        "item-parameter uncertainty withheld: deleting at least one person",
        "cluster leaves the explanatory calibration unidentified")
      sol$cluster_support$correction <- "withheld"
    } else {
      sol$cov_beta <- cov_small
      sol$cov_tau <- B %*% cov_small %*% t(B)
      sol$se_tau <- sqrt(pmax(diag(sol$cov_tau), 0))
      sol$cluster_support$correction <-
        "linearised delete-one-person jackknife"
    }
  } else if (!is.null(sol$cluster_support)) {
    sol$cluster_support$correction <- if (repeated) "CR1" else "none"
  }
  names(sol$beta) <- parameter_names
  dimnames(sol$cov_beta) <- list(parameter_names, parameter_names)
  # Do not turn the curvature at an unfinished optimisation iterate into
  # apparently valid coefficient or threshold inference.  Point estimates
  # remain available to diagnose the failed fit.
  if (!isTRUE(sol$converged)) {
    sol$se_tau[] <- NA_real_
    sol$cov_tau[,] <- NA_real_
    sol$cov_beta[,] <- NA_real_
  }
  thr$tau <- sol$tau
  thr$se <- sol$se_tau
  thr$anchored <- FALSE
  thr$weak <- weak$flag
  thr$se[thr$weak] <- NA_real_
  se <- sqrt(pmax(diag(sol$cov_beta), 0))
  stat <- .wald_ratio(sol$beta, se)
  ref_df <- if (!isTRUE(sol$cluster_inference)) NA_real_ else if (repeated)
    sol$cluster_support$n - 1L else Inf
  coef <- data.frame(term = parameter_names, estimate = sol$beta, se = se,
                     t = stat, df = ref_df,
                     p = 2 * stats::pt(-abs(stat), df = ref_df),
                     stringsAsFactors = FALSE)
  coef$p_adj <- .p_adjust_family(coef$p, method = "holm")
  rownames(coef) <- coef$term
  list(model = "explanatory", thr = thr, cov_tau = sol$cov_tau,
       loglik = sol$loglik, iterations = sol$iterations,
       converged = sol$converged, m = m, anchors = NULL,
       n_parameters = ncol(B), B = B, beta = sol$beta,
       coefficients = coef, cov_beta = sol$cov_beta,
       H_beta = sol$H_beta,
       notes = c(weak$notes, sol$cluster_note),
       cluster_inference = sol$cluster_inference,
       cluster_support = sol$cluster_support)
}

.pcml_nested_test <- function(full, restricted) {
  if (!isTRUE(full$converged) || !isTRUE(restricted$converged))
    stop("both conditional calibrations must converge before comparison")
  sf <- .design_column_scale(full$B)
  Bf <- sweep(full$B, 2L, sf, `/`)
  Br <- sweep(restricted$B, 2L, .design_column_scale(restricted$B), `/`)
  if (nrow(Bf) != nrow(Br))
    stop("the compared models do not describe the same thresholds")
  M <- nrow(Bf)
  S <- cbind(Br, rep(1, M))
  ss <- svd(S)
  rs <- sum(ss$d > max(1e-10, max(ss$d) * 1e-8))
  U <- ss$u[, seq_len(rs), drop = FALSE]
  Portho <- diag(M) - tcrossprod(U)
  A <- Portho %*% Bf
  sa <- svd(t(A))
  r <- sum(sa$d > max(1e-8, max(sa$d) * 1e-8))
  W <- max(0, 2 * (full$loglik - restricted$loglik))
  if (!r) return(list(chisq = W, df = 0L, p = NA_real_,
                      chisq_kent = NA_real_, p_kent = NA_real_,
                      lambda = numeric(0)))
  C <- sa$u[, seq_len(r), drop = FALSE]
  Hinv <- tryCatch(solve(-full$H_beta / outer(sf, sf)),
                   error = function(e) NULL)
  kc <- if (is.null(Hinv))
    list(chisq = NA_real_, p = NA_real_, lambda = numeric(0)) else
    .kent_calibration(W, C, full$cov_beta * outer(sf, sf), Hinv)
  list(chisq = W, df = r, p = stats::pchisq(W, r, lower.tail = FALSE),
       chisq_kent = kc$chisq, p_kent = kc$p, lambda = kc$lambda)
}

.btl_explanatory_design <- function(predictors, formula, objects,
                                    known = objects) {
  if (!is.data.frame(predictors) || !"object" %in% names(predictors))
    stop("`predictors` must be a data frame with an `object` column")
  .check_column_names(predictors)
  if (!inherits(formula, "formula") || length(formula) != 2L)
    stop("`formula` must be one-sided, for example ~ domain + format")
  predictors$object <- .role_text_values(predictors$object)
  if (anyNA(predictors$object) || any(!nzchar(predictors$object)))
    stop("predictor object names must be non-missing and non-empty")
  if (anyDuplicated(predictors$object))
    stop("object predictors need exactly one row per object")
  missing <- setdiff(objects, predictors$object)
  if (length(missing))
    stop("object predictors are missing: ", paste(missing, collapse = ", "))
  # a predictor row for an object the comparisons never mention is an error;
  # a row for an object that WAS compared but was set aside at a response
  # boundary is not, and is simply not part of the design
  extra <- setdiff(predictors$object, known)
  if (length(extra))
    stop("predictor row(s) for object(s) not present in the comparisons: ",
         paste(extra, collapse = ", "))
  meta <- predictors[match(objects, predictors$object), , drop = FALSE]
  for (nm in setdiff(names(meta), "object")) {
    if (is.character(meta[[nm]]) || is.logical(meta[[nm]]))
      meta[[nm]] <- factor(meta[[nm]])
    # a level no calibrated object carries would add an all-zero column
    else if (is.factor(meta[[nm]])) meta[[nm]] <- droplevels(meta[[nm]])
  }
  mf <- tryCatch(stats::model.frame(formula, data = meta,
                                    na.action = stats::na.fail),
    error = function(e) stop("cannot construct the explanatory object model: ",
                             conditionMessage(e), call. = FALSE))
  if (length(attr(attr(mf, "terms"), "offset")))
    stop("formula offsets are not supported by explanatory models; remove offset() terms",
         call. = FALSE)
  # Unused metadata must not constrain the requested model.
  for (nm in setdiff(names(mf), "object"))
    if (is.ordered(mf[[nm]]))
      contrasts(mf[[nm]]) <- .explanatory_ordinal_contrasts(mf[[nm]])
  mm <- tryCatch(stats::model.matrix(formula, data = mf),
    error = function(e) stop("cannot construct the explanatory object model ",
                             "matrix: ", conditionMessage(e), call. = FALSE))
  if (any(!is.finite(mm)))
    stop("the explanatory object predictors produce non-finite model-matrix values")
  mm <- .explanatory_unique_columns(mm)
  B0 <- .explanatory_centre(mm)
  keep <- colSums(B0 != 0) > 0L
  bad <- colnames(mm)[!keep & colnames(mm) != "(Intercept)"]
  if (length(bad))
    stop("predictor(s) have no estimable variation after fixing the scale origin: ",
         paste(bad, collapse = ", "))
  B <- B0[, keep, drop = FALSE]
  if (!ncol(B)) stop("the explanatory formula contains only an intercept")
  if (qr(sweep(B, 2L, .design_column_scale(B), `/`),
         tol = 1e-10)$rank < ncol(B))
    stop("the explanatory object design is not identified")
  rownames(B) <- objects
  list(B = B, offset = stats::setNames(numeric(length(objects)), objects),
       parameter_names = colnames(B), metadata = meta,
       source_predictors = predictors, formula = formula)
}

.btl_explanatory_nested_test <- function(full, restricted) {
  if (!isTRUE(full$converged) || !isTRUE(restricted$converged))
    stop("both comparative judgement models must converge before comparison")
  if (!identical(as.character(full$objects$object),
                 as.character(restricted$objects$object)))
    stop("the compared models do not contain the same objects")
  sf <- .design_column_scale(full$location_design)
  Bf <- sweep(full$location_design, 2L, sf, `/`)
  Br <- sweep(restricted$location_design, 2L,
              .design_column_scale(restricted$location_design), `/`)
  K <- nrow(Bf)
  S <- cbind(Br, rep(1, K))
  ss <- svd(S); rs <- sum(ss$d > max(1e-10, max(ss$d) * 1e-8))
  P <- diag(K) - tcrossprod(ss$u[, seq_len(rs), drop = FALSE])
  A <- P %*% Bf
  sa <- svd(t(A)); r <- sum(sa$d > max(1e-8, max(sa$d) * 1e-8))
  W <- max(0, 2 * (full$loglik - restricted$loglik))
  if (!r) return(list(chisq = W, df = 0L, p = NA_real_,
                      chisq_kent = NA_real_, p_kent = NA_real_,
                      lambda = numeric(0)))
  Cobj <- sa$u[, seq_len(r), drop = FALSE]
  C <- rbind(Cobj,
             matrix(0, nrow(full$sensitivity) - nrow(Cobj), r))
  available <- isTRUE(full$cl$inference_available)
  sp <- c(sf, rep(1, nrow(full$sensitivity) - length(sf)))
  scale_outer <- outer(sp, sp)
  Hinv <- if (available)
    tryCatch(solve(full$sensitivity / scale_outer),
             error = function(e) NULL) else NULL
  kc <- if (is.null(Hinv))
    list(chisq = NA_real_, p = NA_real_, lambda = numeric(0)) else
    .kent_calibration(W, C, full$cov_parameters * scale_outer, Hinv)
  list(chisq = W, df = r,
       p = stats::pchisq(W, r, lower.tail = FALSE),
       chisq_kent = kc$chisq, p_kent = kc$p, lambda = kc$lambda)
}

.btl_explanatory_refit <- function(fit, B, relaxations) {
  B <- .explanatory_unique_columns(B)
  spec <- fit$explanatory$refit_spec
  design <- list(B = B,
                 offset = stats::setNames(numeric(nrow(B)), rownames(B)),
                 parameter_names = colnames(B))
  out <- do.call(btl, c(list(data = spec$data), spec$args,
                        list(.object_design = design)))
  out$reference_fit <- fit$reference_fit
  out$explanatory <- fit$explanatory
  out$explanatory$active_B <- B
  out$explanatory$relaxations <- relaxations
  class(out) <- c("rasch_btl_explanatory", "rasch_btl")
  out
}

#' Fit an explanatory comparative judgement model
#'
#' Constrains Bradley--Terry--Luce object locations to linear functions of
#' observed object characteristics. The formulation applies to dichotomous
#' and ordered comparative judgements; ordered-response thresholds retain the
#' structure selected in \code{thresholds}.
#'
#' @details For objects \eqn{a} and \eqn{b},
#' \deqn{\log\{P(a \succ b)/P(b \succ a)\}=\beta_a-\beta_b,\qquad
#' \beta_i=\mathbf z_i^{T}\boldsymbol\gamma.}
#' The scale origin is fixed at mean object location zero. Numeric predictors
#' are continuous, unordered factors are categorical, and ordered factors use
#' successive contrasts between adjacent levels. Character predictors are
#' converted to unordered factors. Selected
#' interactions may be included in \code{formula}.
#' Design columns are centred and rescaled internally for numerical stability; reported
#' coefficients and standard errors use the supplied predictor units.
#' Coincident coefficient labels receive numeric suffixes; this does not
#' change the predictor design.
#' A free calibration is retained for \code{explanatory_test()}. Standard errors use the same
#' sandwich covariance as \code{btl()}; when judges are identified, coefficient
#' tests use the judge-clustered covariance and a \eqn{t} reference with
#' judge-cluster degrees of freedom. Holm adjustment covers the coefficient
#' family.
#'
#' @inheritParams btl
#' @param predictors Data frame with one row per object, an \code{object}
#'   column, and the predictors named in \code{formula}. Column names must be
#'   unique.
#' @param formula One-sided explanatory formula, including selected
#'   interactions if required. Formula offsets (\code{offset()}) are not supported.
#' @return An object of class \code{"rasch_btl_explanatory"}, inheriting from
#'   \code{"rasch_btl"}.
#' @references Bradley, R. A. and Terry, M. E. (1952). Rank analysis of
#' incomplete block designs: I. The method of paired comparisons. Biometrika,
#' 39, 324--345.
#'
#' Fischer, G. H. (1973). The linear logistic test model as an instrument in
#' educational research. Acta Psychologica, 37, 359--374.
#' @examples
#' set.seed(1)
#' q <- data.frame(object = LETTERS[1:6],
#'                 domain = rep(0:1, each = 3))
#' beta <- setNames(0.8 * q$domain, q$object)
#' pr <- t(combn(q$object, 2))
#' d <- data.frame(a = rep(pr[, 1], each = 20),
#'                 b = rep(pr[, 2], each = 20))
#' p <- plogis(beta[d$a] - beta[d$b])
#' d$winner <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' fit <- btl_explanatory(d, q, ~ domain, "a", "b", winner = "winner")
#' fit$object_coefficients
#' explanatory_test(fit)
#' @seealso \code{\link{btl}}, \code{\link{explanatory_test}},
#'   \code{\link{explanatory_diagnostics}}, and
#'   \code{\link{relax_btl_explanatory}}.
#' @export
btl_explanatory <- function(data, predictors, formula, object_a, object_b,
                            winner = NULL, response = NULL, margin = NULL,
                            judge = NULL, count = NULL, order = NULL,
                            position = FALSE,
                            ties = c("drop", "half", "error"),
                            thresholds = c("free", "pc"), maxit = 60,
                            tol = 1e-8) {
  ties <- match.arg(ties); thresholds <- match.arg(thresholds)
  if (!is.data.frame(data))
    stop("`object_a` and `object_b` must name columns in `data`")
  .check_reshape_column(data, object_a, "object_a")
  .check_reshape_column(data, object_b, "object_b")
  if (!is.data.frame(predictors) || !"object" %in% names(predictors))
    stop("`predictors` must be a data frame with an `object` column")
  .check_column_names(predictors)
  observed_objects <- unique(c(.role_text_values(data[[object_a]]),
                               .role_text_values(data[[object_b]])))
  observed_objects <- observed_objects[!is.na(observed_objects) &
                                         nzchar(observed_objects)]
  predictor_objects <- .role_text_values(predictors$object)
  if (anyNA(predictor_objects) || any(!nzchar(predictor_objects)))
    stop("predictor object names must be non-missing and non-empty")
  unknown <- setdiff(predictor_objects, observed_objects)
  if (length(unknown))
    stop("predictor object(s) are not present in the comparison data: ",
         paste(unknown, collapse = ", "))
  args <- list(object_a = object_a, object_b = object_b, winner = winner,
               response = response, margin = margin, judge = judge,
               count = count, order = order, position = position,
               ties = ties, thresholds = thresholds, maxit = maxit, tol = tol)
  reference <- do.call(btl, c(list(data = data), args))
  if (!isTRUE(reference$converged))
    stop("the unrestricted comparative judgement calibration did not ",
         "converge; an explanatory restriction cannot be assessed against it",
         call. = FALSE)
  usable_objects <- unique(c(
    as.character(reference$observed_comparisons$object_a),
    as.character(reference$observed_comparisons$object_b)))
  usable_objects <- usable_objects[!is.na(usable_objects) & nzchar(usable_objects)]
  # objects set aside at a response boundary are not fitted, so centring the
  # design over them would leave the fitted locations off the sum-zero
  # origin the model and the reference fit both use
  all_objects <- as.character(reference$objects$object)
  calibrated <- all_objects
  if (!is.null(reference$objects$extreme))
    calibrated <- calibrated[!reference$objects$extreme %in% TRUE]
  design <- .btl_explanatory_design(predictors, formula, calibrated,
                                    known = usable_objects)
  fit <- do.call(btl, c(list(data = data), args,
                        list(.object_design = design)))
  fit$reference_fit <- reference
  fit$explanatory <- list(
    formula = formula, formula_text = paste(deparse(formula), collapse = " "),
    metadata = design$metadata, source_predictors = design$source_predictors,
    base_B = design$B, active_B = fit$location_design,
    relaxations = data.frame(), refit_spec = list(data = data, args = args))
  class(fit) <- c("rasch_btl_explanatory", "rasch_btl")
  fit
}

#' Add a fixed object departure to an explanatory comparative judgement model
#'
#' @param fit A fitted object from \code{\link{btl_explanatory}}.
#' @param object Object name.
#' @return A refitted explanatory comparative judgement model.
#' @export
relax_btl_explanatory <- function(fit, object) {
  if (!inherits(fit, "rasch_btl_explanatory"))
    stop("relax_btl_explanatory() needs an explanatory comparative judgement fit")
  if (!isTRUE(fit$converged))
    stop("the explanatory comparative judgement calibration did not converge; it cannot be relaxed",
         call. = FALSE)
  if (!is.atomic(object) || !is.null(dim(object)) || length(object) != 1L ||
      is.na(object))
    stop("`object` must name exactly one object")
  object <- .role_text_values(object)
  if (!nzchar(object)) stop("`object` must name exactly one object")
  obj_tab <- fit$objects
  if ("extreme" %in% names(obj_tab)) {
    at <- match(object, obj_tab$object)
    if (!is.na(at) && isTRUE(obj_tab$extreme[at]))
      .refuse(object, " was set aside at a response boundary; its location ",
              "is an extrapolation for display and cannot be relaxed")
    obj_tab <- obj_tab[!(obj_tab$extreme %in% TRUE), ]
  }
  objects <- as.character(obj_tab$object)
  j <- match(object, objects)
  if (is.na(j)) stop("object not found in the explanatory fit: ", object)
  D <- diag(length(objects))[, j, drop = FALSE]
  D <- (diag(length(objects)) - 1 / length(objects)) %*% D
  rownames(D) <- objects; colnames(D) <- paste0("departure[", object, "]")
  B <- fit$explanatory$active_B
  if (qr(cbind(B, D), tol = 1e-10)$rank == qr(B, tol = 1e-10)$rank)
    stop("that object departure is already represented by the active model")
  rel <- fit$explanatory$relaxations
  rel <- rbind(rel, data.frame(order = nrow(rel) + 1L, object = object,
                               component = "Object location",
                               parameters_added = 1L,
                               stringsAsFactors = FALSE))
  out <- .btl_explanatory_refit(fit, cbind(B, D), rel)
  if (!isTRUE(out$converged))
    stop("the relaxed explanatory comparative judgement calibration did not converge",
         call. = FALSE)
  out
}

#' @export
print.rasch_btl_explanatory <- function(x, ...) {
  cat("Explanatory comparative judgement model\n")
  cat("Formula: ", x$explanatory$formula_text, "\n", sep = "")
  print.rasch_btl(x, ...)
}

.explanatory_attach <- function(out, reference, design, formula, level,
                                relaxations, n_groups_requested,
                                maxit, tol) {
  out$reference_fit <- reference
  out$mc <- reference$mc
  out$explanatory <- list(
    formula = formula, formula_text = paste(deparse(formula), collapse = " "),
    level = level, metadata = design$metadata,
    model_matrix = design$matrix, base_B = design$B,
    active_B = out$est$B, threshold_index = design$threshold_index,
    source_predictors = design$source_predictors,
    relaxations = relaxations)
  out$refit_spec <- list(model = "PCM", n_groups = n_groups_requested,
    anchors = NULL,
    na_codes = reference$refit_spec$na_codes %||% -1,
    key = reference$refit_spec$key,
    pc_components = NULL, maxit = maxit, tol = tol,
    explanatory = TRUE)
  class(out) <- c("rasch_explanatory", "rasch")
  out
}

#' Fit an explanatory Rasch model
#'
#' Fits the linear logistic test model (LLTM) for dichotomous responses or
#' the linear partial credit model (LPCM) for polytomous responses. Item or
#' threshold locations are linear functions of observed predictors. The
#' response model remains Rasch and is estimated by pairwise conditional
#' maximum likelihood.
#'
#' @details
#' For threshold \eqn{k} of item \eqn{i},
#' \deqn{\delta_{ik}=z_{ik}^{T}\gamma.}
#' The adjacent-category log odds are
#' \deqn{\log\{P(X_{ni}=k)/P(X_{ni}=k-1)\}=\theta_n-\delta_{ik}.}
#' The threshold origin is fixed to the same mean-item-location zero used by
#' \code{\link{rasch}}. An intercept therefore sets the arbitrary origin and
#' is not separately estimated. Numeric predictors are continuous, unordered
#' factors are categorical, and ordered factors use successive contrasts
#' between adjacent levels. Character predictors are converted to unordered
#' factors. The reserved factor
#' \code{threshold} identifies the within-item threshold number;
#' \code{threshold_number} supplies its integer value.
#' Design columns are centred and rescaled internally for numerical stability; reported
#' coefficients and standard errors use the supplied predictor units.
#' Coincident coefficient labels receive numeric suffixes; this does not
#' change the predictor design.
#'
#' A free PCM reference is fitted to the same prepared responses and retained
#' on the object. \code{\link{explanatory_test}} applies the first-order Kent
#' calibration required for the pairwise composite likelihood. When an
#' identifier occurs on more than one response row, coefficient covariance is
#' clustered by person. A linearised delete-one-person correction accounts for
#' finite-cluster leverage without refitting the model once per person.
#' Supported repeated-person fits use a \eqn{t} reference
#' with degrees of freedom equal to the number of person clusters contributing
#' conditional information minus one; inference
#' is withheld when the calibration lacks enough independent information.
#' Supported fits without repeated identifiers use the limiting normal
#' reference. Holm adjustment covers the coefficient family.
#' With few persons and unequal numbers of response rows, these approximate
#' tests can still be mildly liberal; the correction does not guarantee nominal
#' coverage in small samples.
#'
#' @param data,items,id,factors,n_groups,na_codes,key,maxit,tol As in
#'   \code{\link{rasch}}.
#' @param predictors Data frame containing an \code{item} column and the
#'   predictors named in \code{formula}. With \code{level = "threshold"}, it
#'   must also contain \code{threshold}, with one row for every fitted item
#'   threshold. Column names must be unique.
#' @param formula One-sided explanatory formula. For example,
#'   \code{~ format + operation + format:operation}. The reserved
#'   \code{threshold} factor permits threshold-specific effects.
#'   Formula offsets (\code{offset()}) are not supported.
#' @param level Whether \code{predictors} contains one row per \code{"item"}
#'   or per \code{"threshold"}. Item rows are expanded over their thresholds.
#' @return An object of class \code{"rasch_explanatory"} inheriting from
#'   \code{"rasch"}. Standard item, person, fit and diagnostic components use
#'   the explanatory thresholds. The \code{explanatory} component contains the
#'   formula, metadata and design matrices; \code{reference_fit} is the free
#'   PCM calibration. \code{est$coefficients} reports the estimates, standard
#'   errors, \eqn{t} statistics, reference degrees of freedom, raw
#'   probabilities and Holm-adjusted probabilities.
#' @references
#' Fischer, G. H. (1973). The linear logistic test model as an instrument in
#' educational research. Acta Psychologica, 37, 359--374.
#'
#' Fischer, G. H. and Ponocny, I. (1994). An extension of the partial credit
#' model with an application to the measurement of change. Psychometrika, 59,
#' 177--192.
#' @examples
#' set.seed(1)
#' q <- data.frame(item = paste0("I", 1:8),
#'                 operation = rep(0:1, each = 4),
#'                 format = rep(c("A", "B"), 4))
#' difficulty <- -1 + 0.7 * q$operation + 0.4 * (q$format == "B")
#' X <- matrix(rbinom(500 * 8, 1,
#'   plogis(outer(rnorm(500), difficulty, "-"))), 500, 8)
#' colnames(X) <- q$item
#' fit <- rasch_explanatory(X, predictors = q,
#'                          formula = ~ operation + format)
#' fit$est$coefficients
#' explanatory_test(fit)
#' @seealso \code{\link{explanatory_test}},
#'   \code{\link{explanatory_diagnostics}}, and
#'   \code{\link{relax_explanatory}}.
#' @export
rasch_explanatory <- function(data, predictors, formula, items = NULL,
                              level = c("item", "threshold"), id = NULL,
                              factors = NULL, n_groups = NULL,
                              na_codes = -1, key = NULL, maxit = 60,
                              tol = 1e-8) {
  level <- match.arg(level)
  reference <- rasch(data, model = "PCM", id = id, factors = factors,
                     items = items, n_groups = n_groups,
                     na_codes = na_codes, key = key, maxit = maxit, tol = tol)
  if (!isTRUE(reference$est$converged))
    stop("the unrestricted Rasch calibration did not converge; an ",
         "explanatory restriction cannot be assessed against it",
         call. = FALSE)
  design <- .explanatory_metadata(predictors, formula, reference$X, level)
  est <- .pcml_design(reference$X, design$B,
                      parameter_names = colnames(design$B),
                      maxit = maxit, tol = tol,
                      cluster = reference$person$id)
  if (!isTRUE(est$converged))
    warning("the explanatory calibration did not converge; estimates, ",
            "diagnostics and probabilities are unreliable", call. = FALSE)
  kind <- if (all(reference$m == 1L)) "LLTM" else "LPCM"
  notes <- c(reference$notes,
             sprintf("%s explanatory calibration: %s", kind,
                     paste(deparse(formula), collapse = " ")),
             est$notes)
  out <- .assemble_fit("PCM", reference$X, est, reference$person$id,
                       reference$factors, n_groups, notes)
  out$explanatory_model <- kind
  .explanatory_attach(out, reference, design, formula, level,
                      relaxations = data.frame(),
                      n_groups_requested = n_groups,
                      maxit = maxit, tol = tol)
}

#' Compare an explanatory model with its free calibration
#'
#' Tests explanatory item, threshold or object restrictions against the
#' corresponding free calibration of the same responses. The inferential
#' result uses the first-order Kent calibration for the fitted likelihood and
#' sandwich covariance. This multivariate comparison is asymptotic; unlike
#' the individual coefficient tests, it has no finite-person-cluster \eqn{t}
#' correction. The calibration coefficient of determination is
#' \deqn{R^2_{cal}=1-\frac{\sum_j(\hat\eta^{free}_j-
#' \hat\eta^{expl}_j-\bar d)^2}{\sum_j(\hat\eta^{free}_j-
#' \bar\eta^{free})^2},}
#' where \eqn{\bar d} removes the arbitrary scale origin. It describes the
#' proportion of variation in the well-determined free threshold calibration
#' (Rasch models) or free object calibration (comparative judgement)
#' reproduced by the explanatory model. It is at most one and may be negative.
#' It is not adjusted for the number of predictors, so with few calibrated
#' parameters it reads above zero even for an uninformative design.
#' \code{r_squared_adj} divides the unexplained proportion by its share of
#' the degrees of freedom, \eqn{1-(1-R^2_{cal})(n-1)/\mathit{df}}, where
#' \eqn{n} counts the calibrated parameters compared and \eqn{\mathit{df}}
#' is \eqn{n} minus the rank of the retained explanatory design with its
#' origin, so exclusions that remove a level's only support reduce it. The
#' correction is exact for independent homoskedastic estimates fitted by
#' least squares, which these calibrations are not, so read it as a
#' descriptive optimism adjustment. Read either beside the test rather than
#' in place of it.
#'
#' @param fit A fitted explanatory Rasch or comparative judgement model.
#' @return A one-row data frame containing the raw and Kent-calibrated
#'   statistics, degrees of freedom and parameter counts. The primary
#'   \code{p} and the retained \code{p_kent} are the Kent-calibrated
#'   probability. \code{p_naive} is the unscaled composite-likelihood
#'   probability and is provided for methodological inspection, not
#'   inference. \code{r_squared} is the calibration coefficient of
#'   determination, \code{r_squared_adj} its degrees-of-freedom-adjusted
#'   counterpart, and \code{r2_basis} names the calibrated parameters used.
#' @export
explanatory_test <- function(fit) {
  if (inherits(fit, "rasch_btl_explanatory")) {
    if (!isTRUE(fit$converged))
      stop("the explanatory comparative judgement calibration did not converge; model comparison is unavailable",
           call. = FALSE)
    if (is.null(fit$reference_fit) || !isTRUE(fit$reference_fit$converged))
      stop("the unrestricted comparative judgement reference fit did not ",
           "converge; model comparison is unavailable", call. = FALSE)
    z <- .btl_explanatory_nested_test(fit$reference_fit, fit)
    free <- fit$reference_fit$objects
    active <- fit$objects
    # an extrapolated boundary row is not a calibrated location; the design
    # rows span the calibrated objects only
    if ("extreme" %in% names(active)) active <- active[!(active$extreme %in% TRUE), ]
    if ("extreme" %in% names(free)) free <- free[!(free$extreme %in% TRUE), ]
    f <- free$location[match(active$object, free$object)]
    a <- active$location
    ok <- is.finite(f) & is.finite(a)
    d <- f[ok] - a[ok]
    den <- sum((f[ok] - mean(f[ok]))^2)
    r2 <- if (sum(ok) > 1L && is.finite(den) && den > 0)
      1 - sum((d - mean(d))^2) / den else NA_real_
    df_sub <- if (sum(ok) > 1L) sum(ok) -
      qr(cbind(1, fit$location_design[ok, , drop = FALSE]))$rank else 0L
    r2_adj <- if (is.finite(r2) && df_sub > 0)
      1 - (1 - r2) * (sum(ok) - 1) / df_sub else NA_real_
    out <- data.frame(
      model = if (nrow(fit$explanatory$relaxations))
        "Partially relaxed explanatory CJ model" else "Explanatory CJ",
      parameters = ncol(fit$location_design),
      free_parameters = ncol(fit$reference_fit$location_design),
      r_squared = r2, r_squared_adj = r2_adj,
      r2_basis = "object calibration",
      chisq = z$chisq, df = z$df, p_naive = z$p,
      chisq_kent = z$chisq_kent, p = z$p_kent, p_kent = z$p_kent,
      stringsAsFactors = FALSE)
    attr(out, "lambda") <- z$lambda
    return(.tag_tables(out))
  }
  if (!inherits(fit, "rasch_explanatory"))
    stop("explanatory_test() needs an explanatory Rasch fit")
  if (!isTRUE(fit$est$converged))
    stop("the explanatory calibration did not converge; model comparison is unavailable",
         call. = FALSE)
  if (is.null(fit$reference_fit) ||
      !isTRUE(fit$reference_fit$est$converged))
    stop("the unrestricted Rasch reference fit did not converge; model ",
         "comparison is unavailable", call. = FALSE)
  z <- .pcml_nested_test(fit$reference_fit$est, fit$est)
  free <- fit$reference_fit$est$thr$tau
  active <- fit$est$thr$tau
  ok <- is.finite(free) & is.finite(active)
  free_weak <- fit$reference_fit$est$thr$weak %||% rep(FALSE, length(free))
  active_weak <- fit$est$thr$weak %||% rep(FALSE, length(active))
  ok <- ok & !free_weak & !active_weak
  d <- free[ok] - active[ok]
  den <- sum((free[ok] - mean(free[ok]))^2)
  r2 <- if (sum(ok) > 1L && is.finite(den) && den > 0)
    1 - sum((d - mean(d))^2) / den else NA_real_
  df_sub <- if (sum(ok) > 1L)
    sum(ok) - qr(cbind(1, fit$est$B[ok, , drop = FALSE]))$rank else 0L
  r2_adj <- if (is.finite(r2) && df_sub > 0)
    1 - (1 - r2) * (sum(ok) - 1) / df_sub else NA_real_
  out <- data.frame(
    model = if (nrow(fit$explanatory$relaxations))
      "Partially relaxed explanatory model" else fit$explanatory_model,
    parameters = fit$est$n_parameters,
    free_parameters = fit$reference_fit$est$n_parameters,
    r_squared = r2, r_squared_adj = r2_adj,
    r2_basis = "threshold calibration",
    chisq = z$chisq, df = z$df, p_naive = z$p,
    chisq_kent = z$chisq_kent, p = z$p_kent, p_kent = z$p_kent,
    stringsAsFactors = FALSE)
  attr(out, "lambda") <- z$lambda
  .tag_tables(out)
}

.explanatory_candidate <- function(fit, item,
                                   component = c("location", "thresholds")) {
  component <- match.arg(component)
  j <- match(item, colnames(fit$X))
  if (is.na(j)) stop("item not found in the explanatory fit: ", item)
  thr <- fit$explanatory$threshold_index %||%
    threshold_index(fit$m)
  rows <- which(thr$item == j); mi <- length(rows); M <- nrow(thr)
  if (component == "location") {
    D <- matrix(0, M, 1L); D[rows, 1L] <- 1
    colnames(D) <- paste0("departure_location[", item, "]")
  } else {
    if (mi < 2L) stop("a dichotomous item has no threshold-structure block")
    D <- matrix(0, M, mi - 1L)
    D[rows, ] <- rbind(diag(mi - 1L), rep(-1, mi - 1L))
    colnames(D) <- paste0("departure_threshold[", item, ",",
                          seq_len(mi - 1L), "]")
  }
  A <- .explanatory_projector(fit$m, thr)
  A %*% D
}

.explanatory_addable <- function(B, D) {
  b_scale <- .design_column_scale(B)
  augmented <- cbind(B, D)
  a_scale <- .design_column_scale(augmented)
  qr(sweep(augmented, 2L, a_scale, `/`), tol = 1e-10)$rank -
    qr(sweep(B, 2L, b_scale, `/`), tol = 1e-10)$rank
}

# Keep a stable subset of a candidate block that adds exactly its remaining
# directions to the active design. A predictor may already span only part of
# a polytomous item's threshold block; appending the whole block would then
# make the refit rank deficient even though a genuine departure remains.
.explanatory_addition <- function(B, D) {
  # Rank decisions are made after column normalisation, while the returned
  # residuals retain the caller's units. This separates a nearly represented
  # added direction before the refitter scales its columns. Keeping the
  # residual in D's units also keeps the reported departure coefficient on
  # the same scale as the requested fixed departure.
  current <- B
  current_rank <- qr(sweep(current, 2L, .design_column_scale(current), `/`),
                     tol = 1e-10)$rank
  stable <- matrix(numeric(0), nrow = nrow(D), ncol = 0L,
                   dimnames = list(rownames(D), NULL))
  raw_map <- matrix(numeric(0), nrow = ncol(D), ncol = 0L)
  keep <- integer(0)
  for (j in seq_len(ncol(D))) {
    d <- D[, j, drop = FALSE]
    ds <- .design_column_scale(d)[1L]
    d_work <- d / ds
    current_scale <- .design_column_scale(current)
    current_work <- sweep(current, 2L, current_scale, `/`)
    q_current <- qr(current_work, tol = 1e-10)
    candidate_rank <- qr(cbind(current_work, d_work),
                         tol = 1e-10)$rank
    if (candidate_rank > current_rank) {
      # This residual spans the same augmented model as d. It is deliberately
      # not normalised here: the coefficient remains the requested
      # departure, while the downstream fit scales the column internally.
      residual <- qr.resid(q_current, d)
      coef <- qr.coef(q_current, d)
      raw <- numeric(ncol(D)); raw[j] <- 1
      if (ncol(raw_map)) {
        n_active <- ncol(B)
        previous <- coef[seq.int(n_active + 1L, length(coef))] /
          current_scale[seq.int(n_active + 1L, length(current_scale))]
        raw <- raw - drop(raw_map %*% previous)
      }
      stable <- cbind(stable, residual)
      raw_map <- cbind(raw_map, raw)
      keep <- c(keep, j)
      current <- cbind(current, residual)
      current_rank <- candidate_rank
    }
  }
  colnames(stable) <- colnames(D)[keep]
  stable <- .explanatory_unique_columns(stable, colnames(B))
  attr(stable, "raw_map") <- raw_map
  stable
}

#' Diagnose fixed departures from an explanatory model
#'
#' Fits each available item-location, polytomous threshold-structure or
#' comparative-judgement object departure separately from the active model.
#' Probabilities use Kent calibration and Holm adjustment over the complete
#' candidate family. The Kent departure tests are first-order asymptotic
#' comparisons and do not use the finite-person-cluster correction applied to
#' individual coefficients.
#' A candidate with a withheld probability, including a refit that errors or
#' fails to converge, remains in that family.
#'
#' @param fit A fitted explanatory Rasch or comparative judgement model.
#' @param p_adjust Multiplicity adjustment over the candidate departures.
#' @return A data frame ordered by adjusted probability. For item fits, a
#'   \code{weak} column marks items whose thresholds the calibration flags
#'   as weakly identified; their probabilities are withheld, since the
#'   departure test rests on the same sparse categories, and a note on the
#'   table records the withholding. The \code{converged} column identifies
#'   candidate refits that converged. Statistics from a failed or
#'   non-convergent candidate are withheld, but it remains in the
#'   multiplicity family.
#' @export
explanatory_diagnostics <- function(fit, p_adjust = "holm") {
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      !is.null(dim(p_adjust)) || !is.null(oldClass(p_adjust)) ||
      is.na(p_adjust))
    stop("`p_adjust` must name one method in stats::p.adjust.methods")
  if (inherits(fit, "rasch_btl_explanatory")) {
    if (!isTRUE(fit$converged))
      stop("the explanatory comparative judgement calibration did not converge; diagnostics are unavailable",
           call. = FALSE)
    if (!p_adjust %in% stats::p.adjust.methods)
      stop("p_adjust must name a method in stats::p.adjust.methods")
    obj_tab <- fit$objects
    if ("extreme" %in% names(obj_tab)) obj_tab <- obj_tab[!(obj_tab$extreme %in% TRUE), ]
    objects <- as.character(obj_tab$object)
    B <- fit$explanatory$active_B
    rows <- list()
    for (j in seq_along(objects)) {
      D <- diag(length(objects))[, j, drop = FALSE]
      D <- (diag(length(objects)) - 1 / length(objects)) %*% D
      rownames(D) <- objects
      colnames(D) <- paste0("departure[", objects[j], "]")
      D <- tryCatch(.explanatory_addition(B, D),
                    error = function(e) e)
      if (inherits(D, "error")) {
        rows[[length(rows) + 1L]] <- data.frame(
          object = objects[j], component = "Object location",
          parameters_added = 1L, departure = NA_real_,
          deviance_reduction = NA_real_, df = NA_integer_, p = NA_real_,
          converged = FALSE, stringsAsFactors = FALSE)
        next
      }
      if (!ncol(D))
        next
      cand <- tryCatch(
        .btl_explanatory_refit(fit, cbind(B, D),
                               fit$explanatory$relaxations),
        error = function(e) e)
      converged <- !inherits(cand, "error") && isTRUE(cand$converged)
      if (inherits(cand, "error")) {
        tst <- list(chisq = NA_real_, df = NA_integer_, p_kent = NA_real_)
      } else if (converged) {
        tst <- tryCatch(.btl_explanatory_nested_test(cand, fit),
                        error = function(e) {
                          converged <<- FALSE
                          list(chisq = NA_real_, df = NA_integer_,
                               p_kent = NA_real_)
                        })
      } else {
        tst <- list(chisq = NA_real_, df = NA_integer_, p_kent = NA_real_)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        object = objects[j], component = "Object location",
        parameters_added = 1L,
        departure = if (converged && !inherits(cand, "error"))
          utils::tail(cand$object_coefficients$estimate, 1L) else NA_real_,
        deviance_reduction = tst$chisq, df = tst$df, p = tst$p_kent,
        converged = converged, stringsAsFactors = FALSE)
    }
    if (!length(rows)) return(.tag_tables(data.frame(
      object = character(0), component = character(0),
      parameters_added = integer(0), departure = numeric(0),
      deviance_reduction = numeric(0), df = integer(0), p = numeric(0),
      converged = logical(0), p_adj = numeric(0))))
    out <- do.call(rbind, rows)
    out$p_adj <- .p_adjust_family(out$p, method = p_adjust)
    out <- out[order(out$p_adj, -out$deviance_reduction), , drop = FALSE]
    rownames(out) <- NULL
    attr(out, "p_adjust") <- p_adjust
    if (any(!out$converged))
      attr(out, "note") <- paste("statistics are withheld for failed or",
        "non-convergent candidate refits; those candidates remain in the",
        "adjustment family")
    return(.tag_tables(out))
  }
  if (!inherits(fit, "rasch_explanatory"))
    stop("explanatory_diagnostics() needs an explanatory Rasch fit")
  if (!isTRUE(fit$est$converged))
    stop("the explanatory calibration did not converge; diagnostics are unavailable",
         call. = FALSE)
  if (!p_adjust %in% stats::p.adjust.methods)
    stop("p_adjust must name a method in stats::p.adjust.methods")
  B <- fit$est$B; rows <- list()
  spec <- fit$refit_spec
  # Candidate construction, refitting and nested testing are deliberately
  # isolated below. One singular candidate must remain an unavailable member
  # of the Holm family rather than aborting all later departures.
  for (item in colnames(fit$X)) for (component in c("location", "thresholds")) {
    ii <- match(item, colnames(fit$X))
    nominal_add <- if (component == "location") 1L else
      max(fit$m[ii] - 1L, 0L)
    if (component == "thresholds" && fit$m[ii] < 2L)
      next
    D_raw <- tryCatch(.explanatory_candidate(fit, item, component),
      error = function(e) e)
    D <- if (inherits(D_raw, "error")) D_raw else
      tryCatch(.explanatory_addition(B, D_raw),
      error = function(e) e)
    if (inherits(D, "error")) {
      add <- nominal_add
      est <- D
    } else {
      add <- ncol(D)
      if (!add) next
      candB <- cbind(B, D)
      est <- tryCatch(.pcml_design(fit$X, candB, colnames(candB),
                                  maxit = spec$maxit, tol = spec$tol,
                                  cluster = fit$person$id),
                      error = function(e) e)
    }
    converged <- !inherits(est, "error") && isTRUE(est$converged)
    if (inherits(est, "error")) {
      tst <- list(chisq = NA_real_, df = NA_integer_, p_kent = NA_real_)
    } else if (converged) {
      tst <- tryCatch(.pcml_nested_test(est, fit$est),
                      error = function(e) {
                        converged <<- FALSE
                        list(chisq = NA_real_, df = NA_integer_,
                             p_kent = NA_real_)
                      })
    } else {
      tst <- list(chisq = NA_real_, df = NA_integer_, p_kent = NA_real_)
    }
    departure <- NA_real_
    if (converged && !inherits(est, "error")) {
      b <- utils::tail(est$beta, ncol(D))
      raw_map <- attr(D, "raw_map")
      b_raw <- if (is.null(raw_map)) b else drop(raw_map %*% b)
      departure <- if (component == "location") unname(b_raw[1L]) else
        max(abs(drop(D_raw %*% b_raw)))
    }
    # a departure test rests on the same sparse categories that made the
    # item's thresholds weak; the probability is withheld there, as the
    # threshold standard errors already are, and the departure stays
    # descriptive
    weak_item <- isTRUE(any(fit$thresholds$weak[fit$thresholds$item == ii]))
    rows[[length(rows) + 1L]] <- data.frame(
      item = item, component = if (component == "location")
        "Item location" else "Threshold structure",
      parameters_added = add, departure = departure,
      deviance_reduction = tst$chisq, df = tst$df,
      p = if (weak_item) NA_real_ else tst$p_kent,
      weak = weak_item, converged = converged,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(.tag_tables(data.frame(
    item = character(0), component = character(0), parameters_added = integer(0),
    departure = numeric(0), deviance_reduction = numeric(0), df = integer(0),
    p = numeric(0), weak = logical(0), converged = logical(0),
    p_adj = numeric(0))))
  out <- do.call(rbind, rows)
  usable <- is.finite(out$p)
  out$p_adj <- NA_real_
  out$p_adj[usable] <- stats::p.adjust(
    out$p[usable], method = p_adjust, n = nrow(out))
  notes <- character(0)
  if (any(out$weak))
    notes <- c(notes, paste("departure probabilities are withheld for",
      "item(s) with weak thresholds:",
      paste(unique(out$item[out$weak]), collapse = ", ")))
  if (any(!out$converged))
    notes <- c(notes, paste("statistics are withheld for failed or",
      "non-convergent candidate refits; those candidates remain in the",
      "adjustment family"))
  if (length(notes)) attr(out, "note") <- paste(notes, collapse = "; ")
  out <- out[order(out$p_adj, -out$deviance_reduction), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "p_adjust") <- p_adjust
  .tag_tables(out)
}

#' Relax a nominated explanatory restriction
#'
#' Adds either one fixed item-location departure or the part of an item's
#' threshold-structure block not already represented by the predictor design,
#' then repeats the complete conditional calibration and downstream Rasch
#' analysis. The departure is fixed rather than random; raw-score sufficiency
#' and the common discrimination remain.
#' Earlier DIF splits and superitem definitions are retained.
#'
#' @param fit A fitted explanatory Rasch model.
#' @param item Item name.
#' @param component Either \code{"location"} or \code{"thresholds"}.
#' @return A partially relaxed \code{"rasch_explanatory"} fit.
#' @export
relax_explanatory <- function(fit, item,
                              component = c("location", "thresholds")) {
  if (!inherits(fit, "rasch_explanatory"))
    stop("relax_explanatory() needs an explanatory Rasch fit")
  if (!isTRUE(fit$est$converged))
    stop("the explanatory calibration did not converge; it cannot be relaxed",
         call. = FALSE)
  if (!is.atomic(item) || !is.null(dim(item)) || length(item) != 1L ||
      is.na(item))
    stop("`item` must name exactly one item")
  item <- .explanatory_match_items(item, colnames(fit$X))
  if (!nzchar(trimws(item))) stop("`item` must name exactly one item")
  component <- match.arg(component)
  D <- .explanatory_addition(
    fit$est$B, .explanatory_candidate(fit, item, component))
  add <- ncol(D)
  if (!add)
    stop("that departure is already represented by the active explanatory model")
  B <- cbind(fit$est$B, D)
  spec <- fit$refit_spec
  est <- .pcml_design(fit$X, B, colnames(B), maxit = spec$maxit,
                      tol = spec$tol, cluster = fit$person$id)
  if (!isTRUE(est$converged))
    stop("the relaxed explanatory calibration did not converge")
  rel <- fit$explanatory$relaxations
  new <- data.frame(order = nrow(rel) + 1L, item = item,
                    component = if (component == "location")
                      "Item location" else "Threshold structure",
                    parameters_added = add, stringsAsFactors = FALSE)
  rel <- rbind(rel, new)
  out <- .assemble_fit("PCM", fit$X, est, fit$person$id, fit$factors,
                       .refit_n_groups(fit),
                       unique(c(fit$notes, est$notes,
                         sprintf("fixed explanatory departure: %s, %s",
                                 item, tolower(new$component)))))
  out$explanatory_model <- fit$explanatory_model
  design <- list(B = fit$explanatory$base_B,
                 matrix = fit$explanatory$model_matrix,
                 metadata = fit$explanatory$metadata,
                 source_predictors = fit$explanatory$source_predictors,
                 threshold_index = fit$explanatory$threshold_index)
  out <- .explanatory_attach(out, fit$reference_fit, design,
                      fit$explanatory$formula, fit$explanatory$level, rel,
                      n_groups_requested = .refit_n_groups(fit),
                      maxit = spec$maxit, tol = spec$tol)
  out$mc <- fit$mc
  # Relaxation changes restrictions, not response columns or their history.
  out$split_map <- fit$split_map
  out$subtest_map <- fit$subtest_map
  out$subtest_binary <- fit$subtest_binary
  out
}

.explanatory_inherit_mc <- function(fit, source, inherit,
                                    exclude = character(0),
                                    person_rows = seq_len(nrow(source))) {
  if (is.null(fit$mc) || is.null(fit$mc$raw)) return(NULL)
  items <- colnames(source)
  old <- unname(inherit[items])
  keep <- items[old %in% colnames(fit$mc$raw) & !items %in% exclude]
  if (!length(keep)) return(NULL)
  raw <- vapply(keep, function(it) {
    old_item <- unname(inherit[it])
    old_score <- fit$X[person_rows, old_item]
    observed <- !is.na(source[, it])
    if (any(observed & (is.na(old_score) | source[, it] != old_score)))
      stop("cannot inherit observed multiple-choice answers for changed ",
           "scores; simulated refits must use inherit_mc = FALSE", call. = FALSE)
    value <- fit$mc$raw[person_rows, old_item]
    value[is.na(source[, it])] <- NA_character_
    value
  }, character(nrow(source)))
  if (!is.matrix(raw)) raw <- matrix(raw, ncol = length(keep))
  colnames(raw) <- keep
  map <- lapply(keep, function(it) fit$mc$map[[unname(inherit[it])]])
  names(map) <- keep
  list(raw = raw, map = map,
       key = vapply(map, .key_label, character(1)))
}

.explanatory_refit_modified <- function(fit, source, inherit = NULL,
                                        location_relaxed = character(0),
                                        fully_relaxed = character(0),
                                        person_rows = NULL, inherit_mc = TRUE) {
  source <- as.matrix(source)
  .check_flag(inherit_mc, "inherit_mc")
  # a source holding a SUBSET of the fitted persons must carry the matching
  # identifiers and person factors: passing the full-length vectors would
  # fail on row alignment, and every caller catching that error would report
  # a failure whose cause is invisible
  if (is.null(person_rows)) person_rows <- seq_len(nrow(source))
  if (length(person_rows) != nrow(source))
    stop("`person_rows` must give one fitted row per row of `source`")
  if (!is.numeric(person_rows) || is.complex(person_rows) ||
      !is.null(dim(person_rows)) || !is.null(oldClass(person_rows)) ||
      any(!is.finite(person_rows)) || any(person_rows != floor(person_rows)) ||
      any(person_rows < 1L | person_rows > nrow(fit$X)))
    stop("`person_rows` must contain valid fitted row indices")
  new_items <- colnames(source)
  old_items <- colnames(fit$X)
  if (is.null(inherit))
    inherit <- stats::setNames(intersect(new_items, old_items),
                               intersect(new_items, old_items))
  if (is.null(names(inherit)) || any(!new_items %in% names(inherit)))
    stop("the explanatory refit needs an inherited predictor source for every item")
  if (any(!unname(inherit[new_items]) %in% old_items))
    stop("an inherited explanatory item is absent from the original fit")
  src <- fit$explanatory$source_predictors
  level <- fit$explanatory$level
  pred <- list()
  m_new <- apply(source, 2L, max, na.rm = TRUE)
  for (it in new_items) {
    old <- unname(inherit[it])
    z <- src[as.character(src$item) == old, , drop = FALSE]
    if (!nrow(z)) stop("predictors are unavailable for inherited item ", old)
    if (level == "item") {
      z <- z[1L, , drop = FALSE]; z$item <- it
    } else {
      kn <- if ("threshold" %in% names(z))
        as.integer(as.character(z$threshold)) else
          as.integer(as.character(z$threshold_number))
      take <- vapply(seq_len(m_new[it]), function(k) {
        hit <- which(kn == k)
        if (length(hit)) hit[1L] else which.max(kn)
      }, integer(1))
      z <- z[take, , drop = FALSE]; z$item <- it
      if ("threshold" %in% names(z)) z$threshold <- seq_len(m_new[it])
      if ("threshold_number" %in% names(z))
        z$threshold_number <- seq_len(m_new[it])
    }
    pred[[it]] <- z
  }
  pred <- do.call(rbind, pred); rownames(pred) <- NULL
  spec <- fit$refit_spec
  out <- rasch_explanatory(source, predictors = pred,
    formula = fit$explanatory$formula, level = level,
    id = fit$person$id[person_rows],
    factors = if (is.null(fit$factors)) NULL else
      fit$factors[person_rows, , drop = FALSE],
    n_groups = .refit_n_groups(fit), maxit = spec$maxit %||% 60,
    tol = spec$tol %||% 1e-8)
  expected_max <- stats::setNames(fit$m[match(inherit[new_items], old_items)],
                                  new_items)
  expected_max[intersect(fully_relaxed, new_items)] <-
    m_new[intersect(fully_relaxed, new_items)]
  .require_fitted_score_structure(out, expected_max, "the explanatory refit")

  # Preserve prior analyst-approved departures on items that survive or are
  # replaced by inherited copies. A split copy receives the source item's
  # departure before its own location is freed.
  rel <- fit$explanatory$relaxations
  if (nrow(rel)) for (r in seq_len(nrow(rel))) {
    targets <- names(inherit)[unname(inherit) == rel$item[r]]
    for (it in targets) {
      component <- if (rel$component[r] == "Item location")
        "location" else "thresholds"
      D <- tryCatch(.explanatory_candidate(out, it, component),
                    error = function(e) NULL)
      if (!is.null(D) && .explanatory_addable(out$est$B, D) > 0L)
        out <- relax_explanatory(out, it, component)
    }
  }
  for (it in unique(location_relaxed)) {
    D <- .explanatory_candidate(out, it, "location")
    if (.explanatory_addable(out$est$B, D) > 0L)
      out <- relax_explanatory(out, it, "location")
  }
  for (it in unique(fully_relaxed)) {
    D <- .explanatory_candidate(out, it, "location")
    if (.explanatory_addable(out$est$B, D) > 0L)
      out <- relax_explanatory(out, it, "location")
    if (out$m[match(it, colnames(out$X))] > 1L) {
      D <- .explanatory_candidate(out, it, "thresholds")
      if (.explanatory_addable(out$est$B, D) > 0L)
        out <- relax_explanatory(out, it, "thresholds")
    }
  }
  # Simulated scores have no observed answer-option identities to inherit.
  out$mc <- if (inherit_mc) .explanatory_inherit_mc(
    fit, source, inherit, exclude = unique(fully_relaxed),
    person_rows = person_rows) else NULL
  out
}

#' @export
print.rasch_explanatory <- function(x, ...) {
  cat(sprintf("rasch %s analysis: %d items, %d persons\n",
              x$explanatory_model, ncol(x$X), nrow(x$X)))
  cat("Formula: ", x$explanatory$formula_text, "\n", sep = "")
  cat(sprintf("Conditional calibration: %d explanatory parameter(s), %d fixed departure(s)\n",
              nrow(x$est$coefficients),
              nrow(x$explanatory$relaxations)))
  if (!isTRUE(x$est$converged)) {
    cat("Free calibration comparison: unavailable because the explanatory calibration did not converge\n")
    return(invisible(x))
  }
  tst <- explanatory_test(x)
  if (tst$df > 0L)
    cat(sprintf("Free calibration comparison: adjusted chi-square %.3f on %d df, p = %s\n",
                tst$chisq_kent, tst$df, .fmt_p(tst$p_kent)))
  invisible(x)
}
