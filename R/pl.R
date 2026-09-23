# ---------------------------------------------------------------------------
# Rank analysis: the Plackett-Luce model for rankings of objects
#
# A ranking of n objects is read as a sequence of choices: the object ranked
# first is chosen from all n, the object ranked second from the remaining
# n - 1, and so on (Luce 1959; Plackett 1975). Each stage is a multinomial
# choice with P(j chosen) = exp(beta_j) / sum over the remaining objects, so
# the model reduces to Bradley-Terry-Luce when every ranking holds two
# objects and the Rasch-type conditional structure of btl() carries over:
# the object locations beta are on a logit scale with a sum-zero origin (or
# anchors), estimated by maximum likelihood with Godambe sandwich standard
# errors clustered by judge when judges are recorded, and fit is read from
# the residuals of the stage choices, pooled by object, judge and ranking.
#
# Rankings may be partial: an object listed in a ranking without a rank is
# in every choice set of that ranking but never chosen (a top-k ranking).
# A ranking that lists only some of the objects is a choice among those.
#
# The sequential reading is not reversal-invariant. Fitting the same
# rankings from the bottom up ("worst first") gives a different likelihood
# and, in general, different locations. The reversal check reports both
# fits, the Vuong statistic for the two non-nested readings and the
# agreement of their locations, so a ranking process that did not proceed
# best-first can be recognised.
# ---------------------------------------------------------------------------

# stage representation: one row per (stage, object) cell
.pl_cells <- function(rk_obj, rk_rank, rk_id, K) {
  # rk_obj: object index; rk_rank: rank or NA (unranked, present); rk_id:
  # ranking index. Returns the cells of the informative stages (stage,
  # obj, y per cell) and the ranking of each stage.
  rid <- split(seq_along(rk_id), rk_id)
  cs <- co <- cy <- cr <- vector("list", length(rid))
  st <- 0L
  for (i in seq_along(rid)) {
    r <- rid[[i]]
    ord <- order(rk_rank[r], na.last = TRUE)
    obj <- rk_obj[r][ord]; rnk <- rk_rank[r][ord]
    n_ranked <- sum(!is.na(rnk)); n <- length(obj)
    n_st <- min(n_ranked, n - 1L)
    if (n_st < 1L) next
    sizes <- n - seq_len(n_st) + 1L
    cs[[i]] <- rep(st + seq_len(n_st), sizes)
    co[[i]] <- unlist(lapply(seq_len(n_st), function(s) obj[s:n]))
    cy[[i]] <- unlist(lapply(sizes, function(z) c(1L, rep(0L, z - 1L))))
    cr[[i]] <- rep(rk_id[r][1L], n_st)
    st <- st + n_st
  }
  list(stage = unlist(cs), obj = unlist(co), y = unlist(cy),
       ranking = unlist(cr), n_stages = st)
}

.pl_prob <- function(beta, cells) {
  e <- exp(beta[cells$obj])
  den <- rowsum(e, cells$stage)[, 1L]
  e / den[cells$stage]
}

.pl_ll <- function(beta, cells) {
  p <- .pl_prob(beta, cells)
  sum(log(p[cells$y == 1L]))
}

.pl_grad <- function(beta, cells, K) {
  p <- .pl_prob(beta, cells)
  g <- rowsum(cells$y - p, cells$obj)
  out <- numeric(K); out[as.integer(rownames(g))] <- g[, 1L]
  out
}

.pl_hess <- function(beta, cells, K) {
  p <- .pl_prob(beta, cells)
  P <- matrix(0, cells$n_stages, K)
  P[cbind(cells$stage, cells$obj)] <- p
  crossprod(P) - diag(colSums(P), K)
}

# score contributions by cluster (a nc x K matrix)
.pl_scores <- function(beta, cells, K, cluster) {
  p <- .pl_prob(beta, cells)
  nc <- max(cluster)
  G <- matrix(0, nc, K)
  s <- rowsum(cells$y - p, (cluster[cells$stage] - 1L) * K + cells$obj)
  idx <- as.integer(rownames(s))
  G[cbind((idx - 1L) %/% K + 1L, (idx - 1L) %% K + 1L)] <- s[, 1L]
  G
}

# Newton-Raphson on the free parameters theta, beta = beta0 + B theta
.pl_newton <- function(cells, K, B, beta0, start, maxit, tol) {
  theta <- start
  ll <- .pl_ll(beta0 + drop(B %*% theta), cells)
  it <- 0L; converged <- FALSE
  n_units <- cells$n_stages
  while (it < maxit) {
    it <- it + 1L
    beta <- beta0 + drop(B %*% theta)
    g <- crossprod(B, .pl_grad(beta, cells, K))
    H <- crossprod(B, .pl_hess(beta, cells, K) %*% B)
    if (max(abs(g)) < 1e-6 * n_units) { converged <- TRUE; break }
    step <- tryCatch(-solve(H, g), error = function(e) NULL)
    if (is.null(step)) step <- -solve(H - diag(1e-6, nrow(H)), g)
    step <- drop(step)
    lam <- 1
    repeat {
      cand <- theta + lam * step
      ll_new <- .pl_ll(beta0 + drop(B %*% cand), cells)
      if (is.finite(ll_new) && ll_new >= ll - 1e-10) break
      lam <- lam / 2
      if (lam < 1e-8) break
    }
    theta <- cand; ll <- ll_new
    if (max(abs(lam * step)) < tol) {
      beta <- beta0 + drop(B %*% theta)
      g <- crossprod(B, .pl_grad(beta, cells, K))
      converged <- max(abs(g)) < 1e-6 * n_units
      break
    }
  }
  beta <- beta0 + drop(B %*% theta)
  list(theta = theta, beta = beta, ll = ll, iterations = it,
       converged = converged)
}

# reachability closure of a directed adjacency matrix (i beats j)
.pl_reach <- function(adj) {
  R <- adj | diag(TRUE, nrow(adj))
  repeat {
    R2 <- (R %*% R) > 0
    if (identical(R2, R)) return(R)
    R <- R2
  }
}

#' Rank analysis with the Plackett-Luce model
#'
#' Calibrates objects from rankings. Each ranking is read as a sequence of
#' choices, the object ranked first chosen from all the objects in the
#' ranking, the object ranked second from those remaining, and so on, with
#' \eqn{P(j \mathrm{\ chosen}) = \exp(\beta_j) / \sum_k \exp(\beta_k)} over
#' the objects still to be placed (Luce 1959; Plackett 1975). A ranking of
#' two objects is a paired comparison, and the model then equals
#' \code{\link{btl}}: the locations are on the same logit scale and, when
#' every ranking is a pair, the two functions agree.
#'
#' The locations are maximum likelihood estimates with a sum-zero origin,
#' or at the values in \code{anchors}. Standard errors are Godambe sandwich
#' errors clustered by judge when a judge column is given and by ranking
#' otherwise, subject to the same conditions as \code{btl()}: with fewer
#' than ten judges, fewer judges than parameters, or fewer than eight
#' effective judges, the clustered covariance is not calibrated and is
#' withheld with a note. \code{se = "model"} instead reports the
#' information-based errors of the Plackett-Luce likelihood, which are
#' valid when the model holds and the rankings are independent, and are
#' available whatever the design.
#'
#' Fit is read from the residuals of the stage choices, \eqn{y - p} for each
#' object in each choice set, pooled by object and by judge into the infit,
#' outfit and standardised fit residuals of \code{btl()}. Each ranking also
#' receives a surprise statistic: its log-likelihood standardised against
#' the mean and variance the fitted model implies for it, so a ranking the
#' consensus makes very unlikely stands out.
#'
#' An object that is chosen at every stage it appears in, or never chosen
#' at any stage with an alternative, has no finite location (the undefeated
#' or winless case of a paired comparison). It is set aside and reported
#' with \code{extreme = TRUE} at an extrapolated location, its score moved
#' half a choice inside the boundary against the calibrated objects. A
#' design in which a group of objects is never ranked below an object
#' outside the group cannot place that group and is refused.
#'
#' Rankings may be partial in two ways. A ranking that lists only some of
#' the objects is a choice among those. An object listed with a missing
#' \code{rank} is present but unranked: it is in every choice set of that
#' ranking and never chosen, which is the top-\eqn{k} design. Tied ranks
#' within a ranking are not modelled; a ranking with a tie is dropped
#' (\code{ties = "drop"}) or refused.
#'
#' \strong{Reversal check.} The sequential reading is not symmetric: read
#' from the bottom up, with the worst object chosen first, the same rankings
#' have a different likelihood and in general different locations. If the
#' judges worked best-first the forward reading should fit at least as well
#' as the reversed one. The check fits both readings to the rankings that
#' are complete and hold at least three objects, reports the Vuong (1989)
#' statistic for the two non-nested models (positive favours best-first)
#' and the correlation and largest difference between the two sets of
#' locations. It is a check on the ranking process, not on the objects.
#'
#' @param data A data frame in long format: one row per object per ranking.
#' @param ranking,object,rank Names of the columns holding the ranking
#'   identifier, the object and the rank (1 = highest). A missing rank
#'   marks an object present in the ranking but unranked.
#' @param judge Optional name of a judge column. Rankings are then clustered
#'   by judge for the sandwich covariance and a judge fit table is reported.
#' @param anchors Optional named numeric vector of fixed object locations.
#' @param ties How to treat a ranking with tied ranks: \code{"drop"} the
#'   ranking with a note, or \code{"error"}.
#' @param se \code{"sandwich"} (clustered Godambe errors, withheld when the
#'   cluster design does not support them) or \code{"model"} (inverse
#'   observed information of the Plackett-Luce likelihood).
#' @param maxit,tol Newton-Raphson controls.
#' @return A \code{"rasch_pl"} object with \code{objects} (location, se,
#'   rankings, stages, chosen, infit, outfit, fit residual, extreme flag),
#'   \code{judges} (when a judge column is given), \code{rankings} (one row
#'   per ranking with its log-likelihood and surprise \code{z}),
#'   \code{reversal} (the reversal check, or \code{NULL} when no ranking is
#'   complete with three or more objects), \code{osi} (object separation),
#'   \code{loglik}, \code{cov_beta}, \code{converged}, \code{iterations},
#'   \code{n_rankings}, \code{n_stages}, \code{se_type}, \code{anchors},
#'   \code{notes} and the \code{call}.
#' @references
#' Luce, R. D. (1959). Individual Choice Behavior. Wiley.
#'
#' Plackett, R. L. (1975). The analysis of permutations. Applied
#' Statistics, 24, 193--202.
#'
#' Hunter, D. R. (2004). MM algorithms for generalized Bradley-Terry
#' models. Annals of Statistics, 32, 384--406.
#'
#' Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
#' non-nested hypotheses. Econometrica, 57, 307--333.
#' @seealso \code{\link{btl}} for paired comparisons, \code{\link{rasch_cj}}
#'   to combine rankings with comparisons and item responses,
#'   \code{\link{plot_pl}}.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1.5, B = -0.5, C = 0, D = 0.5, E = 1.5)
#' rk <- do.call(rbind, lapply(1:60, function(r) {
#'   rem <- names(beta); ord <- character(0)
#'   while (length(rem) > 1) {
#'     pick <- sample(rem, 1, prob = exp(beta[rem]))
#'     ord <- c(ord, pick); rem <- setdiff(rem, pick)
#'   }
#'   data.frame(ranking = r, object = c(ord, rem), rank = 1:5)
#' }))
#' fit <- pl(rk)
#' fit
#' fit$reversal$z
#' @export
pl <- function(data, ranking = "ranking", object = "object", rank = "rank",
               judge = NULL, anchors = NULL, ties = c("drop", "error"),
               se = c("sandwich", "model"), maxit = 100, tol = 1e-8) {
  .check_column_names(data)
  .check_controls(maxit, tol)
  ties <- match.arg(ties); se <- match.arg(se)
  data <- as.data.frame(data)
  for (nm in c("ranking", "object", "rank", "judge")) {
    v <- get(nm, inherits = FALSE)
    if (!is.null(v)) .check_reshape_column(data, v, nm)
  }
  roles <- c(ranking, object, rank, judge)
  if (anyDuplicated(roles))
    stop("ranking role columns must be distinct; repeated: ",
         paste(unique(roles[duplicated(roles)]), collapse = ", "))
  rid <- .role_text_values(data[[ranking]])
  ob <- .role_text_values(data[[object]])
  if (any(!is.na(ob) & !nzchar(ob)))
    stop("blank object identifier(s) in ", object,
         "; a whitespace-only name is not an object")
  rkx <- data[[rank]]
  if (is.complex(rkx)) stop("`", rank, "` is complex; it must be numeric")
  rk <- suppressWarnings(as.numeric(if (is.factor(rkx))
    as.character(rkx) else rkx))
  if (any(is.na(rk) & !is.na(rkx)))
    stop("`", rank, "` has non-numeric value(s); it must be numeric")
  if (any(!is.finite(rk[!is.na(rk)])))
    stop("`", rank, "` must hold finite rank values")
  jd <- if (is.null(judge)) NULL else .role_text_values(data[[judge]])
  if (!is.null(jd) && any(!is.na(jd) & !nzchar(jd)))
    stop("blank judge identifier(s) in ", judge,
         "; a whitespace-only name is not a judge")
  notes <- character(0)

  keep <- !is.na(rid) & !is.na(ob)
  if (!is.null(jd)) keep <- keep & !is.na(jd)
  if (any(!keep)) {
    notes <- c(notes, sprintf("%d row(s) dropped (missing ranking, object or judge)",
                              sum(!keep)))
    rid <- rid[keep]; ob <- ob[keep]; rk <- rk[keep]
    if (!is.null(jd)) jd <- jd[keep]
  }
  if (!length(rid)) stop("no usable rankings")
  # a judge is a property of the ranking, not of its rows
  if (!is.null(jd)) {
    nj <- tapply(jd, rid, function(z) length(unique(z)))
    if (any(nj > 1))
      stop("ranking(s) attributed to more than one judge: ",
           paste(names(nj)[nj > 1], collapse = ", "))
  }
  dup <- tapply(ob, rid, anyDuplicated)
  if (any(dup > 0))
    stop("ranking(s) list an object more than once: ",
         paste(names(dup)[dup > 0], collapse = ", "))
  tied <- tapply(rk, rid, function(z) anyDuplicated(z[!is.na(z)]) > 0)
  if (any(tied)) {
    if (ties == "error")
      stop(sum(tied), " ranking(s) with tied ranks; set ties = 'drop'")
    notes <- c(notes, sprintf("%d ranking(s) with tied ranks dropped", sum(tied)))
    sel <- !(rid %in% names(tied)[tied])
    rid <- rid[sel]; ob <- ob[sel]; rk <- rk[sel]
    if (!is.null(jd)) jd <- jd[sel]
  }
  # a ranking needs at least two objects and one rank to say anything
  n_obj <- tapply(ob, rid, length)
  n_rk <- tapply(rk, rid, function(z) sum(!is.na(z)))
  small <- names(n_obj)[n_obj < 2 | n_rk < 1]
  if (length(small)) {
    notes <- c(notes, sprintf(
      "%d ranking(s) with fewer than two objects or no ranked object dropped",
      length(small)))
    sel <- !(rid %in% small)
    rid <- rid[sel]; ob <- ob[sel]; rk <- rk[sel]
    if (!is.null(jd)) jd <- jd[sel]
  }
  if (!length(rid)) stop("no usable rankings")

  objs <- sort(unique(ob))
  if (!is.null(anchors)) {
    if (!is.numeric(anchors) || is.complex(anchors) || !is.null(dim(anchors)) ||
        !is.null(oldClass(anchors)) || !length(anchors) ||
        is.null(names(anchors)) || anyNA(names(anchors)) ||
        any(!nzchar(trimws(names(anchors)))))
      stop("`anchors` must be a non-empty named numeric vector ",
           "(names = object names)")
    names(anchors) <- .role_text_values(names(anchors))
    if (anyDuplicated(names(anchors)))
      stop("duplicate anchor(s): ",
           paste(unique(names(anchors)[duplicated(names(anchors))]),
                 collapse = ", "))
    if (any(!is.finite(anchors)))
      stop("anchor value(s) must be finite: ",
           paste(names(anchors)[!is.finite(anchors)], collapse = ", "))
    bad <- setdiff(names(anchors), objs)
    if (length(bad))
      stop("`anchors` name(s) do not match any object in the data: ",
           paste(bad, collapse = ", "))
  }
  K <- length(objs)
  ru <- unique(rid); R <- length(ru)
  ri <- match(rid, ru); oi <- match(ob, objs)
  judge_of <- if (is.null(jd)) NULL else jd[match(ru, rid)]

  cells0 <- .pl_cells(oi, rk, ri, K)
  if (!cells0$n_stages) stop("no informative choice in the rankings")

  # extreme objects: chosen whenever they could be, or never chosen. They
  # leave the estimation and are placed afterwards.
  in_set <- rep(TRUE, K)
  cells <- cells0
  extreme <- character(0)
  repeat {
    present <- tabulate(cells$obj, K)
    chosen <- tabulate(cells$obj[cells$y == 1L], K)
    ext <- which(in_set & present > 0 & (chosen == 0 | chosen == present))
    if (!length(ext)) break
    extreme <- c(extreme, objs[ext])
    in_set[ext] <- FALSE
    # drop those objects from every set; a stage whose choice was an
    # extreme object disappears, a set that shrinks below two is over
    drop_stage <- unique(cells$stage[cells$y == 1L & cells$obj %in% ext])
    sel <- !(cells$obj %in% ext) & !(cells$stage %in% drop_stage)
    size <- tabulate(cells$stage[sel], cells$n_stages)
    sel <- sel & size[cells$stage] >= 2L
    kept <- unique(cells$stage[sel])
    cells <- list(stage = match(cells$stage[sel], kept), obj = cells$obj[sel],
                  y = cells$y[sel], ranking = cells$ranking[kept],
                  n_stages = length(kept))
    if (!cells$n_stages) break
  }
  # objects that are listed but never in an informative stage
  never <- objs[in_set & tabulate(cells$obj, K) == 0]
  if (length(never))
    stop("object(s) never in a choice set with an alternative: ",
         paste(never, collapse = ", "))
  fit_idx <- which(in_set)
  if (length(fit_idx) < 2L)
    stop("fewer than two objects remain after setting aside extreme ",
         "objects (", paste(extreme, collapse = ", "), ")")
  if (length(extreme))
    notes <- c(notes, sprintf(
      "%d extreme object(s) set aside (always or never chosen) and reported at extrapolated locations: %s",
      length(extreme), paste(extreme, collapse = ", ")))

  # the design must connect the fitted objects: i beats j when i is chosen
  # from a set holding j; every object must reach every other
  adj <- matrix(FALSE, K, K)
  ch <- cells$obj[cells$y == 1L][cells$stage]
  adj[cbind(ch[cells$y == 0L], cells$obj[cells$y == 0L])] <- TRUE
  adj <- adj[fit_idx, fit_idx, drop = FALSE]
  reach <- .pl_reach(adj)
  if (!all(reach)) {
    scc <- reach & t(reach)
    comp <- match(apply(scc, 1, function(z) paste(which(z), collapse = ",")),
                  unique(apply(scc, 1, function(z) paste(which(z), collapse = ","))))
    # a source group is beaten by nothing outside it
    src <- vapply(unique(comp), function(cc)
      !any(adj[comp != cc, comp == cc]), NA)
    top <- objs[fit_idx][comp == unique(comp)[which(src)[1]]]
    stop("the rankings do not place every object on one scale: ",
         paste(top, collapse = ", "),
         if (length(top) > 1) " are" else " is",
         " never ranked below an object outside that group; add rankings ",
         "that mix the groups, or anchor")
  }

  # map cells to the fitted objects
  fk <- length(fit_idx)
  cells$obj <- match(cells$obj, fit_idx)
  anch_idx <- if (is.null(anchors)) integer(0) else
    match(names(anchors), objs[fit_idx])
  anch_idx <- anch_idx[!is.na(anch_idx)]
  beta0 <- numeric(fk)
  if (length(anch_idx)) {
    beta0[anch_idx] <- anchors[objs[fit_idx][anch_idx]]
    free <- setdiff(seq_len(fk), anch_idx)
    if (!length(free)) stop("every fitted object is anchored; nothing to estimate")
    B <- diag(fk)[, free, drop = FALSE]
  } else {
    B <- rbind(diag(fk - 1L), -1)
  }
  np <- ncol(B)
  fit <- .pl_newton(cells, fk, B, beta0, rep(0, np), maxit, tol)
  if (!fit$converged)
    warning("pl estimation did NOT converge in ", fit$iterations,
            " iterations; estimates and standard errors are unreliable",
            call. = FALSE)
  beta <- fit$beta
  H <- crossprod(B, .pl_hess(beta, cells, fk) %*% B)
  bread <- tryCatch(solve(-H), error = function(e) NULL)

  # covariance: clustered sandwich (judge, else ranking), or model-based
  cluster <- if (is.null(jd)) cells$ranking else judge_of[cells$ranking]
  cl_ids <- match(cluster, unique(cluster))
  nc <- max(cl_ids)
  se_available <- isTRUE(fit$converged) && !is.null(bread)
  se_type <- se
  if (se_available && se == "sandwich") {
    G <- .pl_scores(beta, cells, fk, cl_ids) %*% B
    M <- crossprod(G)
    sv <- svd(M, nu = 0, nv = 0)$d
    score_rank <- if (!length(sv) || max(sv) <= 0) 0L else
      sum(sv > max(sv) * sqrt(.Machine$double.eps))
    nc_eff <- if (is.null(jd)) Inf else {
      shr <- tabulate(cl_ids, nc) / cells$n_stages
      1 / sum(shr^2)
    }
    ok <- score_rank >= np && (is.null(jd) ||
      (nc >= 10L && nc > np && nc_eff >= 8))
    if (ok) {
      cov_th <- bread %*% M %*% bread
    } else {
      se_available <- FALSE
      notes <- c(notes, if (is.null(jd)) sprintf(
        "the ranking-clustered score covariance has rank %d for %d parameters; sandwich standard errors are withheld (se = 'model' gives the information-based errors)",
        score_rank, np) else sprintf(
        "%d judge clusters (%.1f effective) for %d parameters do not support clustered sandwich errors; they are withheld (se = 'model' gives the information-based errors)",
        nc, nc_eff, np))
    }
  } else if (se_available) {
    cov_th <- bread
  }
  cov_beta <- matrix(NA_real_, K, K, dimnames = list(objs, objs))
  se_b <- rep(NA_real_, K)
  loc <- rep(NA_real_, K)
  loc[fit_idx] <- beta
  if (se_available) {
    cb <- B %*% cov_th %*% t(B)
    cov_beta[fit_idx, fit_idx] <- cb
    se_b[fit_idx] <- sqrt(pmax(diag(cb), 0))
  }
  if (length(anch_idx)) {
    cov_beta[fit_idx[anch_idx], ] <- 0; cov_beta[, fit_idx[anch_idx]] <- 0
    se_b[fit_idx[anch_idx]] <- 0
  }

  # fit statistics from the stage residuals
  p <- .pl_prob(beta, cells)
  y <- cells$y
  v <- p * (1 - p)
  z <- (y - p) / sqrt(v)
  c4v <- (1 - 4 * v) / v
  n_cells <- length(y)
  f_cell <- (n_cells - np) / n_cells
  pool <- function(sel) {
    n <- sum(sel)
    if (n < 3)
      return(list(infit_ms = NA_real_, outfit_ms = NA_real_,
                  fit_resid = NA_real_, df = NA_real_, n = n))
    y2 <- sum(z[sel]^2); f <- f_cell * n
    wv <- sum(v[sel])
    infit <- if (wv > 1e-12) sum((y[sel] - p[sel])^2) / (f_cell * wv) else NA_real_
    vv <- sum(c4v[sel])
    fr <- if (vv > 1e-8 && y2 > 0) f * (log(y2) - log(f)) / sqrt(vv) else NA_real_
    list(infit_ms = infit, outfit_ms = y2 / f, fit_resid = fr, df = f, n = n)
  }
  ofit <- lapply(seq_len(fk), function(k) pool(cells$obj == k))
  obj_rankings <- vapply(seq_len(K), function(k)
    length(unique(ri[oi == k])), 0L)
  objects <- data.frame(object = objs, location = loc, se = se_b,
                        rankings = obj_rankings,
                        stages = tabulate(cells0$obj, K),
                        chosen = tabulate(cells0$obj[cells0$y == 1L], K),
                        infit_ms = NA_real_, outfit_ms = NA_real_,
                        fit_resid = NA_real_, df_fit = NA_real_,
                        extreme = !in_set, stringsAsFactors = FALSE)
  objects$infit_ms[fit_idx] <- vapply(ofit, `[[`, 0, "infit_ms")
  objects$outfit_ms[fit_idx] <- vapply(ofit, `[[`, 0, "outfit_ms")
  objects$fit_resid[fit_idx] <- vapply(ofit, `[[`, 0, "fit_resid")
  objects$df_fit[fit_idx] <- vapply(ofit, `[[`, 0, "df")

  # extrapolated locations for the extreme objects: expected number of
  # choices equal to the score moved half a choice inside its range, over
  # the original stages restricted to the calibrated objects
  for (e in which(!in_set)) {
    st_e <- unique(cells0$stage[cells0$obj == e])
    sets <- lapply(st_e, function(s) {
      o <- cells0$obj[cells0$stage == s]
      o[o %in% fit_idx]
    })
    st_ok <- lengths(sets) >= 1L
    if (!any(st_ok)) next
    sets <- sets[st_ok]
    n_e <- length(sets)
    T_e <- sum(cells0$y[cells0$obj == e][st_ok])
    Tstar <- if (T_e <= 0) 0.5 else n_e - 0.5
    g <- function(th) sum(vapply(sets, function(o) {
      d <- c(th, loc[o]); exp(th) / sum(exp(d))
    }, 0)) - Tstar
    lim <- range(beta) + c(-12, 12)
    root <- tryCatch(stats::uniroot(g, lim, tol = 1e-8)$root,
                     error = function(err) NA_real_)
    objects$location[e] <- root
  }

  # per-ranking surprise
  lp <- log(p)
  E_s <- rowsum(p * lp, cells$stage)[, 1L]
  V_s <- rowsum(p * lp^2, cells$stage)[, 1L] - E_s^2
  ll_s <- lp[y == 1L]
  rk_of_stage <- cells$ranking
  ll_r <- rowsum(ll_s, rk_of_stage)[, 1L]
  E_r <- rowsum(E_s, rk_of_stage)[, 1L]
  V_r <- rowsum(V_s, rk_of_stage)[, 1L]
  rr <- as.integer(names(ll_r))
  rankings <- data.frame(ranking = ru, judge = if (is.null(jd)) NA_character_ else
                           judge_of,
                         n_objects = as.integer(tabulate(ri, R)),
                         n_ranked = as.integer(tapply(!is.na(rk), ri, sum)),
                         stages = 0L, loglik = NA_real_, surprise_z = NA_real_,
                         stringsAsFactors = FALSE)
  rankings$stages[rr] <- as.integer(tabulate(rk_of_stage, R))[rr]
  rankings$loglik[rr] <- ll_r
  rankings$surprise_z[rr] <- ifelse(V_r > 1e-10, (E_r - ll_r) / sqrt(V_r), NA_real_)
  rownames(rankings) <- NULL

  judges <- NULL
  if (!is.null(jd)) {
    ju <- unique(judge_of)
    cj <- judge_of[cells$ranking][cells$stage]
    jfit <- lapply(ju, function(j) pool(cj == j))
    judges <- data.frame(judge = ju,
                         rankings = as.integer(table(factor(judge_of, ju))),
                         stages = vapply(ju, function(j)
                           sum(rk_of_stage %in% which(judge_of == j)), 0L),
                         infit_ms = vapply(jfit, `[[`, 0, "infit_ms"),
                         outfit_ms = vapply(jfit, `[[`, 0, "outfit_ms"),
                         fit_resid = vapply(jfit, `[[`, 0, "fit_resid"),
                         df_fit = vapply(jfit, `[[`, 0, "df"),
                         stringsAsFactors = FALSE)
    js <- tapply(rankings$surprise_z, factor(rankings$judge, ju),
                 function(z) mean(z, na.rm = TRUE))
    judges$mean_surprise_z <- as.numeric(js)
    judges <- judges[order(judges$judge), ]
    rownames(judges) <- NULL
  }

  osi <- if (se_available) .psi(objects$location, objects$se, !objects$extreme)
         else list(PSI = NA_real_, separation = NA_real_, strata = NA_real_,
                   var_theta = NA_real_, mean_error_var = NA_real_, n = 0L)

  # reversal check on the complete rankings of three or more objects
  reversal <- NULL
  complete <- which(tabulate(ri, R) >= 3 &
                    tapply(!is.na(rk), ri, sum) == tabulate(ri, R))
  if (length(complete) >= 2L && isTRUE(fit$converged)) {
    selr <- ri %in% complete & oi %in% fit_idx
    if (length(unique(oi[selr])) == fk) {
      fwd <- .pl_cells(match(oi[selr], fit_idx), rk[selr], ri[selr], fk)
      rev <- .pl_cells(match(oi[selr], fit_idx), -rk[selr], ri[selr], fk)
      ok_rev <- fwd$n_stages > 0 && rev$n_stages > 0 &&
        all(tabulate(fwd$obj, fk) > 0) && all(tabulate(rev$obj, fk) > 0)
      if (ok_rev) {
        f_fwd <- if (length(complete) == R) fit else
          .pl_newton(fwd, fk, B, beta0, rep(0, np), maxit, tol)
        f_rev <- .pl_newton(rev, fk, B, beta0, rep(0, np), maxit, tol)
        if (isTRUE(f_fwd$converged) && isTRUE(f_rev$converged)) {
          p_f <- .pl_prob(f_fwd$beta, fwd); p_r <- .pl_prob(f_rev$beta, rev)
          rf <- fwd$ranking; rv <- rev$ranking
          llf <- rowsum(log(p_f[fwd$y == 1L]), rf)[, 1L]
          llr <- rowsum(log(p_r[rev$y == 1L]), rv)[, 1L]
          d <- llf - llr[names(llf)]
          n_d <- length(d)
          zv <- if (n_d > 1 && stats::sd(d) > 0)
            mean(d) / (stats::sd(d) / sqrt(n_d)) else NA_real_
          loc_r <- -f_rev$beta
          loc_r <- loc_r - mean(loc_r) + mean(f_fwd$beta)
          reversal <- list(
            n_rankings = n_d,
            loglik_forward = sum(llf), loglik_reversed = sum(llr),
            z = zv, p = if (is.finite(zv)) 2 * stats::pnorm(-abs(zv)) else NA_real_,
            correlation = if (fk > 2) stats::cor(f_fwd$beta, loc_r) else NA_real_,
            max_abs_difference = max(abs(f_fwd$beta - loc_r)),
            objects = data.frame(object = objs[fit_idx],
                                 forward = f_fwd$beta, reversed = loc_r,
                                 difference = f_fwd$beta - loc_r,
                                 stringsAsFactors = FALSE))
        }
      }
    }
  }

  out <- list(objects = objects, judges = judges, rankings = rankings,
              reversal = reversal, osi = osi, loglik = fit$ll,
              cov_beta = cov_beta, converged = fit$converged,
              iterations = fit$iterations, n_rankings = R,
              n_stages = cells0$n_stages,
              size = range(tabulate(ri, R)),
              clustered = !is.null(jd), se_type = se_type,
              se_available = se_available, n_clusters = nc,
              anchors = anchors, notes = notes, call = match.call())
  out <- .tag_tables(out)
  class(out) <- "rasch_pl"
  out
}

#' @export
print.rasch_pl <- function(x, ...) {
  cat(sprintf("Plackett-Luce rank analysis: %d objects, %d rankings of %s objects%s\n",
              nrow(x$objects), x$n_rankings,
              if (x$size[1] == x$size[2]) x$size[1] else
                paste(x$size, collapse = " to "),
              if (!is.null(x$judges)) sprintf(", %d judges", nrow(x$judges)) else ""))
  se_lab <- if (!x$se_available) "standard errors withheld" else
    if (x$se_type == "model") "information-based SEs" else
    paste0("sandwich SEs clustered by ",
           if (x$clustered) "judge" else "ranking")
  cat(sprintf("Maximum likelihood: %s in %d iterations; %s\n",
              if (x$converged) "converged" else "NOT converged",
              x$iterations, se_lab))
  cat(sprintf("Object separation index %.3f\n", x$osi$PSI))
  if (!is.null(x$anchors))
    cat(sprintf("Anchored at %d object(s) (se = 0): %s\n",
                length(x$anchors), paste(names(x$anchors), collapse = ", ")))
  if (!is.null(x$reversal)) {
    r <- x$reversal
    cat(sprintf(paste0("Reversal check on %d complete rankings: best-first ",
                       "log-likelihood %.2f, worst-first %.2f; Vuong z = %.2f, ",
                       "p = %s; location correlation %.3f\n"),
                r$n_rankings, r$loglik_forward, r$loglik_reversed, r$z,
                .fmt_p(r$p), r$correlation))
  }
  print(.fmt_df(x$objects[, c("object", "location", "se", "rankings",
                              "chosen", "fit_resid", "extreme")]),
        row.names = FALSE)
  if (!is.null(x$judges)) {
    mis <- x$judges[!is.na(x$judges$fit_resid) & abs(x$judges$fit_resid) > 2.5, ]
    cat(sprintf("Judges beyond |fit residual| 2.5: %d%s\n", nrow(mis),
                if (nrow(mis)) paste0(" (", paste(mis$judge, collapse = ", "), ")")
                else ""))
  }
  sur <- x$rankings[!is.na(x$rankings$surprise_z) & x$rankings$surprise_z > 2.5, ]
  cat(sprintf("Rankings beyond surprise z 2.5: %d%s\n", nrow(sur),
              if (nrow(sur) && nrow(sur) <= 10)
                paste0(" (", paste(sur$ranking, collapse = ", "), ")") else ""))
  if (length(x$notes)) cat(sprintf("Notes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

#' Plot Plackett-Luce object locations
#'
#' Caterpillar plot of the object locations with 95 per cent error bars,
#' the objects beyond the fit-residual band marked and extreme objects shown
#' without an interval. The interval uses a t reference with judges minus
#' one degrees of freedom for judge-clustered errors and the normal
#' reference otherwise.
#'
#' @param fit An object from \code{\link{pl}}.
#' @param band Absolute fit-residual value beyond which an object is
#'   highlighted.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' rk <- do.call(rbind, lapply(1:40, function(r) {
#'   rem <- names(beta); ord <- character(0)
#'   while (length(rem) > 1) {
#'     pick <- sample(rem, 1, prob = exp(beta[rem]))
#'     ord <- c(ord, pick); rem <- setdiff(rem, pick)
#'   }
#'   data.frame(ranking = r, object = c(ord, rem), rank = 1:4)
#' }))
#' plot_pl(pl(rk))
#' @export
plot_pl <- function(fit, band = 2.5) {
  if (!inherits(fit, "rasch_pl")) stop("not a rank analysis (pl) fit", call. = FALSE)
  if (!isTRUE(fit$converged))
    stop("the rank calibration did not converge; fitted displays are unavailable",
         call. = FALSE)
  .check_band(band)
  d <- fit$objects[order(fit$objects$location), ]
  k <- nrow(d)
  ref_df <- if (fit$clustered && fit$se_type == "sandwich")
    max(fit$n_clusters - 1L, 1L) else Inf
  critical <- stats::qt(0.975, df = ref_df)
  has_se <- is.finite(d$se) & d$se > 0
  lower <- upper <- rep(NA_real_, k)
  lower[has_se] <- d$location[has_se] - critical * d$se[has_se]
  upper[has_se] <- d$location[has_se] + critical * d$se[has_se]
  xerr <- c(lower[has_se], upper[has_se], d$location)
  xlim <- range(xerr[is.finite(xerr)])
  op <- .rr_canvas(xlim + c(-0.15, 0.15) * diff(xlim), c(0.5, k + 0.5),
                   "Location (logits)", "", grid_y = FALSE, grid_x = TRUE,
                   yaxis = FALSE)
  on.exit(par(op))
  mis <- !is.na(d$fit_resid) & abs(d$fit_resid) > band
  segments(lower[has_se], which(has_se), upper[has_se], which(has_se),
           col = ifelse(mis[has_se], .rr$red, .rr$soft), lwd = 2.2)
  points(d$location, seq_len(k), pch = ifelse(d$extreme, 24, 21), cex = 1.6,
         lwd = 1.2, bg = ifelse(mis, .rr$red, .rr$blue), col = "white")
  text(d$location, seq_len(k), d$object, pos = 3, offset = 0.55, cex = 0.8,
       col = .rr$ink)
  if (any(mis))
    .rr_legend("bottomright", sprintf("|fit residual| > %.1f", band),
               pch = 21, pt.bg = .rr$red, col = "white", pt.cex = 1.4)
  invisible(NULL)
}
