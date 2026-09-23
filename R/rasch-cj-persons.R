# rasch :: persons from item responses and judgements of their work
# ===========================================================================
# The person mode of rasch_cj(). The items are calibrated already (anchored
# thresholds), and each person is located by two kinds of evidence: their
# responses to the items, and judges' comparisons or rankings of their
# work. The persons are the objects the judges compare, so the blocks are
#
#   responses    the partial credit likelihood of each person's responses
#                given the anchored thresholds, unit 1;
#   comparisons  Bradley-Terry in alpha (theta_a - theta_b);
#   rankings     Plackett-Luce in kappa theta.
#
# The response block fixes the origin and the unit, so nothing is
# constrained: theta is one free location per person. Every block is a
# likelihood in a fixed number of parameters (one per person plus a unit per
# judgement frame), so this is not joint maximum likelihood over persons and
# items: the items are held, and each person's information grows with their
# own items and judgements.
#
# A person has a finite location when the evidence points both ways: a
# response score strictly inside its range, or at least one win and one
# loss among the persons still being estimated. A person with evidence in
# one direction only is set aside as extreme and reported at the Warm
# estimate from their responses, as rasch() reports extreme persons. A group
# of judged persons none of whom has responses is not on the test scale and
# is refused.
# ===========================================================================

# Log-likelihood, gradient and curvature of the response block, one value
# per person, from the partial credit model with fixed thresholds.
.cj_person_parts <- function(theta, X, tau_list) {
  N <- length(theta); ll <- g <- h <- numeric(N)
  for (i in seq_along(tau_list)) {
    x <- X[, i]; ok <- which(!is.na(x))
    if (!length(ok)) next
    tau <- tau_list[[i]]; m <- length(tau); cs <- c(0, cumsum(tau)); cat_ <- 0:m
    lp <- outer(theta[ok], cat_) - rep(cs, each = length(ok))
    mx <- lp[cbind(seq_along(ok), max.col(lp, ties.method = "first"))]
    e <- exp(lp - mx); s <- rowSums(e); P <- e / s
    E <- as.vector(P %*% cat_); V <- as.vector(P %*% cat_^2) - E^2
    xo <- x[ok]
    ll[ok] <- ll[ok] + xo * theta[ok] - cs[xo + 1L] - (mx + log(s))
    g[ok] <- g[ok] + xo - E
    h[ok] <- h[ok] - V
  }
  list(ll = ll, g = g, h = h)
}

# Maximum likelihood location of each person from responses alone, for the
# persons whose score is inside its range; Inf or -Inf at an extreme score.
.cj_person_ml <- function(X, tau_list, finite, maxit = 100L, tol = 1e-8) {
  N <- nrow(X); theta <- rep(0, N)
  active <- which(finite)
  for (it in seq_len(maxit)) {
    if (!length(active)) break
    p <- .cj_person_parts(theta[active], X[active, , drop = FALSE], tau_list)
    step <- -p$g / p$h
    step[!is.finite(step)] <- 0
    step <- pmax(pmin(step, 1), -1)
    theta[active] <- theta[active] + step
    active <- active[abs(step) >= tol]
  }
  p <- .cj_person_parts(theta, X, tau_list)
  se <- sqrt(pmax(-1 / p$h, 0))
  raw <- rowSums(X, na.rm = TRUE)
  max_raw <- as.vector((!is.na(X)) %*% vapply(tau_list, length, 1L))
  theta[!finite] <- ifelse(max_raw[!finite] == 0L, NA_real_,
                           ifelse(raw[!finite] >= max_raw[!finite], Inf, -Inf))
  se[!finite] <- NA_real_
  list(theta = theta, se = se, ll = p$ll)
}

# Persons with evidence in both directions among the judgements of those
# still active, given prior evidence up0 (some score above the minimum) and
# down0 (some score below the maximum). Persons lacking either are removed
# with their comparisons and rankings, and the rule is applied again until
# nothing changes.
.cj_evidence_prune <- function(active, win, lose, rk, up0, down0) {
  repeat {
    keep_c <- active[win] & active[lose]
    rk_k <- lapply(rk, function(v) v[active[v]]); rk_k <- rk_k[lengths(rk_k) >= 2L]
    up <- up0; down <- down0
    up[win[keep_c]] <- TRUE; down[lose[keep_c]] <- TRUE
    for (v in rk_k) { up[v[-length(v)]] <- TRUE; down[v[-1L]] <- TRUE }
    ext <- which(active & !(up & down))
    if (!length(ext)) break
    active[ext] <- FALSE
  }
  list(active = active, win = win[keep_c], lose = lose[keep_c], rk = rk_k,
       up = up, down = down, n_dropped = sum(!keep_c))
}

# Resolve the anchors argument to a threshold list named by item.
.cj_anchor_list <- function(anchors) {
  if (inherits(anchors, "rasch_cj")) anchors <- anchors$anchors
  else if (inherits(anchors, "rasch")) {
    if (is.null(anchors$tau_list))
      stop("`anchors` is a fit without item thresholds", call. = FALSE)
    return(anchors$tau_list)
  }
  if (!is.data.frame(anchors) || !all(c("item", "k", "tau") %in% names(anchors)))
    stop("`anchors` must be a rasch() or rasch_cj() fit, or a data frame with ",
         "columns item, k and tau", call. = FALSE)
  item <- .role_text_values(anchors$item)
  k <- suppressWarnings(as.integer(anchors$k)); tau <- as.numeric(anchors$tau)
  if (anyNA(item) || anyNA(k) || any(!is.finite(tau)))
    stop("`anchors` has missing items, threshold numbers or values", call. = FALSE)
  out <- lapply(split(seq_along(item), factor(item, unique(item))), function(r) {
    kk <- k[r]
    if (!setequal(kk, seq_along(kk)))
      stop("`anchors` must give thresholds 1..m of item ", item[r[1]], call. = FALSE)
    tau[r][order(kk)]
  })
  out
}

.cj_persons <- function(data, anchors, id, comparisons, object_a, object_b,
                        winner, threshold_a, threshold_b, rankings, ranking,
                        item, rank, threshold, u_spec, items, na_codes, maxit,
                        tol, call) {
  if (is.null(data))
    stop("the person mode needs response data; a judgement of persons alone ",
         "is btl() or pl()", call. = FALSE)
  if (is.null(anchors))
    stop("the person mode needs `anchors`: the item thresholds from rasch() ",
         "or rasch_cj()", call. = FALSE)
  .check_na_codes(na_codes)
  tau_all <- .cj_anchor_list(anchors)
  notes <- character(0)

  # persons: an id column or vector, else the row names
  D <- if (is.data.frame(data)) data else as.data.frame(data)
  if (!is.null(id)) {
    if (is.character(id) && length(id) == 1L && id %in% names(D)) {
      ids <- .role_text_values(D[[id]]); D <- D[, setdiff(names(D), id), drop = FALSE]
    } else if (length(id) == nrow(D)) {
      ids <- .role_text_values(id)
    } else stop("`id` must name a column of `data` or give one identifier per row",
                call. = FALSE)
  } else {
    rn <- rownames(D)
    ids <- if (is.null(rn) || identical(rn, as.character(seq_len(nrow(D)))))
      sprintf("P%d", seq_len(nrow(D))) else .role_text_values(rn)
  }
  if (anyNA(ids) || any(!nzchar(ids)))
    stop("person identifiers must be non-missing and non-empty", call. = FALSE)
  if (anyDuplicated(ids))
    stop("person identifiers must be unique: ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "), call. = FALSE)
  N <- length(ids)

  # items: those in the data with anchored thresholds; categories must lie
  # within the anchored range, so nothing is rescored
  if (!is.null(items)) D <- D[, items, drop = FALSE]
  item_names <- names(D)
  no_anchor <- setdiff(item_names, names(tau_all))
  if (length(no_anchor))
    stop("item(s) in `data` without anchored thresholds: ",
         paste(no_anchor, collapse = ", "), call. = FALSE)
  unused <- setdiff(names(tau_all), item_names)
  if (length(unused))
    notes <- c(notes, sprintf("%d anchored item(s) not in the data ignored", length(unused)))
  tau_list <- tau_all[item_names]
  m <- vapply(tau_list, length, 1L)
  X <- as.matrix(D)
  code <- matrix(.missing_code_mask(as.vector(X), na_codes), N, ncol(X))
  Xi <- suppressWarnings(apply(X, 2, function(col) as.integer(as.character(col))))
  dim(Xi) <- dim(X); dimnames(Xi) <- list(NULL, item_names)
  Xi[code | (!is.na(Xi) & Xi < 0)] <- NA_integer_
  over <- sweep(Xi, 2, m, ">")
  if (any(over, na.rm = TRUE))
    stop("response(s) above the anchored top category in: ",
         paste(item_names[colSums(over, na.rm = TRUE) > 0], collapse = ", "),
         call. = FALSE)
  X <- Xi
  n_items <- rowSums(!is.na(X))
  raw <- rowSums(X, na.rm = TRUE); raw[n_items == 0L] <- NA_integer_
  max_raw <- as.vector((!is.na(X)) %*% m)
  has_resp <- n_items > 0L

  # judgement sources name persons; a person judged but not in the data has
  # no responses and is placed by the judgements alone
  judged_ids <- character(0)
  if (!is.null(comparisons))
    judged_ids <- c(judged_ids, .role_text_values(comparisons[[object_a]]),
                    .role_text_values(comparisons[[object_b]]))
  if (!is.null(rankings)) judged_ids <- c(judged_ids, .role_text_values(rankings[[item]]))
  judged_ids <- setdiff(unique(judged_ids[!is.na(judged_ids)]), ids)
  if (length(judged_ids)) {
    ids <- c(ids, judged_ids)
    X <- rbind(X, matrix(NA_integer_, length(judged_ids), ncol(X)))
    N <- length(ids)
    n_items <- rowSums(!is.na(X))
    raw <- rowSums(X, na.rm = TRUE); raw[n_items == 0L] <- NA_integer_
    max_raw <- as.vector((!is.na(X)) %*% m)
    has_resp <- n_items > 0L
    notes <- c(notes, sprintf("%d judged person(s) without responses, placed by the judgements",
                              length(judged_ids)))
  }
  cmp <- if (!is.null(comparisons))
    .cj_comparison_keys(comparisons, object_a, object_b, winner, threshold_a,
                        threshold_b, ids, NULL, notes, noun = "persons")
  else list(win = character(0), lose = character(0), n = 0L, notes = notes)
  notes <- cmp$notes
  rkl <- if (!is.null(rankings))
    .cj_ranking_keys(rankings, ranking, item, rank, threshold, ids, NULL, notes,
                     noun = "persons")
  else list(rk = list(), n = 0L, notes = notes)
  notes <- rkl$notes
  has_cmp <- cmp$n > 0L; has_rk <- rkl$n > 0L
  if (!is.null(comparisons) && !has_cmp)
    stop("`comparisons` has no usable rows", call. = FALSE)
  if (!is.null(rankings) && !has_rk)
    stop("`rankings` has no usable rows", call. = FALSE)
  win <- match(cmp$win, ids); lose <- match(cmp$lose, ids)
  rk <- lapply(rkl$rk, function(v) match(v, ids))

  # extreme persons: evidence in one direction only, judged among the
  # persons still in the estimation
  up0 <- has_resp & !is.na(raw) & raw > 0
  down0 <- has_resp & !is.na(raw) & raw < max_raw
  pr <- .cj_evidence_prune(rep(TRUE, N), win, lose, rk, up0, down0)
  in_set <- pr$active
  n_cmp_dropped <- pr$n_dropped
  win <- pr$win; lose <- pr$lose; rk <- pr$rk
  fit_idx <- which(in_set)
  n_fit <- length(fit_idx)
  if (!n_fit) stop("no person has evidence in both directions", call. = FALSE)
  if (any(!in_set))
    notes <- c(notes, sprintf(paste0("%d extreme person(s) set aside (evidence ",
      "in one direction only) and reported at the Warm estimate from their responses"),
      sum(!in_set)))
  if (n_cmp_dropped)
    notes <- c(notes, sprintf("%d comparison(s) involving an extreme person left out",
                              n_cmp_dropped))

  # every judged group must contain a person with responses, or it has no
  # place on the test scale
  fk <- n_fit; pos <- match(seq_len(N), fit_idx)
  adj <- matrix(FALSE, fk, fk)
  adj[cbind(pos[win], pos[lose])] <- TRUE
  for (v in rk) adj[pos[v], pos[v]] <- TRUE
  adj <- adj | t(adj); diag(adj) <- TRUE
  comp <- .cj_components(adj)
  resp_fit <- has_resp[fit_idx]
  orphan <- unique(comp[!resp_fit])
  orphan <- orphan[!vapply(orphan, function(g) any(resp_fit[comp == g]), NA)]
  if (length(orphan))
    stop("judged person(s) with no responses and no judged link to a person ",
         "with responses cannot be placed on the test scale: ",
         paste(ids[fit_idx][comp %in% orphan], collapse = ", "), call. = FALSE)

  # parameters: theta of the fitted persons, then the log units
  frames <- c("responses", if (has_cmp) "comparisons", if (has_rk) "rankings")
  free <- c(if (has_cmp) is.na(u_spec[["comparisons"]]),
            if (has_rk) is.na(u_spec[["rankings"]]))
  names(free) <- setdiff(frames, "responses")
  nu <- sum(free); u_pos <- integer(0)
  if (nu) u_pos <- fk + seq_len(nu)
  names(u_pos) <- names(free)[free]
  unit_of <- function(th, f) if (isTRUE(free[[f]])) exp(th[u_pos[[f]]]) else 1
  W <- matrix(0, fk, fk)
  if (length(win))
    W <- matrix(as.numeric(table(factor(pos[win], levels = seq_len(fk)),
                                 factor(pos[lose], levels = seq_len(fk)))), fk, fk)
  rk_f <- lapply(rk, function(v) pos[v])
  Xf <- X[fit_idx, , drop = FALSE]

  fn <- function(th) {
    theta <- th[seq_len(fk)]
    ll <- sum(.cj_person_parts(theta, Xf, tau_list)$ll)
    if (has_cmp) ll <- ll + .cj_bt_ll(theta, unit_of(th, "comparisons"), W)
    if (has_rk) ll <- ll + .cj_pl_ll(theta, unit_of(th, "rankings"), rk_f)
    ll
  }
  assemble <- function(th) {
    theta <- th[seq_len(fk)]; np <- fk + nu
    g <- numeric(np); H <- matrix(0, np, np)
    p <- .cj_person_parts(theta, Xf, tau_list)
    g[seq_len(fk)] <- p$g; diag(H)[seq_len(fk)] <- p$h
    add <- function(parts, f) {
      jo <- seq_len(fk)
      g[jo] <<- g[jo] + parts$g[jo]
      H[jo, jo] <<- H[jo, jo] + parts$H[jo, jo]
      if (isTRUE(free[[f]])) {
        k <- u_pos[[f]]
        g[k] <<- parts$g[fk + 1L]
        H[jo, k] <<- H[jo, k] + parts$H[jo, fk + 1L]
        H[k, jo] <<- H[k, jo] + parts$H[fk + 1L, jo]
        H[k, k] <<- parts$H[fk + 1L, fk + 1L]
      }
    }
    if (has_cmp) add(.cj_bt_parts(theta, unit_of(th, "comparisons"), W), "comparisons")
    if (has_rk) add(.cj_pl_parts(theta, unit_of(th, "rankings"), rk_f), "rankings")
    list(g = g, H = H)
  }
  gr <- function(th) assemble(th)$g
  he <- function(th) assemble(th)$H

  # separate calibrations: responses alone per person (the supremum, at an
  # infinite location, for a person with an extreme score held by the
  # judgements), and each judgement frame alone over the persons it reaches
  finite_resp <- has_resp & !is.na(raw) & raw > 0 & raw < max_raw
  ml <- .cj_person_ml(X, tau_list, finite_resp, maxit, tol)
  sep <- list(); judged <- list(); blocks <- list()
  none <- rep(FALSE, fk)
  frame_alone <- function(f, win_f, lose_f, rk_f) {
    pr <- .cj_evidence_prune(rep(TRUE, fk), win_f, lose_f, rk_f, none, none)
    jo <- which(pr$active)
    reached <- unique(c(win_f, lose_f, unlist(rk_f)))
    n_ext <- length(setdiff(reached, jo))
    if (n_ext)
      notes <<- c(notes, sprintf(paste0("%d person(s) extreme within the %s ",
        "alone (evidence in one direction only there): no separate location ",
        "from that frame"), n_ext, f))
    if (!length(jo)) return(NULL)
    pos_f <- match(seq_len(fk), jo)
    Wf <- matrix(as.numeric(table(factor(pos_f[pr$win], levels = seq_along(jo)),
                                  factor(pos_f[pr$lose], levels = seq_along(jo)))),
                 length(jo), length(jo))
    rkf <- lapply(pr$rk, function(v) pos_f[v])
    adj_f <- Wf + t(Wf) > 0
    for (v in rkf) adj_f[v, v] <- TRUE
    cp <- .cj_components(adj_f)
    ll <- if (f == "comparisons") function(l) .cj_bt_ll(l, 1, Wf)
      else function(l) .cj_pl_ll(l, 1, rkf)
    parts <- if (f == "comparisons") function(l) .cj_bt_parts(l, 1, Wf)
      else function(l) .cj_pl_parts(l, 1, rkf)
    nj <- length(jo)
    fit <- .cj_fit_alone(ll, function(l) { p <- parts(l)
                           list(g = p$g[seq_len(nj)], H = p$H[seq_len(nj), seq_len(nj), drop = FALSE]) },
                         .cj_basis_groups(cp), maxit = maxit, tol = tol)
    list(fit = fit, jo = jo, comp = cp)
  }
  if (has_cmp) {
    a <- frame_alone("comparisons", pos[win], pos[lose], list())
    if (!is.null(a)) {
      sep$comparisons <- a$fit; judged$comparisons <- a$jo; blocks$comparisons <- a$comp
    }
  }
  if (has_rk) {
    a <- frame_alone("rankings", integer(0), integer(0), rk_f)
    if (!is.null(a)) {
      sep$rankings <- a$fit; judged$rankings <- a$jo; blocks$rankings <- a$comp
    }
  }

  # start at the response estimates where finite, else at the judged
  # position within the frame, else zero
  start <- rep(0, fk)
  fr <- finite_resp[fit_idx]
  start[fr] <- ml$theta[fit_idx][fr]
  for (f in names(sep)) {
    s <- sep[[f]]; jo <- judged[[f]]
    if (!s$converged) next
    fill <- !fr & seq_len(fk) %in% jo
    if (any(fill)) start[fill] <- s$par[match(which(fill), jo)]
  }
  th0 <- c(start, rep(0, nu))
  fit <- .cj_newton(th0, fn, gr, he, maxit, tol)
  th <- fit$par; theta <- th[seq_len(fk)]
  covth <- tryCatch(solve(-fit$H), error = function(e) matrix(NA_real_, fk + nu, fk + nu))
  cov_t <- covth[seq_len(fk), seq_len(fk), drop = FALSE]
  se_t <- sqrt(pmax(diag(cov_t), 0))
  if (!fit$converged) {
    notes <- c(notes, "estimation did not converge; standard errors and tests withheld")
    se_t[] <- NA_real_; cov_t[] <- NA_real_
  }
  for (f in names(free)[free]) {
    if (unit_of(th, f) > 1e3)
      notes <- c(notes, sprintf(paste0("the %s unit ran away: the judgements ",
        "never disagree with an ordering of the persons the responses allow, ",
        "so they are too few per person to give the unit a maximum; ",
        "fix it with units = c(%s = 1)"), f, f))
    if (unit_of(th, f) < 1e-3)
      notes <- c(notes, sprintf(paste0("the %s unit collapsed towards zero: ",
        "the judgements carry no information about the persons, or their ",
        "orientation is reversed (the winner or first rank should be the ",
        "person with the higher location)"), f))
  }
  units_tab <- data.frame(frame = frames, unit = 1, se = NA_real_,
                          estimated = FALSE, stringsAsFactors = FALSE)
  for (f in names(free)) {
    r <- match(f, units_tab$frame)
    units_tab$unit[r] <- unit_of(th, f)
    units_tab$estimated[r] <- isTRUE(free[[f]])
    if (isTRUE(free[[f]]) && fit$converged) {
      k <- u_pos[[f]]
      units_tab$se[r] <- units_tab$unit[r] * sqrt(max(covth[k, k], 0))
    }
  }

  # invariance: each frame's separate locations against the response
  # locations, person by person. There is no likelihood ratio test here:
  # the separate model has a location per person per frame, so its
  # parameters grow with the persons and the ratio is not chi-square.
  persons <- data.frame(person = ids, n_items = n_items, raw = raw,
                        max_raw = max_raw, location = NA_real_, se = NA_real_,
                        location_responses = ml$theta, se_responses = ml$se,
                        stringsAsFactors = FALSE, row.names = NULL)
  persons$location[fit_idx] <- theta; persons$se[fit_idx] <- se_t
  ref_theta <- ml$theta[fit_idx]; ref_se <- ml$se[fit_idx]
  inv_tab <- NULL
  for (f in names(sep)) {
    u <- unit_of(th, f); s <- sep[[f]]; jo <- judged[[f]]; cp <- blocks[[f]]
    nj <- length(jo)
    # centring within each block, over the persons finite in both frames
    both <- is.finite(ref_theta[jo])
    C <- matrix(0, nj, nj)
    for (g in unique(cp)) {
      rows <- which(cp == g & both)
      if (length(rows)) C[rows, rows] <- -1 / length(rows)
    }
    diag(C) <- diag(C) + 1
    l_f <- as.vector(C %*% (s$par / u))
    col <- rep(NA_real_, N); col[fit_idx[jo]] <- l_f
    persons[[paste0("location_", f)]] <- col
    if (fit$converged) {
      r_t <- ref_theta[jo]; r_t[!both] <- 0
      l_r <- as.vector(C %*% r_t)
      diff <- l_f - l_r
      V <- C %*% (diag(ifelse(both, ref_se[jo]^2, 0), nj) + s$cov / u^2) %*% t(C)
      se_d <- sqrt(pmax(diag(V), 0))
      z <- diff / se_d
      # a contrast needs a response location and at least one other such
      # person in the block to centre on
      ok <- both & vapply(cp, function(g) sum(both[cp == g]) > 1L, NA)
      z[!is.finite(z) | !ok] <- NA_real_
      diff[!ok] <- NA_real_; l_r[!ok] <- NA_real_
      p <- 2 * stats::pnorm(-abs(z))
      inv_tab <- rbind(inv_tab, data.frame(
        frame = f, person = ids[fit_idx][jo], reference = l_r, judgements = l_f,
        difference = diff, se = se_d, z = z, p = p,
        p_adj = stats::p.adjust(p, "holm"), stringsAsFactors = FALSE,
        row.names = NULL))
    }
  }
  persons$comparisons <- seq_len(N) %in% fit_idx[judged$comparisons]
  persons$rankings <- seq_len(N) %in% fit_idx[judged$rankings]
  persons$extreme <- !in_set
  # extreme persons at the Warm estimate from their responses, as rasch()
  # reports them
  if (any(!in_set & has_resp)) {
    wle <- .person_estimates(X[!in_set & has_resp, , drop = FALSE], tau_list)
    persons$location[!in_set & has_resp] <- wle$theta
    persons$se[!in_set & has_resp] <- wle$se
  }
  dimnames(cov_t) <- list(ids[fit_idx], ids[fit_idx])

  structure(list(persons = persons, units = units_tab,
                 invariance = list(persons = inv_tab),
                 anchors = data.frame(item = rep(item_names, m),
                                      k = unlist(lapply(m, seq_len), use.names = FALSE),
                                      tau = unlist(tau_list, use.names = FALSE),
                                      stringsAsFactors = FALSE),
                 cov = cov_t, loglik = fit$ll,
                 converged = fit$converged, iterations = fit$iterations,
                 n = c(persons = sum(has_resp), comparisons = cmp$n, rankings = rkl$n),
                 reference = "responses", mode = "persons", notes = notes,
                 call = call),
            class = "rasch_cj")
}
