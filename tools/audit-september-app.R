# Display-only corrections: compare the joint estimator with the audited
# baseline. Run from the repository root, with git history available:
# Rscript tools/audit-september-app.R
pkgload::load_all(quiet = TRUE)
old <- new.env(parent = asNamespace("rasch"))
src <- system2("git", c("show", "59ab6ba:R/rasch-cj.R"), stdout = TRUE)
if (!is.null(attr(src, "status"))) stop("the baseline git revision is unavailable")
eval(parse(text = src), old)
set.seed(9262026)
delta <- setNames(seq(-2, 2, length.out = 6), paste0("I", 1:6))
theta <- runif(400, -3, 3)
X <- sapply(delta, function(d) rbinom(400, 1, plogis(theta - d)))
make_pairs <- function(ids, loc) {
  p <- t(combn(ids, 2)); p <- p[sample(nrow(p), 600, TRUE), ]
  data.frame(object_a = p[, 1], object_b = p[, 2],
    winner = ifelse(runif(nrow(p)) < plogis(loc[p[, 1]] - loc[p[, 2]]), p[, 1], p[, 2]))
}
partial <- make_pairs(names(delta)[4:6], delta)
disconnected <- rbind(partial, make_pairs(names(delta)[1:3], delta))
m <- c(1, 2, 3, 2, 1, 3)
tau <- lapply(seq_along(m), function(i) delta[i] + seq(-.5, .5, length.out = m[i]))
Xp <- sapply(tau, function(t) vapply(theta, function(b) {
  pr <- exp((0:length(t)) * b - c(0, cumsum(t)))
  sample(0:length(t), 1L, prob = pr)
}, 0L))
colnames(Xp) <- names(delta)
obj <- data.frame(item = rep(names(delta), m), k = sequence(m))
keys <- paste(obj$item, obj$k, sep = ":")
ct <- make_pairs(keys, setNames(unlist(tau), keys))
ct$threshold_a <- obj$k[match(ct$object_a, keys)]
ct$threshold_b <- obj$k[match(ct$object_b, keys)]
ct$object_a <- obj$item[match(ct$object_a, keys)]
ct$object_b <- obj$item[match(ct$object_b, keys)]
scenarios <- list(subset = list(X, partial),
  disconnected = list(X, disconnected),
  polytomous_items = list(Xp, partial),
  polytomous_thresholds = list(Xp, ct),
  multiple_tests = list(list(A = X[, 1:3], B = X[, 4:6]), make_pairs(names(delta), delta)))
fields <- c("items", "thresholds", "objects", "units", "invariance", "anchors",
            "cov", "cov_items", "loglik", "loglik_separate", "converged", "iterations", "n")
for (name in names(scenarios)) {
  args <- c(scenarios[[name]], list(units = c(comparisons = 1)))
  a <- do.call(old$rasch_cj, args); b <- do.call(rasch_cj, args)
  inference <- b
  inference$invariance$n_contrasts <- NULL # new reporting metadata only
  stopifnot(a$converged, b$converged, identical(a[fields], inference[fields]))
  loc <- b$frame_locations
  stopifnot(all(is.finite(loc$location)))
  for (f in b$tests) {
    at <- loc[loc$frame == f, ]
    idx <- match(at$item, b$items$item)
    stopifnot(length(idx) > 0, abs(mean(at$location) - mean(b$items$location[idx])) < 1e-10)
  }
  pdf(NULL); plot_cj(b); dev.off()
  cat(name, ": estimates, covariance, tests and likelihood identical; aligned plot OK\n")
}
