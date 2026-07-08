# Case study: wording effects as frame units in survey data
# ===========================================================================
# Self-report scales routinely mix positively worded items with reverse-
# scored negatively worded items, and the reverse-scored items typically
# discriminate less. Under the extended frame of reference model this is a
# difference in the natural unit of two item sets: rho = alpha_set, with the
# wording defining the sets. This script applies the model to the public
# Rosenberg Self-Esteem Scale dataset collected by the Open Source
# Psychometrics Project (about 47,000 respondents; downloaded from the
# source at run time, not redistributed here) and quantifies why the unit
# difference matters: persons with the same raw score receive different
# measures, and group comparisons shift.
#
# Run time: a few minutes (subsampled to 6,000 respondents).
# ===========================================================================
library(rasch)
set.seed(2026)

src <- "http://openpsychometrics.org/_rawdata/RSE.zip"
tmp <- tempfile(fileext = ".zip")
download.file(src, tmp, quiet = TRUE)
unzip(tmp, exdir = dirname(tmp))
d <- read.csv(file.path(dirname(tmp), "RSE", "data.csv"), sep = "\t")

items <- paste0("Q", 1:10)
positive <- c("Q1", "Q2", "Q4", "Q6", "Q7")
negative <- c("Q3", "Q5", "Q8", "Q9", "Q10")   # reverse-scored below
X <- as.matrix(d[, items])
X[X == 0] <- NA                                # 0 = no answer
X[, negative] <- 5 - X[, negative]             # score all items in one direction
keep <- rowSums(is.na(X)) == 0 & d$gender %in% c(1, 2)
df <- data.frame(X[keep, ], gender = c("male", "female")[d$gender[keep]],
                 check.names = FALSE)
df <- df[sample(nrow(df), 6000), ]

# equal-unit Rasch versus wording-set EFRM
f0 <- rasch(df, factors = "gender", items = items)
f1 <- rasch_efrm(df, items = items, groups = rep("all", nrow(df)),
                 item_sets = list(positive = positive, negative = negative))

print(f1$alpha_table, digits = 3)
print(f1$efrm_vs_rasch$unit_tests, digits = 3)
# Note: with a single person group the pairwise 2*delta-ll is invariant to
# the set units by construction; the Wald tests above carry the evidence.

# why it matters -------------------------------------------------------------
ok <- !f0$person$extreme & !f1$person$extreme
t0 <- f0$person$theta[ok]; t1 <- f1$person$theta[ok]
raw <- f0$person$raw[ok]
cat(sprintf("max measure spread within one raw score: %.3f logits\n",
            max(tapply(t1, raw, function(x) diff(range(x))), na.rm = TRUE)))
cat(sprintf("persons moving more than 0.1 logits: %.1f%%\n",
            100 * mean(abs(t1 - t0) > 0.1)))
g <- df$gender[ok]
cat(sprintf("male-female gap: equal units %.3f logits, wording units %.3f logits\n",
            mean(t0[g == "male"]) - mean(t0[g == "female"]),
            mean(t1[g == "male"]) - mean(t1[g == "female"])))

# the visual: a positive and a negative item with their different slopes
op <- par(mfrow = c(1, 2))
plot_icc(f1, "Q6:all")    # positive wording: steeper (larger unit)
plot_icc(f1, "Q9:all")    # negative wording: flatter (smaller unit)
par(op)

# cross-check against free slopes, and sensitivity ---------------------------
# A generalized partial credit model (the polytomous two-parameter model)
# frees one slope per item; its geometric-mean slope ratio between the
# wording sets (about 1.25 here) agrees closely with the single EFRM unit
# ratio (1.266), so the frame structure captures the set-level
# discrimination with one parameter instead of nine. The per-item slopes
# also localise the anomalies: Q8 ("I wish I could have more respect for
# myself"), the scale's well-known ambivalent item, discriminates far below
# the other negatives, and Q4 is the weakest positive. The unit ratio
# remains significantly above one without them (1.115 dropping Q8; 1.196
# also dropping Q4), so the wording effect is real but Q8 inflates it.
if (requireNamespace("mirt", quietly = TRUE)) {
  m <- mirt::mirt(as.data.frame(X[keep, ][sample(sum(keep), 6000), ]), 1,
                  itemtype = "gpcm", verbose = FALSE)
  print(mirt::coef(m, simplify = TRUE)$items[, "a1", drop = FALSE])
}
for (drop in list("Q8", c("Q8", "Q4"))) {
  p2 <- setdiff(positive, drop); n2 <- setdiff(negative, drop)
  f <- rasch_efrm(df, items = c(p2, n2), groups = rep("all", nrow(df)),
                  item_sets = list(positive = p2, negative = n2))
  cat(sprintf("dropping %s: alpha ratio = %.3f\n",
              paste(drop, collapse = "+"),
              f$alpha_table$alpha[f$alpha_table$set == "positive"] /
                f$alpha_table$alpha[f$alpha_table$set == "negative"]))
}
