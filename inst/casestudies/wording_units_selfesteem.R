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
set.seed(25)

src <- "https://openpsychometrics.org/_rawdata/RSE.zip"
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
set.seed(26)
f1 <- rasch_efrm(df, items = items, groups = rep("all", nrow(df)),
                 item_sets = list(positive = positive, negative = negative),
                 se_method = "hybrid", boot_reps = 300)

print(f1$alpha_table, digits = 3)
print(f1$efrm_vs_rasch$unit_tests, digits = 3)
# Note: with a single person group the pairwise 2*delta-ll is invariant to
# the set units by construction; the Wald tests above carry the evidence.
# The two set units are centred, so the second row states the same
# hypothesis as the first: it keeps its own probability and leaves the
# adjusted one to the row it restates.

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

# sensitivity: which items carry the set-unit difference? ---------------------
# The single set unit is an average over the items in a set, so a badly
# behaved item moves it. Refit without the suspects rather than asserting
# what would happen.
unit_ratio <- function(drop = character(), seed = 31L) {
  it <- setdiff(items, drop)
  set.seed(seed)
  f <- rasch_efrm(df[, c(it, "gender")], items = it,
                  groups = rep("all", nrow(df)),
                  item_sets = list(positive = setdiff(positive, drop),
                                   negative = setdiff(negative, drop)),
                  se_method = "hybrid", boot_reps = 300)
  # The two centred set units are one hypothesis stated twice, so the row
  # that restates it withholds its adjusted probability. Take the minimum
  # over the rows that carry one, and report nothing if none does.
  pa <- f$efrm_vs_rasch$unit_tests$p_adj
  c(ratio = unname(f$alpha_table$alpha[f$alpha_table$set == "positive"] /
                   f$alpha_table$alpha[f$alpha_table$set == "negative"]),
    p_adj = if (any(is.finite(pa))) min(pa, na.rm = TRUE) else NA_real_)
}
apriori <- c("Q8", "Q4")   # the usual suspects, named in advance
print(signif(rbind(all_items = unit_ratio(seed = 31),
                   drop_first = unit_ratio(apriori[1], seed = 32),
                   drop_both = unit_ratio(apriori, seed = 33)), 4))

# let the model nominate the suspects ----------------------------------------
# The drops above were chosen a priori, from what is already known about the
# scale. The fitted model can nominate them instead: an item that shares no
# unit with its set misfits within it, and the standardised fit residual
# ranks that misfit. Rank on it rather than threshold on it -- a fixed cut
# states detectability, not magnitude, so at this sample size it selects
# most of the instrument, as the count below shows.
fr <- f1$items[order(abs(f1$items$fit_resid), decreasing = TRUE), ]
fr$item <- sub(":.*$", "", fr$item)     # frame models name items by frame
fr$set <- ifelse(fr$item %in% positive, "positive", "negative")
print(fr[, c("item", "set", "fit_resid", "infit_z")], digits = 3,
      row.names = FALSE)
cat(sprintf("items clearing a fixed |fit_resid| > 2 cut: %d of %d\n",
            sum(abs(fr$fit_resid) > 2, na.rm = TRUE), nrow(fr)))

# The ranking selects, from the data alone, the items named in advance
# above, so the table already printed is also the data-driven drop sequence
# and needs no refitting to reproduce.
cat(sprintf("ranked first and second: %s; named in advance: %s\n",
            paste(fr$item[1:2], collapse = " "),
            paste(apriori, collapse = " ")))

# Read that table as a sensitivity sequence rather than an automatic deletion
# rule. Removing Q8 reduces the ratio from about 1.32 to 1.08, although the
# large sample still gives an adjusted p-value near .015: the two centred
# set units are one hypothesis, so Holm has nothing to adjust it against.
# Removing Q4 as well moves the ratio away from one again. The conclusion is
# therefore that Q8 carries most, but not all, of the original difference and
# that the result is sensitive to the composition of these short wording sets.

# cross-check against free slopes --------------------------------------------
# A generalized partial credit model frees one slope per item, on the same
# respondents, so its per-item slopes are an independent reading of the same
# data. They localise the same anomalies: Q8 ("I wish I could have more
# respect for myself"), the scale's well-known ambivalent item, discriminates
# far below the other negatives, and Q4 is the weakest positive. The two
# orderings printed at the end are the useful comparison -- the frame model's
# fit residuals and the free slopes are computed from different quantities
# and agree on which items are extreme.
#
# Reading the analyses together: with all ten items the positive/negative unit
# ratio is about 1.32. Dropping Q8 reduces it to about 1.08; dropping Q4 as well
# raises it to about 1.19. The set-level effect is therefore carried mainly by
# individual anomalous items rather than by wording alone: a conclusion the
# single-parameter frame model cannot reach on its own, which is why the
# free-slope cross-check belongs here.
#
# The GPCM's geometric-mean slope ratio between the sets lands in the same
# region as the EFRM unit ratio, which is the other useful comparison: one
# parameter per set reproduces what ten free slopes say about the sets on
# average. It is an average, though -- with slopes as spread as Q8's and
# Q6's, a single set unit summarises a heterogeneous set, and the two ratios
# are close rather than equal.
if (requireNamespace("mirt", quietly = TRUE)) {
  m <- mirt::mirt(as.data.frame(df[, items]), 1, itemtype = "gpcm",
                  verbose = FALSE)
  a1 <- mirt::coef(m, simplify = TRUE)$items[, "a1"]
  print(round(a1, 3))
  gm <- function(z) exp(mean(log(z)))
  cat(sprintf("GPCM geometric-mean slope ratio positive/negative: %.3f\n",
              gm(a1[positive]) / gm(a1[negative])))
  cat(sprintf("free slopes, flattest first: %s\n",
              paste(names(sort(a1)), collapse = " ")))
  cat(sprintf("fit residuals, largest first: %s\n",
              paste(fr$item, collapse = " ")))
} else {
  message("install the 'mirt' package for the free-slope cross-check: ",
          "install.packages('mirt')")
}
