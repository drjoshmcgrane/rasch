# Unexpected judgements of one judge, pair by pair

The comparison-level companion of
[`judge_surprise`](https://drjoshmcgrane.github.io/rasch/reference/judge_surprise.md).
Each pair the judge met is oriented to its stronger object (higher
consensus location) and given a standardised residual: `z < 0` means the
stronger object won less than its lead predicts – the judge backed the
underdog. Approximate two-sided normal probabilities are adjusted by
Holm across matchups meeting `min_n`. A surprise is an eligible matchup
with a negative residual whose adjusted probability passes the level
represented by `flag_z`. The fitted model must have converged. An
adequately sampled matchup with unavailable residual inference remains
in the adjustment family. Locations within \\10^{-10}\\ logits are
treated as tied: neither object is called stronger, the two-sided
residual probability remains descriptive and in the Holm family, and
`surprise` is false. A tied row is oriented toward its non-negative
residual for display only.

## Usage

``` r
judge_pair_surprise(fit, judge, min_n = 1L, flag_z = 1.96)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) with
  judges.

- judge:

  The judge to profile.

- min_n:

  Pairs met fewer times are shown but never flagged.

- flag_z:

  Absolute normal-residual threshold defining the familywise flagging
  level; 1.96 corresponds to an adjusted two-sided probability of
  approximately 0.05.

## Value

A list of class `"rasch_btl_judge_pairs"`: `pairs` (per matchup: the
stronger and weaker object and their locations when `tied` is false (at
a tie these columns only orient the row), the location `gap`, tie
indicator `tied`, times met `n`, residual `z`, approximate `p`,
Holm-adjusted `p_adj`, the `net_winner`, and the `surprise` flag);
`all_locations`; the `judge` and settings.

## Examples

``` r
set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
pr <- t(utils::combn(objs, 2))
d <- data.frame(a = rep(pr[, 1], each = 12), b = rep(pr[, 2], each = 12))
d$judge <- sample(paste0("J", 1:5), nrow(d), TRUE)
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
judge_pair_surprise(btl(d, "a", "b", "win", judge = "judge"), "J1")
#> Judge J1: 40 comparisons over 15 matchups
#> No matchup went against the consensus beyond noise.
```
