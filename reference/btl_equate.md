# Equate two paired-comparison calibrations through their common objects

Compares the object locations of a Bradley-Terry-Luce fit with those of
a second fit (or a banked table of object locations), matched by object
name. The use case is standards maintenance and comparative judgement
across panels or years: the same scripts, performances, or products are
judged by different panels – or by one panel in successive years – and a
common set of anchor objects is carried through so that the two rounds
land on a single scale.

## Usage

``` r
btl_equate(fit1, fit2, alpha = 0.05, p_adjust = "holm")
```

## Arguments

- fit1:

  A fitted object from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md): the
  calibration whose scale (origin) the equating targets.

- fit2:

  A second
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) fit,
  or a bank: a data frame with columns `object`, `location`, and
  optionally `se`.

- alpha:

  Significance level for the (multiplicity-adjusted) drift tests.

- p_adjust:

  Multiple-comparison adjustment across the common objects, passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) (default
  `"holm"`).

## Value

A list of class `"rasch_btl_equate"`: the comparison `table` (per common
object: object, both locations and standard errors, their `difference`,
the `shifted_difference` against the estimated origin, the pooled
`se_diff`, `t`, raw and adjusted `p`, and the `drifting` flag); the
estimated `shift` and its `shift_se`; `equated`, the second
calibration's full object table re-expressed on `fit1`'s scale; the
number of common objects `n_common`; `alpha`; `p_adjust`; and `notes`.

## Details

Each calibration is identified by the sum-zero constraint, but the two
constraints are imposed over *different* object sets, so the origins do
not coincide even when the shared objects are unchanged: each scale is
centred on the mean of a different collection. A scale shift between the
two origins is therefore estimated by the precision-weighted mean
difference over the common objects, and each common object is then
tested against the shifted identity line. A flagged object shows drift –
a script the two panels valued differently, or a standard that moved
between years – and weakens the equating link; the surviving objects
carry the second calibration's whole scale onto the first
(`loc2 + shift`).

The single shift presumes the drifting objects are a *minority*. When
most of the common objects have genuinely moved, the precision-weighted
shift is pulled toward the movers and the drift tests can invert – the
stable anchors flag as the apparent drifters. Read wholesale flagging
(several objects, one direction) as a contaminated link, not as evidence
about the individual objects; equate through a vetted anchor subset
instead. The `shift_se` accounts for the covariance of the location
estimates within each calibration (each is sum-zero constrained, so its
locations are not independent).

## References

Bramley, T. (2007). Paired comparison methods. In P. Newton, J. Baird,
H. Goldstein, H. Patrick, & P. Tymms (Eds.), *Techniques for monitoring
the comparability of examination standards* (pp. 246-294). London:
Qualifications and Curriculum Authority.

## Examples

``` r
set.seed(1)
beta <- setNames(seq(-2, 2, length.out = 8), paste0("O", 1:8))
sim <- function(objs) {
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  btl(d, "a", "b", "win")
}
eq <- btl_equate(sim(paste0("O", 1:7)), sim(paste0("O", 2:8)))
eq$table
#>   object location_1      se_1 location_2      se_2 difference
#> 1     O2 -1.0161822 0.1424913 -1.4290667 0.1437563  0.4128844
#> 2     O3 -0.5553217 0.1363539 -1.1019358 0.1420704  0.5466140
#> 3     O4 -0.1387659 0.1454378 -0.5398557 0.1370181  0.4010897
#> 4     O5  0.7692300 0.1424835  0.1589980 0.1227805  0.6102320
#> 5     O6  1.2374851 0.1430865  0.3920252 0.1300981  0.8454599
#> 6     O7  1.8209314 0.1682478  0.9978767 0.1373542  0.8230547
#>   shifted_difference   se_diff           t         p p_adj drifting
#> 1       -0.191440169 0.1997993 -0.95816214 0.3379810     1    FALSE
#> 2       -0.057710563 0.1941698 -0.29721701 0.7663008     1    FALSE
#> 3       -0.203234869 0.1945965 -1.04439116 0.2963044     1    FALSE
#> 4        0.005907421 0.1784021  0.03311296 0.9735845     1    FALSE
#> 5        0.241135288 0.1865883  1.29233859 0.1962399     1    FALSE
#> 6        0.218730044 0.2194493  0.99672260 0.3188992     1    FALSE
```
