# Equate two paired-comparison calibrations through their common objects

Places two Bradley–Terry–Luce calibrations on a common origin using
their shared objects, then tests the shared objects for drift. The
second calibration may be a fitted model or an object bank.

## Usage

``` r
btl_equate(
  fit1,
  fit2,
  alpha = 0.05,
  p_adjust = "holm",
  independent = NULL,
  shift = c("mean", "none")
)
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
  optionally `se`; object names and column names must be unique. Numeric
  fields may be numeric columns, numeric text, or factors with numeric
  labels; other column classes are refused. Locations must be finite.
  Bank-based drift inference with an estimated mean shift requires the
  joint location covariance as a square matrix in
  `attr(fit2, "cov_location")`, ordered like the bank rows (or named by
  object), unless the bank is treated as fixed with zero SEs. Marginal
  standard errors are sufficient with `shift = "none"`. A bank whose
  covariance was estimated from a finite number of independent sampling
  units may carry their residual degrees of freedom in
  `attr(fit2, "df_location")` as one positive numeric value. For a
  polytomous fit the bank must carry `attr(bank, "m")` matching the
  number of fitted score steps.

- alpha:

  Significance level for the (multiplicity-adjusted) drift tests.

- p_adjust:

  Adjustment for the common-object tests, passed to
  [`stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html). The
  default is `"holm"`. A common object remains in the family when its
  drift probability is unavailable.

- independent:

  Whether the calibrations have independent judges and comparisons. For
  two fitted objects the default `NULL` withholds drift tests until
  independence is stated explicitly. Bank tables are treated as
  independent unless `FALSE` is supplied. Dependent calibrations require
  a joint or paired bootstrap for inference.

- shift:

  `"mean"` (default) estimates the origin shift from the common objects;
  `"none"` compares raw locations when both calibrations have already
  been placed on the same externally anchored scale.

## Value

A list of class `"rasch_btl_equate"`: the comparison `table` (per common
object: object, both locations and standard errors, their `difference`,
the `shifted_difference` against the estimated origin, the pooled
`se_diff`, `t`, raw and adjusted `p`, and the `drifting` flag); the
estimated `shift`, its `shift_method` and `shift_se`; `equated`, the
second calibration's full object table re-expressed on `fit1`'s scale;
the number of common objects `n_common`; the number usable for inference
`n_inference`; whether inference was available `inferential`; `alpha`;
`p_adjust`; the requested `shift_setting`; and `notes`. A drift
probability is withheld when its contrast has zero estimated
uncertainty. Such an object remains in the multiplicity family.

## Details

Let \\d_j\\ be the location difference for common object \\j\\ and
\\v_j\\ its marginal variance. With `shift = "mean"`, the origin shift
is the precision-weighted mean \$\$\hat s=\frac{\sum_j d_j/v_j}{\sum_j
1/v_j}.\$\$ If fewer than two common objects have usable variances but
at least two have finite locations, their unweighted mean difference is
returned as a descriptive fallback and recorded in `shift_method`. Every
error built from those precision weights conditions on them as if the
two calibrations' standard errors were known: `shift_se`, each drift
contrast's `se_diff` (and so its probability and `drifting` flag), and
the equated location errors that carry the shift. Weight uncertainty
adds a positive term all of them omit, so they understate uncertainty –
intervals under-cover, drift probabilities run small – when the
calibrations rest on few judges; a judge resample of both calibrations
is the weight-aware alternative. An exact common anchor determines the
shift even when it is the only common object with usable uncertainty.
Each object is tested using its shifted difference \\d_j-\hat s\\. The
covariance calculation retains the dependence induced by the sum-zero
constraints. Drift tests then require independent calibrations and at
least three common objects with usable, positive-semidefinite joint
covariance information. Two common objects identify a descriptive origin
shift, but do not support an object-drift test. With `shift = "none"`,
the origin is fixed before the comparison and each object's variance is
the sum of its two marginal variances; joint covariance information and
a three-object link are unnecessary. One common object is sufficient for
that fixed-origin comparison; estimating a shift still requires at least
two. A judge-clustered ordinary BTL covariance, or a BTL–EFRM covariance
from the judge bootstrap, uses finite judge-cluster degrees of freedom.
A contrast involving only fixed external anchors has exact zero
covariance from that calibration and therefore uses infinite degrees of
freedom even when inference for its estimated objects is unavailable. A
BTL–EFRM location outside the reference set is also limited by the
weakest edge on its strongest supported path to that reference. A
comparison-level parametric-bootstrap BTL–EFRM covariance uses the
asymptotic normal reference instead. Conditional frame errors are
preliminary and do not support drift inference, including comparisons on
a fixed origin. Binary fits have no threshold parameters, so their
recorded threshold structure does not affect compatibility. Polytomous
fits must use the same category scale and threshold structure.

The `equated` table includes uncertainty in the estimated shift. For
independent calibrations, with \\y_j=b_j+\hat s\\,
\$\$\operatorname{Var}(y_j)=\operatorname{Var}(b_j)+
\operatorname{Var}(\hat s)+2\operatorname{Cov}(b_j,\hat s).\$\$ These
SEs are withheld if joint uncertainty is unavailable. A fixed shift
(`shift = "none"` or an exact common anchor) leaves supported original
SEs unchanged. SEs from a conditional frame reference are withheld in
the equated bank. When available, the table carries its full covariance
in `attr(equated, "cov_location")` and conservative finite sampling-unit
degrees of freedom in `attr(equated, "df_location")`. An equated bank is
not independent of either calibration used to construct it.

The common-object set should contain a stable majority. If most common
objects move in the same direction, the estimated shift follows them and
stable objects can appear to drift. In that case, repeat the equating
with a substantively justified anchor set.

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
eq <- btl_equate(sim(paste0("O", 1:7)), sim(paste0("O", 2:8)),
                  independent = TRUE)
eq$table
#>   object location_1      se_1 location_2      se_2 difference
#> 1     O2 -1.0161822 0.1424913 -1.4290667 0.1437563  0.4128844
#> 2     O3 -0.5553217 0.1363539 -1.1019358 0.1420704  0.5466140
#> 3     O4 -0.1387659 0.1454378 -0.5398557 0.1370181  0.4010897
#> 4     O5  0.7692300 0.1424835  0.1589980 0.1227805  0.6102320
#> 5     O6  1.2374851 0.1430865  0.3920252 0.1300981  0.8454599
#> 6     O7  1.8209314 0.1682478  0.9978767 0.1373542  0.8230547
#>   shifted_difference   se_diff           t  df         p p_adj drifting
#> 1       -0.191440169 0.1997993 -0.95816214 Inf 0.3379810     1    FALSE
#> 2       -0.057710563 0.1941698 -0.29721701 Inf 0.7663008     1    FALSE
#> 3       -0.203234869 0.1945965 -1.04439116 Inf 0.2963044     1    FALSE
#> 4        0.005907421 0.1784021  0.03311296 Inf 0.9735845     1    FALSE
#> 5        0.241135288 0.1865883  1.29233859 Inf 0.1962399     1    FALSE
#> 6        0.218730044 0.2194493  0.99672260 Inf 0.3188992     1    FALSE
```
