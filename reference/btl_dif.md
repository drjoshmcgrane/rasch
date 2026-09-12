# DIF analysis for paired comparisons

Tests whether object locations differ across groups of judges. Several
judge factors can be fitted jointly, with optional factor-by-factor
interactions. Uniform DIF is a judge-factor effect; non-uniform DIF is
its interaction with opponent-strength band.

## Usage

``` r
btl_dif(
  fit,
  factors,
  objects = NULL,
  effects = c("main", "factorial"),
  p_adjust = "holm",
  alpha = 0.05,
  flag_logits = 0.5,
  min_n = 20,
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- fit:

  An ordinary paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- factors:

  One judge factor, or a named list containing several. Each factor may
  have one value per comparison row or be a vector named by every judge
  in the fit.

- objects:

  Objects to test; all by default.

- effects:

  `"main"` (default) models several factors additively (each factor's
  main effect and its band interaction); `"factorial"` also crosses the
  factors with one another.

- p_adjust:

  Multiplicity adjustment over all object-by-term tests; the
  resolved-contrast probabilities are adjusted separately in one pool
  over all objects, terms, and contrasts.

- alpha:

  Significance level for adjusted probabilities.

- flag_logits:

  Absolute resolved difference flagged as practically significant.

- min_n:

  Term cells with fewer comparisons involving the object are dropped
  from its resolution, with a note.

- maxit, tol:

  Newton controls for the resolution refits.

## Value

A list of class `"rasch_btl_dif"`: `summary` (one row per object and
group term with the uniform F, adjusted p and partial eta-squared – the
term itself – the non-uniform ones – the term crossed with the opponent
band – plus `uniform_DIF`, `nonuniform_DIF` and `superseded` flags);
`terms` (the full per-object analysis-of-variance table, including its
raw and effective judge support); `levels` (resolved location, SE,
comparison count, judge count and effective judge count per object, term
and complete-design cell); `sizes` (per object, term and marginal or
interaction contrast: difference in logits, judge support for both
sides, SE, t, degrees of freedom, adjusted p, significance and practical
flags); `effects`, `factors`, `alpha`, `p_adjust`, `flag_logits`, and
`notes`. `size_family_n` records the complete planned resolved-contrast
family, including unavailable comparisons. `summary_factors` retains the
factor membership of each displayed term.

## Details

Judges are the independent units. For each object, oriented residuals
are aggregated to one weighted mean per judge and opponent band. A
split-plot analysis then tests judge factors between judges and band
effects within judges. Each factor level requires at least two judges.
Confirmatory Wald tests are available only when the base fit supplies a
valid judge-clustered covariance. The base paired-comparison calibration
must have converged. BTL-EFRM fits are not accepted: the ordinary
residual and resolution models do not contain the fitted panel and set
units.

A significant uniform term is followed by a joint refit in which the
object has one location per cell of the complete judge-factor design.
Main-effect magnitudes average these cells equally over the other
factors. Interaction magnitudes are differences between differences,
with the corresponding higher-order tensor contrast beyond two factors.
A contributing cell needs at least eight effective judges for inference;
otherwise its location and contrasts remain descriptive. Higher-order
terms supersede their component terms. Two-cell contrasts retain the
Welch reference used by the ordinary pairwise comparison. Contrasts
spanning more than two fitted cells use the effective-judge count in
their least-supported cell as a conservative denominator reference.
Models fitted with `order` retain the exposure and carry-over effects in
both the residual analysis and refit. Between-judge tests use HC3
covariance so unequal comparison workloads do not impose equal precision
on judge means. Omnibus probabilities require at least eight judges and
eight effective judges in every factor cell. Holm adjustment is the
default; `"BH"` remains available for false-discovery-rate screening. A
reported object-by-term test remains in the adjustment family when its
probability is unavailable.

Objects are resolved one at a time against the common locations of the
remaining objects. With DIF in several objects, this can induce
compensating apparent DIF in invariant objects (Andrich and Hagquist
2012, 2015). An externally anchored object is not resolved: fixing each
of its copies at the same anchor would define their difference as zero.
Anchors on the other objects are retained in the joint refit. If a
resolved-location covariance is unavailable or not positive
semidefinite, the locations and differences remain descriptive but their
uncertainty and tests are withheld. A note raised by a resolution refit
is reported against that object. Engine notes explaining a statistic
this analysis never reports – the refit's pairwise chi-square
probability and its dependence table's carry-over probability – are not
carried; the base fit retains them.

## References

Andrich, D., & Hagquist, C. (2012). Real and artificial differential
item functioning. *Journal of Educational and Behavioral Statistics*,
37(3), 387-416.

Dittrich, R., Hatzinger, R., & Katzenbeisser, W. (1998). Modelling the
effect of subject-specific covariates in paired comparison studies with
an application to university rankings. *Journal of the Royal Statistical
Society C*, 47(4), 511-525.

MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent
covariance matrix estimators with improved finite sample properties.
*Journal of Econometrics*, 29(3), 305–325.

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 100), b = rep(pr[, 2], each = 100),
                judge = sample(sprintf("J%02d", 1:20), 600, TRUE))
shift <- ifelse(d$judge %in% sprintf("J%02d", 1:10) & d$a == "C", 0.9,
         ifelse(d$judge %in% sprintf("J%02d", 1:10) & d$b == "C", -0.9, 0))
p <- plogis(beta[d$a] - beta[d$b] + shift)
d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
f <- btl(d, "a", "b", winner = "win", judge = "judge")
grp <- setNames(rep(c("g1", "g2"), each = 10), sprintf("J%02d", 1:20))
btl_dif(f, grp, objects = "C")
#> DIF for paired comparisons: 1 factor(s) [group], main effects
#> Residual ANOVA per object and term (uniform = term; non-uniform = term x opponent band)
#>  object  term F_uniform p_uniform_adj uniform_DIF F_nonuniform p_nonuniform_adj
#>       C group    15.645         0.002           *        0.142            0.858
#>  nonuniform_DIF
#>                
#> 
#> Resolved locations (logits; holm over 1 planned comparison(s); practical 0.50)
#>  object  term level_a level_b difference n_judges_a n_judges_b
#>       C group      g2      g1     -1.121         10         10
#>  effective_judges_a effective_judges_b    se      t     df   p_adj significant
#>               9.340              9.326 0.249 -4.508 16.665 < 0.001           *
#>  practical
#>          *
#> Notes: 2 object-term test(s) have 8.0--9.4 effective judges in their smallest cell; see min_effective_judges and interpret these results cautiously; C [group]: level(s) g1, g2 have 8.0--9.4 effective judges; interpret pairwise inference cautiously 
```
