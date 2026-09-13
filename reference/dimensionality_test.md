# Residual-component test of unidimensionality

Estimates each person separately on two item subsets and compares the
two estimates with a per-person t-test (Smith 2002). By default the
subsets are defined by the sign of a residual-component loading (the
first by default; any leading component may be chosen); they can also be
nominated manually (for example, by content). Under unidimensionality
and local independence the two subset estimates are independent given
the person location, so
`t = (theta_A - theta_B) / sqrt(se_A^2 + se_B^2)` is approximately
standard normal and about `alpha` of the tests should reach
significance. Persons with an extreme score on either subset are
excluded (their weighted-likelihood estimates are most biased there).
The proportion of significant tests is reported with a Clopper–Pearson
binomial confidence interval. The interval describes the observed
proportion; it is not a calibrated test of dimensionality. The test
requires a converged calibration and one response row per person.

## Usage

``` r
dimensionality_test(
  fit,
  alpha = 0.05,
  items_positive = NULL,
  items_negative = NULL,
  component = 1,
  min_score_points = 15L,
  B = 0,
  workers = 4L,
  seed = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  with one response row per person. Repeated identifiers are refused
  because the person-level comparisons and their binomial count would
  not be independent. Fully anchored scoring fits require `B = 0`;
  bootstrap refitting is not supported for these fits.

- alpha:

  Nominal significance level for the per-person t-tests.

- items_positive, items_negative:

  Optional character vectors naming the two item subsets; both must be
  given (disjoint, at least two items each), otherwise the sign of a
  residual component defines the split.

- component:

  Which residual principal component's loading sign defines the default
  split (ignored when subsets are named). Default the first component.

- min_score_points:

  Score-point threshold below which the verdict carries a caution.
  Andrich and Marais (2019) recommend subtests of roughly 15 score
  points for stable subtest estimates; shorter subsets (the norm for
  ordinary dichotomous tests) retain the analysis, with a `caution`
  field noting the reduced stability. A quiet verdict under caution is
  inconclusive, not clean: with a four-item subtest the test lacks power
  where nonparametric alternatives still flag.

- B:

  Number of parametric-bootstrap replicates that calibrate the
  proportion of significant tests under the fitted model (see Details).
  The default `0` reports the binomial interval and descriptive reading
  alone; neither split then has an inferential verdict. Each replicate
  refits the calibration, so `B = 200` costs about two hundred fits; the
  bootstrap is available for single-facet fits with a common unit whose
  thresholds were estimated directly.

- workers:

  Number of parallel workers for the bootstrap refits.

- seed:

  Optional integer seed for the bootstrap; the replicates are
  reproducible for a given seed whatever the worker count. The bootstrap
  does not support Box–Muller, including when `seed = NULL`; see
  [`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md).

## Value

A list with the proportion of significant tests, its Clopper–Pearson
confidence interval, the sample sizes (`n` used, `n_excluded_extreme`),
the item split and its source, a `multidimensional` verdict, the
corresponding uncalibrated `binomial_multidimensional` reading, a
`caution` note when the subtests fall short of `min_score_points`, and
`subset_mean_difference`, the mean and standard deviation of the
person-level differences between the two subset estimates, reported
descriptively with the `note` that no test accompanies them (see
Details). With `B > 0` the list also carries `p_boot`, the bootstrap
probability of a proportion at least as large as the observed one under
the fitted unidimensional model; `prop_null`, the mean replicate
proportion (the rate the split produces when nothing is there);
`bootstrap_resolution`, the smallest attainable bootstrap probability;
and `bootstrap`, the replicate proportions with the counts requested,
used, non-converged and failed. When the comparison itself is
unavailable (undefined split, degenerate subsets, too few persons) the
list carries a `note` explaining why and `multidimensional = NA`. Every
result carries `algorithm`, the stamp of the calculation that produced
it; a saved result without the current stamp used a superseded mean or
binomial reference. An analysis file carrying one opens with that result
dropped and a warning, and the rest of the analysis intact.
`verdict_note` explains why a descriptive comparison has no inferential
verdict.

## Details

The binomial reference assumes a per-person null rejection rate of
`alpha`. Unequal targeting and short-subset estimation bias can violate
that assumption even for a split fixed in advance. A split chosen from
the residuals is chosen to make the two subsets disagree, so its
proportion runs above `alpha` under unidimensionality. Package
simulations confirmed that applying the fixed-split binomial rule after
choosing the split from the same residuals is anti-conservative. Without
bootstrap calibration neither split therefore has a binary verdict:
`multidimensional` is `NA`, while the interval and uncalibrated binomial
reading remain available descriptively. `B > 0` supplies a model-based
reference for either split: each replicate draws responses from the
fitted model conditional on every person's raw score and missingness
pattern, refits the calibration, retains a fixed split or repeats a
residual-component split on its own residuals and recomputes the
proportion, so the bootstrap probability `p_boot` carries the same
selection the observed proportion carries. With `B > 0` the verdict is
`p_boot <= alpha`; the binomial interval is still reported, as a
description of the observed proportion rather than a test of it. A
one-sided bootstrap probability cannot be smaller than `1/(B_used + 1)`.
If that floor exceeds `alpha`, the bootstrap has no rejection region and
the verdict is withheld for either split.

The mean difference between the two subset estimates is reported but not
tested. Each subset estimate is a weighted-likelihood estimate on a
short test, and the two subsets differ in difficulty, so their
estimation bias differs systematically: under a perfectly unidimensional
Rasch model the expected difference is non-zero whenever the subsets are
not matched in targeting, and it grows relative to its standard error
with the number of persons. A t-test of that difference therefore tests
the targeting of the split rather than its dimensionality – it rejects
for every sample large enough on unidimensional data – so the difference
is reported as a description of the split and the inference is withheld.
The person-level comparisons also carry this bias; their proportion
needs the bootstrap reference before it can support a dimensionality
verdict.

## References

Smith, E. V. Jr. (2002). Detecting and evaluating the impact of
multidimensionality using item fit statistics and principal component
analysis of residuals. Journal of Applied Measurement, 3(2), 205–231.

Tennant, A., & Pallant, J. F. (2006). Unidimensionality matters! (A tale
of two Smiths?). Rasch Measurement Transactions, 20(1), 1048–1051.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 8)
X <- matrix(rbinom(300 * 8, 1, plogis(outer(rnorm(300), d, "-"))), 300, 8)
colnames(X) <- paste0("I", 1:8)
dimensionality_test(
  rasch(X), items_positive = paste0("I", 1:4),
  items_negative = paste0("I", 5:8))$multidimensional
#> [1] NA
if (FALSE) { # \dontrun{
# A longer run calibrates the data-driven split under the fitted model.
dimensionality_test(rasch(X), B = 99, workers = 1, seed = 1)$p_boot
} # }
```
