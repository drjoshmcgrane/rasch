# Bootstrap fit statistics

Refers fit statistics to replicated datasets fitted in the same way as
the observed data. For a person-by-item Rasch model, the default
generator conditions on each person's observed raw score and missingness
pattern. The person parameter then cancels by sufficiency. Item
parameters and person locations are re-estimated in every replicate. The
original class-interval rule is repeated, including automatic per-item
allocation with missing responses. An automatic count is resolved
against the observed fit and held across replicates, so every replicate
chi-square is scored on the intervals the observed one used. Tied
locations remain in one interval. A result stored by an earlier version
under `theta = "resample"`, `"fixed"` or `"normal"` was built from a
null that mixed interval counts, so it is refused and must be
recomputed; results from the default generator are unchanged. The
generator assumes independent response rows. A fit with repeated person
IDs is therefore refused because this bootstrap does not reproduce
within-person dependence.

## Usage

``` r
fit_bootstrap(
  fit,
  B = 200,
  theta = c("conditional", "resample", "fixed", "normal"),
  workers = 4L,
  seed = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md) or
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).
  Extended-frame, many-facet and explanatory person-by-item fits are not
  supported. Explanatory paired-comparison fits are supported. Fully
  anchored scoring fits are not supported by this refitting procedure.
  An older anchored fit must retain its original anchor settings.

- B:

  Positive whole number of bootstrap replicates.

- theta:

  Generator for a person-by-item fit. `"conditional"` retains each
  observed raw score and missingness pattern. `"resample"` resamples
  estimated person locations, `"fixed"` reuses each location with the
  same response row and missingness pattern, and `"normal"` draws from a
  normal distribution with error-corrected variance. Fixed generation
  requires a finite location for each row with an observed response.
  This argument does not apply to paired comparisons.

- workers:

  Number of parallel bootstrap workers. The default is four, reduced
  when fewer physical cores are available or the R process has a lower
  system limit. Per-replicate seeds and random-number generator settings
  are shared, so results do not depend on the worker count.

- seed:

  Optional non-negative whole-number seed within the integer range. The
  caller's random stream is restored on exit. This calculation does not
  support Box–Muller, including when `seed = NULL`; see
  [`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md).

## Value

An object of class `rasch_fit_bootstrap`. For a person-by-item fit, it
contains `items`, `persons`, `total`, `replicates`, adjustment metadata
and replicate counts, including separate non-convergence and
other-failure counts. For a paired-comparison fit, the corresponding
tables are `pairs`, `objects`, `judges` and `total`. Runs requesting at
least 30 replicates are withheld unless at least 30 and 90 percent of
the requested refits are usable.

## Details

Item chi-squares use the upper tail. Fit residuals, infit and outfit use
equal-tailed probabilities. Holm adjustment is applied separately to
each predeclared item-statistic family; an unavailable item remains in
the multiplicity count. Under the conditional generator, the same
replicates also give person-specific null distributions. Person
probabilities are adjusted with a single-step maximum-statistic
distribution across persons for each statistic. A maximum-statistic
adjustment is withheld for the complete family if any testable member
lacks a usable joint null.

The adjustments describe the fitted global null. They do not guarantee
familywise error control among otherwise fitting items, persons, objects
or judges when another member departs from the model. Each fit statistic
is a separate family. The marginal probabilities are retained.

For a [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
fit, outcomes are generated on the fitted comparison design and the
model is refitted. The result covers the total pairwise chi-square and
pair, object and judge fit. Ordered response thresholds, judge
allocation, counts, anchors, position effects and explanatory object
restrictions are retained. History-dependent effects are generated in
sequence. Half-weighted ties and frame-dependent paired-comparison fits
are refused because the present generator does not reproduce those
models.

Bootstrap probabilities use \\(1+r)/(1+B)\\. Their one-sided resolution
is therefore \\1/(1+B)\\ and their equal-tailed resolution is
\\2/(1+B)\\. For \\L\\ items, Holm-adjusted inference at .05 requires at
least \\20L\\ replicates for the chi-square and \\40L\\ for the
two-sided statistics. Runs of at least 30 replicates also require at
least 30 and 90 percent of the requested refits to be usable. Otherwise
inference is withheld because failed refits can select a non-random
subset of the bootstrap null.

The approach follows Wolfe (2013), who bootstrapped critical values for
Rasch fit statistics by generating from the estimated person and item
parameters and averaging replicate quantiles into cut points. This
implementation conditions on the observed scores instead — generating at
estimated locations carries their estimation error into the null, which
his single 1,000-person design could not surface — and returns
per-statistic probabilities under a declared multiplicity policy rather
than averaged cut points.

The need for a replicated null is not a real-data artefact. Wu and Adams
(2013) derive the mean squares' null variance as about \\2/N\\ but
conclude the standardised statistics have a sample-size-independent
null, reading contrary reports as flawed simulation or as genuine misfit
surfacing in large samples. Under data generated from the model —
neither explanation applies — infit z beyond 1.96 still flags 11.5
percent of correctly fitting items at 250 persons and 68.4 percent at
4,000 with eight items: parameter estimation alone moves the null, and
reproducing that estimation in every replicate is what calibrates it.

## References

Andrich, D. and Marais, I. (2019) *A Course in Rasch Measurement
Theory*. Springer.

Wolfe, E. W. (2013). A bootstrap approach to evaluating person and item
fit to the Rasch model. *Journal of Applied Measurement*, 14(1), 1–9.

Wu, M. and Adams, R. J. (2013). Properties of Rasch residual fit
statistics. *Journal of Applied Measurement*, 14(4), 339–355.

Molenaar, I. W. and Hoijtink, H. (1996). Person-fit test statistics for
the Rasch model. *Applied Measurement in Education*, 9(1), 87–106.

Westfall, P. H. and Young, S. S. (1993). *Resampling-Based Multiple
Testing*. Wiley.

## See also

[`chisq_detail`](https://drjoshmcgrane.github.io/rasch/reference/chisq_detail.md),
[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md),
[`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md).

## Examples

``` r
set.seed(1)
d <- seq(-1.5, 1.5, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
# an exploratory run, kept small for speed: raw probabilities are usable,
# and the warning says what Holm-adjusted flagging at .05 would need
bs <- suppressWarnings(fit_bootstrap(rasch(X), B = 49, seed = 1))
bs$items[c("item", "chisq", "chisq_p_boot", "fit_resid", "fit_resid_p_boot")]
#>  item chisq chisq_p_boot fit_resid fit_resid_p_boot
#>    I1 2.839        0.920    -0.090            0.800
#>    I2 3.867        0.680     0.679            0.640
#>    I3 5.085        0.660     0.313            0.760
#>    I4 6.088        0.500     1.172            0.520
#>    I5 3.258        0.820     0.023            0.800
#>    I6 7.682        0.160    -0.200            0.920
```
