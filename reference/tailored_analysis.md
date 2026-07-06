# Tailored analysis for guessing

Runs the four-step tailored procedure of Andrich, Marais and Humphry
(2012) on a dichotomous analysis. Step 1 is the supplied fit. Step 2
(tailored) sets to missing every observed response whose modelled
probability of success, at the step-1 person and item estimates, is
below `chance`, and re-estimates items and persons. Step 3
(origin-equated) re-analyses the *original* data with the mean location
of the anchor items fixed at their tailored values by average anchoring,
so the two calibrations share an origin. Step 4 (all-anchored) fixes
every item at its tailored difficulty and re-estimates persons on the
original data. Guessing is indicated when difficult items are estimated
harder in the tailored analysis than in the origin-equated one; the
comparison table and
[`plot_equate`](https://drjoshmcgrane.github.io/rmt/reference/plot_equate.md)
on the two calibrations show it directly.

## Usage

``` r
tailored_analysis(fit, chance = 0.25, anchor_items = NULL)
```

## Arguments

- fit:

  A dichotomous fit from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- chance:

  The guessing floor: the probability of success by chance (1/number of
  options; default 0.25).

- anchor_items:

  Items whose mean location fixes the common origin in step 3. The
  default takes the third of the test (at least two items) least
  affected by tailoring – fewest responses removed, ties broken towards
  the easier tailored location – which are the easy items the procedure
  trusts.

## Value

A list of class `"rmt_tailored"`: `tailored`, `origin_equated`, and
`anchored` fits, the comparison `table` (initial, tailored,
origin-equated locations, the tailored-minus-equated `shift`, and its
`z`), the number of responses removed, and the anchor items used.

## References

Waller, M. I. (1989). Modeling guessing behavior: A comparison of two
IRT models. Applied Psychological Measurement, 13, 233-243. Andrich, D.,
Marais, I. and Humphry, S. (2012). Using a theorem by Andersen and the
dichotomous Rasch model to assess the presence of random guessing in
multiple choice items. Journal of Educational and Behavioral Statistics,
37, 417-442.

## Examples

``` r
set.seed(1); N <- 800
d <- seq(-2, 2.5, length.out = 10); th <- rnorm(N)
P <- plogis(outer(th, d, "-"))
P <- 0.25 + 0.75 * P            # uniform guessing floor
X <- matrix(rbinom(N * 10, 1, P), N, 10)
colnames(X) <- paste0("I", 1:10)
ta <- tailored_analysis(rasch(X), chance = 0.25)
ta$table
#>     item    initial    tailored origin_equated removed      shift         z
#> I1    I1 -1.6725135 -1.77025037    -1.77025037       0 0.00000000 0.0000000
#> I2    I2 -1.2209361 -1.35953938    -1.35953938       0 0.00000000 0.0000000
#> I3    I3 -0.7506180 -0.82611031    -0.82611031       7 0.00000000 0.0000000
#> I4    I4 -0.3397033 -0.40510059    -0.40510059      33 0.00000000 0.0000000
#> I5    I5  0.1120400  0.04827271     0.02317712      33 0.02509559 0.2161692
#> I6    I6  0.2715352  0.26473095     0.18275306      90 0.08197789 0.7059928
#> I7    I7  0.5422752  0.51272210     0.45362192      90 0.05910019 0.5167861
#> I8    I8  0.9005334  0.88296336     0.81230152     206 0.07066184 0.5786044
#> I9    I9  1.0370971  1.28376521     0.94890015     206 0.33486506 2.6189468
#> I10  I10  1.1202901  1.36854631     1.03204123     206 0.33650508 2.5838021
```
