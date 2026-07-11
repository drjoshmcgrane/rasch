# Simulate paired-comparison (BTL) data with dial-in misfit

Generates dichotomous or graded paired comparisons from the
Bradley-Terry-Luce model, with optional departures each of which a
paired-comparison diagnostic is built to detect. The result is a data
frame ready for
[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md), with
the truth attached.

## Usage

``` r
simulate_btl(
  n_objects = 8,
  n_judges = 12,
  reps_per_pair = 25,
  model = c("dichotomous", "graded"),
  n_categories = 4,
  object_sd = 1,
  second_attribute = NULL,
  erratic_judges = 0,
  dependence = NULL,
  seed = NULL
)
```

## Arguments

- n_objects, n_judges:

  Objects to scale and judges comparing them.

- reps_per_pair:

  Comparisons made of each object pair.

- model:

  `"dichotomous"` (a winner) or `"graded"` (a rated margin in
  `n_categories` categories).

- n_categories:

  Categories for the graded model.

- object_sd:

  Spread of the object locations (evenly spaced, sum-zero).

- second_attribute:

  `NULL`, or `list(rho=)`: half the judges rank by a second object
  attribute correlated `rho` with the first – genuine
  multidimensionality. Feeds
  [`btl_dimensionality`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md)
  and
  [`btl_transitivity`](https://drjoshmcgrane.github.io/rasch/reference/btl_transitivity.md).

- erratic_judges:

  Proportion of judges who choose at random. Feeds the judge fit
  residual,
  [`btl_transitivity`](https://drjoshmcgrane.github.io/rasch/reference/btl_transitivity.md)
  consistency, and
  [`judge_surprise`](https://drjoshmcgrane.github.io/rasch/reference/judge_surprise.md).

- dependence:

  `NULL`, or `list(exposure=, carry_over=)`: within-judge order effects
  (a seen-before advantage and a pull from the judge's own earlier
  verdicts). Adds an `order` column. Feeds the dependence effects of
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- seed:

  Optional RNG seed.

## Value

A data frame of class `"rasch_sim"`: `object_a`, `object_b`, `winner`
(or `response` when graded), `judge`, and `order` when dependence is
planted; with `attr(x, "truth")`.

## Examples

``` r
d <- simulate_btl(8, 12, erratic_judges = 0.15, seed = 1)
bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
bt$judges          # the erratic judges carry large fit residuals
#>    judge  n  infit_ms outfit_ms  fit_resid df_fit
#> 1     J1 64 1.2295559 1.2860809  1.4008368  63.36
#> 2    J10 61 0.8282835 0.7923607 -1.3444319  60.39
#> 3    J11 41 0.7521572 0.6528004 -1.8596500  40.59
#> 4    J12 62 0.9105348 0.8766737 -0.7595074  61.38
#> 5     J2 50 1.4382514 1.8513748  3.0502717  49.50
#> 6     J3 59 0.9333531 0.8985656 -0.5650125  58.41
#> 7     J4 65 1.0931577 1.1947629  1.1791819  64.35
#> 8     J5 70 1.0077925 1.0865689  0.5635480  69.30
#> 9     J6 60 0.9139414 0.8926580 -0.7310213  59.40
#> 10    J7 54 1.0475564 1.0148149  0.0828951  53.46
#> 11    J8 53 1.0434358 1.0532995  0.2926053  52.47
#> 12    J9 61 0.8308331 0.7986060 -1.4223486  60.39
```
