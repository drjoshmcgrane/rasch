# Simulate person-by-item Rasch data with dial-in misfit

Generates dichotomous or polytomous (partial credit / rating scale) data
from the Rasch model, with optional, individually controllable
departures from it – each of which the package's matching diagnostic is
built to detect. The result is a data frame ready for
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
with the true parameters attached as `attr(x, "truth")`.

## Usage

``` r
simulate_rasch(
  n_persons = 500,
  n_items = 20,
  model = c("dichotomous", "PCM", "RSM"),
  n_categories = 3,
  theta_mean = 0,
  theta_sd = 1,
  theta_dist = "normal",
  difficulty = c(-2.5, 2.5),
  threshold_spread = 1.2,
  discrimination = 1,
  guessing = 0,
  second_dim = NULL,
  dependence = NULL,
  dif = NULL,
  careless = 0,
  response_style = NULL,
  speeded = 0,
  disordered = NULL,
  n_groups = 1,
  missing = 0,
  seed = NULL
)
```

## Arguments

- n_persons, n_items:

  Sample size and test length.

- model:

  `"dichotomous"`, `"PCM"`, or `"RSM"`. Under `"RSM"` every item shares
  one category-threshold pattern (items differ by location only); under
  `"PCM"` each item's threshold spacings and span are drawn afresh, as
  the partial credit model allows.

- n_categories:

  Response categories for polytomous models (\>= 3).

- theta_mean, theta_sd, theta_dist:

  Person distribution: mean, SD, and shape (`"normal"`, `"uniform"`,
  `"skew"`, `"bimodal"`).

- difficulty:

  Two numbers giving the item-location range (evenly spaced), or a
  length-`n_items` vector of locations.

- threshold_spread:

  Half-range of the category thresholds about each item location
  (polytomous).

- discrimination:

  Scalar or length-`n_items`: the slope of each item. Values above 1
  over-discriminate (Guttman-like, negative fit residual); below 1
  under-discriminate (noisy, positive residual). Feeds infit/outfit and
  the item-fit F.

- guessing:

  Scalar or length-`n_items` lower asymptote (dichotomous): low-ability
  persons answer correctly by chance. Feeds
  [`tailored_analysis`](https://drjoshmcgrane.github.io/rasch/reference/tailored_analysis.md).

- second_dim:

  `NULL`, or `list(items=, rho=)`: the named items load on a second
  trait correlated `rho` with the first. Feeds
  [`dimensionality_test`](https://drjoshmcgrane.github.io/rasch/reference/dimensionality_test.md).

- dependence:

  `NULL`, or `list(pairs=, strength=)`: each pair's second item responds
  partly to the first (response dependence). Feeds
  [`residual_correlations`](https://drjoshmcgrane.github.io/rasch/reference/residual_correlations.md)
  /
  [`dependence_magnitude`](https://drjoshmcgrane.github.io/rasch/reference/dependence_magnitude.md).

- dif:

  `NULL`, or `list(items=, uniform=, nonuniform=)`: the named items
  function differently for the last person group – a location shift
  (`uniform`) and/or a slope change (`nonuniform`). Needs
  `n_groups >= 2`. Feeds
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  /
  [`dif_size`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md).

- careless:

  Proportion of persons who answer at random (person misfit; feeds
  person infit/outfit).

- response_style:

  `NULL`, or `list(type=, prop=, strength=)` with `type` `"extreme"` or
  `"middle"`: a proportion `prop` of persons favour the end (or middle)
  categories regardless of the trait, with distortion `strength`
  (default 1.6) on the log-probability scale (polytomous; feeds the
  category diagnostics and person fit).

- speeded:

  Proportion not-reached at the last item: a growing tail of missing
  responses over the final items, as under time pressure (feeds the item
  statistics and the missingness pattern).

- disordered:

  `NULL` or item names/indices given disordered thresholds (polytomous;
  feeds the threshold diagnostics).

- n_groups:

  Number of equal person groups (a `group` factor column is added when
  \> 1, for DIF).

- missing:

  Proportion of responses set missing (completely at random).

- seed:

  Optional RNG seed.

## Value

A data frame of class `"rasch_sim"` (item columns `I01`..., an `id`
column, and a `group` column when grouped), with `attr(x, "truth")`
holding the generating parameters and the planted departures.

## Examples

``` r
# a clean scale with one over-discriminating item and one DIF item
d <- simulate_rasch(400, 12, discrimination = c(3, rep(1, 11)),
                    dif = list(items = "I06", uniform = 1), n_groups = 2,
                    seed = 1)
fit <- rasch(d, id = "id", factors = "group")
fit$items[c("item", "infit_ms", "outfit_ms")]   # item 1 misfits
#>    item  infit_ms outfit_ms
#> 1   I01 1.0129467 0.4265500
#> 2   I02 1.0434044 0.9527288
#> 3   I03 0.9820669 0.8193981
#> 4   I04 1.0541616 1.1210207
#> 5   I05 1.0402727 0.9429850
#> 6   I06 1.0909707 1.1305319
#> 7   I07 1.0264347 0.9907582
#> 8   I08 1.1123775 1.1226994
#> 9   I09 1.0439933 1.0497987
#> 10  I10 1.0452768 0.9465806
#> 11  I11 1.0341611 0.9823771
#> 12  I12 0.9830431 0.7799822
dif_anova(fit)$summary                           # item 6 flags
#>    item  term   F_uniform    p_uniform p_uniform_adj eta2_uniform uniform_DIF
#> 1   I01 group  4.83171573 2.853237e-02  1.151575e-01 1.233110e-02       FALSE
#> 2   I02 group  1.83782166 1.759971e-01  3.519941e-01 4.726448e-03       FALSE
#> 3   I03 group  0.03257713 8.568613e-01  8.930004e-01 8.417155e-05       FALSE
#> 4   I04 group  2.03757803 1.542602e-01  3.519941e-01 5.237484e-03       FALSE
#> 5   I05 group  0.86937693 3.517092e-01  5.935165e-01 2.241417e-03       FALSE
#> 6   I06 group 25.47500950 6.910514e-07  8.292617e-06 6.176134e-02        TRUE
#> 7   I07 group  3.26266938 7.165039e-02  2.149512e-01 8.360188e-03       FALSE
#> 8   I08 group  0.08509096 7.706696e-01  8.930004e-01 2.198249e-04       FALSE
#> 9   I09 group  0.72302977 3.956776e-01  5.935165e-01 1.864810e-03       FALSE
#> 10  I10 group  0.50867110 4.761445e-01  6.348594e-01 1.312670e-03       FALSE
#> 11  I11 group  4.81605772 2.878939e-02  1.151575e-01 1.229163e-02       FALSE
#> 12  I12 group  0.01811637 8.930004e-01  8.930004e-01 4.681013e-05       FALSE
#>    F_nonuniform p_nonuniform p_nonuniform_adj eta2_nonuniform nonuniform_DIF
#> 1    1.54966116   0.17347796        0.5204339    0.0196284713          FALSE
#> 2    0.01997421   0.99983670        0.9998367    0.0002579982          FALSE
#> 3    1.76053942   0.11999503        0.5204339    0.0222401139          FALSE
#> 4    1.71079174   0.13105466        0.5204339    0.0216252638          FALSE
#> 5    1.19165870   0.31252139        0.6250428    0.0151626613          FALSE
#> 6    0.47260534   0.79668141        0.8834756    0.0060689550          FALSE
#> 7    2.56058359   0.02694998        0.3233997    0.0320230728          FALSE
#> 8    0.95357514   0.44624178        0.7649859    0.0121701548          FALSE
#> 9    0.82209172   0.53445762        0.8016864    0.0105097128          FALSE
#> 10   0.51196042   0.76724495        0.8834756    0.0065710119          FALSE
#> 11   0.45468459   0.80985262        0.8834756    0.0058401700          FALSE
#> 12   1.20767050   0.30477298        0.6250428    0.0153632653          FALSE
#>    superseded
#> 1       FALSE
#> 2       FALSE
#> 3       FALSE
#> 4       FALSE
#> 5       FALSE
#> 6       FALSE
#> 7       FALSE
#> 8       FALSE
#> 9       FALSE
#> 10      FALSE
#> 11      FALSE
#> 12      FALSE
```
