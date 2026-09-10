# Simulate person-by-item Rasch data

Generates dichotomous, partial credit, or rating scale data. Optional
arguments introduce item misfit, guessing, multidimensionality, local
dependence, DIF, response styles, or missingness. Generating values are
stored in `attr(x, "truth")`. A positive planted proportion selects at
least one person or response cell when the requested departures can
coexist.

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

- theta_mean, theta_sd:

  Realised sample mean and standard deviation of the person
  distribution.

- theta_dist:

  Shape of the person distribution: `"normal"`, `"uniform"`, `"skew"`,
  or `"bimodal"`.

- difficulty:

  Either the two endpoints of an evenly spaced location range, or one
  location per item.

- threshold_spread:

  Half-range of the category thresholds about each item location
  (polytomous).

- discrimination:

  The item slope, supplied as one value or one per item. One common
  value changes the logit unit without departing from a Rasch model.
  Differences between items are a deliberate slope departure: values
  above the common pattern produce steeper responses and values below it
  produce flatter responses, and require `theta_sd > 0`.

- guessing:

  Scalar or length-`n_items` lower asymptote (dichotomous): low-location
  persons answer correctly by chance. A positive asymptote requires
  `theta_sd > 0`.

- second_dim:

  `NULL`, or `list(items=, rho=)`: the named items load on a second
  trait whose realised sample correlation with the first is `rho`. It
  lies in \[-1, 1); a correlation of 1 would reproduce the primary trait
  rather than plant another dimension. At least three persons are needed
  unless `rho` is -1. Each item is named once, and at least one item
  must remain on the primary trait; moving every item to the second
  trait is still a one-dimensional scale.

- dependence:

  `NULL`, or `list(pairs=, strength=)`: each pair's second item responds
  partly to the first. This departure feeds the residual-dependence
  diagnostics. Each directed pair is listed once.

- dif:

  `NULL`, or `list(items=, uniform=, nonuniform=)`: the named items
  function differently for the last person group: a location shift
  (`uniform`) and/or a slope change (`nonuniform`). Needs
  `n_groups >= 2`; each item is named once and at least one item must
  remain invariant. A common shift or slope change on every item changes
  the group origin or unit rather than defining item DIF. A non-uniform
  effect also requires `theta_sd > 0`.

- careless:

  Proportion of persons who answer at random. At least one person in
  each generated group must retain model-based responses. Careless and
  response-style assignments are disjoint; their requested counts must
  fit.

- response_style:

  `NULL`, or `list(type=, prop=, strength=)` with `type` `"extreme"` or
  `"middle"`: a proportion `prop` of persons favour the end (or middle)
  categories regardless of the trait, with distortion `strength`
  (default 1.6) on the log-probability scale (polytomous). At least one
  person must remain without the style; applying one category weighting
  to everyone merely changes the fitted thresholds.

- speeded:

  Proportion not reached at the last item: a growing tail of missing
  responses over the final items. These cells are kept distinct from any
  completely-at-random missing cells. Persons are selected independently
  of their ability and responses; this does not simulate non-ignorable
  missingness or change the response model.

- disordered:

  `NULL` or item names/indices given disordered thresholds (polytomous;
  feeds the threshold diagnostics).

- n_groups:

  Number of equal person groups (a `group` factor column is added when
  \> 1, for DIF).

- missing:

  Proportion of responses set missing completely at random, drawn from
  cells not already missing through speededness. The requested count
  must fit among those cells and leave at least one observed response.

- seed:

  Optional non-negative whole-number RNG seed. See
  [`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md)
  for generator support.

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
#>  item infit_ms outfit_ms
#>   I01    1.013     0.428
#>   I02    1.043     0.955
#>   I03    0.982     0.822
#>   I04    1.054     1.124
#>   I05    1.040     0.946
#>   I06    1.091     1.134
#>   I07    1.027     0.993
#>   I08    1.113     1.126
#>   I09    1.045     1.052
#>   I10    1.046     0.949
#>   I11    1.036     0.984
#>   I12    0.985     0.781
dif_anova(fit)$summary                           # item 6 flags
#>  item  term F_uniform p_uniform p_uniform_adj eta2_uniform uniform_DIF
#>   I01 group     4.320     0.038         0.805        0.012            
#>   I02 group     1.858     0.174         1.000        0.005            
#>   I03 group     0.033     0.856         1.000        0.000            
#>   I04 group     2.029     0.155         1.000        0.005            
#>   I05 group     0.856     0.356         1.000        0.002            
#>   I06 group    25.647   < 0.001       < 0.001        0.062           *
#>   I07 group     3.186     0.075         1.000        0.008            
#>   I08 group     0.083     0.773         1.000        0.000            
#>   I09 group     0.754     0.386         1.000        0.002            
#>   I10 group     0.520     0.471         1.000        0.001            
#>   I11 group     4.891     0.028         0.620        0.012            
#>   I12 group     0.018     0.895         1.000        0.000            
#>  F_nonuniform p_nonuniform p_nonuniform_adj eta2_nonuniform nonuniform_DIF
#>         1.550        0.173            1.000           0.020               
#>         0.020        1.000            1.000           0.000               
#>         1.761        0.120            1.000           0.022               
#>         1.711        0.131            1.000           0.022               
#>         1.192        0.313            1.000           0.015               
#>         0.473        0.797            1.000           0.006               
#>         2.561        0.027            0.620           0.032               
#>         0.954        0.446            1.000           0.012               
#>         0.822        0.534            1.000           0.011               
#>         0.512        0.767            1.000           0.007               
#>         0.455        0.810            1.000           0.006               
#>         1.208        0.305            1.000           0.015               
#>  superseded
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
#>            
```
