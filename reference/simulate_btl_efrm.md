# Simulate paired-comparison EFRM data with differing frame units

Generates dichotomous paired comparisons whose latent unit differs
across judge-panel by object-set frames – the paired-comparison
extension of the extended frame of reference model (Humphry 2005) fitted
by
[`btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md).
Objects in set `s` have a within-set calibration location `beta`; their
common-scale value is `v = alpha_s beta + kappa_s`. A comparison judged
in panel `g` carries the panel unit `phi_g`: within a set the comparison
logit is `phi_g (beta_a - beta_b)`, across sets it is
`phi_g (v_a - v_b)`. The planted panel units, set units and origins are
recovered by
[`btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md).

## Usage

``` r
simulate_btl_efrm(
  n_objects_per_set = 8,
  n_sets = 2,
  n_judges_per_panel = 6,
  n_panels = 2,
  reps_within = 20,
  reps_cross = 20,
  panel_units = NULL,
  set_units = NULL,
  set_origins = NULL,
  object_sd = 1,
  seed = NULL
)
```

## Arguments

- n_objects_per_set, n_sets:

  Objects in each set and number of sets.

- n_judges_per_panel, n_panels:

  Judges in each panel and number of panels.

- reps_within:

  Replications of each within-set object pair.

- reps_cross:

  Replications of each cross-set object pair.

- panel_units:

  Panel units `phi` (length `n_panels`); the default is all one, and any
  supplied vector is rescaled to geometric mean one.

- set_units:

  Set units `alpha` (length `n_sets`); the default is all one, and
  `alpha_1` is forced to one (the reference set).

- set_origins:

  Set origins `kappa` (length `n_sets`); the default is all zero, and
  `kappa_1` is forced to zero.

- object_sd:

  Spread of the within-set calibration locations.

- seed:

  Optional RNG seed.

## Value

A data frame of class `"rasch_sim"` with columns `object_a`, `object_b`,
`winner`, `judge` and `panel`, and `attr(x, "truth")` holding the
common-scale values `v`, the per-set `beta`, the units `phi`, `alpha`,
`kappa`, and the `object_sets` map to pass to
[`btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md).

## Examples

``` r
d <- simulate_btl_efrm(6, 2, set_units = c(1, 1.4), seed = 1)
bt <- btl_efrm(d, "object_a", "object_b", winner = "winner",
               judge = "judge", panels = "panel",
               object_sets = attr(d, "truth")$object_sets)
bt$alpha_table   # recovers the ~1.4 set unit
#>    set    alpha se_log_alpha       z          p
#> 1 set1 1.000000           NA      NA         NA
#> 2 set2 1.457189    0.1818235 2.07074 0.03838314
```
