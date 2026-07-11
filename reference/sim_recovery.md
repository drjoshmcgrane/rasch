# Parameter recovery of a fit against the simulation truth

Compares the parameters recovered by a fit with the ones a `simulate_*`
function planted (carried on the data as `attr(sim, "truth")`): item
difficulties and person abilities for a Rasch fit, object locations for
a paired-comparison fit, rater severities (with item and person
measures) for a many-facet fit, and the set units for a frames fit.
Locations are mean-centred before comparison, since the model identifies
them only up to an origin.

## Usage

``` r
sim_recovery(fit, sim)
```

## Arguments

- fit:

  A fit of the simulated data
  ([`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md),
  [`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
  or
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)).

- sim:

  The simulated data (from a `simulate_*` function).

## Value

A list of class `"rasch_recovery"`: `summary` (per parameter type: n,
correlation, RMSE, bias) and `pieces` (the true and estimated values
behind each).

## Examples

``` r
d <- simulate_rasch(500, 12, seed = 1)
sim_recovery(rasch(d), d)$summary
#>         parameter   n correlation      rmse         bias
#> 1 item difficulty  12   0.9980002 0.1703475 6.822709e-17
#> 2  person ability 500   0.8177074 0.6756645 5.300462e-03
```
