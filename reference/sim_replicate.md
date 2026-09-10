# Replicate a simulation for Monte Carlo studies

Calls one of the `simulate_*` functions `n` times with successive seeds,
returning the datasets as a list – for power, Type-I, or
parameter-recovery studies.

## Usage

``` r
sim_replicate(FUN, n, ..., seed = NULL)
```

## Arguments

- FUN:

  A simulator, e.g.
  [`simulate_rasch`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md).

- n:

  Number of datasets.

- ...:

  Arguments passed to `FUN` (the same each replicate).

- seed:

  Seed of the first replicate (each subsequent one increments it). See
  [`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md)
  for generator support.

## Value

A list of class `"rasch_sim_batch"`, one simulated dataset per element.

## Examples

``` r
# 8 datasets with a planted DIF item; how often is it flagged?
batch <- sim_replicate(simulate_rasch, 4, n_persons = 300, n_items = 8,
                       dif = list(items = "I05", uniform = 0.8), n_groups = 2,
                       seed = 1)
# sim_apply() is resilient: a replicate the estimator refuses (e.g. a
# small or disconnected draw) contributes NA instead of aborting the run
flagged <- sim_apply(batch, function(d)
  dif_anova(rasch(d, id = "id", factors = "group"))$summary$uniform_DIF[5])
mean(flagged, na.rm = TRUE)
#> [1] 0.5
```
