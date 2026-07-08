# Plot an item's characteristic curves across frames

The signature display of the extended frame of reference model: the
model expected-score curve of one underlying item drawn once per person
group (curves fan with the group units), with observed class-interval
means per group overlaid.

## Usage

``` r
plot_icc_frames(fit, item, n_groups = fit$n_groups, grid = seq(-5, 5, 0.05))
```

## Arguments

- fit:

  A fitted object from
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

- item:

  Underlying item name.

- n_groups:

  Number of class intervals for the observed means.

- grid:

  Logit grid.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
# \donttest{
# see ?rasch_efrm for a complete simulated example
# }
```
