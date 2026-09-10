# Residual map of the leading paired-comparison bimension

Objects placed in the leading bimension plane to show the pattern of
residual comparisons. The coordinates describe the pattern, not its
magnitude or significance; use
[`plot_btl_scree`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_scree.md)
to compare its strength with the conditional reference. Point size grows
with the object's location on the primary scale.

## Usage

``` r
plot_btl_dim_map(x, ...)
```

## Arguments

- x:

  A `"rasch_btl_dim"` object.

- ...:

  Unused.

## Value

Called for its plotting side effect.

## Examples

``` r
# \donttest{
d <- simulate_btl(7, 12, reps_per_pair = 20, seed = 1)
fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
dimensions <- btl_dimensionality(fit, reps = 20)
plot_btl_dim_map(dimensions)

# }
```
