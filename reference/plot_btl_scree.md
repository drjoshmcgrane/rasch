# Scree of paired-comparison residual bimensions

Bimension strengths against the model-simulated noise reference (its
mean and finite-simulation 5 available. A leading bar clearing the band
suggests residual structure beyond the fitted model under that
reference.

## Usage

``` r
plot_btl_scree(x, ...)
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
plot_btl_scree(dimensions)

# }
```
