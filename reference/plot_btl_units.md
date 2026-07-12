# Plot the frame units of a paired-comparison EFRM fit

Caterpillar plot of the estimated units on the log scale: one row per
panel unit `phi_g` and one per set unit `alpha_s`, with 95 per cent
intervals, the reference (unit one) marked, mirroring
[`plot_frames`](https://drjoshmcgrane.github.io/rasch/reference/plot_frames.md)
in the package's house style.

## Usage

``` r
plot_btl_units(fit)
```

## Arguments

- fit:

  A fitted object from
  [`btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md).

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
# \donttest{
# see ?btl_efrm for a complete simulated example
# }
```
