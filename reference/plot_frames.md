# Plot frame units

Caterpillar plot of the frame units `rho_sg = alpha_s phi_g` on the log
scale, grouped by item set and coloured by person group, with 95 per
cent error bars; frames with pooled fit residuals beyond the band are
highlighted.

## Usage

``` r
plot_frames(fit, band = 2.5)
```

## Arguments

- fit:

  A fitted object from
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

- band:

  Pooled fit residual band beyond which a frame is highlighted.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
# \donttest{
# see ?rasch_efrm for a complete simulated example
# }
```
