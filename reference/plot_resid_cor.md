# Plot the residual-correlation heatmap

Cells are coloured by Q3\* – each pair's residual correlation minus the
average off-diagonal correlation – so white marks the value expected
under local independence and warm colour marks dependence. The scale
saturates at a Q3\* of `cap` rather than the +/-1 of an ordinary
correlation: a residual correlation seldom reaches even 0.5 under a
fitting model, and the conventional flag (Q3\* above 0.2; Christensen,
Makransky and Horton 2017, marked on the key) sits well within the
range, so the colour is spent where the values actually discriminate.

## Usage

``` r
plot_resid_cor(fit, cap = 0.5, flag = 0.2)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- cap:

  Q3\* value at which the colour saturates (default 0.5).

- flag:

  Q3\* value marked on the key as the dependence flag (default 0.2).

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
plot_resid_cor(rasch(X))
```
