# Plot a paired-comparison equating comparison

Scatter of the two calibrations' common-object locations with the
shifted identity line, per-object 95 per cent error bars, and a dotted
guide band at the average pooled precision; objects that drift (after
the multiplicity adjustment) are highlighted and labelled. The
counterpart of
[`plot_equate`](https://drjoshmcgrane.github.io/rasch/reference/plot_equate.md)
for Bradley-Terry-Luce scales.

## Usage

``` r
plot_btl_equate(fit1, fit2, ...)
```

## Arguments

- fit1:

  A fitted object from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- fit2:

  A second
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) fit,
  or a bank data frame with columns `object`, `location`, and optionally
  `se`.

- ...:

  Passed to
  [`btl_equate`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)
  (e.g. `alpha`, `p_adjust`).

## Value

Called for its plotting side effect; invisibly the
[`btl_equate`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)
result.

## Examples

``` r
set.seed(1)
beta <- setNames(seq(-2, 2, length.out = 8), paste0("O", 1:8))
sim <- function(objs) {
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  btl(d, "a", "b", "win")
}
plot_btl_equate(sim(paste0("O", 1:7)), sim(paste0("O", 2:8)))
```
