# Targeting plot for a paired-comparison design

The paired-comparison counterpart of a test-information display. Every
object is a dot at its location (x) and its design information (y, the
pooled Fisher information of the comparisons it took part in), the dot
sized by how many comparisons that is. A reference curve, read on the
right axis, traces the information a single *new* comparison would carry
against an opponent at each location – anchored at the centre of the
scale, so it peaks at gap zero and falls away with the gap. The curve is
the visual explanation of why an adaptive design chases near-neighbour
contests: information is bought most cheaply where the two objects are
close, and a well-targeted design lifts the low dots by pairing their
objects against opponents near them.

## Usage

``` r
plot_btl_targeting(fit, grid = NULL)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- grid:

  Optional location grid for the reference curve.

## Value

Called for its plotting side effect; invisibly `NULL`.

## See also

[`btl_information`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md),
[`btl_next_pairs`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
plot_btl_targeting(btl(d, "a", "b", "win"))
```
