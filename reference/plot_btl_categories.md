# Plot graded-comparison category curves

For a graded paired-comparison fit, the probability of each response
category as a function of the location difference `beta_a - beta_b`,
with the symmetric threshold structure marked. The display is the
paired-comparison counterpart of the category probability curves of a
polytomous item.

## Usage

``` r
plot_btl_categories(fit, grid = seq(-4, 4, 0.05))
```

## Arguments

- fit:

  A graded fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) (with
  `response`).

- grid:

  Difference grid, in logits.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
P <- vapply(seq_len(nrow(d)), function(r)
  item_moments(beta[d$a[r]] - beta[d$b[r]], c(-1, 0, 1))$P, numeric(4))
d$grade <- apply(P, 2, function(p) sample(0:3, 1, prob = p))
plot_btl_categories(btl(d, "a", "b", response = "grade"))
```
