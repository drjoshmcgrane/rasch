# Plot a within-judge dependence effect

The graphical display of a paired-comparison dependence effect, the
counterpart of the DIF characteristic curve. For every comparison the
departure of the observed response from what the object locations alone
predict is taken, and the contribution of the *other* dependence effect
is removed (a partial-residual display); these departures are then
averaged in bins of the effect's own history covariate and plotted
against it, with the model's fitted contribution overlaid. Observed
points that rise with the covariate along the fitted line are the effect
the coefficient summarises; a flat, scattered cloud means the estimate
rests on little. Only the informative comparisons (a non-zero covariate:
the two objects' histories differ) carry the effect, and the count in
each bin is printed so a thin exposure tail is visible.

## Usage

``` r
plot_btl_dependence(fit, effect = c("exposure", "carry_over"), bins = 6)
```

## Arguments

- fit:

  An object from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) fitted
  with an `order` column, so `fit$dependence_data` is present.

- effect:

  Which effect to display: `"exposure"` (the seen-before advantage, the
  default) or `"carry_over"` (response dependence).

- bins:

  Number of covariate bins for the continuous carry-over display;
  exposure takes its three natural levels (-1, 0, +1).

## Value

Called for its plotting side effect; invisibly a data frame of the
binned covariate value, observed and fitted departure, and bin count.

## References

Davidson, R. R., & Beaver, R. J. (1977). On extending the Bradley-Terry
model to incorporate within-pair order effects. *Biometrics*, 33(4),
693-702.

## Examples

``` r
set.seed(1)
beta <- c(A = -0.8, B = -0.2, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
d$judge <- sample(sprintf("J%02d", 1:8), nrow(d), TRUE)
d <- d[order(d$judge), ]; d$t <- ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
f <- btl(d, "a", "b", winner = "win", judge = "judge", order = "t")
plot_btl_dependence(f, "carry_over")
```
