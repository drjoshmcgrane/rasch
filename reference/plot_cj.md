# Plot a joint calibration of responses and judgements

Caterpillar plot of the item locations from
[`rasch_cj`](https://drjoshmcgrane.github.io/rasch/reference/rasch_cj.md):
the combined location of each item with its 95 per cent interval, and
beside it the location each frame (the responses, the comparisons, the
rankings, or each test when several are linked) gives the item on its
own, expressed on the reference scale. The spread of the frame markers
around the combined location is the evidence the invariance test
summarises; an item whose location differs between frames at a
Holm-adjusted p below .05 in the per-object invariance table is drawn in
red. A frame that does not reach an item leaves no marker for it.

## Usage

``` r
plot_cj(fit, frames = TRUE)
```

## Arguments

- fit:

  An item-mode object from
  [`rasch_cj`](https://drjoshmcgrane.github.io/rasch/reference/rasch_cj.md).

- frames:

  Logical; draw the separate calibration of each frame beside the
  combined location (the default), or the combined locations alone.

## Value

Called for its plotting side effect; invisibly `NULL`.

## See also

[`rasch_cj`](https://drjoshmcgrane.github.io/rasch/reference/rasch_cj.md),
[`plot_pl`](https://drjoshmcgrane.github.io/rasch/reference/plot_pl.md).

## Examples

``` r
set.seed(1)
delta <- seq(-1.5, 1.5, length.out = 8)
names(delta) <- sprintf("I%02d", 1:8)
theta <- rnorm(300)
X <- sapply(delta, function(d) as.integer(runif(300) < plogis(theta - d)))
pairs <- t(combn(names(delta), 2))[sample(28, 200, replace = TRUE), ]
p_a <- plogis(0.6 * (delta[pairs[, 1]] - delta[pairs[, 2]]))
cj <- data.frame(a = pairs[, 1], b = pairs[, 2],
                 winner = ifelse(runif(200) < p_a, pairs[, 1], pairs[, 2]))
fit <- rasch_cj(X, comparisons = cj, object_a = "a", object_b = "b",
                winner = "winner")
plot_cj(fit)
```
