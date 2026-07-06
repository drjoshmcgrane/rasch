# Plot a person characteristic curve

The person characteristic curve: the probability of success as a
function of item location at the person's estimated measure, with the
person's observed responses overlaid, grouped into item-difficulty
intervals (proportion of maximum score per interval). Erratic responding
(for example lucky guessing on hard items by a low-proficiency person)
shows as observed points far from the curve, complementing the person
fit residual.

## Usage

``` r
plot_pcc(fit, person, n_groups = 5, grid = seq(-5, 5, 0.05))
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- person:

  Row number of the person, or an ID matching `fit$person$id`.

- n_groups:

  Number of item-difficulty intervals for the observed points (capped by
  the number of observed items).

- grid:

  Item-location grid over which to draw the curve.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 12)
X <- matrix(rbinom(300 * 12, 1, plogis(outer(rnorm(300), d, "-"))), 300, 12)
colnames(X) <- paste0("I", 1:12)
plot_pcc(rasch(X), person = 1)
```
