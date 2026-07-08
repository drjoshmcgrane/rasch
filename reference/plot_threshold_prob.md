# Plot threshold probability curves

Conditional probability of success at each threshold,
`P(X = k | X = k - 1 or k)`, a logistic ogive crossing 0.5 at the
threshold location. Disordered thresholds are immediately visible as
out-of-sequence ogives. With `observed = TRUE` the observed conditional
proportions per class interval are overlaid, the direct check on whether
each threshold discriminates (and hence whether collapsing categories
could ever be justified; Andrich and Marais 2019, ch. 22).

## Usage

``` r
plot_threshold_prob(
  fit,
  item,
  grid = seq(-6, 6, 0.05),
  observed = FALSE,
  n_groups = fit$n_groups
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- item:

  Item name or column index.

- grid:

  Logit grid over which to draw the curves.

- observed:

  Overlay the observed conditional threshold proportions per class
  interval.

- n_groups:

  Class intervals for the observed points.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
simP <- function(th, t) { x <- 0:length(t); p <- exp(x * th - c(0, cumsum(t))); p / sum(p) }
th <- rnorm(400)
X <- sapply(1:4, function(i) sapply(th, function(t) sample(0:3, 1, prob = simP(t, c(-1, 0, 1)))))
colnames(X) <- sprintf("P%02d", 1:4)
plot_threshold_prob(rasch(X), "P01")
```
