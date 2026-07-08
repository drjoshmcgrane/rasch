# Plot category probability curves

Plot category probability curves

## Usage

``` r
plot_ccc(
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

  Overlay the observed category proportions per class interval (Andrich
  and Marais 2019, ch. 20).

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
plot_ccc(rasch(X), "P01", observed = TRUE)
```
