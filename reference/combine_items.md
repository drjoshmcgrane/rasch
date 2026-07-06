# Combine items into subtests and re-analyse

Forms one polytomous super-item from each nominated group of items (its
score is the member sum; missing if any member is missing), keeps all
other items as they are, and refits the model with the same settings.
The usual treatment for item pairs flagged by
[`residual_correlations`](https://drjoshmcgrane.github.io/rmt/reference/residual_correlations.md).

## Usage

``` r
combine_items(fit, groups, model = "PCM")
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- groups:

  A list of character vectors, each naming two or more items to combine;
  a single vector is also accepted.

- model:

  Model for the re-analysis; defaults to `"PCM"`, which is almost always
  required because subtests change the maximum scores.

## Value

A new [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md)
fit on the combined structure, with the combinations recorded in its
notes.

## Examples

``` r
set.seed(1); Np <- 500; L <- 8
d <- seq(-2, 2, length.out = L)
X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
colnames(X) <- paste0("I", 1:L)
X[, 5] <- ifelse(runif(Np) < 0.9, X[, 4], X[, 5])   # dependent pair
fit <- rasch(X)
fit2 <- combine_items(fit, list(c("I4", "I5")))
fit2$items$item
#> [1] "I1"    "I2"    "I3"    "I6"    "I7"    "I8"    "I4+I5"
```
