# Plot a Wright map

The conventional vertical person-item map (Wright and Stone 1979): the
person distribution to the left of a shared logit axis and the item
thresholds, labelled by item (and threshold number for polytomous
items), stacked to its right.

## Usage

``` r
plot_wright(fit, bins = 35, xlim = NULL, cex_labels = 0.8)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- bins:

  Number of bins for the person distribution and the threshold label
  rows.

- xlim:

  Optional logit range for the shared scale; persons and thresholds
  outside it are omitted.

- cex_labels:

  Character expansion for the threshold labels.

## Value

Called for its plotting side effect; invisibly `NULL`.

## References

Wright, B. D., & Stone, M. H. (1979). *Best Test Design*. Chicago: MESA
Press.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
plot_wright(rasch(X))
```
