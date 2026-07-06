# Principal components of the residual correlations

The first residual contrast (PC1) carries any second dimension; items
with opposing loadings define the split used by the unidimensionality
t-test. Loadings for the leading components and the eigenvalue table
support inspection beyond the first contrast.

## Usage

``` r
residual_pca(fit, n_components = 10)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- n_components:

  Number of leading components to return loadings and eigenvalue rows
  for (capped at the number of items).

## Value

A list with the residual `eigenvalues`, their `prop`ortions, the
first-contrast `loadings` (sorted), the `loadings_matrix` for the
leading components, the `eigen_table` (component, eigenvalue,
proportion, cumulative), and the `first_eigen`value.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 8)
X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
colnames(X) <- paste0("I", 1:8)
residual_pca(rasch(X))$first_eigen
#> [1] 1.295188
```
