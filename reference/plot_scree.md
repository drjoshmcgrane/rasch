# Scree plot of the residual components with parallel analysis

Eigenvalues of the residual correlation matrix for the leading
components, with a parallel-analysis reference: the mean eigenvalues of
residual-sized random normal matrices sharing the data's missingness
pattern. Observed eigenvalues above the reference suggest structure
beyond chance.

## Usage

``` r
plot_scree(fit, n_components = 10, parallel = TRUE, reps = 20)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- n_components:

  Number of leading components to display.

- parallel:

  Draw the parallel-analysis reference line.

- reps:

  Random replicates for the reference.

## Value

Called for its plotting side effect; invisibly the eigen table.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 8)
X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
colnames(X) <- paste0("I", 1:8)
plot_scree(rasch(X))
```
