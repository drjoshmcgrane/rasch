# Scree plot of the residual components with parallel analysis

Eigenvalues of the residual correlation matrix for the leading
components, with a model-simulated parallel-analysis reference:
responses are simulated from the calibrated model (observed missingness
kept), every person is re-estimated, and the residual eigenvalues
recomputed. Because estimating the person locations couples the
residuals within a person, this reference sits above the classical
random-normal one and is calibrated under the fitted model (Raiche 2005;
Chou & Wang 2010). Observed eigenvalues above the reference suggest
structure beyond what the model itself produces.

## Usage

``` r
plot_scree(fit, n_components = 10, parallel = TRUE, reps = 20)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- n_components:

  Number of leading components to display.

- parallel:

  Draw the parallel-analysis reference line.

- reps:

  Model-simulated replicates for the reference.

## Value

Called for its plotting side effect; invisibly the eigen table.

## References

Raiche, G. (2005). Critical eigenvalue sizes in standardized residual
principal components analysis. *Rasch Measurement Transactions*, 19(1),
1012.

Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in item
response models with principal component analysis on standardized
residuals. *Educational and Psychological Measurement*, 70(5), 717-731.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 8)
X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
colnames(X) <- paste0("I", 1:8)
plot_scree(rasch(X))
```
