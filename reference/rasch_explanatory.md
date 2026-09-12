# Fit an explanatory Rasch model

Fits the linear logistic test model (LLTM) for dichotomous responses or
the linear partial credit model (LPCM) for polytomous responses. Item or
threshold locations are linear functions of observed predictors. The
response model remains Rasch and is estimated by pairwise conditional
maximum likelihood.

## Usage

``` r
rasch_explanatory(
  data,
  predictors,
  formula,
  items = NULL,
  level = c("item", "threshold"),
  id = NULL,
  factors = NULL,
  n_groups = NULL,
  na_codes = -1,
  key = NULL,
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- data, items, id, factors, n_groups, na_codes, key, maxit, tol:

  As in
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- predictors:

  Data frame containing an `item` column and the predictors named in
  `formula`. With `level = "threshold"`, it must also contain
  `threshold`, with one row for every fitted item threshold. Column
  names must be unique.

- formula:

  One-sided explanatory formula. For example,
  `~ format + operation + format:operation`. The reserved `threshold`
  factor permits threshold-specific effects. Formula offsets
  ([`offset()`](https://rdrr.io/r/stats/offset.html)) are not supported.

- level:

  Whether `predictors` contains one row per `"item"` or per
  `"threshold"`. Item rows are expanded over their thresholds.

## Value

An object of class `"rasch_explanatory"` inheriting from `"rasch"`.
Standard item, person, fit and diagnostic components use the explanatory
thresholds. The `explanatory` component contains the formula, metadata
and design matrices; `reference_fit` is the free PCM calibration.
`est$coefficients` reports the estimates, standard errors, \\t\\
statistics, reference degrees of freedom, raw probabilities and
Holm-adjusted probabilities.

## Details

For threshold \\k\\ of item \\i\\,
\$\$\delta\_{ik}=z\_{ik}^{T}\gamma.\$\$ The adjacent-category log odds
are \$\$\log\\P(X\_{ni}=k)/P(X\_{ni}=k-1)\\=\theta_n-\delta\_{ik}.\$\$
The threshold origin is fixed to the same mean-item-location zero used
by [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).
An intercept therefore sets the arbitrary origin and is not separately
estimated. Numeric predictors are continuous, unordered factors are
categorical, and ordered factors use successive contrasts between
adjacent levels. Character predictors are converted to unordered
factors. The reserved factor `threshold` identifies the within-item
threshold number; `threshold_number` supplies its integer value. Design
columns are centred and rescaled internally for numerical stability;
reported coefficients and standard errors use the supplied predictor
units. Coincident coefficient labels receive numeric suffixes; this does
not change the predictor design.

A free PCM reference is fitted to the same prepared responses and
retained on the object.
[`explanatory_test`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_test.md)
applies the first-order Kent calibration required for the pairwise
composite likelihood. When an identifier occurs on more than one
response row, coefficient covariance is clustered by person. A
linearised delete-one-person correction accounts for finite-cluster
leverage without refitting the model once per person. Supported fits use
a \\t\\ reference with degrees of freedom equal to the number of
independent person units contributing conditional information minus one,
whether or not identifiers repeat; inference is withheld when the
calibration lacks enough independent information. Holm adjustment covers
the coefficient family. With few persons and unequal numbers of response
rows, these approximate tests can still be mildly liberal; the
correction does not guarantee nominal coverage in small samples.

## References

Fischer, G. H. (1973). The linear logistic test model as an instrument
in educational research. Acta Psychologica, 37, 359–374.

Fischer, G. H. and Ponocny, I. (1994). An extension of the partial
credit model with an application to the measurement of change.
Psychometrika, 59, 177–192.

## See also

[`explanatory_test`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_test.md),
[`explanatory_diagnostics`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_diagnostics.md),
and
[`relax_explanatory`](https://drjoshmcgrane.github.io/rasch/reference/relax_explanatory.md).

## Examples

``` r
set.seed(1)
q <- data.frame(item = paste0("I", 1:8),
                operation = rep(0:1, each = 4),
                format = rep(c("A", "B"), 4))
difficulty <- -1 + 0.7 * q$operation + 0.4 * (q$format == "B")
X <- matrix(rbinom(500 * 8, 1,
  plogis(outer(rnorm(500), difficulty, "-"))), 500, 8)
colnames(X) <- q$item
fit <- rasch_explanatory(X, predictors = q,
                         formula = ~ operation + format)
fit$est$coefficients
#>       term estimate    se     t  df       p   p_adj
#>  operation    0.697 0.070 9.996 449 < 0.001 < 0.001
#>    formatB    0.423 0.071 5.925 449 < 0.001 < 0.001
explanatory_test(fit)
#>  model parameters free_parameters r_squared r_squared_adj              r2_basis
#>   LLTM          2               7     0.945         0.923 threshold calibration
#>   chisq df p_naive chisq_kent     p p_kent
#>  29.999  5 < 0.001      7.585 0.181  0.181
```
