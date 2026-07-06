# Compare fitted Rasch models

Builds a comparison table for two or more fits from
[`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md),
[`rasch_mfrm`](https://drjoshmcgrane.github.io/rmt/reference/rasch_mfrm.md),
or
[`rasch_efrm`](https://drjoshmcgrane.github.io/rmt/reference/rasch_efrm.md).
For fits of the same response data (identical item columns, maximum
scores, and number of persons) the pairwise conditional log-likelihoods
share their conditional information, and twice the difference from the
reference fit is reported with the difference in parameter counts; this
is descriptive (composite likelihood), and most meaningful for nested
structures such as RSM inside PCM. Across different data preparations
(subtests, splits, facet or frame structures) the likelihoods are not
comparable and the calibration-free columns carry the comparison: total
item-trait chi-square per degree of freedom, item and person fit
residual SDs (ideal 1), PSI, and alpha.

## Usage

``` r
compare_fits(..., reference = 1)
```

## Arguments

- ...:

  Two or more fitted objects, ideally named
  (`compare_fits(PCM = f1, RSM = f2)`).

- reference:

  Index or name of the reference fit for the log-likelihood difference;
  defaults to the first.

## Value

A data frame with one row per fit: label, model, persons, items,
parameters, log-likelihood, comparability with the reference,
`two_delta_ll` and `delta_parameters` (same-data fits only), chi-square
per df, fit residual SDs, PSI, and alpha.

## Examples

``` r
set.seed(1)
simP <- function(th, tau) { x <- 0:length(tau); p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
th <- rnorm(400)
X <- sapply(seq(-1, 1, length.out = 6), function(b)
  sapply(th, function(t) sample(0:3, 1, prob = simP(t, b + c(-0.8, 0, 0.8)))))
colnames(X) <- paste0("R", 1:6)
compare_fits(PCM = rasch(X, model = "PCM"), RSM = rasch(X, model = "RSM"))
#> Model comparison (reference: PCM)
#> 
#>  label model persons items two_delta_ll chisq_per_df item_fit_sd person_fit_sd
#>    PCM   PCM     400     6                     0.829       0.556         0.835
#>    RSM   RSM     400     6      -22.886        0.884       0.600         0.837
#>    PSI alpha
#>  0.745 0.784
#>  0.743 0.784
#> (further columns on the object: loglik, parameters, same_data)
#> 
#> two_delta_ll is a composite (pairwise) likelihood difference against the reference fit, reported only for fits of the same response data; it is descriptive, not chi-square calibrated. Across different data preparations compare chisq_per_df, the fit residual SDs (ideal 1), PSI, and alpha.
```
