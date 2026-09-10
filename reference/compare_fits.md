# Compare fitted Rasch models

Builds a comparison table for two or more fits from
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
[`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md),
or (all together)
[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md). For
fits of the same response data using the same likelihood contributions,
twice the difference from the reference fit is reported with the
difference in parameter counts; this is descriptive (composite
likelihood), and most meaningful for nested structures such as RSM
inside PCM. Paired-comparison data identity is evaluated separately for
each fit against the reference, including sequence fields shared by that
pair. Order values are compared through their within-judge ranks, so
relabelling that preserves order and ties does not change identity.
Adding another model does not change an existing pair's compatibility.

## Usage

``` r
compare_fits(..., reference = 1)
```

## Arguments

- ...:

  Two or more fitted objects, preferably given unique names. Supply
  either all Rasch-family fits or all `btl` fits. For `btl`, fits of the
  same comparison data (same objects, response scale, comparisons, and
  judges) support the likelihood columns – e.g. free versus
  principal-component thresholds, with and without a position effect or
  within-judge dependence. Row order, arbitrary person or judge labels,
  and expansion versus count compression do not change data identity;
  allocation to independent persons or judges does.

- reference:

  Index or name of the reference fit for the log-likelihood difference;
  defaults to the first.

## Value

A data frame with one row per fit: label, model, persons, items (judges,
objects, comparisons for `btl`), parameters, log-likelihood,
`eff_params`, `cl_aic`, `cl_bic`, response-data identity with the
reference (`same_data`), `two_delta_ll` and `delta_parameters` (eligible
same-data comparisons only), chi-square per df, fit residual SDs, PSI,
and alpha (OSI for `btl`).

## Details

The calibrated comparison is carried by the composite-likelihood
information criteria `cl_aic` (Varin and Vidoni 2005) and `cl_bic` (Gao
and Song 2010): \\-2\\cl + c \cdot tr(H^{-1}J)\\, with \\c = 2\\ or
\\\log n\\. Because every response enters every pair its item forms, the
pairwise log-likelihood over-counts the data; the effective parameter
count \\tr(H^{-1}J)\\ from the Godambe matrices – the same quantity
whose eigenvalues calibrate
[`lr_test`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)
– accounts for composite score variability and curvature, which the
nominal parameter count would not. \\n\\ counts independent units:
persons contributing at least one informative pair, or judges for
paired-comparison fits (count-weighted comparisons when unclustered).
Smaller is better; the criteria are valid across models of the same data
whether or not they nest, and are `NA` (with the reason in the printed
note) for EFRM fits, which do not carry Godambe matrices for the
complete linked model. Rasch EFRM log-likelihoods cover only within-set
calibration pairs, so generic likelihood differences involving these
fits are also withheld. Use `fit$efrm_vs_rasch` for the descriptive
group-unit comparison on matched pairs; it does not assess set units.
MFRM fits carry the sensitivity and sandwich covariance for their
structural pairwise calibration and therefore receive the same criteria
as ordinary and explanatory Rasch fits. BTL–EFRM fits also use a
two-stage estimator rather than maximising the combined objective
jointly, so their generic information criteria and likelihood
differences are withheld; use the fit's `equal_unit` component for its
labelled descriptive comparison.

Across different data preparations (subtests, splits, facet or frame
structures), or different allocations of response rows to persons, the
likelihood-based criteria are not comparable and are withheld. The table
retains descriptive context: total trait chi-square per degree of
freedom, calibration and person fit-residual SDs (ideal 1), PSI, and
alpha where applicable (OSI for paired comparisons). Alpha is `NA` when
an MFRM or EFRM item is represented by several response cells. These
columns do not provide a formal selection test across different response
data.

## Examples

``` r
set.seed(1)
simP <- function(th, tau) {
  x <- 0:length(tau)
  p <- exp(x * th - c(0, cumsum(tau)))
  p / sum(p)
}
th <- rnorm(400)
X <- sapply(seq(-1, 1, length.out = 6), function(b)
  sapply(th, function(t)
    sample(0:3, 1, prob = simP(t, b + c(-0.8, 0, 0.8)))))
colnames(X) <- paste0("R", 1:6)
compare_fits(PCM = rasch(X, model = "PCM"),
             RSM = rasch(X, model = "RSM"))
#> Model comparison (reference: PCM)
#> 
#>  label model persons items eff_params   cl_aic   cl_bic two_delta_ll
#>    PCM   PCM     400     6     61.697 8108.103 8353.432             
#>    RSM   RSM     400     6     21.817 8051.228 8137.980      -22.886
#>  chisq_per_df item_fit_sd person_fit_sd   PSI alpha
#>         0.829       0.556         0.835 0.745 0.784
#>         0.884       0.600         0.837 0.743 0.784
#> (further columns on the object: loglik, parameters, same_data)
#> 
#> cl_aic and cl_bic are composite-likelihood information criteria (Varin & Vidoni 2005; Gao & Song 2010): -2 cl penalised by the effective parameter count tr(H^-1 J), which accounts for composite score variability and curvature; smaller is better, valid across models of the same data. Information criteria and two_delta_ll are withheld when the response data, response scale or independent-unit allocation differ from the reference. two_delta_ll is the raw composite difference against the reference, descriptive only. Across different data preparations, chisq_per_df, the fit residual SDs and separation/reliability columns provide descriptive context rather than a formal selection criterion.
```
