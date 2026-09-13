# Write a self-contained HTML report of a Rasch analysis

Writes one HTML file containing the summary statistics, diagnostic
tables, and test-level plots. Images and styles are embedded in the
file. Computed tailored item shifts and externally weighted secondary
person measures can be included with the fitted-model results. For keyed
fits with repeated person IDs, the report explains why distractor
analysis is unavailable; the other model outputs remain available.

## Usage

``` r
report_html(
  fit,
  file,
  title = "Rasch measurement analysis",
  dpi = 150,
  dif = NULL,
  bootstrap = NULL,
  dif_bootstrap = NULL,
  dimensionality = NULL,
  invariance = NULL,
  subtest = NULL,
  tailored = NULL,
  person_weights = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- file:

  Path of the HTML file to write.

- title:

  Report title.

- dpi:

  Resolution of the embedded plots.

- dif, bootstrap:

  Optional computed
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  and
  [`fit_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  results from this fit, exported as run; the DIF table is otherwise
  recomputed at defaults when the fit carries person factors.

- dif_bootstrap:

  Optional
  [`dif_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/dif_bootstrap.md)
  sensitivity analysis from this fit and DIF specification.

- dimensionality:

  Optional computed
  [`plot_scree`](https://drjoshmcgrane.github.io/rasch/reference/plot_scree.md)
  result from this fit. Supplying it keeps the table and figure
  identical to the analysis already run.

- invariance:

  Optional computed
  [`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  result from an EFRM fit.

- subtest:

  Optional computed
  [`dimensionality_test`](https://drjoshmcgrane.github.io/rasch/reference/dimensionality_test.md)
  result from this fit. Supplying it keeps the item split and bootstrap
  calibration used in the analysis.

- tailored:

  Optional computed
  [`tailored_analysis`](https://drjoshmcgrane.github.io/rasch/reference/tailored_analysis.md)
  result from this ordinary dichotomous fit.

- person_weights:

  Optional table returned by
  [`weighted_person_estimates`](https://drjoshmcgrane.github.io/rasch/reference/weighted_person_estimates.md).
  An explicit table takes precedence over a compatible result retained
  by the application.

## Value

Invisibly, `file`.

## Examples

``` r
set.seed(1)
d <- seq(-1.5, 1.5, length.out = 4)
X <- matrix(rbinom(80 * 4, 1, plogis(outer(rnorm(80), d, "-"))), 80, 4)
colnames(X) <- paste0("I", 1:4)
out <- file.path(tempdir(), "report.html")
report_html(rasch(X), out, dpi = 96)
```
