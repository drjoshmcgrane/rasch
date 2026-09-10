# Save the outputs of a Rasch analysis

Writes the summary, estimates, diagnostic tables, person measures, and
model-specific results as CSV. Plots are written as PNG and, optionally,
PDF, together with a plain-text analysis summary. For MFRM and EFRM
fits, item estimates and response-cell diagnostics are saved separately.
Computed tailored item shifts and externally weighted secondary person
measures can be retained with the fitted-model output; the latter
include their resolved weights. For keyed fits with repeated person IDs,
the export records why distractor analysis is unavailable and retains
the other model outputs.

## Usage

``` r
save_outputs(
  fit,
  dir,
  formats = c("png", "pdf"),
  width = 9,
  height = 6,
  dpi = 300,
  item_plots = TRUE,
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

- dir:

  Output directory; created if absent. An existing directory must be
  empty; exports refuse to overwrite files from an earlier analysis.

- formats:

  Plot formats, any of `"png"` and `"pdf"`.

- width, height:

  Plot size in inches.

- dpi:

  PNG resolution.

- item_plots:

  Also write the per-item plot set (one ICC, category curve, threshold
  curve, and frequency chart per item).

- dif:

  Optional
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  result to export as computed — an application analysis carries the DIF
  model the analyst chose, which a default recomputation would silently
  replace. `NULL` computes the default when the fit carries person
  factors.

- bootstrap:

  Optional
  [`fit_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  result from this fit. Its item and person tables, or pair, object and
  judge tables for paired comparisons, join the export with the
  whole-test readings.

- dif_bootstrap:

  Optional
  [`dif_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/dif_bootstrap.md)
  sensitivity analysis from this fit and DIF specification.

- dimensionality:

  Optional computed
  [`plot_scree`](https://drjoshmcgrane.github.io/rasch/reference/plot_scree.md)
  or
  [`btl_dimensionality`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md)
  result from this fit. Supplying it keeps the exported table and figure
  identical to the analysis already run.

- invariance:

  Optional computed
  [`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  result from an EFRM fit.

- subtest:

  Optional computed
  [`dimensionality_test`](https://drjoshmcgrane.github.io/rasch/reference/dimensionality_test.md)
  result from this fit. Supplying it keeps the nominated item split and
  bootstrap calibration used in the analysis.

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

Invisibly, the vector of files written.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(150 * 6, 1, plogis(outer(rnorm(150), d, "-"))), 150, 6)
colnames(X) <- paste0("I", 1:6)
out <- tempfile("rasch-out-")
save_outputs(rasch(X), out, formats = "png", item_plots = FALSE, dpi = 96)
```
