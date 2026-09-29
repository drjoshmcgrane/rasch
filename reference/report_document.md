# Write a Markdown, HTML, Word or PDF analysis report

Renders the active Rasch or paired-comparison fit as one Markdown file,
a self-contained HTML document, an editable Word document, or a PDF. The
report contains the principal estimates, model-specific tables,
diagnostic figures, and software provenance. Complete machine-readable
results remain available from
[`save_outputs`](https://drjoshmcgrane.github.io/rasch/reference/save_outputs.md).
Reports downloaded from the application retain compatible tailored item
shifts and externally weighted secondary person measures. For keyed fits
with repeated person IDs, the report explains why distractor analysis is
unavailable and retains the other model outputs.

## Usage

``` r
report_document(
  fit,
  file,
  format = c("auto", "md", "html", "docx", "pdf"),
  title = "Rasch measurement analysis",
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
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
  [`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md),
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md), or
  [`btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md).

- file:

  Output path ending in `.md`, `.html`, `.docx`, or `.pdf`.

- format:

  Output format. By default it is inferred from `file`. `"md"` writes
  GitHub-flavoured Markdown: the title is the first heading, the tables
  are pipe tables, and the diagnostic figures are left out so the report
  stays one text file that a person or a language model can read as it
  is.

- title:

  Report title.

- dif, bootstrap:

  Optional computed
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  and
  [`fit_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  results from this fit, rendered as run rather than recomputed at
  defaults.

- dif_bootstrap:

  Optional
  [`dif_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/dif_bootstrap.md)
  sensitivity analysis from this fit and DIF specification.

- dimensionality:

  Optional computed
  [`plot_scree`](https://drjoshmcgrane.github.io/rasch/reference/plot_scree.md)
  or
  [`btl_dimensionality`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md)
  result from this fit.

- invariance:

  Optional computed
  [`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  result from an EFRM fit.

- subtest:

  Optional computed
  [`dimensionality_test`](https://drjoshmcgrane.github.io/rasch/reference/dimensionality_test.md)
  result from this fit.

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

Invisibly, the output path.

## Details

Markdown, Word and HTML output require Pandoc, supplied with RStudio and
available through rmarkdown; Markdown needs Pandoc 2.8 or later. PDF
output also requires a LaTeX installation such as TinyTeX.
[`save_outputs`](https://drjoshmcgrane.github.io/rasch/reference/save_outputs.md)
writes the figures a Markdown report leaves out.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- rasch(matrix(rbinom(3000, 1, .5), 300, 10))
report_document(fit, file.path(tempdir(), "analysis.md"))
report_document(fit, file.path(tempdir(), "analysis.docx"))
} # }
```
