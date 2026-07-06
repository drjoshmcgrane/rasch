# Write a self-contained HTML report of a Rasch analysis

Builds a single portable HTML file containing the complete analysis: the
summary statistics, every diagnostic table, and every test-level plot
embedded as an image, styled for reading and sharing. The file has no
external dependencies, so it can be e-mailed or archived as the record
of an analysis.

## Usage

``` r
report_html(fit, file, title = "Rasch measurement analysis", dpi = 150)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- file:

  Path of the HTML file to write.

- title:

  Report title.

- dpi:

  Resolution of the embedded plots.

## Value

Invisibly, `file`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
out <- file.path(tempdir(), "report.html")
report_html(rasch(X), out)
```
