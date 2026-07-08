# Save a kidmap for every person

Writes one kidmap
([`plot_kidmap`](https://drjoshmcgrane.github.io/rasch/reference/plot_kidmap.md))
per person to a single multi-page PDF or a ZIP archive of PNGs, chosen
by the extension of `file`. Persons without a location estimate are
skipped.

## Usage

``` r
save_person_plots(
  fit,
  file,
  persons = NULL,
  level = 0.95,
  width = 8,
  height = 6,
  dpi = 300
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- file:

  Output path ending in `.pdf` (one page per person) or `.zip` (one PNG
  per person).

- persons:

  Row numbers or IDs; all estimated persons by default.

- level:

  Confidence level of the band marking unexpected responses.

- width, height, dpi:

  Device size in inches and PNG resolution.

## Value

Invisibly, the output path.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(60 * 6, 1, plogis(outer(rnorm(60), d, "-"))), 60, 6)
colnames(X) <- paste0("I", 1:6)
f <- rasch(X)
save_person_plots(f, file.path(tempdir(), "kidmaps.pdf"), persons = 1:5)
```
