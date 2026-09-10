# Launch the rasch point-and-click graphical interface

Opens the Shiny application for fitting models and examining their
tables, plots and diagnostics. The R code for each result is available
in the app. Analyses can be saved and reopened, or exported as HTML,
Word or PDF reports.

## Usage

``` r
run_app(...)
```

## Arguments

- ...:

  Passed to
  [`shiny::runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Value

Called for its side effect of launching the app.

## Details

The app's interface packages ('shiny', 'bslib', 'DT', 'bsicons', and
'callr' for cancellable EFRM and fit-bootstrap estimation) are suggested
rather than required by the package. If any are missing, `run_app` lists
them all and, in an interactive session, offers to install them before
launching.

Older saved analyses are checked against the current person-scoring
algorithm. If their scores differ, refit the analysis before reopening
it; the original file is left unchanged. Its source data can be
recovered with `readRDS(file)$data`. Saved EFRM and CJ frame fits
without a current likelihood-check record also require refitting. Their
settings remain in `readRDS(file)$settings`. Superseded DIF results, or
CJ DIF without verified judge-role alignment, are omitted with a
warning. Rerun those analyses before reporting them.

## Examples

``` r
if (interactive()) run_app()
```
