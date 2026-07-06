# Launch the rmt graphical interface

Opens the Shiny application: data upload with ID, person-factor, and
item column nomination; the full analysis; interactive tables and plots;
and one-click export of every table and plot.

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

## Examples

``` r
if (interactive()) run_app()
```
