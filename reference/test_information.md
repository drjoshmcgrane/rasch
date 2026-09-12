# Test information function

Fisher information over a grid of person locations, with the
corresponding standard error of measurement. Ordinary Rasch fits return
one whole-test curve. EFRM fits return one curve per person group and
per item administration pattern actually observed within that group (in
a linking design, persons who took only the core set get a core-only
curve, and the linking subsample gets the pooled one). MFRM fits return
one curve per observed item-by-facet pattern for a person, so ratings
that jointly inform the same person measure are added and mutually
exclusive designs remain separate. Partly answered sets or facet
conditions contribute only their observed items; a missing response is
not treated as an administered item when defining these patterns. Where
item nonresponse leaves nearly every person a pattern of their own, that
is what these fits return: an unanswered item carries no information
about the person who left it, so no pattern is merged into a fuller one
and no curve of theirs is drawn over a design nobody was administered.

## Usage

``` r
test_information(fit, grid = NULL, items = NULL)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- grid:

  Logit grid over which to evaluate the information.

- items:

  Optional item selection: item names or indices. Every design block is
  restricted to the named items, so a restricted person-item map can
  carry the information of its own selection. The `design` labels are
  those of the restricted blocks, and blocks that differ only outside
  the selection are returned once.

## Value

A data frame with `theta`, `info`, and `sem`. For EFRM and MFRM fits it
also contains a `design` column identifying the administrable frame or
facet design.

## Details

For an administrable block \\\mathcal A\\, the information and standard
error of measurement are \$\$I(\theta)=\sum\_{i\in\mathcal A}d_i^2
\operatorname{Var}(X_i\mid\theta),\qquad
\operatorname{SEM}(\theta)=I(\theta)^{-1/2},\$\$ where \\d_i\\ is the
frame unit or discrimination multiplier. For an ordinary Rasch fit,
\\d_i=1\\. Information is returned only for a converged calibration.
Comparative Judgement designs use
[`btl_information`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md)
instead.

## References

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

## See also

[`targeting_table`](https://drjoshmcgrane.github.io/rasch/reference/targeting_table.md)
and
[`plot_tif`](https://drjoshmcgrane.github.io/rasch/reference/plot_tif.md).

## Examples

``` r
set.seed(1)
d <- seq(-1.5, 1.5, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
head(test_information(rasch(X)))
#>   theta       info      sem
#> 1  -6.0 0.02364631 6.503068
#> 2  -5.9 0.02609577 6.190346
#> 3  -5.8 0.02879468 5.893101
#> 4  -5.7 0.03176749 5.610590
#> 5  -5.6 0.03504089 5.342105
#> 6  -5.5 0.03864388 5.086976
```
