# Warm's weighted likelihood estimates by raw score

Computes the weighted likelihood estimate (WLE) of person location for
every possible raw score on a set of items, with standard errors. WLE
estimates are finite at the extreme (zero and maximum) scores, unlike
the maximum likelihood estimate.

## Usage

``` r
person_wle(tau_list, disc = 1)
```

## Arguments

- tau_list:

  List of per-item threshold vectors.

- disc:

  Common discrimination (frame unit) of the items; with a constant
  discrimination the raw score remains sufficient.

## Value

A list with `theta` and `se`, each named by raw score.

## Details

For raw score \\R\\, let \\E(\theta)\\, \\V(\theta)\\, and
\\\mu_3(\theta)\\ be the sums of the item expected scores, variances,
and third central moments. The estimate solves Warm's weighted score
equation \$\$R-E(\theta)+\frac{\mu_3(\theta)}{2V(\theta)}=0.\$\$ When
the equation has several solutions, competing maxima are compared using
log likelihood plus one half log information. Equal maxima use the lower
location; their average need not maximize the weighted likelihood. With
common discrimination \\d\\, its explicit multiplier cancels from this
equation, although the moments are evaluated under \\d\\. The reported
standard error is \$\$\operatorname{SE}(\hat{\theta})=
\\d^2V(\hat{\theta})\\^{-1/2}.\$\$

## References

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. Psychometrika, 54(3), 427–450.

Zhang, J. (2005). Bias correction for the maximum likelihood estimate of
ability. ETS Research Report RR-05-15.

## See also

[`score_table`](https://drjoshmcgrane.github.io/rasch/reference/score_table.md)
and
[`person_extrapolated`](https://drjoshmcgrane.github.io/rasch/reference/person_extrapolated.md).

## Examples

``` r
person_wle(list(c(-1, 0), c(-0.5, 0.5), c(0, 1)))
#> $theta
#>          0          1          2          3          4          5          6 
#> -2.4444062 -1.2371268 -0.5611339  0.0000000  0.5611339  1.2371268  2.4444062 
#> 
#> $se
#>         0         1         2         3         4         5         6 
#> 1.5455221 0.9731988 0.8347598 0.8003924 0.8347598 0.9731988 1.5455221 
#> 
```
