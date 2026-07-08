# Person measures with extrapolated extreme scores

Extreme persons (zero or maximum raw score on their observed items) are
excluded from calibration, but they cannot be left out of group
comparisons; Andrich and Marais (2019, ch. 10) therefore describe an
extrapolated measure for them, continuing the growth of the
score-to-score differences so the last difference is the geometric mean
of its neighbours (see
[`score_table`](https://drjoshmcgrane.github.io/rasch/reference/score_table.md)).
This helper applies the same rule to the person table: for each
missing-data pattern with extreme persons, the score-to-measure
conversion over that pattern's items is extrapolated at its ends, and
the extreme persons receive the extrapolated location with the standard
error \\1/\sqrt{I(\theta)}\\ evaluated there. Non-extreme persons keep
their estimates unchanged. The extrapolation continues the Warm
(weighted likelihood) conversion, matching the package's person
estimates.

## Usage

``` r
person_extrapolated(fit)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  (equal discriminations; not EFRM).

## Value

The fit's person table with two added columns, `theta_extrapolated` and
`se_extrapolated`: equal to `theta` and `se` for non-extreme persons,
extrapolated for extreme persons. Patterns with fewer than three
interior scores cannot be extrapolated and keep their Warm values.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 8)
X <- matrix(rbinom(300 * 8, 1, plogis(outer(rnorm(300, 0, 2), d, "-"))), 300, 8)
colnames(X) <- paste0("I", 1:8)
fit <- rasch(X)
pe <- person_extrapolated(fit)
head(pe[pe$extreme, c("theta", "theta_extrapolated", "se", "se_extrapolated")])
#>        theta theta_extrapolated       se se_extrapolated
#> 1  -3.686068          -3.285814 1.677222        1.440928
#> 4   3.750309           3.362643 1.706623        1.478043
#> 11  3.750309           3.362643 1.706623        1.478043
#> 14 -3.686068          -3.285814 1.677222        1.440928
#> 15  3.750309           3.362643 1.706623        1.478043
#> 24 -3.686068          -3.285814 1.677222        1.440928
```
