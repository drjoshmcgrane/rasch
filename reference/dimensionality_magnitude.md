# Magnitude of multidimensionality from a subtest analysis

Estimates how strongly two or more hypothesised subscales measure
distinct traits, by Andrich's (2016) comparison of two reliability
calculations: one treating all items as independent (which inflates
reliability under multidimensionality) and one on the subtest analysis
in which each subscale is combined into a single polytomous super-item
(which absorbs the unique subscale variance). Under the bifactor
formalisation \\\beta\_{ns} = \beta_n + c\\\beta'\_{ns}\\ (Marais and
Andrich 2008), with \\S\\ subscales of \\K\\ items, \$\$c^2 =
S\\(r_1/r_2 - 1) \frac{SK - 1}{S(K - 1)},\$\$ the latent correlation
between subscales is \\\rho = 1/(1 + c^2)\\, and \\A = S/(S + c^2)\\ is
the proportion of common (non-unique, non-error) variance. Both the
person separation index and coefficient alpha versions are reported
(Andrich and Marais 2019, ch. 24).

## Usage

``` r
dimensionality_magnitude(fit, subtests)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- subtests:

  A list of character vectors assigning *every* item of the fit to one
  subscale (at least two subscales of two or more items). Unequal
  subscale sizes use their mean as \\K\\.

## Value

A list of class `"rmt_dim_magnitude"`: the comparison `table` (rows PSI
and alpha; columns `run1`, `subtest`, `c2`, `c`, `rho`, `A`), the
subtest `refit`, and the design constants `S` and `K`.

## References

Andrich, D. (2016). Components of variance of scales with a bifactor
structure from two calculations of coefficient alpha. Educational
Measurement: Issues and Practice, 35(4), 25-30.

## Examples

``` r
set.seed(1); N <- 500
common <- rnorm(N); u1 <- rnorm(N); u2 <- rnorm(N)
d <- rep(seq(-1, 1, length.out = 5), 2)
X <- sapply(1:10, function(i) rbinom(N, 1,
  plogis(common + 0.8 * (if (i <= 5) u1 else u2) - d[i])))
colnames(X) <- paste0("I", 1:10)
fit <- rasch(X)
dimensionality_magnitude(fit,
  list(paste0("I", 1:5), paste0("I", 6:10)))$table
#>   index      run1   subtest        c2         c       rho         A
#> 1   PSI 0.5770681 0.3690656 1.2680828 1.1260918 0.4409010 0.6119796
#> 2 alpha 0.6618793 0.4969475 0.7467516 0.8641479 0.5724912 0.7281328
```
