# Estimate the magnitude of response dependence between two items

Quantifies how strongly a dependent item's response follows an
independent item's response, in logits, by the resolution method of
Andrich and Kreiner (2010; polytomous generalisation Andrich, Humphry
and Marais 2012), by resolution of the dependent item. The dependent
item is resolved into one item per category of the independent item
(each carrying the responses of the persons who gave that category),
both original items are removed, and the model refitted. Under
dependence of magnitude \\d\\, threshold \\k\\ of the resolved item for
category \\x_i\\ is shifted by \\-d\\ when \\k \le x_i\\ and \\+d\\
otherwise, so each threshold yields \\\hat d_k =
(\hat\delta\_{ji(k)}(x_i = k-1) - \hat\delta\_{ji(k)}(x_i = k))/2\\ and
\\\hat d\\ is their mean (eq. 24.7 of Andrich and Marais 2019). Because
the resolved items are answered by disjoint persons, the estimates are
independent, and the standard error pools the threshold variances:
\\\hat\sigma^2_k = (\hat\sigma^2\_{(k)(k-1)} +
\hat\sigma^2\_{(k)(k)})/4\\, with \\V\[\hat d\] =
\bar{\hat\sigma^2_k}/m\\ (eqs. 24.9-24.11).

## Usage

``` r
dependence_magnitude(fit, dependent, independent)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- dependent, independent:

  Item names or indices: the item hypothesised to depend, and the item
  it depends on. Both must share the same maximum score (the
  formalisation requires it).

## Value

A list of class `"rasch_dependence"`: the estimate `d`, its `se`, `z`
and `p` for the hypothesis \\d = 0\\, the per-threshold table
`thresholds` (columns `k`, `delta_lo`, `delta_hi`, `d_k`, `se_k`), and
the resolved `refit`.

## References

Andrich, D. and Kreiner, S. (2010). Quantifying response dependence
between two dichotomous items using the Rasch model. Applied
Psychological Measurement, 34, 181-192. Andrich, D., Humphry, S. M. and
Marais, I. (2012). Quantifying local, response dependence between two
polytomous items using the Rasch model. Applied Psychological
Measurement, 36, 309-324.

## Examples

``` r
set.seed(1); N <- 700
d0 <- seq(-1.5, 1.5, length.out = 8)
X <- matrix(rbinom(N * 8, 1, plogis(outer(rnorm(N), d0, "-"))), N, 8)
X[, 5] <- ifelse(runif(N) < 0.75, X[, 4], X[, 5])   # I5 follows I4
colnames(X) <- paste0("I", 1:8)
dependence_magnitude(rasch(X), dependent = "I5", independent = "I4")
#> Response dependence of I5 on I4 (Andrich & Kreiner resolution)
#>   d = 1.992 logits (se 0.126), z = 15.81, p = < 0.001
```
