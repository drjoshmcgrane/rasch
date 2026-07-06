# Factorial DIF analysis with Tukey comparisons

Models all nominated person factors jointly: for each item the
standardised residuals are analysed by the full factorial of the person
factors crossed with the trait class interval,
`z ~ (f1 * f2 * ...) * ci`. Terms not involving the class interval are
uniform DIF effects (main effects and factor-by-factor interactions);
terms involving it are non-uniform. Probabilities are adjusted across
items within each term (Benjamini-Hochberg by default). A significant
interaction supersedes the main effects (and lower-order interactions)
of the factors it involves, which is recorded in the `superseded`
column; interpret the highest-order significant terms. Tukey HSD
comparisons are returned for every significant, non-superseded group
term (the cell-mean contrasts for interactions), with Tukey's own
familywise adjustment within each term.

## Usage

``` r
dif_anova_factorial(
  fit,
  factors = NULL,
  n_groups = NULL,
  p_adjust = "BH",
  alpha = 0.05,
  effects = c("factorial", "main"),
  sizes = FALSE,
  id = NULL,
  within = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- factors:

  As in
  [`dif_anova`](https://drjoshmcgrane.github.io/rmt/reference/dif_anova.md);
  at least one factor, usually two or more.

- n_groups:

  Number of trait class intervals. By default set from the smallest
  factor-combination cell so every interval-by-cell count keeps about 30
  expected responses (between 2 and 10 intervals); the value used is
  returned as `n_groups`.

- p_adjust:

  Multiplicity adjustment across items within each term; default `"BH"`.

- alpha:

  Significance level applied to the adjusted probabilities.

- effects:

  `"factorial"` (default) crosses every person factor with every other
  and with the class interval; `"main"` fits the factors additively
  (each factor's main effect and its interaction with the class
  interval, but no factor-by-factor terms).

- sizes:

  Also compute DIF magnitudes in logits
  ([`dif_size`](https://drjoshmcgrane.github.io/rmt/reference/dif_size.md))
  for every significant, non-superseded group term: the item is resolved
  by the term's levels (interaction terms by their cells) and all
  pairwise location differences are returned with Holm familywise
  adjustment and the practical-significance flag. Each size involves a
  re-analysis, so this costs one refit per flagged item-term.

- id, within:

  Person identifier and within-subject factor names for stacked
  repeated-measures designs, as in
  [`dif_anova`](https://drjoshmcgrane.github.io/rmt/reference/dif_anova.md).
  When any factor is within-subject the model becomes a mixed
  (split-plot) analysis of variance – the class interval taken at the
  person level, the within factors carrying a person error stratum – so
  their terms are tested validly. Auto-detected from the fit's person
  identifier.

## Value

A list with `summary`, the compact reading of the analysis (one row per
item and group term with the uniform F, adjusted p, and partial
eta-squared – the term itself – and the non-uniform ones – the term
crossed with class interval – plus `uniform_DIF`, `nonuniform_DIF` and
`superseded` flags); `terms`, the complete per-item analysis of variance
table (term, df, sum of squares, mean square, F, partial eta-squared,
raw and adjusted p, significance, supersession, including the residual
row); and `tukey` (per item, term, and level comparison: difference, 95
per cent interval, and Tukey-adjusted p), plus the `alpha` and
adjustment used. Tukey comparisons are reported for significant,
non-superseded group terms except two-level main effects, where the F
test is already the only comparison. With `sizes = TRUE`, `sizes` holds
the logit DIF magnitudes per item, term, and level pair (two-level main
effects included, since the single difference is exactly the DIF size).

## Details

Sums of squares are sequential (factors in the order given, class
interval last), as is conventional for this residual diagnostic; with
markedly unbalanced groups the term order matters and the
factor-at-a-time
[`dif_anova`](https://drjoshmcgrane.github.io/rmt/reference/dif_anova.md)
is a useful cross-check.

## Examples

``` r
set.seed(1); n <- 800
d <- seq(-1.5, 1.5, length.out = 6)
g1 <- rep(c("a", "b"), each = n / 2)
g2 <- rep(c("x", "y"), times = n / 2)
sh <- matrix(0, n, 6); sh[g1 == "b", 2] <- 0.8
X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
colnames(X) <- paste0("I", 1:6)
fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))
dif_anova_factorial(fit)$terms
#>    item      term  df       sum_sq      mean_sq      F_value            p
#> 1    I1        g1   1 8.230232e-01 8.230232e-01 9.550190e-01 3.287793e-01
#> 2    I1        g2   1 1.360461e+00 1.360461e+00 1.578650e+00 2.093702e-01
#> 3    I1        ci   4 3.246510e+00 8.116275e-01 9.417956e-01 4.391025e-01
#> 4    I1     g1:g2   1 4.485547e-01 4.485547e-01 5.204936e-01 4.708697e-01
#> 5    I1     g1:ci   4 1.623572e+00 4.058930e-01 4.709897e-01 7.570630e-01
#> 6    I1     g2:ci   4 3.921920e-01 9.804799e-02 1.137728e-01 9.776780e-01
#> 7    I1  g1:g2:ci   4 4.632957e+00 1.158239e+00 1.343997e+00 2.520171e-01
#> 8    I1 Residuals 707 6.092836e+02 8.617873e-01           NA           NA
#> 9    I2        g1   1 2.021127e+01 2.021127e+01 2.437812e+01 9.885240e-07
#> 10   I2        g2   1 1.954500e+00 1.954500e+00 2.357448e+00 1.251334e-01
#> 11   I2        ci   4 6.599038e+00 1.649759e+00 1.989881e+00 9.431812e-02
#> 12   I2     g1:g2   1 5.295925e-02 5.295925e-02 6.387757e-02 8.005425e-01
#> 13   I2     g1:ci   4 8.528709e+00 2.132177e+00 2.571757e+00 3.676770e-02
#> 14   I2     g2:ci   4 7.411413e-01 1.852853e-01 2.234846e-01 9.253291e-01
#> 15   I2  g1:g2:ci   4 4.754322e+00 1.188581e+00 1.433624e+00 2.211040e-01
#> 16   I2 Residuals 707 5.861555e+02 8.290743e-01           NA           NA
#> 17   I3        g1   1 2.047455e+00 2.047455e+00 2.306874e+00 1.292499e-01
#> 18   I3        g2   1 3.345232e-01 3.345232e-01 3.769083e-01 5.394600e-01
#> 19   I3        ci   4 6.349981e+00 1.587495e+00 1.788636e+00 1.292444e-01
#> 20   I3     g1:g2   1 2.184829e-02 2.184829e-02 2.461653e-02 8.753712e-01
#> 21   I3     g1:ci   4 3.998454e+00 9.996136e-01 1.126268e+00 3.428888e-01
#> 22   I3     g2:ci   4 2.074367e+00 5.185917e-01 5.842989e-01 6.741067e-01
#> 23   I3  g1:g2:ci   4 5.785853e+00 1.446463e+00 1.629735e+00 1.649018e-01
#> 24   I3 Residuals 707 6.274944e+02 8.875452e-01           NA           NA
#> 25   I4        g1   1 1.478280e+00 1.478280e+00 1.709891e+00 1.914245e-01
#> 26   I4        g2   1 1.665960e+00 1.665960e+00 1.926976e+00 1.655269e-01
#> 27   I4        ci   4 9.397036e+00 2.349259e+00 2.717333e+00 2.888353e-02
#> 28   I4     g1:g2   1 1.101144e+00 1.101144e+00 1.273668e+00 2.594626e-01
#> 29   I4     g1:ci   4 3.068076e+00 7.670189e-01 8.871929e-01 4.710847e-01
#> 30   I4     g2:ci   4 5.837362e+00 1.459341e+00 1.687985e+00 1.508985e-01
#> 31   I4  g1:g2:ci   4 3.609162e+00 9.022905e-01 1.043658e+00 3.837160e-01
#> 32   I4 Residuals 707 6.112339e+02 8.645459e-01           NA           NA
#> 33   I5        g1   1 4.646516e-01 4.646516e-01 5.578515e-01 4.553755e-01
#> 34   I5        g2   1 9.339602e-01 9.339602e-01 1.121294e+00 2.900016e-01
#> 35   I5        ci   4 1.161085e+01 2.902713e+00 3.484940e+00 7.879883e-03
#> 36   I5     g1:g2   1 2.961540e-01 2.961540e-01 3.555566e-01 5.511753e-01
#> 37   I5     g1:ci   4 2.291815e+00 5.729538e-01 6.878769e-01 6.004857e-01
#> 38   I5     g2:ci   4 1.302151e+00 3.255378e-01 3.908342e-01 8.152800e-01
#> 39   I5  g1:g2:ci   4 3.516764e+00 8.791909e-01 1.055539e+00 3.776201e-01
#> 40   I5 Residuals 707 5.888820e+02 8.329307e-01           NA           NA
#> 41   I6        g1   1 5.716892e-01 5.716892e-01 6.714112e-01 4.128362e-01
#> 42   I6        g2   1 3.752787e+00 3.752787e+00 4.407401e+00 3.613727e-02
#> 43   I6        ci   4 1.024388e+00 2.560970e-01 3.007691e-01 8.774804e-01
#> 44   I6     g1:g2   1 3.497013e-04 3.497013e-04 4.107011e-04 9.838371e-01
#> 45   I6     g1:ci   4 1.153657e+01 2.884143e+00 3.387236e+00 9.315719e-03
#> 46   I6     g2:ci   4 4.770759e+00 1.192690e+00 1.400735e+00 2.320386e-01
#> 47   I6  g1:g2:ci   4 4.975015e+00 1.243754e+00 1.460707e+00 2.124431e-01
#> 48   I6 Residuals 707 6.019921e+02 8.514739e-01           NA           NA
#>    eta2_partial        p_adj significant superseded
#> 1  1.348983e-03 4.553755e-01       FALSE      FALSE
#> 2  2.227911e-03 3.140553e-01       FALSE      FALSE
#> 3  5.300164e-03 5.269230e-01       FALSE      FALSE
#> 4  7.356587e-04 9.838371e-01       FALSE      FALSE
#> 5  2.657641e-03 7.570630e-01       FALSE      FALSE
#> 6  6.432795e-04 9.776780e-01       FALSE      FALSE
#> 7  7.546558e-03 3.780257e-01       FALSE      FALSE
#> 8            NA           NA       FALSE      FALSE
#> 9  3.333176e-02 5.931144e-06        TRUE      FALSE
#> 10 3.323357e-03 3.140553e-01       FALSE      FALSE
#> 11 1.113283e-02 1.886362e-01       FALSE      FALSE
#> 12 9.034201e-05 9.838371e-01       FALSE      FALSE
#> 13 1.434158e-02 1.103031e-01       FALSE      FALSE
#> 14 1.262814e-03 9.776780e-01       FALSE      FALSE
#> 15 8.045766e-03 3.780257e-01       FALSE      FALSE
#> 16           NA           NA       FALSE      FALSE
#> 17 3.252293e-03 3.828490e-01       FALSE      FALSE
#> 18 5.328253e-04 5.394600e-01       FALSE      FALSE
#> 19 1.001820e-02 1.938666e-01       FALSE      FALSE
#> 20 3.481708e-05 9.838371e-01       FALSE      FALSE
#> 21 6.331749e-03 6.857776e-01       FALSE      FALSE
#> 22 3.294901e-03 9.776780e-01       FALSE      FALSE
#> 23 9.136322e-03 3.780257e-01       FALSE      FALSE
#> 24           NA           NA       FALSE      FALSE
#> 25 2.412682e-03 3.828490e-01       FALSE      FALSE
#> 26 2.718159e-03 3.140553e-01       FALSE      FALSE
#> 27 1.514110e-02 8.665058e-02       FALSE      FALSE
#> 28 1.798270e-03 9.838371e-01       FALSE      FALSE
#> 29 4.994409e-03 7.066271e-01       FALSE      FALSE
#> 30 9.459785e-03 6.961158e-01       FALSE      FALSE
#> 31 5.870054e-03 3.837160e-01       FALSE      FALSE
#> 32           NA           NA       FALSE      FALSE
#> 33 7.884181e-04 4.553755e-01       FALSE      FALSE
#> 34 1.583477e-03 3.480020e-01       FALSE      FALSE
#> 35 1.933554e-02 4.727930e-02        TRUE      FALSE
#> 36 5.026561e-04 9.838371e-01       FALSE      FALSE
#> 37 3.876720e-03 7.205829e-01       FALSE      FALSE
#> 38 2.206347e-03 9.776780e-01       FALSE      FALSE
#> 39 5.936481e-03 3.837160e-01       FALSE      FALSE
#> 40           NA           NA       FALSE      FALSE
#> 41 9.487613e-04 4.553755e-01       FALSE      FALSE
#> 42 6.195326e-03 2.168236e-01       FALSE      FALSE
#> 43 1.698773e-03 8.774804e-01       FALSE      FALSE
#> 44 5.809064e-07 9.838371e-01       FALSE      FALSE
#> 45 1.880364e-02 5.589431e-02       FALSE      FALSE
#> 46 7.862642e-03 6.961158e-01       FALSE      FALSE
#> 47 8.196516e-03 3.780257e-01       FALSE      FALSE
#> 48           NA           NA       FALSE      FALSE
```
