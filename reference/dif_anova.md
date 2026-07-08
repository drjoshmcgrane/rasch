# Differential item functioning by residual analysis of variance

For each item the standardised residuals are analysed by the nominated
person factor(s) crossed with the trait class interval. A term not
involving the class interval is uniform DIF; a term crossing it is
non-uniform DIF (Hagquist and Marais 2019, ch. 16). With one factor this
is a one-way analysis, `z ~ g * ci`. With several factors they are
modelled jointly – the statistically correct treatment, rather than one
factor at a time – with main effects by default
(`z ~ (f1 + f2 + ...) * ci`); set `effects = "factorial"` to add the
factor-by-factor interactions (`z ~ (f1 * f2 * ...) * ci`). When
interactions are fitted, a significant one supersedes the lower-order
terms built from its variables, recorded in the `superseded` column;
interpret the highest-order significant terms.

## Usage

``` r
dif_anova(
  fit,
  factors = NULL,
  n_groups = NULL,
  p_adjust = "BH",
  alpha = 0.05,
  effects = c("main", "factorial"),
  sizes = FALSE,
  id = NULL,
  within = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- factors:

  A vector (one factor), a data frame of person factors, or a character
  vector naming factor columns nominated in the fit. Defaults to every
  factor stored in the fit.

- n_groups:

  Number of trait class intervals. By default set from the smallest
  factor-combination cell so every interval-by-cell count keeps about 30
  expected responses (between 2 and 10 intervals); the value used is
  returned in `n_groups`.

- p_adjust:

  Multiplicity adjustment across items within each term; default `"BH"`.

- alpha:

  Significance level applied to the adjusted probabilities.

- effects:

  `"main"` (default) models several factors additively (each factor's
  main effect and its class-interval interaction, but no
  factor-by-factor terms); `"factorial"` also crosses the factors with
  each other. Immaterial with a single factor.

- sizes:

  Also compute DIF magnitudes in logits
  ([`dif_size`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md))
  for every significant, non-superseded group term: the item is resolved
  by the term's levels (interaction terms by their cells) and all
  pairwise location differences are returned with Holm familywise
  adjustment and the practical-significance flag. Each size involves a
  re-analysis, so this costs one refit per flagged item-term.

- id, within:

  Person identifier and within-subject factor names for stacked
  repeated-measures designs. When any factor is within-subject the model
  becomes a mixed (split-plot) analysis of variance – the class interval
  taken at the person level, the within factors carrying a person error
  stratum – so their terms are tested validly. Auto-detected from the
  fit's person identifier.

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

Probabilities are adjusted across items within each term
(Benjamini-Hochberg by default). Tukey HSD comparisons are returned for
each significant, non-superseded group term. Sums of squares are
sequential (factors in the order given, class interval last).

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
dif_anova(fit)$summary
#>    item term  F_uniform    p_uniform p_uniform_adj eta2_uniform uniform_DIF
#> 1    I1   g1  0.9538018 3.290852e-01  4.552104e-01 0.0013378171       FALSE
#> 2    I1   g2  1.5766382 2.096575e-01  3.144863e-01 0.0022094868       FALSE
#> 3    I2   g1 24.3495106 1.001290e-06  6.007742e-06 0.0330678710        TRUE
#> 4    I2   g2  2.3546818 1.253516e-01  3.144863e-01 0.0032962362       FALSE
#> 5    I3   g1  2.3019432 1.296561e-01  3.835307e-01 0.0032226473       FALSE
#> 6    I3   g2  0.3761028 5.398922e-01  5.398922e-01 0.0005279554       FALSE
#> 7    I4   g1  1.7072431 1.917653e-01  3.835307e-01 0.0023920775       FALSE
#> 8    I4   g2  1.9239918 1.658513e-01  3.144863e-01 0.0026949533       FALSE
#> 9    I5   g1  0.5582562 4.552104e-01  4.552104e-01 0.0007834534       FALSE
#> 10   I5   g2  1.1221076 2.898243e-01  3.477892e-01 0.0015735139       FALSE
#> 11   I6   g1  0.6706080 4.131137e-01  4.552104e-01 0.0009409790       FALSE
#> 12   I6   g2  4.4021284 3.624594e-02  2.174756e-01 0.0061447729       FALSE
#>    F_nonuniform p_nonuniform p_nonuniform_adj eta2_nonuniform nonuniform_DIF
#> 1     0.4640075  0.762182811       0.76218281    0.0026000063          FALSE
#> 2     0.1169713  0.976504536       0.97650454    0.0006567108          FALSE
#> 3     2.5582202  0.037591977       0.11277593    0.0141683953          FALSE
#> 4     0.2242190  0.924909281       0.97650454    0.0012580724          FALSE
#> 5     1.1295688  0.341325851       0.68265170    0.0063058758          FALSE
#> 6     0.5819722  0.675793611       0.97650454    0.0032588521          FALSE
#> 7     0.8257136  0.508957578       0.69635131    0.0046174207          FALSE
#> 8     1.5816798  0.177323408       0.69062571    0.0088075791          FALSE
#> 9     0.7172319  0.580292761       0.69635131    0.0040132219          FALSE
#> 10    0.3857592  0.818903711       0.97650454    0.0021624998          FALSE
#> 11    3.3737483  0.009529991       0.05717994    0.0186010839          FALSE
#> 12    1.4061183  0.230208570       0.69062571    0.0078376275          FALSE
#>    superseded
#> 1       FALSE
#> 2       FALSE
#> 3       FALSE
#> 4       FALSE
#> 5       FALSE
#> 6       FALSE
#> 7       FALSE
#> 8       FALSE
#> 9       FALSE
#> 10      FALSE
#> 11      FALSE
#> 12      FALSE
```
