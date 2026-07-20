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
  repeated-measures designs, auto-detected from the fit's person
  identifier. Whenever ids repeat, PERSONS are the units of analysis:
  residuals are aggregated to one mean per person (per within-subject
  cell), so duplicated or stacked observations cannot manufacture
  information, and the class interval is taken at the person level.
  Between-person terms are tested with order-invariant Type II sums of
  squares – every term adjusted for every term not containing it, the
  class interval always among them, so entry order cannot decide which
  correlated factor absorbs shared or trait variance. Within-person
  terms are tested on the person-by-cell means through orthonormal
  contrasts with the Greenhouse-Geisser epsilon correction (Maxwell and
  Delaney 2004), so multi-level within factors are valid without the
  sphericity assumption; persons missing any within cell are dropped
  from the within-stratum tests. A factor that varies within persons
  must be declared (or auto-detected) as within-subject; treating it as
  between-subjects is refused. The BH adjustment is applied across items
  separately within each term: each term is read as its own prespecified
  family, not as one pooled screen across all terms.

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
each significant, non-superseded group term. Sums of squares are Type II
(each term adjusted for every term not containing it, the class interval
always among them), so results do not depend on the order factors are
given.

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
#> 1    I1   g1  0.6121334 4.342458e-01  4.342458e-01 0.0008589994       FALSE
#> 2    I1   g2  1.3413518 2.471841e-01  2.966209e-01 0.0018803785       FALSE
#> 3    I2   g1 22.5995342 2.415245e-06  1.449147e-05 0.0307644276        TRUE
#> 4    I2   g2  2.9308828 8.733541e-02  2.620062e-01 0.0040995331       FALSE
#> 5    I3   g1  2.6226689 1.057899e-01  2.632076e-01 0.0036700052       FALSE
#> 6    I3   g2  0.3236896 5.695781e-01  5.695781e-01 0.0004544136       FALSE
#> 7    I4   g1  2.2787190 1.316038e-01  2.632076e-01 0.0031902378       FALSE
#> 8    I4   g2  1.4279983 2.324893e-01  2.966209e-01 0.0020016012       FALSE
#> 9    I5   g1  0.8332442 3.616451e-01  4.342458e-01 0.0011689189       FALSE
#> 10   I5   g2  1.3836032 2.398815e-01  2.966209e-01 0.0019394940       FALSE
#> 11   I6   g1  0.6800322 4.098518e-01  4.342458e-01 0.0009541900       FALSE
#> 12   I6   g2  3.9539511 4.714388e-02  2.620062e-01 0.0055226333       FALSE
#>    F_nonuniform p_nonuniform p_nonuniform_adj eta2_nonuniform nonuniform_DIF
#> 1     0.4714642  0.756715172       0.75671517    0.0026416781          FALSE
#> 2     0.1169713  0.976504536       0.97650454    0.0006567108          FALSE
#> 3     2.5711448  0.036798473       0.11039542    0.0142389570          FALSE
#> 4     0.2242190  0.924909281       0.97650454    0.0012580724          FALSE
#> 5     1.2157323  0.302694888       0.60538978    0.0067836248          FALSE
#> 6     0.5819722  0.675793611       0.97650454    0.0032588521          FALSE
#> 7     0.8410808  0.499308707       0.66662317    0.0047029510          FALSE
#> 8     1.5816798  0.177323408       0.69062571    0.0088075791          FALSE
#> 9     0.7539858  0.555519305       0.66662317    0.0042180085          FALSE
#> 10    0.3857592  0.818903711       0.97650454    0.0021624998          FALSE
#> 11    3.5440568  0.007116541       0.04269925    0.0195217452           TRUE
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
