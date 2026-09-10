# DIF differences between factor levels

Resolves an item by one or more person factors and compares the
resulting locations. Several factors in `by` give pairwise comparisons
between their joint cells and can be used to quantify an interaction.

## Usage

``` r
dif_size(
  fit,
  item,
  by,
  p_adjust = "holm",
  alpha = 0.05,
  flag_logits = 0.5,
  min_n = 20
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md) or
  [`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md).
  EFRM fits are excluded because an ordinary split refit would discard
  their frame units.

- item:

  Item name or index.

- by:

  One or more person-factor names nominated in the fit (several names
  give interaction cells), or a grouping vector/data frame with one
  entry per person.

- p_adjust:

  Adjustment over the pairwise comparisons. The default `"holm"`
  controls familywise error; use `"BH"` only for false-discovery-rate
  screening. `"none"` leaves probabilities unadjusted.

- alpha:

  Significance level for the adjusted probabilities.

- flag_logits:

  Absolute difference flagged as practically significant.

- min_n:

  Levels with fewer distinct responders to the item are dropped (their
  resolved locations would be too unstable to compare), with a note.
  When identifiers repeat, response rows from one person count once
  within each level.

## Value

A list of class `"rasch_dif_size"`. `levels` contains the resolved
location, standard error and sample size for each level. `pairs`
contains logit differences, Wald `t` statistics, confidence intervals,
raw and adjusted probabilities, and practical flags. For dichotomous
items it also contains `ets`, together with the raw and adjusted
probabilities for exceeding the ETS A boundary; for polytomous items it
contains the descriptive `signed_area`. `df` gives the reference degrees
of freedom: infinite for independent response rows and the independent
person-cluster count minus one for a supported repeated-person
calibration. Sampling-uncertainty fields are `NA` when the
resolved-location covariance cannot support Wald inference. The ETS
category is also `NA` if a probability needed to classify the contrast
is unavailable or its standard error is zero.

## Details

Let \\\delta_i\\ contain the resolved locations of item \\i\\, and let
\\\mathbf{c}\_{ab}\\ place 1 on level \\a\\, -1 on level \\b\\, and zero
elsewhere. The reported difference and its standard error are
\$\$\Delta\_{i,ab}=\mathbf{c}\_{ab}^{\mathsf T} \delta_i,\$\$
\$\$\operatorname{SE}(\Delta\_{i,ab})= \sqrt{\mathbf{c}\_{ab}^{\mathsf
T}\mathbf{V}\_i \mathbf{c}\_{ab}},\$\$ where \\\mathbf{V}\_i\\ is the
full covariance of the resolved locations. Wald probabilities are
adjusted over the pairwise family. A comparison with withheld inference
remains in that declared family.

When person identifiers repeat, the resolved-location covariance uses
the person-clustered calibration sandwich and a t reference with the
number of independent person clusters minus one degree of freedom. The
same reference is used for confidence intervals and the ETS
interval-null probability. This permits inference for the logit
difference while allowing response rows from the same person to be
dependent. For a planned within-person question,
[`dif_contrasts`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md)
remains preferable because it tests the nominated contrast directly from
person-level residual contrasts. Inference is withheld if the resolved-
location covariance is unavailable or not positive semidefinite.

## Magnitude conventions

For dichotomous items, `ets` applies the ETS A, B and C rules to the
itemwise comparison. On the logit scale the magnitude cut-points are
\\1/2.35=0.426\\ and \\1.5/2.35=0.638\\. Category A also includes an
adjusted test that is not significant. Category C requires a magnitude
of at least 0.638 and rejection of the interval null \\\|\Delta\|\leq
0.426\\; B is the remainder. Both probabilities are adjusted by
`p_adjust` over the requested pairwise family.

For a partial credit item with \\m_i\\ thresholds, the signed area
between the two expected-score curves has the closed form
\$\$SA\_{ab}=\int\\E_b(X\mid\theta)-E_a(X\mid\theta)\\\\d\theta
=\sum\_{k=1}^{m_i}(\delta\_{iak}-\delta\_{ibk})
=m_i(\beta\_{ia}-\beta\_{ib}).\$\$ This is returned as `signed_area`; a
positive value means that level `a` has the harder resolved item. It is
descriptive and is not given an A/B/C category: score-metric
classifications for polytomous DIF are not interchangeable with a PCM
logit difference. For pooled MFRM items, the areas use the same
precision weight for a facet cell in every group. A comparison is
withheld when the groups do not support the same observed response
categories.

## References

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. Scandinavian Journal of Statistics, 6(2), 65–70.

Zieky, M. (1993). Practical questions in the use of DIF statistics in
item development. In P. W. Holland and H. Wainer (eds), Differential
Item Functioning (pp. 337–364). Erlbaum.

Linacre, J. M. and Wright, B. D. (1989). Mantel-Haenszel DIF and PROX
are equivalent! Rasch Measurement Transactions, 3(2), 51–53.

Cohen, A. S., Kim, S.-H. and Baker, F. B. (1993). Detection of
differential item functioning in the graded response model. Applied
Psychological Measurement, 17(4), 335–350.

Raju, N. S. (1988). The area between two item characteristic curves.
Psychometrika, 53(4), 495–502.

## See also

[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
and
[`dif_contrasts`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md).

## Examples

``` r
set.seed(1); n <- 600
d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
colnames(X) <- paste0("I", 1:8)
fit <- rasch(data.frame(X, grp = g), factors = "grp")
dif_size(fit, "I3", by = "grp")
#> DIF size for I3 by grp (resolved locations, logits)
#>  level location    se weak   n
#>      a   -0.890 0.133    0 300
#>      b    0.018 0.126    0 300
#>  level_a level_b difference    se      t  df       p   p_adj  lower  upper
#>        a       b     -0.907 0.204 -4.441 Inf < 0.001 < 0.001 -1.308 -0.507
#>  significant practical p_beyond_A p_beyond_A_adj ets signed_area
#>            *   >= 0.50      0.009          0.009  C-            
#> p adjusted by holm over 1 pairwise comparison(s); practical criterion 0.50 logits
```
