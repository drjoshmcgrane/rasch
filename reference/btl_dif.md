# DIF analysis for paired comparisons

Tests whether objects function differently for identifiable groups of
judges. One judge factor is analysed on its own; several factors are
modelled jointly – with main effects by default and factor-by-factor
interactions optional – exactly as
[`dif_anova`](https://drjoshmcgrane.github.io/rmt/reference/dif_anova.md)
treats person factors. For each object the standardised residuals of its
comparisons, oriented to the object, are analysed by the judge factor(s)
crossed with opponent-strength bands: a term is uniform DIF, its
crossing with the band non-uniform DIF, and a significant higher-order
group term supersedes the lower-order group terms built from a subset of
its factors. Each term flagged for uniform DIF and not superseded is
then resolved – the object split into one copy per cell of the term's
factors inside a joint refit – and the differences between the resolved
locations reported in logits with judge-clustered Wald tests and the
practical-significance flag, mirroring
[`dif_size`](https://drjoshmcgrane.github.io/rmt/reference/dif_size.md).
Fits with within-judge dependence effects (`order`) keep those effects
in the residual moments and in the refits, so dependence is not mistaken
for judge-group DIF; count-weighted comparisons enter all tests with
their weights.

## Usage

``` r
btl_dif(
  fit,
  factors,
  objects = NULL,
  effects = c("main", "factorial"),
  p_adjust = "BH",
  alpha = 0.05,
  flag_logits = 0.5,
  min_n = 20,
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- fit:

  An object from
  [`btl`](https://drjoshmcgrane.github.io/rmt/reference/btl.md).

- factors:

  A judge factor, or a named list of them, each either one value per row
  of `fit$comparisons` or a vector named by judge.

- objects:

  Objects to test; all by default.

- effects:

  `"main"` (default) models several factors additively (each factor's
  main effect and its band interaction); `"factorial"` also crosses the
  factors with one another.

- p_adjust:

  Multiplicity adjustment across objects within each term; the
  resolved-size probabilities are adjusted in one pool over all objects,
  terms, and cell pairs.

- alpha:

  Significance level for adjusted probabilities.

- flag_logits:

  Absolute resolved difference flagged as practically significant.

- min_n:

  Term cells with fewer comparisons involving the object are dropped
  from its resolution, with a note.

- maxit, tol:

  Newton controls for the resolution refits.

## Value

A list of class `"rmt_btl_dif"`: `summary` (one row per object and group
term with the uniform F, adjusted p and partial eta-squared – the term
itself – the non-uniform ones – the term crossed with the opponent band
– plus `uniform_DIF`, `nonuniform_DIF` and `superseded` flags); `terms`
(the full per-object analysis-of-variance table); `levels` (resolved
location and SE per object, term and cell); `sizes` (per object, term
and cell pair: difference in logits, SE, z, adjusted p, significance and
practical flags); `effects`, `factors`, and `notes`.

## Details

Each object is resolved against the other objects' common locations.
When several objects carry real DIF, resolving them one at a time can
spread a large effect onto clean objects as compensating,
opposite-signed artificial DIF (Andrich & Hagquist 2012, 2015); read
large flags on several objects together with that hazard in mind, and
prefer resolving the largest effect first and re-running.

## References

Andrich, D., & Hagquist, C. (2012). Real and artificial differential
item functioning. *Journal of Educational and Behavioral Statistics*,
37(3), 387-416.

Dittrich, R., Hatzinger, R., & Katzenbeisser, W. (1998). Modelling the
effect of subject-specific covariates in paired comparison studies with
an application to university rankings. *Journal of the Royal Statistical
Society C*, 47(4), 511-525.

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 60), b = rep(pr[, 2], each = 60),
                judge = sample(sprintf("J%02d", 1:12), 360, TRUE))
shift <- ifelse(d$judge %in% sprintf("J%02d", 1:6) & d$a == "C", 0.9,
         ifelse(d$judge %in% sprintf("J%02d", 1:6) & d$b == "C", -0.9, 0))
p <- plogis(beta[d$a] - beta[d$b] + shift)
d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
f <- btl(d, "a", "b", winner = "win", judge = "judge")
grp <- setNames(rep(c("g1", "g2"), each = 6), sprintf("J%02d", 1:12))
btl_dif(f, grp, objects = "C")
#> DIF for paired comparisons: 1 factor(s) [group], main effects
#> Residual ANOVA per object and term (uniform = term; non-uniform = term x opponent band)
#>  object  term F_uniform p_uniform_adj uniform_DIF F_nonuniform p_nonuniform_adj
#>       C group    14.495       < 0.001           *        0.081            0.922
#>  nonuniform_DIF
#>                
#> 
#> Resolved locations (logits; BH over 1 comparison(s); practical 0.50)
#>  object  term level_a level_b difference    se     z   p_adj significant
#>       C group      g1      g2      1.404 0.271 5.186 < 0.001           *
#>  practical
#>          *
```
