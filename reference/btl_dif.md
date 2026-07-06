# DIF analysis for paired comparisons

Tests whether objects function differently for identifiable groups of
judges. For each object: (i) the standardised residuals of its
comparisons, oriented to the object, are analysed by judge group crossed
with opponent-strength bands – a group main effect is uniform DIF and a
group-by-band interaction non-uniform DIF, mirroring
[`dif_anova`](https://drjoshmcgrane.github.io/rmt/reference/dif_anova.md);
and (ii) the object is resolved into one copy per judge group inside a
joint refit and the differences between the resolved locations are
reported in logits with judge-clustered Wald tests, familywise
adjustment, and the practical-significance flag, mirroring
[`dif_size`](https://drjoshmcgrane.github.io/rmt/reference/dif_size.md).

## Usage

``` r
btl_dif(
  fit,
  groups,
  objects = NULL,
  p_adjust = "holm",
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

- groups:

  Judge grouping: one value per row of `fit$comparisons`, or a vector
  named by judge.

- objects:

  Objects to test; all by default.

- p_adjust:

  Familywise adjustment over all pairwise location comparisons; the
  ANOVA probabilities are adjusted across objects by Benjamini-Hochberg
  within each term, as in
  [`dif_anova`](https://drjoshmcgrane.github.io/rmt/reference/dif_anova.md).

- alpha:

  Significance level for adjusted probabilities.

- flag_logits:

  Absolute resolved difference flagged as practically significant.

- min_n:

  Group levels with fewer comparisons involving the object are dropped
  from its resolution, with a note.

- maxit, tol:

  Newton controls for the resolution refits.

## Value

A list of class `"rmt_btl_dif"`: `anova` (per object: uniform and
non-uniform F, raw and adjusted p, flags), `levels` (resolved location
and SE per object and group), `sizes` (per object and group pair:
difference in logits, SE, z, adjusted p, significance and practical
flags), and `notes`.

## References

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
btl_dif(f, groups = grp, objects = "C")
#> DIF for paired comparisons: 1 object(s) by judge group
#> Residual ANOVA (uniform = group; non-uniform = group x opponent band; BH across objects)
#>  object   n F_uniform p_uniform_adj uniform_DIF F_nonuniform p_nonuniform_adj
#>       C 180    14.495       < 0.001           *        0.081            0.922
#>  nonuniform_DIF
#>                
#> 
#> Resolved locations (logits; holm over 1 comparison(s); practical 0.50)
#>  object level_a level_b difference    se     z   p_adj significant practical
#>       C      g1      g2      1.404 0.271 5.186 < 0.001           *         *
```
