# Rank analysis with the Plackett-Luce model

Calibrates objects from rankings. Each ranking is read as a sequence of
choices, the object ranked first chosen from all the objects in the
ranking, the object ranked second from those remaining, and so on, with
\\P(j \mathrm{\\ chosen}) = \exp(\beta_j) / \sum_k \exp(\beta_k)\\ over
the objects still to be placed (Luce 1959; Plackett 1975). A ranking of
two objects is a paired comparison, and the model then equals
[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md): the
locations are on the same logit scale and, when every ranking is a pair,
the two functions agree.

## Usage

``` r
pl(
  data,
  ranking = "ranking",
  object = "object",
  rank = "rank",
  judge = NULL,
  anchors = NULL,
  ties = c("drop", "error"),
  se = c("sandwich", "model"),
  split = c("first", "half"),
  maxit = 100,
  tol = 1e-08
)
```

## Arguments

- data:

  A data frame in long format: one row per object per ranking.

- ranking, object, rank:

  Names of the columns holding the ranking identifier, the object and
  the rank (1 = highest). A missing rank marks an object present in the
  ranking but unranked.

- judge:

  Optional name of a judge column. Rankings are then clustered by judge
  for the sandwich covariance and a judge fit table is reported.

- anchors:

  Optional named numeric vector of fixed object locations.

- ties:

  How to treat a ranking with tied ranks: `"drop"` the ranking with a
  note, or `"error"`.

- se:

  `"sandwich"` (clustered Godambe errors, withheld when the cluster
  design does not support them) or `"model"` (inverse observed
  information of the Plackett-Luce likelihood).

- split:

  How the choices are grouped for the invariance check: `"first"` (the
  first choice of each ranking against the later choices) or `"half"`
  (the early half of each ranking against the late half).

- maxit, tol:

  Newton-Raphson controls.

## Value

A `"rasch_pl"` object with `objects` (location, se, rankings, stages,
chosen, infit, outfit, fit residual, extreme flag), `judges` (when a
judge column is given), `rankings` (one row per ranking with its
log-likelihood and surprise `z`), `reversal` (the reversal check, or
`NULL` when no ranking is complete with three or more objects),
`invariance` (the invariance check: `groups`, `lr`, `df`, `p` and the
per-object contrasts in `objects`, or `NULL` when withheld), `osi`
(object separation), `loglik`, `cov_beta`, `converged`, `iterations`,
`n_rankings`, `n_stages`, `se_type`, `anchors`, `notes` and the `call`.

## Details

The locations are maximum likelihood estimates with a sum-zero origin,
or at the values in `anchors`. Standard errors are Godambe sandwich
errors clustered by judge when a judge column is given and by ranking
otherwise, subject to the same conditions as
[`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md): with
fewer than ten judges, fewer than eight effective judges, no residual
effective cluster degrees of freedom, or a rank-deficient score
covariance, errors are withheld with a note. Judge-clustered errors
include the \\G/(G-1)\\ correction, where \\G\\ is the judge count.
`se = "model"` instead reports the information-based errors of the
Plackett-Luce likelihood, which are valid when the model holds and the
rankings are independent, and are available whatever the design.

Fit is read from the residuals of the stage choices, \\y - p\\ for each
object in each choice set, pooled by object and by judge into the infit,
outfit and standardised fit residuals of
[`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md). Each
ranking also receives a surprise statistic: its log-likelihood
standardised against the mean and variance the fitted model implies for
it, so a ranking the consensus makes very unlikely stands out.

An object that is chosen at every stage it appears in, or never chosen
at any stage with an alternative, and is not anchored, has no finite
location (the undefeated or winless case of a paired comparison). It is
set aside and reported with `extreme = TRUE` at an extrapolated
location, its score moved half a choice inside the boundary against the
calibrated objects. A design in which a group of objects is never ranked
below an object outside the group cannot place that group and is
refused.

Rankings may be partial in two ways. A ranking that lists only some of
the objects is a choice among those. An object listed with a missing
`rank` is present but unranked: it is in every choice set of that
ranking and never chosen, which is the top-\\k\\ design. Tied ranks
within a ranking are not modelled; a ranking with a tie is dropped
(`ties = "drop"`) or refused.

**Reversal check.** The sequential reading is not symmetric: read from
the bottom up, with the worst object chosen first, the same rankings
have a different likelihood and in general different locations. If the
judges worked best-first the forward reading should fit at least as well
as the reversed one. The check fits both readings to the rankings that
are complete and hold at least three objects, reports the Vuong (1989)
statistic for the two non-nested models (positive favours best-first)
and the correlation and largest difference between the two sets of
locations. It is a check on the ranking process, not on the objects.
When judges are recorded, the likelihood differences are clustered by
judge and referred to a t distribution on judges minus one degrees of
freedom. Inference is withheld without sufficient cluster support. The
reference is asymptotic and assumes the two readings are
distinguishable; it does not test distinguishability at their common
equal-worth model.

**Invariance check.** Luce's choice axiom, which the model rests on,
says that once an object is chosen the remaining objects compete on the
same scale as before. The locations estimated from the first choice of
each ranking should then agree with the locations estimated from the
later choices, in the way item locations should agree across class
intervals in Andersen's (1973) test. The check fits the model to the two
groups of choices (`split = "first"`: first choice against later
choices; `split = "half"`: the early against the late half of each
ranking), each with its own extreme objects set aside and its own
connectivity required, and reports the likelihood ratio against the
pooled fit on the stages both groups kept, with degrees of freedom the
free parameters of the groups less those of the pooled fit. Each object
calibrated in both groups also receives a Wald contrast of its two
locations, the groups centred on the objects they share unless both are
on the anchor scale, with Holm-adjusted p values. The contrasts use the
requested covariance as the main fit, including the covariance between
early and later estimates from the same ranking or judge. Each group
must support that covariance; ordinary errors never replace withheld
sandwich errors. Judge-clustered contrasts use a t reference, and the
likelihood-ratio probability is withheld because it assumes independent
rankings. With `se = "model"`, independence is assumed and the
likelihood-ratio reference is chi-square. The check is withheld with a
note when every ranking has a single informative choice (a design of
pairs), when a group cannot calibrate the objects on its own, or when
fewer than two objects are calibrated in both groups. It answers a
different question from the reversal check: a set of rankings may fit
better reversed and still be invariant across positions, or fail
invariance in either direction. With few distinct rankings, the
contrasts can be conservative or unavailable even when the main
calibration has usable standard errors.

## References

Andersen, E. B. (1973). A goodness of fit test for the Rasch model.
Psychometrika, 38, 123–140.

Luce, R. D. (1959). Individual Choice Behavior. Wiley.

Plackett, R. L. (1975). The analysis of permutations. Applied
Statistics, 24, 193–202.

Hunter, D. R. (2004). MM algorithms for generalized Bradley-Terry
models. Annals of Statistics, 32, 384–406.

Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
non-nested hypotheses. Econometrica, 57, 307–333.

## See also

[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) for
paired comparisons,
[`rasch_cj`](https://drjoshmcgrane.github.io/rasch/reference/rasch_cj.md)
to combine rankings with comparisons and item responses,
[`plot_pl`](https://drjoshmcgrane.github.io/rasch/reference/plot_pl.md).

## Examples

``` r
set.seed(1)
beta <- c(A = -1.5, B = -0.5, C = 0, D = 0.5, E = 1.5)
rk <- do.call(rbind, lapply(1:60, function(r) {
  rem <- names(beta); ord <- character(0)
  while (length(rem) > 1) {
    pick <- sample(rem, 1, prob = exp(beta[rem]))
    ord <- c(ord, pick); rem <- setdiff(rem, pick)
  }
  data.frame(ranking = r, object = c(ord, rem), rank = 1:5)
}))
fit <- pl(rk)
fit
#> Plackett-Luce rank analysis: 5 objects, 60 rankings of 5 objects
#> Maximum likelihood: converged in 5 iterations; sandwich SEs clustered by ranking
#> Object separation index 0.979
#> Reversal check on 60 complete rankings: best-first log-likelihood -217.06, worst-first -226.30; Vuong z = 1.77, p = 0.076; location correlation 0.989
#> Invariance check (first choice vs later choices): LR = 5.24 on 4 df, p = 0.264; objects moving (Holm p < 0.05): none
#>  object location    se rankings chosen fit_resid extreme
#>       A   -1.617 0.201       60     19    -0.311        
#>       B   -0.491 0.165       60     48     0.435        
#>       C    0.061 0.136       60     56    -0.791        
#>       D    0.568 0.173       60     57     1.068        
#>       E    1.479 0.159       60     60    -0.273        
#> Rankings beyond surprise z 2.5: 1 (28)
fit$reversal$z
#> [1] 1.773556
fit$invariance$objects
#>  object  first  later difference    se      z     p p_adj
#>       A -1.963 -1.659     -0.303 0.853 -0.356 0.722 1.000
#>       B -0.353 -0.571      0.218 0.467  0.467 0.641 1.000
#>       C -0.017  0.031     -0.047 0.414 -0.115 0.909 1.000
#>       D  0.928  0.290      0.638 0.349  1.831 0.067 0.335
#>       E  1.405  1.910     -0.506 0.495 -1.022 0.307 1.000
```
