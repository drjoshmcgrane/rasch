# Fit a many-facet Rasch model

Estimates the many-facet Rasch model (Linacre 1989) for long-format data
in which each row is one scored response carrying a person, an item, a
score, and one or more facet levels (for example the rater). Every
item-by-facet combination becomes a virtual item whose thresholds are
the item's thresholds shifted by the facet severities, and the whole
structure is estimated in one pass of the pairwise conditional
likelihood, in which the person parameter cancels. Facet severities are
reported with standard errors and pooled fit statistics; the returned
object is also a full
[`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md) fit at
the virtual-item level, so every diagnostic table and plot in the
package applies to it.

## Usage

``` r
rasch_mfrm(
  data,
  person,
  item = NULL,
  score = NULL,
  facets,
  items = NULL,
  n_groups = NULL,
  adjust_N = NA,
  na_codes = -1,
  interaction = NULL,
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- data:

  Long-format data frame.

- person:

  Name of the person identifier column.

- item:

  Name of the item column.

- score:

  Name of the integer score column (categories from 0; gaps are
  collapsed per item with a note).

- facets:

  Character vector naming one or more facet columns (for example a rater
  column).

- items:

  Optional character vector of item score columns for data in wide
  format: one row per person-by-facet combination (for example one row
  per script per rater) with one column per item or criterion. The long
  form (`item` + `score`) remains available for data where the facet
  varies within items.

- n_groups:

  Number of class intervals for the item-trait chi-square; `NULL` (the
  default) applies the class-interval rule of Andrich and Marais (2019,
  ch. 15) (at least 50 non-extreme persons per interval, at most 10
  intervals, at least 2).

- adjust_N:

  Optional reference sample size for the chi-square.

- na_codes:

  Score values to read as missing (default `-1`); any negative score is
  also treated as missing.

- interaction:

  Optional name of one facet to interact with the items (interactive
  facet mode). Adds item-by-facet terms `gamma[item, level]` with double
  sum-to-zero constraints on top of the additive severities, so each
  level may be more or less severe on particular items; estimates are
  returned in `interaction_effects`. The interactive model remains in
  the Rasch class (all discriminations equal one and the parameters are
  additive), but a significant interaction qualifies specific
  objectivity in practice: comparisons of the interacting facet's levels
  become item-dependent, which is itself the substantive finding.

- maxit, tol:

  Newton-Raphson iteration cap and convergence tolerance.

## Value

An object of classes `"rasch_mfrm"` and `"rasch"`. In addition to every
component of a
[`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md) fit
(computed over the virtual items), it carries `facet_effects` (per
facet: level, severity, standard error, observation count, pooled fit),
`item_effects` (underlying item locations and pooled fit),
`item_thresholds` (the structural `delta_ik` with standard errors), and
`facet_spec`. Two fit residuals are reported per facet level and per
underlying item. `fit_resid` is the facet-margin statistic of the
published three-facet fit tables (Andrich and Marais 2019, ch. 26 and
app. C), the mean of the constituent virtual items' fit residuals; it
weighs each virtual item equally, so an erratic level shows the average
of its per-item misfit. `fit_resid_pooled` is the log-of-mean-square
statistic summed over the margin's observed cells of non-extreme
persons, with its degrees of freedom in `df_fit`; it weighs each
response equally and is the more powerful statistic when misfit is
spread evenly over the level's cells.

## Examples

``` r
set.seed(1)
simP <- function(th, tau) { x <- 0:length(tau); p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
persons <- sprintf("P%03d", 1:120); raters <- paste0("R", 1:4)
th <- setNames(rnorm(120, 0, 1.3), persons)
rho <- setNames(c(-0.6, -0.2, 0.2, 0.6), raters)
tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4))
d <- expand.grid(person = persons, item = names(tau), rater = raters,
                 stringsAsFactors = FALSE)
d$score <- mapply(function(p, i, r)
  sample(0:2, 1, prob = simP(th[p], tau[[i]] + rho[r])), d$person, d$item, d$rater)
fit <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                  facets = "rater")
fit$facet_effects$rater
#>    level   severity         se   n  infit_ms outfit_ms    fit_resid
#> R1    R1 -0.6741794 0.08866656 354 1.0748808 1.0555789  0.352198069
#> R2    R2 -0.1984195 0.07726631 354 1.0559605 1.0389336  0.337676916
#> R3    R3  0.1850839 0.07784122 354 1.0009494 0.9972743 -0.007416974
#> R4    R4  0.6875150 0.07293532 354 0.9666351 0.9568306 -0.336987303
#>    fit_resid_pooled df_fit
#> R1      0.736766065  322.5
#> R2      0.548246161  322.5
#> R3      0.008326603  322.5
#> R4     -0.508911000  322.5
```
