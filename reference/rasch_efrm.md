# Fit the extended frame of reference model

Fits Humphry's extended frame of reference model, in which the unit can
differ across item-set by person-group frames. For item \\i\\ in set
\\s\\ and person \\n\\ in group \\g\\,
\$\$P(X\_{ni}=x)=\frac{\exp\\\rho\_{sg}\[x\theta_n-
\sum\_{k=1}^{x}\delta\_{ik}\]\\}
{\sum\_{y=0}^{m_i}\exp\\\rho\_{sg}\[y\theta_n-
\sum\_{k=1}^{y}\delta\_{ik}\]\\},\qquad \rho\_{sg}=\alpha_s\phi_g.\$\$

## Usage

``` r
rasch_efrm(
  data,
  item_sets,
  groups,
  id = NULL,
  factors = NULL,
  items = NULL,
  n_groups = NULL,
  na_codes = -1,
  maxit = 50,
  tol = 1e-07,
  min_link_persons = 30,
  se_method = c("hybrid", "bootstrap"),
  boot_reps = NULL,
  progress = NULL,
  cancel = NULL,
  workers = 4L,
  seed = NULL
)
```

## Arguments

- data:

  Persons-by-items data (matrix or data frame, like
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)),
  plus a person-group column.

- item_sets:

  A named list mapping set names to item-column names, or a named
  character vector mapping every analysed item exactly once to a set.
  The vector cannot name items outside the analysis. Items not mentioned
  form their own set `"(rest)"` when a list is given.

- groups:

  Name of the person-group column in `data`, or a vector with one entry
  per person. Several columns define crossed group cells. Their units
  are returned in `phi_table`; `phi_factorial` and `phi_factorial_tests`
  contain the GLS factorial decomposition and omnibus Wald tests, each
  referred to \\F(q, B-q)\\ when the cell-unit covariance came from the
  full bootstrap and to \\\chi^2(q)\\ when it is analytic. Raw
  probabilities are retained in `p`; decisions use `p_adj`,
  Holm-adjusted across the factorial terms. Structurally unidentified
  units are refused. Very imprecise but identified units are retained
  with a warning.

- id:

  Person identifier, either a column name or one value per row. EFRM
  data require one response row per person, so identifiers must be
  unique when supplied. Missing or blank identifiers are treated as
  different unknown persons.

- factors, items, n_groups, na_codes:

  As in
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- maxit, tol:

  Outer iteration cap and convergence tolerance of the bilinear pairwise
  stage.

- min_link_persons:

  Minimum number of common persons required for a set pair to contribute
  to the unit linking.

- se_method:

  `"hybrid"` (sandwich + linking bootstrap + delta propagation; default)
  or `"bootstrap"` (full person bootstrap of all stages).

- boot_reps:

  Bootstrap replicates; defaults to 300 for the linking bootstrap and
  200 for the full bootstrap. Use zero to omit set-link uncertainty; in
  a multi-set fit, common-unit item and threshold standard errors are
  then unavailable. Otherwise at least 30 replicates are required. A
  bootstrap covariance is reported only when more than half of the
  requested replicates are usable. Inference is returned only when at
  least 30 replicates succeed, a majority of those requested, and the
  requested count exceeds the number of independent directions in the
  largest covariance block used by the fit. The fit stops if the linking
  covariance cannot meet that rule; an unsuccessful full bootstrap falls
  back to hybrid standard errors with a warning and retains its
  replicate accounting.

- progress:

  Optional function called as `progress(stage, current, total)` during
  long uncertainty calculations. It is intended for interfaces and batch
  logging and does not alter estimation.

- cancel:

  Optional zero-argument function checked between bootstrap batches.
  Returning `TRUE` stops with a `rasch_cancelled` condition. A serial
  fit uses one replicate per batch.

- workers:

  Number of parallel bootstrap workers. The default is four, reduced
  when fewer physical cores are available or the R process has a lower
  system limit. Random samples are generated before distribution, so a
  fixed seed gives the same result for any worker count. Every worker
  holds its own copy of the bootstrap state.

- seed:

  Optional bootstrap seed. The caller's random-number state is restored
  when estimation finishes; see
  [`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md)
  for generator support.

## Value

An object of classes `"rasch_efrm"` and `"rasch"`. Model-specific
components include `frames`, `phi_table`, `alpha_table`, `set_table`,
common-unit item and threshold tables, group-specific `score_curves`
(expected weighted sufficient score and conditional standard error by
person location and exact observed-item pattern, one row block per
design and labelled as in
[`test_information`](https://drjoshmcgrane.github.io/rasch/reference/test_information.md)),
`efrm_vs_rasch`, and `linking`, the person support used for unit
inference in `unit_support`, and the active covariance blocks in
`unit_cov`. For a full-bootstrap fit, all blocks in `unit_cov` are
calculated from the same usable person resamples; otherwise they are the
analytic within-frame and, when requested, hybrid linking covariances
used by the reported tests. With several item sets and `boot_reps = 0`,
`cov_delta` and the corresponding common-unit standard errors are
unavailable because set-link uncertainty has not been estimated. The
requested, usable and failed uncertainty replicates used by the returned
uncertainty method are reported as `boot_reps_requested`,
`boot_reps_used` and `boot_reps_failed`; the hybrid set-link counts are
repeated inside `linking`. When a full bootstrap was requested, its
requested, attempted, usable and failed counts are retained separately
in the corresponding `full_boot_reps_*` components, including when the
fit falls back to hybrid standard errors. See the extended frame of
reference vignette for their interpretation. If the within-frame
calibration does not converge, its covariance blocks, standard errors
and all later inferential probabilities are withheld. Failure of only a
set link does not invalidate the already converged within-frame
calibration or group-unit estimates, but common-unit item, frame and
person uncertainty is withheld because it depends on that link,
including the standard errors in `score_curves`.

## Details

The partial credit model holds within each frame in its natural unit.
\\\phi_g\\ and \\\alpha_s\\ are unit *ratios* in the sense of Humphry
and Andrich (2008, eq. 15): each is the common reference unit over the
frame's own unit. The identification constraints set the geometric mean
of the group units and of the set units to one; no observed group or set
is the reference level. A value above one therefore denotes a finer
natural unit than the corresponding geometric-mean unit and steeper
curves on the common scale. Ratios between two observed levels are
obtained directly, for example as \\\alpha_s/\alpha_t\\. Person-group
ratios \\\phi_g\\ are identified from common item thresholds across
groups. Item sets partition the items, so set ratios \\\alpha_s\\ are
identified instead from persons observed in more than one set. The
set-linking graph and the group-by-set frame graph must each connect to
a common scale. Direct overlap between every pair of item sets is not
required: sets can be linked through intermediate sets. Pairs without
enough informative common persons contribute no edge; the remaining
graph must still connect all sets. A bootstrap replicate is unusable if
any supported link fails numerically or does not converge, even when
other links still connect the sets.

Set units use a semiparametric likelihood for persons observed in each
linked pair of sets. For sets \\a\\ and \\b\\, it maximises
\$\$\prod_n\int P(X\_{na}\mid u)P(X\_{nb}\mid ru+c)\\dF\_{g(n)}(u),\$\$
where the masses of each observed group's \\F_g\\, the scale ratio \\r\\
and the offset \\c\\ are estimated jointly on a fixed grid. This avoids
prescribing a normal or common person distribution across groups. A link
whose scale or offset reaches the numerical search boundary is refused
rather than reported as a finite estimate. So is a link whose grid
truncates the person distribution: persons with a finite location in at
least one set may place no more than two per cent of their posterior
mass on the ends of the grid. Persons at the same extreme tail in both
sets remain in the likelihood but do not count towards this check or the
link's inferential support; their likelihood rises towards one end of
any finite grid and carries no information about the link. Opposing
extremes retain a finite compromise location and do contribute. The
conditional thresholds and group units are held fixed in this step; only
\\r\\, \\c\\, and the nuisance masses are estimated. The linked
parameters are then
\$\$\delta\_{ik}=\widetilde\delta\_{ik}/\alpha_s+\mu_s, \qquad
\rho\_{sg}=\alpha_s\phi_g.\$\$ Score moments supply starting values and
screen weak links. Response patterns must span a score range of at least
four within a set. Overlapping item sets are not permitted. The public
convergence flag covers the conditional calibration, the set-link
transformation and its nonparametric nuisance masses; `stage1_converged`
records the conditional stage separately. The conditional stage requires
a small score and negative curvature of the exact likelihood Hessian in
all free directions. A stationary point with flat or positive curvature
is refused, including in bootstrap refits. This checks an identified
local maximum, not a global maximum.

The empirical Godambe covariance from the conditional stage requires at
least ten informative persons, at least eight effective persons, more
effective persons than fitted stage-one directions, and full rank in the
projected score covariance. If these conditions fail, point estimates
can be returned with `boot_reps = 0`, but unit uncertainty is withheld.
A hybrid or full-bootstrap fit is refused because its linking and unit
covariance would otherwise inherit an unsupported stage-one covariance.

The hybrid covariance combines the pairwise Godambe covariance with a
person bootstrap for set linking. Each replicate jointly redraws the
within-frame thresholds and group units, then rebuilds the link. The
joint draws retain covariance among common-scale thresholds, set units
and group units. With `se_method = "bootstrap"`, the complete model is
refitted to each person resample. Refits that do not converge or have
unidentified group units are discarded and counted as failed replicates.

The `efrm_vs_rasch` component records the within-frame composite
log-likelihood comparison between group-dependent and equal group units.
This difference is descriptive and contains no information about set
units, which are identified at the linking stage. The accompanying Wald
omnibus tests provide inference for the group- and set-unit families.
Their probabilities are Holm-adjusted as one omnibus family; the
individual unit contrasts form a second Holm-adjusted follow-up family.
Each family counts the distinct hypotheses it declares, available or
not: the units are centred, so with two groups (or two sets) the two
reported rows are one hypothesis stated twice and count once, as the
omnibus rank already does. The second row of such a pair keeps its
estimate and unadjusted probability, but its adjusted probability and
flag are withheld, so one difference is not reported as two deviating
units. Beyond two groups (or two sets) no two reported rows are the same
hypothesis, so each stays a member: that family is then one larger than
its free dimension, which leaves the adjustment conservative rather than
liberal. A bootstrap standard error is a standard deviation over the
retained replicates, so its contrast is referred to \\t(B-1)\\; an
analytic standard error keeps the normal reference. The reference is
reported as `df`. An omnibus Wald test on an estimated (bootstrap)
covariance is referred to \\F(q, B-q)\\ on the same grounds, reported as
`df`, `df2` and `f`, so that a one-dimensional omnibus and its unit
contrast report the same probability; an analytic covariance keeps the
chi-square reference. This holds for every omnibus Wald test the fit
reports, including the crossed group-unit decomposition in
`phi_factorial_tests`, so one printed fit never refers two tests on the
same draws to two different references. Unit estimates are retained for
sparse designs, but probabilities require at least 50 persons or
effective persons in every group. Set-unit inference requires at least
50 informative common persons on the strongest bottleneck path from
every set to the first set, which is used only as the support graph's
bookkeeping root. Thus a weak upstream link limits a terminal set, while
a weak redundant edge does not suppress a stronger route. Group-unit and
dependent set-unit probabilities are withheld when any group unit has a
reported standard error above 5 log units. The estimates and covariance
remain descriptive. This check uses the returned uncertainty method,
including the full bootstrap when available.

The model assumes that an item retains its location and discrimination
across the frames in which it appears, apart from the frame unit.
[`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
examines this assumption by separate frame calibrations. Misfit
concentrated within one item set can also distort its estimated unit;
inspect item fit and targeting before interpreting unit differences.
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
and
[`resolve_frames`](https://drjoshmcgrane.github.io/rasch/reference/resolve_frames.md)
provide refitted sensitivity analyses.

The dichotomous model follows Humphry (2005) and Humphry and Andrich
(2008). The polytomous, multigroup and crossed-frame forms are
extensions implemented in this package. The discrete nonparametric
margin follows the Rasch estimation approach of Follmann (1988); its use
for linked item-set units is an extension implemented here.

## References

Andrich, D. (1982). An extension of the Rasch model for ratings
providing both location and dispersion parameters. Psychometrika, 47(1),
105–113.

Andrich, D. and Luo, G. (2003). Conditional pairwise estimation in the
Rasch model for ordered response categories using principal components.
Journal of Applied Measurement, 4(3), 205–221.

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

Follmann, D. (1988). Consistent estimation in the Rasch model based on
nonparametric margins. Psychometrika, 53, 553–562.
[doi:10.1007/BF02294407](https://doi.org/10.1007/BF02294407)

Humphry, S. M. (2005). Maintaining a Common Arbitrary Unit in Social
Measurement. PhD thesis, Murdoch University.

Humphry, S. M. (2010). Modeling the effects of person group factors on
discrimination. Educational and Psychological Measurement, 70(2),
215–231.

Humphry, S. M. (2012). Item set discrimination and the unit in the Rasch
model. Journal of Applied Measurement, 13(2), 165–180.

Montuoro, P. and Humphry, S. M. (2024). Modeling the effect of reading
item clarity on item discrimination. Journal of Applied Measurement,
24(3/4), 121–132.

Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
Rasch model. Journal of Applied Measurement, 9(3), 249–264.

## See also

[`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md),
which tests the item invariance this model assumes rather than imposing
it, and
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md),
which removes an item the test flags and refits. Also
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
[`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
[`test_information`](https://drjoshmcgrane.github.io/rasch/reference/test_information.md),
and
[`simulate_efrm`](https://drjoshmcgrane.github.io/rasch/reference/simulate_efrm.md).

## Examples

``` r
# \donttest{
set.seed(1); Np <- 400
simP <- function(th, tau, r) { x <- 0:length(tau)
  p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p) }
grp <- rep(c("A", "B"), each = Np / 2)
phi <- c(A = 0.8, B = 1.25)
d <- seq(-1.5, 1.5, length.out = 10)
theta <- rnorm(Np)
X <- sapply(seq_along(d), function(i) sapply(seq_len(Np), function(n)
  sample(0:1, 1, prob = simP(theta[n], d[i], phi[grp[n]]))))
colnames(X) <- sprintf("I%02d", seq_along(d))
fit <- rasch_efrm(data.frame(X, grp = grp), item_sets = list(core = colnames(X)),
                  groups = "grp")
fit$phi_table
#>  group   phi se_log_phi
#>      A 0.842      0.042
#>      B 1.187      0.042
# }
```
