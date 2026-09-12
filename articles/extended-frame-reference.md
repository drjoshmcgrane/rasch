# The Extended Frame of Reference model

``` r

library(rasch)
```

## Frames and units

Rasch’s criterion of invariant comparison is stated within a frame of
reference (Rasch 1961). The extended frame of reference model allows the
unit to differ over linked combinations of item sets and person groups
(Humphry 2005; Humphry and Andrich 2008). For item \\i\\ in set \\s\\
and person \\n\\ in group \\g\\,

\\ P(X\_{ni}=x)=
\frac{\exp\\\rho\_{sg}\[x\theta_n-\sum\_{k=1}^{x}\delta\_{ik}\]\\}
{\sum\_{y=0}^{m_i}\exp\\\rho\_{sg}\[y\theta_n-
\sum\_{k=1}^{y}\delta\_{ik}\]\\},\qquad \rho\_{sg}=\alpha_s\phi_g . \\

Following Humphry and Andrich (2008, eq. 15), the multiplicative
parameter of a frame is the ratio of the common reference unit to the
frame’s own unit: \\\alpha_s\\ carries that ratio for item set \\s\\,
and \\\phi_g\\ for person group \\g\\. Their geometric means are fixed
at one; no observed set or group is the reference. A set with \\\alpha_s
= 1.3\\ has a smaller natural unit than the geometric-mean set unit, so
its response curves are steeper on the common scale. Two observed sets
can be compared directly by the ratio \\\alpha_s/\alpha_t\\. The partial
credit model holds within each frame in its natural unit; dichotomous
items are the one-threshold case.

Use this model when variation in the unit is part of the measurement
account. Poor fit to an equal-unit model is not by itself evidence for a
frame-dependent unit.

## Fit the model

``` r

d <- simulate_efrm(
  n_per_group = 150,
  items_per_set = 8,
  set_unit_ratio = 1.30,
  group_unit_ratio = 1.10,
  seed = 25
)
truth <- attr(d, "truth")
d$site <- rep(c("A", "B"), length.out = nrow(d))
d
#> Simulated efrm data: 300 persons, 2 sets x 2 groups, 16 items (2 categories)
#> Generating values and departures:
#>   - set-unit ratio 1.30 across 2 sets
#>   - group-unit ratio 1.10 across 2 groups
```

The simulated data carry the units the fit below has to recover: a
set-unit ratio of 1.30 and a group-unit ratio of 1.10.

``` r

fit <- rasch_efrm(
  d,
  item_sets = truth$item_sets,
  groups = "group",
  id = "id",
  factors = "site",
  boot_reps = 50,
  workers = 1,
  seed = 25
)
fit
#> rasch extended frame of reference analysis: 16 items in 2 set(s) x 2 group(s) = 4 frames, 300 persons
#> Within-frame pairwise conditional ML: converged in 10 iterations
#> PSI 0.781, separation quality: reasonable
#> 
#> Person group units (phi):
#>  group   phi se_log_phi
#>     g1 0.897      0.040
#>     g2 1.115      0.040
#> 
#> Item set units (alpha) and locations:
#>   set alpha se_log_alpha     mu n_items
#>  set1 0.840        0.043  0.081       8
#>  set2 1.190        0.043 -0.081       8
#> 
#> Equal-unit comparison: 2(ll_EFRM - ll_equal) = 20.620 with 1 extra unit parameter(s)
#> (composite likelihood: descriptive; informative for group units (phi))
#> Omnibus Wald tests of equal units (Holm-adjusted family):
#>               term df df2   wald      f       p   p_adj significant
#>  group units (phi)  1 Inf  7.442  7.442   0.006   0.006           *
#>  set units (alpha)  1  49 16.034 16.034 < 0.001 < 0.001           *
#> Holm-adjusted exploratory unit contrasts (H0: unit = 1):
#>        parameter estimate    se      z  df       p   p_adj significant
#>      log phi[g1]   -0.109 0.040 -2.728 Inf   0.006   0.006           *
#>      log phi[g2]    0.109 0.040  2.728 Inf   0.006                    
#>  log alpha[set1]   -0.174 0.043 -4.004  49 < 0.001 < 0.001           *
#>  log alpha[set2]    0.174 0.043  4.004  49 < 0.001                    
#> (the units are centred: the second row of a two-row family restates the first, so its adjustment is withheld)
#> 
#> Notes: person measures use the weighted score; per-group score curves replace the raw-score table (see score_curves); a universal raw-score conversion is not defined across the expanded frame response cells; use score_curves and design-specific information; two centred unit rows are one hypothesis: the adjusted probability and flag are reported on the first row of the pair and withheld on the second, which restates it
```

The principal model-specific tables are:

``` r

fit$phi_table       # person-group units
#>  group   phi se_log_phi
#>     g1 0.897      0.040
#>     g2 1.115      0.040
fit$alpha_table     # item-set unit ratios (reference unit / set unit)
#>   set alpha se_log_alpha
#>  set1 0.840        0.043
#>  set2 1.190        0.043
fit$set_table       # linked set locations
#>   set     mu alpha n_items
#>  set1  0.081 0.840       8
#>  set2 -0.081 1.190       8
fit$frames          # complete frame units
#>   set group n_persons n_items alpha   phi   rho se_log_rho origin infit_ms
#>  set1    g1       150       8 0.840 0.897 0.754      0.062  0.081    1.031
#>  set2    g1       150       8 1.190 0.897 1.067      0.056 -0.081    1.048
#>  set1    g2       150       8 0.840 1.115 0.937      0.056  0.081    1.063
#>  set2    g2       150       8 1.190 1.115 1.327      0.062 -0.081    1.012
#>  outfit_ms fit_resid n_responses
#>      0.989     0.185        1168
#>      0.988    -0.316        1168
#>      1.008     0.422        1184
#>      0.865    -1.361        1184
fit$linking         # set-linking design
#> $phi_edges
#> $phi_edges[[1]]
#> [1] 2 1
#> 
#> $phi_edges[[2]]
#> [1] 2 1
#> 
#> 
#> $alpha_edges
#>  set_a set_b   n log_slope converged edge_mass    loglik
#>   set1  set2 294     0.348         1     0.000 -3522.116
#> 
#> $boot_reps_requested
#> [1] 50
#> 
#> $boot_reps_used
#> [1] 50
#> 
#> $boot_reps_failed
#> [1] 0
#> 
#> $alpha_method
#> [1] "finite-grid semiparametric maximum likelihood"
#> 
#> $alpha_grid
#>  lower  upper points 
#>     -8      8     61
```

The frame-defining group is part of the model and is not retested as
DIF. Other person factors can be examined with
[`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md).
Its bootstrap conditions on each person’s subtotal within each observed
item set, then repeats the Extended Frames fit and the complete DIF
analysis. Its familywise probabilities refer to the fitted global
invariant null.

``` r

dif <- dif_anova(fit)
dif_bootstrap(fit, dif, B = 999, seed = 2026)$summary
```

Person-group units are estimated from common threshold patterns across
groups. Set units are estimated from persons observed in more than one
set. For a linked pair, the natural coordinates satisfy \\u_b=r u_a+c\\.
The package estimates \\r\\, \\c\\, and masses for the nuisance person
distribution on a finite grid by semiparametric likelihood:

\\ \prod_n \int P(X\_{na}\mid u)P(X\_{nb}\mid ru+c)\\dF\_{g(n)}(u). \\

The masses are estimated jointly with the link, separately for each
observed person group, so neither a normal shape nor a common
distribution across groups is prescribed. The link holds the conditional
item thresholds and person-group units fixed while estimating only
\\r\\, \\c\\, and the nuisance masses. The common-scale parameters are
then

\\ \delta\_{ik}=\widetilde\delta\_{ik}/\alpha_s+\mu_s, \qquad
\rho\_{sg}=\alpha_s\phi_g. \\

The item thresholds and person-group units remain pairwise conditional
estimates; the mixing distribution is used only for the relative set
unit, which is not identified by that conditional stage. The discrete
nonparametric margin follows Follmann (1988); its use for linked
item-set units is an extension implemented here.

The group-by-set frame graph and the set-linking graph must each connect
to a common scale. A linking response pattern must span a score range of
at least four within a set. Overlapping item sets are not permitted.

Several group factors may define crossed frames:

``` r

fit_crossed <- rasch_efrm(
  data,
  item_sets = item_sets,
  groups = c("language", "cohort"),
  id = "id"
)
fit_crossed$phi_factorial
fit_crossed$phi_factorial_tests
```

The factorial tables give the GLS decomposition of the estimated log
group units. Use `p_adj`, which applies Holm’s correction across the
factorial terms, for decisions. The tables do not replace the
frame-level estimates.

## Uncertainty and equal-unit comparisons

The default hybrid covariance combines the pairwise Godambe covariance
with a person bootstrap for set linking. Within-frame thresholds and
group units are jointly redrawn, and the link is rebuilt from each draw.
The resulting covariance is retained for the common-scale thresholds and
complete frame units. Set `se_method = "bootstrap"` to refit the
complete model to each person resample.

The linking calculations use a compiled numerical kernel. Bootstrap
replicates run on four workers by default, or fewer where the system
limit is lower. A fixed `seed` gives the same samples and results for
every worker count without altering the caller’s random-number stream.
In the Shiny application the fit runs in a separate process. The analyst
can select up to four workers, follow completed batches, or cancel the
process without replacing the current analysis. R scripts can supply
`progress` and `cancel` callbacks for the same controls.

The fit records the requested, usable and failed bootstrap counts. The
requested number must exceed the free directions in the largest
covariance block used by the fit. The stored bootstrap row contains
several constrained blocks and is not treated as one covariance. Once
estimation begins, the linking covariance needs at least 30 usable
replicates and a majority of those requested; otherwise the fit stops. A
numerical failure or non-convergence on any supported link invalidates
the whole replicate. Pairs without enough informative common persons can
be omitted only if the remaining links connect all sets. Full-bootstrap
refits with unidentified group units are also discarded. The conditional
calibration checks the exact likelihood curvature as well as the score.
Stationary saddles and flat solutions are refused; the same check
applies to bootstrap refits. Saved frame fits without a current
likelihood-check record must be refitted before reopening in the app. So
must an extended frame fit whose stored score curves predate the shared
design enumeration, since their designs and labels are no longer the
ones
[`test_information()`](https://drjoshmcgrane.github.io/rasch/reference/test_information.md)
and the curve plots would give. Neither refusal alters the source file.
The alpha–phi cross-covariance has the same usable-draw rule and is
withheld if it cannot be estimated. A full bootstrap that falls below
its minimum returns hybrid standard errors instead, retaining its
attempted, usable and failed counts in `full_boot_reps_*`.

``` r

fit$efrm_vs_rasch$unit_omnibus
#>               term df df2   wald      f       p   p_adj significant
#>  group units (phi)  1 Inf  7.442  7.442   0.006   0.006           *
#>  set units (alpha)  1  49 16.034 16.034 < 0.001 < 0.001           *
fit$efrm_vs_rasch$unit_tests
#>        parameter estimate    se      z  df       p   p_adj significant
#>      log phi[g1]   -0.109 0.040 -2.728 Inf   0.006   0.006           *
#>      log phi[g2]    0.109 0.040  2.728 Inf   0.006                    
#>  log alpha[set1]   -0.174 0.043 -4.004  49 < 0.001 < 0.001           *
#>  log alpha[set2]    0.174 0.043  4.004  49 < 0.001
```

The omnibus Wald tests assess the set- and group-unit families. Their
probabilities are adjusted together by Holm’s method. Individual unit
contrasts form a separate Holm-adjusted follow-up family. The raw
composite-likelihood difference is descriptive and compares the
group-unit stage only. Set units are identified in the linking stage and
are assessed by their Wald test.
[`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
withholds generic likelihood differences involving EFRM: its calibration
likelihood uses within-set pairs, whereas ordinary PCM also uses
cross-set pairs. The matched group-unit comparison is retained in
`fit$efrm_vs_rasch`. Probabilities require at least 50 persons or
effective persons in every group and a path between each pair of sets
whose links each have at least 50 informative common persons. A person
at the same extreme tail in both sets remains in the likelihood but does
not count as link support. Sparse designs retain the unit estimates
without an inferential probability.

Two-set simulations under normal, bimodal and deliberately different
group distributions gave set-unit bias within 0.004 log-units,
empirical-to-reported SE ratios from 0.97 to 1.05, and null rejection
from 4.0 to 5.0 per cent for the hybrid method. The complete person
bootstrap was mildly conservative in the corresponding null design. Sets
with only three dichotomous items were not identified; four-item sets
were recovered without bias in adequately supported samples.

A separate null study checked the complete decision families. With 200
persons per group and six items per set, 489 of 500 fits were analysed.
The Holm familywise rejection rates were 3.5 per cent for the omnibus
tests and 1.6 per cent for the individual unit contrasts. A boundary
design with 100 persons per group and four items per set refused 493 of
1,000 fits and did not converge in another 24. Among the 483 analysed
fits, the corresponding rates were 5.2 and 2.5 per cent. A design should
not be used merely because it is formally identifiable.

## Item invariance across frames

The fitted EFRM gives an item one location across its frames, scaled by
the frame unit.
[`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
examines that restriction by calibrating each frame separately. Each
compared set and frame needs at least 50 distinct persons contributing
an informative item pair: not both zero or both at their item maxima
after any item removal or category recoding in that frame’s separate
calibration. Items with weak separate-frame standard errors are excluded
from the comparison.

For frame \\f\\, the common-scale item location is
\\\hat\delta\_{if}^{\*}=\hat\delta\_{if}/\hat\rho_f\\. The separate
calibrations have independent origins, so differences are centred over
the common thresholds before testing. Let \\w_i=m_i/\sum_jm_j\\, where
\\m_i\\ is the number of thresholds for item \\i\\. The default
covariance conditions on the fitted frame units:

\\ C\\V_1/\hat\rho_1^2+V_2/\hat\rho_2^2\\C^{\mathsf T}, \qquad
C=I-\mathbf 1\mathbf w^{\mathsf T} . \\

``` r

inv <- frame_invariance(fit, se_method = "conditional")
inv$summary
inv$locations
inv$discrimination
```

A whole-person bootstrap within person group retains each sampled row’s
item-set response pattern and includes uncertainty in the fitted units.
A replicate is used only if it retains the same item-comparison family,
and therefore the same centred frame origin, as the observed analysis:

``` r

inv_boot <- frame_invariance(
  fit,
  se_method = "bootstrap",
  boot_reps = 300,
  seed = 1
)
```

The conditional method returns raw and Holm-adjusted probabilities for
the location comparisons. Its standardised-infit and slope-ratio columns
are descriptive. The bootstrap returns discrimination probabilities as
well and applies Holm adjustment to the combined location and
discrimination family.

In 300 null simulations with 500 persons per frame, the bootstrap SE
ratios were 1.00 for locations and 1.03 for log discrimination ratios;
95 per cent coverage was 0.948 and 0.954, and combined Holm familywise
error was 3.0 per cent. A separate 2,000-replicate check found that
probabilities based only on the conditional standardised-infit
comparison gave 7.1 per cent combined familywise error. The conditional
method therefore leaves those probabilities undefined. At the same
sample size, bootstrap power was 96.3 per cent for two items shifted by
one logit, but 9.6 per cent for two items made 1.5 times as steep.
Discrimination departures require appreciably more information.

The location contrasts are relative to the mean difference of the common
items. If a small number of items is shifted, centring also gives the
remaining items non-zero relative contrasts. The table identifies the
pattern of departures; item content or external anchors are needed to
decide which items provide the defensible reference.

## Refit after a departure

[`resolve_frames()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_frames.md)
gives an item a separate location in each frame. The item continues to
contribute to person measurement but no longer links those frames.
[`drop_items()`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
removes it. Both functions refit the model and update the frame, item
and person estimates. Each set must retain common items linking its
groups’ origins; a unit link through another set is not sufficient.

``` r

resolved <- resolve_frames(fit, "S1I02")
removed <- drop_items(fit, "S1I02")
```

The remaining common items and linked sets must still identify the frame
units. The functions preserve the original crossed factors, uncertainty
method and fitting controls, and return a result only when the revised
model converges.

## References

Follmann, D. (1988). Consistent estimation in the Rasch model based on
nonparametric margins. *Psychometrika*, 53, 553–562.
<https://doi.org/10.1007/BF02294407>

Humphry, S. M. (2005). *Maintaining a Common Arbitrary Unit in Social
Measurement*. PhD thesis, Murdoch University.

Humphry, S. M., and Andrich, D. (2008). Understanding the unit in the
Rasch model. *Journal of Applied Measurement*, 9(3), 249–264.

Rasch, G. (1961). On general laws and the meaning of measurement in
psychology. In *Proceedings of the Fourth Berkeley Symposium on
Mathematical Statistics and Probability* (Vol. 4, pp. 321–333).
University of California Press.
