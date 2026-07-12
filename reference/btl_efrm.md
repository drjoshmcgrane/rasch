# Fit the extended frame of reference model for paired comparisons

Estimates object locations from paired comparisons when the unit of the
latent scale differs across frames – judge-panel by object-set cells –
an extension of Humphry's (2005) extended frame of reference model to
the Bradley-Terry-Luce family. Objects are partitioned into sets and
judges into panels. Each object has a common-scale value
`v_k = alpha_s beta_k + kappa_s`, with `beta_k` the within-set
calibration location, `alpha_s > 0` the set unit and `kappa_s` the set
origin; a comparison judged in panel `g` carries the panel unit `phi_g`.
Within a set the comparison logit is `phi_g (beta_a - beta_b)` (the set
origin cancels and the set unit is confounded with the spread of `beta`,
so neither is identified within a set); across sets it is
`phi_g (v_a - v_b)`, which places the sets on one common scale.

## Usage

``` r
btl_efrm(
  data,
  object_a,
  object_b,
  winner,
  judge,
  panels,
  object_sets,
  response = NULL,
  ties = c("drop", "error"),
  min_link = 20,
  se_method = c("bootstrap", "conditional"),
  boot_reps = 60,
  maxit = 60,
  tol = 1e-08
)
```

## Arguments

- data:

  A data frame with one comparison per row.

- object_a, object_b:

  Names of the columns holding the two compared objects.

- winner:

  Name of the column holding the winner; its value must equal the row's
  `object_a` or `object_b` entry, `"tie"` or `"draw"` marks a tie, and
  anything else is treated as missing.

- judge:

  Name of the judge column (clusters the stage-one standard errors and
  defines the panels when `panels` is a judge attribute).

- panels:

  Either the name of a judge-attribute column in `data` or a named
  vector mapping judge to panel.

- object_sets:

  A named list mapping set names to character vectors of object names;
  every compared object must belong to exactly one set.

- response:

  Not supported: this first implementation fits dichotomous winner data
  only. Supplying it raises an informative error.

- ties:

  `"drop"` (default, removed with a note) or `"error"`.

- min_link:

  Minimum number of cross-set comparisons a set pair must supply to be
  used for linking; sets not reachable from the reference set through
  sufficient cross-set pairs raise an error.

- se_method:

  `"bootstrap"` (the default): a parametric bootstrap of the ENTIRE
  two-stage pipeline – winners resampled from the fitted probabilities,
  both stages refitted `boot_reps` times – so the reported standard
  errors carry every source of sampling variability, including the
  stage-one uncertainty that flows into the linking. `"conditional"`:
  the fast analytic errors, exact for `beta` and `phi` but conditional
  on stage one for `alpha` and `kappa`, which they therefore UNDERSTATE
  (verified by simulation); for quick inspection only.

- boot_reps:

  Bootstrap replicates for `se_method = "bootstrap"`.

- maxit, tol:

  Newton iteration cap and convergence tolerance.

## Value

An object of class `"rasch_btl_efrm"`: `objects` (object, set,
`beta_set` and its standard error, common-scale `v` and its standard
error), `phi_table` (panel units with Wald tests against `log phi = 0`),
`alpha_table` and `kappa_table` (set units and origins with Wald tests;
the reference set carries `alpha = 1`, `kappa = 0` with no standard
error), `frames` (panel by set: unit `rho = phi alpha`, comparison
count, pooled fit residual), `equal_unit` (the descriptive single-unit
comparison), `n_cross` (cross-set comparison counts per set pair),
`notes` and `converged`.

## Details

One modelling convention deserves plain statement. Rewritten with a
single latent value per object, the model says within-set-`s` contests
are judged at discrimination `phi_g / alpha_s` while every cross-set
contest is judged at exactly `phi_g`: the common scale is *defined* as
the scale of between-frame judgement. That is an assumption about how
discriminal dispersion behaves when unlike objects meet (a Thurstonian
question with no assumption-free answer), not a consequence of the
theory; the per-frame fit residuals in `frames` are where a violation
would show.

Estimation follows
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
in two conditional stages. In stage one the within-set comparisons of
each set, pooled over panels, fit the bilinear model
`logit = rho_{gs} (b_a - b_b)` with `b` sum-zero and the most-used panel
fixed at `rho = 1`; the ratios `rho_{gs} = phi_g / phi_{ref(s)}`
estimate the panel units up to each set's reference panel and are
reconciled across sets by a precision-weighted least squares over the
panel-by-set linking graph, then normalised to geometric mean one. In
stage two the cross-set comparisons are a low-dimensional maximum
likelihood in `(log alpha, kappa)` for the non-reference sets, with
`alpha_1 = 1` and `kappa_1 = 0` fixing the reference set (the first,
alphabetically).

The set units are identified here WITHIN the conditional framework: the
cross-set linking uses only comparison outcomes and makes no
distributional assumption about the objects. This is the substantive
difference from the persons-by-items EFRM, whose item-set units can only
be identified from the person side – that is, from the distribution of
the persons over the linked sets. The paired-comparison design supplies
its own conditional link, so no such distributional step is needed. See
Humphry (2005) and Humphry and Andrich (2008) for the theory of the
unit, and Thurstone (1927) and David (1988) for the
varying-discriminal-dispersion lineage from which the frame-dependent
unit descends. The paired-comparison form is this package's extension of
Humphry's model.

The ESTIMATES are always the staged conditional estimator: the
within-frame calibrations are invariant to the linking data (the
frame-of-reference property), a deliberate trade of some efficiency for
invariance, exactly as anchored equating trades. Inference defaults to a
parametric bootstrap of the whole pipeline (`se_method = "bootstrap"`):
winners are resampled from the fitted probabilities and both stages
refitted, so the standard errors of `phi`, `alpha`, `kappa`, `beta` and
the common-scale `v` carry every stage's sampling variability – verified
by simulation to restore nominal coverage where the conditional errors
cover at a third of their nominal rate on linked designs. The
`"conditional"` option keeps the fast analytic errors (judge-clustered
sandwich for stage one, inverse observed information for stage two,
conditional on stage one) for quick inspection; its `alpha` and `kappa`
errors understate, and the fit says so.

Two honesty notes on the bootstrap. It is model-based: replicates are
drawn as independent Bernoulli outcomes at the fitted probabilities,
which is self-consistent (the model has no judge parameter) but does not
carry extra-model dependence within judges; the conditional stage-one
errors are judge-clustered and guard against exactly that, so when the
two disagree materially the larger is the cautious choice. And a
parameter that reaches its boundary in some replicates (a set unit
driven to zero when a resampled within-set order flips against the
cross-set evidence, the signature of a two-object set with a near-even
internal pair) has no normal sampling distribution: its standard error
is reported as `NA` and a note names the parameter and the boundary
count rather than manufacturing a number. Relatedly, a set whose
within-set contests are all near-even (or all one-sided) carries no
stable information about the panel-unit ratios; such sets are screened
out of the `phi` reconciliation, refit with the panel units held at the
reconciled `phi` (which the frame model says apply to them regardless),
and named in a note.

A single set (`S = 1`) reduces the model to panel units alone; stage two
is skipped and the print states the panel-units model. When additionally
`G = 1` the fit reduces exactly to
[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) on the
same data. The equal-unit (single-unit) comparison refits plain
[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) on all
comparisons pooled and reports the descriptive composite log-likelihood
difference against the frames model; because that comparison is a
composite likelihood, the inference on the units is carried by the Wald
tests on `log phi_g` and `log alpha_s` in `phi_table` and `alpha_table`.

## References

Bradley, R. A. and Terry, M. E. (1952). Rank analysis of incomplete
block designs: I. The method of paired comparisons. Biometrika, 39,
324-345.

David, H. A. (1988). The Method of Paired Comparisons (2nd ed.).
Griffin.

Humphry, S. M. (2005). Maintaining a common arbitrary unit in social
measurement. PhD thesis, Murdoch University.

Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
Rasch model. Journal of Applied Measurement, 9(3), 249-264.

Luce, R. D. (1959). Individual Choice Behavior. Wiley.

Thurstone, L. L. (1927). A law of comparative judgment. Psychological
Review, 34, 273-286.

## Examples

``` r
# \donttest{
d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                       set_units = c(1, 1.4), set_origins = c(0, 0.8),
                       seed = 1)
fit <- btl_efrm(d, "object_a", "object_b", winner = "winner",
                judge = "judge", panels = "panel",
                object_sets = attr(d, "truth")$object_sets)
fit$alpha_table
#>    set    alpha se_log_alpha        z          p
#> 1 set1 1.000000           NA       NA         NA
#> 2 set2 1.606727    0.1908155 2.485119 0.01295081
# }
```
