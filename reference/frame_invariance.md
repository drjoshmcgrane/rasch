# Test item invariance across frames

Calibrates each frame separately and compares the locations and
discriminations of items administered in more than one frame.

## Usage

``` r
frame_invariance(
  fit,
  alpha = 0.05,
  adjust = c("holm", "none"),
  se_method = c("conditional", "bootstrap"),
  boot_reps = 200,
  seed = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

- alpha:

  Significance level used for flags.

- adjust:

  Either `"holm"` or `"none"`. Both raw and adjusted probabilities are
  returned.

- se_method:

  `"conditional"` treats the estimated frame units as fixed;
  `"bootstrap"` refits the complete analysis to whole-person resamples
  within group, preserving each person's item-set response pattern. As
  in
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md),
  one response row is required per person.

- boot_reps:

  Number of bootstrap replicates. At least 30 are required. At least 90
  per cent, and no fewer than 30, must yield the complete set of
  comparisons.

- seed:

  Optional bootstrap seed. See
  [`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md)
  for generator support.

## Value

An object of class `"rasch_frame_invariance"`. The `locations` and
`discrimination` tables contain the pairwise item comparisons; `summary`
contains set-level RMSD and RMSE summaries. Under the conditional
method, discrimination `p`, `p_adj`, and `flagged` are `NA`. `excluded`
lists items dropped or rescored by a separate calibration, whose
observed category structures differed between calibrations, or whose
separate-frame estimate was weakly determined. The remaining components
record the multiplicity and uncertainty settings, including the
algorithm identifier, declared comparison-family size `family_n`, and
the requested, usable, non-converged and other-failure bootstrap counts.
`bootstrap_stratified` records whether persons were resampled within
group rather than globally.

## Details

Let \\\hat\delta\_{if}\\ be the location of item \\i\\ from a separate
calibration of frame \\f\\, and let \\\hat\rho_f\\ be that frame's unit
from the fitted EFRM. The common-scale location is
\\\hat\delta\_{if}^{\*}=\hat\delta\_{if}/\hat\rho_f\\. Because each
separate calibration has its own origin, pairwise differences are
centred over the common thresholds before testing.

The conditional method treats the fitted frame units as fixed. Let
\\w_i=m_i/\sum_jm_j\\, where \\m_i\\ is the number of thresholds for
item \\i\\. If \\C=I-\mathbf{1}\mathbf{w}^{\mathsf T}\\ centres the
common items on the threshold-weighted origin, the covariance of the
location differences is
\$\$C\\V_1/\hat\rho_1^2+V_2/\hat\rho_2^2\\C^{\mathsf T}.\$\$ This is
fast and conditions on the estimated units. The discrimination table
gives the difference between the two standardised infit statistics,
divided by \\\sqrt{2}\\, together with fitted slopes and their ratio.
These quantities are descriptive under the conditional method; it does
not report discrimination probabilities.

With `se_method = "bootstrap"`, whole persons are resampled within their
observed group, retaining each sampled person's item-set response
pattern, and the EFRM and separate frame calibrations are refitted. A
replicate is usable only when it retains the observed set of item
comparisons, so every centred difference has the same frame origin.
Location tests then use the empirical covariance of the centred
differences. The discrimination test uses the bootstrap standard error
of the log slope ratio. Both are standard deviations over the usable
replicates, so their statistics are referred to \\t(B-1)\\ rather than
the normal, where \\B\\ is the number of usable replicates. This
includes uncertainty in the fitted frame units but is more
computationally demanding.

Raw and Holm-adjusted probabilities are reported. With conditional
uncertainty, Holm adjustment covers the location comparisons. With
bootstrap uncertainty, it covers the combined family of location and
discrimination comparisons. An unavailable comparison remains in the
applicable family. A discrimination probability is unavailable when
either separate-frame slope is on its imposed estimation boundary; the
ratio remains descriptive and the comparison remains in the Holm family.
The summary gives the root mean squared location difference and root
mean squared standard error for each set and frame pair. Items from
different sets cannot be compared because the sets partition the items.
Location differences are relative to the mean difference of the common
items. Concentrated DIF can therefore produce non-zero centred contrasts
for items that were not themselves shifted. The table identifies the
pattern of relative departures; item content or external anchors are
needed to determine which items provide the defensible reference. A
compared set-by-frame cell must contain at least 50 distinct persons
contributing an informative item pair. A pair is informative unless both
responses are zero or both are at their item maxima, using the retained
items and recoded categories of that frame's separate calibration. Items
with weakly determined standard errors in either separate calibration
are listed in `excluded` rather than tested.

A flagged item may be resolved with
[`resolve_frames`](https://drjoshmcgrane.github.io/rasch/reference/resolve_frames.md)
when it remains useful within frames, or removed with
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
when it fits poorly more generally. Either change requires a refit. The
invariance tests require a converged frame calibration.

## References

Humphry, S. M. (2005). *Maintaining a Common Arbitrary Unit in Social
Measurement*. PhD thesis, Murdoch University.

## See also

[`resolve_frames`](https://drjoshmcgrane.github.io/rasch/reference/resolve_frames.md)
to give a flagged item a location per frame,
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
to remove it altogether, and
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
for the model whose assumption is tested.

## Examples

``` r
d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 1,
                   n_groups = 2, group_unit_ratio = 1.4, seed = 2)
tr <- attr(d, "truth")
fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                  id = "id", boot_reps = 0)
frame_invariance(fit)
#> Item invariance across frames (each frame calibrated separately)
#> 
#> Uncertainty: conditional on the fitted frame units 
#> 
#>   set frame_1 frame_2 n_items n_excluded  rmsd  rmse ratio n_location
#>  set1      g1      g2       8          0 0.199 0.200 0.994          0
#>  n_discrimination
#>                  
#> 
#> rmsd/rmse above 1 indicates item behaviour the frame units do not account for
#> 
#> No available item-location comparison differs across frames at alpha = 0.05 (Holm-adjusted).
#> 
#> The discrimination comparisons are descriptive:
#>   set frame_1 frame_2  item infit_1 infit_2 infit_z disc_1 disc_2 disc_ratio
#>  set1      g1      g2 S1I01   1.147   1.113   0.630  1.021  1.088      1.066
#>  set1      g1      g2 S1I02   1.035   1.023   0.172  1.267  1.241      0.980
#>  set1      g1      g2 S1I03   1.063   1.107  -0.356  1.212  1.131      0.933
#>  set1      g1      g2 S1I04   0.996   1.037  -0.464  1.448  1.266      0.875
#>  set1      g1      g2 S1I05   1.015   1.087  -0.763  1.374  1.173      0.854
#>  set1      g1      g2 S1I06   1.081   1.130  -0.357  1.175  1.110      0.945
#>  set1      g1      g2 S1I07   1.131   0.977   1.639  1.070  1.382      1.292
#>  set1      g1      g2 S1I08   1.069   1.118  -0.434  1.121  1.079      0.962
#>  disc_boundary
#>               
#>               
#>               
#>               
#>               
#>               
#>               
#>               
#> Use se_method = "bootstrap" for discrimination probabilities.
```
