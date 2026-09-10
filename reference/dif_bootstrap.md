# Bootstrap sensitivity analysis for DIF

Repeats a residual DIF analysis under the fitted invariant model. Rasch,
explanatory Rasch and Multiple Ratings fits condition on each person's
raw score. Extended Frames fits condition on the person's subtotal
within each observed item set. Comparative Judgement fits draw outcomes
from the fitted comparison model. Every replicate refits the model and
repeats the complete DIF design. The ordinary adjusted analysis remains
the primary DIF result.

## Usage

``` r
dif_bootstrap(fit, dif = NULL, B = 999, workers = 4L, seed = NULL)
```

## Arguments

- fit:

  A fitted Rasch, Multiple Ratings, Extended Frames, explanatory Rasch
  or ordinary Comparative Judgement model. Fully anchored scoring fits
  are not supported by this refitting procedure.

- dif:

  A current
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  or
  [`btl_dif`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  result from `fit`. For person-by-item models it may be omitted and the
  default DIF analysis is then computed. Comparative Judgement requires
  an explicit result because its judge factors have no default.

- B:

  Number of bootstrap replicates.

- workers:

  Number of parallel workers. The default is four, subject to the limits
  reported by the operating system and job scheduler.

- seed:

  Optional non-negative whole-number seed.

## Value

An object of class `"rasch_dif_bootstrap"`. Its `summary` and `terms`
tables add marginal and familywise bootstrap probabilities to the
observed DIF analysis. `replicates` contains the complete replicated F
statistics and their F-reference probabilities; the requested, usable,
non-converged and failed counts are recorded separately.

## Details

Responses are drawn from the conditional distribution \$\$P(\mathbf
X_p=\mathbf x\mid R_p=r_p,\boldsymbol\delta),\$\$ where \\R_p\\ is
person \\p\\'s observed raw score over that person's observed items. The
person parameter cancels by sufficiency. Thus the generator retains the
observed score, booklet or missingness pattern and repeated-person row
structure without drawing an ability distribution. For Extended Frames,
the same calculation is applied within each item set: the frame unit is
common to the set, so conditioning on the set subtotal cancels the
person parameter. This conditions on more information than a single
total when a person sees several sets, but remains an exact null
conditional distribution. For paired comparisons, responses are drawn
from the fitted category probabilities (and generated sequentially when
history effects were fitted), retaining judges and the comparison
design. Half-weighted ties and undefeated or winless objects are refused
because they do not supply fitted outcome probabilities for the required
null.

A replicate contributes only when its calibration converges and every
member of the declared item- or object-by-DIF-term family is estimable.
Marginal probabilities compare each observed term's F-reference
probability with the corresponding replicated probabilities. This puts
refits whose numerator or denominator degrees of freedom differ on the
same tail scale. Familywise probabilities use the single-step
distribution of the smallest term-wise F-reference probability in each
replicate. The same transformation is applied to the observed and
replicated statistics. These adjustments describe the fitted global
invariant null. They do not guarantee familywise error control among
otherwise invariant items or objects when another member has DIF,
because that departure can affect the fitted calibration and matching
scores. For `B >= 30`, at least 90 per cent of the requested replicates
and no fewer than 30 must be usable; a smaller exploratory run must
retain a majority.

BTL-EFRM and explanatory Comparative Judgement fits are refused because
[`btl_dif`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
does not define judge-group DIF after those structural restrictions. A
bootstrap cannot supply an estimand that the fitted model does not
define.

## References

Andrich, D. and Marais, I. (2019). *A Course in Rasch Measurement
Theory*. Springer.

Westfall, P. H. and Young, S. S. (1993). *Resampling-Based Multiple
Testing: Examples and Methods for p-Value Adjustment*. Wiley.

## See also

[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md),
[`btl_dif`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md),
[`fit_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md),
[`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md)

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(300 * 6, 1, .5), 300, 6)
colnames(X) <- paste0("I", 1:6)
g <- factor(rep(c("A", "B"), each = 150))
fit <- rasch(X, factors = data.frame(group = g))
d <- dif_anova(fit)
# Small only to keep the example quick; use substantially more replicates
# for inferential work.
db <- suppressWarnings(dif_bootstrap(fit, d, B = 3, seed = 1))
head(db$summary)
#>  item  term F_uniform p_uniform p_uniform_adj eta2_uniform uniform_DIF
#>    I1 group     0.368     0.544         1.000        0.001            
#>    I2 group     2.584     0.109         1.000        0.009            
#>    I3 group     0.308     0.579         1.000        0.001            
#>    I4 group     0.057     0.812         1.000        0.000            
#>    I5 group     0.081     0.776         1.000        0.000            
#>    I6 group     0.245     0.621         1.000        0.001            
#>  F_nonuniform p_nonuniform p_nonuniform_adj eta2_nonuniform nonuniform_DIF
#>         0.447        0.774            1.000           0.006               
#>         0.427        0.789            1.000           0.006               
#>         0.706        0.588            1.000           0.010               
#>         0.344        0.848            1.000           0.005               
#>         0.452        0.771            1.000           0.006               
#>         1.229        0.299            1.000           0.017               
#>  superseded p_uniform_boot_adj p_uniform_boot p_nonuniform_boot_adj
#>                          1.000          0.750                 1.000
#>                          0.750          0.250                 1.000
#>                          1.000          0.250                 1.000
#>                          1.000          1.000                 1.000
#>                          1.000          0.750                 1.000
#>                          1.000          0.750                 1.000
#>  p_nonuniform_boot n_boot_nonuniform n_boot_uniform uniform_DIF_boot
#>              1.000                 3              3                 
#>              0.500                 3              3                 
#>              0.750                 3              3                 
#>              1.000                 3              3                 
#>              1.000                 3              3                 
#>              0.750                 3              3                 
#>  nonuniform_DIF_boot
#>                     
#>                     
#>                     
#>                     
#>                     
#>                     
```
