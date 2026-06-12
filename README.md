# RaschR

A pairwise conditional Rasch measurement engine in R, built entirely from
published measurement theory. Item estimation is pairwise conditional maximum
likelihood (Andrich & Luo 2003; Zwinderman 1995); person locations use Warm's
weighted likelihood; the diagnostic suite covers fit, reliability, targeting,
dimensionality, local dependence, and DIF over multiple person factors; and a
modern Shiny interface exposes the lot with one-click export of every table
and plot.

## Install

```r
# install.packages("remotes")
remotes::install_github("drjoshmcgrane/RaschR")
```

The analysis engine is base R only (`stats`, `graphics`, `grDevices`,
`utils`). The Shiny GUI additionally needs `shiny`, `bslib`, and `DT`.

## Quick start

```r
library(RaschR)

# data frame with an ID column, item columns (names carry through), and
# person factors for DIF
fit <- rasch(responses, model = "PCM",
             id = "person_id", factors = c("gender", "site"))

summary(fit)         # full test-of-fit report
fit$items            # location + SE, fit residual, infit/outfit, chi-square
fit$person           # WLE + SE per person, with ID, factors, raw score, person fit
fit$thresholds       # every threshold with its standard error
score_table(fit)     # raw score to measure conversion
fit$psi; fit$psi_noext; fit$alpha   # PSI (with/without extremes), Cronbach's alpha

plot_icc(fit, "Q05")                  # ICC with observed class-interval means
plot_icc(fit, "Q05", group = "gender")  # ICC by group: the graphical DIF display
plot_ccc(fit, "Q05")                  # category probability curves
plot_threshold_prob(fit, "Q05")       # threshold probability ogives
plot_pimap(fit)                       # person-item threshold distribution
plot_threshold_map(fit)               # all thresholds by item
plot_tcc(fit); plot_tif(fit)          # test characteristic and information curves
plot_item_map(fit); plot_person_fit(fit)
plot_resid_cor(fit); plot_pca(fit)

dimensionality_test(fit)   # residual-PCA contrast + person t-test
residual_correlations(fit) # local dependence
dif_anova(fit)             # uniform and non-uniform DIF per factor (BH-adjusted)
dif_anova_factorial(fit)   # all factors jointly: factorial ANOVA per item,
                           # interaction-precedence, Tukey HSD on group terms

# local dependence remedy: combine dependent items into a subtest and refit
fit2 <- combine_items(fit, list(c("Q04", "Q05")))

# DIF remedy: split the offending item by the factor and refit; each group
# receives its own location and the distance between them is the DIF size
fit3 <- split_items(fit, "Q05", by = "gender")

# equating: anchor nominated thresholds (or item means, with k = NA) at
# fixed values; everything else is estimated on the anchored scale and
# person measures become comparable across separately analysed datasets
fit_eq <- rasch(responses, anchors = data.frame(item = c("Q01", "Q10"),
                                                k = c(1, NA), tau = c(-1.2, 0.9)))

# equating check: compare two calibrations through their common items
eq <- equate_tests(fit, fit_eq)   # or a bank: data.frame(item, location, se)
eq$table                          # t-tests against the shifted identity line
plot_equate(fit, fit_eq)

# many-facet Rasch model for rated, long-format data (one row per response)
mf <- rasch_mfrm(ratings, person = "person", item = "criterion",
                 score = "score", facets = "rater")
mf$facet_effects$rater    # severities with SEs and pooled fit per rater
plot_facets(mf)           # severity caterpillar plot

# interactive facet mode: lets a rater be severe on particular items
mfi <- rasch_mfrm(ratings, person = "person", item = "criterion",
                  score = "score", facets = "rater", interaction = "rater")
mfi$interaction_effects   # item-by-rater gamma with SEs

# multiple choice: score raw responses against a key, keep them for
# distractor analysis (rest-measure based, with miskey flagging)
mc <- rasch(responses_raw, key = c(Q1 = "A", Q2 = "C", Q3 = "B"))
distractor_analysis(mc)    # per option: n, proportion, mean location, point-biserial
plot_distractors(mc, "Q2") # option curves over class intervals

# extended frame of reference model: the unit of the scale differs across
# item-set by person-group frames (rho = alpha_set x phi_group)
ef <- rasch_efrm(responses, item_sets = list(numeracy = num_items,
                                             literacy = lit_items),
                 groups = "year_group")
ef$phi_table; ef$alpha_table   # group and set units with SEs
ef$frames                      # one row per frame: unit, origin, pooled fit
plot_frames(ef)                # unit caterpillar
plot_icc_frames(ef, "Q07")     # fanned ICCs: parallel within, fanned across

# compare model fits: same-data fits get the pairwise conditional
# log-likelihood difference (descriptive); all fits get calibration-free
# fit descriptors (chi-square/df, fit residual SDs, PSI, alpha)
compare_fits(PCM = fit, RSM = rasch(responses, model = "RSM"))

save_outputs(fit, "results/")   # every table (CSV), every plot (PNG + PDF), summary.txt
```

`responses` is a persons-by-items matrix, or a data frame that also carries an
ID column and person factors. Item scores are integers starting at 0; missing
data is allowed (estimation is pairwise, person measures use each person's
observed items). Empty categories are collapsed and constant items dropped,
with notes recorded on the fit.

## GUI

```r
RaschR::run_app()
```

A bslib (Bootstrap 5) interface: upload a CSV/TSV, nominate the ID variable,
any number of person factors, and the item columns, then run. Long-format
(rated) data is supported too: nominate the person, item, score, and facet
columns (optionally with an item-by-facet interaction) and the app fits the
many-facet model, with a Facets tab for severities, caterpillar plots, and
interaction effects. Anchor values can be uploaded as a CSV for equating
(blank `k` rows anchor item means), an Equating tab compares the current
calibration against an uploaded reference with drift tests and a plot,
locally dependent items can be combined into a subtest from the Local
dependence tab, and DIF can be resolved by splitting the flagged item from
the DIF tab. A third data layout fits the extended frame of reference model:
nominate the person-group column and an item-set map (CSV upload or
inference from item-name prefixes), and a Frames tab reports the unit
tables, the frame caterpillar plot, fanned cross-frame ICCs, and the
equal-unit comparison. The other tabs cover the summary
(PSI, alpha, power of fit, targeting, score table), item statistics with
per-item curves, person estimates, test-level plots, multi-factor DIF with
ICC-by-group plots, dimensionality, and local dependence. Every plot has PNG
and PDF download buttons, every table a CSV button, and the Export tab
downloads the entire analysis as a ZIP archive.

## What it implements

- Pairwise conditional maximum likelihood (Andrich & Luo 2003; Zwinderman
  1995): the person parameter cancels within every item pair and the
  conditional likelihood is maximised by Newton-Raphson (`pcml`). Standard
  errors come from a Godambe sandwich estimator, which corrects the
  over-optimism of the naive pairwise information.
- Partial credit (PCM) and rating scale (RSM, via a constrained design
  matrix) models; dichotomous data is the special case.
- Warm (1989) weighted likelihood person estimates, finite at extreme scores,
  computed per missing-data pattern; extremes flagged.
- Item and person fit residuals (Wilson-Hilferty standardised mean squares
  with a degrees-of-freedom correction for the estimated person locations),
  infit and outfit; item-trait interaction chi-square over class intervals
  with Benjamini-Hochberg adjustment and an optional sample-size adjustment.
- Person Separation Index with and without extreme persons, Cronbach's alpha,
  targeting, the power-of-test-of-fit assessment, the score-to-measure table,
  and the test information function.
- Threshold and category diagnostics (ordering, reversals, never-modal
  categories, category frequencies).
- Residual-PCA dimensionality test (Smith person t-test), residual-correlation
  local dependence, and DIF by two-way residual ANOVA over any number of
  person factors.
- Anchored estimation for test equating: individual thresholds, or item mean
  locations with free thresholds (average anchoring, `k = NA`), are held at
  known values and the rest of the calibration, including all person
  measures, lands on the anchored scale.
- Common-item equating tests (`equate_tests`, `plot_equate`): two
  calibrations (or a calibration against an item bank) are compared through
  their common items with t-tests against the shifted identity line;
  drifting items are flagged.
- The many-facet Rasch model (Linacre 1989) by the same pairwise conditional
  likelihood: item-by-facet combinations are calibrated jointly with the
  facet severities entering the design matrix; severities are reported with
  sandwich standard errors and pooled fit statistics, and the result is a
  full diagnostic object. Interactive facet mode (`interaction =`) adds
  item-by-facet terms so a level may be more or less severe on particular
  items.
- Structure amendments by re-analysis: subtests (`combine_items`) merge
  locally dependent items into one polytomous super-item; item splitting
  (`split_items`) resolves DIF by giving each group level its own copy of
  the offending item.
- The Guttman scalogram (`guttman_table`, `plot_guttman`) with the
  coefficient of reproducibility.
- Multiple-choice support: raw responses scored 0/1 against a key
  (`rasch(..., key = )`), with the raw options retained for
  `distractor_analysis` (count, proportion, mean rest-measure location, and
  point-biserial per option, miskeys flagged) and `plot_distractors` option
  curves. The rest measure excludes the analysed item, so the keyed option
  cannot credit its own takers.
- The extended frame of reference model (`rasch_efrm`; Humphry 2005;
  Humphry & Andrich 2008) — to our knowledge its first software
  implementation. Frames are item-set by person-group cells with units
  `rho = alpha_set x phi_group`. Person-group units come from person-free
  within-frame pairwise conditioning (a set taken by two groups shows the
  same threshold pattern at two scales); item-set units come from persons
  common to the sets via error-corrected true-score variance ratios,
  reconciled over the linking graph. Person measures use weighted-score
  sufficiency, and per-group score curves replace the raw-score table.
  Includes frame fit pooling, an equal-unit model comparison, the frame-unit
  caterpillar plot, and fanned cross-frame ICCs. Humphry (2005) states the
  model for dichotomous responses; the polytomous form fitted here (the
  frame unit multiplying the whole exponent over the item's partial-credit
  thresholds) is this package's extension of that statement, characterised
  by preserving within-frame partial-credit structure — hence the validity
  of the pairwise conditional cancellation — and weighted-score
  sufficiency, and reducing exactly to the dichotomous model and the
  ordinary PCM in the corresponding special cases.
- `save_outputs()` writes the complete analysis (all tables, all plots, a
  text summary) to a folder in one call.

## Methodological references

Andrich & Luo (2003) and Zwinderman (1995) conditional pairwise estimation;
Warm (1989) weighted likelihood; Andrich item-trait interaction chi-square;
Smith (2002) residual-component dimensionality; Hagquist & Andrich (2017)
DIF by residual ANOVA; Christensen, Makransky & Horton (2017) residual
correlations for local dependence; Linacre (1989) many-facet Rasch
measurement; Humphry (2005) and Humphry & Andrich (2008) on the unit in
Rasch measurement and the extended frame of reference; Godambe sandwich
covariance for composite likelihoods.

## Validation

Every estimation and diagnostic component is validated against simulated data
with known parameters in `tests/testthat`: parameter recovery for dichotomous,
PCM, and RSM data; sandwich standard errors checked against empirical sampling
variability; DIF detected on planted items only; dimensionality verdicts on
one- and two-dimensional data; missing-data person estimation; and the full
export. `R CMD check` runs clean. The fit residual is an information-weighted
standardised mean square (Wilson and Hilferty transformation), approximately
N(0,1) under fit: negative values indicate overfit, positive values underfit.
