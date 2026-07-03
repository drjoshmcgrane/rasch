# rmt

A pairwise conditional Rasch measurement engine in R, built entirely from
published measurement theory. Item estimation is pairwise conditional maximum
likelihood (Andrich & Luo 2003; Zwinderman 1995); person locations use Warm's
weighted likelihood; the diagnostic suite follows the conventions set out
in Andrich & Marais (2019) —
the log-of-mean-square fit residual with apportioned degrees of freedom,
automatically sized class intervals, the item-trait chi-square with its
per-interval detail table, the ANOVA item-fit F — and covers reliability, targeting,
dimensionality (with magnitude estimation), local dependence (with magnitude
estimation), guessing (tailored analysis), and DIF over multiple person
factors; and a modern Shiny interface exposes the lot with one-click export
of every table and plot.

## Install

```r
# install.packages("remotes")
remotes::install_github("drjoshmcgrane/rmt")
```

The analysis engine is base R only (`stats`, `graphics`, `grDevices`,
`utils`). The Shiny GUI additionally needs `shiny`, `bslib`, and `DT`.

## Quick start

```r
library(rmt)

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

chisq_detail(fit, "Q05")   # per-class-interval chi-square detail table
score_table(fit, method = "mle", extremes = "extrapolated")
                           # complete-data estimates with geometric
                           # extreme-score extrapolation, frequencies, cum %
ctt_table(fit)             # traditional statistics: facility, item-total r,
                           # discrimination index, alpha, classical SEM
lr_test(fit)               # PCM vs rating parameterisation: the raw
                           # composite chi-square + Kent-calibrated version

dimensionality_test(fit)   # residual-PCA contrast + person t-test + paired t
dimensionality_magnitude(fit, list(setA, setB))
                           # Andrich (2016): c, latent subscale correlation
                           # rho = 1/(1+c^2), common-variance proportion A
residual_correlations(fit) # local dependence
dependence_magnitude(fit, dependent = "Q05", independent = "Q04")
                           # Andrich & Kreiner d in logits, with SE and test
dif_anova(fit)             # uniform and non-uniform DIF per factor: the full
                           # two-way table (group, class interval, interaction)
                           # with partial eta-squared effect sizes; BH, Holm,
                           # or Bonferroni adjustment across items
dif_anova_factorial(fit, sizes = TRUE)
                           # all factors jointly: x-way factorial ANOVA per
                           # item with interactions, interaction-precedence,
                           # Tukey HSD post-hocs on significant group terms
                           # (interaction cells included), and DIF magnitudes
                           # in logits for every significant term
dif_size(fit, "Q05", by = "gender")
                           # DIF magnitude on the measurement scale: resolved
                           # locations per level (or interaction cell), all
                           # pairwise differences with Holm familywise
                           # adjustment and a practical-significance flag
                           # (default 0.5 logits)
tailored_analysis(fit, chance = 0.25)
                           # four-step guessing procedure (dichotomous MC)

# local dependence remedy: combine dependent items into a subtest and refit,
# then screen the super-items' spread against the binomial bound
fit2 <- combine_items(fit, list(c("Q04", "Q05")))
spread_test(fit2)          # spread below the LUB indicates dependence

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

# double keying: several correct options, all scoring 1
mc2 <- rasch(responses_raw, key = c(Q1 = "A/C", Q2 = "C", Q3 = "B"))

# polytomous option scoring (Andrich & Styles 2011): informative
# distractors receive partial credit and the item is fitted as polytomous;
# distractor_rescore() proposes the scoring from the rest-measure evidence
prop <- distractor_rescore(mc)
prop$option_scores         # item, option, score - review, edit, refit:
mc3 <- rasch(responses_raw, key = prop$option_scores)

# extended frame of reference model: the unit of the scale differs across
# item-set by person-group frames (rho = alpha_set x phi_group)
ef <- rasch_efrm(responses, item_sets = list(numeracy = num_items,
                                             literacy = lit_items),
                 groups = "year_group")
ef$phi_table; ef$alpha_table   # group and set units with SEs
ef$frames                      # one row per frame: unit, origin, pooled fit
plot_frames(ef)                # unit caterpillar
plot_icc_frames(ef, "Q07")     # fanned ICCs: parallel within, fanned across

# paired comparisons: the Bradley-Terry-Luce model, the conditional
# (person-free) form of the dichotomous Rasch model, estimated by the same
# conventions -- sum-zero identification, sandwich SEs (judge-clustered),
# fit residuals for objects and judges, pairwise goodness of fit
bt <- btl(comparisons, object_a = "left", object_b = "right",
          winner = "preferred", judge = "judge")
bt$objects; bt$judges     # locations + fit; erratic judges flagged
plot_btl(bt)              # object caterpillar

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
rmt::run_app()
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
- An optional alternative estimator (`pcml_pc`, or `rasch(pc_components =)`
  end to end, with ranks from "location only" through "equal spread" (the
  dispersion model of Andrich 1982) to the full four-component model), for
  use when some
  categories are sparsely populated, reparameterises each item's thresholds
  as Andrich's (1978, 1985) orthogonal-polynomial principal components —
  location, spread, skewness, and kurtosis (Pedler 1987). For items with up
  to 3 thresholds this recovers the same thresholds and log-likelihood as
  the free partial credit model (the default, via `pcml`); the component
  family stops at the quartic (kurtosis) term, so items with more thresholds
  are necessarily smoothed to that four-parameter polynomial trend across
  categories, which can stabilise thresholds when some categories are
  sparsely observed.
- Warm (1989) weighted likelihood person estimates, finite at extreme scores,
  computed per missing-data pattern; extremes flagged. The score-to-measure
  table also offers plain MLE and the geometric extrapolation of the
  extreme-score measures described by Andrich & Marais (2019, ch. 10) (the
  top difference d solves b = sqrt(ad) on the
  two preceding differences), verified against the worked example there;
  `person_extrapolated()` applies the same rule to
  the person table, per missing-data pattern, so extreme persons can enter
  group comparisons.
- The fit residual exactly as derived in Andrich & Marais (2019, ch. 23):
  squared standardised residuals summed
  over the observed cells of non-extreme persons, compared with equally
  apportioned model-testing degrees of freedom, symmetrised by the
  log-of-mean-square transform with model-based variance sum(C4/V^2 - 1).
  The untransformed natural residual, its degrees of freedom, and
  infit/outfit (Wilson-Hilferty z) are reported alongside, with distribution
  summaries (mean, SD, skewness, kurtosis), fit-location correlations, and
  the cell df factor, as in the summary blocks of Andrich & Marais (2019,
  app. C).
- The item-trait interaction chi-square over automatically sized class
  intervals (at least 50 non-extreme persons per interval, at most 10, at
  least 2, by default; Andrich & Marais 2019, ch. 15), per-item degrees of freedom, Benjamini-Hochberg and Bonferroni
  adjustments, an optional sample-size adjustment, the per-item
  class-interval detail printout (`chisq_detail`: size, location max/mean,
  residual, chi-square component, OM/EV, effect size ES, per-category
  OBS.P/EST.P/OBS.T), and the class-interval ANOVA item-fit F.
- Person Separation Index with and without extreme persons, the item
  separation index, Cronbach's alpha (flagged not-applicable with missing
  data), targeting, the power-of-test-of-fit assessment, the
  score-to-measure table with complete-responder frequencies, the test
  information function, and classical-test-theory companion statistics
  (`ctt_table`: facility, item-total and item-rest correlations, the
  upper-lower discrimination index, alpha-if-deleted, classical SEM).
- A likelihood-ratio test of the partial credit against the rating
  parameterisation (`lr_test`), reporting both the conventional raw
  pairwise-composite chi-square and a first-order calibrated version (Kent
  1982; Varin, Reid & Firth 2011) whose eigenvalue adjustment comes from
  the same Godambe matrices as the sandwich standard errors — simulation
  shows the raw test is severely anticonservative while the calibrated one
  holds size.
- Threshold and category diagnostics (ordering, reversals, never-modal
  categories, category frequencies).
- Residual-PCA dimensionality test (Smith person t-test, plus the paired
  t-test of subset means), residual-correlation local dependence, and a
  complete DIF procedure extending the single-factor residual ANOVA of
  Hagquist & Andrich (2017): the full two-way table per item (group, class
  interval, group-by-interval) with partial eta-squared effect sizes and a
  choice of false-discovery-rate or familywise adjustment across items;
  with several person factors, an x-way factorial per item including all
  factor-by-factor interactions, where a significant interaction supersedes
  the terms it contains, significant interactions receive Tukey HSD
  post-hoc comparisons over their cells, and significant main effects with
  more than two levels receive Tukey post-hocs over their levels (Tukey's
  own familywise control within each family); and DIF magnitude on the
  logit scale (`dif_size`) by the resolved-item method of Andrich & Marais
  (2019, ch. 16), with the full sandwich covariance between resolved
  locations, Holm familywise adjustment over the pairwise comparisons, and
  a practical-significance criterion (default half a logit) so statistical
  and practical DIF are judged separately.
- Magnitude estimation for both violations of independence, not just flags:
  `dependence_magnitude` implements the Andrich & Kreiner (2010) resolution
  method (polytomous form Andrich, Humphry & Marais 2012) with the pooled
  threshold-variance standard error; `spread_test` screens subtest spread
  parameters against Andrich's (1985) binomial least upper bounds; and
  `dimensionality_magnitude` implements Andrich's (2016) two-calculation
  reliability comparison, returning
  the latent subscale correlation rho = 1/(1+c^2) and the common-variance
  proportion A — verified against the worked example in Andrich & Marais
  (2019, Table 24.4).
- Tailored analysis for guessing (`tailored_analysis`): the four-step
  procedure of Andrich, Marais & Humphry (2012) — initial, Waller-cutoff
  tailored, average-anchored origin-equated, and all-anchored person
  re-estimation — with the guessing signature summarised per item.
- Racked and stacked reshaping of repeated measurements (`rack_data`,
  `stack_data`) for change-in-items and change-in-persons designs, the
  latter enabling DIF-over-time with time as a person factor.
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
  sandwich standard errors and pooled fit statistics — the
  log-of-mean-square fit residual summed over the level's cells with its
  degrees of freedom, alongside the facet-margin mean of the virtual items'
  fits reported in the three-facet tables of Andrich & Marais (2019) — and
  the result is a
  full diagnostic object. Interactive facet mode (`interaction =`) adds
  item-by-facet terms so a level may be more or less severe on particular
  items.
- Structure amendments by re-analysis: subtests (`combine_items`) merge
  locally dependent items into one polytomous super-item; item splitting
  (`split_items`) resolves DIF by giving each group level its own copy of
  the offending item.
- The Guttman scalogram (`guttman_table`, `plot_guttman`) with the
  coefficient of reproducibility.
- The Bradley-Terry-Luce model for paired comparisons (`btl`,
  `plot_btl`) as a member of the same family: within an item pair, given
  one correct response, the Rasch probability that it was the easier item
  is exactly of BTL form (Andrich 1978), and the package's pairwise
  conditional estimation maximises precisely such a likelihood — verified
  by a test in which `btl()` on the pair-conditional comparisons extracted
  from Rasch data reproduces `pcml()`'s item locations to solver
  tolerance. Estimated by the package's standard conventions (conditional
  ML, sum-zero identification, Godambe sandwich standard errors clustered
  by judge when judges are identified), with log-of-mean-square fit
  residuals for
  objects and judges (an erratic judge flags exactly as an erratic person
  does), the classical pairwise goodness-of-fit chi-square, an object
  separation index, tie handling, extreme-object removal, and a
  connectivity check on the comparison graph.
- Multiple-choice support: raw
  responses scored against a key (`rasch(..., key = )`) as 0/1, with double
  keying (`"A/C"` credits both options), or with full polytomous option
  scoring (a data frame `item, option, score`) so that informative
  distractors receive partial credit and the item is fitted as polytomous
  (Andrich & Styles 2011); `distractor_rescore()` proposes such a scoring
  from the rest-measure evidence (separation of a distractor's takers from
  the uncredited ones), for substantive review before refitting. Raw
  options are retained for `distractor_analysis` (count, proportion, mean
  rest-measure location, and point-biserial per option, miskeys flagged)
  and `plot_distractors` option curves, which mark credited distractors.
  The rest measure excludes the analysed item, so the keyed option cannot
  credit its own takers.
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

## Measurement-theoretic status

Everything in the package that claims to be Rasch measurement is
conditionally estimated, sufficiency-respecting, and invariance-preserving:
item comparisons are person-free by pairwise conditioning, and the partial
credit, rating scale, and additive many-facet models are members of the
Rasch class. The diagnostics (fit residuals, item-trait chi-square,
residual components, DIF analysis of variance, threshold ordering) exist to
police the theory's requirements, and the structural remedies (subtests,
item splitting, anchoring) are orthodox practice that restore rather than
parameterise away invariance.

Departures from the classical model are deliberate and labelled. Warm's
weighted likelihood adds a penalty beyond the conditional likelihood (this
is what makes extreme-score estimates finite). Cronbach's alpha and the
distractor point-biserials are classical-test-theory companions, reported
as descriptives only. Interactive facet mode remains in the Rasch class but
a significant item-by-facet interaction qualifies specific objectivity in
practice: comparisons of the interacting facet's levels become
item-dependent. The extended frame of reference model is strictly Rasch
within every frame; across frames it is an argued extension of the theory
of the unit (Humphry 2005; Humphry & Andrich 2008) whose status the
literature still debates, and its item-set units are necessarily identified
from the person side (we show conditional identification is impossible), a
departure from purely distribution-free comparison that belongs to the
model, not the implementation. The free-slope model used in the case study
is a diagnostic of the frame restriction, never a measurement model.

## Case study: wording effects as frame units

`inst/casestudies/wording_units_selfesteem.R` applies the extended frame of
reference model to the public Rosenberg Self-Esteem Scale dataset of the
Open Source Psychometrics Project (downloaded from the source at run time),
treating positively and negatively worded items as two item sets. The
positively worded items carry a unit about 27 per cent larger than the
reverse-scored negative items (alpha ratio 1.266, 95 per cent CI 1.246 to
1.286; Wald p < 1e-89), and it matters: persons with identical raw scores
differ by up to 0.85 logits once the wording units are modelled, over half
of all respondents move by more than 0.1 logits, and the male-female gap is
understated by about 10 per cent under equal units. A free-slope
generalized partial credit model agrees at the set level (geometric-mean
slope ratio about 1.25), and sensitivity analysis localises part of the
effect: the unit ratio is 1.115 without the scale's well-known ambivalent
item Q8 (1.196 also dropping the weakest positive item), still
significantly above one in every configuration. With a single person group
the pairwise equal-unit comparison is invariant to set units by
construction, so the Wald tests on the log units carry the evidence -- the
fit object reports both and says which applies.

## Methodological references

Andrich & Luo (2003) and Zwinderman (1995) conditional pairwise estimation;
Andrich (1978, 1985) and Pedler (1987) principal-components reparameterisation
of thresholds; Warm (1989) weighted likelihood; Andrich item-trait interaction
chi-square;
Smith (2002) residual-component dimensionality; Andrich (2016)
multidimensionality magnitude from two calculations of reliability;
Andrich & Kreiner (2010) and Andrich, Humphry & Marais (2012) response
dependence magnitude; Andrich (1985) spread-parameter bounds; Waller (1989)
and Andrich, Marais & Humphry (2012) tailored analysis of guessing;
Hagquist & Andrich (2017)
DIF by residual ANOVA; Christensen, Makransky & Horton (2017) residual
correlations for local dependence; Linacre (1989) many-facet Rasch
measurement; Humphry (2005) and Humphry & Andrich (2008) on the unit in
Rasch measurement and the extended frame of reference; Godambe sandwich
covariance for composite likelihoods; Kent (1982) and Varin, Reid & Firth
(2011) composite likelihood-ratio calibration; Andrich & Marais (2019) for
the output conventions followed throughout.

## Validation

Every estimation and diagnostic component is validated against simulated data
with known parameters in `tests/testthat`: parameter recovery for dichotomous,
PCM, and RSM data; the principal-components reparameterisation checked against
the free partial credit model for items with few thresholds, against known
locations at reduced rank, and against its own rank ceiling (kurtosis
identified from 5 thresholds but not at exactly 4); sandwich standard errors
checked against empirical sampling
variability; DIF detected on planted items only; dimensionality verdicts on
one- and two-dimensional data; missing-data person estimation; and the full
export. The published conventions are validated directly: the fit-residual
degrees of freedom reproduce the (N-1)(I-1) apportionment on complete
data, the geometric extreme-score extrapolation reproduces the worked
example in Andrich & Marais (2019, ch. 10; 4.762), the dependence
magnitude d recovers planted
simulation values, the dimensionality block reproduces the algebra of
Andrich (2016), and the calibrated likelihood-ratio test holds size on
rating-scale-true data where the raw composite statistic rejects wildly.
`R CMD check` runs clean. The fit residual is the log-of-mean-square
statistic of Andrich & Marais (2019, ch. 23), approximately N(0,1) under
fit: negative values indicate
overfit (too-Guttman responding), positive values underfit (erratic
responding).
