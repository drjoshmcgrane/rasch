# rmt 1.1.0

* `dif_anova()` is now the single DIF analysis-of-variance function. One
  factor is analysed one-way; several factors are modelled jointly -- the
  statistically correct treatment -- with main effects by default and
  factor-by-factor interactions optional (`effects = "factorial"`). It
  handles within-subject factors (repeated-measures / mixed ANOVA) and
  returns a classed object with a `summary`, the full `terms` table, and
  Tukey comparisons. The separate `dif_anova_factorial()` is removed.
* `resolve_dif()` resolves DIF iteratively by item splitting, largest
  effect first, to clear artificial DIF.

# rmt 1.0.0

First stable release. The package delivers a complete Rasch Measurement
Theory workflow built entirely from published measurement theory, with a
pairwise conditional estimation core and a Shiny interface that exposes
every analysis with reproducing R code attached to every output.

## Models

* `rasch()`: dichotomous and polytomous (partial credit and rating scale)
  Rasch models by pairwise conditional maximum likelihood (Andrich & Luo
  2003; Zwinderman 1995), with Godambe sandwich standard errors, sum-zero
  identification, Warm (1989) weighted likelihood person estimation, and
  the principal-component threshold parameterisation as an estimation
  option.
* `rasch_mfrm()`: many-facet models with additive facet severities or
  item-by-facet interactions, accepting wide (one column per item) or long
  data.
* `rasch_efrm()`: the extended frame-of-reference model (Humphry) with
  frame units estimated by within-frame pairwise conditional likelihood.
* `btl()`: paired comparisons as the conditional form of the dichotomous
  Rasch model (Bradley & Terry 1952; Andrich 1978), including the
  adjacent-categories graded extension (Tutz 1986; Agresti 1992) with
  symmetric thresholds -- the Davidson (1970) ties model is its
  three-category case -- winner-plus-margin data entry, judge-clustered
  errors, and within-judge exposure and carry-over effects estimated
  jointly with the locations.

## Diagnostics

* The test-of-fit suite follows Andrich & Marais (2019): log-of-mean-square
  fit residuals with apportioned degrees of freedom, the item-trait
  chi-square on automatically sized, tie-preserving class intervals with
  per-interval detail, the ANOVA item fit, and person fit.
* Unidimensionality: residual principal components with parallel analysis,
  the Smith (2002) subset t-test, and magnitude estimation.
* Local dependence: Yen's (1984) Q3 with the Christensen, Makransky &
  Horton (2017) flagging convention, response-dependence magnitude
  (Andrich & Kreiner), subtest formation, and the spread test.
* Guessing: tailored analysis with origin-equated comparison calibrations.
* DIF: two-way residual ANOVA per factor, the joint factorial analysis
  with a compact uniform/non-uniform summary, resolved DIF magnitudes in
  logits with familywise control and a practical-significance criterion,
  planned contrasts derived from the factor structure (with
  repeated-measures support via person-level residual scores), and class
  intervals sized to the cells each analysis actually uses. For paired
  comparisons, `btl_dif()` tests object-by-judge-group interaction by the
  same two routes.
* Equating and invariance: common-item comparison with drift flags, and
  calibration anchoring.

## Display and reporting

* A complete base-graphics plot suite: expected value curves with observed
  class-interval means, category and threshold probability curves, the
  person-item threshold distribution, Wright maps, kidmaps (Wright, Mead &
  Ludlow 1980), object characteristic curves for paired comparisons,
  threshold and item maps, test characteristic and information curves,
  residual displays, and Guttman scalograms -- all with adjustable class
  intervals and scale ranges, and batch export to multi-page PDF or ZIP at
  publication resolution.
* `fit_summary_table()` and `targeting_table()` return the headline
  statistics as tidy tables; `save_outputs()` writes every table and plot
  to disk; `report_html()` produces a single-file HTML report.
* `run_app()` launches the Shiny interface: model-adaptive navigation,
  reproducing R code beneath every output, and one-click export.
