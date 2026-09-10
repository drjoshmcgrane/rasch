# Simulation validation battery

Version note: development builds 1.12.1 through 1.14.2 were never
released; their content ships publicly as the 1.12.0 CRAN update. The
version numbers and commit hashes below refer to those internal
development states and remain the accurate provenance record.

Scripts and results for the release validation batteries summarised in the
plant-and-detect vignette ("Validation studies"). Nothing here ships in the
CRAN package (see `.Rbuildignore`); everything runs from the package root
with `Rscript`, loading the in-tree package via `pkgload::load_all(".")`.

## Layout

- `tests/testthat/test-format.R` checks that missing DIF probabilities and
  zero standard errors cannot produce an ETS A/B classification. It also
  checks `dif_size()` with a degenerate resolved covariance. Classification
  is unchanged for supported tests; this is a reporting conformance check,
  not a new calibration study.
- `studies/frame-link-safety.R` checks complete and quasi-complete CJ frame
  separation, finite near-separated links, iteration-limit accounting, and
  refusal of unverified saved frame calibrations. Project checks cover base,
  history, kept and nested-result fits, with source-file preservation and an
  app restore regression. `results/frame-link-safety.csv` records the runner,
  source-tree and test hashes. These are fixed-data conformance checks, not
  estimates of coverage or error rates.
- `studies/frame-likelihood-curvature.R` checks stationary saddles in response
  EFRM and CJ frame calibration against independently evaluated likelihoods
  and numerical Hessians. It checks scale invariance of the curvature guard,
  public refusal, fixed-unit estimation, bootstrap exclusion and covariance
  reconstructed from accepted draws. Results in
  `results/frame-likelihood-curvature.csv` identify the source tree, runner
  and regression-file hashes. These are deterministic numerical checks, not
  coverage or Type I error studies. Historical frame validation predates this
  guard and does not measure its refusal rate; saved analyses need refitting
  to apply it.
  `test-frame-refits.R` also checks a threshold-only flat direction after
  resolving every item in a set. Retaining one common item restores the
  origin link when another set already identifies the group-unit ratio.
- `studies/efrm-full-identification.R` checks full-bootstrap group-unit
  identification on five fixed datasets: balanced binary and polytomous
  scales, groups of three and twelve persons alongside a group of 120,
  and a two-set, two-group design. The recording wrappers do not alter
  solver output. Unidentified draws must be rejected, and every reported
  covariance block must match the accepted draws. Results are in
  `results/efrm-full-identification.csv`, with execution-state hashes:
  190 of 200 draws were retained and ten unidentified draws were rejected.
  The rerun with the likelihood-curvature guard retained the same draws;
  one rejected draw also failed the convergence check.
  All covariance reconstructions matched exactly. The two small-group
  cases still withhold group-unit probabilities.
  These are numerical and accounting checks, not coverage or Type I error
  studies. Historical full-bootstrap results predate this rejection rule;
  analyses that included unidentified unit draws need refitting. The
  separate `test-efrm-full-bootstrap-identification.R` checks controlled
  failures, the usable-draw floor and the hybrid fallback, without requiring
  analytic standard errors from an otherwise identified bootstrap refit.
- `tests/testthat/test-btl-pair-tie-invariance.R` checks that tied CJ object
  locations cannot acquire a directional surprise through relabelling.
  Tied matchups retain their descriptive probabilities and their place in
  the planned Holm family. The report-template regression in
  `test-export-safety-and-external-dif.R` checks that externally specified
  DIF and its bootstrap sensitivity analysis appear in rendered reports.
- `tests/testthat/test-efrm-link-graph.R` checks indirect set linking on a
  three-set chain, insufficient direct overlap, disconnected designs and
  numerical-link refusal. A failed supported edge invalidates the whole
  bootstrap replicate, even when two other edges connect all three sets.
  Fault-injection checks verify numerical failure, non-convergence, usable
  counts, the minimum-draw refusal and covariance computed from accepted
  draws only. The internal moment-link fallback follows the same rule.
  Supported direct-link estimates are unchanged.
  A fixed-seed case exercises both hybrid and full-bootstrap uncertainty
  with 40 requested replicates and retains failed-replicate accounting.
  `studies/efrm-link-graph-conformance.R` adds eight end-to-end cases: two
  sets, a complete three-set design, pairwise booklets and a three-set chain,
  each at unit ratios 1 and 1.3 with 600 persons and five items per set.
  All fits converged; 465 of 480 requested hybrid draws were usable and 15
  failed. Results are in `results/efrm-link-graph-conformance.csv`, with
  execution-state source and script hashes. Subsequent help-text edits change
  the R-tree hash but not the estimators.
  These are graph, numerical and accounting checks, not Type I error or
  coverage studies. Historical multi-set hybrid results predate the corrected
  failure handling and do not establish its calibration; affected analyses
  need refitting. Two-set scenarios already rejected a failed sole link.
  `test-structural-factor-roles.R` checks factor identity through structural
  refits, including the frame-invariance bootstrap. The app simulation-bundle
  test checks that exported code recreates the responses and truth metadata.
- `tests/testthat/test-factor-item-name-disambiguation.R` checks explicit
  item selectors with separately supplied, same-named person factors against
  matrix input. `test-compare-btl-pairwise-identity.R` checks that increasing
  within-judge sequence relabelling preserves model comparability while real
  order changes do not. These change input handling, not estimation.
  `test-export-safety-and-external-dif.R` checks exact DIF table exports and
  refusal of non-empty output folders without changing existing files.
- `tests/testthat/test-dif-factor-handoff.R` checks automatic follow-ups
  against directly supplied factor values and refuses obsolete saved tables
  before report creation. The app-project tests check omission of those
  results and CJ DIF whose judge role cannot be verified, while preserving
  the data and fitted models. These are compatibility checks; affected DIF
  analyses must be rerun. `test-sim-recovery-boundaries.R` checks BTL recovery
  against pre-boundary observations, including refusal of changed source
  rows. Extreme objects remain excluded from recovery summaries. Neither
  correction changes the fitted estimators or simulation generators.
- `tests/testthat/test-parallel-rng.R` checks bootstrap worker parity under
  Mersenne-Twister, L'Ecuyer-CMRG, Box-Muller and the older Rounding sampler.
  It requires an installed package and actual socket workers; a serial
  fallback cannot count as a passing parallel check. A public item-fit
  bootstrap comparison also verifies replicate arrays, item and person
  results, and preservation of the caller's RNG state. This corrects seeded
  reproducibility under non-default generators, not the bootstrap null or
  estimator. Earlier default-generator calibration studies remain applicable;
  earlier non-default parallel runs need rerunning to reproduce their serial
  counterparts.
- `studies/joint-dif-wle-maxima.R` checks the revised incomplete-panel
  adjustment and competing WLE maxima. `SV_PART=dif`, `public`, and `wle`
  select its components; an unset `SV_PART` runs all three. Results are in
  `joint-dif-residual-model.csv`, `joint-dif-public-screen.csv`, and
  `wle-separated-bank-maxima.csv`. The principal known-residual null uses
  3,000 replicates: uniform/non-uniform rejection is 5.7/5.3%, with MC SEs
  0.42/0.41 percentage points. Secondary null cells use 300 replicates;
  two effect sizes per alternative use 200. The separate 200-dataset
  end-to-end Rasch check gave 3.0% rejection for the null A effect on an
  item carrying B DIF (MC SE 1.21 percentage points). These are marginal
  primary-test checks, not all-item familywise or bootstrap calibration.
  All requested replicates were usable. No dense-grid comparison point
  exceeded the selected WLE objective across all scores in 100 dichotomous
  and PCM banks with easy/hard gaps up to 16 logits. This is numerical
  conformance, not a confidence-interval coverage study or a proof about
  arbitrary banks. Rows identify the script hash and execution-state R-tree
  hash; later documentation and app changes need not share that tree hash.
  Earlier incomplete-panel validation used marginal centering and does not
  validate the revised joint adjustment. Complete-panel references are
  unchanged. The new regression files also check zero-unit BTL-EFRM refusal,
  pairwise comparison identity, stable explanatory candidates, app source
  changes, saved-result migration and conditional-bootstrap replay.
- `tests/testthat/test-btl-dimensionality-availability.R` checks that
  incomplete pair coverage and shared-order confounding withhold the whole
  CJ dimensionality reference, not just its categorical flag. It covers
  ordinary and extended-frame comparisons, printed output, plots and stale
  saved results. Descriptive decompositions and completed draws are retained.
  The conditional simulation method is unchanged; previously displayed
  probabilities from these unsupported designs should not be interpreted.
  `test-btl-dimensionality-project-migration.R` checks authenticated omission
  of the affected saved references without losing data, models or history.
  Supported older references remain usable; tampered bundles and results
  bound to a different fit are refused.
- `tests/testthat/test-explanatory-coefficient-labels.R` checks coincident
  model-matrix labels and fixed-departure names against equivalent models
  with unambiguous predictor names. Rasch and CJ estimates, diagnostics and
  nested-model comparisons agree. This changes naming, not estimation.
- `tests/testthat/test-app-cj-dif-code.R` executes the app's emitted CJ DIF
  code after judge metadata edits and after reopening a saved analysis.
  It must reproduce the completed run's factor maps and summary, even when
  the calibration's source data contain earlier judge metadata.
  Reopened DIF runs retain the corrected assignments without changing the
  calibration data; a new upload can replace the restored metadata.
- `tests/testthat/test-legacy-refit-restrictions.R` checks that structural
  refits, tailored analysis and bootstrap procedures cannot release anchor or
  principal-component restrictions when replay settings are unavailable.
  It includes average anchoring, unchanged supported refits and recovery on
  the identified origin. These are model-preservation checks; estimators and
  generators are unchanged. Analyses previously recalibrated from such older
  constrained fits should be rerun from the source data and original settings.
- `tests/testthat/test-lr-fitted-restrictions.R` checks that missing replay
  settings cannot remove an older fit's anchor or principal-component
  restrictions from the `lr_test()` eligibility check. Unrestricted legacy
  fits retain the same comparison. The comparison method and its simulation
  calibration are unchanged.
- `tests/testthat/test-app-export-selected-dif.R` exercises the app's actual
  report and archive handlers when the requested factorial DIF analysis is
  unavailable. No default main-effects result may replace it. Valid analyses
  retain their selected specification. The same check applies to custom
  dimensionality requests and rejects stale results after an unsuccessful run.
  Project files can still save settings without a computed result. These are
  export conformance checks.
- `tests/testthat/test-recovery-estimands.R` checks response-scale agreement
  before Rasch, EFRM and MFRM recovery comparisons. Its sparse PCM example
  verifies refusal after category compression. BTL-EFRM checks align a valid
  change of reference origin and refuse a reference with a generating unit
  other than one. Generators and estimators are unchanged; these checks
  prevent comparison of different estimands. Earlier recovery summaries with
  compressed categories or a changed reference must be reviewed separately.
  Older MFRM simulations without a recorded category count retain an explicit
  note that response-scale equivalence cannot be verified.
- `tests/testthat/test-varying-person-units.R` checks mixed-category EFRM
  and externally weighted person scoring against their score equations and
  SE formulas after changes of unit. It includes missing responses, extreme
  scores, unequal and zero external weights, and equal-weight reduction to
  EFRM scoring. `test-weighted-project-migration.R` checks authenticated
  migration of both earlier weighted solvers and rejects changed results.
  These numerical conformance checks do not change the estimators. Historical
  calibration studies using ordinary units remain applicable; person scoring
  from unusually rescaled calibrations should be rerun.
- `tests/testthat/test-distractor-report-refusal.R` checks that a valid keyed
  calibration with repeated IDs can export CSVs, figures and reports while
  recording the intentional refusal of distractor analysis. Unexpected
  calculation failures still stop the export. This changes reporting, not
  the repeated-person inference policy.
  `test-report-probability-format.R` checks the actual document formatter's
  handling of very small adjusted BTL-EFRM unit probabilities.
- `tests/testthat/test-btl-information-tails.R` checks dichotomous CJ
  information against category probabilities in both pair orientations,
  its nonzero tail values, targeting-plot refusal without graphics changes,
  and a polytomous case whose information does not peak at a zero gap.
  These are information and display conformance checks; CJ estimation and
  its historical calibration studies are unchanged.
- `tests/testthat/test-mfrm-factor-collapse.R` compares named-column and
  data-frame person factors for long and wide MFRM data, including partial
  missingness, row reordering, conflicting values and DIF results.
  `test-person-scoring-units.R` checks equivalent common-unit WLE and MLE
  locations and SEs, extreme-score extrapolation and the WLE score equation
  after changes of unit. These are conformance checks, not new calibration
  studies. Historical studies with fully observed person factors and ordinary
  measurement units do not need rerunning for these corrections.
- `tests/testthat/test-scored-missing-codes.R` compares keyed scoring with
  manually prepared data when raw missing codes overlap generated scores.
  It covers binary and polytomous keys, mixed keyed/numeric data, anchors,
  explanatory models, and item removal, splitting and superitem refits.
  An EFRM refit also checks preservation of renumbered scores. Comparisons
  include response matrices, likelihoods, person estimates and calibration
  covariance. `test-combine-model-request.R` checks model selection for
  superitem refits, including refusal of unsupported explanatory RSM requests.
  These are conformance checks, not new Type I error estimates. Historical
  studies using numeric scores and the default negative missing code are
  unaffected; analyses with overlapping raw missing codes need to be rerun.
- `tests/testthat/test-interval-policy.R` checks that same-data bootstrap
  refits reproduce item chi-squares under both automatic and explicit interval
  rules. Its missing-response design has six intervals for six items and two
  for two less-exposed items. It also checks tied-score refits, plots, exports
  and saved-result handling. These are conformance checks, not Type I error
  estimates. Historical Rasch item-fit bootstrap results for missing responses
  with automatic intervals are superseded by the study below; their
  original files and provenance are retained. The CJ bootstrap is unchanged.
- `studies/item-fit-bootstrap-intervals.R` repeats the 100-dataset linked-booklet
  null with the historical seeds and 600 bootstrap draws per dataset. A
  second 100-dataset null uses eight items, with two answered by 120 of 400
  persons, and 399 draws per dataset. The summary, per-dataset accounting
  and item-level probabilities are stored under the same stem in `results/`.
  Each row identifies the study, harness and package source state. Marginal
  rates are means of dataset-level item proportions, with Monte Carlo SEs
  calculated across datasets. A partially unavailable family is omitted from
  that metric's denominator and counted separately, not treated as no rejection.
  After completion, `Rscript tools/simval/check-item-fit-intervals.R` independently
  reconstructs the rates and Monte Carlo SEs, checks bootstrap resolution,
  and verifies the accounting and source hashes.
- `tests/testthat/test-tailored-fixed-calibration.R` checks the final tailored
  scoring fit's anchor contract and refusal of unsupported structural, item-fit,
  DIF and dimensionality refits. Older fits without refit metadata receive the
  same protection; observed diagnostics and person scoring remain available.
  These are conformance checks, not a new simulation study. The tailored
  item-shift bootstrap and its historical validation results are unchanged.
- `tests/testthat/test-explanatory-formula-offsets.R` checks explicit refusal
  of unsupported formula offsets in item- and threshold-level Rasch and CJ
  models. Ordinary predictors named `offset` remain valid. This is a formula
  conformance check; supported estimators and historical simulations are unchanged.
- `tests/testthat/test-compare-ic.R` checks that generic EFRM likelihood
  differences are withheld despite identical response data, while the matched
  group-unit comparison is retained. CJ comparisons distinguish declared
  response supports but retain same-scale threshold-model comparisons.
  These are conformance checks; estimators and historical simulation results
  are unchanged.
- `tests/testthat/test-explanatory-mc-refits.R` checks keyed-response alignment
  after person subsetting and reordering, and keeps observed answer options
  out of simulated refits. With one blank row among 200 respondents, the
  keyed and numeric-score scree analyses both use 20 of 20 replicates and
  produce identical reference curves and probabilities under the same seed.
  `test-explanatory-item-names.R` checks exact and unambiguous
  whitespace-normalised item matching, including threshold-level metadata
  and item relaxation. These are conformance checks, not coverage studies;
  historical simulation results were not regenerated.
- `tests/testthat/test-explanatory-unused-predictors.R` checks that unused
  ordinal metadata leaves Rasch and CJ estimates, covariance and diagnostic
  probabilities unchanged. Selected ordinal predictors retain adjacent
  contrasts and must have sufficient observed levels. These are conformance
  checks; no estimator or simulation results changed.
- `tests/testthat/test-explanatory-offsets.R` checks additive offsets of
  +/-1e12 in explanatory Rasch and CJ models. Estimates, covariance, likelihoods
  and diagnostic probabilities must remain unchanged for `~ x`. Separate
  checks retain the weighted threshold origin, refuse genuinely constant
  columns and preserve interaction-only formulas. These are conformance
  checks, not additional simulations of inferential calibration.
- `tests/testthat/test-explanatory-units.R` (under the package root) checks
  predictor-unit invariance for LLTM, LPCM and dichotomous and graded CJ.
  Rescalings from 1e-12 to 1e9 must preserve locations, likelihoods, covariance,
  adjusted coefficient probabilities, Kent comparisons and composite
  information criteria, with coefficient uncertainty transformed back to the
  supplied units. The LPCM case includes repeated-person covariance; graded
  CJ includes judge clustering and a position effect. Deliberately unfinished
  fits must still withhold inference. These are numerical conformance checks,
  not new coverage studies; the historical simulations below were not rerun.
- `results/equating-shift-covariance.csv` checks CJ equated-location uncertainty
  against the full linear covariance transformation in 20 pairs of fitted
  calibrations, including a non-common object (maximum difference 1.7e-16).
  Twenty fixed-bank links also reproduce the shared shift variance exactly.
  The script is `studies/equating-shift-covariance.R`. These are numerical
  conformance checks with fitted weights held fixed, not coverage studies.
- `results/dimensionality-matched-sample.csv` checks reliability sample
  handling in 100 datasets: complete null and bifactor data, random missing
  responses, ability-related missingness, and partial-credit data with missing
  responses. Both PSI values agree with a separate calculation on the shared
  complete-response sample to within 1.2e-16. Complete-data results are unchanged.
  This is numerical conformance, not a study of magnitude bias or coverage.
  The reproducible script is `studies/dimensionality-matched-sample.R`.
- `round1/` — the August 2026 battery (rasch 1.14.2): one directory per
  model/diagnostic family (`dichotomous`, `pcm_rsm`, `mfrm`, `efrm`, `btl`,
  `btl-efrm`, `equating`, `dif`, `dimensionality`, `mc`, `tailored`), each
  containing the study scripts with their seeds inline. `dm_*.R` are the
  follow-up studies that diagnosed and verified the `dependence_magnitude`
  standard-error correction (null calibration, covariance diagnosis at
  1,200 replicates, post-fix confirmation).
- `results/round1-checks.csv` — the structured result of every round-1
  check: area, condition, metric, observed, expected, pass/fail, notes
  (including Monte Carlo error where computed), and provenance columns.
  231 of 234 checks passed. The three failures are retained deliberately:
  two are the `dependence_magnitude` null-calibration failure that led to
  the 1.14.2 standard-error correction (those rows describe PRE-fix
  behaviour; the post-fix confirmations are `round1/dm_post_fix.R` and the
  1,200-replicate diagnosis `round1/dm_worker.R`), and one is a
  mis-designed power scenario in the battery itself (non-uniform DIF
  planted on an item whose discrimination crossover falls outside the
  observed trait range -- undetectable by construction; the redesigned
  check appears as a passing row). Round 1 was executed 2026-08-10/11
  against the 1.14.1 release-candidate R code (commits 433e2a5..326f0b9),
  before the 1.14.2 fix; the `area` column maps to `round1/`
  subdirectories as follows: dichotomous-core -> `dichotomous/`,
  PCM/RSM -> `pcm_rsm/`, MFRM -> `mfrm/`, EFRM -> `efrm/`,
  BTL/Recovery -> `btl/`, BTL-EFRM -> `btl-efrm/`,
  equate_tests -> `equating/`, DIF -> `dif/`,
  residual_correlations/dependence/dimensionality -> `dimensionality/`.
- `results/dependence-magnitude-fix.csv` — the full lifecycle of the one
  round-1 defect as structured rows: pre-fix Type I 7.5% (1,200
  replicates), 5.4% with the corrected formula on the same draws, 4.8%
  post-fix on fresh seeds, 3.5% partial credit; each row names its script
  and the package state it ran against.
- `parse_check.R` — parses every script in the battery and fails loudly on
  any syntax error; run it (and a smoke subset) before trusting the
  reproducibility claim.
- `results/scree-conditional-reference.csv` records the score-conditional
  scree study. Across five supported null designs, familywise rejection over
  ten displayed components was 2.0--5.5%, and 40/1,000 (4.0%) overall. Power
  for a planted second dimension was 98.5--100% in the complete designs and
  71.5% in the booklet design. Every supported dataset and all 100,000 inner
  draws were analysed. A sparse-PCM stress cell analysed 178/200 datasets
  (2 fits refused and 20 reference analyses unavailable); its conditional
  familywise rate was 2.8%, with 8,751/8,900 inner draws used. The study is
  `studies/scree-conditional-reference.R`; its rows carry script hash
  `7c42ef04e0a145c398cfa46a977c0819` and R-tree hash `17a08ef1fa4c`.
- `results/dimensionality-bootstrap.csv` separates a content split fixed in
  advance from the residual-derived split and its parametric-bootstrap
  calibration. Across 100 null datasets, the uncalibrated automatic rule
  rejected 15% of dichotomous samples and 41% of PCM samples; bootstrap
  rejection was 5% and 9% (MCSE 2.2 and 2.9 points). All three procedures
  detected the planted PCM second dimension. All 300 outer datasets were
  analysed; 29,697 of 29,700 inner refits converged and none otherwise
  failed. These cells use one response row per person; the public procedure
  now refuses repeated identifiers because its person comparisons and
  score-conditional null generator do not model within-person dependence.
  The study and R-tree hashes are
  `abae964531ee3fc544ecc0ba226c89d6` and `17a08ef1fa4c`.
- `results/repeated-id-item-fit-guard.csv` checks the inferential boundary for
  stacked response data. Exact row duplication leaves threshold estimates and
  their person-clustered covariance unchanged. The ordinary item-fit
  probabilities are withheld, descriptive statistics remain available, and
  the independent-row item-fit bootstrap is refused in every replicate.
- `results/btl-dimensionality-reference.csv` checks the pooled observed-minus-
  expected point residual and its finite simulated upper-tail decision.
  Binary, ordered-response and fitted-position null cells rejected 2.5--3.5%
  with 20 reference draws and 2.5--2.8% with 200 (1,000 datasets per cell;
  MCSE 0.49--0.58 percentage points). These designs were conservative.
  At 200 draws, power for a directed-cycle departure of 0.75 logits was
  5.4%, 97.6% and 7.8%, respectively; at 1.50 logits it was 82.6%, 100%
  and 79.0% (500 datasets per cell). Thus the binary diagnostic had little
  power for the weaker departure in this design. The BTL-EFRM conditional
  reference rejected 2.8% and 3.6% (250 datasets per cell; MCSE 1.04 and
  1.18 points). This arm assumes independent outcomes within its fitted
  design; it does not validate general within-judge dependence. All 12,500
  analyses completed, with no refusals, non-convergence, errors or withheld
  references. Data draws were paired across reference sizes, so these are
  not 12,500 independent datasets. The study is
  `studies/btl-dimensionality-reference.R`; its rows carry script hash
  `4ca05b42c9b81f77fa975d489e719ee2` and R-tree hash `bceec205f5ba`.
  `results/btl-dimensionality-reference-pre-pooled-expected.csv` preserves
  the superseded results; they are not evidence for the current statistic.

## Known limitations surfaced by the battery

- Small-sample undercoverage for extreme items: with 20-100 persons per
  form, item-threshold coverage runs 0.88-0.94 (SE ratios up to ~1.5) for
  items placed 3 logits from the person mean. The fully crossed reference
  arm shows the same undercoverage as the linked-booklet arm at every
  sample size, and both are nominal by 600 persons -- this is
  small-sample behaviour of the pairwise-conditional standard errors for
  poorly-informed items, not a cost of the booklet structure, and it is
  compounded at the smallest sizes by conditioning on the 4-20% of fits
  the identification guard refuses (`structural-missingness.csv`,
  weak_links rows).
- Pooled threshold rows in the sparse-category scenarios mix items of
  very different precision; their emp_sd/mean_se exceeds 1 by Jensen's
  inequality even when every item calibrates. The per-item rows carry the
  calibration-relevant figures.

Exact provenance state of the result files: the three studies regenerated
in the final release round (`btl-clustered.csv`, `btl-share-sweep.csv`,
`efrm-fix-sweep.csv`) carry R-tree hash `4a0e2bf8b357`, matching commit
`8de917b` exactly; `lr-smalln-topup.csv` carries the dirty-tree hash of
its execution state (`42945dd80663`): the `lr_test` estimator code it
exercised was unchanged relative to the release tree, and the hash
difference comes from documentation comments and unrelated BTL source
edits in flight at execution time. Files produced before the tree-hash
mechanism existed carry none. Documentation-only edits under `R/` move the tree hash of
later commits without touching any estimator, so a result file's hash
identifies the sources it ran against, not necessarily HEAD.

The accounting fields added to `btl-clustered.csv`,
`coherence-fixes.csv` and `tailored-bootstrap.csv` on 2026-09-03 were
reconstructed from their recorded denominators without rerunning the fitted
models. They distinguish guarded inferences, unavailable individual metrics
and failed fits; the estimates and Monte Carlo results are unchanged. Their
script hashes continue to identify the scripts that produced those results.
The short commit identifier in `misfit-repair.csv` is stored with a `git:`
prefix so CSV type inference cannot mistake its hexadecimal notation for a
number; this likewise leaves the study results unchanged.

Anchored calibration is checked in `studies/anchored-estimation.R`, with six
rows in `results/anchored-estimation.csv`. Mixed-score PCM with individual
threshold and item-location anchors gave maximum absolute item bias 0.022
logits, mean itemwise empirical SD/mean SE 1.003 and coverage 0.946 over 250
replicates. Two disconnected
four-item blocks, identified by one true anchor in each, gave maximum absolute
item bias 0.013, ratio 1.040 and coverage 0.935; one replicate had a weak-item
SE withheld. Anchored BTL gave maximum absolute object bias 0.018, ratio 1.014
and coverage 0.943 over 500 replicates. No fit failed or did not converge, and
the anchor constraints held to numerical precision. These results are
conditional on the supplied anchor
values: they do not include uncertainty from an earlier calibration. The
script and R-tree hashes are `b1fd106131141c5b76e966cad99bd69c` and
`ae784f3ee517`.

The current-estimator studies run on 2026-08-21 carry R-tree hash
`360e3609691b`: `alpha-correction-limits.csv`,
`alpha-npml-coverage.csv`, `mfrm-pooled-dif.csv`,
`tailored-bootstrap-topup.csv`, `btl-equating-clustered.csv`, and
`cross-package-validation.csv`. Their script hashes identify the exact study
files used. `btl-efrm-current.csv` and `btl-efrm-bias-sweep.csv` were rerun on
2026-08-27 after the reconciled-panel refit; both carry R-tree hash
`dc534f567649`. Their current results are described below.

The structural-null, structural-alternative and repeated-measures DIF
bootstrap studies run on 2026-09-01 and 2026-09-02 used the earlier raw-F
marginal reference and its floor on the adjusted result. Their R-tree hash is
`067abd31b102`. They remain in the repository as an execution record, but
their bootstrap rates are not presented as validation of the current
reference-probability calculation. The current conformance and paired
transition results are identified below by their own R-tree and script hashes.

The DIF bootstrap implementation is checked through two independent public
refit studies. `studies/dif-bootstrap-public-conformance.R` covers ordinary
dichotomous, PCM and RSM fits; its 49 generated datasets preserved every
person score and missingness pattern and reproduced the retained F statistics,
reference probabilities and bootstrap adjustments exactly. The structural
extension is in `studies/dif-bootstrap-model-conformance.R`, with results in
`results/dif-bootstrap-model-conformance.csv`. Across 19 replicates for each
of explanatory Rasch, Multiple Ratings, Extended Frames and Comparative
Judgement, all refits were usable. Person scores, item-set subtotals, or the
paired-comparison design were preserved as appropriate, and the independently
orchestrated F and reference-probability matrices agreed exactly with
`dif_bootstrap()`. The result rows record the executed R-tree and script hashes.

`studies/dif-bootstrap-reference-transition.R` compares the earlier raw-F
calculation with the current reference-probability calculation on the same
observed analyses and the same retained conditional draws. Its incomplete
four-occasion design is the setting in which Greenhouse--Geisser references
and sparse class intervals make raw F values least safely comparable. The
result table reports current marginal calibration, current and earlier
familywise rejection, and the frequency with which the two calculations
change a decision.

The earlier null-calibration screen of the structural extensions is in
`studies/dif-bootstrap-structural-null.R`; the two elevated 50-dataset cells
were followed by 150 fresh datasets in
`studies/dif-bootstrap-structural-followup.R`. Combined minimum-p FWER was
9/170 (5.29%) for Extended Frames and 13/200 (6.50%) for ordinary Comparative
Judgement; the corresponding Monte Carlo standard errors were 1.72 and 1.74
percentage points. The clustered-judge stress screen was 1/50. Extended
Frames had a material refusal rate in this small linked design: in the
follow-up, 12 observed links failed, one fit did not converge, and 11 of 137
otherwise fitted analyses retained fewer than 90 of 99 complete-family
refits. Performance is conditional on the 126 analysed datasets. The result
file includes the attempted, usable, non-converged and other-failure counts,
including work done inside a refused bootstrap. These numerical bootstrap
rates describe the pre-reference-probability implementation and are retained
as historical results, not current calibration claims.

The corresponding earlier partial-alternative screen is in
`studies/dif-bootstrap-structural-alternative.R` for explanatory Rasch,
Multiple Ratings, Extended Frames and Comparative Judgement. Each condition
plants one uniform-DIF member at two magnitudes and records affected-member
power separately from familywise error among the remaining invariant members.
This distinction is required because the bootstrap adjustment describes the
fitted global invariant null rather than strong control under every partial
alternative. Its bootstrap rates are likewise historical.

The earlier mixed-design screen is in
`studies/dif-bootstrap-repeated.R`. It crosses a between-person group with
four observations of a within-person occasion factor. Under an exact balanced
Rasch null, primary and bootstrap familywise error were 6.86% and 4.86% over
350 datasets. They were 2.86% and 2.29% in a local-dependence stress condition,
and 4.57% and 4.00% with unequal groups and group-dependent panel loss. With
uniform occasion DIF of 0.70 and 1.20 logits, primary/bootstrap power was
48.8/41.6% and 98.4/96.0%; familywise error among the other items was no more
than 4.0%. Adjusted power for non-uniform occasion DIF was weak: 2.4/2.4% for
a slope increment of 0.80 and 5.6/4.8% for an increment of 1.50. This was not
a bootstrap-specific loss; the primary residual analysis was also insensitive
after adjustment over the complete 42-test family. Every conditional refit was
usable in all seven conditions. Those exact bootstrap rates predate the
reference-probability marginal calculation; the current paired transition
study replaces them as release evidence for the revised algorithm.

The explanatory-model study run on 2026-08-23 is in
`studies/explanatory-models.R`, with results in
`results/explanatory-models.csv`. It covers LLTM and LPCM coefficients,
dichotomous and ordered comparative judgement, judge-clustered covariance,
Kent-adjusted comparisons with free calibrations, and Holm-adjusted fixed
departure diagnostics. Each principal condition used 1,000 replicates; the
diagnostic conditions used 300. The result rows carry script hash
`e72c33409a15c6ecc04bc5fb91f413ca` and R-tree hash `dd5f1154d99f`.

The edge-case extension is in `studies/explanatory-edge-cases.R`, with results
in `results/explanatory-edge-cases.csv`. Four LLTM/LPCM conditions cover 300
to 2,000 persons, mixed maximum scores, and four-category items. Across 1,000
replicates per condition, coefficient bias was at most 0.0049 logits,
empirical SD/mean SE was 0.993--1.026, coverage was 0.942--0.954, and
Kent-adjusted null rejection was 4.3--5.8%. The unscaled probability rejected
98.6--100% and is retained only as `p_naive`. No fit was refused or failed to
converge. Mean calibration R-squared was 0.872--0.994 under the correctly
specified generating models. The rows carry script hash
`42edd1a6a6fda75a19006c01953c2dae` and R-tree hash `954bcba447fe`.

The compiled EFRM linking kernel was checked against the retained R
implementation with the same seed and 30 hybrid bootstrap replicates. Across
set units, their standard errors, origins, edge likelihoods and thresholds,
the largest absolute difference was 1.30e-11; every convergence flag agreed.
The study is `studies/efrm-cpp-parity.R` and its six result rows are in
`results/efrm-cpp-parity.csv`. They carry script hash
`1a6468c78f37e505d3c01870f8a1c693` and R-tree hash `2ff015352a0d`.

Parallel EFRM bootstrap execution was then checked from a fresh isolated
installation. Serial, two-worker and default four-worker hybrid fits used the
same 300 pre-generated resamples; the full-bootstrap comparison used 30. Every
checked estimate and convergence flag was identical. On the executing machine,
two and four workers reduced the hybrid elapsed time from 86.0 seconds to 52.4
and 33.6 seconds; two workers reduced the full-bootstrap time from 8.6 to 6.4
seconds. These timings are contextual, not general performance guarantees.
The study is `studies/efrm-parallel-parity.R`, with results in
`results/efrm-parallel-parity.csv`; the rows carry script hash
`64574bd4400ad20f64b2746b0b1453ed` and R-tree hash `f25f578582a2`.

The BTL-EFRM judge bootstrap was checked separately from a fresh isolated
installation. Serial and default four-worker fits used the same 200 judge
resamples and agreed exactly after excluding the recorded worker count. On the
executing machine, elapsed time fell from 17.34 to 5.47 seconds (3.17 times
faster). The timing is machine-specific. The study is
`studies/btl-efrm-parallel-parity.R`, with results in
`results/btl-efrm-parallel-parity.csv`; the row carries script hash
`db8d37e3a52d75b7a342fa6e43d7f8f8` and R-tree hash `eba5f653300c`.

`results/item-fit-bootstrap.csv` records the 2026-08-30 execution of the
item-fit study. Its performance rows are conditional on analysed replicates,
but that execution combined refusals and other failures in the
`n_nonconv` field and did not retain inner-bootstrap counts. The current
`studies/item-fit-bootstrap.R` distinguishes refusals, non-convergence and
other errors, carries `B`, `B_used`, `B_nonconverged` and `B_errors`, and
treats an entirely unavailable adjusted family as unavailable rather than as
no rejection. The existing CSV retains its original script hash and is not
presented as a rerun of that accounting revision.

The missing-response rerun is in `results/item-fit-bootstrap-intervals.csv`,
with all 200 attempts and 2,300 item rows in the corresponding `-attempts`
and `-items` files. In the 100-dataset booklet condition, item-wise chi-square
rejection was 5.33% (MCSE 0.53 percentage points), and Holm familywise rejection
was 7% for chi-square and 6% for fit residuals (MCSE 2.56 and 2.39 points).
The 100-dataset unequal-exposure condition gave 3.63% item-wise chi-square
rejection (MCSE 0.69 points) and 3% for both Holm families (MCSE 1.71 points).
The exact 95% familywise intervals are 2.86--13.89% and 2.23--12.60% for the
booklets, and 0.62--8.52% for each unequal-exposure family. These studies do
not establish tight error-rate bounds beyond the designs tested.
All 99,900 bootstrap refits were usable, with no outer refusals,
non-convergence, other errors or unavailable metrics. Same-data refits
reproduced the observed item chi-squares exactly in all 200 datasets.
The independent checker verified the complete accounting, per-item-to-summary
calculations, Holm resolution, and study, harness and R-source hashes.

The score-conditional person-fit and fitted-design comparative-judgement
bootstraps are checked in `studies/person-cj-fit-bootstrap.R`. The supported
four-category follow-ups are `studies/person-fit-bootstrap-pcm-topup.R` and
`studies/cj-fit-bootstrap-polytomous-topup.R`; their result files have the
same stems under `results/`. At B = 99, dichotomous person-fit marginal error
was 2.41% and leave-one-out maxT familywise error was 5.0%. In the
well-targeted four-category PCM follow-up they were 3.81% and 6.0%. Each used
100 datasets without an outer refusal or an unusable inner refit. The broad
240-person PCM stress condition in the main study refused 71/100 datasets
because too many bootstrap samples lost sparse categories; its 29 conditional
results are not used to claim calibration. The current comparative-judgement
rerun gave total-test, pair-family, object-family and judge-family error of
3%, 3%, 1% and 2% for dichotomous responses and 3%, 8%, 8% and 3% for
polytomous responses. The two 8% estimates have MCSE 2.7 points. Neither CJ
design was refused. The current model studies carry R-tree hash
`17a08ef1fa4c`; their exact package and study hashes are in the result rows.
Outer errors and inner-bootstrap attempts, usable refits, non-convergence and
other failures are recorded in separate columns; the inner counts reconcile
to the attempted total in every row.

The corrected algorithm is also checked directly in
`studies/maxt-exchangeability.R`. With 20,000
iid global-null experiments, 99 reference rows and a family of ten, the
leave-one-out familywise rates were 4.955% and 5.105% for centred upper- and
two-sided statistics, and 4.895% and 5.110% for their studentised counterparts
(Monte Carlo standard errors 0.153--0.156 percentage points). The structured
rows are in `results/maxt-exchangeability.csv`.

History-dependent paired-comparison fit bootstrap is checked in
`studies/cj-fit-bootstrap-history.R`, with the fresh-seed dichotomous top-up in
`studies/cj-fit-bootstrap-history-dich-topup.R`. Exposure 0.5 and carry-over
0.4 were fitted jointly at 30 judges. Across 1,000 dichotomous datasets, total,
pair-family, object-family and judge-family error was 4.8%, 3.9%, 4.8% and
5.5%; marginal judge error was 4.28%. Across 200 four-category datasets the
corresponding rates were 3.0%, 3.5%, 6.0%, 6.0% and 3.65%. All 238,800 inner
refits were usable, with no outer refusal, non-convergence or other error. The
first 200 dichotomous datasets put judge-family error at 9.0%; the 800 fresh
seeds gave 4.625%, leaving 5.5% over the predeclared 1,000-dataset program.
The object- and judge-family rows use the superseded standardisation; the
total, pair-family and marginal-judge rows are unaffected.
The result files use the study names under `results/`; both carry R-tree hash
`daddaa5c5222` and their script hashes match the files.

Second-round rows carry exact provenance automatically: `sv_row()` stamps
each row with the study script (`options(simval.script = ...)`), the
package git SHA (`+dirty` when the R code differs from HEAD), and the run
date. `round1-checks.csv` predates this mechanism; its rows carry the
executed date, the commit range, and a per-row `script_dir` mapping to the
`round1/` directory that produced them (exact per-script attribution was
not recorded by the round-1 orchestration and is documented at directory
granularity deliberately).
- `studies/` — the second-round studies (custom Wald/contrast tests,
  clustered comparative-judgement inference, tailored bootstrap
  calibration, `lr_test` size and power, equating multiplicity, structural
  missingness, WLE coverage), each writing its result table to `results/`.
  The `custom-wald-tests` and `btl-clustered` tables reflect re-runs after
  the 1.14.2 statistical corrections (their pre-fix rows motivated those
  corrections and are superseded); `tailored-bootstrap` documents the
  bootstrap's strong conservatism at feasible replicate counts, which is
  the procedure's documented design.
- `results/efrm-fix-sweep.csv` and `results/efrm-fix/` — covariance checks
  for the superseded moment-based set link. They established the need to
  propagate within-frame calibration uncertainty, which the current
  semiparametric link retains, but they do not validate its point estimator.
- `results/round2-followups.csv` — the adjudication trail for every
  round-2 suspect: the MFRM q=25 exoneration (600 fixed-truth replicates
  per cell), the btl_efrm origin-test correction (8.5% pre-fix to 5.5%
  post-fix at 400 replicates each), the PCM item-level guard validation,
  and the concentration guard's field behaviour.
- `results/audit3-fixes.csv` — checks the corrections from the third audit.
  Under a 10:90 nuisance-cell imbalance, the repeated-measures DIF follow-up
  held its 5% size (5.25%, 2,000 replicates) when it retained the equal-cell
  estimand; the superseded person-frequency shortcut rejected every dataset
  because it targeted a different contrast. The corresponding mixed
  time-by-group interaction had 5.4% size. BTL-DIF pairwise size was 5.5%
  with eight judges per level and 4.83% with ten. In the five-level design
  with four judges per level, the superseded global degrees of freedom gave
  9.45% pairwise rejection; current public inference is withheld in that
  region. At the binomial spread boundary, the raw estimate fell below the
  bound in 47.1% of datasets, while the current one-sided test rejected 5.2%
  (1,000 replicates). The planted-dependence condition had 100% power.
- `results/audit3-btl-imbalance-topup.csv` — fresh-seed adjudication of the
  mildly imbalanced BTL-DIF cell. With ten raw and 9.31 effective judges per
  level, the current pair-specific reference rejected exactly 5.0% over 2,000
  datasets. A diagnostic minimum-cell reference rejected 3.55% and was not
  adopted for the two-cell comparison. There were no refused or non-converged
  fits.
- `results/btl-dif-multicell-df.csv` — the separate four-cell case. In 500
  balanced 2 by 2 null datasets with 12 judges per cell, the weakest-cell
  reference rejected 3.6% and gave 0.964 interval coverage. The superseded
  pooled-count extension rejected 6.2% and covered 0.938. Multi-cell
  contrasts therefore use the conservative weakest-cell rule while the
  validated two-cell Welch rule is unchanged.
- `results/btl-dif-hc3.csv` — null calibration of the between-judge residual
  test under a fourfold variance ratio. With 8 versus 16 judges, the classical
  equal-variance test rejected 11.72% or 1.78%, depending on which group was
  less precise; HC3 gave 5.57% and 4.04%. With eight judges per group, HC3
  gave 4.13% (10,000 replicates per condition).
- `results/dif-hc3.csv` — ordinary between-person DIF under balanced groups,
  a 1:4 ability imbalance, and unequal observations per person. The adopted
  hybrid (HC3 for uniform terms, residual ANOVA for class-interval
  interactions) gave 4.0%, 6.4%, and 4.6% Holm familywise rejection. Applying
  HC3 to every term gave 22.0% under the ability imbalance and was rejected.
  At a planted 0.6-logit shift, hybrid power was 29.6% against 23.6% for the
  classical analysis (500 replicates per condition).
- `results/dif-hc3-multilevel.csv` — a three-level factor with group sizes in
  the ratio 1:2:3 and different group locations. Hybrid Holm familywise
  rejection was 6.2%, compared with 6.4% for the classical analysis and 20.2%
  when HC3 was also applied to the class-interval interaction. At a planted
  0.6-logit shift, hybrid power was 36.6%, compared with 32.8% classically
  (500 replicates per condition).
- `results/dif-hc3-homoskedastic.csv` — balanced group-by-interval cells with
  independent normal errors and a common variance. Hybrid Holm familywise
  rejection remained between 4.3% and 5.2%. With ten observations per cell,
  HC3 reduced planted-item power from 21.7% to 18.8% for a two-level factor
  and from 23.5% to 21.9% for a three-level factor. The difference was 0.9
  percentage points with 30 observations per cell; fixed-effect power then
  approached its ceiling (5,000 replicates per condition).
- `results/dif-hc3-homoskedastic-local-power.csv` — the same comparison with
  effects reduced as cell size increased, preventing the larger designs from
  reaching the power ceiling. For two levels, the classical power advantage
  declined from 3.10 percentage points at ten observations per cell to 0.84,
  0.42, and 0.16 points at 30, 75, and 150. For three levels it declined from
  1.72 points at ten per cell to 0.54 at 50. This confirms a real but diminishing
  efficiency cost for HC3 when the classical assumptions hold exactly (5,000
  paired replicates per condition).
- `results/dif-conditional-bootstrap.csv` — a conditional Rasch null bootstrap
  that preserves each raw score and refits the model 199 times per dataset.
  Its single-step minimum-p rows already used replicated reference
  probabilities and remain applicable: familywise rejection was 6.33% for
  both procedures in the balanced null and 7.0% versus 6.0% with 1:4 group
  sizes and a 0.8-logit ability difference. Adjusted power was 29.7% versus
  24.3% for a 0.6-logit uniform shift and 5.7% versus 4.3% for a centred 0.7
  slope departure (hybrid versus bootstrap; 300 datasets per condition).
  Its item-wise bootstrap rows used the superseded raw-F calculation and are
  retained only as historical output.
- `results/dif-conditional-bootstrap-extended.csv` — the same comparison for
  four-category PCM and RSM data, a three-level group, and two correlated
  person factors. For a response vector \(x\) with raw score \(r\), the
  polytomous sampler draws from
  \(P(X=x\mid r) \propto \exp\{-\sum_i\sum_{k=1}^{x_i}\tau_{ik}\}\), and every
  draw is checked against the conditioned score before the model is refitted.
  Global-null rejection was broadly consistent with 5% in all designs. The
  bootstrap was usually a little more conservative and a little less
  powerful than the hybrid analysis. It did not remove artificial flags on
  invariant items when another item truly had DIF (100 datasets and 99
  bootstrap refits per condition). These statements concern its unchanged
  minimum-p rows; its marginal bootstrap rows used raw F and are historical.
- `results/dif-conditional-bootstrap-confirm.csv` — fresh-seed confirmation
  with 200 datasets and 199 bootstrap refits. Under the imbalanced global
  null, Holm familywise rejection was 3.0% versus 2.5% for the PCM and 4.0%
  versus 2.5% for the RSM (hybrid versus bootstrap). With a 1.4 slope
  departure on one item, false flags among the other five items occurred in
  19.0% versus 13.5% of PCM datasets and 14.5% versus 11.0% of RSM datasets.
  The bootstrap attenuates score contamination but does not solve it. These
  adjusted rows used the same minimum-reference-p calculation as the current
  public result. The file's marginal rows use the earlier raw-F comparison
  and are not current validation evidence.
- `results/dif-bootstrap-public-conformance.csv` — direct conformance of the
  exported `dif_bootstrap()` result to a separately orchestrated refit loop.
  The check covers a correlated multifactor design, a mixed repeated-measures
  design, PCM and RSM. It verifies the complete F matrix, marginal empirical
  reference probabilities, minimum-p familywise probabilities, raw-score
  preservation and missingness preservation. All discrepancies and failure
  counts were zero over 49 conditional refits per design. This connects the
  public implementation to the operating-characteristic studies above; it is
  not a further estimate of size or power.
- `results/dif-bootstrap-structural-null.csv` and
  `results/dif-bootstrap-structural-followup.csv` — historical global-null
  calibration from the pre-reference-probability implementation
  for explanatory Rasch, Multiple Ratings, Extended Frames and Comparative
  Judgement, including a judge-heterogeneity stress condition and a separate
  top-up of the two initially elevated cells. The follow-up separates losses
  at the observed fit, observed DIF analysis and bootstrap stages.
- `results/dif-bootstrap-structural-alternative.csv` — affected-member power
  and invariant-member familywise error under two uniform-DIF magnitudes in
  each supported structural model family (100 attempted datasets per cell,
  B = 99). This is also a historical pre-change screen. For effects of 0.5
  and 0.9 logits, primary/bootstrap power was
  5/8% and 51/49% for explanatory Rasch, 23/23% and 78/74% for Multiple
  Ratings, and 7/8% and 15/17% for Comparative Judgement. Corresponding
  invariant-member FWER was 3/2% and 8/5%, 10/9% and 23/20%, and 5/6% and
  7/10%. Extended Frames analysed 83 and 80 datasets; conditional power was
  6.0/8.4% and 31.3/25.0%, with invariant-item FWER 0/0% and 6.3/3.8%.
  Observed-fit/DIF refusals accounted for 9 and 11 datasets, and bootstrap
  failures for 8 and 9. The study confirms that the fitted-global-null
  bootstrap is a sensitivity analysis: it can reduce artificial flags, but
  does not provide strong control after one member departs.
- `results/dif-bootstrap-repeated.csv` — familywise error and affected-item
  power from the earlier raw-F-floor implementation for a four-occasion mixed
  DIF design. It includes an exact balanced
  null, a local-dependence stress condition, group-dependent panel loss, and
  uniform and non-uniform alternatives at two magnitudes. Each row records
  Greenhouse--Geisser epsilon and complete inner-refit accounting.
- `results/dif-bootstrap-reference-transition.csv` — current paired
  transition check for an incomplete four-occasion null design. Both
  calculations use each observed dataset's same 99 conditional refits. It
  reports marginal and familywise rejection under the reference-probability
  calculation, the corresponding historical calculations, and changed
  decisions. This is the current release evidence for the algorithm change.
- `results/dif-score-purification.csv` — an initial comparison of leave-one-out,
  fixed-anchor, and iterative matching scores. Leave-one-out testing was
  rejected because null familywise error reached 28--34%. Re-estimating the
  full analysis from a five-item anchor scale also lost too much uniform-DIF
  power. This screen motivated the staged comparison below (100 datasets per
  condition).
- `results/dif-score-purification-refined.csv` — 500-replicate comparison of
  anchor-based class intervals, full anchor recalibration, strongest-item
  iteration, and the public split-and-refit workflow. The unmodified hybrid
  gave global-null Holm familywise error of 4.0% for the PCM, 3.8% for the RSM,
  and 6.0% with two correlated person factors. Preselecting a five-item anchor
  scale was liberal (7.2--10.6%) and is not a valid default.

  For a uniform 0.6-logit shift, `resolve_dif()` split the planted item in
  84.4% of PCM and 86.6% of RSM datasets, against initial detection of 85.0%
  and 87.4%. It split an invariant item in 0.6% and 1.8%, and the final
  invariant-item familywise rates were 4.0% and 5.6%. The existing
  split-and-refit procedure therefore supplies an effective purification step
  for uniform DIF.

  For a centred 1.4 slope departure, correct-term power was 89.0% for the PCM
  and 86.8% for the RSM, while familywise flags among invariant items rose to
  14.8% and 16.8%. A strongest-item, one-at-a-time procedure selected the
  planted item first in 97.6% and 96.0% of datasets; it selected an invariant
  item first in 0.2% and 1.0%, and ever excluded one in 3.4% and 5.2%. After
  anchor recalibration, remaining false flags occurred in 1.6% and 2.0%.
  Retesting the selected item on the short anchor scale needlessly reduced
  power. With two correlated person factors, non-target factor error remained
  controlled (1.6%) but correct non-uniform power was only 12.0%. This is a
  power limit, not a calibration defect. No method from this study has been
  installed as an automatic non-uniform-DIF remedy.
- `results/item-fit-hc3.csv` — sensitivity study for class-interval item fit.
  HC3 was rejected: item-wise null rejection ranged from 21.9% to 48.3% over
  8--30 items, against 5.6--17.0% for the conventional ANOVA. The
  conventional ANOVA remained approximate (Holm familywise rejection 7.5%
  at 30 items and 11.0--31.5% at 8--15 items). The item-trait test was
  calibrated from ten items onward (4.0--7.0% familywise) but not with eight
  items (12.0--17.0%). These results support the short-test qualification in
  `?rasch`, not an HC3 replacement (200 replicates per condition).
- `results/item-fit-interval-count.csv` — reducing the requested number of
  class intervals did not repair short-test calibration and generally reduced
  power. The interval-count change was therefore rejected (100 replicates per
  condition).
- `results/btl-cluster-jackknife.csv` — CR1 versus delete-one-judge covariance
  for the core BTL fit. CR1 Type I was 5.4% with ten balanced judges and 4.2%
  with one of twenty judges carrying 20% of the work; jackknife rates were
  5.6% and 4.0%. The concentrated design below the public effective-judge
  guard gave 6.4% for CR1 and 5.4% for the jackknife, both within Monte Carlo
  uncertainty of 5%. CR1 therefore remains the default (500 replicates per
  condition).
- `results/btl-equating-clustered.csv` — common-object drift under two
  independent 12-judge panels. Welch--Satterthwaite probabilities with Holm
  adjustment gave 4.4% familywise rejection over 1,000 null replicates; the
  superseded normal reference gave 7.4% on the same fitted samples.
- `results/pcml-cluster-cr1-screen.csv` — seed-paired comparison of the
  repeated-person PCML covariance as fitted with the same covariance
  multiplied by G/(G-1), retaining the Welch--Satterthwaite reference.
  Across 10, 12, 20 and 50 person clusters, the current dichotomous Holm
  familywise rates were 4.6%, 4.6%, 5.8% and 5.0%; the scaled rates were
  3.8%, 2.0%, 5.2% and 4.6%. A five-item PCM cell at 20 clusters gave 6.3%
  and 5.5% (Monte Carlo SE about one percentage point). The correction
  slightly improved pooled unadjusted rejection in some cells, but did not
  improve the declared Holm family consistently and was markedly conservative
  at 12 clusters. It was therefore not added to ordinary PCML (500 replicates
  per condition; four dichotomous datasets were withheld by the support guard,
  and five PCM datasets had an unavailable item contrast).
- `results/btl-btm-agreement.csv` — `btl()` against `sirt::btm` on shared
  dichotomous comparisons with eps = 0 and fixed home advantage: identical
  likelihood, mean max-difference 1.1e-15 and worst 3.1e-15 over 23 clean
  replicates. Replicates with a fully extreme object are characterised
  separately: btl sets the unidentified object aside and reports an
  extrapolated boundary location (score half a point inside the boundary,
  SE withheld), whereas btm keeps it, diverging at eps = 0 and shrinking
  it under its default eps = 0.3. The two policies agreed in direction in
  every case, differing by 1.6 logits on average in the flat region of
  the likelihood — a policy difference, not an estimation difference.
- `results/equating-holm-refresh.csv` — the null familywise cells of the
  equating study re-run after `equate_tests()` moved from BH to Holm
  adjustment, mirroring the original design: 4.8-5.0% at 3, 5, and 10
  anchors (2,000 replicates each, no refusals). The null rows of the
  round-2 equating-multiplicity table are BH-era and superseded by this
  table; its drift-power rows also predate the switch and read at most
  slightly high for the current function.
- `results/explanatory-r2-adjusted.csv` — sampling behaviour of the
  calibration R-squared reported by `explanatory_test()`. Under
  uninformative designs the raw coefficient averaged 0.170 (12 items) and
  0.085 (24 items) while the rank-based adjusted coefficient centred at
  -0.015 and -0.002; a true 12-item design gave 0.948 raw and 0.936
adjusted (300 replicates per condition).

- `results/person-external-weights.csv` — externally weighted person
  estimation with fixed generating calibrations. The 18 conditions cover
  equal, moderate, strong and zero weights; dichotomous and partial credit
  items; three person locations; and differing model units. Across 5,000
  persons per condition, absolute bias was at most 0.016 logits,
  empirical SD/mean reported SE was 0.940--1.004, and 95% coverage was
  0.941--0.978. No estimate was refused. The study is
  `studies/person-external-weights.R`.
- `results/comparison-validation.csv` — the release gate for the app's
  automatic model-comparison cards: every comparison surface validated
  under its null model and at least two departure magnitudes. Citation
  rows point surfaces already validated elsewhere (`lr_test`, the EFRM,
  BTL-EFRM, and MFRM omnibus tests) at their provenance CSVs. New cells
  (400 replicates per selection condition): CL-AIC null false selection
  5.2% for PCM-vs-RSM and 4.5% for free-vs-two-component thresholds
  (matching the theoretical multi-parameter AIC rates), rising to ~17%
  for the one-parameter comparative-judgement threshold comparison
  (the familiar P(chi-sq_1 > 2) = 15.7% AIC property); detection 95-100%
  at the stronger departures; CL-BIC selects the smaller model almost
  always under nulls with little power against mild departures. Effect
  tests: position/exposure nulls 5.8%/5.9% (800 replicates); carry-over
  8.3% at 14 judges falling to 5.3% at 30 judges, which is why its
  probability is now withheld below 30 judges; power at 0.6 logits was
  62/39/77%.
  The EFRM log set-unit bias is +0.0036 against TAM's +0.0008 for
  dichotomous data and +0.0035 against +0.0020 for polytomous data. In a
  crossed two-set by two-group design, EFRM log-alpha bias is +0.0141 and
  log-phi bias +0.0004. The person-group unit is +0.003 against a per-group
  `lme4::glmer` coefficient-slope anchor at −0.003. `rasch_mfrm` and
  `tam.mml.mfr` item and rater estimates correlate above 0.9999.
- `harness.R` — shared reporting helpers for the second-round studies: one
  row per scenario with bias, empirical SD, mean reported SE, SE ratio,
  95% coverage, Type I / familywise error or power, refusal and
  convergence rates, and the Monte Carlo standard error of each rate.
- `results/cross-package-validation.csv` — estimates checked against
  independent implementations on shared datasets (25 replicates per
  cell). Identical-likelihood comparators agree to solver precision:
  `btl` dichotomous vs `BradleyTerry2::BTm` to 2.6e-9 logits, polytomous
  comparative judgement vs constraint-matched `VGAM::vglm(acat)` to
  ~1e-7 with log-likelihoods equal to 5e-12. Estimator-variant
  comparators agree to the expected order: `sirt::rasch.pairwise` mean
  max-difference 0.02 logits, `eRm` full conditional ML 0.07
  (dichotomous) and 0.15 (PCM thresholds, worst 0.34 at sparse
  extremes), with equal truth RMSE. Current EFRM comparisons are summarised
  above. The `btl_efrm` panel-unit ratio is unbiased within Monte Carlo error
  and tracks a
  per-panel intercept-free adjacent-category anchor.
- `results/alpha-correction.csv`, `alpha-correction-designs.csv`,
  `alpha-n-sweep.csv` and `studies/alpha-bootstrap-pointest.R` evaluate the
  superseded moment-based set link. They are retained to document why it was
  replaced and must not be cited as evidence for the current estimator. The
  scratch script that produced `alpha-n-sweep.csv` was not retained, so that
  file is a historical, non-reproducible record; its recorded path and hash
  identify the missing script rather than an executable repository study.
- `results/alpha-correction-limits.csv` — point-estimator stress tests for the
  current finite-grid semiparametric link. The study covers targeting, unit
  ratios, short sets, small samples, heavy-tailed and bimodal populations,
  missingness, guessing and within-set discrimination departures. Sets with
  fewer than four score steps are refused by design. Absolute bias is at most
  0.022 under the model. At 80 persons, 11% of datasets are refused and 2%
  do not converge; the 41- and 101-point grid results are effectively equal.
- `results/alpha-npml-coverage.csv` — sampling calibration of the current
  estimator under normal, wide bimodal and deliberately different group
  distributions. Raw marginal hybrid set-unit Type I is 4.0--5%, SE ratios
  are 0.97--1.05, and coverage is 0.927--0.960. Common-scale item SE ratios are
  0.97--1.04. The complete bootstrap is mildly conservative under the null
  (2.5% rejection, 0.975 coverage) and calibrated under the planted ratio
  (SE ratio 0.99, coverage 0.938).
- `results/efrm-unit-multiplicity-supported.csv` — complete-family null
  calibration for the EFRM decision policy, with 200 persons per group and
  six items per set. Of 500 attempted fits, 489 were analysed and 11 were
  refused. At least one raw omnibus probability was below 0.05 in 10.4% of
  analysed fits. Holm familywise rejection was 3.5% for the omnibus family
  and 1.6% for the separate individual-contrast family. The script is
  `studies/efrm-unit-multiplicity-supported.R`; the rows carry script hash
  `6341acc22941fb7041fb08ec8f8bef0c` and R-tree hash `7356d398664f`.
- `results/frame-unit-multiplicity.csv` — the same complete-family check for
  BTL-EFRM, and an EFRM boundary design. The BTL-EFRM rows predate the
  reconciled-panel refit and are superseded; they are retained as provenance,
  not as evidence for the current estimator. The four-item,
  100-person-per-group EFRM boundary attempted 1,000 fits: 493
  were refused, 24 did not converge, and 483 were analysed. Conditional on
  analysis, its Holm rates were 5.2% and 2.5%. The script is
  `studies/frame-unit-multiplicity.R`; the rows carry script hash
  `45abbc7d7938f9fc2416d2dbbf356eeb` and the same R-tree hash.
- `results/audit-adjusted-dependence.csv` — null familywise error and power
  for the crossed-EFRM and BTL dependence decisions, the finite-object
  correlation check, and the affected dimensionality-power cell. Crossed
  EFRM rejected 6.7% in the first 1,000 null fits; the independent top-up in
  `results/crossed-efrm-factorial-topup.csv` gave 5.55% over 2,000 fresh fits,
  with marginal rates 5.25--5.75% and no refusals or non-convergence. The
  pooled familywise rate is 5.93% over 3,000 fits. BTL dependence familywise
  rejection was 5.9% over 1,000 fits. The simulator reproduced requested
  finite-object correlations to 3.9e-16. Those retained rows carry script
  hash `3f1af2fcb6098b5996bda3e0987763f6`; the independent top-up carries
  `70e9130564c43819760d330e925f53b1`. Both identify R-tree hash
  `8a6cc825b06a`. The dimensionality cell was rerun using the pooled-expected
  residual and explicit conditional-independence opt-in: power was 83%
  over 100 analyses (MCSE 3.76 percentage points), with no refusals,
  non-convergence, errors or withheld references. This replaces the earlier
  87% claim. The new row from `studies/audit-adjusted-dependence.R` carries
  script hash `d8c6861e23c1e23c0970a08fb8f6a43f` and R-tree hash
  `bceec205f5ba`; the unrelated rows were preserved, not rerun.
- `results/btl-efrm-current.csv` — current judge- and independent-outcome
  bootstrap calibration for BTL-EFRM, rerun after the reconciled-panel refit.
  Over 300 null fits, raw marginal judge-bootstrap Type I was 3.7% for panel
  units, 7.3% for set units and 5.0% for origins. The corresponding
  independent-outcome rates were 4.3%, 6.7% and 3.7%. Set-unit coverage was
  0.890 with the judge bootstrap and 0.933 with the independent-outcome
  bootstrap; coverage for the other units was 0.923--0.963. This design has
  six judges per panel and lies in the documented caution band.
- `results/btl-efrm-multiplicity-current.csv` — 1,000 current-estimator null
  fits in that caution-band design. Raw marginal Type I was 3.3% for panel
  units, 6.7% for set units and 6.0% for origins. Holm familywise error was
  3.9% across the three omnibus decisions and 3.0% across the individual
  follow-ups. Set-unit coverage was 0.900; no fit was refused or failed to
  converge. The script is `studies/btl-efrm-multiplicity-current.R`; its rows
  carry script hash `a3635aff9d2c9ed79f437686651f3dc0` and R-tree hash
  `8b2530afb990`.
- `results/btl-efrm-supported-topup.csv` — the supported-design top-up with
  12 judges per panel and the public default of 200 resamples. Across 500
  null fits, raw set-unit Type I was 4.6%, the empirical-SD/mean-SE ratio was
  0.992 and coverage was 0.934. The Holm-adjusted set-unit rate was 2.0% in
  the omnibus family and 1.6% in the follow-up family. There were no refusals
  or non-convergences. The script is
  `studies/btl-efrm-supported-topup.R`; its rows carry script hash
  `1f58bf65d9a4d192773c296596a78135` and the same R-tree hash.
- `results/btl-efrm-bias-sweep.csv` — finite-sample attenuation from the
  staged BTL-EFRM set link, rerun after the reconciled-panel refit. Log
  set-unit bias decreases from −0.106 at 10 repetitions per pair to −0.041
  at 20, −0.016 at 50 and −0.006 at 100 (500 datasets per cell). The
  caution-band bootstrap study above shows that this finite-sample
  attenuation can also reduce Wald coverage.
- `results/coherence-fixes.csv` — direct checks of the repaired BTL-EFRM
  fit and the multifactor DIF estimands. Across three information levels,
  reported BTL-EFRM log likelihoods agreed with likelihoods reconstructed
  from every stored fitted probability to numerical precision. Mean object
  RMSE fell from 0.401 to 0.232 and 0.123 as repetitions increased. In a
  correlated 3:1:1:3 two-factor design, ordinary DIF marginal magnitudes had
  biases −0.021 and +0.017 logits, and paired-comparison DIF biases +0.026
  and +0.032; 95% coverage was 0.92--0.96. The study is
  `studies/coherence-fixes.R`; its 29 rows carry script hash
  `25e3f34eddee8bb9d71f77d3c79e8a31` and R-tree hash `dc534f567649`.
- `results/mfrm-pooled-dif.csv` — null calibration after putting uniform and
  non-uniform item tests in the one multiplicity family used for decisions.
  Familywise rejection is 4.7% with balanced raters and 3.8% when the second
  group has only two raters (1,000 attempted datasets per cell).
- `results/tailored-bootstrap-topup.csv` — full automatic anchor selection
  repeated inside each of 399 person-bootstrap draws. Clean-item familywise
  error is 2.5% at guessing 0.15 and 0% at 0.30. At least one of two planted
  hard items is detected in 17.5% and 26.3% of datasets respectively, showing
  that calibration is conservative and power is limited for this design.
- `results/tailored-bootstrap.csv` — current-tree rerun with complete outer
  and inner accounting. Under no guessing, 0/54 analysed datasets had a
  Holm-adjusted item flag (exact 95% interval 0--6.6%); one of 55 datasets was
  refused by the 90% usable-refit rule. In all, 21,802/21,945 null refits were
  usable, with no non-convergence and 143 other failures. All 44 datasets and
  all 19,980 inner refits in the power grid were usable. Its per-cell counts
  are too small for precise power estimates; the 80-replicate top-up above is
  used for those. The current result carries script hash
  `55f5ce4436029ac4cbfb6f853336b04f` and R-tree hash `3ae2582c9cf7`.
- `results/cross-package-diagnostics.csv` — the diagnostics checked
  against independent implementations, at the level each comparison
  supports. Formula parity: alpha vs `psych::alpha` to 4e-15. Value
  level: infit/outfit vs `eRm::itemfit` r 0.97/0.99 once the mean-square
  convention is aligned (our E[z^2] divisor sits 8-10% above eRm's /n),
  with 25/25 same-direction flags on a planted over-discriminating item;
  PSI vs `eRm::SepRel` to 0.004 under matched conventions; q3/q3_star vs
  TAM Q3/aQ3 r 0.92 native, 0.97 at shared person estimates (top-1
  localisation of a planted dependent pair is weak for both: 2/25 vs
  8/25 — the flag rule, not top rank, is the operative criterion).
  Rank level: person fit vs PerFit lzstar/U3 |Spearman| 0.97-0.98, 87%
  worst-decile overlap, equal planted-careless detection (0.76 vs 0.77).
  Decision level: uniform DIF detection 84/88/80% (rasch/Waldtest/MH,
  BH-aligned) with comparable false-positive rates; dimensionality —
  ours gave no false flags in the sampled null datasets and 67% power at a balanced
  planted second dimension where DETECT flags 100% (also with clean
  nulls) and the quasi-exact Tmd/T11 tests pair full power with 13-20%
  null false-flag rates. Citation rows name the diagnostics with no
  external parallel (dependence_magnitude, spread_test, the tailored
  bootstrap, the comparative judgement family, equate_tests), which
  remain simulation-validated only.
- `results/humphry-item-side.csv` — is the item-side variance-ratio
  argument (Humphry 2005, eq. 2.27) biased? Not detectably. Over 15 design
  cells (8 to 40 items, 200 to 5,000 persons per frame, planted unit ratio
  1.30, 30 replicates each) the raw SD ratio's pooled log bias is +0.00008
  with a standard error of 0.0025, and no cell reaches two standard
  errors. The attenuation mechanism is real but arithmetic puts it out of
  reach: item standard errors contribute one to two per cent of the
  observed item variance, and the two frames' error variances differ in
  the offsetting direction. Correcting is pointless here — subtracting
  `mean(se^2)`, or the covariance-correct `tr(V)/(K-1)` that accounts for
  the sum-zero constraint, moves the estimate by about 0.004. Contrast the
  person side, where the error share exceeds half and the naive
  construction biased the ratio by five per cent.
- `results/common-item-channel.csv` — when two frames share items, three
  channels estimate the unit ratio from the same data: the bilinear
  conditional fit `rasch_efrm` uses for phi, the SD ratio of two
  within-frame calibrations, and their regression slope. Conditional ML
  and the SD ratio are indistinguishable — both unbiased, with sampling
  standard deviations agreeing to three decimals across all six cells. The
  slope channel is attenuated as errors-in-variables predicts (−0.013 to
  −0.034, worsening as the item grid densifies and the true item variance
  falls), which is the control confirming the comparison measures what it
  claims. Implication for design: a bridge design putting some items in
  both set contexts would estimate alpha on an unbiased channel needing no
  correction and no new estimator, at a log-ratio SD of 0.087 against
  about 0.105 for the person-side link at equal budget, improving to 0.024
  at 40 items. Untested caveat: the study gave each frame independent
  persons; a bridge design shares persons across contexts, so the two
  calibrations' errors correlate. That caveat is settled in
  `bridge-item-design.csv`, and the bridge idea does not survive it.
- `results/bridge-item-design.csv` — a bridge design for item-set units
  requires the same person to answer an item in both set contexts, which
  means re-administration. Sharing persons is harmless on its own: with
  responses conditionally independent given theta the SD ratio stays
  unbiased and its sampling standard deviation is no worse than under
  independent persons (0.067 against 0.087 at 8 items and 250 persons),
  since correlated errors partly cancel in a ratio. Carry-over between the
  two administrations is what breaks it: half a logit inflates the unit
  ratio by about 12 per cent and a full logit by 17, in every cell, with no
  decay in sample size or item count — a person repeating their first
  answer makes the second administration look more consistent, which the
  calibration reads as a larger unit. Since item sets are defined by item
  properties, an item belongs to one set, so a bridge means literal
  re-administration and conditional independence is not credible. The
  person-side link is therefore the practical route to item-set units. The
  dependence would at least be detectable: paired
  administrations show a large Q3 in `residual_correlations`, though
  detection only tells you to abandon the bridge.
- `results/humphry-pgd-replication.csv` — a simulation of the
  person-group-discrimination design in Humphry (2005, ch. 4): 12 common
  items linking Year 5 and Year 7, calibrated separately, the unit ratio
  read off the ratio of the common items' location standard deviations.
  Planted 1.306 (his phi_5 = 0.875, phi_7 = 1.143); recovered 1.303 to
  1.311 in every cell, bias within 0.004 of zero at 200 replicates. The
  ability gap between year levels, swept from 0 to 1 logit, changes
  nothing, so differential targeting is not a threat to the estimator.
  Sampling SD is 0.04 at his N of 980 and 0.02 at 5,000, which reconciles
  his own two figures: the full-population 1.22 and the sample 1.30 differ
  by 0.064 in the log, about 1.6 SD of their difference — sampling
  variation, not a discrepancy needing explanation.
- `results/humphry-pgd-misfit.csv` — the exposure that design does carry.
  Item-level departures on the common items enter the dispersion as if
  they were unit differences: two of twelve items with a one-logit uniform
  DIF shift move the recovered ratio from 1.306 to 1.403, and four items
  with discrimination multiplied by 1.5 take it to 1.442. The direction is
  not fixed — four items with moderate alternating-sign DIF pull it down
  to 1.289, their shifts partly cancelling within the spread. Robust
  dispersion measures do not rescue it: the median absolute deviation is
  worse than the standard deviation in five of six cells (1.712 against
  1.442 under differential discrimination), because a MAD over 12 items is
  itself a poor dispersion estimator and that noise outweighs its
  contamination resistance. The discipline the method needs is screening
  the common items for DIF and misfit before computing the ratio, which is
  what his own RMSD 0.24 against RMSE 0.12 diagnostic was detecting.
- `results/channel-head-to-head.csv` — a historical comparison of Humphry's
  item-side estimator and the superseded score-moment person-side link on
  identical simulated data (same persons, same items, two frames differing
  only in unit, conditional independence given theta). Both are unbiased
  everywhere, so the comparison is of efficiency: the item-side
  channel is 1.96 times more precise at 8 items and 980 persons, 1.46
  times at 12 items, and 1.15 to 1.33 times by 20 to 40 items. Adding
  persons sharpens item locations and helps the item-side channel; adding
  items sharpens person estimates and helps the score-moment channel. This
  comparison does not validate the current semiparametric link.
- `results/misfit-both-channels.csv` — a historical comparison of ordinary
  item misfit in the item-side and superseded score-moment channels. Four of
  twelve items with discrimination doubled or halved attenuate the recovered
  ratio by 2.3 per cent (item side) or 1.0 per cent (score-moment side), and
  scattering every item's discrimination log-normally costs about 1 per cent.
  It does not validate the current semiparametric link; the corresponding
  current-estimator departures are in `alpha-correction-limits.csv`.
- `results/humphry-isd-replication.csv` — historical comparison of the
  earlier item-set estimators on Humphry's item-set discrimination
  study replicated on its own design (4 sets of 10 items spanning -4 to 4,
  N = 1000, planted ISDs 0.604/0.906/1.209/1.511). ISD is estimated
  person-side in the thesis -- "a matrix of log ratios of standard
  deviations for common persons across the sets", corrected by equation
  2.29, var(WLE) minus the mean squared standard error -- which is the
  construction this package replaced. On his design the uncorrected ratio
  is attenuated to 2.089 against a planted 2.502 end-to-end; equation 2.29
  overshoots to 2.698 (+7.8 per cent); the superseded score-moment
  correction lands at 2.483 (-0.8 per cent). Note that the product
  constraint fixes the mean of log alpha, so mean bias is zero by
  construction for every estimator and only the SPREAD can be wrong --
  which is why the end-to-end ratio is the discriminating statistic.
  Caveat: the per-set SDs here (0.74/1.15/1.57/1.99) sit about 12 per cent
  below his Table 3.10 (0.84/1.30/1.74/2.14). The thesis is internally
  inconsistent about the person spread (his expected SDs imply 1.51, the
  text reports a generated 1.76) and he used RUMM2020's WLEs, so a design
  detail differs; the ordering of the three estimators is unaffected since
  all three run on identical data. These results do not validate the current
  semiparametric estimator.
- `results/pgd-ours-vs-his.csv` — the like-for-like comparison on a
  common-item linking design. Given common items and disjoint person
  groups this package does not use its person-side link at all: it
  estimates the person-group unit by conditional ML on the bilinear
  threshold structure. Run against Humphry's SD ratio on his own design
  (12 WALNA common items, N = 980 per year, ability gap 0.5, planted phi
  ratio 1.306), the two are indistinguishable on clean data: 1.306
  against 1.308, sampling SD 0.040 against 0.041. Under item-level
  departures conditional ML is consistently the worse of the two, by 3 to
  5 percentage points (1.458 against 1.406 with two DIF items; 1.480
  against 1.460 with four differentially discriminating items). Weighting
  by information is the reason: an item whose discrimination is inflated
  in one frame carries more information there, so conditional ML leans on
  the items that mislead it, where a dispersion weights items by squared
  distance from the mean. The conclusion is symmetrical. On a common-item
  design his estimator is at least as good as ours, and our advantage is
  confined to designs where item sets partition the items and his channel
  does not exist. It also settles a tempting change: switching phi to a
  dispersion ratio would buy a few points under contamination that
  screening should remove anyway, at the cost of the conditional standard
  errors an SD ratio cannot provide.
- `results/alpha-set-misfit.csv` — what actually threatens an item-set
  unit, and the largest effect measured anywhere in this battery. Sets
  partition the items, so no item appears in two sets and DIF across sets
  is undefined; the hazard is misfit concentrated in ONE set, which the
  other carries nothing to offset. With 8 items per set and a planted
  ratio of 1.40: two over-discriminating items in set 1 recover 1.17, two
  under-discriminating items recover 1.73, and four over-discriminating
  items recover 1.02 — a real 40 per cent unit difference read as none.
  The same misfit spread evenly across both sets cancels almost exactly
  (1.41 against 1.42 clean). The sign follows the mechanism: an
  over-discriminating item disperses its own set's person estimates, which
  reads as a larger unit for that set. The under-discriminating case is
  the practical one — it is the wording case study's Q8, an ambivalent
  item filed among the negatively worded ones, whose removal moved that
  analysis from 1.24 to 1.03. Set membership is thus the most
  consequential decision a user makes here, it is a falsifiable hypothesis
  rather than a given, and the diagnostic is ordinary item fit within each
  set.
- `results/misfit-repair.csv` — does the diagnose-and-drop workflow put a
  planted unit ratio back? Plant misfit, read the diagnostic, drop what it
  flags with `drop_items()`, refit, and compare against dropping the
  planted items regardless of what was flagged. Planted ratio 1.40,
  N = 500 per group, 200 replicates.

  Dropping is always a sufficient cure: the oracle recovers 1.406 to 1.430
  in every cell, against clean references of 1.408 and 1.423. Nothing is
  lost by removing an item, so the binding constraint is never the repair —
  it is whether the diagnostic finds the item. Recovery is not a simple
  function of sensitivity across the whole table, because the two damage
  directions push the ratio opposite ways and cancel when pooled; read the
  cells within a direction, where the cell with the lowest sensitivity is
  also the cell whose repair leaves the most damage behind.

  Sensitivity ranges from 91 to 22 per cent across the departures. Items
  with DIF across person frames are found 91 per cent of the time and the
  loop closes completely (1.448 damaged, 1.413 repaired, 1.406 oracle).
  Items that merely discriminate differently across frames are found 40
  per cent of the time and the loop half closes (1.530, 1.479, 1.406).
  Under-discriminating items concentrated in one item set are found 22 per
  cent of the time, the item fit test flags nothing at all in 114 of 200
  replicates, and the repair is nearly worthless (1.694, 1.638, 1.430).

  The multiplicity adjustment is the wrong instrument for this job.
  Screening ten items with Holm costs 8 to 42 points of sensitivity,
  and the loose screens recover more: `|infit z| > 2` lifts the
  under-discriminating cell from 1.638 to 1.486 (sensitivity 22 to 95 per
  cent), and unadjusted probabilities lift the differential-discrimination
  cell from 1.479 to 1.433 (40 to 82 per cent). The loose screen is not
  free. Where misfit is strong it over-flags: `|infit z| > 2` flags 17 per
  cent of sound items in the over-discriminating cell and `drop_items()`
  refuses 48 of 200 drops for emptying a set, so the surviving mean rests
  on 152 replicates and is not comparable with the rest.
  `fit-residual-screens.csv` later replaced this screen with the
  standardised fit residual, which detects as much without the
  over-flagging; read that entry before adopting `|infit z| > 2`.

  Screening should therefore be separated from confirming, which
  `frame_invariance()` did not allow at the time of this run: it flagged on
  Holm-adjusted probabilities only. Its `adjust` argument now offers both.

  One cell resists every screen. Four of ten items breaking invariance
  leaves 1.663 damaged against 1.410 oracle, and the best screen reaches
  only 1.602 while flagging a quarter of the sound items. Two of ten and
  four of ten are the only contamination levels this study runs, so where
  between them screening stops substituting for a coherent item set is not
  established here -- only that at two of ten it substitutes partially and
  at four of ten it does not.
- `results/fit-residual-screens.csv` — which standardised fit statistic
  should drive a screen. Both the frame comparison and the item-set screen
  used `infit_z`, the cube-root standardisation of a mean square; the
  package also computes the log-transformed fit residual
  `f (log y2 - log f) / sqrt(v)`. The answer differs by which comparison is
  meant, so the two were tested separately.

  Across frames, `infit_z` should stay. It detects two items discriminating
  1.8 times as steeply in 67, 95 and 100 per cent of replicates at 500,
  1,000 and 2,000 persons per frame against the fit residual's 50, 85 and
  99, at type I rates of 3.0, 4.3 and 5.5 per cent against 2.0, 2.8 and
  4.0. The extra power is largely the difference between a test at nominal
  size and a conservative one, and on ranking — whether the planted items
  are the largest departures, which is what a screen depends on — they are
  indistinguishable (80.5 against 78.5 per cent at 500 persons). An earlier
  single dataset suggested the fit residual ranked better; it does not
  generalise.

  Within a set the fit residual is the statistic to use, though the cut
  this study used does not generalise -- see
  `fit-residual-threshold-n.csv`, which supersedes the threshold advice
  below while leaving the choice between statistics standing. At the 500
  persons measured here, the question is not which test is more powerful
  but which is better calibrated at the conventional cut of 2. Against a planted 1.40 distorted to 1.69 by two
  under-discriminating items, `|fit_resid| > 2` recovered 1.439 where
  `|infit z| > 2` reached 1.486 and the chi-square test reached 1.638, with
  an oracle of 1.430. In the over-discriminating cell it recovered 1.376
  against 1.321 and 1.337, oracle 1.430. It detects as much as `infit_z`
  (88 and 93 per cent against 86 and 95) while flagging 6 and 2 per cent of
  sound items against 17 and 4, and that gap decides the repair:
  `drop_items()` refused 48 of 200 drops under `infit_z` for emptying a set
  against 4 under the fit residual.
- `results/fit-residual-threshold-n.csv` — why a fixed cut on a fit
  statistic does not survive a real sample size, and what to do instead.
  Run against the Rosenberg Self-Esteem data behind the wording case study
  (6,000 respondents, ten items), the `|fit_resid| > 2` screen recommended
  by `fit-residual-screens.csv` selects seven of the ten items, `|infit z|
  > 2` selects eight, and the chi-square test selects all ten, so
  `drop_items()` refuses every one for emptying a set. The recommendation
  did not survive contact with the data it was written for.

  The tempting explanation — real items never fit exactly while simulated
  ones do — is wrong. The cut degrades with N even when the sound items
  are generated from the model exactly: they clear it 0.8 times out of 14
  at 500 persons, 3.1 at 2,000 and 7.0 at 6,000. Two of eight items
  discriminating twice as steeply forces the fitted model to a compromise
  under which the rest genuinely depart, and that departure is fixed in
  size, so only its detectability grows. A fixed threshold on any fit
  statistic is a statement about power, not magnitude. Carried through to
  the repair it is worse than useless: the drop it implies is refused in
  61 per cent of replicates at 2,000 persons and 100 per cent at 6,000.

  Ranking survives what thresholding does not. The two planted items are
  the two largest departures in 100 per cent of replicates at 2,000
  persons and above when the others fit exactly, and 75 to 78 per cent
  when every item carries its own small slope departure. On the
  self-esteem data the ranking puts Q8 first at 22.5 and Q4 second at
  12.4, and a free-slope model fitted to the same respondents ranks the
  same two lowest — the two orderings agree on the extreme three items.
  Dropping Q8 alone moves the unit ratio from 1.24 to 1.03 at p = 0.63.

  The stopping rule has to be the unit test rather than the ratio, and this
  dataset shows why: dropping the second-ranked item as well does not
  improve on 1.03 but returns a significant difference of 1.16, still
  favouring the positive set. Once the sets no longer differ in unit there is nothing left
  for a further drop to explain, so the advice is to order by
  `abs(fit_resid)`, drop the largest, refit, and stop when the unit test
  goes quiet. `inst/casestudies/wording_units_selfesteem.R` runs exactly
  this sequence, so the figures above are reproducible rather than prose.
- `results/resolve-versus-drop.csv` — when an item breaks frame invariance,
  is it better removed or resolved? `misfit-repair.csv` showed that dropping
  a flagged item restores a planted unit ratio; it did not ask what the
  repair costs. Dropping takes the item out of every frame, so it stops
  contributing to any person's measure, including in the frames it behaved
  well in. Resolving gives it a location per frame: it stops linking the
  frames, which is what the diagnosis found wrong with it, and goes on
  measuring the person inside their own frame.

  One item shifted 1.2 logits in frame 2, twelve items, 500 persons per
  frame, planted group-unit ratio 1.40, 200 replicates.

  On the unit the two remedies are indistinguishable: 1.400 dropped against
  1.401 resolved, from 1.386 damaged and 1.399 clean. Both are unbiased, and
  the damage itself is modest — one differentially functioning item in twelve
  moves the ratio by about one per cent, so the repair is not what
  distinguishes them.

  The cost does. Dropping raises the mean person standard error from 0.7268
  to 0.7592, a 4.5 per cent loss of precision paid by every respondent, and
  lowers the correlation between the person estimates and the locations that
  generated them from 0.8501 to 0.8374. Resolving leaves both where the
  clean analysis had them, 0.7267 and 0.8498 — indistinguishable from never
  having had the problem. It buys this with one parameter per extra frame
  and by removing the item from the link, so the group units then rest on
  the items that remain common; `resolve_frames()` refuses when the remaining
  frame graph or information matrix no longer identifies those units.

  So the two remedies are not a trade-off on this evidence: where an item
  measures well within each frame and only its comparability fails,
  resolving dominates. Dropping earns its place where the item is a poor
  measure wherever it appears, which is a different diagnosis than the one
  `frame_invariance()` makes.
- `results/chained-linking.csv` — does a unit ratio recover when the two
  frames share no items at all? A vertical design rarely gives every pair of
  year levels a common block: year 3 and year 5 share one anchor, year 5 and
  year 7 share a different one, and years 3 and 7 share nothing. The model
  accepts this, because identification of the person-group units runs over
  the CONNECTED COMPONENTS of the graph whose edges are group pairs sharing
  at least two items within one set. A chain is enough. Breaking it is
  refused, naming the components: "relative units (phi) are unidentified
  between: year3+year5 | year7".

  Accepting a design is not recovering from it, and the year 3 to year 7
  comparison is made entirely through year 5. It recovers anyway. Against
  planted ratios of 1.250 directly linked and 1.562 chained, at 300, 700 and
  2,000 persons per year the direct ratio came back at 1.238, 1.262 and
  1.249, and the chained one at 1.536, 1.573 and 1.567. No bias worth the
  name in either, at any of the three sizes.

  What the chain costs is precision, not accuracy. The chained ratio's
  empirical standard deviation runs about 1.4 times the direct one at every
  sample size — 0.171 against 0.119, 0.115 against 0.080, and 0.068 against
  0.051. Contrast standard errors use \eqn{c^T V c}. They track the empirical
  spread closely from 700 persons per year; at 300, the indirect-link mean SE
  is conservative (0.374 against 0.171).
- `results/frame-invariance-conditional-topup.csv` — the 2,000-replicate
  null check that led to conditional discrimination probabilities being
  withdrawn. The standardised-infit comparison produced 7.1% combined Holm
  familywise error, with rejection strongly dependent on item position. The
  descriptive infit and fitted-slope columns remain useful, but they are not
  treated as tests.
- `results/frame-invariance-bootstrap.csv` — calibration of the replacement
  bootstrap inference at 500 persons per frame. Across 300 null replicates,
  empirical-SD/mean-SE ratios were 1.002 for locations and 1.034 for log
  discrimination ratios, coverage was 0.948 and 0.954, and the combined Holm
  familywise error was 3.0%. At two planted items, power was 96.3% for a
  one-logit location shift and 9.6% for a 1.5-fold discrimination change over
  120 replicates. The latter comparison is valid but weak at this design.
  This was a supported complete-category design rather than a screen of
  changing comparison sets; the probability and SE calculations are
  unchanged when the comparison family is retained. A sparse polytomous
  regression test separately fixes the boundary: a
  resample that gains or loses a comparison is counted as an other failure
  rather than being centred on a different set of items.
  The subsequent support-guard correction counts distinct persons with an
  informative conditional pair in each set and frame. Its regression tests
  (`test-frame-invariance-support.R`) cover mixed categories, missing responses,
  repeated IDs, and extreme-response padding for dichotomous and polytomous
  data, including frames missing a lowest or highest category. Support is
  checked after the separate calibration's item removal and category recoding;
  the minimum is not reapplied to select bootstrap replicates.
  These are conformance checks; the simulation results above have not
  been rerun under the revised guard.
- `results/frame-invariance-power.csv` — the conditional location study,
  with 1,000 null and 400 departure replicates at 500, 1,000 and 2,000
  persons per frame. Empirical-SD/mean-SE ratios were 0.908--0.916, coverage
  was 0.966--0.968, and Holm familywise error was 1.9--2.6%. Holm-adjusted
  power for each of two one-logit shifts was 95.1% at 500 persons and 100% at
  1,000 and 2,000. Because separate frame origins centre the common items,
  unshifted items can carry non-zero relative contrasts when a few items
  move; the study records that rate rather than treating it as an ordinary
  false-positive rate.
- `sha-map-2026-08-16.txt` — commit-ID map (old, new) from the 2026-08-16
  message-only history rewrite. `package_sha` values stamped in result
  tables before that date are pre-rewrite IDs; look them up in the first
  column to find the corresponding commit in the current history. File
  contents were untouched by the rewrite, and the `r_tree_md5` content
  hashes remain directly verifiable.

## Conventions

### Principal-component kurtosis (September 2026)

`studies/pc-kurtosis.R` and `results/pc-kurtosis.csv` check the corrected
fourth component against the published adjacent-threshold polynomials.
Across 100 five-category datasets, the maximum difference from unrestricted
PCM was 3.21e-8 logits for thresholds, 1.20e-9 for their covariance and
1.28e-11 for log likelihood. All four components were identified.

Fixed-truth recovery used 100 datasets at each of four, five and six
thresholds, with 1,200 persons and six items. Mean item-wise 95% kurtosis
coverage was 95.5%, 94.0% and 95.2% (Monte Carlo SE 0.74, 0.99 and 0.83
percentage points). The I3 empirical-SD/mean-SE ratios were 1.06, 1.09 and
0.83. All 300 fits converged with available estimates. These limited
recovery checks are not a principal null-size study. Earlier PC checks
through three thresholds did not exercise kurtosis and do not validate
the previous fourth-component formula. Result hashes identify the code
loaded for this run, not later unrelated edits in the working tree.

### CJ panel covariance (September 2026)

`studies/btl-panel-covariance.R` and `results/btl-panel-covariance.csv`
compare conditional panel-unit variances with independently stacked
stage-one judge scores in 50 fitted designs. These cover shared, partially
shared and disjoint judge pools, with positive and negative cross-set
dependence. All comparisons agreed within 1.4e-17; there were no refusals,
non-convergences or errors. Point estimates were unchanged to numerical
precision. The former SE divided by the corrected SE averaged .81 under
shared positive dependence, 1.29 under shared negative dependence, and 1.00
for disjoint pools. This is a covariance-conformance check, not a coverage
or null-rejection study; the default judge bootstrap is unchanged.

### Repeated-person explanatory calibration (September 2026)

`results/explanatory-repeated-person.csv` records 4,500 runs with the
linearised delete-one-person covariance. Balanced null conditions rejected
5.2%, 3.4% and 6.0% of tests (1,000, 500 and 500 runs). With twenty persons
contributing alternating one and four rows, rejection was 6.35% over 2,000
runs (Monte Carlo SE 0.55 percentage points), and 95% coverage was 93.65%.
The correction therefore remains mildly liberal in this small, unbalanced
condition. The 500-run coefficient-power condition gave 70.6% power.

`results/explanatory-cluster-jackknife.csv` compares covariance methods on
the same draws. The linearised and full delete-one-person methods both
rejected 5.0% in its independent 500-run unbalanced null condition; their
empirical-SD/mean-SE ratios were 1.028 and 1.026. This supports the numerical
approximation, not a claim of exact small-sample calibration.

Earlier tailored-bootstrap results predate retention of no-tailoring draws
as zero shifts. They remain a record of the previous algorithm, not validation
of the revised bootstrap. Focused regression checks cover zero draws, retained
item structure and non-convergence accounting; updated operating-characteristic
studies are still needed for the revision.
`results/tailored-zero-draws.csv` checks that boundary directly: 130 of 200
draws at 35 persons and 95 of 200 at 150 persons required no tailoring and
were retained as zero changes. Total usable counts were 193 and 198; the
remaining 7 and 2 were other failures. These are accounting checks, not
Type I error estimates. The main and top-up scripts now seed the inner bootstrap
explicitly, so future results do not depend on the worker's starting RNG state;
their earlier results predate that reproducibility fix.

### Structural refit scoring (September 2026)

`results/structural-score-preservation.csv` checks the structural-refit guard
against explicit PCM and EFRM refits over 80 randomised input designs. All
40 control designs retained their scores and remained usable. All 40 boundary
designs lost a category during the explicit refit and were refused by the
structural wrapper. There were no setup errors or disagreements. This is an
algorithm-conformance check, not a Type I error or coverage study. Regression
tests also cover explanatory item dropping, dependence resolution, superitem
formation and a supported EFRM single-group boundary where no rescoring occurs.

### CJ information and comparison context (September 2026)

`results/btl-information-context.csv` checks design information against
numerical likelihood curvature in 120 fits: dichotomous and four-category
comparisons, each with true position effects of -1.2, zero and 1.2. The
largest absolute difference was 6.2e-8. Omitting the fitted position effect
overstated total information by 16.9--18.1% on average in the non-zero
conditions. This checks the information calculation, not Type I error or
coverage. Regression tests separately check history-conditioned information,
observed-opponent response curves and position-adjusted recommendations.

### Partial frame and facet information (September 2026)

`results/structural-information-patterns.csv` checks 50 EFRM and 50 MFRM
response-pattern designs against independently computed numerical likelihood
curvature. Each uses eight response cells with mixed category counts and
40 missingness patterns; EFRM also uses unequal set units. Every curve retains
exactly its observed items. The largest absolute information discrepancy was
2.1e-7. These are fixed-calibration conformance checks, not coverage or Type I
error studies. Regression tests also check SEM, item selection, plotted expected
totals and agreement with an estimated EFRM's stored score curves.
Person-item map tests also check that selecting a subgroup excludes response
patterns found only outside that subgroup, without changing the calibration.

### General reporting

- True generating parameters are held fixed across replicates within a
  scenario; only responses (and persons, where stated) are redrawn.
  Redrawing the truth each replicate folds between-replicate truth
  variation into the empirical SD and fakes a calibration failure.
- Principal null-calibration claims target >= 1,000 replicates (Monte
  Carlo SE ~0.7 percentage points at a true 5% rate); power and secondary
  conditions use fewer, with the Monte Carlo error reported alongside.
- Every study records refusals and non-convergences rather than silently
  dropping them.
