# Changelog

## rasch 1.12.1

- Automatic DIF resolution withholds its verdict when the final
  assessment could not estimate a single item-term test, instead of
  reporting that no DIF remains. Where only some tests were estimable
  the stop reason says how many were not, and the counts of remaining
  and non-uniform DIF are reported as unavailable rather than as zero.
  Saved resolutions from earlier versions are dropped when a project is
  reopened.

- Planned DIF contrasts withhold only the affected item when a split
  refit is refused, and say why, instead of abandoning the whole call.
  The refusal reason is reported in place of an unrelated
  response-category note.

- A linear trend across an ordered factor is no longer withheld because
  the middle level carries the rounding residue that
  [`contr.poly()`](https://rdrr.io/r/stats/contrast.html) leaves in
  place of an exact zero. Contrast weights below a relative tolerance
  are zeroed, so a level with no weight is no longer treated as a
  required cell. Saved planned contrasts from earlier versions are
  dropped on reopening.

- The notes recording which requested DIF tests were not estimable reach
  every surface that reports the analysis: the console, the HTML report,
  the Word and PDF report, and the saved summary, not only the app.

- Item and person tables in the app colour and embolden cells by their
  own value again. A conditional style was passed to the table as a
  function where an expression was expected, so every flag was decided
  from the row index: whole columns read as misfitting or significant,
  and residual and DIF magnitude flags never appeared.

- Switching an analysis to Comparative Judgement clears the Rasch
  results that do not carry over, so a CJ analysis can be saved after
  planned contrasts, automatic DIF resolution or a tailored analysis was
  run on the same data.

- Uploaded equating references, Comparative Judgement banks and item
  panel maps are checked when they are supplied rather than when the
  analysis is saved. An unusable upload is refused at the point of
  upload and leaves the rest of the analysis saveable.

- The app’s dependence and spread panels print the note explaining why a
  statistic was withheld, as the console does, rather than an
  unexplained NA.

- [`sim_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/sim_recovery.md)
  accepts a fit made without an identifier column, matching persons on
  their responses rather than on labels the fit never held, and accepts
  a fit whose estimator dropped an item with no variation. Both cases
  previously failed with a claim that the fit did not belong to the
  simulation. A person-to-group allocation that differs from the planted
  one is still refused.

- Test characteristic, test information and targeting displays enumerate
  the administration designs actually observed. Item-level missingness
  no longer produces one curve per person, and a design that no person
  was administered is never formed. Where more designs exist than the
  palette can name, the legend counts them instead of cycling colours.

- `test_information(items = )` labels each curve with the items that
  produced it, and no longer returns a curve twice.

- Extended frame score curves use the same design enumeration and labels
  as
  [`test_information()`](https://drjoshmcgrane.github.io/rasch/reference/test_information.md)
  and the curve plots, so a fit describes its designs one way. A design
  label always names its item set, and a set or item name containing the
  separator is quoted. Saved extended frame analyses holding the earlier
  curves must be refitted before reopening; their source files remain
  unchanged.

- The EFRM unit-test family counts the hypotheses that are free. With
  two groups or two sets the second reported coordinate restates the
  first, and now leaves its adjusted probability and flag to the row it
  restates instead of being counted twice.

- `frame_invariance(se_method = "bootstrap")` refers its Wald statistics
  to a t distribution on the bootstrap degrees of freedom rather than to
  the normal. A 30-draw bootstrap rejected an invariant item about 6.6
  per cent of the time at the nominal 5 per cent.

- The two-panel Comparative Judgement frame family counts the single
  panel-unit hypothesis once, so a panel-unit difference is reported
  with the evidence it has and the step-down multiplier for the scale
  rows is the number of distinct questions.

- [`btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)
  records that its precision-weighted shift standard error treats the
  estimated weights as fixed, which understates the shift’s uncertainty
  in small panels.

- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  no longer reports the judge-clustered chi-square boilerplate as a
  finding about a particular object.

- Explanatory coefficient tests refer to a t distribution on the
  residual degrees of freedom for a fit without repeated persons. The
  normal reference rejected a true null about 8.6 per cent of the time
  at the nominal 5 per cent with the smallest supported person count.

- The item-fit bootstrap keeps the observed class-interval count in
  every replicate under the `"resample"`, `"normal"` and `"fixed"`
  generators, so the null and the observed statistic use the same number
  of intervals. Probabilities no longer depend on which side of a
  50-person boundary the sample happens to sit.

- [`dimensionality_test()`](https://drjoshmcgrane.github.io/rasch/reference/dimensionality_test.md)
  reports the difference between the two subsets’ mean estimates as a
  description and no longer as a paired t-test. The test rejected on
  every unidimensional replicate whenever the subsets differed in
  difficulty, which is targeting rather than dimensionality. Saved
  person-subset tests are dropped when a project is reopened.

- Printing a dependence result saved by an earlier version reports its
  statistic again instead of failing, and a withheld item-trait
  probability prints as unavailable rather than as an empty value, in
  the console and in the saved summary alike.

- Notes containing parentheses or brackets appear as prose in the Word,
  PDF and HTML reports. They were escaped in a way the report renderer
  read as mathematics.

- The simulation-validation studies run from the tree as committed: the
  conditional DIF bootstrap confirmation pins the committed helper, the
  item-fit interval checker verifies against the current sources, and
  the README describes the recorded dimensionality guard results and how
  to reproduce them.

- DIF magnitudes leave the ETS category unavailable when a required test
  probability is unavailable or its standard error is zero, rather than
  treating the missing test as evidence for category A or B.

- EFRM conditional calibration and CJ frame panel-ratio fits check the
  exact likelihood curvature as well as the score. Stationary saddles
  can no longer supply converged fits, unit inference or bootstrap
  covariance. Previously saved frame analyses without a current
  likelihood-check record must be refitted before reopening; their
  source files remain unchanged. Resolving every common item in a set is
  also refused when it leaves the groups’ relative origins unidentified,
  even if another set links their units.

- CJ frame linking checks partly separated outcomes, including
  comparisons with zero loading along a separating direction. Such links
  have no finite estimate. An exhausted iteration limit can no longer
  certify convergence.

- The installed Shiny app resolves the internal helpers needed for item
  tables, saved CJ analyses and simulation recovery, including
  standalone launches.

- CJ dimensionality compares observed and fitted expected points within
  each object pair, retaining category thresholds and fitted position or
  history effects. Simulated histories use their own fitted
  expectations. Older saved dimensionality results require
  recalculation; fitted models are unchanged.

- EFRM bootstrap replicates are discarded when any supported set link
  fails numerically or does not converge, including redundant links.
  Failed draws cannot enter the covariance or count as usable
  replicates. Full-bootstrap refits with unidentified group units are
  also discarded.

- CSV, HTML, Word and PDF exports include supplied DIF results based on
  external person factors. Bulk exports require a new or empty folder,
  preventing files from separate analyses being mixed.

- CJ pair-surprise diagnostics do not assign a stronger or weaker
  direction to tied object locations. Relabelling objects no longer
  changes these flags.

- Explicit item selectors can distinguish item columns from separately
  supplied person factors with the same names in Rasch and EFRM fits.

- CJ model comparisons recognise unchanged within-judge sequences after
  relabelling their order values.

- EFRM item sets can link through intermediate sets without direct
  overlap between every pair. Insufficient overlap is distinguished from
  a numerical link failure; the retained linking graph must still
  connect all sets.

- Structural Rasch and EFRM refits preserve external factors whose names
  match item names, including frame-invariance bootstrap refits.

- Downloaded simulation scripts load `rasch` before recreating the data.

- EFRM DIF bootstraps preserve frame and person-factor names, including
  an ordinary factor named `group`.

- Calculations that restore random-number streams refuse the Box-Muller
  normal generator before changing the stream. R does not expose its
  cached normal value for restoration. The default Inversion generator
  is unaffected;
  [`?rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md)
  describes the supported behaviour.

- BTL-EFRM documents its different treatment of response-boundary
  objects and distinguishes separated outcomes from disconnected
  comparison designs.

- Saved CJ DIF results without verified judge-role alignment are omitted
  on reopening. Reports refuse older automatic DIF magnitudes that may
  have used stale factor values; the analysis must be rerun.

- Bootstrap workers retain the coordinator’s random-number generator
  settings, including non-default generators, so a fixed seed gives the
  same draws in serial and parallel execution. Installations with source
  references retained also select the coordinator’s package library
  before loading the worker namespace.

- BTL simulation recovery checks the original comparisons before extreme
  objects were set aside, then assesses the calibrated objects only.

- Principal-component calibration uses the corrected kurtosis
  polynomial. The full four-component model matches unrestricted PCM
  through four thresholds; longer scales retain a restricted polynomial
  structure. Saved fits that used the earlier kurtosis polynomial
  require a refit.

- Automatic DIF magnitudes retain the factor values supplied to that
  analysis, including replacements for factors stored in the fit.

- BTL-EFRM panel-unit reconciliation retains covariance between object
  sets assessed by the same judges. Point-estimation weights are
  unchanged.

- Saved fits record their person-scoring algorithm. Older fits are
  checked before reopening; materially changed scores require a refit
  rather than silently restoring superseded estimates and diagnostics.

- CJ DIF in the app uses the judge role from the fitted analysis, not a
  subsequently changed column selector. Unavailable overall fit
  probabilities are labelled explicitly in summary tables.

- Incomplete repeated-measures DIF fits occasion and person factors
  jointly, with equal total weight per person and person-cluster CR3
  covariance. Marginal occasion centering no longer introduces group
  effects when occasion coverage and person-factor composition differ.
  Complete-panel analyses retain their existing covariance references.

- WLE scoring compares competing maxima for separated item banks instead
  of accepting the first score-equation root. This applies to ordinary
  and extended-frame scoring, including externally weighted measures.
  Saved weighted tables from earlier solvers are authenticated and
  recomputed.

- BTL-EFRM refuses a linking unit at the zero boundary rather than
  replacing it with one. The same check applies inside bootstrap refits;
  genuinely indistinguishable within-set locations retain their
  placement convention.

- CJ model comparisons check each fit against its reference separately.
  Adding another fit cannot erase a known difference in comparison
  order.

- Judge-clustered CJ dimensionality is descriptive by default.
  `independent_comparisons = TRUE` requests the conditional sensitivity
  reference; it does not account for general within-judge dependence.
  Unsupported saved references and superseded mixed-panel DIF results
  are omitted on reopening, without removing the source data or fitted
  models.

- Explanatory diagnostics use stable added directions and retain failed
  candidates in the Holm family, allowing other departures to be
  reported.

- Loading another dataset in the app clears the active analysis and its
  derived results. Explicitly kept model snapshots remain available.
  HTML reports label an unavailable total-fit probability explicitly.

- CJ dimensionality withholds probabilities and reference bands wherever
  its design checks withhold inference, including incomplete pair
  coverage in BTL-EFRM. The observed residual decomposition remains
  available. Saved analyses omit superseded unsupported references after
  integrity checks, retaining the data, fitted models and analysis
  history.

- Explanatory Rasch and CJ fits accept distinct predictors whose
  generated coefficient labels coincide. Labels receive unique suffixes
  without changing the design, estimates or model comparisons. Fixed
  departures also avoid existing coefficient names.

- R code for a completed CJ DIF analysis retains the judge-factor maps
  used in that run, including metadata corrected after calibration.
  Reopened analyses retain those assignments for subsequent DIF runs;
  new metadata can still replace them.

- Structural refits and bootstrap procedures refuse older anchored or
  component-constrained fits when their original refit settings are
  missing. They no longer release fitted restrictions silently.
  Simulation recovery also retains an anchored fit’s origin when those
  settings are unavailable.

- [`lr_test()`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)
  checks restrictions recorded by the fitted estimator as well as its
  refit settings. Older anchored or principal-component fits no longer
  bypass the requirement for an unrestricted PCM comparison.

- App reports and result archives stop when a selected DIF analysis or
  custom dimensionality test is unavailable, rather than substituting a
  default analysis. Analysis files can still save settings without a
  computed result.

- [`sim_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/sim_recovery.md)
  refuses comparisons after fitted categories have been removed or
  merged. Paired-comparison frame recovery aligns generating origins to
  the fitted reference set and refuses a changed reference whose
  generating unit is not one. The model help identifies the alphabetical
  reference-set rule explicitly.

- EFRM and externally weighted person scoring retain accuracy under
  changes of measurement unit. The score equations and SE formulas are
  unchanged. Saved weighted tables from earlier solvers are
  authenticated and recomputed when an analysis is loaded.

- Keyed analyses with repeated person IDs now export the available
  results and an explanation that distractor analysis is unavailable,
  instead of stopping the export or leaving an empty report section.

- BTL-EFRM unit tables in rendered reports show very small adjusted
  probabilities as `< 0.001`, not zero.

- Dichotomous CJ information remains symmetric when a pair is reversed,
  including in the far tails. Targeting plots refuse a grid whose
  reference information is numerically unavailable. The help now
  distinguishes the zero-gap information peak for dichotomous
  comparisons without a position effect from the more general
  ordered-comparison case.

- MFRM person factors supplied as a data frame retain each person’s
  observed value when other rows are missing. Results no longer depend
  on which response row comes first; conflicting observed values are
  still refused.

- Common-unit person scoring scales its root-search tolerance with the
  unit and calculates SEs without squaring the discrimination. Large
  changes of unit no longer produce inaccurate WLEs, MLEs or
  extrapolated-score SEs.

- Missing-data codes apply to raw responses, not scores generated by a
  key or earlier category renumbering. Keyed scoring and structural
  refits no longer delete valid scores that coincide with an original
  missing code. Affected analyses should be rerun from the raw data.

- [`combine_items()`](https://drjoshmcgrane.github.io/rasch/reference/combine_items.md)
  validates its requested model for every fit. Explanatory fits refuse
  an RSM request rather than silently returning a different model; their
  default refit still retains predictors and frees new superitems.

- Tailored analysis retains keyed scoring and DIF-split records in its
  returned fits. The tailored option data exclude censored responses,
  keeping subsequent distractor analysis and refits on the retained
  data.

- The RSM fit returned by
  [`lr_test()`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)
  retains keyed scoring, DIF splits and superitem definitions for
  subsequent analyses.

- Alpha if deleted checks each reduced scale’s covariance matrix
  separately. Missing covariances involving a deleted item no longer
  suppress an otherwise estimable reduced-scale alpha; unsupported
  reductions remain unavailable.

- Saved analyses authenticate and recompute weighted person tables from
  the earlier solver, preserving the data, fits and requested weights.
  Current tables retain strict integrity and reproduction checks.
  Available-case CTT summaries with missing responses withhold SEM
  rather than combine pairwise alpha with a complete-case score
  variance.

- Relaxing explanatory restrictions preserves earlier DIF-split and
  superitem records for subsequent anchor checks and dependence
  diagnostics. Weighted person estimation rescales weights within each
  response pattern, preventing missing estimates when only very
  small-weight items were answered.

- Refits preserve automatic and explicitly requested class-interval
  rules. Item-fit bootstraps retain per-item allocation with missing
  responses; earlier Rasch bootstrap results must be recomputed before
  reuse. Category plots, plot exports and downstream refits accept fits
  whose tied locations produce one interval. Fixed tailored fits are
  labelled as person scoring against a fixed item calibration. The
  missing-response validation now includes 200 datasets and 99,900
  usable bootstrap refits under the corrected rule.

- The final tailored scoring fit records its fixed thresholds
  explicitly. Downstream item changes and refit-based bootstraps refuse
  fully anchored scoring fits, including older saved fits, instead of
  releasing their anchors.

- Explanatory Rasch and CJ formulas now refuse unsupported
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms instead of
  silently omitting their fixed contributions.

- Model comparisons withhold generic EFRM likelihood differences, which
  need not use the same item pairs. The dedicated group-unit comparison
  is unchanged. CJ comparisons now require the same declared response
  scale, including endpoint categories absent from the observed
  responses.

- Explanatory refits align keyed responses with retained persons.
  Simulated scree and DIF refits no longer inherit observed answer
  options; blank response rows no longer make every keyed-response scree
  replicate fail. Predictor matching and item relaxation preserve exact
  item names, including leading or trailing spaces, and refuse ambiguous
  whitespace-normalised matches.

- Unused ordinal columns in predictor tables no longer prevent
  explanatory Rasch or CJ models from fitting. Ordinal contrasts are
  applied to the predictors included in the formula.

- Explanatory Rasch and CJ models use rescaled design columns for
  estimation and numerical checks. Changing predictor units no longer
  causes premature convergence or prevents Newton updates. Coefficients,
  covariance and model comparisons retain the supplied predictor units
  and equivalent restrictions. Stable column centring also prevents
  large constant offsets from making a varying predictor appear
  unidentifiable.

- Frame-invariance tests count distinct persons contributing informative
  item pairs within each set and frame. Adding all-zero or all-maximum
  responses no longer enables inference in an otherwise undersupported
  frame. Support uses the separate calibration’s retained items and
  category structure, including when a frame lacks an endpoint category.

- Equating retains the shift fixed by an exact common anchor when other
  common items or objects have unavailable SEs. CJ equated-location SEs
  now include shift uncertainty and its covariance with the reference
  locations; unavailable joint uncertainty is reported rather than
  ignored.

- Dimensionality magnitudes compare reliability on matched
  complete-response rows, retaining both calibrations. Missing subscale
  scores no longer change the PSI comparison sample. Tables report rows
  used and excluded; saved app results from the earlier calculation must
  be recomputed.

- Named-column arguments in reshaping, CJ, and Multiple Ratings require
  character strings. Numeric selectors can no longer pass a name check
  and then read a different column by position. Names such as `"1"`
  remain valid.

- Repeated-measures DIF reports unavailable within-person tests as `NA`
  when missing factorial cells confound their adjusted mean with the
  between-person design. These tests remain in the multiplicity family.
  Missing within-person combinations no longer abort the analysis, and
  item-specific missing levels no longer change the test’s requested
  design.

- Group-restricted person-item maps now show only the EFRM or MFRM
  response patterns observed in the selected subgroup. Information
  curves from other subgroups no longer accompany that group’s person
  distribution.

- EFRM and MFRM information and test-characteristic curves now retain
  each observed item pattern. Partly answered sets or facet conditions
  no longer contribute unobserved items, overstate information, or
  understate SEM.

- Crossed WrightMap panels keep distinct groups separate when their
  combined labels coincide. Supplied person-panel estimates must be real
  numbers. Large crossed-factor designs also retain exact cell
  identities when their potential combinations exceed integer precision
  in floating-point numbers. Grouping keys also treat column names such
  as `collapse` as data, so these facet names cannot merge separate MFRM
  administration designs.

- CJ design information now retains fitted position and history effects.
  Position-adjusted pair recommendations use the displayed presentation
  order. History-dependent fits no longer offer context-free
  recommendations or a single new-comparison information curve. Their
  response plots also retain these effects, comparing observed means
  with fitted expectations for each opponent and, when selected, judge
  group. Graded category curves include position effects and explicitly
  label the full linear-predictor axis for history-dependent fits.

- Structural refits now check the scoring structure after calibration,
  not just the categories present beforehand. Dropping or combining
  items cannot silently rescore retained items when a category loses
  conditional information. EFRM checks this within each frame;
  explanatory refits also preserve inherited item scoring.
  Dependence-magnitude refits apply the same check to every resolved
  copy before comparing their thresholds.

- Tailored bootstraps retain no-tailoring resamples as zero shifts and
  count non-convergence at each calibration stage separately from other
  failures. Saved results from the earlier bootstrap are omitted on
  restore. The app starts at 999 draws and displays warnings when too
  few usable draws make Holm-adjusted significance unattainable.

- Rasch simulation keeps numbered groups in numerical order, so DIF
  targets the intended last group when there are ten or more groups.
  Not-reached responses are documented as independently selected missing
  tails, not incorrectly labelled as non-ignorable departures from the
  response model.

- DIF resolution no longer reports zero remaining DIF when its final
  diagnostic assessment fails; it reports the failed assessment instead.
  Item splitting also checks the fitted category structure: a
  within-group category observed only in extreme patterns can no longer
  be merged and presented as a comparable split calibration.

- The conditional DIF bootstrap now compares the observed and replicated
  F-reference probabilities for both its marginal and familywise
  results. It previously used raw F values for the marginal result even
  when sparse refits changed the term’s degrees of freedom, making
  unlike statistics comparable. Results record the revised algorithm and
  saved analyses drop bootstrap DIF output from the superseded
  calculation.

- Tailored analyses now carry fitted-model and result signatures. Saved
  analyses omit older unauthenticated tailored tables, and current
  tailored item shifts are retained in app reports and result archives.
  Generated app code preserves whether anchors were selected
  automatically or supplied, so bootstrap reruns use the same
  anchor-selection procedure.

- Frame-invariance bootstrap replicates now retain the exact observed
  item comparison family. A category that was absent from an observed
  frame can no longer reappear in a resample and silently change the
  items used to centre that replicate’s frame origin. The result records
  the revised algorithm, and saved pre-correction frame-invariance
  results are omitted on restore.

- Multi-set EFRM fits with `boot_reps = 0` now withhold common-unit item
  and threshold standard errors. Those errors require the omitted
  set-link uncertainty; reporting the conditional calibration component
  alone made them too small. Their set-unit inference flag now also
  records that the covariance is unavailable, rather than describing
  sample support alone. Single-set analytic errors are unchanged.
  Frame-unit plots also omit confidence intervals whenever the fitted
  unit families do not meet their inferential support conditions; a
  finite descriptive standard error is no longer drawn as an inferential
  interval.

- Calibration-sandwich support and row-based diagnostic support are now
  distinguished. A repeated person row that contributes a fitted
  residual but no conditional item pair still withholds row-independent
  fit probabilities and null generators, without incorrectly changing
  the calibration covariance. Residual DIF procedures use the same
  row-based distinction. Their paired contrasts now withhold inference
  when a required between-person stratum cannot supply a variance,
  rather than dropping it and changing the planned equal-stratum
  estimand.

- BTL sandwich inference now checks the empirical score-covariance rank
  for both independent comparisons and judge clusters, rather than
  relying on the number of sampling units. Repeated identical profiles
  can no longer produce zero standard errors and highly significant
  probabilities when the sandwich has no variation in the fitted
  directions. With identified judges, the row-based pair-fit chi-square
  remains descriptive and its probability is withheld; the same applies
  to BTL–EFRM, where judges are required. Object-location plots use the
  finite judge reference carried by the fitted covariance. BTL DIF
  follow-up adjustment retains every contrast opened by a significant
  term even when no resolved cell meets its support threshold. In a
  one-panel BTL–EFRM fit, the structurally fixed panel unit is no longer
  counted as a hypothesis in the Holm family for the estimated set units
  and origins.

- Explanatory Rasch coefficient tests now retain repeated-person cluster
  support from the conditional calibration. Supported repeated-person
  fits use a linearised delete-one-person covariance correction and a t
  reference with degrees of freedom equal to the number of person
  clusters contributing conditional information minus one; unsupported
  fits withhold coefficient inference, while supported fits without
  repeated identifiers retain the limiting normal reference. The
  multivariate Kent model comparison remains a first-order asymptotic
  test.

- Pairwise conditional Rasch calibrations now withhold item-parameter
  covariance and standard errors when the independent or effective
  person count, parameter count, or empirical score rank cannot support
  them. The same conditions apply to person clusters when IDs repeat. A
  single shared ID can therefore no longer produce false zero standard
  errors that equating would treat as exact; genuinely fixed anchor
  locations remain exact. Unit and effective-unit counts use informative
  conditional item-pair contributions, so appending empty or
  deterministic response rows cannot manufacture inferential support.
  EFRM conditional-stage covariance now applies the same person-support
  and projected-rank checks; unsupported fits retain descriptive point
  estimates only and cannot enter hybrid or full-bootstrap unit
  inference. EFRM also refuses repeated supplied person identifiers: its
  set-link likelihood and person bootstrap require one response row per
  person. Missing identifiers remain separate unknown persons.
  Common-item equating from supported clustered calibrations now uses
  contrast-specific Welch–Satterthwaite degrees of freedom rather than a
  normal reference. Resolved DIF magnitudes likewise use the independent
  person-cluster count minus one for their probabilities, confidence
  intervals and ETS interval-null decisions. Response-dependence
  magnitudes and the spread-parameter screen use the same finite-cluster
  reference. Ordinary and paired-comparison equating plots use the
  corresponding contrast-specific critical values for their interval
  displays. MFRM fits now retain the structural calibration’s Godambe
  ingredients, so
  [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  reports CL-AIC and CL-BIC for comparable MFRM structures rather than
  withholding them with the two-stage EFRM criteria.

- Fixed-origin equating now refuses contradictory exact anchors instead
  of averaging them into a spurious shift and apparent drift in the
  remaining items or objects. A location-anchored polytomous item also
  retains its exact zero item-location standard error when sparse
  categories make its free threshold-spread standard errors unavailable.

- Missing-response codes are now matched before numeric conversion in
  ordinary Rasch and many-facet analyses. Text labels and numerically
  equivalent forms, such as `"09"` and `9`, therefore behave
  consistently. Superitem refits also discard a stored missing code when
  it would become a valid summed score.

- Paired-comparison order diagnostics now assess variation between
  judges, not repeated presentations within a judge. Position-only
  effects are not treated as sequential dependence, and simulated
  comparisons use independent judge orders so planted dependence remains
  assessable. Static first-position effects are labelled separately from
  within-judge history effects in printed, downloaded, and app
  summaries.

- Fixed-origin equating can use one common item or object; estimating a
  shift still requires at least two. A normal-parametric fit bootstrap
  now uses a point mass when corrected person variance is zero and
  refuses an unavailable variance instead of silently changing to an
  empirical resampling method.

- Legacy paired-comparison dependence results without adjusted
  probabilities now report that Holm probabilities are unavailable in
  direct printing, plots, and the app. The same rule applies to
  residual-dimensionality printing. Raw probabilities are not
  substituted for adjusted inference.

- Weighted person scoring now refuses set maps that become duplicated
  after item-name matching, rather than silently replacing one
  assignment. Saved application results are authenticated against the
  active calibration and retained in HTML, Word and PDF reports and the
  results archive as externally weighted secondary person measures. The
  archive also carries the resolved response-cell weights needed to
  interpret and reproduce them, and a successful switch to Comparative
  Judgement clears an earlier Rasch weighted table with the rest of the
  calibration-specific results. A table returned by
  [`weighted_person_estimates()`](https://drjoshmcgrane.github.io/rasch/reference/weighted_person_estimates.md)
  can be supplied explicitly through `person_weights` when saving
  outputs or writing a report; it is reproduced against the fitted model
  before use.

- [`simulate_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_mfrm.md)
  now treats `item_sd` as the realised sample standard deviation of its
  item difficulties, consistent with the argument name and with the
  other simulator scale controls. A value of zero still gives equal item
  difficulties.

- Simulation controls now refuse departures that the generated design
  cannot distinguish: DIF must leave an invariant item, EFRM drift must
  leave one in each affected set, a second attribute cannot duplicate
  the first exactly, and random or halo contamination cannot replace
  every identifying group, panel, attribute camp or sampling unit. A
  response style must leave an unstyled person because a category
  weighting applied to everyone is absorbed by the fitted thresholds.
  Halo needs differing item difficulties, while guessing, heterogeneous
  slopes and non-uniform DIF need person variation. Common item
  discrimination is recorded as a change of logit unit. Recovery output
  labels any deliberate departure outside the fitted model as a
  descriptive comparison rather than an unbiased recovery target.
  Partial-credit threshold structures are therefore identified as
  unsupported when the same simulated responses are fitted with a
  rating-scale model, and MFRM recovery checks every fitted facet
  allocation carried by the simulation. The application also requires
  slope departures to leave one reference item. A nonzero
  local-dependence request must contain an item pair, and scale controls
  described as realised standard deviations are refused when their
  spread cannot be represented numerically.

- Fit summaries and the app now distinguish ordinary pairwise
  conditional calibration, response-cell calibration, and semiparametric
  EFRM set linking.

- Generic model comparison now withholds information criteria and
  likelihood differences for BTL–EFRM fits. Their combined objective is
  evaluated at a two-stage estimate rather than maximised jointly; the
  fit’s `equal_unit` component retains the corresponding labelled
  descriptive comparison.

- Explanatory diagnostics and fixed relaxations now add only the part of
  a polytomous threshold-departure block not already represented by the
  active predictor design.

- BTL–EFRM judge-bootstrap degrees of freedom now follow each set’s
  supported path to the reference set. A strong terminal link can no
  longer mask a weak link earlier in the chain, including when the
  fitted locations are used for common-object equating. Set pairs below
  `min_link` are now omitted from the linking fit, as documented, rather
  than entering the likelihood after being excluded only from the
  connectivity check. Unit plots use these row-specific references for
  their intervals and omit an interval when inference is unavailable.

- EFRM set-unit support now follows the strongest bottleneck path
  through the linking graph. A strong terminal edge can no longer lend
  support across a weak upstream link, and a weak redundant edge no
  longer suppresses a set that has a stronger route.

- EFRM score curves are now keyed to the exact observed-item pattern
  within each group. A partial booklet form no longer receives the
  expected weighted score and information of every item in a set merely
  because it contains one item from that set.

- EFRM now refuses repeated non-missing person identifiers because its
  set-link likelihood and person bootstrap require one response row per
  person. Missing or blank identifiers remain distinct unknown persons.

- [`simulate_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_efrm.md)
  now refuses a zero person spread when it generates more than one item
  set, because relative set units then have no shared person variation
  from which to identify their scale. It also refuses non-zero item
  drift with one person group, where the shift is an ordinary
  item-location change rather than a violation of invariance across
  frames.

- Item-fit and dimensionality bootstraps now reject a sparse polytomous
  replicate if it loses a category and would therefore be fitted on a
  shorter item scale. Such draws are included in the existing failure
  accounting.

- A full EFRM person bootstrap now supplies every covariance exposed in
  `unit_cov`, including the joint within-frame threshold/group-unit
  block and the set-unit/group-unit cross-block. These matrices
  therefore agree with the bootstrap standard errors and frame-unit
  tests returned by the same fit.

- Single-value display arguments no longer flatten one-cell matrices:
  [`plot_facets()`](https://drjoshmcgrane.github.io/rasch/reference/plot_facets.md)
  requires a plain facet name,
  [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  requires a plain reference name, and report titles require plain
  character scalars. The app also withholds legacy BTL dependence
  decisions when an adjusted probability is absent, rather than
  highlighting the raw probability.

- [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md)
  now guarantees that an item requested through `disordered` has an
  adjacent threshold reversal, including PCM simulations with randomly
  varying threshold spans. Generating threshold lists are also named by
  item.
  [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md)
  and
  [`simulate_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_efrm.md)
  now refuse `missing = 1`, which would remove every response and leave
  no model to fit.

- The two equating functions now refuse a fitted model from the other
  response family directly, rather than passing it to data-frame
  coercion as though it were an item or object bank.

- [`split_items()`](https://drjoshmcgrane.github.io/rasch/reference/split_items.md)
  now requires a plain grouping vector. A one-column matrix can no
  longer be flattened silently into a person-factor allocation.

- Tailored-analysis anchor items now require unique, non-missing item
  names in a plain vector. Matrix-shaped or repeated anchors are refused
  before the recalibrations begin.

- A data-driven dimensionality bootstrap now withholds its verdict when
  failed refits leave a probability floor above the requested alpha
  level. Printed, app and exported results distinguish that case from an
  uncalibrated residual split and retain the bootstrap resolution and
  replicate accounting.

- Traditional SEM is now withheld when coefficient alpha is negative.
  The tailored-bootstrap resolution warning also covers a Holm
  probability floor exactly equal to .05 and gives the first usable
  replicate count.

- Item-curve overlays now retain each fitted item’s class intervals
  under missing data, honour `observed = FALSE` with group displays, and
  reject repeated selectors. Plot palettes cycle for larger category and
  group sets, and residual components beyond the tenth can be selected
  directly.

- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
  now refuses non-data-frame input at entry with a direct message
  instead of failing later through a long-format subscript or an
  internal assertion.
  [`rack_data()`](https://drjoshmcgrane.github.io/rasch/reference/rack_data.md)
  and
  [`stack_data()`](https://drjoshmcgrane.github.io/rasch/reference/rack_data.md)
  now check source column names before coercion, so duplicate columns
  cannot be silently renamed and reshaped as different items. The
  built-in person-item and Wright targeting plots now also refuse a
  failed calibration or EFRM link instead of drawing its last numerical
  iterate as an estimate. Direct BTL scale, category, characteristic,
  dependence and frame-unit plots apply the same convergence rule. BTL
  transitivity and residual-dimension plots now give a direct class
  error when passed the wrong result object. Complete threshold,
  item-fit, person-fit, person diagnostic, fit-residual, facet and
  frame-unit plots follow the same rule. Complete exports and reports
  reject a failed calibration before creating output, rather than
  stopping later with a partial result directory. Rest-measure
  distractor analyses and plots also refuse a failed calibration.
  [`plot_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/plot_recovery.md)
  now gives a direct class error for an unrelated object.

- A paired-comparison solution rejected by the late quasi-separation
  check now withdraws its complete covariance-based inference state.
  Composite information criteria, object and graded-threshold standard
  errors, threshold components and the pairwise fit probability are no
  longer reported from a calibration that the estimator has marked as
  non-convergent. Composite information criteria for ordinary Rasch fits
  likewise require a finite positive-semidefinite Godambe covariance and
  an identified positive- definite sensitivity matrix; an absolute trace
  can no longer turn invalid ingredients into an apparently usable
  penalty. A non-converged BTL-EFRM fit now likewise retains its final
  estimates only for diagnosis and withholds standard errors and
  inferential probabilities. The same rule now applies to ordinary
  response-data and paired-comparison optimisers that stop before
  convergence: locations and residual patterns remain for diagnosis,
  while standard errors, separation reliability and probabilities are
  unavailable. The exported low-level
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  and
  [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md)
  entry points now warn on the same condition and withhold their
  covariance-based uncertainty. Explanatory fits and their
  model-comparison helper now also require the unrestricted reference
  calibration to have converged; a simpler restricted model can no
  longer be assessed against a failed reference fit.

- Rasch simulation now refuses duplicate dependence pairs, repeated item
  selectors and non-scalar response-style types. Zero-sized DIF and MFRM
  interaction requests, and zero-strength dependence, are no longer
  reported as planted departures. Zero-size EFRM item drift is treated
  in the same way. Response-style/careless and halo/erratic-rater
  proportions are refused when they cannot be realised without one
  mechanism erasing the other. Completely-at-random missing cells are
  now sampled only from responses not already removed by planted
  speededness, so both requested counts remain observable; incompatible
  proportions are refused.

- Simulation recovery now verifies that the fitted responses and
  allocations belong to the supplied simulation before comparing
  parameters. Row and item order remain immaterial. EFRM sets, person
  groups and paired-comparison panels are matched by membership, so
  renamed frames retain their recovery summaries; a different partition
  is refused.

- Fixed departures generated by the app for explanatory models are now
  orthogonal to the nominated predictor design. The requested magnitude
  can no longer be partly, or in a saturated design wholly, absorbed by
  the explanatory coefficients.

- [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  now treats row order, item-column order, arbitrary person or judge
  labels, and count compression as changes of presentation rather than
  changes of data. It still distinguishes different person or judge
  allocations, and a half-scored tie remains one comparison rather than
  two independent opposing judgements.

- Simulation recovery now retains the identified origin for externally
  anchored Rasch and paired-comparison fits.
  [`distractor_rescore()`](https://drjoshmcgrane.github.io/rasch/reference/distractor_rescore.md)
  refuses to propose scores when no respondent selected the full-credit
  option. Export functions also reject an invalid fitted object before
  creating a partial output directory.

- Person scoring, test information and model-curve plots now follow an
  externally anchored calibration origin instead of searching or drawing
  on a fixed zero-centred interval. Model curves are withheld when
  calibration or EFRM linking did not converge.

- When an ordinary Rasch analysis contains repeated person identifiers,
  its pairwise item-calibration sandwich now sums score contributions
  within person before estimating uncertainty. Point estimates still use
  every response row; duplicating an occasion no longer makes item
  standard errors spuriously smaller. This covariance is also used by
  explanatory refits and resolved DIF magnitudes. CL-BIC counts
  contributing persons rather than response occasions. Model comparisons
  now include the independent-person allocation in their data check and
  withhold composite information criteria when it differs from the
  reference fit. Unique-identifier analyses are unchanged. Item-trait
  and class-interval ANOVA probabilities are now withheld when
  identifiers repeat because their reference distributions count
  response rows rather than independent persons. The ordinary item-fit
  bootstrap is refused for the same reason; descriptive fit statistics
  and repeated-measures DIF remain available.

- The threshold diagnostic now determines exactly which partial-credit
  categories are modal over a non-zero interval. It no longer depends on
  a finite plotting grid that could miss a very narrow category region.
  Thresholds that are equal by construction now leave only the extreme
  categories modal; their crossing points differ by rounding, and that
  gap was previously reported as an interval.

- The semiparametric EFRM set link now judges grid truncation from the
  persons with a finite location in at least one set. Persons at the
  same extreme tail in both sets remain in the likelihood but no longer
  count towards the two per cent edge-mass rule, since the masses
  explaining them sit on the end of any finite grid. The rule refused
  about one person resample in twelve on an ordinary design, so a
  bootstrap
  [`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  rarely reached its usable-replicate minimum and its usable replicates
  under-represented the resamples richest in extreme scores. The rule is
  now documented.

- Paired-comparison EFRM simulation now requires plain vectors for panel
  units, set units and set origins instead of silently flattening
  matrices.

- Anchor-table columns must now contain one plain value per row.
  Matrix-valued data-frame columns can no longer be flattened or
  recycled into a different set of item or threshold anchors.
  Paired-comparison object anchors apply the same plain-vector rule.

- Externally imposed person-estimation weights and item-set maps now
  reject matrices instead of flattening them into a different weighting
  scheme. Set-map shape is checked before label canonicalisation.
  BTL-EFRM panel-unit reconciliation likewise refuses asymmetric
  within-set covariance matrices rather than silently symmetrising their
  precision weights.

- Numeric columns in external calibration banks must be vectors.
  Matrix-valued data-frame columns can no longer be flattened and
  misaligned with their item or object names.

- User-defined DIF contrast weights must be plain numeric vectors;
  shaped weights can no longer be flattened into a different contrast.

- When dichotomous paired-comparison ties are divided equally between
  the two outcomes without a judge column, the two half rows now remain
  one sampling unit in the sandwich covariance. Count-compressed ties
  retain their original replication count in that calculation.

- Item selectors in ordinary and EFRM fits, EFRM item-set maps, and
  paired-comparison margin columns now reject matrix-valued inputs
  rather than flattening them into a different fitted design or response
  scale. Paired-comparison EFRM panel and object-set maps apply the same
  rule. Multiple Ratings facet, item, interaction and factor selectors
  are likewise required to be plain vectors. By-value person identifiers
  in ordinary and EFRM fits must also be plain vectors. By-value Rasch,
  EFRM and Comparative Judgement person or judge factors and EFRM frame
  groups now apply the same rule.

- EFRM and paired-comparison EFRM omnibus Wald tests now refuse
  materially asymmetric covariances. A singular covariance is used only
  when the effect lies wholly in its estimable subspace; an effect in a
  discarded direction is no longer silently ignored. The Kent
  composite-likelihood calibration now applies the same covariance
  validity gate to its Godambe covariance and inverse sensitivity.

- Extrapolated extreme-score person locations are now reported only when
  the adjacent finite-score estimates support a stable geometric
  continuation. Otherwise
  [`score_table()`](https://drjoshmcgrane.github.io/rasch/reference/score_table.md)
  refuses the extrapolation and
  [`person_extrapolated()`](https://drjoshmcgrane.github.io/rasch/reference/person_extrapolated.md)
  retains the fitted Warm estimate.

- Uncertainty-weighted adaptive paired-comparison recommendations now
  require an aligned, finite, symmetric positive-semidefinite object
  covariance (or finite object standard errors when no covariance is
  stored). Information- only recommendations remain available with
  `weight_se = FALSE`.

- Covariance gates and public scalar controls now reject malformed
  shaped, classed or complex inputs consistently. Manual dimensionality
  subsets and DIF design selectors also refuse duplicate, overlapping or
  ambiguous specifications before fitting.

- Fit-derived person, targeting, information, scalogram and Wright-map
  outputs now verify EFRM link convergence from the model’s per-edge
  status as well as its overall convergence flag. This also protects
  restored or modified fits whose two records do not agree.

- Tailored-analysis bootstrap runs can take an explicit seed without
  changing the caller’s random-number state. The app now validates its
  replicate count and seed and includes that seed in the displayed R
  code.

- Judge-level unexpected-object and unexpected-pair flags now use
  Holm-adjusted approximate normal probabilities over the eligible
  family. Their earlier per-row z cut-off did not control familywise
  error.

- Paired-comparison DIF now retains anchored, failed or otherwise
  unavailable planned follow-up contrasts in the Holm family size. A
  failed resolution can no longer make the surviving adjusted
  probabilities smaller.

- Response-dependence magnitudes, resolved DIF contrasts, both equating
  procedures and conditional frame-invariance comparisons now refuse or
  withhold Wald inference when their fitted covariance is unavailable,
  asymmetric or materially indefinite. Point estimates remain available
  for description where the comparison itself remains defined. Multiple
  Ratings DIF also withholds its covariance-weighted pooled point
  estimate in this case, because the invalid covariance would otherwise
  determine the pooling weights themselves.

- EFRM corrected-moment set linking again runs its requested bootstrap.
  Its non-iterative edges record no optimiser status; those structural
  missing values were incorrectly treated as non-convergence, leaving
  unit standard errors unavailable.

- Crossed EFRM group-unit decompositions now withhold their GLS table
  when either the unit covariance or the derived coefficient covariance
  is asymmetric or materially indefinite.

- Residual PCA plots and biplots now use a finite neutral plotting range
  when their displayed loading vectors are identically zero.

- [`plot_facets()`](https://drjoshmcgrane.github.io/rasch/reference/plot_facets.md)
  now refuses non-character facet selectors rather than allowing a
  factor code to be interpreted as a list position. Item selectors
  likewise refuse matrices and other shaped objects rather than treating
  a one-cell object as a scalar item index. Multiple-choice item subsets
  also reject duplicate canonical names. Batch item and person plot
  exports now validate their complete selection before opening an output
  file; fractional row numbers can no longer be truncated to a different
  person.

- The application DIF summary now distinguishes unavailable uniform and
  non-uniform tests from non-significant adjusted results.

- Available-case coefficient alpha now requires the pairwise covariance
  to be positive semidefinite at a scale-relative numerical tolerance.
  The former absolute tolerance could accept materially indefinite
  binary-item covariance matrices.

- Hybrid EFRM set linking no longer falls back silently to a person-only
  bootstrap when the joint stage-one covariance cannot be used.
  Requested unit inference is refused because that fallback omits
  calibration uncertainty; a descriptive fit remains available with
  `boot_reps = 0`.

- Equating now validates bank covariance symmetry, positive
  semidefiniteness and agreement with stated marginal standard errors at
  the covariance’s own numerical scale. Precision-weighted links and
  pooled Multiple Ratings DIF magnitudes no longer use an absolute
  variance floor.

- Whitespace-only DIF levels are treated as missing factor metadata
  rather than as a person or judge group. EFRM now fits the same trimmed
  frame-group labels that it validates, so visually identical labels
  cannot define different frame units. Whitespace-only data, item and
  factor names are refused.

- Custom DIF contrast weights are normalised by their direction rather
  than their absolute magnitude; rescaling every weight cannot make a
  valid contrast disappear.

- Structural item changes now reject missing and blank item names, and
  direct DIF splits treat blank group labels as missing. Externally
  weighted person estimates identify extreme response patterns from
  their responses rather than an absolute weighted-score tolerance, and
  no longer merge distinct weighted totals for computational caching.

- Numeric `NaN` values supplied as BTL judge-group labels, including
  grouped ICC displays, are now treated as missing rather than as a
  literal `"NaN"` group. Missing and whitespace-only BTL anchor names
  are also refused explicitly.

- Externally weighted person estimates now normalise relative weights
  without underflow when all positive weights are close to the smallest
  finite double. Exact fitted item names take precedence when item
  weights or item-set maps are matched, so literal leading or trailing
  spaces in a data-column name are not discarded; padded selector input
  remains available when unambiguous. Native EFRM person estimates
  likewise identify extreme responses from the response pattern rather
  than an absolute weighted-score tolerance, and no longer merge
  distinct weighted scores for computational caching.

- MFRM marginal fit summaries return an unavailable value when every
  contributing response-cell residual is unavailable, rather than
  emitting a warning and storing `NaN`. MFRM simulation recovery now
  identifies the planted rater facet from the strongest level-label
  match, using its name only to resolve a tie, instead of assuming that
  it is the first fitted facet.
  [`sim_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/sim_recovery.md)
  also refuses a non-convergent fit or a model family that does not
  match the generating simulation, rather than comparing coincident
  parameter labels. Replicated simulations also retain a valid seed at
  the upper integer boundary instead of overflowing during seed
  incrementation.

- The paired-comparison fit bootstrap now gives a deliberate refusal
  when a compressed comparison count exceeds the integer range accepted
  by R’s multinomial generator.

- In the app, planned Multiple Ratings DIF contrasts now select
  underlying items rather than internal item-by-facet response cells.
  Saved-analysis notices now also name the superseded derived result
  actually omitted, including frame-invariance results, rather than
  always referring to a fit bootstrap.

- Parallel bootstrap requests from an installed package now start the
  requested socket workers. An incorrect namespace-name check had made
  EFRM, BTL–EFRM, DIF, dimensionality and fit bootstraps run serially
  while still recording the requested worker count. The application now
  refuses fractional or unavailable worker selections instead of
  truncating them silently.

- Joint bootstrap references now require every transformed statistic or
  fitted parameter in a replicate to be finite. Infinite log
  mean-squares, frame units, or replicated eigenvalues are treated as
  failed draws rather than entering maximum-statistic or covariance
  calculations. EFRM NPML compression also retains the full precision of
  weighted sufficient statistics instead of pooling values that agree
  only to 12 significant digits.

- Kent-adjusted model comparisons and unit-family Wald tests now
  withhold inference when their estimated covariance is singular in a
  tested direction or materially indefinite. BTL–EFRM panel-unit
  reconciliation no longer adds a fixed variance floor to its precision
  weights. Test-information curves are withheld for an unconverged
  calibration, and zero-degree-of-freedom model summaries no longer
  display an infinite or undefined chi-square ratio. Targeting summaries
  and WrightMap displays now likewise refuse a non-convergent
  calibration. Item-trait detail retains an interval with unavailable or
  zero model variance but does not manufacture an infinite residual from
  it. Fully fixed item locations are excluded from item separation
  reliability, rather than being counted as estimates with zero standard
  error.

- PSI bands are now described as person separation quality rather than
  as the power of fit tests. The old `power_of_fit` component remains as
  an alias for saved-object compatibility. EFRM documentation now states
  the identification conventions for ordinary and paired-comparison
  frame models exactly.

- Person-subset dimensionality tests and simulated residual references
  now refuse repeated person identifiers. Their current null generators
  treat response rows as independent and cannot preserve within-person
  dependence. The fixed-split interval is described as Clopper–Pearson
  rather than as an exact test of dimensionality.

- [`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  with bootstrap uncertainty now requires the standard 90% usable-refit
  rate and reports usable, non-converged and other-failure counts. A
  contrast with zero bootstrap uncertainty is retained descriptively,
  but its Wald probability is withheld. Every observed frame in a
  compared item set must also yield a usable separate calibration;
  unavailable frames are no longer silently omitted. An item dropped or
  rescored by one separate calibration is recorded as an unavailable
  comparison and retained in the multiplicity family. The procedure now
  gives a clear refusal when no frame pair retains two items with a
  comparable category structure.

- Schema-2 application projects made with the earlier text signature now
  reach the intended maximum-statistic-result migration. Multi-design
  test information plots show standard-error-of-measurement curves as
  well as information curves.

- Wald follow-up tables now withhold probabilities when a contrast has
  zero estimated uncertainty. Positive finite standard errors remain
  testable, so rescaling a predictor or contrast cannot change
  inference. This applies consistently to DIF, frame units, explanatory
  coefficients, superitem spread and both ordinary and Comparative
  Judgement equating.

- Missing numeric identifiers, including `NaN`, are no longer converted
  into literal factor levels. Unknown person identifiers remain separate
  in DIF and tailored-bootstrap calculations, and malformed `NaN` anchor
  indices are refused.

- A Comparative Judgement anchor whose object has no usable comparison
  after data preparation is now refused instead of being silently
  removed from the fitted anchor set.

- Paired-comparison equating now refuses malformed bank `df_location`
  and polytomous score-scale attributes rather than silently replacing
  the former with an asymptotic reference or coercing the latter from
  another class.

- Saved schema-2 application analyses carrying frame-invariance results
  from the earlier comparison rules now reopen after that derived result
  is authenticated and omitted; the application asks the analyst to
  recompute it.

- A tailored-analysis bootstrap that cannot retain enough refits now
  signals the standard structured refusal, including requested, usable,
  non-converged and other-failure counts.

- [`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  now retains excluded comparisons in its Holm family. Bootstrap
  discrimination probabilities are withheld when either fitted slope is
  on its imposed boundary; the ratio remains descriptive.

- [`dimensionality_test()`](https://drjoshmcgrane.github.io/rasch/reference/dimensionality_test.md)
  documents that its default split, chosen from the residuals, makes the
  ordinary fixed-split binomial rule anti-conservative. The earlier
  documentation claimed the procedure held an exact null. A new `B`
  argument calibrates the split by a parametric bootstrap under the
  fitted model (score-conditional replicates, each refitted and split
  afresh from its own residuals), giving `p_boot`, the mean replicate
  proportion `prop_null`, and a verdict at the chosen `alpha` level. An
  automatic split without bootstrap calibration is now descriptive and
  has no binary verdict. The application can run the bootstrap and
  carries the exact item split and result into saved analyses, tables
  and reports. Refactored
  [`residual_pca()`](https://drjoshmcgrane.github.io/rasch/reference/residual_pca.md)
  and the fit-bootstrap refit so the replicates travel the observed
  computation.

- A category observed only in extreme response patterns (every answered
  item at its floor, or every one at its ceiling) carries no pairwise
  conditional information, and its partial-credit threshold diverged
  during estimation until the projected information matrix was reported
  singular.
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  now merges such a category with its neighbour during data preparation,
  with a note, and drops an item with no conditional information at all.
  An item with a free threshold in that position is refused when it is
  anchored; an item whose every threshold is anchored keeps its coding,
  since nothing is estimated for it.
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  names the item and category instead of reporting the singular matrix.

- Tailored analysis step 3 now uses average item anchoring (RUMM2030’s
  convention): the retained items are constrained to the mean location
  of the step-2 calibration rather than each being fixed individually,
  so the step-3 standard errors reflect the uncertainty of the
  re-estimated locations.
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  and
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  accept `anchors` with `average = TRUE`, and the application offers the
  option. Anchor documentation now distinguishes a change of origin from
  the extra relative constraints imposed by several fixed parameters,
  and states that standard errors are conditional on the supplied anchor
  values.
  [`drop_items()`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
  refuses to remove an externally anchored item because that would
  silently change the fitted scale.

- [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md),
  [`report_document()`](https://drjoshmcgrane.github.io/rasch/reference/report_document.md)
  and
  [`save_outputs()`](https://drjoshmcgrane.github.io/rasch/reference/save_outputs.md)
  now say why a scree plot is absent: a residual-PCA refusal or a
  reference that could not be completed is warned about and, in the
  reports, noted in place of the plot, where the precomputed-scree
  change had dropped it without a word. The document template wraps a
  long refusal so it stays inside the figure.
  [`spread_test()`](https://drjoshmcgrane.github.io/rasch/reference/spread_test.md)
  documents that Andrich’s tabulated bounds are exact for two- and
  three-item subtests and conservative for larger ones.
  [`distractor_analysis()`](https://drjoshmcgrane.github.io/rasch/reference/distractor_analysis.md)
  no longer warns when nobody chose the keyed option.

- Bootstrap maximum-statistic probabilities now standardise each
  simulated row against the other rows. This gives the simulated and
  observed statistics the same external standardisation and removes the
  anti-conservative shrinkage caused by using a null row in its own mean
  and standard deviation. Constant statistics remain in the declared
  family. Saved projects omit fit-bootstrap results made with the
  earlier adjustment and ask the analyst to recompute them.

- Residual scree references and tailored-analysis bootstraps now require
  the same usable-refit proportion as the fitted-model bootstrap: at
  least 90% for 30 or more requested refits. Their results report
  requested, usable, non-converged and other-failure counts. Scree and
  paired-comparison dimensionality functions accept a seed and restore
  the caller’s random- number state.

- The application keeps the computed dimensionality and EFRM invariance
  analyses with the saved project and passes them to reports and
  complete exports. The selected EFRM uncertainty method is therefore
  not replaced by a conditional rerun. Complete exports and editable
  reports now include the dimensionality tables and figures for Rasch
  and Comparative Judgement analyses.

- Simulated response styles now preserve declared local dependence by
  rebuilding dependent responses from the styled source responses, and
  their log-probability tilt remains stable at large finite strengths.
  [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) now
  refuses count-compressed rows with judgment order: the row does not
  retain the sequence needed to estimate exposure and carry-over
  effects. The dimensionality diagnostic retains the same guard for
  older saved fits.

- Scree parallel analysis now draws from the exact Rasch distribution
  conditional on each person’s observed score and missingness pattern.
  Its shaded reference band runs from the simulated mean to a
  finite-simulation familywise 5% upper critical curve. Each simulated
  maximum is standardised against the other draws, so its reference
  matches the externally standardised observed value. The returned table
  gives marginal and adjusted upper-tail probabilities for every
  displayed component. Saved result and application-project signatures
  use platform-independent binary serialisation; analysis files carrying
  the earlier text signatures remain valid across Windows and Unix line
  endings.

- Crossed Extended Frames groups retain their original column names in
  fitted objects and factorial unit tables. Frame-defining groups are
  therefore excluded correctly from a subsequent DIF analysis even when
  their names contain spaces, reserved words or punctuation.

- Explanatory models in the application now preserve exact predictor
  names and distinguish literal colons from interactions. Ordered
  predictors keep unique adjacent-contrast names. Comparative Judgement
  DIF display settings are kept outside the signed result, and saved
  primary DIF analyses are checked against the active fit whether or not
  a bootstrap was requested.

- Fit-bootstrap and primary DIF results now carry their own integrity
  fingerprints. Before an analysis is restored, reported or exported,
  the package checks the fitted-model identity, replicate accounting,
  declared family, adjusted probabilities and the agreement between
  detailed and summary tables. Edited or incomplete results are refused.

- DIF bootstrap results are now checked before they are restored,
  reported or exported. The replicate matrices, accounting totals,
  declared test family, primary DIF tables and derived probabilities
  must agree; an incomplete or inconsistent result is refused before any
  output files are written.

- [`dif_bootstrap()`](https://drjoshmcgrane.github.io/rasch/reference/dif_bootstrap.md)
  adds an optional sensitivity analysis for DIF. It conditions on person
  scores for Rasch, explanatory and Multiple Ratings models, on item-set
  subtotals for Extended Frames, and on the fitted outcome model for
  Comparative Judgement. Each replicate refits the model and repeats the
  declared DIF analysis. Adjusted probabilities use a single-step
  minimum-p reference over the complete family;
  [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  or
  [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  remains primary. They describe the fitted global invariant null, not
  strong familywise control when another member has DIF. The application
  runs the analysis in the background and retains it with its primary
  DIF specification in saved projects, reports and exports. Validation
  now separates power for the affected member from familywise error
  among invariant members under partial alternatives in every supported
  model family.

- Fitted-model provenance now includes the complete fit, including
  person identities, comparative-judgement thresholds and
  sequential-dependence estimates. Application project files use a
  sealed schema that detects a changed dataset, fit, history or saved
  result before restoration, and structurally validates every stored fit
  and history entry. Schema-1 files are structurally checked and
  upgraded when opened; saved bootstrap or DIF results that predate the
  current fingerprints are omitted with a warning, while the data, fit
  and history are retained. Saving the analysis again records the new
  integrity information.

- Fit bootstraps now withhold inference when fewer than 90% of 30 or
  more requested refits are usable. This prevents sparse polytomous
  analyses from reporting a null distribution selectively thinned by
  failed category structures. Conditional generation is evaluated on the
  log scale, so finite extreme thresholds no longer overflow. Printed
  bootstrap objects give a short summary rather than dumping the stored
  replicate matrices.

- DIF results carry the identity of their fitted model and cannot be
  exported or reported with a different calibration. Saved application
  projects now validate every fitted-model history entry, not only the
  active one. A saved DIF bootstrap must also match its primary DIF
  specification and restored controls. Reports include the corresponding
  DIF magnitude estimates.

- Adjusted probabilities for DIF magnitudes and contrasts,
  paired-comparison DIF magnitudes, explanatory coefficients and
  diagnostics, MFRM interaction follow-ups, and tailored analysis now
  retain the complete declared family when one probability is
  unavailable.

- The paired-comparison application example separates sequential
  dependence from residual dimensionality. Its main dataset has planted
  exposure, null carry-over and null judge factors; a second dataset
  supplies the dimensionality example.

- Paired-comparison fitted-design bootstraps with exposure and
  carry-over have been calibrated for dichotomous and four-category
  responses at 30 judges. Total and adjusted pair, object and judge
  error remained near 5%; all 238,800 refits were usable.

- Fit inference now retains its predeclared family when an item or
  frame-unit probability is unavailable. BTL-EFRM omnibus tests are
  withheld rather than reduced when a requested covariance coordinate is
  missing. Fixed-location item bootstraps keep each location with its
  own response row, and whole-test fit-residual nulls use one fixed item
  set. Saved application bootstraps are checked against the active
  restored fit. Parallel bootstraps loaded from a development source
  tree run serially instead of finding another installed package on
  socket workers.

- [`fit_bootstrap()`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  now calibrates person fit under its default score-conditional null. It
  reports marginal probabilities and a maximum-statistic adjustment
  across persons for each fit statistic. Paired-comparison fits are also
  supported: outcomes are generated on the fitted comparison design, the
  model is refitted, and bootstrap probabilities are returned for the
  total pairwise chi-square and pair, object and judge fit. Ordered
  thresholds and sequential history effects are reproduced. The
  application runs either fit bootstrap in a cancellable background
  process and carries the result into its tables, saved analyses and
  reports.

- Fit bootstraps now refuse a non-converged observed fit, reject
  non-converged refits, and report non-convergence separately from other
  failed replicates. A maximum-statistic adjustment is withheld for the
  complete predeclared family when any member lacks a usable joint null.
  Bootstrap results carry their model identity and can only be exported
  with the fit that produced them. The application shows unavailable
  inference neutrally. Documentation now states that adjustment is
  within each statistic under the fitted global null; it is not strong
  familywise control for otherwise fitting members when another member
  misfits.

- Bootstrap accounting now follows the quantities actually used. EFRM
  and BTL-EFRM size their rank requirement from the free directions in
  the largest covariance block, rather than the concatenated row
  retained from each replicate. Item-fit bootstrap probabilities are
  withheld when too few refits or too few values of a particular
  statistic remain; the item table records the usable count for every
  statistic, and Holm adjustment retains the full family when one result
  is unavailable. Seeds must be non-negative whole numbers. The
  asymptotic item-fit probabilities remain available as descriptive
  diagnostics but are labelled as approximate in summaries, reports and
  the application. A calibrated bootstrap probability is the operative
  result after
  [`fit_bootstrap()`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  is run.

- A person-item map now refuses an EFRM item-by-frame cell and
  person-group selection with no common response cell. It no longer
  replaces the requested item silently with all items from the selected
  group.

- Release review, second pass.
  [`fit_bootstrap()`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  warns when the requested replicates cannot reach Holm-adjusted
  significance at .05 for the item count — the floor is 2L/(B + 1)
  two-sided, so a run whose adjusted probabilities could never flag
  anything says so before it spends the time. A restricted person-item
  map carries the information of its own selection:
  [`test_information()`](https://drjoshmcgrane.github.io/rasch/reference/test_information.md)
  gains an `items` argument and
  [`plot_pimap()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pimap.md)
  passes its item subset through, where the curve previously described
  the whole instrument over a subset’s distributions. A person group is
  addressed unambiguously:
  [`plot_pimap()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pimap.md)
  accepts the qualified `"factor: level"` form and refuses a bare level
  that two factors share, naming the qualified candidates, where the
  first factor was silently taken; the application’s selector now offers
  the qualified form throughout, and its generated code quotes names by
  [`deparse()`](https://rdrr.io/r/base/deparse.html) so a level carrying
  a quote still yields runnable code. A saved analysis carries its
  bootstrap: the project stores and restores
  [`fit_bootstrap()`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  results, and reports and complete exports render the analysis as run —
  the application’s configured DIF model and any bootstrap null — rather
  than recomputing defaults
  ([`save_outputs()`](https://drjoshmcgrane.github.io/rasch/reference/save_outputs.md),
  [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md)
  and
  [`report_document()`](https://drjoshmcgrane.github.io/rasch/reference/report_document.md)
  gain `dif` and `bootstrap` arguments), and every bootstrap probability
  a report renders is the Holm-adjusted one. Bootstrap covariance
  acceptance is dimension-aware: successes must also exceed the free
  directions in the largest covariance block used by the fit, since the
  stored row contains several separate constrained blocks. The alpha–phi
  cross-covariance requires 30 joint draws and a majority, and is
  withheld with a warning below that. On an extended-frame map
  restricted to one group, the information curves are that group’s
  designs alone, and
  [`test_information()`](https://drjoshmcgrane.github.io/rasch/reference/test_information.md)
  refuses fractional indices and accepts an extended-frame fit’s
  underlying item names.

- [`fit_bootstrap()`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  generates its replicates conditionally by default: each person’s
  responses are drawn from the Rasch conditional distribution given
  their observed raw score over their own observed items, with
  sufficiency cancelling the person parameter. Release review caught two
  defects in the ability-resampling default this replaces. Resampled
  estimates carry their estimation error, and the standardised
  statistics feel the inflated spread as the sample grows — infit z
  rejected 14.4% of correctly fitting items at 4,000 persons, against
  2.1% conditionally, with the chi-square at 5.2% and the fit residual
  at 3.7%. And abilities drawn independently of each person’s
  missingness sever any tie between who answers and what is missing: a
  linked-booklet design’s observed 1.29-logit group difference fell to
  0.12 in resampled replicates and is reproduced at 1.28 conditionally,
  exactly, because the scores that carry it are held; the recorded
  booklet scenario holds the item-level rate at 1.9% and Holm familywise
  error at 3% for the chi-square and 6% for the fit residual. The
  ability-sampling schemes remain by name. Review also tightened the
  arithmetic around the resolution floor — a two-sided bootstrap
  probability cannot fall below 2/(B + 1), so Holm across L items needs
  B of at least 40L to flag at .05, which the documentation, the
  application’s guidance and its adjusted fit-residual display now state
  and follow — and two operational faults: the shared bootstrap harness
  now degrades to a serial run with a warning when its sockets cannot be
  opened (every draw is made before dispatch, so the result is
  identical), where it previously errored inside documented examples;
  and bootstrap acceptance again honours the documented contract of at
  least 30 successful replicates and a majority of those requested,
  where a bare majority — 16 of 30 — had been allowed to price a
  covariance.

- [`fit_bootstrap()`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  gives the item fit statistics a null distribution they can be referred
  to. Every one of them is computed at estimated person locations and
  referred to a distribution derived as though those locations were
  known, and each is miscalibrated by an amount that moves with the
  sample. The item-trait chi-square forms its class intervals by
  selecting on the estimates, so the expected interval means carry a
  bias that does not shrink while the intervals grow: at eight items a
  correctly fitting item is rejected 8.9% of the time at 250 persons and
  99.9% at 4,000, and the whole-test total rejects every correctly
  fitting instrument from 1,000 persons on. The fit residual fails the
  other way, its null SD running from 0.74 at 250 persons to 1.00 at
  4,000 with a mean drifting to -0.36, so the conventional 2.5 cut flags
  under 1% of correctly fitting items and means something different at
  each sample size. Infit z beyond 1.96 flags 11.5% at 250 persons and
  68.4% at 4,000. The bootstrap generates from the fitted model and
  sends every replicate down the same road as the observed data, so the
  same bias enters the null and cancels: under the score-conditional
  default, item-level error across all five sample sizes is 0.039-0.051
  for the chi-square, 0.039-0.054 for the fit residual, 0.033-0.050 for
  infit z and 0.037-0.055 for outfit z, and the whole-test total holds
  at 0.020-0.080. Power is retained in full — an item generated at
  discrimination 2.5 is detected 1.000 of the time at 2,000 persons
  while the clean items beside it come back from 0.723 to 0.037. One set
  of replicates serves every statistic, since a replicate must refit the
  whole model in any case, and the whole-test readings include the mean
  and SD of the item fit residuals: the convention reads that SD against
  1, and it is nowhere near 1 until the sample runs to a few thousand.

- `adjust_N` is removed from
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
  [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
  [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  and
  [`rasch_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_explanatory.md),
  and from the application’s run settings. Rescaling every item’s
  chi-square to a reference sample size applies one constant across
  items, so it can move how many items are flagged but not which: at
  2,000 persons its power on a planted 2.5-slope item equalled its error
  rate on the fitting items beside it at every reference value tried,
  and at a reference of 500 it silenced that item in all 120 replicates.
  The justification offered for it — that a fitting item’s chi-square is
  unchanged by sample size while a misfitting item’s grows — does not
  hold: under a model true by construction the mean chi-square runs 4.0,
  7.9, 11.5, 18.2 and 30.7 for 250 to 4,000 persons on fixed degrees of
  freedom.
  [`chisq_detail()`](https://drjoshmcgrane.github.io/rasch/reference/chisq_detail.md)
  loses `adjust_factor`, `chisq_unadjusted` and the `chisq_adjusted`
  interval column with it; its intervals now sum to the item’s
  chi-square directly. Saved application projects carrying the setting
  still open, the value being ignored.

- The application’s Summary panel renders after a comparative judgement
  fit. Its headline tiles and test-of-fit table sat inside the panel
  shown only for the other models, where their own
  show-for-comparative-judgement condition contradicted the ancestor’s
  and the panel stayed empty; the override history moves out of that
  panel with them, matching the object estimates it already reported.

- [`plot_btl_dependence()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dependence.md)
  includes a fitted first-position advantage in its partial-residual
  baseline. The baseline summed only the exposure and carry-over terms,
  and the expected-score map is nonlinear, so on a fit made with
  `position = TRUE` every observed and fitted departure in the plot was
  displaced; the estimates were never affected. Found by review; the
  regression test pins the corrected baseline’s partial residuals at
  zero under the fitted model.

- The application’s Items page can bootstrap the fit statistics on
  request: one button runs
  [`fit_bootstrap()`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  on the active fit and swaps the chi-square and fit-residual
  probabilities in the table, its CSV, the headline tile and the note
  for the calibrated ones, labelled as such and cleared by the next
  estimation. The default 999 replicates clears the familywise
  resolution floor — the smallest Holm-adjusted probability is the item
  count over B + 1, so flagging at .05 needs at least 20 replicates per
  item, which the control’s help states.

- The application’s Wright map controls sit above the plot. They
  rendered in the card footer, where choosing the WrightMap renderer
  revealed its panel and display controls two screen-heights down; a
  control that changes what is drawn belongs where the change is made.
  The data preview also shows 25 rows rather than 10 — the sidebar it
  sits beside is longer than ten rows ever were.

- The application’s panels start closed, except where the accordion is
  the page: Local dependence, Trait dimensionality and DIF open their
  first panel, since everything those pages show lives inside one.
  Elsewhere twelve accordions opened a panel on arrival — ten naming one
  and two relying on the default, which opens the first — so every page
  presented a choice already made. `More` also loses the only route to
  reopening a saved analysis: **Open saved analysis** now sits beside
  **Upload data** on the welcome screen, since reopening one is a way to
  start.

- The workflow and comparative judgement vignettes demonstrate the
  application over the analysis they set out in code, rather than a
  first panel or two: the workflow vignette walks Data, Summary, Items,
  the class-interval detail, Persons, Targeting, local dependence, DIF,
  the R code disclosure and export, and the comparative judgement
  vignette gains the section it never had, covering what changes when
  the roles are two objects and an observed preference. The data
  structures vignette shows an explanatory design at both levels with
  the same predictors, so the difference between an item-level and a
  threshold-level file is visible as the shape it is, with the number of
  rows a threshold design needs derived from the data and the three
  errors that refuse one named.

- [`plot_pimap()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pimap.md)
  can restrict either side of the map: `group` takes one level of a
  fitted person factor, and `items` takes item names or one item-set
  name of an extended-frame fit, whose virtual item-by-group cells are
  matched through their underlying items. The selection is named in the
  legend, so a restricted map cannot be read as the whole instrument.
  The application offers both as person-group and item-set controls
  beside the person-item map.

- The application colours its headline tiles by the side of the model
  they describe: persons blue and items amber, the colours the Wright
  and person-item maps already use. The data page follows, with rows
  blue and columns amber. Tiles that report a check keep their own
  reading, now green when it passes and red when it wants attention:
  item and person misfit, disordered thresholds, separation, coefficient
  alpha, the item-trait probability and the explanatory comparison.
  Amber is left to mean the item side alone, so the two readings cannot
  be confused.

- The application’s comparative judgement example is dichotomous. It
  also carried an ordered preference and a margin, which invited a
  polytomous fit on its ten judges; a polytomous fit with the
  judging-order covariates estimates as many parameters as there are
  judge clusters, so the example withheld its own cluster-robust
  inference. The ordered response is demonstrated on designs that can
  carry it.

- The application’s post-estimation notification lists its notes. They
  were joined on a newline, which the notification renders as ordinary
  whitespace, so a fit reporting more than one note read as a single
  run-on sentence – withheld cluster-robust inference and a withheld
  carry-over probability ran together into one unreadable line.

- A judge factor named by judge may name judges the fit set aside – one
  who only ever tied, or whose rows were dropped – since a map built
  from the source data necessarily does. Those entries are ignored and
  reported in the notes rather than refused, in
  [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  and
  [`plot_btl_icc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_icc.md)
  alike; every fitted judge must still carry an entry, which is what
  catches a mistyped judge name.

- Bootstrap covariance estimation requires at least 30 usable replicates
  and a majority of those requested. The absolute floor prevents a
  covariance from being reported from 16 draws when `boot_reps = 30`;
  the rank rule may raise the minimum for a larger covariance block.

- [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md)
  now labels unnamed response matrices consistently, and direct
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  calls reject an empty anchor table. Available-case item-rest
  correlations exclude respondents with no observed rest score.

- Frame-invariance bootstraps keep singleton person-group strata in
  their original group. Paired-comparison order effects are refused when
  the comparison design confounds them with the object locations.

- Automatic DIF follow-ups now use the adjustment method and
  significance level supplied to
  [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md).
  [`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
  also checks item selectors before fitting the contrast family.

- A saved app analysis now restores its data roles, estimation controls
  and embedded supporting data, so it can be re-estimated after it is
  reopened. Simulation recovery covers paired-comparison Extended
  Frames, including object locations, panel and set units, and set
  origins.

- Shiny background fits are tied to the data and analysis that launched
  them. A completed EFRM or paired-comparison frame fit is discarded if
  that context has changed, and opening a saved analysis cancels work
  still in progress. A new Comparative Judgement fit clears earlier
  requested DIF and frame results, while reopening a saved fit retains
  the results stored with it.

- The Shiny application can simulate ordinary and explanatory Rasch,
  Comparative Judgement, Multiple Ratings, Extended Frames and paired-
  comparison Extended Frames data. Model parameters can be varied and
  model-specific departures planted.
  [`simulate_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_efrm.md)
  now supports item drift, careless response and missingness, while
  [`simulate_btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl_efrm.md)
  can plant erratic judges. A positive planted proportion affects at
  least one observation whenever the requested departures can coexist;
  incompatible mixtures are diagnosed. Rasch simulation truth records
  speeded persons and the cells selected for missingness.

- [`weighted_person_estimates()`](https://drjoshmcgrane.github.io/rasch/reference/weighted_person_estimates.md)
  provides supplementary person measures from externally imposed item or
  item-set weights. It leaves the fitted calibration and the ordinary
  measures used for fit, reliability, targeting and DIF unchanged. Its
  Warm correction and sandwich standard error now both use the
  variability of the weighted score; the earlier correction treated
  weights as replicated observations. A new first-listed vignette gives
  the data structure required by every model family.

- Simulation truth now retains person identifiers, so recovery matches
  person estimates by ID after rows are reordered. A planted second
  trait has its requested realised correlation. Paired-comparison
  simulators balance comparisons across the declared judges and refuse
  designs with too few comparisons. MFRM interaction probabilities
  require adequate support in every item-by-facet cell rather than only
  at the pooled facet level.

- EFRM convergence now covers the fitted nonparametric masses as well as
  the set transformation. Linking uncertainty requires at least 30
  usable resamples and more than half of those requested; the requested,
  usable and failed counts are retained in the fit. The full-bootstrap
  attempt remains recorded when too few refits succeed and hybrid
  standard errors are returned. BTL-EFRM applies the same
  usable-resample rule and accounting. Fixed-iteration EFRM linking no
  longer recalculates an unused likelihood at every intermediate mass
  update.

- Automatic DIF resolution now acts on the adjusted omnibus result.
  Pairwise contrasts describe the location of a multifactor effect but
  do not impose a second significance test. Thin or incompatible factor
  cells are still refused. A prior manual split now counts as one
  resolved source item rather than several anchors, and its provenance
  survives later item dropping or subtest formation. The reported
  residual count is the number of source items, not the number of
  item-by-factor terms.

- BTL-DIF retains the validated Welch reference for two-cell contrasts
  and uses the least-supported effective-judge cell as a conservative
  reference for contrasts spanning more cells. In 500 balanced four-cell
  null fits, the new reference rejected 3.6%, against 6.2% for the
  superseded pooled-count rule. BTL-EFRM no longer counts a set unit
  fixed after an identification failure as an estimated parameter.

- The application uses Holm-adjusted probabilities for its DIF, frame
  invariance, threshold-spread, explanatory and supplementary item-fit
  decisions. BTL judge-group factors that vary within judge are refused
  rather than reduced to the first comparison row.

- Subtests must be formed before DIF splitting. A group-specific split
  copy cannot be combined because it does not provide a common item
  across groups.

- BTL-EFRM now refits every object set after its panel units are
  reconciled. The reported object locations, expected probabilities,
  composite likelihood and equal-unit comparison therefore come from the
  same fitted parameters. Previously, stable sets retained their
  independently optimised locations while their probabilities were
  evaluated at the reconciled units.

- Multifactor DIF magnitudes now match the adjusted omnibus estimand.
  Ordinary and paired-comparison main effects average complete factor
  cells equally over nuisance factors; interactions use differences
  between differences. DIF contrasts are withheld when the compared
  groups do not share an observed score structure. ETS classifications
  use the adjusted probabilities for both significance and departure
  beyond category A.

- Structural refits preserve score categories. Subtests, DIF splits and
  EFRM frame resolutions are refused when a required score category is
  absent, rather than allowing ordinary data preparation to renumber the
  scores. Resampling refits treat category loss as an unsuccessful
  replicate. Supplied app anchors and scoring keys now fail closed when
  malformed, inapplicable or unmatched.

- Crossed EFRM factorial tests and BTL order-effect tests now retain raw
  probabilities but use Holm-adjusted probabilities for decisions. In
  fresh null simulations, crossed-EFRM familywise rejection was 5.55%
  over 2,000 fits (5.93% pooled over 3,000); the three marginal rates
  were 5.25–5.75%. BTL familywise rejection was 5.9% over 1,000 fits
  with position, exposure and carry-over fitted together.

- [`simulate_btl()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl.md)
  now constructs a second object attribute with the requested realised
  correlation, rather than obtaining that correlation only in
  expectation. The largest error over 1,200 short and long object sets
  was 3.9e-16; dimensionality power was 87% in the re-run of the strong
  design. Structured simulator options are checked before their
  components are used.

- DIF splitting omits factor levels for which an item has fewer than two
  observed categories, and refuses a split unless at least two levels
  remain. EFRM refits and wide MFRM conversion use collision-free
  internal names; named BTL-EFRM panel maps take precedence over
  data-column names. Numeric `NaN` responses are refused rather than
  treated as ordinary missing data.

- EFRM and BTL-EFRM retain raw unit-test probabilities but use
  Holm-adjusted probabilities for decisions. Omnibus tests are adjusted
  across the reported unit families. BTL-EFRM follow-up contrasts form
  one family across panel units, set units and set origins, rather than
  three separately adjusted tables. The application uses the same
  adjusted probabilities. In null simulations, EFRM omnibus-family
  rejection was 3.5% and follow-up-family rejection was 1.6% among 489
  analysed fits. After the reconciled-panel BTL-EFRM refit, the
  corresponding rates were 3.9% and 3.0% over 1,000 null fits in the
  six-judge-per-panel caution design. A 500-fit supported-design top-up
  gave 4.6% raw set-unit rejection and 0.934 interval coverage.

- Character role arguments are resolved from matching column names
  rather than their length. Named judge maps are matched by judge even
  when their length happens to equal the number of comparisons. Missing
  DIF identifiers remain as separate analysis units, conflicting
  external factor columns are refused, and scores outside R’s integer
  range no longer become missing on coercion. Paired-comparison equating
  now permits its documented descriptive two-object link while
  continuing to require three objects for drift tests.

- Common-item and common-object equating now label an unweighted
  descriptive shift when fewer than two common items or objects have
  usable variances. The functions no longer describe that fallback as
  precision-weighted or say that observations used in it were excluded
  from the shift.

- Every location axis is labelled at whole logits when the span allows,
  through one shared tick rule; the previous defaults could leave an
  axis extreme between labels. The person-item map begins and ends its
  proportion axis on labelled ticks, closes the axis corner, and drops
  its unused bottom margin.

- The Wright map and the person-item map mark the person and threshold
  means with dashed lines in their distributions’ colours.

- [`plot_person_fit()`](https://drjoshmcgrane.github.io/rasch/reference/plot_person_fit.md)
  and
  [`plot_item_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_item_map.md)
  can display the standardised infit or outfit in place of the fit
  residual, under the same +/-2.5 band, and annotate the flagged count
  with its percentage. The person table reports `infit_z` alongside
  `outfit_z`.

- Input and selection boundaries are hardened package-wide, closing two
  further review rounds. Every estimator and the DIF family validate
  their arguments through shared checks: iteration caps, tolerances,
  class-interval counts, reference sample sizes, significance levels,
  adjustment methods, practical thresholds, linking minima, component
  and replication counts, and plot bins all reject fractional,
  non-finite, or out-of-range values instead of silently truncating.
  Selection can no longer alter the analysis silently: EFRM refuses
  misspelled or missized id, factor, and item inputs exactly as ordinary
  Rasch does, matrix input honours the items argument in both,
  multiple-choice keys refuse fractional option scores and duplicate
  item entries, replication counts must be whole and finite, explanatory
  threshold labels and predictor rows are validated, duplicate column
  names are refused in direct estimation, and unknown names in anchor,
  distractor, and object-set requests are errors rather than silent
  drops. Pooled MFRM items flow through the DIF follow-ups under the
  stricter item resolver.

- Nine defect families from an adversarial review are corrected.
  Wide-format many-facet data now scores factor columns by their labels,
  where reordered factor levels previously shifted item locations. Item,
  threshold, and anchor indices are validated: a fractional or unknown
  index errors instead of silently truncating to a different item, and
  the numeric fitting controls reject fractional interval counts and
  non-finite references. Duplicate named mappings – comparison anchors,
  item-set maps, judge-panel maps, and named judge factors – are refused
  instead of silently taking one of the conflicting values. A boundary
  comparison object keeps its extrapolated location when a dependence
  effect was dropped after its removal, reports count-weighted
  comparisons, stays inside the plotted range, and is refused by with
  its calibrated companions unaffected. Explanatory departure
  probabilities are withheld for items whose thresholds the calibration
  marks as weak, and the judge diagnostics report count-weighted
  comparison totals.

- A row that is dropped can no longer define the analysis it is dropped
  from. The paired-comparison response scale is derived after zero-count
  and unusable rows are removed, where a single zero-count row could
  take an otherwise identical model from three categories to six; margin
  levels and the presence of ties follow the kept rows for the same
  reason. An ordered factor’s declared levels remain the stated scale
  whatever the weights are, so an empty declared extreme is still
  refused as the identifiability question it is. A tie carries no
  margin, so its margin value no longer opens win and loss categories
  nothing was judged in, which had left both extremes empty and stopped
  the fit; and a response column that is entirely missing reports no
  usable comparisons rather than failing on an impossible vector length.

- A missing person identifier is unknown, not shared. The DIF family no
  longer reads missing identifiers as repeats – which declared a
  repeated-measures design and changed every test, and in
  [`dif_size()`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
  withheld every standard error and Wald test in the table – and the
  tailored bootstrap resamples each unidentified row as its own person
  rather than clustering them into one. A genuinely repeated design
  still withholds its Wald inference.

- Item and object banks are read through their labels. A factor column
  of locations, standard errors or maximum scores was passed to
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html), which returns
  level codes, so a bank of -2.5, 0.25 and 4.0 equated as 1, 2 and 3.
  Numeric text remains accepted; invalid text and non-numeric classes
  are refused. An attached joint covariance now completes individual
  missing standard errors rather than doing so only when the bank
  omitted the whole column. Bank, predictor, key and anchor tables also
  require unique column names.

- A judging sequence must order the comparisons it describes. Values
  repeating within a judge left the order of those comparisons to the
  row order of the data, so the same data read in a different order
  carried different exposure and carry-over covariates; repeated and
  non-finite sequence values are now refused. A retained non-tie margin
  must be an ordered factor or a finite positive magnitude; zero denotes
  a tie. Margins on ties and excluded rows do not define the response
  scale. Logical, complex and Date columns are refused. A replication
  count must be real, since coercing a complex one discards its
  imaginary part in silence.

- A paired-comparison call must state one outcome. Supplying both
  `winner` and `response` fitted the response and ignored the winner, so
  changing every winner left the fit identical; the combination is now
  refused, in
  [`btl_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/btl_explanatory.md)
  too. Every column-role argument names exactly one existing column
  before it is dereferenced, in the comparison, frame-adjusted
  comparison, explanatory comparison and many-facet entry points alike,
  and in the many-facet case through the wide entry as well as the long
  one, which needs at least one item column. Many-facet person, item,
  score, facet and person-factor roles must also be distinct.

- Paired-comparison equating uses the documented precision-weighted
  origin shift when two common objects have usable standard errors.
  Three remain necessary for object-level drift inference, but that
  inferential threshold no longer changes the descriptive link
  estimator.

- Named item-set, judge-panel, judge-factor and Wright-map assignments
  must cover their fitted units exactly. Missing entries no longer
  discard units, and extra entries no longer pass as unnoticed spelling
  errors. Empty or blank panels are refused. Repeated-measure person and
  time roles must be distinct.

- Simulation counts, parameters and seeds are read only as plain numeric
  values. Factors and classed or complex vectors can no longer be
  interpreted through their storage codes, and replicate seeds cannot
  overflow the integer range. A zero-effect BTL dependence specification
  is treated as no planted dependence.

- The low-level
  [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)
  and
  [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md)
  estimators now refuse negative, non-consecutive, constant and
  all-missing item-score columns. Their pair tables require at least two
  observed categories numbered from zero; negative values were
  previously omitted by
  [`tabulate()`](https://rdrr.io/r/base/tabulate.html) rather than
  diagnosed. The main
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  entry point also refuses non-finite scores rather than treating them
  as ordinary non-numeric missing entries.

- MFRM and EFRM exports and app displays draw their observed residual
  scree without requesting the parallel reference that is unavailable
  for their virtual response-cell designs.

- A multiple-choice scoring table refuses a missing or blank item name
  in every key form – the option/score table, the item/key table and the
  named vector – as well as a missing or blank option, and a key refuses
  a blank value; both had scored their item zero throughout and then
  dropped it as constant under a misleading message. A saved analysis
  file’s schema is read as stored rather than coerced, so “1”, TRUE, 1.5
  and a factor’s level code are no longer accepted as schema 1. An
  analysis file that declares a model type that is not one character –
  several, missing, or a factor whose integer code would select a branch
  by level order – is reported as unsupported instead of failing inside
  a length-one condition.

- External person factors supplied as a data frame are checked for
  constancy within person exactly as named columns are, instead of the
  first row’s value being kept. A frame group given by value no longer
  removes a person factor whose name matches one of the group labels,
  and a named panel map must give one stated panel for every judge.

- Set names are trimmed and validated in both frame families: a
  whitespace-only set name, or a blank set in the item-to-set map, is
  refused rather than fitted.

- An empty frame definition is refused rather than dropped: an item set
  or object set naming no identifier would have been fitted away,
  answering a different question from the one asked. A paired-comparison
  DIF call needs at least one judge factor.

- Administration patterns are built with explicit person-by-set
  dimensions. A single respondent – one person in a group, or a
  one-person fit – simplified to a vector that was then read transposed,
  so the design could name sets or facet cells the person never saw.

- A person factor may not take a name the fitted person table generates
  for itself: a factor called `class_interval` silently replaced the
  intervals every fit statistic is computed over, and one called `theta`
  stopped the fit. The generated names are reserved centrally. A frame
  design refuses repeated group columns, and the generated name for
  items no set lists is refused when a nominated set already uses it.

- The standalone HTML report names the estimator the fit used, as the
  fit summary table does, instead of always reporting pairwise
  conditional estimation.

- Names are carried as values rather than parsed out of labels: an
  item-set name containing the label separator keeps its items in the
  score curves; a predictor level no item or object carries is dropped
  before the design is built, where it added an all-zero column and made
  an identified model look rank-deficient; and a resolved comparison
  copy whose generated name already belongs to another object is refused
  with its magnitude withheld, rather than the two silently merging. A
  report table longer than its display limit carries the omission note
  as a caption, so it still renders as a table.

- EFRM score curves are keyed by the administration as well as the
  group: people in one group who sat different item sets have different
  maximum scores and different expected totals, and previously shared
  one curve. The table gains `design` and `n_persons` columns.

- A response style redraws from the probabilities the response was drawn
  under, so planted local dependence survives it, and a style of zero
  strength or zero prevalence is no longer recorded as an active effect.
  In a chain of dependence pairs the middle item’s expectation now
  includes its own carry-over, where the residual it passed on otherwise
  carried the first pair’s shift as a systematic mean into the second.

- An object set aside at a response boundary keeps its predictor row, so
  an explanatory comparison model no longer fails when one object always
  wins or always loses; a predictor row for an object the comparisons
  never mention is still an error.

- [`plot_pcc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pcc.md)
  and
  [`plot_kidmap()`](https://drjoshmcgrane.github.io/rasch/reference/plot_kidmap.md)
  refuse a person identifier that appears in more than one row, naming
  the rows, rather than silently drawing the first. Repeated identifiers
  are the ordinary case in stacked and racked longitudinal data.
  [`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
  validates the repeated-measures identifier where it is described
  rather than failing later on length.

- [`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  computes the covariance of the centred location differences as C{V1 +
  V2}C’, not C{V1 + V2}C. The centring matrix is not symmetric when the
  compared items differ in maximum score, so the standard errors,
  statistics, p values, rmse and ratio were wrong in that case: a
  polytomous item’s standard error was inflated (22% in a five-category
  example) and every dichotomous item’s deflated, flagging short items
  and hiding long ones. Equal maximum scores were unaffected.

- Derived fits keep the controls they were built from.
  [`lr_test()`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)‘s
  rating-scale refit carries the reference sample size and the person
  factors, so both parameterisations’ item-trait statistics are on one
  scale and the refit can be used for follow-up analyses; the fourth
  step of
  [`tailored_analysis()`](https://drjoshmcgrane.github.io/rasch/reference/tailored_analysis.md)
  carries the reference sample size its three siblings use; and a
  subtest total is no longer read as missing when it happens to equal a
  missing-data code, which had deleted complete responses and renumbered
  the categories that remained.

- [`chisq_detail()`](https://drjoshmcgrane.github.io/rasch/reference/chisq_detail.md)
  reports the rescaled per-interval components beside the standardised
  ones, so the detail and the item table reconcile under `adjust_N`, and
  gives the adjustment factor and the unadjusted total.

- Item and object drift are refused for explanatory calibrations, where
  a location is a function of its predictors: a drifted item is smeared
  over every item sharing its design cell and the standard errors belong
  to the design coefficients. The dimensionality magnitude is refused
  for the same reason, since its subtest refit frees every superitem. An
  explanatory comparison design now centres on the objects actually
  calibrated, so its locations sit on the model’s own origin.

- An inestimable DIF analysis is refused rather than returned malformed,
  [`plot_frames()`](https://drjoshmcgrane.github.io/rasch/reference/plot_frames.md)
  draws without intervals when the units carry no standard error instead
  of failing on an empty range, the frame-invariance summary counts item
  comparisons and items separately, and the observed points of the item
  displays follow the fit’s own per-item class intervals, so a graphical
  fit check shows the intervals of the test it illustrates.

- Simulation plants what it records, in three further cases: speededness
  needs a not-reached tail, a dependence source may not be regenerated
  as a later target, and a bias planted on a rater who answers at random
  is refused. Frame group units follow their group labels rather than
  the sorted level order, which attached the wrong unit to each group
  from ten groups upward.
  [`sim_apply()`](https://drjoshmcgrane.github.io/rasch/reference/sim_apply.md)
  requires one atomic scalar per replicate.

- Reports and exports check their arguments before they write: the
  output path, title, plot dimensions and resolution are validated
  first, so a bad size can no longer leave a populated folder that reads
  as a complete export. A plot archive is written fresh rather than
  appended to, and the report closes only the devices it opened. Report
  tables state how many rows were omitted and print small probabilities
  in the package’s own vocabulary rather than as an impossible zero.

- [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  treats presentation as data and the position covariate as a model
  term: paired-comparison fits differing only in which object was
  presented first, or in the judging sequence, are no longer reported as
  the same data, while a plain fit and a position-effect fit of the same
  comparisons still are.

- Four procedures no longer relax an explanatory restriction in silence.
  [`lr_test()`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)
  refuses an explanatory fit, whose rating re-parameterisation would
  drop the design and leave the two models not nested;
  [`tailored_analysis()`](https://drjoshmcgrane.github.io/rasch/reference/tailored_analysis.md)
  refuses one, because the tailored recalibration would differ from the
  original by the design as well as the tailoring; and
  [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  refuses an explanatory comparison fit, whose resolved copies would
  either be forced equal by the design or estimated without it. The
  parallel scree reference now analyses each simulated draw under the
  model that was fitted, where an explanatory calibration was previously
  compared against an unrestricted refit.

- Simulation plants what it records. An item named as the second element
  of several dependence pairs now carries every one of them, where a
  later pair regenerated the item and erased the earlier dependence; a
  dichotomous item cannot be disordered, so the request is refused with
  a warning instead of recorded as planted; and a set or group unit
  ratio must be 1 when there is only one set or group, since a ratio
  between frames cannot be planted in a single frame.

- A failed plot export no longer closes the caller’s graphics device: it
  closes only a device it opened itself. Item names that sanitise to the
  same file stem now keep separate files, where the later plot silently
  overwrote the earlier and the export still looked complete.

- [`plot_pcc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pcc.md)
  draws the person characteristic curve from the fitted model’s own
  expectations. For a polytomous, many-facet, or explanatory fit it
  previously drew a dichotomous logistic curve, which ignores the
  thresholds and the frame’s units; the curve is now the expected
  proportion of the maximum score across the fitted item locations, and
  a dichotomous fit with a common discrimination keeps its exact
  logistic form.

- Selection can no longer alter a multiplicity family in silence.
  Duplicate items, objects, contrast names, and contrast cells are
  refused, because a repeated hypothesis quietly changes the Holm
  adjustment; a grouping given to the DIF family must name fitted
  factors or supply one value per person, rather than being recycled
  into a grouping the fit never contained; dimensionality subsets must
  be free of duplicates and disjoint; and a person-factor frame must
  carry unique, non-empty names whichever input branch assembles it.

- Sustained adversarial review closed the remaining input and display
  boundaries. Role selection is unambiguous in
  [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  and
  [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  alike: a vector as long as the person count is read by value even when
  its labels collide with column names, and a vector whose values match
  a data column exactly is refused, with `items` named as the
  resolution, where an item whose responses happened to agree with a
  role vector was previously dropped from the fit without notice. A
  repeated-measures identifier must carry one entry per fitted row: a
  short vector was recycled, understating the sampling units and
  changing every test in the DIF table. A person with a missing or blank
  frame group is refused rather than expanded into rows missing from
  every set. Reshaping refuses missing or blank person identifiers and
  occasions and requires one existing column name for the person and
  occasion; paired comparisons, frame-adjusted comparisons, and
  many-facet data refuse blank objects, judges, panels, persons, items,
  and facet levels, since a whitespace label is not a level but was
  calibrated as one; and every stated flag must be TRUE or FALSE.
  Display controls are checked before anything is drawn: an evaluation
  grid needs at least two finite locations, class intervals and bins
  must be whole numbers, limits must be two finite ascending values,
  label sizes and colour caps must be positive and finite,
  single-person, single-item, and single-facet displays take exactly one
  name, and limits admitting no thresholds report an empty range instead
  of drawing an empty panel. Batch plot exports and the HTML report name
  the plots they could not draw, where a file was previously written
  with the failures omitted; a batch in which nothing could be drawn is
  now an error rather than a returned path to an archive that was never
  created, and the archive is confirmed on disk before the path is
  returned. Export device dimensions are validated before any device is
  opened.

- The observed points of a paired-comparison ICC no longer abort the
  display when every comparator falls below the informativeness
  threshold; the model curve and the omission note draw on their own.

- A package-wide sweep for reporting-table misalignments found and
  corrected three defects. A numeric anchor index now resolves against
  the data as supplied, where previously a dropped constant item shifted
  the anchor to the wrong item. The class-interval detail refuses an
  item whose only responders are extreme, instead of failing obscurely
  and aborting exports.
  [`dependence_magnitude()`](https://drjoshmcgrane.github.io/rasch/reference/dependence_magnitude.md)
  withholds its standard error and probability when the resolved
  thresholds are weakly identified, reporting the magnitude
  descriptively, as the thresholds themselves already were.

- An undefeated or winless comparison object is now reported in the
  object table at an extrapolated location, its score moved half a point
  inside the boundary against the calibrated scale, with
  `extreme = TRUE` and its standard error withheld – the reporting
  practice already used for extreme person measures. The row takes no
  part in estimation, inference, or equating. Validated against
  [`sirt::btm`](https://rdrr.io/pkg/sirt/man/btm.html): identical
  likelihood to machine precision on clean replicates, with the boundary
  policies agreeing in direction.

- The heaviest examples are smaller, and CRAN runs fewer scenario test
  blocks, keeping the check well inside the incoming pretest budget.

## rasch 1.12.0

CRAN release: 2026-08-24

### Models and inference

- [`wright_map()`](https://drjoshmcgrane.github.io/rasch/reference/wright_map.md)
  sends fitted person and item estimates to `WrightMap`. It supports
  several person distributions and the item-panel layout introduced in
  WrightMap 1.5, including the person-group and item-set structure of
  EFRM fits.

- [`rasch_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_explanatory.md)
  fits the linear logistic test model and linear partial credit model
  from continuous, categorical or ordinal item- or threshold-level
  predictors. Formulae may include selected interactions.
  [`explanatory_test()`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_test.md)
  compares the restrictions with a free calibration using the Kent
  adjustment;
  [`explanatory_diagnostics()`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_diagnostics.md)
  and
  [`relax_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/relax_explanatory.md)
  support fixed item and threshold departures. Refitted departures
  propagate to item and person estimates and are retained through item
  deletion, DIF splitting, superitem construction and
  response-dependence resolution. Keyed option responses remain
  available after item deletion, splitting and fixed-departure refits.

- [`btl_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/btl_explanatory.md)
  applies a fixed explanatory design to object locations in dichotomous
  or ordered comparative judgements. It supports the same model
  comparison, Holm-adjusted diagnostics and fixed departures while
  retaining the nominated ordered-response threshold structure.

- A worked case study on the documentation site uses the verbal
  aggression data to develop and check an explanatory partial credit
  model.

- [`explanatory_test()`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_test.md)
  now places the Kent-calibrated probability in both `p` and `p_kent`.
  The unscaled composite-likelihood probability is named `p_naive` so it
  cannot be mistaken for the inferential result. The table also reports
  calibration R-squared, with an adjusted counterpart whose null
  expectation is near zero, against the free threshold or object
  calibration.

- Pairwise conditional calibrations now use the remaining Newton move as
  a second convergence check. This prevents numerical false refusals at
  large sample sizes without changing the estimates.

- Holm adjustment for item-fit statistics now excludes items whose tests
  are unavailable. Their probabilities remain `NA`.

- Sparse-unit safeguards now use the sampling units that inform each
  test. MFRM interaction tests use the least-supported item-by-level
  cell; EFRM unit tests require adequate persons on every group or set
  link; BTL-EFRM judge bootstraps require adequate effective judges in
  every panel or link; and frame-invariance tests exclude weak frame
  calibrations.

- BTL-EFRM judge bootstraps now distinguish refit errors from
  non-convergence and report the underlying worker error when parallel
  refits fail.

- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
  supports several facets and an optional item-by-facet interaction.
  Omnibus and cell follow-up tests use the fitted joint covariance.

- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  supports crossed person-group factors and reports their GLS factorial
  decomposition. Set-unit linking uses a finite-grid semiparametric
  likelihood, with a separate nuisance distribution for each observed
  person group. Hybrid standard errors retain the joint uncertainty of
  the within-frame calibration and set link; full person-bootstrap
  inference remains available. The convergence flag covers both
  estimation stages, and non-converged links are excluded from bootstrap
  covariance calculations. EFRM data require one response row per
  person.

- The repeated semiparametric linking calculations in the EFRM bootstrap
  now use a compiled numerical kernel. Bootstrap replicates can also be
  distributed over a reproducible, cross-platform worker cluster. The
  Shiny application runs EFRM fits in a background process, defaults to
  four workers where the system permits, records the bootstrap seed,
  reports progress and permits the fit to be cancelled without retaining
  a partial result.

- BTL-EFRM judge bootstraps likewise default to four workers where
  available. A fixed seed gives the same result for any worker count.
  The application runs these fits in the background and supports
  progress reporting and cancellation.

- [`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  compares item locations and discrimination across separately
  calibrated frames. The conditional method tests locations and reports
  discrimination descriptively. The whole-person bootstrap within person
  group preserves item-set response patterns and provides inference for
  both, with one combined Holm family.

- MFRM and EFRM summaries report item estimates separately from the
  item-by-facet or item-by-frame response cells used in estimation.
  Coefficient alpha is not reported for the expanded response-cell
  matrix. EFRM DIF tests pool residual evidence by item and exclude the
  person factors that define the frames.

- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md),
  [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  and
  [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  add ordered paired comparisons, judge-clustered inference,
  judge-factor DIF, linked object sets and judge panels.
  Paired-comparison diagnostics now include equating, transitivity,
  residual dimensions, design information and adaptive pair selection.

- Carry-over probabilities are withheld below 30 judges.
  [`btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)
  uses Welch–Satterthwaite degrees of freedom when fitted calibrations
  have a finite number of judge clusters. Conditional BTL-EFRM unit
  probabilities are withheld; the application defaults to the judge
  bootstrap. BTL-EFRM judge bootstraps use finite-judge references,
  whereas its independent-outcome parametric bootstrap uses normal and
  chi-square references.

### Differential item functioning

- Confirmatory multiplicity defaults are now consistently Holm
  familywise adjustments across item fit, DIF, equating and the
  application. BH remains available where false-discovery-rate screening
  is explicitly requested. BTL DIF uses HC3 covariance for unequal judge
  workloads and withholds omnibus probabilities below eight judges or
  eight effective judges in a factor cell.

- Item-fit documentation now distinguishes the principal item-trait test
  from the supplementary class-interval ANOVA and notes the limits of
  both in short administrations. HC3 was evaluated for item fit and was
  not adopted.

- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  fits several person factors jointly using Type II sums of squares.
  Repeated measurements use the person as the sampling unit and separate
  between- and within-person error strata. Multiplicity adjustment
  covers the complete family of uniform and non-uniform DIF tests rather
  than treating each term as a separate family;
  [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  follows the same rule. Uniform between-person terms now use HC3
  covariance. Class-interval interactions retain the residual-ANOVA
  reference used for non-uniform DIF.

- [`dif_contrasts()`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md)
  and
  [`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
  provide planned and post-hoc logit contrasts, including simple effects
  and difference-in-differences for interactions. MFRM follow-ups pool
  the fitted facet cells of an underlying item; resolved EFRM follow-ups
  are withheld because an ordinary split would discard the frame units.
  The residual-mean Tukey table has been removed from
  [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md);
  [`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
  is the supported follow-up for multilevel terms.

- Repeated-measures DIF follow-ups use the full design-cell weights in
  their person-level tests. Reported resolved estimates and
  probabilities therefore address the same marginal contrast when
  nuisance factors are imbalanced.

- [`dif_size()`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
  reports resolved pairwise logit differences. Dichotomous items receive
  the itemwise ETS A/B/C classification. Polytomous items report the PCM
  signed expected-score area descriptively, without importing an
  incompatible score-metric classification.

- [`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
  splits confirmed DIF items iteratively while retaining a minimum
  anchor set. Automatic splitting is restricted to uniform DIF;
  non-uniform DIF remains visible for item review. MFRM residuals can be
  pooled to their source items, and EFRM factors that do not define
  frames can be tested.

- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  retains anchors and fitted dependence terms in its resolution refit.
  Resolved pairwise inference is withheld unless each factor cell has at
  least eight effective judges; pairwise degrees of freedom use the two
  cells’ effective counts, and the pairwise table reports the raw and
  effective support for both cells. BTL-EFRM fits require a
  frame-specific analysis rather than the equal-unit resolution model.

### Diagnostics and model changes

- Identification checks now cover item, facet, frame and
  paired-comparison graphs, rank, separation and sparse categories.
  Unidentified estimates are refused; identified but weak estimates are
  marked or have inference withheld.
- [`dependence_magnitude()`](https://drjoshmcgrane.github.io/rasch/reference/dependence_magnitude.md)
  uses the joint covariance of resolved thresholds. Equating tests
  require independent calibrations and the covariance of banked
  locations.
- [`spread_test()`](https://drjoshmcgrane.github.io/rasch/reference/spread_test.md)
  applies the binomial least-upper-bound only to superitems formed
  entirely from dichotomous components. It now distinguishes a point
  estimate below the bound from adjusted one-sided evidence of
  dependence. Its significance level and multiplicity adjustment are
  available in the application. The component structure is retained
  through subsequent item splits and removals.
- The tailored-analysis bootstrap resamples complete persons, including
  all rows of a repeated-measures record.
- [`drop_items()`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md),
  [`resolve_frames()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_frames.md),
  DIF splitting and superitem construction refit the active model and
  update downstream item and person estimates. Refit specifications
  retain anchors, keyed scoring, threshold constraints, factors and
  frame-linking controls; a non-converged downstream calibration is not
  returned as a completed analysis.
- Classical whole-test statistics and the Guttman scalogram are withheld
  when an item is represented by several facet or frame response cells.
  They remain available for a one-cell-per-item reduction.
- MFRM characteristic and information curves now combine facet
  conditions administered to the same person. Distinct rating designs
  receive separate curves rather than being added into a test no person
  received.
- Automatic model comparisons are available for the main model families.
  Structural changes are accompanied by before-and-after item and person
  summaries.

### Application and documentation

- The Shiny application uses responsive control and result columns,
  compact explainers for outputs and options, scalable plots and
  downloadable tables. Plot controls sit below the plot, and related
  item curves may be overlaid.
- Analyses can be saved as `.rasch` projects and reopened. Reports can
  be produced as self-contained HTML, Word or PDF documents; the R code
  for each displayed result is available in the application.
- The application covers the extended model suite, including model
  comparison, DIF follow-ups, frame-invariance checks and refitted
  structural changes.
- The manuals and vignettes have been revised to state the fitted
  models, estimands, identification requirements and uncertainty methods
  directly.
- The shipped EFRM and BTL-EFRM case studies now use the current linking
  and uncertainty methods.
- [`plot_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_scree.md)
  and
  [`plot_btl_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_scree.md)
  label their component axes at whole components only, instead of
  overprinting the default axis.

## rasch 1.11.7

CRAN release: 2026-07-30

- [`print()`](https://rdrr.io/r/base/print.html) preserves the reference
  distribution used by saved BTL fits: current and transitional results
  are labelled `t`, while older results without cluster degrees of
  freedom retain their original `z` label.

## rasch 1.11.6

- BTL print methods read both current and earlier dependence-statistic
  columns. Current clustered statistics are labelled `t`, in accordance
  with their t reference distribution.

## rasch 1.11.5

- Ordered paired-comparison margins must be numeric or an ordered
  factor; unordered categorical margins are rejected.
- Invalid graded responses and malformed secondary-dimension
  specifications now produce errors rather than being coerced or
  dropped.
- Clustered dependence and position statistics are labelled `t`.

## rasch 1.11.4

- Score validation reads factor labels and rejects non-integer,
  non-numeric, or non-finite values.
- Ordered BTL responses require an ordered factor or integer scores.
- Clustered dependence and position tests use a t reference with `G - 1`
  degrees of freedom.
- [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md)
  validates secondary-trait correlations and item sets.

## rasch 1.11.3

- All estimators share the same integer-score validation.
- [`item_moments()`](https://drjoshmcgrane.github.io/rasch/reference/item_moments.md)
  uses a log-sum-exp calculation for wide category ranges.
- BTL object-separation reliability is withheld when the clustered
  covariance is rank deficient.
- [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  requires each judge to belong to one panel.
- Undefined reliability and omnibus statistics are reported as `NA`.
- Secondary-trait simulation retains its requested mean and standard
  deviation, and MLE score-table calculations use the correct
  common-unit score equation.

## rasch 1.11.2

- The BTL dimensionality reference simulates count-weighted data at the
  unordered-pair level and returns the leading strength from each
  replicate.
- Judge-bootstrap EFRM inference requires more than one judge per panel
  and notes panels with fewer than five judges.
- BTL fits report when the number of judge clusters cannot support a
  full-rank clustered covariance.

## rasch 1.11.1

- BTL dimensionality reference draws now reproduce count-weighted
  binomial or multinomial sampling.
- Pairwise chi-square degrees of freedom include position and dependence
  parameters; untestable designs return `NA`.
- Judge-clustered covariance uses the CR1 small-sample factor.
- `btl_efrm(se_method = "judge_bootstrap")` resamples judges within
  panels and refits both stages.
- The interpretation of
  [`btl_next_pairs()`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)
  as a ranking rule is stated in its documentation.

## rasch 1.11.0

- Threshold and item-location covariance is transformed consistently
  after recentering items with different maximum scores.
- Warm WLE uses the common-discrimination score equation in which the
  common discrimination cancels.
- Untestable item-trait statistics return `NA`, and the total test
  includes testable items only.
- Equating drift tests account for the estimated origin shift and joint
  covariance of common-item locations.
- Judge-clustered inference requires more than one judge and reports a
  caution for fewer than ten clusters.
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  uses the judge as the sampling unit and count-weighted opponent bands.
- [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  and [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md)
  warn when estimation has not converged.

## rasch 1.10.3

- [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  reports unknown `id`, `factors`, and `items` columns as errors.
- Fractional scores are rejected rather than truncated.
- MFRM rows with missing design identifiers are omitted with a note.
- [`equate_tests()`](https://drjoshmcgrane.github.io/rasch/reference/equate_tests.md)
  excludes common items without usable locations or standard errors from
  weighted linking and drift inference.
- [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md)
  escapes data-derived labels and notes.

## rasch 1.10.2

- The Shiny application adds consistent hover labels to scalograms,
  residual heatmaps, equating plots, tailored analysis, and
  paired-comparison displays.

## rasch 1.10.1

- Person- and item-fit plots in the Shiny application show identifiers,
  locations, and fit residuals on hover.

## rasch 1.10.0

- [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  adds composite-likelihood AIC and BIC based on the Godambe effective
  parameter count.
- Model comparison now covers BTL position, threshold, and dependence
  specifications. The Shiny comparison page supports these fits.
- Optional cross-package tests compare results with `eRm`, `sirt`, and
  `psychotools`.

## rasch 1.9.3

- Pairwise estimation now checks that the observed item graph, together
  with any anchors, identifies a common scale.
- Thresholds adjoining critically sparse categories are marked weak.
  Their threshold and item-location standard errors are withheld.

## rasch 1.9.2

- Case study: `inst/casestudies/party_blocs_crisis.R` applies
  [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  to the Tuebingen 2009 party-preference data.
- Sets without stable panel-ratio information are excluded from the unit
  reconciliation and refitted at the reconciled panel units.
- Boundary-unstable bootstrap parameters receive `NA` standard errors
  with the number of boundary replicates reported.
- BTL-EFRM convergence is assessed by the gradient per comparison.

## rasch 1.9.1

- The Shiny Frames page supports BTL-EFRM panel and object-set
  definitions, unit tables, frame fit, and unit plots.

## rasch 1.9.0

- [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  fits the paired-comparison extension of the extended frame of
  reference model, with panel units, object-set units, and set origins.
- Bootstrap standard errors refit both stages. Conditional standard
  errors remain available for descriptive work.
- [`plot_btl_units()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_units.md)
  and
  [`simulate_btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl_efrm.md)
  support display and simulation.

## rasch 1.8.0

- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) adds
  anchored estimation and a first-position effect.
- [`btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)
  and
  [`plot_btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_equate.md)
  provide common-object linking and drift tests for paired-comparison
  calibrations.
- [`btl_information()`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md),
  [`plot_btl_targeting()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_targeting.md),
  and
  [`btl_next_pairs()`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)
  provide design information and greedy next-pair selection.
- Count-weighted BTL sandwich covariance now reproduces expanded-data
  covariance. Equating uses each calibration’s stored covariance.
- The Shiny Targeting and Equating pages support paired comparisons.

## rasch 1.7.1

- Rasch simulation layers now retain all previously specified DIF,
  dimensionality, response-style, dependence, and guessing terms.
- PCM simulation uses item-specific threshold structures and applies
  stricter validation to DIF, guessing, and disordered-threshold
  specifications.
- BTL dimensionality references include fitted within-judge dependence.
- MFRM simulation and recovery use the recorded item and rater
  parameters consistently.
- The Simulate page links recovery output to the current generated
  dataset.
- A simulation vignette and package logo were added.

## rasch 1.7.0

- Simulation functions add population-distribution controls, response
  styles, speededness, and MFRM halo effects.
- [`sim_replicate()`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md),
  [`sim_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/sim_recovery.md),
  and
  [`plot_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/plot_recovery.md)
  support repeated simulation and parameter-recovery summaries.
- The Shiny Simulate page exposes the population controls and recovery
  output.

## rasch 1.6.1

- The Shiny application adds a Simulate page for the four data layouts
  and their model-departure controls.

## rasch 1.6.0

- [`simulate_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_efrm.md)
  gains `n_categories` for partial credit items within frames, with
  planted thresholds recorded in the truth attribute.
- [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md),
  [`simulate_btl()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl.md),
  [`simulate_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_mfrm.md),
  and
  [`simulate_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_efrm.md)
  generate data from the package’s model families and can introduce
  nominated departures. Generating parameters are stored in the returned
  data.

## rasch 1.5.1

- [`plot_btl_judge_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_judge_map.md)
  now displays individual matchups.
  [`judge_pair_surprise()`](https://drjoshmcgrane.github.io/rasch/reference/judge_pair_surprise.md)
  returns the corresponding residuals.

## rasch 1.5.0

- [`judge_surprise()`](https://drjoshmcgrane.github.io/rasch/reference/judge_surprise.md)
  and
  [`plot_btl_judge_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_judge_map.md)
  compare a judge’s object-level preferences with the consensus object
  scale. The display is available from the Shiny Judge fit page.

## rasch 1.4.0

- [`btl_transitivity()`](https://drjoshmcgrane.github.io/rasch/reference/btl_transitivity.md)
  reports circular triads and Kendall’s consistency coefficient for
  suitable paired-comparison designs.
- [`btl_dimensionality()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md)
  decomposes the skew-symmetric residual preference matrix and compares
  its leading component with a model-based reference.
- New plots display BTL transitivity, scree, and residual maps.

## rasch 1.3.1

- The package was renamed from its development name, `rmt`, to `rasch`.
  Result classes use the `rasch_` prefix.

## rasch 1.3.0

- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) adds
  count-weighted exposure and carry-over effects, separation handling,
  and
  [`plot_btl_dependence()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dependence.md).
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  carries fitted dependence effects into its residual analysis and
  handles aggregated comparison counts.
- Mixed-design DIF terms are evaluated in their corresponding error
  strata. Repeated-person logit contrasts are provided by
  [`dif_contrasts()`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md).
- Residual parallel analysis uses data simulated from the fitted model.
- The Shiny application retains the settings used by each DIF analysis
  and applies consistent fit flags.

## rasch 1.2.0

- [`plot_pca_biplot()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pca_biplot.md)
  draws the item loadings on the first two residual principal components
  on equal axes.
- [`residual_correlations()`](https://drjoshmcgrane.github.io/rasch/reference/residual_correlations.md)
  now also returns the adjusted-Q3 `star_matrix` and
  [`plot_resid_cor()`](https://drjoshmcgrane.github.io/rasch/reference/plot_resid_cor.md)
  can draw raw Q3 or adjusted Q3\*.
- The Shiny trait and local-dependence pages pair tables with their
  plots and allow the original data to be restored after restructuring.
- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  is now the single DIF analysis-of-variance function. One factor is
  analysed one-way; several factors are fitted jointly. It supports
  repeated-measures and mixed designs.
- [`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
  resolves DIF iteratively by item splitting.

## rasch 1.0.0

First stable release.

### Models

- [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  fits dichotomous, partial credit, and rating scale models by pairwise
  conditional maximum likelihood, with Warm WLE person estimates.
- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
  fits additive and item-by-facet many-facet models.
- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  fits the extended frame of reference model.
- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) fits
  dichotomous and ordered paired-comparison models.

### Diagnostics

- Item and person fit, category functioning, targeting, reliability,
  test information, residual dimensionality, and local dependence.
- DIF analysis, item splitting, tailored analysis, common-item equating,
  and anchored calibration.
- BTL judge fit, invariance, and paired-comparison diagnostics.

### Display and reporting

- Base-graphics functions cover item, person, threshold, targeting,
  information, residual, and paired-comparison displays.
- [`fit_summary_table()`](https://drjoshmcgrane.github.io/rasch/reference/fit_summary_table.md)
  and
  [`targeting_table()`](https://drjoshmcgrane.github.io/rasch/reference/targeting_table.md)
  return the headline statistics;
  [`save_outputs()`](https://drjoshmcgrane.github.io/rasch/reference/save_outputs.md)
  and
  [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md)
  export results.
- [`run_app()`](https://drjoshmcgrane.github.io/rasch/reference/run_app.md)
  launches the Shiny interface and shows the R call corresponding to
  each analysis.
