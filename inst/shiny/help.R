# Canonical Shiny explainers
# ---------------------------------------------------------------------------
# Each substantive result card is keyed by its Shiny output id.  The app uses
# these short definitions in click/focus popovers; longer derivations and
# qualifications belong in the manual and vignettes.  Keeping the copy here
# prevents the same statistic acquiring different definitions on different
# pages.

APP_HELP <- c(
  # Headline metrics -------------------------------------------------------
  metric_persons = paste(
    "Number of persons represented in the active analysis. After a structural",
    "change, this refers to the current fitted response matrix."
  ),
  metric_items = paste(
    "Number of items represented in the active analysis. A superitem reduces",
    "this count; a DIF split can increase it."
  ),
  metric_cells = paste(
    "Number of calibrated response columns. These are items in an ordinary",
    "fit and item-by-frame or item-by-facet cells in the corresponding",
    "structural model."
  ),
  metric_psi = paste(
    "Person Separation Index: the proportion of observed person-location",
    "variance not attributed to measurement error. Larger values indicate",
    "more reproducible separation of persons."
  ),
  metric_alpha = paste(
    "Coefficient alpha describes raw-score internal consistency and is",
    "reported beside the model-based PSI for an ordinary administered item",
    "matrix. It is not applicable when an item is represented by several",
    "frame or facet response cells."
  ),
  metric_item_trait = paste(
    "Approximate asymptotic probability for the overall item- or",
    "response-cell-trait chi-square. It treats estimated person locations as",
    "known and can be miscalibrated; use the item-fit bootstrap where available.",
    "It is withheld when person IDs repeat because the reference counts rows,",
    "not independent persons."
  ),
  metric_power = paste(
    "Describes person separation from the PSI in broad bands.",
    "It is not the statistical power of a fit test, which also depends on the",
    "sample, targeting, test statistic and departure being tested."
  ),
  metric_objects = "Number of objects represented in the active comparative judgement analysis.",
  metric_comparisons = paste(
    "Number of usable comparative judgements contributing to the active fit,",
    "including count weights where supplied."
  ),
  metric_judges = "Number of nominated judges contributing usable comparisons to the active fit.",
  metric_osi = paste(
    "Object Separation Index: the proportion of observed object-location",
    "variance not attributed to measurement error. Larger values indicate a",
    "more reproducible object ordering."
  ),
  metric_pair_fit = paste(
    "The overall pairwise fit probability. After a fit bootstrap, it uses the",
    "fitted-design null. Small values indicate departure from the fitted",
    "object-pair expectations. Without a bootstrap, it is unavailable when",
    "judges contribute repeated comparisons. BTL–EFRM retains the statistic",
    "descriptively but does not report a whole-fit probability."
  ),
  metric_item_misfit = paste(
    "Before bootstrapping, this is the number of items with an approximate",
    "asymptotic Holm probability below .05. After bootstrapping, it uses the",
    "calibrated bootstrap probability. Adjustment is within that statistic",
    "under the fitted global null. Item-fit inference is unavailable when",
    "person IDs repeat."
  ),
  metric_disordered = paste(
    "Number of polytomous items whose estimated thresholds are not ordered.",
    "Inspect category use and threshold uncertainty before considering rescoring."
  ),
  metric_extreme = paste(
    "Number of persons with a zero or maximum score on all observed items.",
    "Their finite locations require an extrapolated scoring method."
  ),
  metric_person_misfit = paste(
    "Before bootstrapping, the number of persons with an absolute fit residual",
    "above the working value of 2.5. After bootstrapping, it counts adjusted",
    "fit-residual probabilities below .05 under the fitted global null."
  ),
  metric_expl_model = paste(
    "LLTM denotes a dichotomous explanatory model; LPCM denotes its",
    "polytomous threshold extension. Explanatory CJ models object locations."
  ),
  metric_expl_terms = paste(
    "Number of estimable predictor coefficients after fixing the scale origin."
  ),
  metric_expl_relax = paste(
    "Number of fixed item, threshold or object departures in the active model."
  ),
  metric_expl_test = paste(
    "Kent-adjusted comparison of the active explanatory restrictions with a",
    "free calibration of the same responses."
  ),

  # Summary ---------------------------------------------------------------
  fitsum_tbl = paste(
    "Summarises overall Rasch model fit, including the item-trait",
    "chi-square and the distribution of calibration-cell and person fit",
    "residuals. The test examines whether fitted orderings remain invariant",
    "over class intervals."
  ),
  targeting_tbl = paste(
    "Summarises the alignment of person locations and calibration thresholds,",
    "together with separation and reliability. Targeting is strongest when",
    "the thresholds cover the part of the scale occupied by the persons."
  ),
  btl_fitsum_tbl = paste(
    "Summarises fit of the comparative judgement model. Pairwise chi-square",
    "compares observed and expected responses for each object pair; the",
    "object separation index describes reproducibility of the object ordering.",
    "A row-based probability is unavailable when judges contribute repeated",
    "comparisons; BTL–EFRM retains the chi-square only as a descriptive summary."
  ),
  change_est_tbl = paste(
    "Compares original and active item or object locations after a structural",
    "change. The active estimates are used by every subsequent result."
  ),
  change_person_tbl = paste(
    "Compares original and active person locations after a structural change,",
    "showing how the revised calibration propagates into person measurement."
  ),
  expl_test_tbl = paste(
    "Compares the active explanatory model with a free calibration of the same",
    "responses. Use the first-order, asymptotic Kent-adjusted probability for",
    "this multivariate comparison."
  ),
  expl_coef_tbl = paste(
    "Reports predictor effects on item, threshold or object location in logits, with",
    "standard errors and Holm-adjusted probabilities. Repeated-person or",
    "judge-clustered fits use their supported cluster degrees of freedom."
  ),
  expl_diag_tbl = paste(
    "Tests each available fixed item, threshold or object departure separately from",
    "the active model. Holm adjustment covers the complete displayed family."
  ),
  expl_calibration = paste(
    "Compares explanatory and freely estimated locations on the same logit",
    "scale. Points far from the diagonal identify the largest restrictions."
  ),
  expl_relax_tbl = paste(
    "Lists the fixed departures in the active explanatory calibration."
  ),

  # Item, object and person tables ----------------------------------------
  items_tbl = paste(
    "Reports item fit for ordinary calibrations and response-cell fit for",
    "Extended Frames or Multiple Ratings. Select a row to inspect its curves;",
    "marked cells identify working criteria rather than a single decision rule.",
    "Class-interval probabilities are approximate, particularly when fewer",
    "than ten responses locate each person."
  ),
  structural_items_tbl = paste(
    "Reports item locations on the common scale for an Extended Frames or",
    "Multiple Ratings fit. Item-by-frame or item-by-facet response-cell fit is",
    "reported separately because those cells are not additional items."
  ),
  distractor_tbl = paste(
    "Reports how each response option relates to person location and item score.",
    "Options intended to represent increasing performance should generally be",
    "chosen by progressively more able respondents."
  ),
  person_tbl = paste(
    "Reports person locations, standard errors, fit and score information.",
    "Extreme scores and sparse response patterns can carry limited information",
    "and are identified in the table."
  ),
  person_weight_tbl = paste(
    "Reports supplementary person locations from externally imposed relative",
    "item or item-set weights. Sandwich standard errors allow for the weighted",
    "score; the fitted calibration and ordinary Rasch estimates do not change."
  ),
  btl_obj_tbl = paste(
    "Reports each object's estimated location, standard error and fit over its",
    "comparisons. Select a row to inspect the object characteristic curve."
  ),
  btl_thr_tbl = paste(
    "Reports the ordered response thresholds for a polytomous comparative",
    "judgement model. Thresholds locate transitions between adjacent response",
    "categories on the comparison scale."
  ),
  btl_comp_tbl = paste(
    "Reports the principal threshold components used by a constrained",
    "polytomous comparative judgement model. These components provide a more",
    "parsimonious alternative to freely estimated thresholds."
  ),
  btl_pairs_tbl = paste(
    "Compares observed and expected responses for each object pair. Large",
    "departures identify pairs that are not well represented by the fitted",
    "one-dimensional ordering."
  ),
  btl_judges_tbl = paste(
    "Summarises each judge's fit over their comparisons. Large positive fit",
    "residuals indicate response patterns that are less consistent with the",
    "common object ordering."
  ),
  btl_trans_judges_tbl = paste(
    "Summarises circular preferences within each judge's comparisons.",
    "Consistency equals one for a transitive ordering and zero at the",
    "random-tournament benchmark."
  ),
  chisq_int_tbl = paste(
    "Shows the selected item's observed and expected mean score in each",
    "person class interval, with its contribution to the item chi-square."
  ),
  chisq_cat_tbl = paste(
    "Shows observed and expected response-category counts by class interval.",
    "These cells locate departures contributing to the item's overall test."
  ),
  ctt_tbl = paste(
    "Reports conventional item statistics alongside the Rasch estimates,",
    "including endorsement, item-total association and alpha if removed.",
    "It is withheld when an item is represented by several response cells."
  ),
  rescore_tbl = paste(
    "Summarises the evidence used to propose ordered partial-credit scores for",
    "multiple-choice options. Review proposed scores before reanalysis."
  ),

  # Targeting and equating -------------------------------------------------
  btl_info_tbl = paste(
    "Reports the Fisher information supplied by each object's observed",
    "comparisons, retaining fitted position, history and frame effects.",
    "This describes the design; it is not the fitted standard error."
  ),
  btl_next_tbl = paste(
    "Ranks candidate comparisons by the information expected from one",
    "additional judgement. Priority weighting favours the pairs expected to",
    "reduce total location uncertainty most. The stronger object is presented",
    "first. History-dependent fits need a specified judge and history."
  ),
  btl_eq_tbl = paste(
    "Places two comparative judgement calibrations on a common origin using",
    "shared objects. The shift and object differences show agreement between",
    "the calibrations after linking. The shift is precision-weighted when",
    "usable standard errors are available; otherwise it is an unweighted",
    "descriptive mean. An exact common anchor fixes the shift."
  ),
  eq_tbl = paste(
    "Compares common-item locations after the selected origin alignment.",
    "The mean shift is precision-weighted when usable standard errors are",
    "available and otherwise unweighted and descriptive. An exact common",
    "anchor fixes the shift. Where joint",
    "uncertainty is available, adjusted t tests identify item drift and report",
    "the applicable finite-cluster or limiting reference degrees of freedom."
  ),
  eq_plot = paste(
    "Plots the common-item calibrations against the aligned identity line.",
    "Items departing from the line have shifted after origin alignment."
  ),

  # DIF -------------------------------------------------------------------
  dif_tbl = paste(
    "Summarises uniform and non-uniform DIF for the selected item and factor",
    "terms. Holm-adjusted probabilities control the complete item-by-term",
    "family by default. Uniform between-person terms use HC3 covariance;",
    "class-interval interactions retain the residual-ANOVA reference for",
    "complete panels. Incomplete panels use joint adjustment and",
    "person-cluster CR3 covariance for all between-person terms."
  ),
  dif_full_tbl = paste(
    "Contains the complete item-by-term DIF results. Uniform DIF is associated",
    "with a factor effect; non-uniform DIF with a factor-by-class-interval",
    "interaction."
  ),
  dif_boot_tbl = paste(
    "Repeats the complete DIF analysis under the fitted invariant model while",
    "holding each response row's score and observed-item pattern fixed;",
    "Extended Frames instead holds its item-set subtotals fixed.",
    "Bootstrap probabilities use single-step minimum-p familywise adjustment",
    "under the fitted global invariant null. They are a sensitivity analysis",
    "beside the primary residual ANOVA, not strong control after another item",
    "has DIF."
  ),
  bdif_anova_tbl = paste(
    "Summarises object DIF across judge factors. A factor effect indicates a",
    "location difference between judge groups; an opponent-band interaction",
    "indicates non-uniform DIF. HC3 inference allows judge workloads to differ;",
    "a factor cell needs at least eight judges and eight effective judges."
  ),
  bdif_boot_tbl = paste(
    "Draws outcomes from the fitted comparison model, retaining judges,",
    "comparison design, response categories and fitted history effects.",
    "Single-step minimum-p probabilities cover the complete object-by-term",
    "family under the fitted global invariant null. They supplement the",
    "primary residual analysis and do not provide strong control after another",
    "object has DIF."
  ),
  bdif_sizes_tbl = paste(
    "Reports resolved object-location differences between judge-factor levels",
    "in logits, with the number and effective number of judges supporting each",
    "level. Pairwise inference requires at least eight effective judges in each",
    "level; otherwise the differences remain descriptive. Adjusted probabilities",
    "refer to the full displayed contrast family."
  ),
  dif_posthoc_tbl = paste(
    "Resolves the selected DIF term into adjusted comparisons on the logit",
    "scale. Interaction rows are contrasts of contrasts; main-effect rows are",
    "equal-cell marginal comparisons over the other fitted factors. Repeated",
    "designs use the same cell weights in their person-level tests."
  ),
  dif_size_tbl = paste(
    "Compares resolved interaction cells pairwise. Use these rows to locate the",
    "pattern after reading the interaction contrasts above. ETS letters apply",
    "to dichotomous items. Polytomous items report the PCM signed",
    "expected-score area descriptively. Wald tests report their applicable",
    "person-cluster degrees of freedom."
  ),
  resolve_tbl = paste(
    "Records each automatic item split, its triggering term and effect size.",
    "The resulting fit is used by subsequent results."
  ),
  contr_tbl = paste(
    "Reports planned one-degree-of-freedom contrasts derived from the factor",
    "structure: two-level differences, ordered trends, nominal comparisons and",
    "factor-pair interactions. Probabilities are adjusted over this planned family."
  ),

  # Facets and frames ------------------------------------------------------
  facet_tbl = paste(
    "Reports the severity, uncertainty and fit of each facet level. Positive",
    "values denote greater severity on the common logit scale."
  ),
  facet_int_omnibus = paste(
    "Tests the complete item-by-facet interaction family. This omnibus result",
    "should be read before examining individual interaction cells. Inference",
    "uses the least-supported item-by-level cell."
  ),
  facet_int_tbl = paste(
    "Reports item-by-facet departures from the additive many-facet model.",
    "Cell probabilities are adjusted over the exploratory follow-up family."
  ),
  phi_tbl = paste(
    "Reports the relative measurement unit for each person group in an EFRM.",
    "The group units have geometric mean one, so no observed group is the",
    "reference. Intervals and tests are calculated on the log-unit scale."
  ),
  alpha_tbl = paste(
    "Reports each item set's unit ratio: the reference unit over the set's",
    "own, so a value above one means the finer unit and steeper curves.",
    "The link uses persons observed in more than one set, estimating their",
    "distribution on a finite grid within each person group.",
    "The set units have geometric mean one; set origins are separate."
  ),
  frame_tbl = paste(
    "Reports the unit and fit for each observed item-set by person-group frame.",
    "The frame unit combines the relevant set and group units."
  ),
  btlef_phi_tbl = paste(
    "Reports panel units in the comparative judgement frame model. Their",
    "geometric mean is one; adjusted tests compare each panel with equal unit."
  ),
  btlef_units_tbl = paste(
    "Reports set units and origins for linked comparative judgement sets.",
    "Units describe scale changes; origins describe translations between sets.",
    "The first set fixes unit one and origin zero.",
    "The adjusted probabilities form one Holm family across panel units, set units and origins."
  ),
  btlef_frames_tbl = paste(
    "Reports the unit, comparison count and fit of each panel-by-set frame.",
    "Sparse frames should be interpreted alongside their uncertainty and the",
    "linking design."
  ),
  btlef_cmp_tbl = paste(
    "Compares the comparative judgement frame model with its equal-unit fit on",
    "the same comparisons. The likelihood difference is descriptive."
  ),
  btlef_omnibus_tbl = paste(
    "Jointly tests whether each family of panel units, set units or set origins",
    "can be replaced by its equal-unit restriction. Decisions use Holm-adjusted",
    "probabilities across these omnibus tests. Judge-bootstrap tests need",
    "six judges and 5.5 effective judges per panel, and eight per set link."
  ),
  efrm_cmp_tbl = paste(
    "Compares group-dependent and equal group units using the same within-frame",
    "conditional information. It does not test item-set units; the likelihood",
    "difference is descriptive."
  ),
  efrm_omnibus_tbl = paste(
    "Jointly tests the equal-unit restriction for the group and item-set unit",
    "families. Decisions use Holm-adjusted probabilities across the omnibus",
    "tests. They provide the inferential model comparison when",
    "every group and set link has at least 50 contributing persons."
  ),

  # Dimensionality and local dependence ----------------------------------
  pc_tbl = paste(
    "Reports residual principal-component coordinates and their explained",
    "variance. Components describe structure remaining after the Rasch trait",
    "has been fitted."
  ),
  loadings_tbl = paste(
    "Reports item loadings on the selected residual component. Items with",
    "opposing large loadings define the strongest residual contrast."
  ),
  eigen_tbl = paste(
    "Reports the observed residual eigenvalues and the variance represented by",
    "each component. The scree plot adds the fitted-model simulation reference."
  ),
  dm_tbl = paste(
    "Compares reliability when items are treated separately with a refit in",
    "which each item subset is a super-item. The resulting coefficients estimate",
    "unique loading, correlation and common variance across the subsets.",
    "Both use complete-response rows; PSI also requires usable estimates in",
    "both fits. n and n_excluded count response rows used and excluded."
  ),
  dep_tbl = paste(
    "Reports threshold differences after the dependent item is resolved by the",
    "independent item's categories and the model is refitted. Their half-range",
    "estimates the overall dependence magnitude in logits. The Wald test uses",
    "person-cluster degrees of freedom when response rows repeat."
  ),
  spread_tbl = paste(
    "Compares each recorded superitem's threshold spread with the binomial",
    "bound. The one-sided adjusted probability tests whether the spread is",
    "below that bound and reports the applicable person-cluster degrees of",
    "freedom. The bound is available only for dichotomous-item superitems."
  ),
  cormat_q3_tbl = paste(
    "Shows Yen's Q3 residual correlations between items. Large positive values",
    "indicate response dependence remaining after conditioning on the Rasch",
    "trait."
  ),
  cormat_q3s_tbl = paste(
    "Shows adjusted Q3 residual correlations after removing their overall",
    "mean. The adjustment makes unusually dependent item pairs easier to",
    "identify."
  ),
  btl_bimensions_tbl = paste(
    "Reports rotational dimensions of the observed-minus-expected pair logits.",
    "Judge-clustered results are descriptive by default;",
    "the simulated reference requires conditionally independent comparisons."
  ),
  btl_trans_tbl = paste(
    "Summarises circular triads in the observed comparisons. A high loop rate",
    "indicates departures from a single transitive object ordering."
  ),
  btl_dep_tbl = paste(
    "Reports fitted within-judge history effects. Exposure represents having",
    "seen an object before; carry-over represents the influence of earlier",
    "verdicts involving that object."
  ),
  btl_dep_comps = paste(
    "Lists the comparisons that inform the selected dependence effect, with",
    "their history covariates and fitted contribution."
  ),

  # Other analyses ---------------------------------------------------------
  guess_tbl = paste(
    "Compares the original and tailored item calibrations after responses below",
    "the nominated chance probability are removed. Shifts are descriptive unless",
    "bootstrap inference was requested."
  ),
  frame_inv_loc_tbl = paste(
    "Compares item locations from separate frame calibrations after conversion",
    "to the common unit and centring over common items. The selected uncertainty",
    "method is shown with the result."
  ),
  frame_inv_summary_tbl = paste(
    "Summarises each set and frame comparison. RMSD is the observed spread of",
    "the location differences; RMSE is their expected sampling spread. The",
    "counts give the items tested, excluded and flagged. Discrimination counts",
    "are available with bootstrap uncertainty. Each compared frame requires",
    "at least 50 informative persons; weak item estimates are excluded."
  ),
  frame_inv_disc_tbl = paste(
    "Compares item discrimination across separate frame calibrations. Slopes",
    "show direction and relative size; a boundary flag marks estimates at the",
    "fitted limit. Conditional results are descriptive; bootstrap results include",
    "probabilities and Holm adjustment."
  ),
  cmp_tbl = paste(
    "Compares fits retained during this session. Information criteria are",
    "comparable only for models using the same observations and response",
    "scale; other rows support descriptive comparison. Generic EFRM",
    "information criteria and likelihood differences are withheld; use its",
    "dedicated model comparison. Smaller CL-AIC or CL-BIC is preferred."
  ),
  sim_recovery_tbl = paste(
    "Compares generating and fitted parameters for the current simulated data.",
    "Locations are aligned to the model's identifying origin. If the fitted",
    "model does not represent a planted departure, the comparison is marked",
    "as descriptive rather than parameter recovery."
  ),
  sim_preview = paste(
    "Shows the first rows of the simulated dataset currently loaded for analysis.",
    "Download the data alone as CSV or with its generating call, truth and",
    "explanatory metadata in the reproducibility bundle."
  ),
  preview = "Shows the first rows and current column roles of the dataset to be analysed.",

  # Item explorer ---------------------------------------------------------
  icc = paste(
    "Shows expected item score over person location, with observed class-interval",
    "means when selected. Overlays use a common proportional score scale."
  ),
  ccc = paste(
    "Shows the probability of each response category over person location.",
    "Threshold markers locate transitions between adjacent categories."
  ),
  tpc = paste(
    "Shows the probability of passing each item threshold over person location,",
    "with observed class-interval proportions when selected."
  ),
  cfreq = "Shows the observed number of responses in each category for the selected item.",

  # Plots: core Rasch ------------------------------------------------------
  wright = paste(
    "Places person locations and calibration thresholds on the same logit scale.",
    "Extended Frames and Multiple Ratings use response-cell thresholds.",
    "Good targeting places thresholds across the range occupied by the persons.",
    "The optional WrightMap renderer can separate person distributions and",
    "item sets or facets into panels; both default to a single panel."
  ),
  pim_p = paste(
    "Compares the distributions of person locations and calibration thresholds on",
    "their common logit scale. Dashed lines mark the means. Test information",
    "can be added on the right-hand scale."
  ),
  thrmap = paste(
    "Displays every fitted threshold on the common logit scale. Extended Frames",
    "and Multiple Ratings show the response cells used in estimation."
  ),
  tcc = paste(
    "Shows the modelled expected total score over the latent trait. Extended",
    "Frames and Multiple Ratings draw a curve for each administrable design."
  ),
  tif = paste(
    "Shows information over the latent trait for each administrable design.",
    "Higher information indicates greater precision; a single design also",
    "shows its conditional standard error on the right-hand scale."
  ),
  imap = paste(
    "Plots fitted locations against fit residuals. Extended Frames and Multiple",
    "Ratings show response cells; points outside the band warrant inspection."
  ),
  pfit = paste(
    "Plots person locations against person fit residuals. Select or hover over",
    "unusual points to identify response patterns for further review."
  ),
  rdist_i = paste(
    "Shows the distribution of item or response-cell fit residuals against the",
    "standard normal reference expected under adequate fit."
  ),
  rdist_p = paste(
    "Shows the distribution of person fit residuals against the standard normal",
    "reference expected under adequate fit."
  ),
  kidmap = paste(
    "Places the thresholds the person achieved on one side of their",
    "estimated location and the thresholds not achieved on the other,",
    "inside the person's confidence band. Achieved thresholds above the",
    "band are unexpected successes; unachieved thresholds below it are",
    "unexpected failures."
  ),
  guttman = paste(
    "Orders items by difficulty and compares the person's responses with a",
    "deterministic Guttman pattern. The Rasch model remains probabilistic; this",
    "display is a descriptive response-pattern check for an ordinary",
    "dichotomous item matrix."
  ),
  scree = paste(
    "Shows residual eigenvalues against the score-conditional null reference",
    "band. Red points exceed the familywise 5% limit after adjustment across",
    "the displayed components."
  ),
  pca_biplot = paste(
    "Displays items and persons in the selected residual-component space.",
    "Separation of item clusters can indicate a secondary response structure."
  ),
  rcor_q3 = paste(
    "Visualises the Q3 residual-correlation matrix. Strong positive cells identify",
    "item pairs with more shared response variation than the Rasch trait explains."
  ),
  rcor_q3s = paste(
    "Visualises the mean-adjusted Q3 residual-correlation matrix. The adjustment",
    "centres the matrix so unusually dependent pairs stand out."
  ),
  distractor_plot = paste(
    "Shows selection of each response option over person location. Ordered scoring",
    "is supported when higher-scoring options become more likely at higher locations."
  ),
  guess_plot = paste(
    "Compares item locations before and after the tailored analysis on a common",
    "origin. Departures from the identity line show calibration shifts."
  ),

  # Plots: facets, frames and comparative judgement -----------------------
  facet_plot = paste(
    "Displays facet-level severities with confidence intervals. Positive values",
    "denote greater severity and zero is the average facet level."
  ),
  frame_plot = paste(
    "Displays the estimated units of the observed EFRM frames with confidence",
    "intervals. Unit one represents the equal-frame Rasch model."
  ),
  frame_icc = paste(
    "Shows the selected item's model curve and observed class-interval means",
    "within each frame. Differences between frames reflect their estimated units."
  ),
  btl_plot = paste(
    "Displays comparative judgement object locations with confidence intervals.",
    "Objects farther to the right have a greater modelled probability of being",
    "preferred."
  ),
  btl_occ = paste(
    "Shows the selected object's expected response over opponent location, with",
    "observed means for opponents supported by enough comparisons. Position-",
    "and history-adjusted fits use model expectations for those comparisons."
  ),
  btl_cats = paste(
    "Shows the probabilities of the ordered comparative judgement response",
    "categories. Position effects are included. For history-dependent fits,",
    "the axis is the full linear predictor, including position and history terms."
  ),
  btl_judge_map = paste(
    "Shows one judge's observed comparison responses against their modelled",
    "expectations. Red matchups favour the weaker object and pass the Holm",
    "familywise rule. Tied locations have no directional flag."
  ),
  btl_judge_consist = paste(
    "Displays each judge's transitivity consistency against their number of",
    "informative triples. Sparse judges carry less stable consistency estimates."
  ),
  btl_involve_plot = paste(
    "Shows the selected judge's largest object-pair departures. It identifies",
    "which comparisons contribute most to that judge's fit result."
  ),
  btl_targeting_plot = paste(
    "Plots object location against design information, with point size showing",
    "comparison count. Equal-unit fits also show the information expected from",
    "one new comparison; frame and history-dependent fits have no single",
    "reference curve. A position-adjusted curve presents the reference object first."
  ),
  btl_eq_plot = paste(
    "Compares linked object locations from two comparative judgement calibrations.",
    "Agreement with the shifted identity line indicates stable common-object linking."
  ),
  btl_scree = paste(
    "Shows residual bimension strengths. A 5% reference band is shown only",
    "when conditional independence is assumed and the design supports it;",
    "a leading value above that band suggests structured preference cycles."
  ),
  btl_dim_map = paste(
    "Maps the pattern of residual comparisons, not their magnitude or significance.",
    "Use the scree plot to compare strength with the conditional reference."
  ),
  btl_dep_plot = paste(
    "Shows observed residual departure over the selected history covariate, with",
    "the fitted dependence effect and the number of comparisons in each bin.",
    "Probabilities are Holm-adjusted across the fitted dependence effects."
  ),
  bdif_occ = paste(
    "Shows the selected object's model curve and observed responses by judge-factor",
    "level. Separation between groups is the graphical display of object DIF."
  ),
  dif_icc = paste(
    "Shows the selected item's characteristic curve and observed class-interval",
    "means by factor level. Separation between groups is the graphical display of DIF."
  ),
  btlef_units_plot = paste(
    "Displays panel and set units with confidence intervals on the log scale.",
    "Intervals use the reported reference degrees of freedom and are omitted",
    "where inference is unavailable. Zero corresponds to a unit of one."
  )
)

app_help <- function(id, fallback = NULL) {
  hit <- unname(APP_HELP[id])
  if (length(hit) && !is.na(hit) && nzchar(hit)) return(hit)
  if (!is.null(fallback) && length(fallback) && nzchar(fallback)) return(fallback)
  NULL
}
