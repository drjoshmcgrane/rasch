# Package index

## Estimation

Fit Rasch models and estimate item and person parameters.

- [`rasch()`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  : Fit a Rasch model
- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
  : Fit a many-facet Rasch model
- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
  : Fit the extended frame of reference model
- [`rasch_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/rasch_explanatory.md)
  : Fit an explanatory Rasch model
- [`btl_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/btl_explanatory.md)
  : Fit an explanatory comparative judgement model
- [`explanatory_test()`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_test.md)
  : Compare an explanatory model with its free calibration
- [`explanatory_diagnostics()`](https://drjoshmcgrane.github.io/rasch/reference/explanatory_diagnostics.md)
  : Diagnose fixed departures from an explanatory model
- [`relax_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/relax_explanatory.md)
  : Relax a nominated explanatory restriction
- [`relax_btl_explanatory()`](https://drjoshmcgrane.github.io/rasch/reference/relax_btl_explanatory.md)
  : Add a fixed object departure to an explanatory comparative judgement
  model
- [`pcml()`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md) :
  Estimate Rasch thresholds by pairwise conditional maximum likelihood
- [`pcml_pc()`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md)
  : Estimate Rasch thresholds using a principal-component
  parameterisation
- [`threshold_index()`](https://drjoshmcgrane.github.io/rasch/reference/threshold_index.md)
  : Enumerate item-category thresholds
- [`item_moments()`](https://drjoshmcgrane.github.io/rasch/reference/item_moments.md)
  : Category-score moments for a polytomous item

## Test of fit and comparison

Examine fit, targeting and reliability, or compare fitted models.

- [`fit_summary_table()`](https://drjoshmcgrane.github.io/rasch/reference/fit_summary_table.md)
  : Test-of-fit summary as a table
- [`targeting_table()`](https://drjoshmcgrane.github.io/rasch/reference/targeting_table.md)
  : Targeting and reliability summary as a table
- [`chisq_detail()`](https://drjoshmcgrane.github.io/rasch/reference/chisq_detail.md)
  : Class-interval detail for one item's chi-square test of fit
- [`fit_bootstrap()`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md)
  : Bootstrap fit statistics
- [`test_information()`](https://drjoshmcgrane.github.io/rasch/reference/test_information.md)
  : Test information function
- [`lr_test()`](https://drjoshmcgrane.github.io/rasch/reference/lr_test.md)
  : Compare the partial credit and rating scale models
- [`compare_fits()`](https://drjoshmcgrane.github.io/rasch/reference/compare_fits.md)
  : Compare fitted Rasch models
- [`guttman_table()`](https://drjoshmcgrane.github.io/rasch/reference/guttman_table.md)
  : Guttman-ordered response matrix and reproducibility

## Persons

- [`person_wle()`](https://drjoshmcgrane.github.io/rasch/reference/person_wle.md)
  : Warm's weighted likelihood estimates by raw score
- [`weighted_person_estimates()`](https://drjoshmcgrane.github.io/rasch/reference/weighted_person_estimates.md)
  : Person estimates with externally imposed weights
- [`person_extrapolated()`](https://drjoshmcgrane.github.io/rasch/reference/person_extrapolated.md)
  : Person measures with extrapolated extreme scores
- [`score_table()`](https://drjoshmcgrane.github.io/rasch/reference/score_table.md)
  : Raw score to measure conversion table

## Invariance and DIF

Examine DIF across groups or occasions and equate calibrations.

- [`dif_anova()`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  : Differential item functioning by residual analysis of variance
- [`dif_bootstrap()`](https://drjoshmcgrane.github.io/rasch/reference/dif_bootstrap.md)
  : Bootstrap sensitivity analysis for DIF
- [`dif_contrasts()`](https://drjoshmcgrane.github.io/rasch/reference/dif_contrasts.md)
  : Planned DIF contrasts
- [`dif_posthoc()`](https://drjoshmcgrane.github.io/rasch/reference/dif_posthoc.md)
  : Pairwise follow-up comparisons for a DIF term
- [`dif_size()`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
  : DIF differences between factor levels
- [`split_items()`](https://drjoshmcgrane.github.io/rasch/reference/split_items.md)
  : Split items by a person factor to resolve DIF
- [`drop_items()`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
  : Drop items and refit
- [`resolve_frames()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_frames.md)
  : Resolve items that do not hold across frames
- [`frame_invariance()`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md)
  : Test item invariance across frames
- [`resolve_dif()`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
  : Resolve differential item functioning by iterative item splitting
- [`equate_tests()`](https://drjoshmcgrane.github.io/rasch/reference/equate_tests.md)
  : Equate two test calibrations through their common items
- [`tailored_analysis()`](https://drjoshmcgrane.github.io/rasch/reference/tailored_analysis.md)
  : Tailored analysis for guessing

## Independence and dimensionality

Examine dimensionality and local response dependence.

- [`residual_correlations()`](https://drjoshmcgrane.github.io/rasch/reference/residual_correlations.md)
  : Residual correlations for local dependence (Yen's Q3)
- [`residual_pca()`](https://drjoshmcgrane.github.io/rasch/reference/residual_pca.md)
  : Principal components of the residual correlations
- [`dimensionality_test()`](https://drjoshmcgrane.github.io/rasch/reference/dimensionality_test.md)
  : Residual-component test of unidimensionality
- [`dimensionality_magnitude()`](https://drjoshmcgrane.github.io/rasch/reference/dimensionality_magnitude.md)
  : Magnitude of multidimensionality from a subtest analysis
- [`dependence_magnitude()`](https://drjoshmcgrane.github.io/rasch/reference/dependence_magnitude.md)
  : Estimate the magnitude of response dependence between two items
- [`spread_test()`](https://drjoshmcgrane.github.io/rasch/reference/spread_test.md)
  : Spread-parameter test for dependence within subtests
- [`combine_items()`](https://drjoshmcgrane.github.io/rasch/reference/combine_items.md)
  : Combine items into subtests and re-analyse
- [`rack_data()`](https://drjoshmcgrane.github.io/rasch/reference/rack_data.md)
  [`stack_data()`](https://drjoshmcgrane.github.io/rasch/reference/rack_data.md)
  : Reshape repeated measurements for racked or stacked analysis

## Multiple choice and traditional statistics

- [`distractor_analysis()`](https://drjoshmcgrane.github.io/rasch/reference/distractor_analysis.md)
  : Distractor analysis for multiple-choice items
- [`distractor_rescore()`](https://drjoshmcgrane.github.io/rasch/reference/distractor_rescore.md)
  : Propose polytomous option scores from the distractor evidence
- [`ctt_table()`](https://drjoshmcgrane.github.io/rasch/reference/ctt_table.md)
  : Traditional (classical test theory) statistics

## Paired comparisons

Fit and examine models for dichotomous or ordered paired comparisons.

- [`btl()`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) :
  Fit comparative judgement models to paired comparisons
- [`btl_dif()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dif.md)
  : DIF analysis for paired comparisons
- [`btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
  : Fit the extended frame of reference model for paired comparisons
- [`btl_transitivity()`](https://drjoshmcgrane.github.io/rasch/reference/btl_transitivity.md)
  : Transitivity of paired comparisons
- [`btl_dimensionality()`](https://drjoshmcgrane.github.io/rasch/reference/btl_dimensionality.md)
  : Experimental residual dimensionality of paired comparisons
- [`btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/btl_equate.md)
  : Equate two paired-comparison calibrations through their common
  objects
- [`btl_information()`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md)
  : Information and targeting of a paired-comparison design
- [`btl_next_pairs()`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)
  : Recommend the next informative comparisons (adaptive step)
- [`judge_surprise()`](https://drjoshmcgrane.github.io/rasch/reference/judge_surprise.md)
  : Unexpected judgements of one judge
- [`judge_pair_surprise()`](https://drjoshmcgrane.github.io/rasch/reference/judge_pair_surprise.md)
  : Unexpected judgements of one judge, pair by pair
- [`plot_btl()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl.md)
  : Plot Bradley-Terry-Luce object locations
- [`plot_btl_categories()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_categories.md)
  : Plot polytomous-comparison category curves
- [`plot_btl_icc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_icc.md)
  : Plot an object characteristic curve
- [`plot_btl_dependence()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dependence.md)
  : Plot a within-judge dependence effect
- [`plot_btl_transitivity()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_transitivity.md)
  : Consistency plot for paired-comparison transitivity
- [`plot_btl_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_scree.md)
  : Scree of paired-comparison residual bimensions
- [`plot_btl_dim_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_dim_map.md)
  : Residual map of the leading paired-comparison bimension
- [`plot_btl_judge_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_judge_map.md)
  : Unexpected-judgement map for one judge (pair level)
- [`plot_btl_units()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_units.md)
  : Plot the frame units of a paired-comparison EFRM fit
- [`plot_btl_equate()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_equate.md)
  : Plot a paired-comparison equating comparison
- [`plot_btl_targeting()`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_targeting.md)
  : Targeting plot for a paired-comparison design

## Plots

Plot fitted models and diagnostic results.

- [`wright_map()`](https://drjoshmcgrane.github.io/rasch/reference/wright_map.md)
  : Draw a Wright map with WrightMap
- [`plot_catfreq()`](https://drjoshmcgrane.github.io/rasch/reference/plot_catfreq.md)
  : Plot category frequencies
- [`plot_ccc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_ccc.md)
  : Plot category probability curves
- [`plot_distractors()`](https://drjoshmcgrane.github.io/rasch/reference/plot_distractors.md)
  : Plot multiple-choice option curves
- [`plot_equate()`](https://drjoshmcgrane.github.io/rasch/reference/plot_equate.md)
  : Plot a test-equating comparison
- [`plot_facets()`](https://drjoshmcgrane.github.io/rasch/reference/plot_facets.md)
  : Plot facet severities
- [`plot_frames()`](https://drjoshmcgrane.github.io/rasch/reference/plot_frames.md)
  : Plot frame units
- [`plot_guttman()`](https://drjoshmcgrane.github.io/rasch/reference/plot_guttman.md)
  : Plot the Guttman scalogram
- [`plot_icc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_icc.md)
  : Plot an item characteristic curve
- [`plot_icc_frames()`](https://drjoshmcgrane.github.io/rasch/reference/plot_icc_frames.md)
  : Plot an item's characteristic curves across frames
- [`plot_item_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_item_map.md)
  : Plot the item map (location against fit residual)
- [`plot_kidmap()`](https://drjoshmcgrane.github.io/rasch/reference/plot_kidmap.md)
  : Plot a kidmap
- [`plot_pca()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pca.md)
  : Plot residual principal-component loadings
- [`plot_pca_biplot()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pca_biplot.md)
  : Biplot of the first two residual components
- [`plot_pcc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pcc.md)
  : Plot a person characteristic curve
- [`plot_person_fit()`](https://drjoshmcgrane.github.io/rasch/reference/plot_person_fit.md)
  : Plot person fit
- [`plot_pimap()`](https://drjoshmcgrane.github.io/rasch/reference/plot_pimap.md)
  : Plot the person-item threshold distribution
- [`plot_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/plot_recovery.md)
  : Plot fitted against generating parameters
- [`plot_resid_cor()`](https://drjoshmcgrane.github.io/rasch/reference/plot_resid_cor.md)
  : Plot the residual-correlation heatmap
- [`plot_resid_dist()`](https://drjoshmcgrane.github.io/rasch/reference/plot_resid_dist.md)
  : Plot the fit residual distribution
- [`plot_scree()`](https://drjoshmcgrane.github.io/rasch/reference/plot_scree.md)
  : Scree plot of the residual components with parallel analysis
- [`plot_tcc()`](https://drjoshmcgrane.github.io/rasch/reference/plot_tcc.md)
  : Plot the test characteristic curve
- [`plot_threshold_map()`](https://drjoshmcgrane.github.io/rasch/reference/plot_threshold_map.md)
  : Plot the threshold map
- [`plot_threshold_prob()`](https://drjoshmcgrane.github.io/rasch/reference/plot_threshold_prob.md)
  : Plot threshold probability curves
- [`plot_tif()`](https://drjoshmcgrane.github.io/rasch/reference/plot_tif.md)
  : Plot the test information function
- [`plot_wright()`](https://drjoshmcgrane.github.io/rasch/reference/plot_wright.md)
  : Plot a Wright map

## Data simulation

Simulate data from the fitted model families.

- [`simulate_rasch()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md)
  : Simulate person-by-item Rasch data
- [`simulate_btl()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl.md)
  : Simulate paired-comparison data
- [`simulate_mfrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_mfrm.md)
  : Simulate many-facet Rasch data
- [`simulate_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_efrm.md)
  : Simulate extended frame-of-reference data with differing units
- [`simulate_btl_efrm()`](https://drjoshmcgrane.github.io/rasch/reference/simulate_btl_efrm.md)
  : Simulate paired-comparison EFRM data with differing frame units
- [`sim_replicate()`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md)
  : Replicate a simulation for Monte Carlo studies
- [`sim_apply()`](https://drjoshmcgrane.github.io/rasch/reference/sim_apply.md)
  : Apply a statistic across a simulation batch
- [`sim_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/sim_recovery.md)
  : Compare fitted and generating parameters
- [`plot_recovery()`](https://drjoshmcgrane.github.io/rasch/reference/plot_recovery.md)
  : Plot fitted against generating parameters
- [`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md)
  : Random-number generation

## Export and interface

Export results or launch the Shiny application.

- [`save_outputs()`](https://drjoshmcgrane.github.io/rasch/reference/save_outputs.md)
  : Save the outputs of a Rasch analysis
- [`save_item_plots()`](https://drjoshmcgrane.github.io/rasch/reference/save_item_plots.md)
  : Save a plot for every item
- [`save_person_plots()`](https://drjoshmcgrane.github.io/rasch/reference/save_person_plots.md)
  : Save a kidmap for every person
- [`report_html()`](https://drjoshmcgrane.github.io/rasch/reference/report_html.md)
  : Write a self-contained HTML report of a Rasch analysis
- [`report_document()`](https://drjoshmcgrane.github.io/rasch/reference/report_document.md)
  : Write an editable or print-ready analysis report
- [`run_app()`](https://drjoshmcgrane.github.io/rasch/reference/run_app.md)
  : Launch the rasch point-and-click graphical interface
