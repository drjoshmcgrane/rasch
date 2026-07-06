# Package index

## Estimation

Pairwise conditional estimation of the Rasch model family.

- [`rasch()`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md) :
  Fit and diagnose a Rasch model by pairwise conditional estimation
- [`rasch_mfrm()`](https://drjoshmcgrane.github.io/rmt/reference/rasch_mfrm.md)
  : Fit a many-facet Rasch model
- [`rasch_efrm()`](https://drjoshmcgrane.github.io/rmt/reference/rasch_efrm.md)
  : Fit the extended frame of reference model
- [`pcml()`](https://drjoshmcgrane.github.io/rmt/reference/pcml.md) :
  Estimate Rasch thresholds by pairwise conditional maximum likelihood
- [`pcml_pc()`](https://drjoshmcgrane.github.io/rmt/reference/pcml_pc.md)
  : Estimate Rasch thresholds via the Andrich principal-components
  reparameterisation
- [`threshold_index()`](https://drjoshmcgrane.github.io/rmt/reference/threshold_index.md)
  : Enumerate item-category thresholds
- [`item_moments()`](https://drjoshmcgrane.github.io/rmt/reference/item_moments.md)
  : Category-score moments for a polytomous item

## Test of fit and comparison

- [`fit_summary_table()`](https://drjoshmcgrane.github.io/rmt/reference/fit_summary_table.md)
  : Test-of-fit summary as a table
- [`targeting_table()`](https://drjoshmcgrane.github.io/rmt/reference/targeting_table.md)
  : Targeting and reliability summary as a table
- [`chisq_detail()`](https://drjoshmcgrane.github.io/rmt/reference/chisq_detail.md)
  : Class-interval detail for one item's chi-square test of fit
- [`test_information()`](https://drjoshmcgrane.github.io/rmt/reference/test_information.md)
  : Test information function
- [`lr_test()`](https://drjoshmcgrane.github.io/rmt/reference/lr_test.md)
  : Likelihood-ratio test of the partial credit against the rating scale
  model
- [`compare_fits()`](https://drjoshmcgrane.github.io/rmt/reference/compare_fits.md)
  : Compare fitted Rasch models
- [`guttman_table()`](https://drjoshmcgrane.github.io/rmt/reference/guttman_table.md)
  : Guttman-ordered response matrix and reproducibility

## Persons

- [`person_wle()`](https://drjoshmcgrane.github.io/rmt/reference/person_wle.md)
  : Warm's weighted likelihood estimates by raw score
- [`person_extrapolated()`](https://drjoshmcgrane.github.io/rmt/reference/person_extrapolated.md)
  : Person measures with extrapolated extreme scores
- [`score_table()`](https://drjoshmcgrane.github.io/rmt/reference/score_table.md)
  : Raw score to measure conversion table

## Invariance and DIF

- [`dif_anova()`](https://drjoshmcgrane.github.io/rmt/reference/dif_anova.md)
  : Differential item functioning by two-way residual ANOVA
- [`dif_anova_factorial()`](https://drjoshmcgrane.github.io/rmt/reference/dif_anova_factorial.md)
  : Factorial DIF analysis with Tukey comparisons
- [`dif_contrasts()`](https://drjoshmcgrane.github.io/rmt/reference/dif_contrasts.md)
  : Planned DIF contrasts derived from the factor structure
- [`dif_size()`](https://drjoshmcgrane.github.io/rmt/reference/dif_size.md)
  : DIF magnitude in logits with pairwise comparisons
- [`split_items()`](https://drjoshmcgrane.github.io/rmt/reference/split_items.md)
  : Split items by a person factor to resolve DIF
- [`resolve_dif()`](https://drjoshmcgrane.github.io/rmt/reference/resolve_dif.md)
  : Resolve differential item functioning by iterative item splitting
- [`equate_tests()`](https://drjoshmcgrane.github.io/rmt/reference/equate_tests.md)
  : Equate two test calibrations through their common items
- [`tailored_analysis()`](https://drjoshmcgrane.github.io/rmt/reference/tailored_analysis.md)
  : Tailored analysis for guessing

## Independence and dimensionality

- [`residual_correlations()`](https://drjoshmcgrane.github.io/rmt/reference/residual_correlations.md)
  : Residual correlations for local dependence (Yen's Q3)
- [`residual_pca()`](https://drjoshmcgrane.github.io/rmt/reference/residual_pca.md)
  : Principal components of the residual correlations
- [`dimensionality_test()`](https://drjoshmcgrane.github.io/rmt/reference/dimensionality_test.md)
  : Residual-component test of unidimensionality
- [`dimensionality_magnitude()`](https://drjoshmcgrane.github.io/rmt/reference/dimensionality_magnitude.md)
  : Magnitude of multidimensionality from a subtest analysis
- [`dependence_magnitude()`](https://drjoshmcgrane.github.io/rmt/reference/dependence_magnitude.md)
  : Estimate the magnitude of response dependence between two items
- [`spread_test()`](https://drjoshmcgrane.github.io/rmt/reference/spread_test.md)
  : Spread-parameter test for dependence within subtests
- [`combine_items()`](https://drjoshmcgrane.github.io/rmt/reference/combine_items.md)
  : Combine items into subtests and re-analyse
- [`rack_data()`](https://drjoshmcgrane.github.io/rmt/reference/rack_data.md)
  [`stack_data()`](https://drjoshmcgrane.github.io/rmt/reference/rack_data.md)
  : Reshape repeated measurements for racked or stacked analysis

## Multiple choice and traditional statistics

- [`distractor_analysis()`](https://drjoshmcgrane.github.io/rmt/reference/distractor_analysis.md)
  : Distractor analysis for multiple-choice items
- [`distractor_rescore()`](https://drjoshmcgrane.github.io/rmt/reference/distractor_rescore.md)
  : Propose polytomous option scores from the distractor evidence
- [`ctt_table()`](https://drjoshmcgrane.github.io/rmt/reference/ctt_table.md)
  : Traditional (classical test theory) statistics

## Paired comparisons

- [`btl()`](https://drjoshmcgrane.github.io/rmt/reference/btl.md) : Fit
  the Bradley-Terry-Luce model to paired comparisons
- [`btl_dif()`](https://drjoshmcgrane.github.io/rmt/reference/btl_dif.md)
  : DIF analysis for paired comparisons
- [`plot_btl()`](https://drjoshmcgrane.github.io/rmt/reference/plot_btl.md)
  : Plot Bradley-Terry-Luce object locations
- [`plot_btl_categories()`](https://drjoshmcgrane.github.io/rmt/reference/plot_btl_categories.md)
  : Plot graded-comparison category curves
- [`plot_btl_icc()`](https://drjoshmcgrane.github.io/rmt/reference/plot_btl_icc.md)
  : Plot an object characteristic curve

## Plots

- [`plot_catfreq()`](https://drjoshmcgrane.github.io/rmt/reference/plot_catfreq.md)
  : Plot category frequencies
- [`plot_ccc()`](https://drjoshmcgrane.github.io/rmt/reference/plot_ccc.md)
  : Plot category probability curves
- [`plot_distractors()`](https://drjoshmcgrane.github.io/rmt/reference/plot_distractors.md)
  : Plot multiple-choice option curves
- [`plot_equate()`](https://drjoshmcgrane.github.io/rmt/reference/plot_equate.md)
  : Plot a test-equating comparison
- [`plot_facets()`](https://drjoshmcgrane.github.io/rmt/reference/plot_facets.md)
  : Plot facet severities
- [`plot_frames()`](https://drjoshmcgrane.github.io/rmt/reference/plot_frames.md)
  : Plot frame units
- [`plot_guttman()`](https://drjoshmcgrane.github.io/rmt/reference/plot_guttman.md)
  : Plot the Guttman scalogram
- [`plot_icc()`](https://drjoshmcgrane.github.io/rmt/reference/plot_icc.md)
  : Plot an item characteristic curve
- [`plot_icc_frames()`](https://drjoshmcgrane.github.io/rmt/reference/plot_icc_frames.md)
  : Plot an item's characteristic curves across frames
- [`plot_item_map()`](https://drjoshmcgrane.github.io/rmt/reference/plot_item_map.md)
  : Plot the item map (location against fit residual)
- [`plot_kidmap()`](https://drjoshmcgrane.github.io/rmt/reference/plot_kidmap.md)
  : Plot a kidmap
- [`plot_pca()`](https://drjoshmcgrane.github.io/rmt/reference/plot_pca.md)
  : Plot residual principal-component loadings
- [`plot_pcc()`](https://drjoshmcgrane.github.io/rmt/reference/plot_pcc.md)
  : Plot a person characteristic curve
- [`plot_person_fit()`](https://drjoshmcgrane.github.io/rmt/reference/plot_person_fit.md)
  : Plot person fit
- [`plot_pimap()`](https://drjoshmcgrane.github.io/rmt/reference/plot_pimap.md)
  : Plot the person-item threshold distribution
- [`plot_resid_cor()`](https://drjoshmcgrane.github.io/rmt/reference/plot_resid_cor.md)
  : Plot the residual-correlation heatmap
- [`plot_resid_dist()`](https://drjoshmcgrane.github.io/rmt/reference/plot_resid_dist.md)
  : Plot the fit residual distribution
- [`plot_scree()`](https://drjoshmcgrane.github.io/rmt/reference/plot_scree.md)
  : Scree plot of the residual components with parallel analysis
- [`plot_tcc()`](https://drjoshmcgrane.github.io/rmt/reference/plot_tcc.md)
  : Plot the test characteristic curve
- [`plot_threshold_map()`](https://drjoshmcgrane.github.io/rmt/reference/plot_threshold_map.md)
  : Plot the threshold map
- [`plot_threshold_prob()`](https://drjoshmcgrane.github.io/rmt/reference/plot_threshold_prob.md)
  : Plot threshold probability curves
- [`plot_tif()`](https://drjoshmcgrane.github.io/rmt/reference/plot_tif.md)
  : Plot the test information function
- [`plot_wright()`](https://drjoshmcgrane.github.io/rmt/reference/plot_wright.md)
  : Plot a Wright map

## Export and interface

- [`save_outputs()`](https://drjoshmcgrane.github.io/rmt/reference/save_outputs.md)
  : Save every output of a Rasch analysis to a folder
- [`save_item_plots()`](https://drjoshmcgrane.github.io/rmt/reference/save_item_plots.md)
  : Save a plot for every item
- [`save_person_plots()`](https://drjoshmcgrane.github.io/rmt/reference/save_person_plots.md)
  : Save a kidmap for every person
- [`report_html()`](https://drjoshmcgrane.github.io/rmt/reference/report_html.md)
  : Write a self-contained HTML report of a Rasch analysis
- [`run_app()`](https://drjoshmcgrane.github.io/rmt/reference/run_app.md)
  : Launch the rmt graphical interface
