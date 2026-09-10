# rasch: Models and Diagnostics for Rasch Measurement Theory

Fits and evaluates models within Rasch Measurement Theory. The package
includes models for item responses, explanatory item and threshold
structures, multiple ratings, linked frames of reference, and
comparative judgement of paired stimuli, together with functions for
fit, targeting, reliability, dimensionality, local dependence,
differential item functioning, equating, and simulation. The suite
follows Rasch (1960) and Andrich and Marais (2019); the frame of
reference models follow Humphry (2005), where the model is introduced
and named, and Humphry and Andrich (2008); the dichotomous comparative
judgement model is the conditional form of the dichotomous Rasch model
(Andrich 1978), and the polytomous model is its adjacent-categories
extension (Tutz 1986).

## Rasch models

For an item with ordered scores \\x=0,\ldots,m_i\\, the partial credit
model has adjacent-category log odds
\$\$\log\\P(X\_{ni}=x)/P(X\_{ni}=x-1)\\=\theta_n-\delta\_{ix}.\$\$
Dichotomous items have one threshold. The rating scale model imposes a
common threshold structure across items. The total score is sufficient
for \\\theta_n\\, so item parameters can be estimated without specifying
a population distribution for person locations. The diagnostic functions
examine whether comparisons remain invariant across persons, items,
groups, occasions, raters, and other parts of the measurement design.

Use [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
for the dichotomous, partial credit, and rating scale models;
[`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md)
for additive item, rater, and other facet effects;
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md)
when the unit varies across linked frames; and
[`rasch_explanatory`](https://drjoshmcgrane.github.io/rasch/reference/rasch_explanatory.md)
for the linear logistic test and linear partial credit models.
[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) fits
comparative judgement models for dichotomous and polytomous paired
comparisons;
[`btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)
fits the linked-frame extension; and
[`btl_explanatory`](https://drjoshmcgrane.github.io/rasch/reference/btl_explanatory.md)
constrains object locations by observed characteristics.

## Graphical interface

[`run_app`](https://drjoshmcgrane.github.io/rasch/reference/run_app.md)
launches the package's Shiny application. It supports data import, model
fitting, diagnostics, plots, saved analysis projects and reports. The
corresponding R code is shown for each result.

## Reproducible simulations

Simulators and bootstrap procedures accept a random seed. See
[`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md)
for random-number generator support.

## References

Rasch, G. (1960). Probabilistic Models for Some Intelligence and
Attainment Tests. Copenhagen: Danish Institute for Educational Research.
(Expanded edition, 1980, Chicago: University of Chicago Press.)

Rasch, G. (1961). On general laws and the meaning of measurement in
psychology. In Proceedings of the Fourth Berkeley Symposium on
Mathematical Statistics and Probability (Vol. 4, pp. 321–334). Berkeley:
University of California Press.

Rasch, G. (1977). On specific objectivity: An attempt at formalizing the
request for generality and validity of scientific statements. Danish
Yearbook of Philosophy, 14, 58–94.

Andrich, D. (1978). Relationships between the Thurstone and Rasch
approaches to item scaling. Applied Psychological Measurement, 2(3),
451–462.

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

Fischer, G. H. (1973). The linear logistic test model as an instrument
in educational research. Acta Psychologica, 37(6), 359–374.

Fischer, G. H. and Ponocny, I. (1994). An extension of the partial
credit model with an application to the measurement of change.
Psychometrika, 59(2), 177–192.

Humphry, S. M. (2005). Maintaining a Common Arbitrary Unit in Social
Measurement. PhD thesis, Murdoch University.

Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
Rasch model. Journal of Applied Measurement, 9(3), 249–264.

Tutz, G. (1986). Bradley-Terry-Luce models with an ordered response.
Journal of Mathematical Psychology, 30(3), 306–316.

## See also

Useful links:

- <https://drjoshmcgrane.github.io/rasch/>

- <https://github.com/drjoshmcgrane/rasch>

- Report bugs at <https://github.com/drjoshmcgrane/rasch/issues>

## Author

**Maintainer**: Josh McGrane <drjoshmcgrane@gmail.com>

Authors:

- Josh McGrane <drjoshmcgrane@gmail.com>
