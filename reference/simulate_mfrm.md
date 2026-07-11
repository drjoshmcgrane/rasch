# Simulate many-facet (rated) data with dial-in misfit

Generates ratings from the many-facet Rasch model (Linacre 1989): every
rater rates every person on every item, from person ability, item
difficulty, and rater severity. Departures each feed an MFRM diagnostic.

## Usage

``` r
simulate_mfrm(
  n_persons = 80,
  n_items = 5,
  n_raters = 6,
  n_categories = 4,
  theta_sd = 1.2,
  item_sd = 1,
  rater_severity_sd = 0.6,
  erratic_raters = 0,
  interaction = NULL,
  halo = 0,
  seed = NULL
)
```

## Arguments

- n_persons, n_items, n_raters:

  Facet sizes (fully crossed).

- n_categories:

  Rating categories.

- theta_sd, item_sd:

  Spread of person ability and item difficulty.

- rater_severity_sd:

  Spread of rater severities (the core facet; recovered in
  `facet_effects`).

- erratic_raters:

  Proportion of raters who rate at random (feeds the rater fit
  residual).

- interaction:

  `NULL`, or `list(rater=, item=, bias=)`: one rater is unusually harsh
  (positive) or lenient (negative) on one item. Feeds the item-by-rater
  interaction (fit with `interaction = `).

- halo:

  Proportion of raters showing a halo effect: they rate by the person's
  overall level and barely differentiate items (feeds the rater fit
  residual and the item-by-rater interaction).

- seed:

  Optional RNG seed.

## Value

A long data frame of class `"rasch_sim"` (`person`, `item`, `rater`,
`score`) ready for
[`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
with the truth attached.

## Examples

``` r
d <- simulate_mfrm(60, 5, 6, rater_severity_sd = 0.8, seed = 1)
mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                 facets = "rater")
cor(mf$facet_effects$rater$measure, attr(d, "truth")$severity)  # recovered
#> Error in cor(mf$facet_effects$rater$measure, attr(d, "truth")$severity): 'x' must be numeric
```
