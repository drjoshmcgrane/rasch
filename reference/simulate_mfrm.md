# Simulate many-facet Rasch data

Generates fully crossed ratings from a many-facet Rasch model (Linacre
1989), with optional erratic raters, item-by-rater interaction, or halo.
A positive rater proportion selects at least one rater when the
requested departures can coexist.

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

- theta_sd:

  Realised sample standard deviation of person ability.

- item_sd:

  Realised sample standard deviation of the deterministic item
  difficulties; zero gives equal item difficulties.

- rater_severity_sd:

  Realised sample standard deviation of rater severities (the core
  facet; recovered in `facet_effects`).

- erratic_raters:

  Proportion of raters who rate at random (feeds the rater fit
  residual). Erratic and halo raters are disjoint, and together must
  leave at least one ordinary rater.

- interaction:

  `NULL`, or `list(rater=, item=, bias=)`: one rater is unusually harsh
  (positive) or lenient (negative) on one item. Feeds the item-by-rater
  interaction (fit with `interaction = `).

- halo:

  Proportion of raters showing a halo effect: they rate by the person's
  overall level and barely differentiate items (feeds the rater fit
  residual and the item-by-rater interaction). Its requested count must
  fit among the non-erratic raters and leave an ordinary rater. A
  positive halo proportion requires `item_sd > 0`.

- seed:

  Optional non-negative whole-number RNG seed. See
  [`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md)
  for generator support.

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
cor(mf$facet_effects$rater$severity, attr(d, "truth")$severity)  # recovered
#> [1] 0.9911336
```
