# Simulate paired-comparison data

Generates dichotomous or ordered paired comparisons from the
Bradley–Terry–Luce model. Optional arguments introduce a second object
attribute, erratic judges, or within-judge dependence. Generating values
are stored in `attr(x, "truth")`.

## Usage

``` r
simulate_btl(
  n_objects = 8,
  n_judges = 12,
  reps_per_pair = 25,
  model = c("dichotomous", "polytomous", "graded"),
  n_categories = 4,
  object_sd = 1,
  second_attribute = NULL,
  erratic_judges = 0,
  dependence = NULL,
  seed = NULL,
  object_locations = NULL
)
```

## Arguments

- n_objects, n_judges:

  Objects to scale and judges comparing them. Every judge is allocated
  at least one comparison; the simulator refuses a design with fewer
  comparisons than judges.

- reps_per_pair:

  Comparisons made of each object pair.

- model:

  `"dichotomous"` (a winner) or `"polytomous"` (a rated margin in
  `n_categories` categories; an earlier development-era value `"graded"`
  is accepted as an alias).

- n_categories:

  Categories for the polytomous model.

- object_sd:

  Realised sample standard deviation of the object locations (evenly
  spaced and sum-zero).

- second_attribute:

  `NULL`, or `list(rho=)`: half the judges rank by a second object
  attribute whose realised correlation with the first is `rho`. It lies
  in \[-1, 1); at 1 the attributes are identical and no second attribute
  is planted. This introduces residual dimensionality and possible
  intransitivity.

- erratic_judges:

  Proportion of judges who choose at random. At least one judge must
  retain model-based comparisons, in each camp when a second attribute
  is generated.

- dependence:

  `NULL`, or `list(exposure=, carry_over=)`: within-judge order effects
  (a seen-before advantage and a pull from the judge's own earlier
  verdicts). Adds an `order` column. Feeds the dependence effects fitted
  by [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- seed:

  Optional non-negative whole-number RNG seed. See
  [`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md)
  for generator support.

- object_locations:

  Optional numeric vector of generated object locations. It must have
  length `n_objects`; names, when supplied, must identify the generated
  objects. Values are centred to identify the origin and take precedence
  over `object_sd`.

## Value

A data frame of class `"rasch_sim"`: `object_a`, `object_b`, `winner`
(or `response` when polytomous), `judge`, and `order` when dependence is
planted; with `attr(x, "truth")`.

## Examples

``` r
d <- simulate_btl(8, 12, erratic_judges = 0.15, seed = 1)
bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
bt$judges          # the erratic judges carry large fit residuals
#>  judge  n infit_ms outfit_ms fit_resid df_fit
#>     J1 59    0.834     0.724    -1.783 58.410
#>    J10 58    1.362     1.516     2.989 57.420
#>    J11 58    1.021     0.983    -0.117 57.420
#>    J12 58    0.902     0.828    -1.179 57.420
#>     J2 58    0.865     0.836    -1.170 57.420
#>     J3 58    0.958     0.945    -0.369 57.420
#>     J4 59    0.995     0.931    -0.471 58.410
#>     J5 58    1.371     1.528     2.592 57.420
#>     J6 58    1.018     1.057     0.400 57.420
#>     J7 59    0.796     0.876    -0.807 58.410
#>     J8 58    0.812     0.835    -1.226 57.420
#>     J9 59    1.161     1.133     0.724 58.410
```
