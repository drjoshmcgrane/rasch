# Simulate extended frame-of-reference data with differing units

Generates data whose latent unit differs across item-set by person-group
frames (Humphry 2005): a person in group g responding to an item in set
s does so at the frame unit rho = alpha_set \* phi_group scaling the
whole exponent. The planted set- and group-unit ratios are recovered by
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).
A positive careless-response or missingness proportion selects at least
one person or response cell.

## Usage

``` r
simulate_efrm(
  n_per_group = 300,
  items_per_set = 8,
  n_sets = 2,
  n_groups = 2,
  set_unit_ratio = 1.3,
  group_unit_ratio = 1,
  n_categories = 2,
  theta_sd = 1.3,
  seed = NULL,
  item_drift = NULL,
  careless = 0,
  missing = 0
)
```

## Arguments

- n_per_group:

  Persons in each group.

- items_per_set:

  Items in each set.

- n_sets, n_groups:

  Numbers of item sets and person groups.

- set_unit_ratio, group_unit_ratio:

  Ratio of the last generated set or group unit to the first.
  Intermediate units are geometrically spaced; 1 gives equal units and
  hence an ordinary Rasch fit.

- n_categories:

  Response categories per item: 2 (the default) gives dichotomous items;
  larger values give partial credit items whose evenly spaced thresholds
  are centred on the item locations, with the frame unit scaling the
  whole exponent as in the dichotomous case.

- theta_sd:

  Realised sample standard deviation of person ability. This must be
  positive when more than one item set is generated because relative set
  units are identified from person variation shared across sets.

- seed:

  Optional non-negative whole-number RNG seed. See
  [`rasch_rng`](https://drjoshmcgrane.github.io/rasch/reference/rasch_rng.md)
  for generator support.

- item_drift:

  Optional `list(items=, group=, shift=)`. The named item or items move
  by `shift` logits in one generated person group, violating item
  invariance across frames. At least one item in every affected set must
  remain invariant: shifting a set's items together is a frame-origin
  difference, not item drift within that set. A non-zero drift requires
  at least two person groups.

- careless:

  Proportion of persons whose complete response vectors are replaced by
  random category choices. At least one person in every group must
  retain model-based responses.

- missing:

  Proportion of response cells set missing completely at random after
  the responses are generated. It must leave at least one observed
  response.

## Value

A wide data frame of class `"rasch_sim"`, containing an ID, item
columns, and group. Its truth attribute contains the item-set map
required by
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

## Examples

``` r
d <- simulate_efrm(200, 6, set_unit_ratio = 1.3, seed = 1)
tr <- attr(d, "truth")
ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                 id = "id", boot_reps = 0)    # point estimates only
ef$alpha_table   # planted ratio 1.3, recovered within small-sample noise
#>   set alpha se_log_alpha
#>  set1 0.850             
#>  set2 1.176             
```
