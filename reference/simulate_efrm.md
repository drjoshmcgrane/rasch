# Simulate extended frame-of-reference data with differing units

Generates data whose latent unit differs across item-set by person-group
frames (Humphry 2005): a person in group g responding to an item in set
s does so at the frame unit rho = alpha_set \* phi_group scaling the
whole exponent. The planted set- and group-unit ratios are recovered by
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

## Usage

``` r
simulate_efrm(
  n_per_group = 300,
  items_per_set = 8,
  n_sets = 2,
  n_groups = 2,
  set_unit_ratio = 1.3,
  group_unit_ratio = 1,
  theta_sd = 1.3,
  seed = NULL
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

  Geometric span of the set and group units across their levels (1 =
  equal units, i.e. an ordinary Rasch fit).

- theta_sd:

  Spread of person ability.

- seed:

  Optional RNG seed.

## Value

A wide data frame of class `"rasch_sim"` (`id`, item columns, `group`)
with `attr(x, "truth")$item_sets` the set map to pass to
[`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

## Examples

``` r
d <- simulate_efrm(300, 8, set_unit_ratio = 1.3, seed = 1)
tr <- attr(d, "truth")
ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group")
ef$alpha_table   # recovers the ~1.3 set-unit ratio
#>    set     alpha se_log_alpha
#> 1 set1 0.8913905   0.04276361
#> 2 set2 1.1218428   0.04276361
```
