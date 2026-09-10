# Resolve items that do not hold across frames

Gives each named item a separate location in every frame in which it was
administered, then refits the EFRM.

## Usage

``` r
resolve_frames(fit, items, boot_reps = NULL)
```

## Arguments

- fit:

  A fitted object from
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md).

- items:

  Item names to resolve.

- boot_reps:

  Bootstrap replicates for the refit. The default retains the fitted
  specification; a number overrides it.

## Value

A refitted object of class `"rasch_efrm"`, carrying a note for each item
resolved. The resolved versions appear in the item table as
`"item (frame)"`.

## Details

A resolved item continues to contribute to person measurement within
each frame but no longer constrains the link between those frames. Its
versions are named `"item (frame)"`. The remaining common items and the
linked set design must still identify the frame units; otherwise the
refit is refused by the model's connectivity and rank checks. Each set
must also retain links between its groups of item versions to identify
their relative origins. A unit link supplied by another set cannot
replace these origin links.

Resolve an item when its within-frame measurement remains defensible but
its cross-frame location does not. This refit does not estimate a
separate discrimination and therefore does not resolve a
discrimination-only flag from
[`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md).
Review or remove such an item instead. Use
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
when the item should no longer contribute to measurement. Each resolved
version must observe every score category of the source item. If it does
not, the refit is refused rather than renumbering that frame's scores.

## See also

[`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md),
which identifies the items to resolve;
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md),
which removes an item instead; and
[`split_items`](https://drjoshmcgrane.github.io/rasch/reference/split_items.md),
the equivalent for an ordinary fit.

## Examples

``` r
d <- simulate_efrm(n_per_group = 200, items_per_set = 6, n_sets = 2,
                   n_groups = 2, set_unit_ratio = 1.3, seed = 5)
tr <- attr(d, "truth")
fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                  id = "id", boot_reps = 0)
fit2 <- resolve_frames(fit, "S1I02", boot_reps = 0)
grep("S1I02", fit2$items$item, value = TRUE)
#> [1] "S1I02 (g1):g1" "S1I02 (g2):g2"
```
