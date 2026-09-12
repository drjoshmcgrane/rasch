# Compare fitted and generating parameters

Compares fitted parameters with the generating values from a
`simulate_*` function (carried on the data as `attr(sim, "truth")`):
item difficulties and person abilities for a Rasch fit, object locations
for a paired-comparison fit, rater severities (with item and person
measures) for a many-facet fit, set and group units for a Rasch frames
fit, and common object locations, panel and set units, and set origins
for a paired-comparison frames fit. Locations are mean-centred where the
model identifies them only up to an origin. An externally anchored Rasch
or paired-comparison fit retains its identified origin, so recovery and
bias are reported on the anchored scale. The fit must be from these
simulated responses and model family, and must have converged. Row and
item order do not matter; response values and the person, judge, facet
and frame allocations do. A fit given no `id` column labels its persons
by row number, which are not identifiers: the responses are still
verified, but person recovery and the EFRM group units are withheld with
a note rather than compared against an unverified person allocation; the
EFRM set units are still reported, with the note recording that the
fitted group allocation is unchecked. Pass `id =` when fitting to
recover them. An item the estimator dropped, such as one everyone
answered identically, is named in the note and left out of the
comparison. Recovery is unavailable when fitting removes or merges
generating response categories, because the fitted locations then
describe a different scale. For a many-facet simulation, the planted
rater facet must be identifiable uniquely by its name or level labels.
EFRM set parameters are matched by their item or object membership, not
by the spelling of the set labels; a different fitted partition is
refused. Paired-comparison frames align generating origins to the fitted
reference set. That set must have a generating unit of one: fixing
another unit to one imposes a different cross-set restriction, so
recovery is refused. A common generating discrimination changes only the
logit unit, so ordinary Rasch recovery uses the equivalently rescaled
item thresholds and person locations. When the generator includes a
departure that the fitted model does not represent, the comparisons are
labelled as descriptive rather than as recovery of a single correctly
specified target. This includes, for example, heterogeneous item
discriminations in an ordinary Rasch fit.

## Usage

``` r
sim_recovery(fit, sim)
```

## Arguments

- fit:

  A fit of the simulated data
  ([`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md),
  [`rasch_mfrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_mfrm.md),
  [`rasch_efrm`](https://drjoshmcgrane.github.io/rasch/reference/rasch_efrm.md),
  or
  [`btl_efrm`](https://drjoshmcgrane.github.io/rasch/reference/btl_efrm.md)).

- sim:

  The simulated data (from a `simulate_*` function).

## Value

A list of class `"rasch_recovery"`: `summary` (per parameter type: n,
correlation, RMSE, bias) and `pieces` (the true and estimated values
behind each). `note` identifies unrepresented generating departures, an
unverifiable original response scale, generated items the fit does not
estimate, and comparisons withheld because the fit carries no person
identifiers.

## Examples

``` r
d <- simulate_rasch(500, 12, seed = 1)
sim_recovery(rasch(d, id = "id"), d)$summary
#>         parameter   n correlation      rmse bias
#> 1 item difficulty  12   0.9980002 0.1703475   NA
#> 2  person ability 500   0.8177074 0.6756437   NA
```
