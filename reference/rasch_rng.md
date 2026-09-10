# Random-number generation

Supplying a `seed` makes a simulation or bootstrap reproducible and
restores the caller's random-number stream on exit. Bootstrap methods
that assign seeds to individual replicates also restore those local
streams.

## Details

These operations do not support R's Box–Muller normal generator: its
cached normal value is not part of `.Random.seed`, so restoring that
vector would change subsequent draws. They refuse before changing the
stream. Direct `simulate_*` calls with `seed = NULL` can still use
Box–Muller.
[`sim_replicate`](https://drjoshmcgrane.github.io/rasch/reference/sim_replicate.md)
assigns replicate seeds even when its own `seed` is `NULL`.

The default Inversion generator is supported. To select it explicitly,
use `RNGkind(normal.kind = "Inversion")` before setting the seed for the
analysis. Changing the generator starts a different normal stream; it
does not recover a previous Box–Muller stream.

## See also

[`Random`](https://rdrr.io/r/base/Random.html),
[`simulate_rasch`](https://drjoshmcgrane.github.io/rasch/reference/simulate_rasch.md),
[`fit_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/fit_bootstrap.md),
[`dif_bootstrap`](https://drjoshmcgrane.github.io/rasch/reference/dif_bootstrap.md).
