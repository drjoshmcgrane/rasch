# Consistency plot for paired-comparison transitivity

With `by = "judge"` (the default when judges exist), plots each judge's
consistency – one minus the circular-triad rate over the chance rate –
as a dot against the chance line at zero: the individual-judge lens, a
judge-fit analogue. With `by = "object"`, plots each object's
circular-triad involvement instead: the structural lens, showing which
objects sit in the most contradictions.

## Usage

``` r
plot_btl_transitivity(x, by = c("auto", "judge", "object"), ...)
```

## Arguments

- x:

  A `"rasch_btl_transitivity"` object.

- by:

  `"auto"` (judges if present, else objects), `"judge"`, or `"object"`.

- ...:

  Unused.

## Value

Called for its plotting side effect.
