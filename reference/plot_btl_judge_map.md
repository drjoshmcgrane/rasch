# Unexpected-judgement map for one judge (pair level)

The judge counterpart of the kidmap, drawn matchup by matchup. Each pair
the judge met is a segment on the consensus location axis, spanning its
two objects, positioned horizontally by its standardised residual. For a
non-tied pair, a negative residual means the stronger object
under-performed. A tied pair has no directional interpretation and is
oriented to a non-negative residual for display only. A filled dot marks
the object the judge's verdict favoured and a hollow dot the other. The
rug marks every object's location.

## Usage

``` r
plot_btl_judge_map(fit, judge, min_n = 1L, flag_z = 1.96, ...)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md) with
  judges.

- judge:

  The judge to map.

- min_n, flag_z:

  Passed to
  [`judge_pair_surprise`](https://drjoshmcgrane.github.io/rasch/reference/judge_pair_surprise.md).

- ...:

  Unused.

## Value

Called for its plotting side effect; invisibly the
`rasch_btl_judge_pairs` object.

## Examples

``` r
# \donttest{
d <- simulate_btl(6, 10, reps_per_pair = 20, seed = 1)
fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
plot_btl_judge_map(fit, judge = "J1")

# }
```
