# Transitivity of paired comparisons

The single-dimension analogue for paired comparisons of the
unidimensionality question. A Bradley-Terry-Luce scale implies that
preferences stack into one consistent order: if A beats B and B beats C
then A should beat C. A *circular triad* (A beats B, B beats C, C beats
A) is a local contradiction, like rock-paper-scissors. A few are
sampling noise; many, systematically, mean the comparisons are not being
driven by a single attribute. The rate of circular triads is compared
with the value expected from pure guessing (one quarter of triples),
and, when every pair has been compared, Kendall's coefficient of
consistency is reported (Kendall & Babington Smith 1940). With judges,
each judge's own consistency is reported too, flagging judges whose
choices approach chance.

## Usage

``` r
btl_transitivity(fit, min_triples = 5L)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- min_triples:

  A judge is reported only if this many complete triples (all three
  pairs judged) are available.

## Value

A list of class `"rasch_btl_transitivity"`: `summary` (one row: objects,
pairs compared, complete triples, circular triads, the circular rate,
the chance rate 0.25, the consistency index `1 - rate/0.25`, and
Kendall's `zeta` when the design is a complete round-robin); `objects`
(each object's circular-triad involvement); `judges` (per-judge
consistency, when judges exist); and `notes`.

## References

Kendall, M. G., & Babington Smith, B. (1940). On the method of paired
comparisons. *Biometrika*, 31(3/4), 324-345.

## Examples

``` r
set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
pr <- t(utils::combn(objs, 2))
d <- data.frame(a = rep(pr[, 1], each = 20), b = rep(pr[, 2], each = 20))
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
btl_transitivity(btl(d, "a", "b", "win"))
#> Paired-comparison transitivity: 6 objects, 20 complete triples
#> Circular triads: 0 (0.0% of triples; chance 25%) -> consistency 1.00
#> Kendall coefficient of consistency (complete design): 1.000
```
