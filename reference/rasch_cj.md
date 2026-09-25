# Calibrate items from responses and comparative judgements together

Fits one set of item thresholds to persons' responses and to judges'
comparisons or rankings of the same items, giving each judgement frame
its own unit relative to the response frame, and tests whether the
frames agree about the items.

## Usage

``` r
rasch_cj(
  data,
  comparisons = NULL,
  object_a = "object_a",
  object_b = "object_b",
  winner = "winner",
  threshold_a = "threshold_a",
  threshold_b = "threshold_b",
  rankings = NULL,
  ranking = "ranking",
  item = "item",
  rank = "rank",
  threshold = "threshold",
  units = c(comparisons = NA, rankings = NA),
  items = NULL,
  na_codes = -1,
  objects = c("items", "persons"),
  anchors = NULL,
  id = NULL,
  maxit = 200,
  tol = 1e-08
)
```

## Arguments

- data:

  Persons-by-items response data, dichotomous or polytomous, as for
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
  or `NULL` to combine comparisons and rankings without responses. In
  the person mode, responses to the anchored items, coded from 0 to the
  item's top category. A named list of such tables gives several tests,
  the first the reference; the judgements must link them.

- comparisons:

  Optional data frame of paired comparisons, one row each.

- object_a, object_b, winner:

  Column names in `comparisons`, as in
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md). The
  winner is the object judged to have the higher location, for
  difficulty the harder one; it is given as the item name, or as
  `"item:k"` when both sides are thresholds of the same item. Winners
  equal to `"tie"` or `"draw"` are dropped.

- threshold_a, threshold_b:

  Optional column names in `comparisons` giving the threshold number of
  each side, `NA` for the item's location. Used only when the columns
  are present.

- rankings:

  Optional data frame of rankings in long form, one row per object
  placed: a ranking identifier, the item, and its rank within that
  ranking. Rankings may cover any subset of two or more objects.

- ranking, item, rank:

  Column names in `rankings`. Rank 1 is the object with the highest
  location.

- threshold:

  Optional column name in `rankings` giving the threshold number of each
  row, `NA` for the item's location. Used only when the column is
  present.

- units:

  Named vector giving the unit of the `comparisons` and `rankings`
  frames, and of each test after the first when `data` is a list: `NA`
  (the default) to estimate, or `1` to fix at the reference unit.

- items:

  Optional item columns to analyse, as in
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- na_codes:

  As in
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- objects:

  What the judges compare: `"items"` (the default), or `"persons"` to
  measure the persons from their responses and judgements of their work.

- anchors:

  Person mode only: the calibrated item thresholds, as a
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  fit, an item-mode `rasch_cj` fit, or a data frame with columns `item`,
  `k` and `tau` giving every threshold of every item in `data`; with
  several tests, one such object covering every test's items, or a list
  of them named by test.

- id:

  Person mode only: the name of a column of `data` holding the person
  identifiers the judgement tables use, or a vector of them, one per
  row; with several tests, a column name found in each table or a list
  of vectors named by test. Defaults to the row names.

- maxit, tol:

  Newton iteration cap and convergence tolerance on the parameter scale.

## Value

An object of class `"rasch_cj"` with components `items` (item, location,
se, and the separate calibration of each item-level object from each
frame on the reference scale), `thresholds` (item, k, threshold, se),
`objects` (the judged objects and the frames reaching them), `units`
(frame, unit, se, and whether it was estimated), `invariance` (a list
with the likelihood ratio test and the per-object table comparing each
judgement frame with the reference frame), `anchors` (a data frame ready
for
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)),
`cov` (covariance of the thresholds), `cov_items` (covariance of the
item locations), `loglik`, `converged`, `iterations`, frame sizes in
`n`, the `reference` frame, the `tests` named in `data`, and `notes`.
The per-object invariance table has an `against` column naming the frame
each row is compared with.

In the person mode, `mode` is `"persons"` and the object holds `persons`
(person, n_items, raw, max_raw, the combined location and se, the
location and se from responses alone, `Inf` or `-Inf` at an extreme
score, the centred location from each judgement frame alone on the test
scale, whether each frame reaches the person, and whether the person was
set aside as extreme), `units`, `tests` (with several tests: test, unit,
se_unit, shift, se_shift), `invariance` (a list with the per-person
table `persons`), `anchors` (the thresholds used, by test when there are
several), `cov` (covariance of the estimated locations), `loglik`,
`converged`, `iterations`, `n` and `notes`.

## Details

The response frame contributes the conditional likelihood of each
person's pattern given the raw score under the partial credit model, so
person parameters are eliminated. Unlike the pairwise likelihood in
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
this conditions on the whole response score. A dichotomous item has one
threshold, its location. Comparisons contribute the Bradley-Terry
likelihood \\P(a \succ b) = \mathrm{logistic}\\\alpha (\lambda_a -
\lambda_b)\\\\, in which a judge's own location cancels. Rankings
contribute the Plackett-Luce likelihood, a ranking of \\n\\ objects
being \\n - 1\\ sequential choices with \\P(\text{next} = a) \propto
\exp(\kappa \lambda_a)\\ over the objects still unplaced; with \\n = 2\\
it is Bradley-Terry. The three blocks share the thresholds and are
summed, so the estimator is a full likelihood: standard errors are from
the inverse observed information and no weighting of the sources is
required. This assumes independent response rows, independent
comparisons or rankings within each source, and independent sources.
Repeated judgements from the same judge are not cluster-adjusted here.

What a judge compares is an object. By default an object is an item and
its location \\\lambda\\ is the mean of the item's thresholds, so
judgements inform the location of a polytomous item and leave the
spacing of its thresholds to the responses. A threshold column in the
judgement data makes the object a single threshold of the item,
\\\lambda\\ being that threshold: a judge then says that reaching
category \\k\\ of one item is harder than reaching category \\l\\ of
another. Both kinds of object may appear in one source, and two
thresholds of one item may be compared with each other.

\\\alpha\\ and \\\kappa\\ are the units of the judgement frames relative
to the response frame (Humphry and Andrich 2008). A unit of 0.5 means
the judges separate the objects at half the discrimination the test
does. They are estimated on the log scale; `units` fixes either at one
instead, which asserts that the frame shares the test's unit.

The model assumes each object has one location across the frames it
appears in. `invariance` tests that assumption in two ways. The
likelihood ratio compares the fitted model with separate locations per
frame (under which the units are absorbed into the locations and drop
out). Each separate calibration has one origin constraint per connected
block. The degrees of freedom are the sum of their free parameter counts
minus the joint model's count, including its estimated units. The object
table places each frame's separate calibration on the reference scale by
dividing by the fitted unit, centres it within each block, and tests
each object's difference from the reference calibration by a Wald
statistic. Its covariance includes uncertainty in the fitted units and
their shared-data covariance with the separate calibrations, using the
joint information and the delta method. Probabilities are Holm-adjusted
across objects within each frame comparison. An object that fails is one
the judges and the test-takers disagree about, and combining the sources
moves its estimate toward whichever source has the more information.

Items are calibrated when at least one frame reaches them, and a
disconnected item graph is refused. An item every person answered the
same way carries no conditional information; its responses are left out
and it is located, as a single threshold, by the judgements, which must
reach it.

The fitted thresholds are returned as an anchor table in the `anchors`
component. Passing it to
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md) runs
item and person fit, targeting and the other response diagnostics on the
combined calibration; `rasch` keeps at least one item free, so pass the
rows of every item but one and that item is re-estimated on the anchored
scale.

With `data = NULL` and both judgement sources supplied, the comparisons
become the reference frame and the rankings' unit is relative to theirs;
the objects are the items the two sources name, and threshold columns
are refused because no responses define the thresholds.

Ties are dropped with a note; judge clustering and judge fit are not
provided.

**Several tests.** Two tests of one construct with no item in common
cannot be equated from their responses, and nothing says their responses
are in the same unit. Give `data` as a named list, one response table
per test, and the judgements link them: each test is its own conditional
block with its own unit relative to the first, the reference, whose unit
is 1. The judgements must reach items of every test, since only they
connect the tests. The origin of each block is free in the conditional
likelihood, so the tests differ by a unit only, and the joint fit places
all the items on the reference scale. An item in more than one test is
one item, and must have the same categories in each. The invariance
table compares each test after the first with each judgement frame, and
each judgement frame with the reference test, and the likelihood ratio
test counts the free locations of every frame against the joint
parameters as before. Fix a test's unit at 1 by naming it in `units`.

**Measuring persons.** With `objects = "persons"` the judges compare the
persons' work rather than the items, and the function locates each
person from their responses and from those judgements together. The
items must be calibrated already: `anchors` gives their thresholds, from
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
from an item-mode `rasch_cj` fit, or as a data frame. The response block
is then the partial credit likelihood of each person's responses given
the anchored thresholds, one location per person, and the judgement
blocks are as above with the persons as objects and units relative to
the test. The items are held, so a person's information grows with their
own items and judgements, and the standard errors are from the observed
information. The judgement tables name persons by the row names of
`data`, by the `id` column or vector, or by `"P1"`, `"P2"`, ... when
there are none; a person judged but absent from `data` has no responses
and is placed by the judgements. Threshold columns have no meaning here
and must be omitted or missing.

A person has a finite location when the evidence points both ways: a
score inside its range, or at least one win and one loss (a rank above
someone and a rank below someone) among the persons being estimated. A
person lacking either is set aside as extreme and reported at the Warm
estimate from their responses, as `rasch` reports extreme persons, or
with no location if they have no responses; a judged group none of whose
members has responses is not on the test scale and is refused. The
invariance table compares each frame's separate placement of every
person, divided by the fitted unit and centred within each connected
block of the design over the persons both frames locate, with the
response placement, by a Wald contrast with Holm adjustment. The
contrast covariance includes estimation of the test shifts and all
linking units. The anchored item calibration is treated as fixed. There
is no likelihood ratio test in this mode: the separate model has a
location per person per frame, so its parameters grow with the persons
and the ratio is not chi-square. A person who fails is one whose judged
work does not match their responses.

The unit of a judgement frame is estimated alongside one location per
person, so it carries the incidental-parameter bias of joint estimation
rather than the consistency of the item mode, where persons are
conditioned out. In simulations with 300 persons, 8 to 40 dichotomous
items and 4 to 20 comparisons per person the unit was 5 to 10 percent
high and its standard error understated, while the person locations and
their standard errors were calibrated. Fix the unit with `units` when it
is known. Too few judgements per person leave the unit without a
maximum, because some ordering of the persons the responses allow agrees
with every judgement; the fit then reports a runaway unit.

Persons from several tests are measured together in the same way: a
named list of response tables, each with its own anchors (one anchor set
covering every item, or a list of sets named by test), and the
judgements compare persons across the tests. A person of a test after
the first responds at their location times that test's unit, less its
origin shift, on the test's own scale; the judgements identify both, and
the `tests` table reports them. A test after the first needs at least
one judged person with a score inside its range for its origin and two
for its unit; fix the unit at 1 by naming the test in `units` when there
is only one. Response-only locations are mapped to the reference scale
by the fitted unit and shift, and the unit carries the same upward bias
as the judgement units. Person identifiers must be unique across the
tests.

## References

Bradley, R. A. and Terry, M. E. (1952). Rank analysis of incomplete
block designs: I. The method of paired comparisons. Biometrika, 39,
324–345.

Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
Rasch model. Journal of Applied Measurement, 9(3), 249–264.

Plackett, R. L. (1975). The analysis of permutations. Applied
Statistics, 24(2), 193–202.

## See also

[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
[`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md),
[`frame_invariance`](https://drjoshmcgrane.github.io/rasch/reference/frame_invariance.md).

## Examples

``` r
set.seed(1)
delta <- seq(-1.5, 1.5, length.out = 8)
names(delta) <- sprintf("I%02d", 1:8)
theta <- rnorm(300)
X <- sapply(delta, function(d) as.integer(runif(300) < plogis(theta - d)))
pairs <- t(combn(names(delta), 2))[sample(28, 200, replace = TRUE), ]
p_a <- plogis(0.6 * (delta[pairs[, 1]] - delta[pairs[, 2]]))
cj <- data.frame(a = pairs[, 1], b = pairs[, 2],
                 winner = ifelse(runif(200) < p_a, pairs[, 1], pairs[, 2]))
fit <- rasch_cj(X, comparisons = cj, object_a = "a", object_b = "b",
                winner = "winner")
fit$units
#>         frame      unit        se estimated
#> 1   responses 1.0000000        NA     FALSE
#> 2 comparisons 0.5839283 0.1244835      TRUE
fit$invariance$lr
#>   statistic df         p
#> 1   9.22088  6 0.1615317

# persons from their responses and judgements of their work
pairs <- cbind(sample(300, 1200, replace = TRUE), sample(300, 1200, replace = TRUE))
pairs <- pairs[pairs[, 1] != pairs[, 2], ]
p_a <- plogis(0.8 * (theta[pairs[, 1]] - theta[pairs[, 2]]))
work <- data.frame(a = paste0("P", pairs[, 1]), b = paste0("P", pairs[, 2]))
work$winner <- ifelse(runif(nrow(work)) < p_a, work$a, work$b)
pfit <- rasch_cj(X, comparisons = work, object_a = "a", object_b = "b",
                 winner = "winner", objects = "persons", anchors = fit)
head(pfit$persons)
#>   person n_items raw max_raw    location        se location_responses
#> 1     P1       8   5       8 -0.18971656 0.5746338           0.619502
#> 2     P2       8   3       8 -1.03789001 0.7098991          -0.612274
#> 3     P3       8   3       8 -0.64292531 0.5860613          -0.612274
#> 4     P4       8   6       8  2.60722071 0.7949869           1.306549
#> 5     P5       8   5       8 -0.02951153 0.6177500           0.619502
#> 6     P6       8   1       8 -0.21241340 0.4897669          -2.253734
#>   se_responses location_comparisons comparisons rankings extreme
#> 1    0.7970050           -1.3648429        TRUE    FALSE   FALSE
#> 2    0.7998178                   NA       FALSE    FALSE   FALSE
#> 3    0.7998178           -0.2871821        TRUE    FALSE   FALSE
#> 4    0.8735421                   NA       FALSE    FALSE   FALSE
#> 5    0.7970050           -1.7776361        TRUE    FALSE   FALSE
#> 6    1.1138402            1.0997898        TRUE    FALSE   FALSE

# two tests with no item in common, linked by comparisons of their items
delta2 <- 1.5 * delta; names(delta2) <- sprintf("J%02d", 1:8)
theta2 <- 1.5 * rnorm(300)
X2 <- sapply(delta2, function(d) as.integer(runif(300) < plogis(theta2 - d)))
both <- c(delta, delta / 1)
names(both) <- c(names(delta), names(delta2))
pairs <- t(combn(names(both), 2))[sample(120, 400, replace = TRUE), ]
p_a <- plogis(0.6 * (both[pairs[, 1]] - both[pairs[, 2]]))
cj2 <- data.frame(a = pairs[, 1], b = pairs[, 2],
                  winner = ifelse(runif(400) < p_a, pairs[, 1], pairs[, 2]))
two <- rasch_cj(list(first = X, second = X2), comparisons = cj2,
                object_a = "a", object_b = "b", winner = "winner")
two$units
#>         frame      unit        se estimated
#> 1       first 1.0000000        NA     FALSE
#> 2      second 1.5033540 0.4573696      TRUE
#> 3 comparisons 0.5161414 0.1139780      TRUE
```
