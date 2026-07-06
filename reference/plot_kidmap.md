# Plot a kidmap

The person diagnostic map (Wright, Mead and Ludlow 1980): item
thresholds the person achieved to the right of a vertical logit axis and
thresholds not achieved to the left, with the person's location drawn as
a dashed line inside its confidence band. Achieved thresholds above the
band (unexpected successes) and unachieved thresholds below it
(unexpected failures) are highlighted; a clean response pattern shows
achieved thresholds below the band and unachieved ones above it.

## Usage

``` r
plot_kidmap(
  fit,
  person,
  level = 0.95,
  bins = 35,
  xlim = NULL,
  cex_labels = 0.8
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- person:

  Row number of the person, or an ID matching `fit$person$id`.

- level:

  Confidence level of the band around the person location used to mark
  unexpected responses.

- bins:

  Number of vertical bins used to stack the threshold labels.

- xlim:

  Optional logit range; thresholds outside it are omitted.

- cex_labels:

  Character expansion for the threshold labels.

## Value

Called for its plotting side effect; invisibly `NULL`.

## References

Wright, B. D., Mead, R. J., & Ludlow, L. H. (1980). *KIDMAP:
person-by-item interaction mapping* (Research Memorandum No. 29).
Chicago: University of Chicago, MESA Psychometric Laboratory.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 12)
X <- matrix(rbinom(300 * 12, 1, plogis(outer(rnorm(300), d, "-"))), 300, 12)
colnames(X) <- paste0("I", 1:12)
plot_kidmap(rasch(X), person = 1)
```
