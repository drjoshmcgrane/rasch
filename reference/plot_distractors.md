# Plot multiple-choice option curves

The proportion choosing each response option across class intervals of
the rest measure (the person estimate from the other items), with the
keyed option drawn solid and bold. The keyed option should rise with the
trait and every distractor should fall; a rising distractor is the
graphical signature of a miskey or an ambiguous option.

## Usage

``` r
plot_distractors(fit, item, n_groups = fit$n_groups)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md) run
  with a `key`.

- item:

  Keyed item name.

- n_groups:

  Number of class intervals.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1); Np <- 400
th <- rnorm(Np)
raw <- sapply(seq(-1, 1, length.out = 6), function(d) {
  ok <- rbinom(Np, 1, plogis(th - d))
  ifelse(ok == 1, "A", sample(c("B", "C", "D"), Np, replace = TRUE))
})
colnames(raw) <- paste0("M", 1:6)
fit <- rasch(raw, key = setNames(rep("A", 6), colnames(raw)))
plot_distractors(fit, "M3")
```
