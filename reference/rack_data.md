# Reshape repeated measurements for racked or stacked analysis

Repeated measurements (the same persons and items at two or more time
points) enter a Rasch analysis in one of two designs (Andrich & Marais
2019, ch. 26). *Racking* keeps one row per person and duplicates the
items per time point (columns `item@time`), so change over time shows in
the item estimates. *Stacking* keeps one column per item and duplicates
the persons per time point (rows `person@time`), so change shows in the
person estimates and DIF of items over time can be examined with `time`
as a person factor.

## Usage

``` r
rack_data(data, person, time, items)

stack_data(data, person, time, items)
```

## Arguments

- data:

  A long data frame with one measurement per row.

- person, time:

  Names of the person and time-point columns.

- items:

  Character vector naming the item columns.

## Value

`rack_data`: a wide data frame with one row per person and
`length(items) * n_times` item columns. `stack_data`: a data frame with
one row per person-time (`id` column), the original item columns, and
`time` as a factor column for DIF analysis.

## Examples

``` r
d <- data.frame(pid = rep(1:100, 2), t = rep(1:2, each = 100),
                Q1 = rbinom(200, 1, 0.6), Q2 = rbinom(200, 1, 0.5))
racked <- rack_data(d, person = "pid", time = "t", items = c("Q1", "Q2"))
names(racked)
#> [1] "id"   "Q1@1" "Q2@1" "Q1@2" "Q2@2"
stacked <- stack_data(d, person = "pid", time = "t", items = c("Q1", "Q2"))
head(stacked)
#>    id time Q1 Q2
#> 1 1@1    1  0  1
#> 2 2@1    1  0  0
#> 3 3@1    1  1  0
#> 4 4@1    1  0  1
#> 5 5@1    1  1  0
#> 6 6@1    1  1  1
```
