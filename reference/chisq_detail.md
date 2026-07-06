# Class-interval detail for one item's chi-square test of fit

The per-class-interval breakdown behind an item's item-trait chi-square,
as dissected in Andrich and Marais (2019, ch. 13): for every class
interval the size, the maximum and mean person location, the
standardised residual between observed and expected interval means, its
squared chi-square component, the observed and expected means (OM, EV),
the sample-size-free effect size ES = (OM - EV)/sqrt(mean V), and per
response category the observed proportion (OBS.P), the mean model
probability (EST.P), and the observed conditional threshold proportion
(OBS.T), the proportion scoring k among those scoring k - 1 or k.

## Usage

``` r
chisq_detail(fit, item)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- item:

  Item name or index.

## Value

A list with `item`, `location`, the `intervals` data frame, the
`categories` data frame, the whole-sample observed mean `ave`, and the
item's total `chisq`, `df`, and `p`. Intervals with fewer than 2
responders are shown but carry no chi-square contribution
(`used = FALSE`), matching the item-trait computation.

## Examples

``` r
set.seed(1)
d <- seq(-1.5, 1.5, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
chisq_detail(rasch(X), "I3")$intervals
#>   interval   n   theta_max  theta_mean   obs_mean exp_value   residual
#> 1        1  52 -1.61901055 -1.61901055 0.07692308 0.1962929 -2.1671803
#> 2        2  67 -0.73551099 -0.73551099 0.38805970 0.3714206  0.2818740
#> 3        3 118  0.01710555  0.01710555 0.54237288 0.5563815 -0.3062983
#> 4        4  78  0.75817394  0.75817394 0.78205128 0.7246324  1.1352375
#> 5        5  43  1.61003743  1.61003743 0.93023256 0.8604966  1.3198477
#>        chisq          es used
#> 1 4.69667044 -0.30053383 TRUE
#> 2 0.07945293  0.03443639 TRUE
#> 3 0.09381862 -0.02819704 TRUE
#> 4 1.28876421  0.12854034 TRUE
#> 5 1.74199785  0.20127488 TRUE
```
