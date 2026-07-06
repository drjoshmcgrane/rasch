# Category-score moments for a polytomous item

Category-score moments for a polytomous item

## Usage

``` r
item_moments(theta, tau_i, disc = 1)
```

## Arguments

- theta:

  Person location, in logits.

- tau_i:

  Numeric vector of the item's threshold parameters.

- disc:

  Discrimination (frame unit) multiplier on the exponent; 1 for the
  ordinary Rasch model.

## Value

A list with category probabilities `P`, expected score `E`, variance
`V`, and third and fourth central moments `mu3`, `mu4`.

## Examples

``` r
item_moments(0.5, c(-1, 0, 1))
#> $P
#> [1] 0.0576288 0.2582744 0.4258225 0.2582744
#> 
#> $E
#> [1] 1.884742
#> 
#> $V
#> [1] 0.7337796
#> 
#> $mu3
#> [1] -0.2057782
#> 
#> $mu4
#> [1] 1.285077
#> 
```
