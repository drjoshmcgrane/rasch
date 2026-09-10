# Diagnose fixed departures from an explanatory model

Fits each available item-location, polytomous threshold-structure or
comparative-judgement object departure separately from the active model.
Probabilities use Kent calibration and Holm adjustment over the complete
candidate family. The Kent departure tests are first-order asymptotic
comparisons and do not use the finite-person-cluster correction applied
to individual coefficients. A candidate with a withheld probability,
including a refit that errors or fails to converge, remains in that
family.

## Usage

``` r
explanatory_diagnostics(fit, p_adjust = "holm")
```

## Arguments

- fit:

  A fitted explanatory Rasch or comparative judgement model.

- p_adjust:

  Multiplicity adjustment over the candidate departures.

## Value

A data frame ordered by adjusted probability. For item fits, a `weak`
column marks items whose thresholds the calibration flags as weakly
identified; their probabilities are withheld, since the departure test
rests on the same sparse categories, and a note on the table records the
withholding. The `converged` column identifies candidate refits that
converged. Statistics from a failed or non-convergent candidate are
withheld, but it remains in the multiplicity family.
