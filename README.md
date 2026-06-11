# RaschR

A pairwise Rasch measurement engine in R, built from published measurement
theory. Estimation is pairwise conditional; person locations use Warm's
weighted likelihood; a full diagnostic suite covers fit, dimensionality,
local dependence, and DIF.

## Install

```r
# install.packages("remotes")
remotes::install_github("drjoshmcgrane/RaschR")
```

Dependencies are base R only (`stats`, `graphics`, `grDevices`). The Shiny GUI
additionally needs the `shiny` package.

## Quick start

```r
library(RaschR)

fit <- rasch(responses, model = "PCM", solver = "LS", n_groups = 10)

fit$items            # locations, infit/outfit, fit residual, item-trait chi-square
fit$psi              # person separation index and separation ratio
fit$thresholds_diag  # ordering, reversals, never-modal categories per item

plot_icc(fit, "Q05")        # ICC with observed class-interval means
plot_ccc(fit, "Q05")        # category probability curves
plot_pimap(fit)             # person-item threshold map
plot_resid_cor(fit)         # residual-correlation heatmap (blues to reds)

dimensionality_test(fit)              # residual-PCA contrast + person t-test
residual_correlations(fit)            # local dependence
dif_anova(fit, group = sex_vector)    # uniform and non-uniform DIF
```

`responses` is a persons-by-items matrix or data frame of integer scores
starting at 0. Missing data is allowed (pairwise deletion in estimation).

## GUI

```r
shiny::runApp(system.file("shiny", "app.R", package = "RaschR"))
```

The app opens on built-in demo data and accepts a CSV upload (items as columns,
an optional grouping column for DIF).

## What it implements

- Pairwise conditional estimation. Dichotomous data has a clean principal
  eigenvector solution (`est_eigen_dich`) and a Choppin/least-squares solution;
  polytomous data uses inverse-variance weighted pairwise comparisons solved by
  least squares (`solve_LS`) or reciprocal averaging (`solve_RA`), which
  coincide once the comparison design is structurally incomplete.
- Partial credit (PCM) and rating scale (RSM, fitted under constraint) models.
- Warm's weighted likelihood person estimates, finite at extreme scores.
- Standardised residuals, infit and outfit mean squares, standardised fit
  residual; item-trait interaction chi-square over class intervals with
  Bonferroni adjustment and an optional sample-size adjustment.
- Person Separation Index and separation ratio.
- Threshold and category diagnostics (ordering, reversals, never-modal cats).
- Residual-PCA dimensionality test (Smith person t-test), residual-correlation
  local dependence, DIF by two-way residual ANOVA.

## Documentation

Function documentation is written as `roxygen2` comments in the source. To
generate the manual pages and refresh `NAMESPACE`:

```r
devtools::document()   # or roxygen2::roxygenise()
```

After that, `R CMD check --as-cran` runs clean. Every documented function
carries a runnable example.

## Methodological references

Choppin (1985) pairwise estimation; Andrich & Luo (1993) pairwise algorithm;
Warm (1989) weighted likelihood; Smith (2002) residual-component dimensionality;
Andrich item-trait interaction chi-square.

## Status and known limits

Every estimation and diagnostic component is validated against simulated data
with known parameters in `tests/testthat`. Two items remain before this is
calibrated rather than merely correct:

1. The standardised fit residual is from the correct family (a normalising
   transform of the residual mean square) but its exact scaling constant is a
   placeholder until checked against output from established Rasch software.
2. The item-trait chi-square shows the genuine sample-size sensitivity of the
   statistic (the total/df ratio rises with N even on model-true data); this is
   a property of the statistic, managed via class-interval count, Bonferroni,
   and the sample-size adjustment, not a defect.

A single exported analysis from established Rasch software would let both the
fit residual constant and the item/threshold values be confirmed against ground
truth.
