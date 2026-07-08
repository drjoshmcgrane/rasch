# Likelihood-ratio test of the partial credit against the rating scale model

A likelihood-ratio test in the tradition of Andersen (1973): an
unrestricted (partial credit) analysis is compared with the rating
re-parameterisation of the same model on the same data. Twice the
difference in the pairwise conditional log-likelihoods is referred to a
chi-square on the difference in the number of threshold parameters. A
non-significant outcome supports adopting the simpler rating
parameterisation.

## Usage

``` r
lr_test(fit, maxit = 60, tol = 1e-08)
```

## Arguments

- fit:

  A `"PCM"` fit from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  with equal maximum scores across items (the rating parameterisation
  requires them).

- maxit, tol:

  Passed to the rating-scale refit.

## Value

A list of class `"rasch_lr"`: raw `chisq`, `df`, `p` (the conventional
display); adjusted `chisq_adj`, `p_adj`, and the eigenvalues `lambda`;
the two log-likelihoods; and the rating-scale refit (`fit_rsm`).

## Details

The likelihood here is the pairwise composite likelihood, not a full
likelihood, and twice its difference is not chi-square distributed: each
response enters every pair its item forms, so the raw statistic is
inflated. Two statistics are therefore reported. `chisq` is the raw
composite value with its naive `p`, the conventional display. The
limiting law of the raw statistic is \\\sum_j \lambda_j \chi^2_1\\ (Kent
1982; Varin, Reid and Firth 2011) with \\\lambda_j\\ the eigenvalues of
\\(C'H^{-1}C)^{-1}\\C'H^{-1}JH^{-1}C\\ over the \\r\\ constrained
directions \\C\\ (the part of the partial-credit threshold space outside
the rating subspace), estimated from the same Godambe \\H\\ and \\J\\
matrices that supply the sandwich standard errors; matching the mean
gives `chisq_adj` \\= r W / \sum_j \lambda_j\\ on \\r\\ degrees of
freedom. Use `p_adj` for inference; the naive `p` is severely
anticonservative and kept only for comparability with conventional
software displays.

## References

Kent, J. T. (1982). Robust properties of likelihood ratio tests.
Biometrika, 69, 19-27. Varin, C., Reid, N. and Firth, D. (2011). An
overview of composite likelihood methods. Statistica Sinica, 21, 5-42.

## Examples

``` r
set.seed(1)
tau <- c(-0.7, 0.7)
X <- sapply(seq(-1, 1, length.out = 6), function(d) vapply(rnorm(300),
  function(b) sample(0:2, 1, prob = item_moments(b, tau + d)$P), 0L))
colnames(X) <- paste0("Q", 1:6)
lr_test(rasch(X, model = "PCM"))
#> Likelihood-ratio test: partial credit vs rating parameterisation
#>   Raw composite chi-square 8.741 on 5 df, p = 0.120 (conventional display; anticonservative)
#>   Adjusted chi-square 2.062 on 5 df, p = 0.840 (Kent 1982 first-order calibration)
#>   log-likelihood (pairwise composite): PCM -2814.271, RSM -2818.642
```
