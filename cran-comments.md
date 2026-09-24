# CRAN comments for rasch 1.13.0

## Summary

This feature update adds a Plackett-Luce model for rankings (`pl()`), a
joint calibration of item responses, paired comparisons and rankings
(`rasch_cj()`), including persons measured from their responses and from
judgements of their work and tests of one construct with no item in common
linked through the judgements, differential bundle and test functioning
(`dif_anova(bundles = )`, `dtf()`, `plot_dtf()`), and a conditional Wald
test of differential item functioning (`dif_wald()`). The DIF analysis of
a split fit forms its class intervals on one score-to-measure mapping for
every group. No existing estimator is changed. The application and report
corrections are described in NEWS.

The new calibrations were audited before submission with null simulation
screens and regression tests, held in the repository and not in the
package.

## Check time

CRAN runs representative end-to-end workflows and focused regression tests.
The complete test suite, including the tests of the new calibrations,
remains enabled with `NOT_CRAN=true` locally and in continuous integration
on Windows, macOS and Linux. Statistical validation studies remain in the
repository.

Three vignettes use recorded model or bootstrap calculations. The analysis
code, datasets, seeds and replication counts are retained; tables and
figures are rebuilt from those results. The regeneration script executes
the vignette code and records source and result hashes. Source builds and
CI verify these records. The records were regenerated for this version;
all three saved result files reproduced byte-for-byte. One longer bootstrap
example is marked as optional. The package's statistical computations and
default replication counts have not been reduced to shorten checks. CRAN's
two-core limit is respected.

The installed size is 8.1 MB, as the accepted 1.12.1 was 8.0 MB; the
`doc` directory holds the eight vignettes.

## Local checks

`R CMD check --as-cran --timings` on macOS Sequoia 15.6, R 4.6.1
(aarch64-apple-darwin23): 0 errors, 0 warnings, 0 notes, in 325 seconds
including network checks. The five-second example threshold was set
explicitly; the slowest example was `resolve_dif` at 3.477 s, and the
`--run-donttest` pass was also OK. The CRAN test selection passed 743
expectations, with 28 longer tests skipped on CRAN. All eight vignettes and
the PDF and HTML manuals passed.

The complete test suite ran with `NOT_CRAN=true` on the same machine with
no failures; its five skips are parallel-worker integration tests that
need an installed package namespace.
