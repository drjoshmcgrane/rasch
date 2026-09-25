# CRAN comments for rasch 1.14.0

## Summary

This feature update follows the accepted 1.12.1. It adds a Plackett-Luce
model for rankings (`pl()`), a joint calibration of item responses, paired
comparisons and rankings (`rasch_cj()`, with its calibration map
`plot_cj()`), including persons measured from their responses and from
judgements of their work and tests of one construct with no item in common
linked through the judgements, differential bundle and test functioning
(`dif_anova(bundles = )`, `dtf()`, `plot_dtf()`), and a conditional Wald
test of differential item functioning (`dif_wald()`). The DIF analysis of
a split fit forms its class intervals on one score-to-measure mapping for
every group. The Shiny application gains the same estimators: a rankings
page, a joint calibration of judgements uploaded beside the responses that
anchors the response analysis, and the conditional Wald tests, item
bundles and differential test functioning in its DIF panel. No existing
estimator is changed. The application and report corrections are
described in NEWS.

The new calibrations were audited before submission with null simulation
screens and regression tests, held in the repository and not in the
package.

## Check time

CRAN runs representative end-to-end workflows and focused regression tests.
The complete test suite, including the tests of the new calibrations and
of the application pages, remains enabled with `NOT_CRAN=true` locally and
in continuous integration on Windows, macOS and Linux. Statistical
validation studies remain in the repository.

Three vignettes use recorded model or bootstrap calculations. The analysis
code, datasets, seeds and replication counts are retained; tables and
figures are rebuilt from those results. The regeneration script executes
the vignette code and records source and result hashes. Source builds and
CI verify these records. The records were regenerated for this version;
all three saved result files reproduced byte-for-byte. One longer bootstrap
example is marked as optional. The package's statistical computations and
default replication counts have not been reduced to shorten checks. CRAN's
two-core limit is respected.

The installed size is 8.2 MB, as the accepted 1.12.1 was 8.0 MB; the
`doc` directory holds the eight vignettes.

## Local checks

`R CMD check --as-cran --timings` on macOS Sequoia 15.6, R 4.6.1
(aarch64-apple-darwin23): 0 errors, 0 warnings, 0 notes, in 341 seconds
including network checks. The five-second example threshold was set
explicitly; the slowest example was `plot_scree` at 1.8 s, and the
`--run-donttest` pass was also OK. The CRAN test selection passed 743
expectations, with 28 longer tests skipped on CRAN. All eight vignettes and
the PDF and HTML manuals passed.

The complete test suite ran with `NOT_CRAN=true` on the same machine with
no failures; its five skips are parallel-worker integration tests that
need an installed package namespace.

## win-builder

The submitted tarball was uploaded to win-builder for R-release and
R-devel on 25 September 2026. Both returned 0 errors, 0 warnings and
0 notes on Windows Server 2022:

* R-release (R 4.6.1): installation 58 seconds; check 572 seconds.
* R-devel (2026-09-21 r90579): installation 61 seconds; check 654 seconds.

An earlier build of this version, differing only in the `simulate_mfrm`
example, passed R-release (R 4.6.1) with 0 errors, 0 warnings and 0 notes
(installation 61 seconds; check 578 seconds) and drew one NOTE on R-devel:
that example ran in 10.19 seconds against the ten-second Windows
threshold, having taken 9.9 seconds in the same check of 1.13.0. The
example now simulates four raters instead of six; it runs in 1.0 second
locally and 4.0 to 4.5 seconds on win-builder, and still recovers the
rater severities. The `resolve_dif` example, which drew the same NOTE in
an earlier upload of 1.13.0, ran in 6.7 to 7.4 seconds on win-builder.
