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
bundles and differential test functioning in its DIF panel. The analysis
report can be written as Markdown, which the application offers in place
of its HTML export. No existing estimator is changed. The corrections to
the DIF tests, to the score-to-measure table and the information curves
of a split fit, and to the application and the report are described in
NEWS.

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
(aarch64-apple-darwin23): 0 errors, 0 warnings, 0 notes, in 312 seconds
including network checks. The five-second example threshold was set
explicitly; the slowest example was `plot_scree` at 1.9 s, and the
`--run-donttest` pass was also OK. The CRAN test selection passed 743
expectations, with 28 longer tests skipped on CRAN. All eight vignettes and
the PDF and HTML manuals passed.

The complete test suite ran with `NOT_CRAN=true` on the same machine with
no failures (11,976 expectations); its five skips are parallel-worker
integration tests that need an installed package namespace.

## win-builder

The submitted tarball was uploaded to win-builder on 1 October 2026 and
returned 0 errors, 0 warnings and 0 notes on both platforms, Windows
Server 2022: R-release (R 4.6.1) installation 63 seconds, check 612
seconds; R-devel (2026-09-30 r90605) installation 62 seconds, check 654
seconds. The slowest example on either platform was `plot_scree` at
8.4 seconds against the ten-second Windows threshold, with `resolve_dif`
next at 7.6 seconds; the CRAN test selection passed on both.

Earlier builds of this version were checked there as the work went in.
The build of 25 September, before the Markdown report and the
application and DIF corrections, returned 0 errors, 0 warnings and 0
notes on both (R-release check 572 seconds; R-devel check 654 seconds),
as did the build of 30 September, before the correction of the
information curves of a split fit (R-release check 593 seconds; R-devel
check 657 seconds). The build of 29 September passed R-devel (check 696
seconds) while the R-release checker was ending every upload of that
afternoon after seven seconds, at the CRAN incoming feasibility step,
including a control upload of a build that had passed that morning; the
service was at fault and the checks of 30 September and 1 October
completed normally. A build before those, differing from the 25
September build only in the `simulate_mfrm` example, drew one NOTE on
R-devel: that example ran in 10.19 seconds against the ten-second
Windows threshold, having taken 9.9 seconds in the same check of 1.13.0.
The example now simulates four raters instead of six; it runs in 1.0
second locally and 4.3 seconds on win-builder R-devel, and still
recovers the rater severities.
