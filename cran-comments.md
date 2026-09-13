# CRAN comments for rasch 1.12.1

## Summary

This maintenance update corrects estimation, diagnostic, plotting and saved
analysis defects described in NEWS. It also adds app simulation controls,
supplementary weighted person estimates and a data-structures vignette.

The author's name has been updated from Josh McGrane to Joshua A McGrane.
This is the same author and maintainer; the email address is unchanged.

## Check time

CRAN runs representative end-to-end workflows and focused regression tests.
The complete test suite remains enabled with `NOT_CRAN=true`, including in
continuous integration on Windows, macOS and Linux. Statistical validation
studies remain in the repository.

Three vignettes use recorded model or bootstrap calculations. The analysis
code, datasets, seeds and replication counts are retained; tables and
figures are rebuilt from those results. The regeneration script executes
the vignette code and records source and result hashes. Source builds and
CI verify these records. One longer bootstrap example is marked as optional.
The package's statistical computations and default replication counts have
not been reduced to shorten checks. CRAN's two-core limit is respected.

## Verification

`R CMD check --as-cran --timings` on macOS, R 4.6.1
(aarch64-apple-darwin23): 0 errors, 0 warnings, 0 notes. The CRAN test
selection passed 743 expectations, with 28 longer tests skipped on CRAN.
All eight vignettes and the PDF and HTML manuals passed.

On the same machine, building the source package fell from 180 to 60 seconds
and the full check from 473 to 315 seconds. Within the check, tests fell
from 77 to 47 seconds and vignette rebuilding from 108 to 27 seconds.
Overall check times include variable network checks.

The recorded-vignette tests and scree redraw also passed under R 4.5.1.
This tarball was uploaded to win-builder for R-release and R-devel on
13 September 2026. Both returned 0 errors, 0 warnings and 0 notes on
Windows Server 2022:

* R-release (R 4.6.1): installation 36 seconds; check 371 seconds.
* R-devel (2026-09-12 r90533): installation 35 seconds; check 406 seconds.

These times are from the result notifications; both logs have been verified.
