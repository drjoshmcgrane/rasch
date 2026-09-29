# App follow-up audit — 26 September 2026

Baseline: `59ab6ba`, following the app additions in `b097983`.
The earlier estimator audit remains in `audit-september-additions.md`.
This follow-up concerns the new app displays, generated code and joint map.

## Corrections

| Issue | Correction |
| --- | --- |
| Withheld invariance tests appeared as zero moving objects, with a green indicator. | Ranking and joint-calibration displays distinguish unavailable from nonsignificant results. Text and joint-summary CSVs report tested and withheld contrasts. Printed summaries use the same rule. |
| Separate calibration markers used unrelated origins. | `plot_cj()` uses new `frame_locations` metadata, aligning each connected block to its combined origin without changing within-block differences. Response frames retain markers for items without judgements. Raw separate estimates remain in `items`, with their origins explained in the help. Old fits without the metadata warn and show combined locations only. |
| A significant reversal test was always coloured red. | The indicator is adverse only when worst-first fits significantly better. Best-first evidence is not marked as a reversal problem. |
| Repeated bundle names silently replaced earlier definitions. | The parser rejects duplicate names, including duplicates differing only in surrounding whitespace. |
| Backslashes in bundle names broke generated R code. | Non-syntactic names use quoted, escaped string literals. Tests round-trip backslashes, quotes, backticks, newlines and reserved words. |
| DTF help misdescribed unsigned functioning and confidence bands. | The help distinguishes cancellation across person locations from cancellation between items, and describes zero exclusion by a pointwise band. |

The follow-up sweep found two related paths: failed separate calibrations
could still supply plot markers, and a failed frame could vanish from a
three-frame contrast table while other tests remained. Failed-frame markers
are now omitted. `invariance$n_contrasts` records the planned contrast count,
so omitted tests count as withheld rather than disappearing from the
denominator. The plot help also distinguishes a significant threshold
contrast from a difference in the item mean.

## Statistical scope

These changes correct reporting and plotting, not estimation. For plotting,
each separate frame location receives a block-specific additive constant:
the mean combined location minus the mean separate location over that
block's reached objects. A response frame uses its informative item means;
a judgement frame uses its judged objects, including thresholds where
applicable. Separate origins cannot be compared as absolute differences.

`tools/audit-september-app.R` compares the current estimator with `59ab6ba`
in five reproducible designs: a judged subset, disconnected judgement
blocks, polytomous item judgements, polytomous threshold judgements, and
multiple response tests. Existing estimates, covariance matrices,
likelihoods, convergence results and invariance tables are exactly
identical in all five. The new reporting metadata is excluded from that
identity comparison. This is a regression check, not a new simulation
claim about Type I error or coverage.

The map was also rendered and inspected for the multiple-test design.

## Verification

The integration block passed 1,014 expectations with zero failures, errors
or warnings. It included the three app integration files, the new audit
tests, `plot-cj`, `shiny-help`, `rasch-cj`, `pl`, `dtf` and the earlier
`audit-new-calibrations` tests. The final accounting follow-up
(`latest-app-audit`, `plot-cj`, `shiny-help`, `rasch-cj`) passed 472
expectations, also with zero failures, errors or warnings. This includes
failed separate frames and a failed combined fit. Runs overlap and their
counts should not be added.

Reproduce the combined targeted block from the package root:

```r
testthat::test_local(filter = paste0(
  "^(latest-app-audit|app-joint|app-rankings|app-dif-extensions|plot-cj|",
  "shiny-help|rasch-cj|pl|dtf|audit-new-calibrations)$"))
```

```sh
Rscript tools/audit-september-app.R
```

Help pages were regenerated from source. R source and app files parse, and
`git diff --check` is clean. The full suite, `R CMD check`, release build,
commit and push are not part of this work. This is a targeted follow-up,
not a claim that the entire package is free of defects.
