# September additions: audit and validation

Reviewed on 24 September 2026, against `de67772`. Scope: the nine commits
after `1b4b492`, covering rankings, joint response–judgement calibration,
multiple tests and person calibration, DIF Wald tests, bundle and test
functioning, common intervals after splits, and app/report integration.
The corrections described here are working-tree changes, not a release.

## Corrections

| Area | Reproduced problem | Correction and evidence |
| --- | --- | --- |
| Joint calibration convergence | A zero score at a stationary minimum or along an unidentified direction was accepted as convergence. Separate failed fits could supply covariance and likelihood-ratio inference. | Require finite derivatives and negative-definite likelihood curvature. Withhold failed-fit inference. Tests distinguish a minimum, flat direction and identified maximum. |
| Joint calibration contrasts | Contrasts divided separate estimates by a unit estimated from the same data but omitted its uncertainty and shared-data covariance. | Combine direct and linking influences before forming covariance, in item and person modes. An independent matrix calculation checks the cross terms; null screens cover fixed and estimated units. |
| Disconnected response booklets | Separate response calibration imposed one origin on disconnected item components, although each component needs its own origin. | Component-specific constraints and parameter counts. A two-booklet example verifies LR degrees of freedom of 2 rather than 3. |
| Scores, thresholds and anchors | Fractional person scores or threshold numbers could be truncated; factors could be read as integer level codes; failed fits could supply anchors. | Parse numeric labels explicitly, reject malformed scores and thresholds, and require eligible converged anchor fits. Preserve literal item names containing a colon and refuse ambiguous threshold keys. |
| Shared response items | The same named item could be recoded from different raw category mappings in different tests and then treated as identical. | Require matching category counts and raw-to-model score mappings. |
| Ranking covariance | Follow-up contrasts could replace withheld sandwich errors with model errors, and omitted covariance between early and later estimates from the same ranking/judge. | Retain the requested covariance, enforce subgroup support, and use differences of cluster influences. Judge covariance uses the `G/(G−1)` correction and contrasts use finite-cluster t degrees of freedom. An independent score/Hessian calculation checks the result. |
| Ranking follow-ups | A chi-square LR reference remained available for dependent judges. Reversal comparisons used the wrong anchor orientation and could compare different retained data. | Withhold judge-clustered LR probabilities; fit both reversal models to the same complete rankings, reverse fixed anchors, and cluster likelihood differences. Apply support checks to the retained data. |
| Ranking numerics | Large shifts of fixed anchors overflowed probabilities; extreme fixed anchors could be discarded. Required role columns could be NULL. | Stable log-softmax, origin-aware starts, retained fixed anchors and explicit column validation. A 1,000-logit origin shift preserves estimates, uncertainty and likelihood after translation. |
| Differential test functioning | Indefinite covariance or weak categories could produce inference; repeated-person degrees of freedom and unavailable members of adjustment families were lost. Redundant bundle contrasts inflated homogeneity degrees of freedom. | Apply covariance/category/support checks, complete-family Holm adjustment and supported t/F references; use the rank of bundle contrast covariance. A bundle with two invariant members and one split member has one independent contrast. |
| DIF Wald edge case | No factor levels meeting `min_n` caused a malformed data frame instead of an unavailable result. | Retain the unavailable item-factor rows with zero contrast degrees of freedom. |
| Saved split-item DIF | Saved results could retain the earlier group-specific interval mapping despite the new common-interval algorithm. | Stamp primary DIF, bootstrap and resolution results. Authenticate saved projects, omit affected obsolete results, and retain source data and fitted models. |
| DTF plots | Custom titles/limits duplicated arguments; unavailable bands could exclude the actual curve from the range. | Merge plot overrides and include finite curve values in automatic limits. |

The regression file is `tests/testthat/test-audit-new-calibrations.R`.
Existing tests were updated only for intentional result metadata, clearer
validation messages, and an anchor fixture that previously used a failed
free-unit fit. The replacement fixture fixes the unit and asserts convergence.

## Covariance calculation

For independent source `f`, let `I_f` be its natural-parameter information,
`J_f` the Jacobian into the joint parameterisation, `V` the joint covariance,
`D_f` the direct influence of a separate fit, and `G` the contrast derivative
through the estimated link. The implemented contrast covariance is

```
A_f = D_f + G V J_f'
Var(contrast) = sum_f A_f I_f A_f'
```

Taking the cross-product after adding the two influences retains their
covariance. Adding their variances separately does not. For ranking contrasts,
the corresponding empirical calculation uses the difference of the two
within-cluster influences before its cross-product. This follows the
cluster-summed estimating-function construction documented by
[sandwich](https://sandwich.r-forge.r-project.org/reference/vcovCL.html).

## Checks

The original new-feature tests passed before the audit (509 expectations).
Independent edge cases nevertheless reproduced failures, so that initial
pass was not treated as evidence that these paths were sound.

Final targeted runs were:

| Test block | Passed expectations | Failures / errors / warnings |
| --- | ---: | --- |
| New models, DIF/DTF and independent regression cases | 571 | 0 / 0 / 0 |
| Saved projects, structural refits and independent regression cases | 420 | 0 / 0 / 0 |
| DIF integration, project migrations and app regression files | 479 | 0 / 0 / 0 |

The blocks overlap; these counts must not be added as distinct tests.
An earlier integration run also exercised DIF bootstrap, subtests and report
formatting. Its only failing expectation was the old exact field list for
`resolve_dif`; that list was updated for the new interval stamp and the
structural-refit block rerun successfully.

All 43 R source files parsed. A global-symbol scan of the five new model
files found no unresolved symbols, and `git diff --check` was clean. The
full package suite and `R CMD check` are not part of this run.

The five help-page examples for `pl`, `rasch_cj`, `dif_wald`, `dtf` and
`plot_dtf` executed successfully. Local elapsed times were 0.13, 1.92, 1.95,
0.45 and 0.51 seconds respectively. These are local measurements, not a
guarantee about CRAN hardware.

## Null screens

Eight scenarios used 400 attempts each. Probabilities below are conditional
on the prespecified contrast being available; unavailable attempts are
retained in the accounting. MCSE is in percentage points and SD/SE means
empirical SD divided by mean reported SE.

| Scenario | Analysed / attempted | Unavailable | Rejection at 5% | MCSE (pp) | SD/SE |
| --- | ---: | ---: | ---: | ---: | ---: |
| PL, independent rankings, model covariance | 400 / 400 | 0 | 5.25% | 1.12 | 1.005 |
| PL, independent rankings assigned to 20 judges | 400 / 400 | 0 | 4.75% | 1.06 | 1.033 |
| PL, 20 judges, six identical rankings per judge | 302 / 400 | 98 | 0.66% | 0.47 | 0.682 |
| PL, 50 judges, three identical rankings per judge | 393 / 400 | 7 | 1.02% | 0.51 | 0.926 |
| Joint item calibration, fixed judgement unit | 400 / 400 | 0 | 4.25% | 1.01 | 0.946 |
| Joint item calibration, estimated judgement unit | 400 / 400 | 0 | 5.00% | 1.09 | 0.981 |
| DIF Wald, prespecified item | 400 / 400 | 0 | 4.00% | 0.98 | 0.963 |
| DTF, prespecified mean shift | 400 / 400 | 0 | 4.00% | 0.98 | 0.984 |

Thus 3,095 of 3,200 attempts yielded the requested contrast. The independent
designs are compatible with the nominal rate at this simulation precision.
The duplicated-ranking designs are markedly conservative, not well calibrated
to 5%. They contain only 20 or 50 independent rankings; sparse first-choice
calibrations can lose objects or give weak Wald approximations. This is
evidence of limited information and power, not a reason to tune the correction
to these cells. The manual now states that position contrasts may be
conservative or unavailable despite usable main-model standard errors.

The table is from `output/september-audit-final-screen/summary.rds` only;
earlier development screens are superseded. Every final replicate file has
400 entries of three values, including explicit NA entries for unavailable
contrasts. The completed run exited successfully, and the saved launch hash
of the study script matches its final file. R-source hashes describe the
launch state; subsequent help-text edits changed source comments, not the
estimators used by these screens.

## Reproduction and scope

Run from the package root:

```r
testthat::test_local(filter = "^(audit-new-calibrations|rasch-cj|pl|dtf|dif-wald|dif-split-intervals)$")
```

```sh
Rscript tools/audit-september-additions.R 400 output/september-audit
```

The simulation script records attempted, analysed and unavailable counts,
binomial intervals, Monte Carlo errors, bias and empirical SD / mean SE.
Each saved summary includes source hashes captured at launch and the R
session information. Hashes identify the executed files, not a later HEAD.
The RDS files remain local under `output/`; the script and this record are
excluded from the CRAN tarball.

The screens use one prespecified contrast per replication and unadjusted
5% probabilities. They assess marginal size, not familywise error or power.
They do not establish calibration for every item count, missingness pattern,
polytomous design or linking arrangement.

Remaining assumptions are explicit in the help pages:

- `rasch_cj` assumes independent sources and independent observations within
  each source. It does not provide judge-clustered joint inference.
- Person-mode item anchors are treated as fixed. Its free-unit incidental-
  parameter bias is not removed by correcting contrast covariance; the new
  null screens do not establish person-mode interval coverage.
- The PL reversal reference assumes distinguishable models. There is no
  distinguishability test at the overlapping equal-worth model. Its regular
  Vuong probability is not validated for that case.
- DTF uncertainty is conditional on the selected splits, anchors and
  averaging locations. It does not include model selection or population
  sampling uncertainty in those locations.

These are limits on interpretation, not claims of universal validation.
