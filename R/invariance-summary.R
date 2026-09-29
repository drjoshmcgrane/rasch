# Shared reporting for per-object invariance tests. An unavailable test is
# not a nonsignificant result; counts refer only to finite adjusted p values.
.invariance_summary <- function(tab, labels, n_expected = NULL) {
  tested <- is.finite(tab$p_adj)
  n_tested <- sum(tested)
  total <- if (is.null(n_expected)) length(tested) else max(n_expected, length(tested))
  n_withheld <- total - n_tested
  moved <- unique(labels[tested & tab$p_adj < 0.05])
  coverage <- if (!n_tested) {
    if (n_withheld) sprintf("no estimable contrasts; %d withheld", n_withheld)
    else "no estimable contrasts"
  } else
    sprintf("%d contrasts tested%s", n_tested,
            if (n_withheld) sprintf("; %d withheld", n_withheld) else "")
  list(tested = n_tested, withheld = n_withheld, moved = moved,
       n_moved = if (n_tested) length(moved) else NA_integer_,
       coverage = coverage,
       text = if (!n_tested) sprintf("unavailable (%s)", coverage) else
         sprintf("%s (%s)", if (length(moved)) paste(moved, collapse = ", ")
                 else "none", coverage),
       status = if (length(moved)) "bad" else
         if (!n_tested || n_withheld) "neutral" else "good")
}
