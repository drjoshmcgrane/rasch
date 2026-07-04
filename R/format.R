# rmt :: shared display formatting
# ===========================================================================
# One formatting vocabulary for every surface that prints results (console
# methods, the HTML report, and the app): probabilities below the display
# resolution read "< 0.001" rather than "0.000"; numeric columns carry
# fixed decimals; integer-valued columns print as integers; logical flags
# print as "*" or blank.
# ===========================================================================

.fmt_p <- function(p, digits = 3) {
  lim <- 10^-digits
  out <- ifelse(is.finite(p) & p < lim, paste0("< ", format(lim, scientific = FALSE)),
                ifelse(is.finite(p), sprintf(paste0("%.", digits, "f"), p), ""))
  out
}

.is_pcol <- function(nm)
  grepl("^p$|^p_|_p$|^prob$|^p\\.", nm)

# Format a data frame for display: fixed decimals, "< 0.001" probabilities,
# clean integers, starred logicals. Returns a character data frame.
.fmt_df <- function(d, digits = 3) {
  out <- d
  for (j in seq_along(d)) {
    v <- d[[j]]; nm <- names(d)[j]
    if (is.logical(v)) out[[j]] <- ifelse(is.na(v), "", ifelse(v, "*", ""))
    else if (is.numeric(v)) {
      if (.is_pcol(nm)) out[[j]] <- .fmt_p(v, digits)
      else if (all(is.na(v) | v == round(v))) out[[j]] <-
          ifelse(is.na(v), "", format(v, scientific = FALSE, trim = TRUE))
      else out[[j]] <- ifelse(is.na(v), "",
                              sprintf(paste0("%.", digits, "f"), v))
    }
  }
  out
}
