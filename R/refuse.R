# rasch :: telling a refusal apart from a fault
# ===========================================================================
# Some analyses cannot be run on some designs, and saying so is a result
# rather than a failure. Residual principal components are undefined when
# item columns share no respondents, because a correlation between them does
# not exist to be decomposed; a factor that defines a frame cannot also be
# tested for differential functioning within that frame, because it takes one
# level among the persons who respond there. Both are conclusions the package
# draws deliberately, and both are worth reading.
#
# A caller cannot tell either apart from a genuine fault when both arrive as
# a bare error. The application showed these considered refusals in the same
# red as a crash, which reads as the package having broken. Signalling them
# with their own condition class lets a caller answer each in the right
# voice, and leaves anything unclassed to be treated as the fault it is.
# ===========================================================================

# Signal a deliberate refusal. Inherits from "error", so an unprepared caller
# still stops and every existing expectation on the message still holds.
.refuse <- function(...) {
  stop(structure(
    class = c("rasch_refusal", "error", "condition"),
    list(message = paste0(...), call = NULL)))
}

# NULL is a recorded request for automatic intervals, not missing metadata.
# Older fits record only the realised count, which may be below the public
# minimum because tied locations are never split across intervals.
.refit_n_groups <- function(fit) {
  if ("n_groups" %in% names(fit$refit_spec)) return(fit$refit_spec$n_groups)
  ng <- fit$n_groups
  if (length(ng) == 1L && is.finite(ng) && ng >= 2L) ng else NULL
}

# Restrictions in the fitted estimator remain authoritative when older
# objects lack the optional settings used to replay a calibration.
.has_calibration_anchors <- function(fit) {
  NROW(fit$refit_spec$anchors) > 0L || NROW(fit$est$anchors) > 0L ||
    any(fit$est$thr$anchored %in% TRUE)
}

.has_pc_thresholds <- function(fit) {
  !is.null(fit$refit_spec$pc_components) ||
    !is.null(fit$est$n_components)
}

# Fully fixed scoring fits have no item-estimation stage to repeat. In
# particular, older tailored results have n_parameters = 0 but no refit_spec;
# absence of that metadata must not be interpreted as an unanchored model.
.require_refittable_calibration <- function(fit) {
  if (inherits(fit, "rasch") &&
      (isTRUE(fit$refit_spec$fixed_calibration) ||
       isTRUE(fit$est$n_parameters == 0L)))
    .refuse("this fully anchored calibration is a scoring fit; downstream ",
            "recalibration and bootstrap procedures are not supported. ",
            "Use the original calibration for a new analysis")
  if (inherits(fit, "rasch")) {
    if (.has_calibration_anchors(fit) &&
        NROW(fit$refit_spec$anchors) == 0L)
      .refuse("the fitted calibration has anchors but its saved anchor ",
              "settings are unavailable; refit from the source data with ",
              "the original anchors before recalibration or bootstrapping")
    if (.has_pc_thresholds(fit) &&
        is.null(fit$refit_spec$pc_components))
      .refuse("the fitted calibration has principal-component threshold ",
              "constraints but its saved component settings are unavailable; ",
              "refit from the source data with the original pc_components ",
              "before recalibration or bootstrapping")
  }
  invisible(TRUE)
}

# Dichotomous ETS classification on the log-odds scale. `p` and `p_beyond`
# are the probabilities used for decisions; callers that report a family of
# comparisons pass their familywise-adjusted values here.
ETS_DELTA_PER_LOGIT <- 2.35

.ets_p_beyond <- function(difference, se, df = Inf) {
  a_cut <- 1.0 / ETS_DELTA_PER_LOGIT      # 0.43 as published
  d <- abs(difference)
  # C additionally requires the magnitude to be significantly beyond the A
  # ceiling. Shervish's interval-null p-value retains both tails. A finite
  # cluster reference must be used here as well as for the ordinary Wald
  # probability; otherwise the ETS verdict can be more liberal than the test
  # on which it is reported.
  out <- stats::pt((-a_cut - d) / se, df = df) +
    stats::pt((a_cut - d) / se, df = df)
  out[!is.finite(difference) | !is.finite(se) | se <= 0 |
        is.na(df) | df <= 0] <- NA_real_
  out
}

.ets_category <- function(difference, se, p, alpha = 0.05,
                          p_beyond = NULL) {
  a_cut <- 1.0 / ETS_DELTA_PER_LOGIT      # 0.43 as published
  c_cut <- 1.5 / ETS_DELTA_PER_LOGIT      # 0.64 as published
  d <- abs(difference)
  sig <- is.finite(p) & p < alpha
  if (is.null(p_beyond)) p_beyond <- .ets_p_beyond(difference, se)
  beyond <- is.finite(se) & se > 0 & is.finite(p_beyond) & p_beyond < alpha
  out <- ifelse(!sig | d <= a_cut, "A",
                ifelse(d >= c_cut & beyond, "C", "B"))
  # An unavailable test is not a nonsignificant test. The interval-null
  # probability is required only when it distinguishes B from C.
  unavailable <- !is.finite(difference) | !is.finite(se) | se <= 0 |
    !is.finite(p) | (sig & d >= c_cut & !is.finite(p_beyond))
  out[unavailable] <- NA_character_
  sign_c <- ifelse(!is.finite(difference) | out == "A", "",
                   ifelse(difference > 0, "+", "-"))
  ifelse(is.na(out), NA_character_, paste0(out, sign_c))
}

# EFRM has no scalar link-convergence flag. Each estimated connection between
# item sets carries its own status; all must converge before a result that
# depends on the linked scale is reported. A one-set model has no link edges.
.efrm_link_converged <- function(fit) {
  if (!inherits(fit, "rasch_efrm")) return(TRUE)
  edges <- fit$linking$alpha_edges
  n_sets <- if (is.data.frame(fit$alpha_table))
    nrow(fit$alpha_table) else NA_integer_
  if (is.null(edges) || !NROW(edges))
    return(is.finite(n_sets) && n_sets <= 1L)
  is.data.frame(edges) && "converged" %in% names(edges) &&
    all(edges$converged %in% TRUE)
}
