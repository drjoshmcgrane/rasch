# Fresh-seed confirmation of the extended conditional-bootstrap cells selected
# by the 100 x 99 screen. The helper definitions are loaded only after
# verifying their exact content hash. In addition to each
# rate, this study records the paired current-minus-bootstrap difference and
# its Monte Carlo standard error.
#
# The pinned hash is the helper as committed. Its reference-probability fix
# postdates the recorded results, whose rows name the earlier helper they ran
# against (813c7f6af9cbb06cac9a2b8c24ab07e5), so a rerun is distinguishable
# from them.

helper <- "tools/simval/studies/dif-conditional-bootstrap-extended.R"
helper_md5 <- "4624788f06cc401781b2cc774a246b8e"
if (!identical(unname(tools::md5sum(helper)), helper_md5))
  stop("extended-study helper hash does not match the committed version")
ex <- parse(file = helper)
for (i in seq_len(19L)) eval(ex[[i]], envir = environment())

NREP <- as.integer(Sys.getenv("SV_REPS", "200"))
NBOOT <- as.integer(Sys.getenv("SV_BOOT", "199"))
NCORE <- max(1L, as.integer(Sys.getenv("SV_CORES", "1")))

scenarios <- list(
  scenario("PCM", "imbalanced"),
  scenario("PCM", "imbalanced", "uniform", 0.6),
  scenario("PCM", "imbalanced", "nonuniform", 1.4),
  scenario("RSM", "imbalanced"),
  scenario("RSM", "imbalanced", "uniform", 0.6),
  scenario("RSM", "imbalanced", "nonuniform", 1.4),
  scenario("dichotomous", "multifactor", "nonuniform", 1.4))

rate_row <- function(sc, z, refused, quantity, col, field, note) {
  n <- nrow(z)
  val <- z[, col]
  args <- list(study = "dif-conditional-bootstrap-confirm",
    scenario = sc$name, quantity = quantity, n_reps = n,
    n_attempted = NREP, n_refused = refused, n_nonconv = 0L,
    effect = max(abs(sc$uniform), abs(sc$nonuniform)),
    notes = paste0(sc$model, "; ", sc$design, "; ", NBOOT,
      " conditional refits; helper md5 ", helper_md5, "; ", note))
  args[[field]] <- mean(val, na.rm = TRUE)
  if (field == "type1") args$mc_override <- list(
    type1 = stats::sd(val, na.rm = TRUE) / sqrt(sum(is.finite(val))))
  do.call(sv_row, args)
}

paired_row <- function(sc, z, refused, quantity, current, bootstrap, note) {
  dif <- z[, current] - z[, bootstrap]
  n <- sum(is.finite(dif))
  sv_row(study = "dif-conditional-bootstrap-confirm", scenario = sc$name,
    quantity = quantity, n_reps = n, n_attempted = NREP,
    n_refused = refused, n_nonconv = 0L,
    effect = max(abs(sc$uniform), abs(sc$nonuniform)),
    bias = mean(dif, na.rm = TRUE),
    mean_se = stats::sd(dif, na.rm = TRUE) / sqrt(n), se_ratio = NA_real_,
    notes = paste0("paired rate difference current minus bootstrap; mean_se ",
      "is its Monte Carlo SE; ", sc$model, "; ", sc$design, "; ", NBOOT,
      " conditional refits; helper md5 ", helper_md5, "; ", note))
}

confirm_condition <- function(sc) {
  seeds <- sample.int(.Machine$integer.max, NREP)
  z <- parallel::mclapply(seeds, one_replicate, sc = sc, mc.cores = NCORE,
                          mc.preschedule = TRUE)
  refused <- sum(vapply(z, is.null, TRUE))
  z <- z[!vapply(z, is.null, TRUE)]
  if (!length(z)) stop("no replicate completed in ", sc$name)
  z <- do.call(rbind, z)
  null <- sc$uniform == 0 && sc$nonuniform == 0
  if (null) rbind(
    rate_row(sc, z, refused, "current uniform item-wise Type I",
             "current_raw_uniform", "type1", "hybrid reference"),
    rate_row(sc, z, refused, "bootstrap uniform item-wise Type I",
             "bootstrap_raw_uniform", "type1", "empirical tail"),
    paired_row(sc, z, refused, "paired uniform Type I difference",
               "current_raw_uniform", "bootstrap_raw_uniform", "item-wise"),
    rate_row(sc, z, refused, "current non-uniform item-wise Type I",
             "current_raw_nonuniform", "type1", "residual-ANOVA reference"),
    rate_row(sc, z, refused, "bootstrap non-uniform item-wise Type I",
             "bootstrap_raw_nonuniform", "type1", "empirical tail"),
    paired_row(sc, z, refused, "paired non-uniform Type I difference",
               "current_raw_nonuniform", "bootstrap_raw_nonuniform", "item-wise"),
    rate_row(sc, z, refused, "current Holm FWER", "current_fwer",
             "familywise", "complete item-by-term family"),
    rate_row(sc, z, refused, "bootstrap minimum-p FWER", "bootstrap_fwer",
             "familywise", "complete item-by-term family"),
    paired_row(sc, z, refused, "paired FWER difference", "current_fwer",
               "bootstrap_fwer", "complete item-by-term family"))
  else rbind(
    rate_row(sc, z, refused, "current target power", "current_target",
             "power", paste("I3", sc$target_kind, "group term")),
    rate_row(sc, z, refused, "bootstrap target power", "bootstrap_target",
             "power", paste("I3", sc$target_kind, "group term")),
    paired_row(sc, z, refused, "paired target-power difference",
               "current_target", "bootstrap_target", "target term"),
    rate_row(sc, z, refused, "current unaffected-item FWER",
             "current_unaffected_fwer", "familywise", "I1-I2 and I4-I6"),
    rate_row(sc, z, refused, "bootstrap unaffected-item FWER",
             "bootstrap_unaffected_fwer", "familywise", "I1-I2 and I4-I6"),
    paired_row(sc, z, refused, "paired unaffected-item FWER difference",
               "current_unaffected_fwer", "bootstrap_unaffected_fwer",
               "I1-I2 and I4-I6"),
    if (sc$design == "multifactor")
      rate_row(sc, z, refused, "current non-target-factor FWER",
               "current_other_factor_fwer", "familywise", "all region terms")
      else NULL,
    if (sc$design == "multifactor")
      rate_row(sc, z, refused, "bootstrap non-target-factor FWER",
               "bootstrap_other_factor_fwer", "familywise", "all region terms")
      else NULL,
    if (sc$design == "multifactor")
      paired_row(sc, z, refused, "paired non-target-factor FWER difference",
                  "current_other_factor_fwer", "bootstrap_other_factor_fwer",
                  "all region terms") else NULL)
}

set.seed(8.43e7)
rows <- list()
for (i in seq_along(scenarios)) {
  message(sprintf("[%d/%d] %s", i, length(scenarios), scenarios[[i]]$name))
  rows[[i]] <- confirm_condition(scenarios[[i]])
  sv_write(do.call(rbind, rows), "dif-conditional-bootstrap-confirm")
}
