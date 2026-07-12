# Case study: party blocs and crisis concern as frames in paired comparisons
# ===========================================================================
# When respondents choose between political parties two at a time, the
# Bradley-Terry-Luce scale assumes every contest is judged in the same unit.
# The extended frame of reference model for paired comparisons relaxes
# exactly that: judge panels may judge with different discernment (panel
# units phi), and object sets may be calibrated in a different unit and
# origin than the common scale (set units alpha, origins kappa), with the
# linking identified from cross-set comparison outcomes alone.
#
# This script applies the model to the GermanParties2009 data shipped with
# the 'psychotools' package: 192 respondents in Tuebingen, June 2009, three
# months before the federal election, each choosing within all 15 pairs of
# five parties and the option of abstaining (Strobl, Wickelmaier & Zeileis,
# 2011, J. Educational and Behavioral Statistics 36, 135-153). Two frame
# questions are posed. Do respondents who feel personally affected by the
# 2008-9 economic crisis judge party contests with different discernment
# than those who do not (panels)? And do contests inside an ideological
# bloc -- left: Die Linke, Die Gruenen, SPD; right: CDU/CSU, FDP -- run in a
# different unit than the cross-bloc contests that define the common scale
# (sets)?
#
# The data also make the study a small design clinic: the right bloc has
# only one internal pair (CDU/CSU vs FDP) and the sample splits it nearly
# evenly, so that set carries no stable information about the panel-unit
# ratios and its own unit is weakly identified. The fit detects both,
# screens the set out of the phi reconciliation, reports a boundary-unstable
# alpha as NA, and says so in notes -- the honest answers to questions this
# design cannot support, and a lesson for comparative-judgment designs:
# set units want three or more objects per set with decisive internal
# contests.
#
# Run time: a few seconds.
# ===========================================================================
library(rasch)
if (!requireNamespace("psychotools", quietly = TRUE))
  stop("this case study needs the 'psychotools' package for its data: ",
       "install.packages('psychotools')")
data("GermanParties2009", package = "psychotools")
gp <- GermanParties2009

# one comparison per row; winner carries the choice, judge the respondent
m <- as.matrix(gp$preference)          # 192 x 15, +1 first named option wins
pair_names <- strsplit(colnames(m), ":", fixed = TRUE)
d6 <- do.call(rbind, lapply(seq_along(pair_names), function(j) {
  a <- pair_names[[j]][1]; b <- pair_names[[j]][2]
  data.frame(object_a = a, object_b = b,
             winner   = ifelse(m[, j] == 1, a, b),
             judge    = sprintf("R%03d", seq_len(nrow(m))),
             crisis   = as.character(gp$crisis),
             gender   = as.character(gp$gender),
             educ     = ifelse(gp$education == "5", "university", "school"),
             stringsAsFactors = FALSE)
}))
d <- d6[d6$object_a != "none" & d6$object_b != "none", ]   # parties only

# the equal-unit baseline ----------------------------------------------------
# One scale, judge-clustered errors. The aggregate preferences are perfectly
# transitive, and with abstention kept in as a sixth option it splits the
# scale: every party except Die Linke is preferred to not voting at all.
bt <- btl(d, "object_a", "object_b", "winner", judge = "judge")
print(bt$objects[order(-bt$objects$location),
                 c("object", "location", "se", "fit_resid")], digits = 3)
print(btl_transitivity(bt)$summary, digits = 3)
bt6 <- btl(d6, "object_a", "object_b", "winner", judge = "judge")
print(bt6$objects[order(-bt6$objects$location), c("object", "location", "se")],
      digits = 3)

# the frames: crisis-concern panels x ideological-bloc sets ------------------
blocs <- list(left = c("Linke", "Gruene", "SPD"), right = c("CDU/CSU", "FDP"))
crisis_of <- setNames(d$crisis[!duplicated(d$judge)],
                      d$judge[!duplicated(d$judge)])
set.seed(2009)
f <- btl_efrm(d, "object_a", "object_b", "winner", "judge",
              panels = crisis_of, object_sets = blocs, boot_reps = 80)
print(f)
print(f$frames, digits = 3)
print(f$objects[order(-f$objects$v), c("object", "set", "v", "se_v")],
      digits = 3)

# What the frames say. The common-scale order reproduces the equal-unit
# order (Gruene > SPD > CDU/CSU > FDP > Linke). Crisis-affected respondents
# judge party contests with a smaller unit than the unaffected (phi 0.87
# versus 1.15, a discernment ratio of 1.33) -- less decisively, not more --
# but the difference does not reach significance by either standard error
# (bootstrap z = 1.58; judge-clustered conditional z = 1.26; see below).
# The right bloc's origin is firmly below the left's on the common scale
# (kappa = -0.35, z = -5.6): in this university-town sample the left bloc
# is preferred wholesale, consistent with the equal-unit locations. The
# right bloc's unit alpha, by contrast, is weakly identified -- its only
# internal pair (CDU/CSU vs FDP, 55:45) provides almost no internal spread
# to compare against the cross-bloc scale -- and the fit says so twice: the
# set is screened out of the panel-unit reconciliation, and the bootstrap
# reports the alpha SE as NA after log(alpha) reached the boundary in a
# minority of replicates. The per-frame fit residuals are all |z| < 1, so
# the cross-frame convention (cross-set contests judged at phi_g on the
# common scale) is consistent with these data.

# bootstrap versus judge-clustered errors ------------------------------------
# The parametric bootstrap resamples comparisons independently at the fitted
# probabilities -- self-consistent, since the model has no judge parameter --
# while the conditional stage-one errors are judge-clustered and guard
# against extra-model dependence within a respondent's ten choices. They
# disagree modestly here (se log phi 0.089 versus 0.112); the cautious
# reading takes the larger, and the phi conclusion is the same under either.
fc <- btl_efrm(d, "object_a", "object_b", "winner", "judge",
               panels = crisis_of, object_sets = blocs,
               se_method = "conditional")
cat(sprintf("se(log phi): bootstrap %.3f, judge-clustered conditional %.3f\n",
            f$phi_table$se_log_phi[1], fc$phi_table$se_log_phi[1]))
cat(sprintf("phi ratio %.2f: z = %.2f (bootstrap), %.2f (conditional)\n",
            max(f$phi_table$phi) / min(f$phi_table$phi),
            abs(f$phi_table$z[1]), abs(fc$phi_table$z[1])))

# sensitivity: other panel definitions ---------------------------------------
# Neither gender nor education (university versus school-level, the sample's
# only broad split) defines panels with distinguishable units: the frame
# structure is not an artefact of how the panels are drawn, and the crisis
# contrast above is the largest of the three.
for (pv in c("gender", "educ")) {
  pmap <- setNames(d[[pv]][!duplicated(d$judge)], d$judge[!duplicated(d$judge)])
  pmap <- pmap[!is.na(pmap)]
  ds <- d[d$judge %in% names(pmap), ]
  set.seed(2009)
  fs <- btl_efrm(ds, "object_a", "object_b", "winner", "judge",
                 panels = pmap, object_sets = blocs, boot_reps = 80)
  cat(sprintf("%s: phi = %s; |z| = %.2f\n", pv,
              paste(sprintf("%.2f (%s)", fs$phi_table$phi,
                            fs$phi_table$panel), collapse = ", "),
              abs(fs$phi_table$z[1])))
}

# the visual: one unit per frame ---------------------------------------------
plot_btl_units(f)
