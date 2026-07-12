# Calibration battery for the package's inferential core.
#
# Every reported standard error is checked against the empirical sampling
# variability of its estimator on model-true data (the check that exposed the
# count-weighted BTL sandwich inflation). Run on demand:
#   Rscript tools/calibration.R
# Pass bands: SE-ratio in ~[0.85, 1.20] for the replicated paths; anchored
# translation equivalence to numerical precision.
#
# DESIGN RULE (learned the hard way): hold the TRUE parameters FIXED across
# replicates -- redrawing the truth each replicate folds between-replicate
# truth variation into the "empirical SD" and fakes a calibration failure.

suppressWarnings(pkgload::load_all("/Users/josh/Documents/Claude_Code/rmt", quiet=TRUE))
ok <- function(l,v) cat(sprintf("%-64s %s\n", l, v))

## A. pcml dichotomous item-location SEs, 150 reps (band should be tight)
L <- 10; Np <- 500
dtrue <- scale(seq(-2,2,length.out=L), scale=FALSE)[,1]
est <- ses <- matrix(NA, 150, L)
for (r in 1:150) { set.seed(5000+r)
  th <- rnorm(Np, 0, 1.3)
  X <- matrix(rbinom(Np*L,1,plogis(outer(th,dtrue,"-"))), Np, L)
  f <- pcml(X); est[r,] <- f$thr$tau; ses[r,] <- f$thr$se }
rA <- colMeans(ses)/apply(est,2,sd)
ok("A pcml dichot SE ratio per item (min..max)", sprintf("%.2f .. %.2f", min(rA), max(rA)))

## B. PCM threshold SEs, 100 reps
simP <- function(t, tau) { x<-0:length(tau); p<-exp(x*t-c(0,cumsum(tau))); p/sum(p) }
tt <- list(c(-1.2,-0.1), c(-0.6,0.5), c(-0.9,0.2), c(0.1,1.1), c(-0.4,0.8), c(-1.0,0.6))
tt <- lapply(tt, function(t) t - mean(unlist(tt)))
estB <- sesB <- matrix(NA, 100, sum(lengths(tt)))
for (r in 1:100) { set.seed(6000+r)
  th <- rnorm(600, 0, 1.4)
  X <- sapply(seq_along(tt), function(i) vapply(th, function(t)
    sample(0:2, 1, prob=simP(t, tt[[i]])), 0L))
  colnames(X) <- sprintf("P%d", 1:6)
  f <- rasch(X, model="PCM")
  estB[r,] <- f$thresholds$tau; sesB[r,] <- f$thresholds$se }
rB <- colMeans(sesB)/apply(estB,2,sd)
ok("B PCM threshold SE ratio (min..max)", sprintf("%.2f .. %.2f", min(rB), max(rB)))

## C. anchored rasch at distant origins
set.seed(3); th <- rnorm(600,0,1.3)
X <- matrix(rbinom(600*8,1,plogis(outer(th, seq(-1.5,1.5,length.out=8), "-"))),600,8)
colnames(X) <- sprintf("I%d",1:8)
free <- rasch(X)
fv <- free$items$location
for (delta in c(0, 6)) {
  anc <- data.frame(item=c("I1","I8"), k=1, tau=fv[c(1,8)] + delta)
  fa <- rasch(X, anchors=anc)
  errI <- max(abs(fa$items$location - (fv + delta)))
  errP <- max(abs(fa$person$theta - (free$person$theta + delta)), na.rm=TRUE)
  ok(sprintf("C anchored rasch +%d: conv, item err, person err", delta),
     sprintf("%s, %.2e, %.2e", fa$est$converged, errI, errP))
}

## D. MFRM severity SE calibration, fixed truth, 120 reps
simP2 <- function(t, tau) { x<-0:length(tau); p<-exp(x*t-c(0,cumsum(tau))); p/sum(p) }
lam <- c(R1=-0.7, R2=-0.2, R3=0.3, R4=0.6)
del <- c(I1=-0.8, I2=-0.2, I3=0.3, I4=0.7)
base_tau <- c(-1.1, 0, 1.1)
estD <- sesD <- matrix(NA, 120, 4)
for (r in 1:120) { set.seed(7000+r)
  th <- rnorm(70, 0, 1.2)
  g <- expand.grid(p=1:70, i=1:4, rt=1:4)
  sc <- vapply(seq_len(nrow(g)), function(k)
    sample(0:3, 1, prob=simP2(th[g$p[k]], base_tau + del[g$i[k]] + lam[g$rt[k]])), 0L)
  d <- data.frame(person=sprintf("P%03d",g$p), item=names(del)[g$i],
                  rater=names(lam)[g$rt], score=sc)
  mf <- rasch_mfrm(d, person="person", item="item", score="score", facets="rater")
  fe <- mf$facet_effects$rater
  estD[r,] <- fe$severity[match(names(lam), fe$level)]
  sesD[r,] <- fe$se[match(names(lam), fe$level)]
}
rD <- colMeans(sesD)/apply(estD,2,sd)
ok("D MFRM severity SE ratio, fixed truth (min..max)",
   sprintf("%.2f .. %.2f", min(rD), max(rD)))

## E. dif_size resolved-difference SE, 100 reps (planted 0.8 uniform DIF)
diffs <- sesE <- numeric(100)
for (r in 1:100) { set.seed(8000+r)
  d <- simulate_rasch(500, 10, dif=list(items="I05", uniform=0.8), n_groups=2, seed=8000+r)
  f <- rasch(d, id="id", factors="group")
  ds <- dif_size(f, "I05", by="group")
  diffs[r] <- ds$pairs$difference; sesE[r] <- ds$pairs$se }
ok("E dif_size: empSD vs mean reported SE", sprintf("%.3f vs %.3f (ratio %.2f)",
   sd(diffs), mean(sesE), mean(sesE)/sd(diffs)))
ok("E dif_size: mean recovered diff (planted 0.8)", sprintf("%.3f", mean(abs(diffs))))

## F. person WLE SEs at fixed true theta, complete data
set.seed(9); dv <- seq(-2,2,length.out=12)
for (th0 in c(0, 1.5)) {
  ests <- numeric(400); se_rep <- NA
  # item params treated as known-ish: use one large calibration, then repeated persons
  Xc <- matrix(rbinom(1500*12,1,plogis(outer(rnorm(1500,0,1.4),dv,"-"))),1500,12)
  colnames(Xc) <- sprintf("I%02d",1:12)
  fc <- rasch(Xc)
  st <- person_wle(fc$tau_list)
  for (r in 1:400) { set.seed(9000+r)
    x <- rbinom(12,1,plogis(th0-dv)); s <- sum(x)
    ests[r] <- st$theta[as.character(s)] }
  # reported SE at the MEAN raw score
  s_typ <- as.character(round(mean(vapply(1:400, function(r){set.seed(9000+r); sum(rbinom(12,1,plogis(th0-dv)))},0))))
  ok(sprintf("F WLE at theta=%.1f: empSD vs reported SE", th0),
     sprintf("%.3f vs %.3f", sd(ests, na.rm=TRUE), st$se[s_typ]))
}
