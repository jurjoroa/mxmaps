# ---------------------------------------------------------------------------
# Purpose: repair BROKEN adjacency pairs one at a time in a large state, by a
#          directed local site move, accepting only measured improvement.
#
#          WHY THIS EXISTS SEPARATELY FROM optimise_state.R.
#          optimise_state.R probes random local perturbations. That found the
#          optimum for Baja California (6 cells, 8 pairs, 75% -> 100% in 40
#          tries) and is hopeless for Oaxaca: 570 cells is 1140 degrees of
#          freedom and 375 broken pairs, so a random proposal almost never
#          lands on a defect. The eight states with n > 100 hold 4362 of the
#          6288 national pairs -- 69% -- so a method that only works on small
#          states does not matter.
#
#          THE TWO CHANGES THAT MAKE n=570 AFFORDABLE.
#          1 WARM START. A site move perturbs the optimal weights slightly, so
#            Newton reconverges in 1-2 iterations from the previous w instead of
#            4-5 from w = 0. Cold-starting every candidate was most of the cost.
#          2 DIRECTED PROPOSALS. Enumerate the broken pairs and, for each, try
#            moving those two sites toward each other. 375 aimed attempts beat
#            thousands of blind ones, and each has a concrete hypothesis: this
#            pair is broken because these two cells do not reach each other.
#
#          THE ACCEPTANCE RULE IS THE WHOLE SAFETY ARGUMENT. A proposal is kept
#          only if the TOTAL kept-pair count rises and every constraint still
#          holds. unpin_sites.R computed a plausible descent direction over all
#          2469 sites at once, applied it unconditionally, and lost 1.36
#          cohesion points before starving the solver. Repairing one pair can
#          break a neighbouring one, so nothing may be trusted pairwise -- only
#          the global count decides.
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: repair-state v1
# Date:    28-August-2026
#
# Usage:   Rscript repair_state.R <state> <passes> [out_sites.csv]
#          REP_LEASH=<f>   max site displacement, in state cell widths (0.35)
#          REP_TRIALS=<n>  proposals attempted per broken pair (4)
#          REP_SEED=<n>    RNG seed (1)
#          REP_KNN=<k>     clip-set size during the search (120)
#          REP_SITES=<csv> start from these sites instead of true centroids
# ---------------------------------------------------------------------------

Sys.setenv(V11_LIB = "1")
if (!nzchar(Sys.getenv("V11_GAMMA_TAB")))
  Sys.setenv(V11_GAMMA_TAB = "data/gamma_by_state.csv")
suppressPackageStartupMessages(suppressMessages(
  source("/Users/jorgeroa/Documents/GitHub/mxmaps/ot_municipios_v11_2.R")))
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript repair_state.R <state> <passes> [out.csv]")
S <- sprintf("%02d", as.integer(args[1]))
PASSES <- as.integer(args[2])
OUT <- if (length(args) >= 3) args[3] else NA_character_
.n <- function(nm, d) { v <- Sys.getenv(nm); if (nzchar(v)) as.numeric(v) else d }
LEASH_F <- .n("REP_LEASH", 0.35)
TRIALS <- as.integer(.n("REP_TRIALS", 4))
set.seed(as.integer(.n("REP_SEED", 1)))
TINY <- 0.25
MIN_EQ <- 0.45

muns <- mxmun_sf[mxmun_sf$state_code == S, ]
n <- nrow(muns)
state_clip <- st_make_valid(state_sf$geometry[state_sf$state_code == S])
A_dom <- poly_area(state_clip)
P0 <- st_coordinates(st_transform(
  st_centroid(st_transform(muns, 6372), of_largest_polygon = TRUE), 4326))[, 1:2,
                                                                          drop = FALSE]
dup <- duplicated(round(P0, 10))
if (any(dup)) P0[dup, ] <- P0[dup, ] + 1e-7

gam <- state_areas$gamma[state_areas$state_code == S]
ta <- vapply(seq_len(n), function(i) poly_area(muns$geometry[i]), numeric(1))
tgt <- pmax(ta, 1e-12)^gam
tgt <- tgt / sum(tgt) * A_dom
share <- A_dom / n
for (r in 1:3) {
  tgt <- pmax(pmin(tgt, 2.2 * share), 0.45 * share)
  tgt <- tgt / sum(tgt) * A_dom
}
bb <- st_bbox(state_clip)
pad <- 0.6 * max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"])
sq <- matrix(c(bb["xmin"] - pad, bb["ymin"] - pad, bb["xmax"] + pad, bb["ymin"] - pad,
               bb["xmax"] + pad, bb["ymax"] + pad, bb["xmin"] - pad, bb["ymax"] + pad),
             ncol = 2, byrow = TRUE)

pm <- readRDS("data/true_pairs_meta.rds")
ps <- pm[pm$state_code == S, ]
li <- match(mxmun_sf$region[ps$i], muns$region)
lj <- match(mxmun_sf$region[ps$j], muns$region)
ok <- !is.na(li) & !is.na(lj)
PI <- li[ok]; PJ <- lj[ok]
NP <- length(PI)
cw <- sqrt(A_dom / n)
LEASH <- LEASH_F * cw

# REP_KNN sizes the clip set used DURING THE SEARCH. build_pd() cost scales
# with it, and the search only has to RANK candidates -- the winner is re-solved
# by the real engine at its own KNN before anything ships. Oaxaca at KNN=120
# costs ~12s per evaluation, which is 5 hours for one pass over 375 broken
# pairs; at 40 it is affordable. A too-small clip set makes cells too large and
# can overstate adjacency, which is exactly why the final verification is not
# optional.
SKNN <- as.integer(.n("REP_KNN", 120))
mk_knn <- function(P) {
  dm <- as.matrix(dist(P))
  kk <- min(SKNN, n - 1L)
  lapply(seq_len(n), function(i) order(dm[i, ])[2:(kk + 1)])
}

# Warm-started damped Newton. Ranking only -- a winner is re-solved by the real
# engine through V11_SITES before it is believed.
solve_w <- function(P, knn, w0 = NULL, maxit = 60, tol = 0.01) {
  dm <- as.matrix(dist(P))
  w <- if (is.null(w0)) rep(0, n) else w0 - mean(w0)
  b <- build_pd(P, w, knn, sq, state_clip)
  for (it in seq_len(maxit)) {
    ar <- b$ar
    live <- ar > 0
    if (all(live) && max(abs(ar / tgt - 1)) < tol) break
    res <- tgt - ar
    L <- matrix(0, n, n)
    for (i in seq_len(n)) {
      fl <- b$faces[[i]]
      if (!length(fl)) next
      js <- as.integer(names(fl))
      cij <- fl / (2 * pmax(dm[i, js], 1e-12))
      L[i, js] <- L[i, js] - cij
      L[i, i] <- L[i, i] + sum(cij)
    }
    L <- (L + t(L)) / 2
    dead <- which(!live | diag(L) <= 0)
    if (length(dead)) {
      md <- mean(diag(L)[diag(L) > 0], na.rm = TRUE)
      for (i in dead) L[i, i] <- max(L[i, i], if (is.finite(md)) md else 1)
    }
    reg <- mean(abs(diag(L))) * 1e-8
    dw <- tryCatch(solve(L + matrix(1 / n, n, n) + diag(reg, n), res),
                   error = function(e) NULL)
    if (is.null(dw) || any(!is.finite(dw))) dw <- res / pmax(diag(L), 1e-12)
    dw <- dw - mean(dw)
    rn <- sqrt(sum((ar - tgt)^2))
    tau <- 1
    stepped <- FALSE
    for (bt in 1:25) {
      wn <- w + tau * dw
      wn <- wn - mean(wn)
      bn <- build_pd(P, wn, knn, sq, state_clip)
      if (all(bn$ar > 0) && sqrt(sum((bn$ar - tgt)^2)) < rn) {
        w <- wn; b <- bn; stepped <- TRUE; break
      }
      tau <- tau / 2
    }
    if (!stepped) break
  }
  list(b = b, w = w)
}

evaluate <- function(P, w0 = NULL) {
  knn <- mk_knn(P)
  r <- solve_w(P, knn, w0)
  b <- r$b
  if (any(b$ar <= 0)) return(NULL)
  if (any(vapply(b$cells, is.null, TRUE))) return(NULL)
  g <- st_transform(st_sfc(b$cells, crs = 4326), 6372)
  a <- as.numeric(st_area(g))
  r_eq <- a / mean(a)
  if (min(r_eq) < MIN_EQ || any(r_eq < TINY)) return(NULL)
  nb <- st_is_within_distance(g, g, dist = 1, sparse = TRUE)
  surv <- mapply(function(i, j) j %in% nb[[i]], PI, PJ)
  list(kept = sum(surv), surv = surv, w = r$w, eq_min = min(r_eq))
}

# Writing sites is a FUNCTION so it can checkpoint after every pass. The first
# version only wrote at the very end, and when a 3-hour wave was interrupted it
# lost every accepted repair -- Edomex had already reached +6.92 pp and Chiapas
# +5.94, and both were unrecoverable because only the scores had been printed.
# A long search must be restartable from its own last good state.
write_sites <- function(P, path) {
  allP <- st_coordinates(st_transform(
    st_centroid(st_transform(mxmun_sf, 6372), of_largest_polygon = TRUE),
    4326))[, 1:2, drop = FALSE]
  sd_ <- Sys.getenv("REP_SITES")
  if (nzchar(sd_) && file.exists(sd_)) {
    sv <- read.csv(sd_, colClasses = c(region = "character"))
    mm <- match(mxmun_sf$region, sv$region)
    good <- !is.na(mm)
    allP[good, ] <- cbind(sv$lon[mm[good]], sv$lat[mm[good]])
  }
  allP[match(muns$region, mxmun_sf$region), ] <- P
  tmp <- paste0(path, ".tmp")
  write.csv(data.frame(region = mxmun_sf$region, lon = allP[, 1], lat = allP[, 2]),
            tmp, row.names = FALSE)
  file.rename(tmp, path)   # atomic, so a kill mid-write cannot corrupt it
}

# Resume: REP_SITES seeds the SEARCH as well as preserving other states, so a
# killed run continues from its last checkpoint instead of restarting. P0 stays
# the TRUE CENTROIDS regardless -- it is the leash anchor, and measuring
# displacement from a previous checkpoint instead would let the leash walk.
P_start <- P0
{
  sd0 <- Sys.getenv("REP_SITES")
  if (nzchar(sd0) && file.exists(sd0)) {
    sv0 <- read.csv(sd0, colClasses = c(region = "character"))
    m0 <- match(muns$region, sv0$region)
    if (!any(is.na(m0))) {
      P_start <- cbind(sv0$lon[m0], sv0$lat[m0])
      moved0 <- sum(sqrt(rowSums((P_start - P0)^2)) > 1e-12)
      cat(sprintf("resuming from %s (%d/%d sites already moved)\n",
                  basename(sd0), moved0, n))
    }
  }
}

cur <- evaluate(P_start)
if (is.null(cur)) stop("baseline violates a constraint")
cat(sprintf("state %s: n=%d gamma=%.2f | baseline kept %d/%d = %.2f%% | eq_min %.3f\n",
            S, n, gam, cur$kept, NP, 100 * cur$kept / NP, cur$eq_min))
start_kept <- cur$kept
P <- P_start
t0 <- Sys.time()

for (pass in seq_len(PASSES)) {
  broken <- which(!cur$surv)
  if (!length(broken)) break
  fixed <- 0
  # hardest first: a pair whose endpoints have the most other broken pairs is
  # the one most likely to be sitting in a genuinely bad local configuration
  deg_broken <- tabulate(c(PI[broken], PJ[broken]), nbins = n)
  broken <- broken[order(-(deg_broken[PI[broken]] + deg_broken[PJ[broken]]))]
  att <- 0
  for (bi in broken) {
    i <- PI[bi]; j <- PJ[bi]
    att <- att + 1
    if (att %% 25 == 0)
      cat(sprintf("    [pass %d] %d/%d attempted, +%d, kept %d (%.0f min)\n",
                  pass, att, length(broken), fixed, cur$kept,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    if (cur$surv[bi]) next            # an earlier repair may have fixed it
    for (tr in seq_len(TRIALS)) {
      Q <- P
      d <- Q[j, ] - Q[i, ]
      L2 <- sqrt(sum(d^2))
      if (L2 < 1e-12) next
      u <- d / L2
      m <- LEASH * runif(1, 0.15, 1.0) * switch(as.character(tr %% 4 + 1),
             "1" = 1, "2" = 0.5, "3" = 0.75, "4" = 0.3)
      jit <- LEASH * 0.15
      # move both endpoints toward each other, with a little jitter so repeated
      # trials on the same pair are not identical
      Q[i, ] <- Q[i, ] + m * u + rnorm(2, 0, jit)
      Q[j, ] <- Q[j, ] - m * u + rnorm(2, 0, jit)
      dd <- sqrt((Q[, 1] - P0[, 1])^2 + (Q[, 2] - P0[, 2])^2)
      ov <- dd > LEASH
      if (any(ov)) {
        sc <- LEASH / dd[ov]
        Q[ov, 1] <- P0[ov, 1] + (Q[ov, 1] - P0[ov, 1]) * sc
        Q[ov, 2] <- P0[ov, 2] + (Q[ov, 2] - P0[ov, 2]) * sc
      }
      s <- evaluate(Q, cur$w)
      if (!is.null(s) && s$kept > cur$kept) {
        P <- Q; cur <- s; fixed <- fixed + 1
        break
      }
    }
  }
  cat(sprintf("pass %d: +%d pairs -> kept %d/%d = %.2f%% (%.0f min elapsed)\n",
              pass, fixed, cur$kept, NP, 100 * cur$kept / NP,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (!is.na(OUT) && cur$kept > start_kept) {
    write_sites(P, OUT)
    cat(sprintf("  checkpointed -> %s\n", basename(OUT)))
  }
  if (fixed == 0) break
}

cat(sprintf("FINAL state %s: %d -> %d of %d pairs | %.2f%% -> %.2f%% (%+.2f pp)\n",
            S, start_kept, cur$kept, NP, 100 * start_kept / NP,
            100 * cur$kept / NP, 100 * (cur$kept - start_kept) / NP))

if (!is.na(OUT) && cur$kept > start_kept) {
  write_sites(P, OUT)
  cat("Wrote", OUT, "\n")
}
