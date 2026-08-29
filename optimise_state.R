# ---------------------------------------------------------------------------
# Purpose: search SITE POSITIONS for one state, keeping the prescribed-area
#          guarantee, and accept a move only if it measurably improves.
#
#          WHY SITE POSITIONS, AND WHY THIS IS NOT unpin_sites.R AGAIN.
#          Aurenhammer-Hoffmann-Aronov: for any distinct sites and any positive
#          target areas summing to the domain, weights EXIST making the power
#          diagram's areas exactly the targets. So moving a site costs nothing
#          in area accuracy -- only displacement. unpin_sites.R exploited that
#          with a blind gradient step over every broken pair and LOST 1.36
#          cohesion points, then starved sites into solver failure. The fix is
#          not a better direction, it is a better acceptance rule: propose,
#          re-solve, and KEEP ONLY IF THE SCORE IMPROVED and every constraint
#          still holds. A hill-climb that rejects regressions cannot make the
#          layout worse than its starting point, by construction.
#
#          TWO OBJECTIVES.
#          cohesion -- maximise kept true intra-state adjacencies. Aimed at
#            small states sitting below their k=6 ceiling with no high-degree
#            excuse: Baja California keeps 6 of 8 pairs with a 100% ceiling.
#          lex      -- cohesion first, then shape among equals. Prefer this.
#          shape    -- maximise mean isoperimetric quotient (a regular hexagon
#            is 0.907, a square 0.785) subject to cohesion not dropping. Aimed
#            at Oaxaca, where 146 of 570 cells (25.6%) are quadrilaterals and
#            median elongation is 1.72. NOTE what this cannot do: a power cell
#            has exactly as many sides as it has neighbours, and Oaxaca's
#            municipios have median degree 5, so the side COUNT is the dual of
#            the real adjacency graph and is not a free parameter. Only the
#            regularity of the polygon given that count is.
#
#          THE SEARCH SOLVER IS A DUPLICATE, ON PURPOSE AND ONLY FOR RANKING.
#          The damped Newton below re-implements the engine's loop compactly so
#          hundreds of candidates can be scored in ONE session (a process launch
#          costs ~165s of setup for under a second of solving). build_pd() and
#          clip_lab() are SOURCED from the engine via V11_LIB, never copied. Any
#          winner must then be re-solved by the real unmodified engine through
#          V11_SITES before it is believed -- this script's numbers rank
#          candidates, they do not ship.
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: optimise-state v1
# Date:    28-August-2026
#
# Usage:   Rscript optimise_state.R <state> <objective> <iters> [out_sites.csv]
#          OPT_LEASH=<f>   max site displacement, in state cell widths (0.50)
#          OPT_SEED=<n>    RNG seed (1)
# ---------------------------------------------------------------------------

Sys.setenv(V11_LIB = "1")
if (!nzchar(Sys.getenv("V11_GAMMA_TAB")))
  Sys.setenv(V11_GAMMA_TAB = "data/gamma_by_state.csv")
suppressPackageStartupMessages(suppressMessages(
  source("/Users/jorgeroa/Documents/GitHub/mxmaps/ot_municipios_v11_2.R")))
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("usage: Rscript optimise_state.R <state> <cohesion|shape|lex> <iters> [out.csv]")
S <- sprintf("%02d", as.integer(args[1]))
OBJ <- args[2]
ITERS <- as.integer(args[3])
OUT <- if (length(args) >= 4) args[4] else NA_character_
LEASH_F <- { v <- Sys.getenv("OPT_LEASH"); if (nzchar(v)) as.numeric(v) else 0.50 }
set.seed({ v <- Sys.getenv("OPT_SEED"); if (nzchar(v)) as.integer(v) else 1L })

TINY <- 0.25
MIN_EQ <- 0.45

# ---- per-state setup, mirroring the engine exactly ------------------------
muns <- mxmun_sf[mxmun_sf$state_code == S, ]
n <- nrow(muns)
sp <- state_sf$geometry[state_sf$state_code == S]
state_clip <- st_make_valid(sp)
A_dom <- poly_area(state_clip)
P0 <- st_coordinates(st_transform(
  st_centroid(st_transform(muns, 6372), of_largest_polygon = TRUE), 4326))[, 1:2,
                                                                          drop = FALSE]
dup <- duplicated(round(P0, 10))
if (any(dup)) P0[dup, ] <- P0[dup, ] + 1e-7

gam <- state_areas$gamma[state_areas$state_code == S]
ta <- vapply(seq_len(n), function(i) poly_area(muns$geometry[i]), numeric(1))
wt <- pmax(ta, 1e-12)^gam
tgt <- wt / sum(wt) * A_dom
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

# true intra-state pairs for this state, from the shared cache
pm <- readRDS("data/true_pairs_meta.rds")
ps <- pm[pm$state_code == S, ]
li <- match(mxmun_sf$region[ps$i], muns$region)
lj <- match(mxmun_sf$region[ps$j], muns$region)
keep <- !is.na(li) & !is.na(lj)
PI <- li[keep]; PJ <- lj[keep]
cat(sprintf("state %s: n=%d gamma=%.2f | %d true intra-state pairs\n",
            S, n, gam, length(PI)))

cw_deg <- sqrt(A_dom / n)          # cell width in the units P is expressed in
LEASH <- LEASH_F * cw_deg

mk_knn <- function(P) {
  dm <- as.matrix(dist(P))
  kk <- min(120L, n - 1L)
  if (n > 1) lapply(seq_len(n), function(i) order(dm[i, ])[2:(kk + 1)]) else list(integer(0))
}

# ---- damped Newton on the OT dual (see header: ranking only) --------------
solve_w <- function(P, knn, maxit = 120, tol = 0.01) {
  dm <- as.matrix(dist(P))
  w <- rep(0, n)
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
  b
}

score <- function(P) {
  knn <- mk_knn(P)
  b <- solve_w(P, knn)
  if (any(b$ar <= 0)) return(NULL)                 # a starved site is a reject
  g <- st_sfc(b$cells[!vapply(b$cells, is.null, TRUE)], crs = 4326)
  if (length(g) != n) return(NULL)
  gm <- st_transform(g, 6372)
  a <- as.numeric(st_area(gm))
  r_eq <- a / mean(a)
  if (min(r_eq) < MIN_EQ || any(r_eq < TINY)) return(NULL)   # legibility floor
  nb <- st_is_within_distance(gm, gm, dist = 1, sparse = TRUE)
  kept <- sum(mapply(function(i, j) j %in% nb[[i]], PI, PJ))
  per <- vapply(seq_len(n), function(i) {
    co <- st_coordinates(gm[i]); co <- co[!is.na(co[, 1]), , drop = FALSE]
    x <- co[, 1]; y <- co[, 2]
    sum(sqrt((x - c(x[-1], x[1]))^2 + (y - c(y[-1], y[1]))^2))
  }, 0)
  list(kept = kept, coh = 100 * kept / length(PI),
       iq = mean(4 * pi * a / per^2), eq_min = min(r_eq))
}

# OPT_SEED_SITES starts the search from an existing site table instead of the
# true centroids. This is what makes the objectives composable in PHASES, which
# is the only way lexicographic optimisation actually works: a one-pass greedy
# rule that wants "more pairs, else rounder" gets stuck, because reaching the
# better pair count usually requires passing through a configuration that
# improves neither. Measured on Baja California: one-pass lex stalled at 7/8
# with IQ 0.345, while cohesion-only reached 8/8 at IQ 0.175. Run phase 1
# (cohesion) to completion, then phase 2 (shape) seeded from its winner.
P_start <- P0
seed_f <- Sys.getenv("OPT_SEED_SITES")
if (nzchar(seed_f)) {
  if (!file.exists(seed_f)) stop("OPT_SEED_SITES not found: ", seed_f)
  sv <- read.csv(seed_f, colClasses = c(region = "character"))
  m <- match(muns$region, sv$region)
  if (any(is.na(m))) stop("seed sites do not cover state ", S)
  P_start <- cbind(sv$lon[m], sv$lat[m])
  cat(sprintf("seeded from %s\n", basename(seed_f)))
}

base <- score(P_start)
if (is.null(base)) stop("baseline itself violates a constraint -- check setup")
cat(sprintf("baseline: kept %d/%d = %.2f%% | IQ %.4f | eq_min %.3f\n",
            base$kept, length(PI), base$coh, base$iq, base$eq_min))

better <- function(new, old) {
  if (is.null(new)) return(FALSE)
  if (OBJ == "cohesion") return(new$kept > old$kept)
  if (OBJ == "shape")
    # IQ must rise and cohesion must NOT fall
    return(new$iq > old$iq + 1e-6 && new$kept >= old$kept)
  # "lex": cohesion first, then shape among equals. A pure cohesion objective
  # is blind to what it spends: on Baja California it found 8/8 pairs (up from
  # 6/8) while driving mean IQ from 0.276 to 0.175, a 37% loss of shape
  # regularity nobody asked for. Lexicographic ordering takes the same
  # adjacency win and then keeps improving roundness at no cost to it.
  if (new$kept > old$kept) return(TRUE)
  new$kept == old$kept && new$iq > old$iq + 1e-6
}

P <- P_start
best <- base
acc <- 0
for (it in seq_len(ITERS)) {
  # perturb a random subset; smaller states get whole-vector proposals
  # Proposal size has to shrink with n. At n=570 a 5% proposal perturbs ~28
  # sites at once and nearly every such move breaks something, so the search
  # never accepts. Large states get a LOCAL proposal: one seed municipio plus a
  # couple of its nearest neighbours, which is the scale at which a single
  # broken pair can actually be repaired.
  if (n <= 12) {
    k <- max(1L, rbinom(1, n, 0.6))
    idx <- sample.int(n, k)
  } else if (n <= 60) {
    k <- sample(1:3, 1)
    idx <- sample.int(n, k)
  } else {
    seed_i <- sample.int(n, 1)
    nn <- order((P[, 1] - P[seed_i, 1])^2 + (P[, 2] - P[seed_i, 2])^2)
    idx <- nn[seq_len(sample(1:3, 1))]
  }
  k <- length(idx)
  step <- LEASH * runif(1, 0.05, 0.45)
  Q <- P
  ang <- runif(k, 0, 2 * pi)
  Q[idx, 1] <- Q[idx, 1] + step * cos(ang)
  Q[idx, 2] <- Q[idx, 2] + step * sin(ang)
  d <- sqrt((Q[, 1] - P0[, 1])^2 + (Q[, 2] - P0[, 2])^2)
  over <- d > LEASH
  if (any(over)) {
    sc <- LEASH / d[over]
    Q[over, 1] <- P0[over, 1] + (Q[over, 1] - P0[over, 1]) * sc
    Q[over, 2] <- P0[over, 2] + (Q[over, 2] - P0[over, 2]) * sc
  }
  s <- score(Q)
  if (better(s, best)) {
    P <- Q; best <- s; acc <- acc + 1
    cat(sprintf("  it %4d ACCEPT kept %d/%d = %.2f%% | IQ %.4f | eq_min %.3f\n",
                it, s$kept, length(PI), s$coh, s$iq, s$eq_min))
  }
}
cat(sprintf("\n%d/%d proposals accepted\n", acc, ITERS))
cat(sprintf("FINAL   kept %d/%d = %.2f%% (was %.2f%%) | IQ %.4f (was %.4f)\n",
            best$kept, length(PI), best$coh, base$coh, best$iq, base$iq))

if (!is.na(OUT) && acc > 0) {
  # write the FULL national site table so the engine can be run over all states
  # with only this one changed
  allP <- st_coordinates(st_transform(
    st_centroid(st_transform(mxmun_sf, 6372), of_largest_polygon = TRUE),
    4326))[, 1:2, drop = FALSE]
  m <- match(muns$region, mxmun_sf$region)
  allP[m, ] <- P
  write.csv(data.frame(region = mxmun_sf$region, lon = allP[, 1], lat = allP[, 2]),
            OUT, row.names = FALSE)
  cat("Wrote", OUT, "\n")
}
