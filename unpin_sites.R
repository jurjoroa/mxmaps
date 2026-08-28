# ---------------------------------------------------------------------------
# Purpose: compute ONE round of adjacency-directed site movement for the v11
#          optimal-transport engine, writing a region,lon,lat csv it can read
#          back through V11_SITES.
#
#          WHY MOVING SITES IS FREE, AND WHY v8 STILL FAILED DOING IT.
#          Aurenhammer-Hoffmann-Aronov: for ANY distinct sites and ANY positive
#          target areas summing to the domain, weights EXIST making the power
#          diagram's cell areas exactly the targets. So a site move costs
#          nothing in area accuracy -- the Newton solve re-establishes exact
#          areas afterwards. What it costs is DISPLACEMENT from the true
#          centroid, and that is the only thing this script has to ration.
#          v8's CVT moved sites to equalise AREAS, which is the weights' job,
#          and paid for it in the adjacency the true centroids had for free
#          (Delaunay of true centroids alone reproduces 83.1% of real
#          adjacency). Here the move is directed at ADJACENCY and leashed.
#
#          WHAT IT AIMS AT. Measured over v11.2's 6288 pairs:
#            - 37.3% of broken pairs have a third municipio's territory on the
#              straight centroid-to-centroid line; the modelled counterfactual
#              cure is 78.5% -> 85.6%.
#            - a cell gets ~6 usable faces whether its municipio has 4 true
#              neighbours or 22, and 11.5% of all faces are spent on
#              cross-state contact. Capacity is MISALLOCATED, not scarce:
#              1897 surplus faces against 1493 shortfalls.
#          So: pull broken true-adjacent pairs together, and push apart cells
#          that touch when their municipios do NOT, which frees a face.
#
#          THE LEASH IS ASYMMETRIC ON PURPOSE. Border cells are displaced 34.2
#          km against 45.2 km for interior cells yet survive 11 pp worse at
#          matched degree -- they are sitting on unused positional slack. They
#          get a longer leash. Border membership is computed from each
#          municipio's own same-state shared perimeter, NOT from st_union +
#          st_boundary, which over-flags 2339 of 2469 municipios on union
#          slivers.
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: unpin v1
# Date:    27-August-2026
#
# Usage:   Rscript unpin_sites.R <artifact.rds> <out_sites.csv> [in_sites.csv]
#          UNPIN_STEP=<f>        step per round, in state cell widths (0.25)
#          UNPIN_LEASH_INT=<f>   interior leash, in state cell widths (0.60)
#          UNPIN_LEASH_BDY=<f>   state-boundary leash, same units (1.00)
#          UNPIN_REPEL=<f>       weight on false-adjacency repulsion (0.30)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
source("hexlayout_common.R")

MX_CRS <- 6372
.env <- function(nm, d) {
  v <- Sys.getenv(nm)
  if (nzchar(v)) as.numeric(v) else d
}
STEP <- .env("UNPIN_STEP", 0.25)
LEASH_INT <- .env("UNPIN_LEASH_INT", 0.60)
LEASH_BDY <- .env("UNPIN_LEASH_BDY", 1.00)
REPEL <- .env("UNPIN_REPEL", 0.30)
BDY_TOL <- 0.05  # >5% of perimeter not shared with a same-state neighbour

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("usage: Rscript unpin_sites.R <artifact.rds> <out.csv> [in.csv]")
art <- args[1]
out_csv <- args[2]
in_csv <- if (length(args) >= 3) args[3] else NA_character_

pm <- readRDS("data/true_pairs_meta.rds")
mx <- build_mxmun_sf()
mx_m <- st_transform(mx, MX_CRS)
n <- nrow(mx)

# True centroids: the anchor the leash is measured from, in metres.
cen <- st_coordinates(st_centroid(mx_m, of_largest_polygon = TRUE))[, 1:2]

# Current sites. Round 1 starts at the true centroids (the v11.2 pin).
cur <- cen
if (!is.na(in_csv) && file.exists(in_csv)) {
  sv <- read.csv(in_csv, colClasses = c(region = "character"))
  pts <- st_transform(st_as_sf(sv, coords = c("lon", "lat"), crs = 4326),
                      MX_CRS)
  m <- match(mx$region, sv$region)
  ok <- !is.na(m)
  cc <- st_coordinates(pts)
  cur[ok, ] <- cc[m[ok], ]
  cat(sprintf("starting from %s (%d site(s))\n", in_csv, sum(ok)))
} else {
  cat("starting from the true centroids\n")
}

# ---- state cell width, the unit every distance here is expressed in --------
sz <- state_size_table(mx)$per_state
cw <- setNames(sz$cs_km * 1000, sz$state_code)  # metres
cw_i <- cw[mx$state_code]

# ---- border municipios, from own shared perimeter (see header) -------------
per_km <- as.numeric(st_length(st_boundary(mx_m))) / 1000
shared_km <- rep(0, n)
for (k in seq_len(nrow(pm))) {
  shared_km[pm$i[k]] <- shared_km[pm$i[k]] + pm$border_km[k]
  shared_km[pm$j[k]] <- shared_km[pm$j[k]] + pm$border_km[k]
}
bdy_share <- 1 - pmin(shared_km / pmax(per_km, 1e-9), 1)
is_bdy <- bdy_share > BDY_TOL
cat(sprintf("state-boundary municipios: %d (%.1f%%)\n",
            sum(is_bdy), 100 * mean(is_bdy)))

# ---- which true pairs are broken in this artifact? -------------------------
hex <- readRDS(art)
hex_m <- st_transform(hex, MX_CRS)
inter <- st_is_within_distance(hex_m, hex_m, dist = 1, sparse = TRUE)
h_of <- match(mx$region, hex$region)
pa <- h_of[pm$i]
pb <- h_of[pm$j]
ok <- !is.na(pa) & !is.na(pb)
pm$surv <- NA
pm$surv[ok] <- mapply(function(x, y) y %in% inter[[x]], pa[ok], pb[ok])
cat(sprintf("artifact %s: cohesion %.2f%% (%d/%d), %d broken\n",
            basename(art), 100 * mean(pm$surv, na.rm = TRUE),
            sum(pm$surv, na.rm = TRUE), sum(!is.na(pm$surv)),
            sum(!pm$surv, na.rm = TRUE)))

# ---- false adjacency: cells touch, municipios do not, same state ----------
true_key <- paste(pmin(pm$i, pm$j), pmax(pm$i, pm$j))
m_of <- match(hex$region, mx$region)   # hex row -> municipio row
false_pairs <- do.call(rbind, lapply(seq_len(nrow(hex)), function(a) {
  bs <- inter[[a]]
  bs <- bs[bs > a]
  if (!length(bs)) return(NULL)
  ia <- m_of[a]
  ib <- m_of[bs]
  same <- mx$state_code[ia] == mx$state_code[ib]
  bs <- bs[same]
  ib <- ib[same]
  if (!length(bs)) return(NULL)
  k <- paste(pmin(ia, ib), pmax(ia, ib))
  keep <- !(k %in% true_key)
  if (!any(keep)) return(NULL)
  data.frame(i = ia, j = ib[keep])
}))
cat(sprintf("false adjacencies (cells touch, municipios do not): %d\n",
            if (is.null(false_pairs)) 0L else nrow(false_pairs)))

# ---- accumulate a unit direction per municipio ----------------------------
vx <- rep(0, n)
vy <- rep(0, n)
add_pull <- function(i, j, wgt) {
  dx <- cur[j, 1] - cur[i, 1]
  dy <- cur[j, 2] - cur[i, 2]
  d <- sqrt(dx^2 + dy^2)
  if (d < 1e-9) return(invisible(NULL))
  vx[i] <<- vx[i] + wgt * dx / d
  vy[i] <<- vy[i] + wgt * dy / d
}
brk <- pm[!is.na(pm$surv) & !pm$surv, ]
for (k in seq_len(nrow(brk))) {
  add_pull(brk$i[k], brk$j[k], 1)
  add_pull(brk$j[k], brk$i[k], 1)
}
if (!is.null(false_pairs)) {
  for (k in seq_len(nrow(false_pairs))) {
    add_pull(false_pairs$i[k], false_pairs$j[k], -REPEL)
    add_pull(false_pairs$j[k], false_pairs$i[k], -REPEL)
  }
}

# Direction only -- magnitude comes from STEP * cell width, so one municipio
# with 12 broken pairs does not take a 12x longer step than one with 1.
vnorm <- sqrt(vx^2 + vy^2)
mv <- vnorm > 1e-12
ux <- ifelse(mv, vx / pmax(vnorm, 1e-12), 0)
uy <- ifelse(mv, vy / pmax(vnorm, 1e-12), 0)

new <- cur
new[, 1] <- cur[, 1] + STEP * cw_i * ux
new[, 2] <- cur[, 2] + STEP * cw_i * uy

# ---- leash: hard cap on displacement from the TRUE centroid ---------------
leash <- ifelse(is_bdy, LEASH_BDY, LEASH_INT) * cw_i
dx <- new[, 1] - cen[, 1]
dy <- new[, 2] - cen[, 2]
dd <- sqrt(dx^2 + dy^2)
over <- dd > leash
if (any(over)) {
  sc <- leash[over] / dd[over]
  new[over, 1] <- cen[over, 1] + dx[over] * sc
  new[over, 2] <- cen[over, 2] + dy[over] * sc
}

cat(sprintf("moved %d/%d sites | %d hit the leash | mean move %.2f km | mean |P-C| %.2f km\n",
            sum(mv), n, sum(over), mean(sqrt((new[, 1] - cur[, 1])^2 +
                                             (new[, 2] - cur[, 2])^2)) / 1000,
            mean(sqrt((new[, 1] - cen[, 1])^2 +
                      (new[, 2] - cen[, 2])^2)) / 1000))

ll <- st_coordinates(st_transform(
  st_as_sf(data.frame(x = new[, 1], y = new[, 2]),
           coords = c("x", "y"), crs = MX_CRS), 4326))
write.csv(data.frame(region = mx$region, lon = ll[, 1], lat = ll[, 2]),
          out_csv, row.names = FALSE)
cat("Wrote", out_csv, "\n")
