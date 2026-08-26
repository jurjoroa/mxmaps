# Verify the two properties v9 claims by construction:
#   1. state outline preserved EXACTLY (union of cells == true state polygon)
#   2. adjacency preserved EXACTLY (every true intra-state neighbour pair still
#      touches) -- plus geometry validity and absence of overlaps.
# Usage: Rscript verify_v9.R <rds_path>
suppressPackageStartupMessages({ library(sf); library(dplyr); library(mxmaps) })
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
args <- commandArgs(trailingOnly = TRUE)
rds <- if (length(args) >= 1) args[1] else "data/mxmunicipio_equalarea_v9.rds"
hex <- readRDS(rds)
cat("Verifying:", rds, "\n\n")

data("mxmunicipio.map")
df <- mxmunicipio.map
regions <- unique(df$region)
polys <- lapply(regions, function(r) {
  sub <- df[df$region == r, ]; groups <- unique(sub$group)
  rings <- lapply(groups, function(g) {
    pts <- sub[sub$group == g, ]; pts <- pts[order(pts$order), ]
    coords <- cbind(pts$long, pts$lat)
    if (nrow(coords) < 4) return(NULL)
    if (!all(coords[1,] == coords[nrow(coords),])) coords <- rbind(coords, coords[1,])
    list(coords = coords, hole = any(pts$hole))
  })
  rings <- Filter(Negate(is.null), rings)
  if (length(rings) == 0) return(st_polygon())
  st_polygon(c(lapply(rings[!sapply(rings,`[[`,"hole")],`[[`,"coords"),
               lapply(rings[ sapply(rings,`[[`,"hole")],`[[`,"coords")))
})
mx <- st_make_valid(st_sf(region = regions, geometry = st_sfc(polys, crs = 4326)))
mx$state_code <- substr(mx$region, 1, 2)

adj_all <- st_intersects(mx, mx)

for (s in sort(unique(hex$state_code))) {
  h <- st_transform(hex[hex$state_code == s, ], 6372)
  t <- st_transform(mx[mx$state_code == s, ], 6372)

  # --- 1. outline
  uh <- st_union(st_geometry(h)); ut <- st_union(st_geometry(t))
  sym <- suppressWarnings(st_sym_difference(uh, ut))
  lost <- if (length(sym) == 0) 0 else sum(as.numeric(st_area(sym)))
  outline_err <- 100 * lost / as.numeric(st_area(ut))

  # --- 2. adjacency
  ti <- match(t$region, mx$region)
  tp <- do.call(rbind, lapply(seq_len(nrow(t)), function(i) {
    js <- adj_all[[ti[i]]]
    js <- js[js != ti[i] & mx$state_code[js] == s]
    if (!length(js)) return(NULL)
    data.frame(a = i, b = match(mx$region[js], t$region))
  }))
  tp <- tp[!is.na(tp$b) & tp$a < tp$b, ]
  hi <- match(t$region, h$region)
  touch <- st_is_within_distance(h, h, dist = 1, sparse = TRUE)
  kept <- sum(vapply(seq_len(nrow(tp)), function(k) {
    a <- hi[tp$a[k]]; b <- hi[tp$b[k]]
    as.integer(!is.na(a) && !is.na(b) && b %in% touch[[a]])
  }, 0L))

  # --- 3. validity + overlaps
  n_invalid <- sum(!st_is_valid(h))
  ov <- st_intersection(st_make_valid(h), st_make_valid(h))
  ov <- ov[ov$region != ov$region.1, ]
  ov_area <- if (nrow(ov) == 0) 0 else
    sum(as.numeric(st_area(ov[st_dimension(ov) == 2, ]))) / 2
  ov_pct <- 100 * ov_area / as.numeric(st_area(ut))

  cat(sprintf(
    "state %s n=%3d | outline err %6.4f%% | adjacency %4d/%4d = %5.1f%% | invalid %d | overlap %6.4f%%\n",
    s, nrow(h), outline_err, kept, nrow(tp), 100 * kept / max(nrow(tp), 1),
    n_invalid, ov_pct))
}
