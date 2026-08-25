suppressPackageStartupMessages({
  library(sf); library(dplyr); library(mxmaps)
})
sf_use_s2(FALSE)

# Per-state evaluation of a hex artifact: displacement (km), area evenness
# (vs equal share and vs a gamma=0.35 true-area reference), cohesion (%),
# shape regularity (IQ; regular hexagon = 0.91).
# Usage: Rscript evaluate_states.R <rds_path>

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
args <- commandArgs(trailingOnly = TRUE)
rds <- if (length(args) >= 1) args[1] else "data/mxmunicipio_hex_sf_v8_2.rds"
hex <- readRDS(rds)
cat("Evaluating:", rds, "\n\n")

data("mxmunicipio.map")
df <- mxmunicipio.map
regions <- unique(df$region)
polys <- lapply(regions, function(r) {
  sub <- df[df$region == r, ]
  groups <- unique(sub$group)
  rings <- lapply(groups, function(g) {
    pts <- sub[sub$group == g, ]
    pts <- pts[order(pts$order), ]
    coords <- cbind(pts$long, pts$lat)
    if (nrow(coords) < 4) return(NULL)
    if (!all(coords[1,] == coords[nrow(coords),])) coords <- rbind(coords, coords[1,])
    list(coords = coords, hole = any(pts$hole))
  })
  rings <- Filter(Negate(is.null), rings)
  if (length(rings) == 0) return(st_polygon())
  exterior <- lapply(rings[!sapply(rings, `[[`, "hole")], `[[`, "coords")
  holes    <- lapply(rings[sapply(rings,  `[[`, "hole")], `[[`, "coords")
  st_polygon(c(exterior, holes))
})
mxmun_sf <- st_make_valid(st_sf(region = regions,
                                geometry = st_sfc(polys, crs = 4326)))
mxmun_sf$state_code <- substr(mxmun_sf$region, 1, 2)

true_km2 <- mxmun_sf |>
  st_transform(6372) |>
  mutate(true_km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry() |>
  select(region, true_km2)

adj_all <- st_intersects(mxmun_sf, mxmun_sf)
true_pairs <- do.call(rbind, lapply(seq_along(regions), function(i) {
  js <- adj_all[[i]]
  js <- js[js > i & mxmun_sf$state_code[js] == mxmun_sf$state_code[i]]
  if (length(js) == 0) return(NULL)
  data.frame(i = i, j = js)
}))

hex_m <- st_transform(hex, 6372)
areas <- as.numeric(st_area(hex_m))
pos <- st_coordinates(st_centroid(hex_m))
tkm <- true_km2$true_km2[match(hex$region, true_km2$region)]
inter <- st_is_within_distance(hex_m, hex_m, dist = 1, sparse = TRUE)
pos_nat <- match(hex$region, mxmun_sf$region)
tcn <- st_coordinates(st_centroid(st_transform(mxmun_sf, 6372),
                                  of_largest_polygon = TRUE))

out <- lapply(sort(unique(hex$state_code)), function(s) {
  idx <- which(hex$state_code == s)
  n <- length(idx)
  share_s <- sum(areas[idx]) / n
  r_eq <- areas[idx] / share_s
  r_tg <- areas[idx] / (sum(areas[idx]) * tkm[idx]^0.35 / sum(tkm[idx]^0.35))
  dx <- pos[idx, 1] - tcn[pos_nat[idx], 1]
  dy <- pos[idx, 2] - tcn[pos_nat[idx], 2]
  disp <- sqrt(dx^2 + dy^2) / 1000
  pres <- 0L; tot <- 0L
  for (a in idx) {
    na_row <- pos_nat[a]
    js <- true_pairs$i[true_pairs$j == na_row]
    js <- c(js, true_pairs$j[true_pairs$i == na_row])
    js <- match(mxmun_sf$region[js], hex$region)
    js <- js[!is.na(js)]
    for (b in js) {
      if (b <= a) next
      tot <- tot + 1L
      if (length(inter[[a]][inter[[a]] == b]) > 0) pres <- pres + 1L
    }
  }
  iq <- vapply(idx, function(i) {
    co <- st_coordinates(hex_m$geometry[i])
    co <- co[!is.na(co[, 1]), ]
    if (nrow(co) < 3) return(NA_real_)
    x <- co[, 1]; y <- co[, 2]
    P2 <- sum(sqrt((x - c(x[-1], x[1]))^2 + (y - c(y[-1], y[1]))^2))
    4 * pi * areas[i] / (P2^2)
  }, 0)
  data.frame(state_code = s, n = n,
             disp_km = round(mean(disp)),
             eq_p5 = round(quantile(r_eq, .05), 2),
             eq_p95 = round(quantile(r_eq, .95), 2),
             tg_p5 = round(quantile(r_tg, .05), 2),
             tg_p95 = round(quantile(r_tg, .95), 2),
             cohesion = round(100 * pres / max(tot, 1)),
             iq = round(mean(iq, na.rm = TRUE), 2))
})
out_df <- do.call(rbind, out)
print(out_df, row.names = FALSE)
cat("\nNational cohesion (weighted by pairs):",
    round(sum(out_df$cohesion * out_df$n) / sum(out_df$n), 1), "%\n")
