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

state_km <- true_km2 |>
  dplyr::mutate(state_code = substr(region, 1, 2)) |>
  dplyr::group_by(state_code) |>
  dplyr::summarise(state_km2 = sum(true_km2), n = dplyr::n(), .groups = "drop") |>
  dplyr::mutate(cs_km = sqrt(state_km2 / n / 0.866))
cs_nat <- sqrt(sum(state_km$state_km2) / sum(state_km$n) / 0.866)
gam_tab <- state_km |>
  dplyr::mutate(gamma_auto = pmax(pmin(0.18 * (cs_km / cs_nat - 0.75), 0.45), 0))

out <- lapply(sort(unique(hex$state_code)), function(s) {
  idx <- which(hex$state_code == s)
  n <- length(idx)
  share_s <- sum(areas[idx]) / n
  r_eq <- areas[idx] / share_s
  # The area-vs-target reference used to hardcode gamma = 0.35, but the
  # generators set gamma PER STATE (the derived rule plus, historically, 12
  # hand overrides), so for 31 of 32 states this column was scoring against a
  # gamma the layout never used. Derive the same per-state rule here. This
  # still cannot see the 12 hand overrides -- the real fix is for generators to
  # ship their targets alongside the artifact.
  gam_s <- gam_tab$gamma_auto[match(s, gam_tab$state_code)]
  if (is.na(gam_s)) gam_s <- 0.35
  r_tg <- areas[idx] / (sum(areas[idx]) * tkm[idx]^gam_s / sum(tkm[idx]^gam_s))
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
             eq_min = round(min(r_eq), 3),
             n_tiny = sum(r_eq < 0.25),
             pairs = tot,
             cohesion = round(100 * pres / max(tot, 1)),
             pres_n = pres,
             iq = round(mean(iq, na.rm = TRUE), 2))
})
out_df <- do.call(rbind, out)
print(out_df[, setdiff(names(out_df), "pres_n")], row.names = FALSE)

# National cohesion, two ways. The pair-weighted figure is the real one: it is
# the fraction of ALL true intra-state municipio adjacencies that survive as
# touching cells. The cell-weighted figure is retained only for comparability
# with the v8.3-v8.5 commit messages, which quoted it while labelling it
# "weighted by pairs" -- dense states carry more pairs per cell, so it reads
# optimistically.
cat("\nNational cohesion (pair-weighted, exact):",
    round(100 * sum(out_df$pres_n) / sum(out_df$pairs), 1), "%\n")
cat("National cohesion (cell-weighted, legacy):",
    round(sum(out_df$cohesion * out_df$n) / sum(out_df$n), 1), "%\n")

# LEGIBILITY TAIL. The p5/p95 area-ratio columns cannot see a handful of
# near-invisible cells, and a version can buy cohesion by creating them --
# which is exactly what v8.8 did (38 cells under 0.25x their state mean against
# 2 for v8.6). A municipio rendered at 1/10 of its state's mean cell area is
# not legible, and legibility is the entire point of the layout, so this is a
# first-class number and not a footnote.
cat(sprintf("Legibility tail: %d cells < 0.25x their state mean, %d < 0.10x, smallest %.4fx\n",
            sum(out_df$n_tiny),
            sum(vapply(sort(unique(hex$state_code)), function(s) {
              idx <- which(hex$state_code == s)
              a <- areas[idx]; sum(a / mean(a) < 0.10)
            }, 0L)),
            min(out_df$eq_min)))
