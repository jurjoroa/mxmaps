library(sf)
library(geogrid)
library(ggplot2)
library(dplyr)
library(scales)
library(mxmaps)

# --- 1. Reconstruct sf from the fortified mxmunicipio.map -------------------
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
mxmun_sf <- st_sf(region = regions, geometry = st_sfc(polys, crs = 4326))
mxmun_sf$state_code <- substr(mxmun_sf$region, 1, 2)

# --- 2. Hex grid assignment (per state, chunked for large states) ------------
# Large states (Oaxaca=570, Puebla=217, Veracruz=212) are split into
# geographic k-means chunks of ≤60 each — keeps each assignment problem
# tractable (O(n³) → ~100x faster than solving all at once).

run_geogrid <- function(sub) {
  grid <- calculate_grid(sub, learning_rate = 0.03, grid_type = "hexagonal", seed = 1)
  res  <- assign_polygons(sub, grid)
  res[, c("region", "geometry")]
}

geogrid_chunked <- function(sub, max_n = 60) {
  n <- nrow(sub)
  if (n <= max_n) return(run_geogrid(sub))
  cents <- st_coordinates(st_centroid(sub))
  k <- ceiling(n / max_n)
  sub$chunk <- kmeans(cents, centers = k, nstart = 5)$cluster
  chunks <- lapply(unique(sub$chunk), function(ck) run_geogrid(sub[sub$chunk == ck, ]))
  do.call(rbind, chunks)
}

cat("Computing hex layout (~60 seconds)...\n")
states  <- unique(mxmun_sf$state_code)
results <- lapply(states, function(s) {
  sub <- mxmun_sf[mxmun_sf$state_code == s, ]
  cat("State", s, "- n:", nrow(sub), "\n")
  geogrid_chunked(sub, max_n = 60)
})
hex_sf <- do.call(rbind, results)

# Save layout so future runs skip the 53-second compute
saveRDS(hex_sf, "data/mxmunicipio_hex_sf.rds")
# hex_sf <- readRDS("data/mxmunicipio_hex_sf.rds")

# --- 3. Join census data and plot -------------------------------------------
data("df_mxmunicipio_2020")

hex_data <- hex_sf |>
  left_join(
    df_mxmunicipio_2020 |>
      select(region, state_abbr, indigenous_language, pop, afromexican) |>
      mutate(value = indigenous_language / pop),
    by = "region"
  )

ggplot(hex_data) +
  geom_sf(aes(fill = value), color = "dark grey", linewidth = 0.05) +
  scale_fill_distiller(
    "Share",
    palette   = "YlGnBu",
    direction = 1,
    labels    = percent_format(accuracy = 1),
    na.value  = "grey80"
  ) +
  coord_sf() +
  theme_void() +
  ggtitle(
    "Share of population speaking an indigenous language",
    subtitle = "Mexico — 2020 census, by municipio (hexagonal cartogram)"
  )

# --- swap variable: Afro-Mexican share ---------------------------------------
# hex_data$value <- hex_data$afromexican / hex_data$pop
# then re-run the ggplot block above
