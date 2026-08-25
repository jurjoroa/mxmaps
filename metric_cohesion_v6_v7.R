suppressPackageStartupMessages({
  library(sf); library(dplyr); library(mxmaps)
})
sf_use_s2(FALSE)

# Territorial cohesion metric: share of true-adjacent municipio pairs (queen
# contiguity, same state) whose hex cells are also adjacent. v6 vs v7.

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

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

adj_all <- st_intersects(mxmun_sf, mxmun_sf)
true_pairs <- do.call(rbind, lapply(seq_along(regions), function(i) {
  js <- adj_all[[i]]
  js <- js[js > i & mxmun_sf$state_code[js] == mxmun_sf$state_code[i]]
  if (length(js) == 0) return(NULL)
  data.frame(i = i, j = js)
}))
cat(sprintf("True same-state adjacent pairs: %d\n", nrow(true_pairs)))

cohesion <- function(hex, label) {
  hex_m <- st_transform(hex, 6372)
  inter <- st_is_within_distance(hex_m, hex_m, dist = 1e-6, sparse = TRUE)
  pos <- match(hex$region, mxmun_sf$region)
  ok <- pres <- 0L
  for (k in seq_len(nrow(true_pairs))) {
    a <- pos[true_pairs$i[k]]; b <- pos[true_pairs$j[k]]
    if (is.na(a) || is.na(b)) next
    ok <- ok + 1L
    if (length(inter[[a]][inter[[a]] == b]) > 0) pres <- pres + 1L
  }
  cat(sprintf("%s cohesion: %d/%d = %.1f%%\n", label, pres, ok, 100 * pres / ok))
  per_state <- tapply(seq_len(ok), rep(1, ok), length)  # placeholder
  invisible(pres / ok)
}
v6 <- readRDS("data/mxmunicipio_hex_sf_v6.rds")
v7 <- readRDS("data/mxmunicipio_hex_sf_v7.rds")
cohesion(v6, "v6")
cohesion(v7, "v7")
