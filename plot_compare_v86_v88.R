# ---------------------------------------------------------------------------
# Purpose: three-panel comparison so the anchor-floor change is visible:
#          v8.6 (anchor 0 for 79% of municipios) vs v8.8 (anchor floored at
#          0.30) vs true territory, coloured by inhabitants per 100k.
# Author:  Claude (Opus 5) for Jorge Roa
# Version: v8.8
# Date:    26-August-2026
# Usage:   Rscript plot_compare_v86_v88.R [state_code ...]
# Output:  state_pages/cmp88_<code>_<abbr>.png + cmp88_national.png
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales)
  library(patchwork); library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
dir.create("state_pages", showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
ONLY <- if (length(args)) sprintf("%02d", as.integer(args)) else NULL

a86 <- readRDS("data/mxmunicipio_hex_sf_v8_6.rds")
a88 <- readRDS("data/mxmunicipio_hex_sf_v8_8.rds")
borders <- readRDS("data/mxstate_borders.rds")

data("mxmunicipio.map"); df <- mxmunicipio.map
rg <- unique(df$region)
pl <- lapply(rg, function(r){s<-df[df$region==r,]; g<-unique(s$group)
  rr<-lapply(g,function(gg){p<-s[s$group==gg,];p<-p[order(p$order),];co<-cbind(p$long,p$lat)
    if(nrow(co)<4) return(NULL); if(!all(co[1,]==co[nrow(co),])) co<-rbind(co,co[1,])
    list(coords=co,hole=any(p$hole))}); rr<-Filter(Negate(is.null),rr)
  if(!length(rr)) return(st_polygon())
  st_polygon(c(lapply(rr[!sapply(rr,`[[`,"hole")],`[[`,"coords"),
               lapply(rr[ sapply(rr,`[[`,"hole")],`[[`,"coords")))})
real <- st_make_valid(st_sf(region=rg, geometry=st_sfc(pl, crs=4326)))
real$state_code <- substr(real$region,1,2)
state_sf <- real |> group_by(state_code) |>
  summarise(geometry=st_union(geometry), .groups="drop") |> st_make_valid()

data("df_mxmunicipio_2020")
cs <- df_mxmunicipio_2020 |> select(region, state_code, pop)
natp <- sum(cs$pop)
sp <- cs |> group_by(state_code) |> summarise(spop=sum(pop), .groups="drop")
cs <- cs |> left_join(sp, by="state_code") |>
  mutate(p_state = pop/spop*1e5, p_nat = pop/natp*1e5) |>
  select(region, p_state, p_nat)
lk <- df_mxmunicipio_2020 |> select(state_code, state_name, state_abbr) |> distinct()

a86 <- a86 |> left_join(cs, by="region")
a88 <- a88 |> left_join(cs, by="region")
realj <- real |> left_join(cs, by="region")

pan <- function(d, xlim, ylim, title, fv, is_real) {
  p <- ggplot() +
    scale_fill_distiller(palette="YlGnBu", direction=1, name="Inhabitants\nper 100k",
                         labels=comma_format(accuracy=1), na.value="grey90")
  if (is_real) p <- p + geom_sf(data=d, aes(fill=.data[[fv]]), color="grey30", linewidth=0.12)
  else p <- p + geom_sf(data=d, aes(fill=.data[[fv]]), color="grey40", linewidth=0.08) +
    geom_sf(data=borders, fill=NA, color="black", linewidth=0.45)
  p + coord_sf(xlim=xlim, ylim=ylim, expand=FALSE) + theme_void(base_size=11) +
    theme(plot.title=element_text(face="bold", size=12), legend.position="right") +
    labs(title=title)
}

if (is.null(ONLY)) {
  xl <- c(-118.5,-86.5); yl <- c(14,32.8)
  p <- pan(a86, xl, yl, "v8.6 — anchor 0 for 79% of municipios", "p_nat", FALSE) +
       pan(a88, xl, yl, "v8.8 — anchor floored at 0.30", "p_nat", FALSE) +
       pan(realj, xl, yl, "True municipio territory", "p_nat", TRUE) +
       plot_layout(guides="collect", nrow=1) +
    plot_annotation(
      title="Mexico — effect of flooring the position anchor, inhabitants per 100k",
      subtitle="v8.6 cohesion 73.2%, mean displacement 41 km  |  v8.8 cohesion 78.1%, mean displacement 22 km",
      theme=theme(plot.title=element_text(face="bold", size=17),
                  plot.subtitle=element_text(color="grey40", size=11)))
  ggsave("state_pages/cmp88_national.png", p, width=22, height=7, dpi=115, bg="white")
  cat("Wrote state_pages/cmp88_national.png\n")
}

codes <- if (is.null(ONLY)) sort(unique(real$state_code)) else ONLY
for (s in codes) {
  b <- st_bbox(state_sf$geometry[state_sf$state_code==s])
  px <- 0.15*(b["xmax"]-b["xmin"]); py <- 0.15*(b["ymax"]-b["ymin"])
  xl <- c(b["xmin"]-px, b["xmax"]+px); yl <- c(b["ymin"]-py, b["ymax"]+py)
  nm <- lk$state_name[lk$state_code==s][1]; ab <- lk$state_abbr[lk$state_code==s][1]
  p <- pan(filter(a86, state_code==s), xl, yl, "v8.6 — unanchored", "p_state", FALSE) +
       pan(filter(a88, state_code==s), xl, yl, "v8.8 — anchored", "p_state", FALSE) +
       pan(filter(realj, state_code==s), xl, yl, "True territory", "p_state", TRUE) +
       plot_layout(guides="collect", nrow=1) +
    plot_annotation(title=sprintf("%s (%s) — inhabitants per 100k, anchor effect", nm, ab),
      subtitle="Same cells and same areas; only the pull of each cell toward its municipio's true position differs",
      theme=theme(plot.title=element_text(face="bold", size=16),
                  plot.subtitle=element_text(color="grey40", size=10)))
  f <- sprintf("state_pages/cmp88_%s_%s.png", s, ab)
  ggsave(f, p, width=20, height=7.5, dpi=115, bg="white")
  cat("Wrote", f, "\n")
}
