library(sf)
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)
library(mxmaps)

sf_use_s2(FALSE)

# Load both layouts for side-by-side comparison
hex_v1 <- readRDS("data/mxmunicipio_hex_sf.rds")
hex_v2 <- readRDS("data/mxmunicipio_hex_sf_v2.rds")

data("df_mxmunicipio_2020")

census <- df_mxmunicipio_2020 |>
  mutate(
    indigenous_share  = indigenous_language / pop,
    afromexican_share = afromexican / pop,
    female_share      = pop_female / pop
  )

mx_hex_map <- function(layout, fill_var, title, palette = "YlGnBu",
                      label_fmt = percent_format(accuracy = 1)) {
  ggplot(layout |> left_join(census, by = "region")) +
    geom_sf(aes(fill = .data[[fill_var]]),
            color = "white", linewidth = 0.05) +
    scale_fill_distiller(palette = palette, direction = 1,
                         labels = label_fmt, na.value = "grey80") +
    coord_sf() +
    theme_void(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey40", size = 10),
      legend.title  = element_blank(),
      legend.key.height = unit(0.8, "cm")
    ) +
    labs(title = title)
}

# --- Side-by-side v1 vs v2 (indigenous-language share) ---------------------
p1 <- mx_hex_map(hex_v1, "indigenous_share",
                 title = "v1 — per-state geogrid (overlapping hexes)")
p2 <- mx_hex_map(hex_v2, "indigenous_share",
                 title = "v2 — Mexico-wide hex lattice (no overlaps)")

comparison <- p1 + p2 +
  plot_annotation(
    title    = "Indigenous-language share, by municipio (2020)",
    subtitle = "v1 has between-state hex overlaps; v2 snaps every municipio to a unique cell on a shared lattice",
    theme    = theme(plot.title    = element_text(face = "bold"),
                     plot.subtitle = element_text(color = "grey40"))
  )

ggsave("v1_vs_v2_comparison.png", comparison,
       width = 14, height = 6, dpi = 150)
cat("Saved v1_vs_v2_comparison.png\n")

# --- v2 panel of four indicators -------------------------------------------
quad <- (mx_hex_map(hex_v2, "pop",
                   title = "Total population", palette = "Blues",
                   label_fmt = label_comma()) +
         mx_hex_map(hex_v2, "indigenous_share",
                   title = "Indigenous-language share")) /
        (mx_hex_map(hex_v2, "afromexican_share",
                   title = "Afro-Mexican share", palette = "YlOrRd") +
         mx_hex_map(hex_v2, "female_share",
                   title = "Female share", palette = "RdPu"))

ggsave("v2_indicators_quad.png", quad,
       width = 14, height = 12, dpi = 150)
cat("Saved v2_indicators_quad.png\n")
