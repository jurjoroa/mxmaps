suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales); library(patchwork)
  library(mxmaps)
})
sf_use_s2(FALSE)

# v3 plotting template: hex cells filled with data, state borders colored
# distinctly so adjacent states are easy to tell apart.

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

hex     <- readRDS("data/mxmunicipio_hex_sf_v3.rds")
borders <- readRDS("data/mxstate_borders.rds")

data("df_mxmunicipio_2020")
census <- df_mxmunicipio_2020 |>
  select(region, indigenous_language, pop, afromexican, pop_female) |>
  mutate(
    indigenous_share  = indigenous_language / pop,
    afromexican_share = afromexican / pop,
    female_share      = pop_female / pop
  )
hex_data <- hex |> left_join(census, by = "region")

# State labels: human-readable names from df_mxmunicipio_2020 (one row per state)
state_lookup <- df_mxmunicipio_2020 |>
  select(state_code, state_name, state_abbr) |>
  distinct()
borders <- borders |> left_join(state_lookup, by = "state_code")
suppressWarnings({
  borders_centroids <- st_centroid(borders) |>
    mutate(x = st_coordinates(geometry)[, 1],
           y = st_coordinates(geometry)[, 2])
})

# 32-color qualitative palette — high-contrast, no near-duplicates between
# typical neighbors. Built by hand from Polychrome / Set3 / Paired families.
state_colors <- c(
  "01" = "#E41A1C", "02" = "#377EB8", "03" = "#4DAF4A", "04" = "#984EA3",
  "05" = "#FF7F00", "06" = "#FFD92F", "07" = "#A65628", "08" = "#F781BF",
  "09" = "#1B9E77", "10" = "#7570B3", "11" = "#E7298A", "12" = "#66A61E",
  "13" = "#E6AB02", "14" = "#A6761D", "15" = "#666666", "16" = "#1F78B4",
  "17" = "#33A02C", "18" = "#FB9A99", "19" = "#FDBF6F", "20" = "#CAB2D6",
  "21" = "#B15928", "22" = "#8DD3C7", "23" = "#80B1D3", "24" = "#FCCDE5",
  "25" = "#BC80BD", "26" = "#CCEBC5", "27" = "#FFED6F", "28" = "#D95F02",
  "29" = "#7298c0", "30" = "#dd5050", "31" = "#5b9d6f", "32" = "#9467bd"
)

base_theme <- theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15, margin = margin(b = 4)),
        plot.subtitle = element_text(color = "grey40", size = 10, margin = margin(b = 8)),
        legend.title = element_blank(),
        legend.key.height = unit(1.5, "cm"))

# ---- Plot 1: data + colored state borders ---------------------------------
p1 <- ggplot() +
  geom_sf(data = hex_data, aes(fill = indigenous_share),
          color = "grey60", linewidth = 0.04) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1,
                       labels = percent_format(accuracy = 1),
                       na.value = "grey90") +
  geom_sf(data = borders, aes(color = state_code),
          fill = NA, linewidth = 0.85) +
  scale_color_manual(values = state_colors, guide = "none") +
  geom_text(data = borders_centroids, aes(x = x, y = y, label = state_abbr),
            size = 2.6, color = "black", fontface = "bold") +
  coord_sf() + base_theme +
  labs(title = "Share speaking an indigenous language",
       subtitle = "v3 hex grid — each state outlined in its own color")

# ---- Plot 2: cell edges colored by state (alternative styling) ------------
hex_data$state_color <- state_colors[hex_data$state_code]
p2 <- ggplot() +
  geom_sf(data = hex_data,
          aes(fill = indigenous_share, color = state_code),
          linewidth = 0.18) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1,
                       labels = percent_format(accuracy = 1),
                       na.value = "grey90") +
  scale_color_manual(values = state_colors, guide = "none") +
  geom_sf(data = borders, fill = NA, color = "black", linewidth = 0.35) +
  coord_sf() + base_theme +
  labs(title = "Same data, cell edges colored by state",
       subtitle = "Every hex's outline matches its state color — overlap with thin black state border")

ggsave("v3_colored_borders.png", p1,
       width = 12, height = 9, dpi = 150, bg = "white")
ggsave("v3_colored_borders_alt.png", p2,
       width = 12, height = 9, dpi = 150, bg = "white")
ggsave("v3_colored_borders_compare.png", p1 / p2,
       width = 13, height = 17, dpi = 130, bg = "white")

cat("Wrote: v3_colored_borders.png, v3_colored_borders_alt.png, v3_colored_borders_compare.png\n")
