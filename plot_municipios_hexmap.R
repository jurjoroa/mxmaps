library(sf)
library(ggplot2)
library(dplyr)
library(scales)
library(mxmaps)

# Load precomputed hex layout (run hexamap_municipios.R once to generate)
hex_sf <- readRDS("data/mxmunicipio_hex_sf.rds")

data("df_mxmunicipio_2020")

census <- df_mxmunicipio_2020 |>
  mutate(
    indigenous_share = indigenous_language / pop,
    afromexican_share = afromexican / pop
  )

hex_data <- hex_sf |>
  left_join(census, by = "region")

# helper so each map is one function call
mx_hex_map <- function(data, fill_var, title, subtitle = "Mexico — 2020 census, by municipio",
                       palette = "YlGnBu", label_fmt = percent_format(accuracy = 1)) {
  ggplot(data) +
    geom_sf(aes(fill = .data[[fill_var]]), color = "dark grey", linewidth = 0.05) +
    scale_fill_distiller(palette = palette, direction = 1,
                         labels = label_fmt, na.value = "grey80") +
    coord_sf() +
    theme_void(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", margin = margin(b = 4)),
      plot.subtitle = element_text(color = "grey40", margin = margin(b = 8)),
      legend.title  = element_blank(),
      legend.key.height = unit(1.5, "cm")
    ) +
    labs(title = title, subtitle = subtitle)
}

# --- Map 1: Total population -------------------------------------------------
mx_hex_map(hex_data, "pop",
           title   = "Total population",
           palette = "Blues",
           label_fmt = label_comma())

# --- Map 2: Indigenous-language speakers -------------------------------------
mx_hex_map(hex_data, "indigenous_share",
           title   = "Share speaking an indigenous language",
           palette = "YlGnBu")

# --- Map 3: Afro-Mexican population ------------------------------------------
mx_hex_map(hex_data, "afromexican_share",
           title   = "Share self-identifying as Afro-Mexican",
           palette = "YlOrRd")

# --- Map 4: Female share of population ---------------------------------------
hex_data <- hex_data |> mutate(female_share = pop_female / pop)

mx_hex_map(hex_data, "female_share",
           title   = "Female share of population",
           palette = "RdPu")
