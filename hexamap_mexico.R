library(mxmaps)
library(ggplot2)
library(dplyr)
library(scales)

# --- map geometry (fortified data.frame: long, lat, group, region) -----------
data("mxhexbin.map")

# --- state labels (centroid coords baked into the package) -------------------
state_labels <- data.frame(
  region     = c("03","02","25","26","18","10","12","16","06","14","32","08",
                 "17","15","11","01","05","20","21","09","22","24","19","07",
                 "29","13","30","28","27","04","23","31"),
  state_abbr = c("BCS","BC","SIN","SON","NAY","DGO","GRO","MICH","COL","JAL",
                 "ZAC","CHIH","MOR","MEX","GTO","AGS","COAH","OAX","PUE",
                 "CDMX","QRO","SLP","NL","CHPS","TLAX","HGO","VER","TAM",
                 "TAB","CAMP","QROO","YUC"),
  long = c(-107.424,-107.424,-104.328,-104.328,-101.232,-101.232,-98.136,
           -98.136,-98.136,-98.136,-98.136,-98.136,-95.040,-95.040,-95.040,
           -95.040,-95.040,-91.944,-91.944,-91.944,-91.944,-91.944,-91.944,
           -88.848,-88.848,-88.848,-88.848,-88.848,-85.752,-82.656,
           -79.560,-79.560),
  lat  = c(11.387,14.962,13.175,16.750,11.387,14.962,-1.125,2.450,6.025,
            9.600,13.175,16.750,0.662,4.237,7.812,11.387,14.962,-4.700,
           -1.125,2.450,6.025,9.600,13.175,-2.913,0.662,4.237,7.812,
            11.387,-1.125,0.662,-1.125,2.450)
)

# --- data --------------------------------------------------------------------
data("df_mxstate_2020")

plot_df <- df_mxstate_2020 |>
  mutate(value = indigenous_language / pop)   # share speaking indigenous language

# join value onto hex geometry by the INEGI region code
map_data <- mxhexbin.map |>
  left_join(plot_df |> select(region, value), by = "region")

# --- plot (mirrors what MXHexBinChoropleth$render() produces) ----------------
ggplot(map_data, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = value), color = "dark grey", linewidth = 0.2) +
  geom_text(
    data     = state_labels,
    aes(long, lat, label = state_abbr, group = NULL),
    color    = "black",
    size     = 4.5
  ) +
  scale_fill_continuous(
    "Share",
    low      = "#eff3ff",
    high     = "#084594",
    na.value = "black",
    labels   = percent_format(accuracy = 1)
  ) +
  coord_quickmap() +
  theme_void() +
  ggtitle("Share of population speaking an indigenous language (2020)")
