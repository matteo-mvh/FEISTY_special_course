# ----------------------------------------
# Global plots of FEISTY output (Robinson projection)
# ----------------------------------------

# Load libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)

# Load the saved FEISTY output
load("data/Global_fish_biomass_No_Fishing.RData")  # loads 'out'

# Make sure numeric columns are numeric
out$lon <- as.numeric(out$lon)
out$lat <- as.numeric(out$lat)
out$totB_all <- as.numeric(out$totB_all)
out$carbon_inject <- as.numeric(out$carbon_inject)

# Convert to sf object
points_sf <- st_as_sf(out, coords = c("lon", "lat"), crs = 4326)

# Get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Function to plot global points with Robinson projection
plot_global <- function(data_sf, value_col, title, palette="viridis") {
  ggplot() +
    geom_sf(data = world, fill = "gray90", color = "gray50") +
    geom_sf(data = data_sf, aes(color = !!sym(value_col)), size = 1) +
    scale_color_viridis_c(option = "viridis") +
    coord_sf(crs = "+proj=robin") +
    theme_minimal() +
    labs(color = title, title = title)
}

# Plot total biomass
plot_global(points_sf, "totB_all", "Total Fish Biomass (g/m²)")

# Plot carbon injection
plot_global(points_sf, "carbon_inject", "Carbon Injection (g C/m²)")

