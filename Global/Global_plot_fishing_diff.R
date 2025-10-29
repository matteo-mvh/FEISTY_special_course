# ----------------------------------------
# Plot differences between No Fishing and other fishing scenarios
# ----------------------------------------

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)

# Define fishing scenarios
fishing_scenarios <- c("No_Fishing", "Demersal", "Forage_Fish", "Large_Pelagic")

# Load No Fishing baseline
load(paste0("data/Global_fish_biomass_", fishing_scenarios[1], ".RData"))
baseline <- out
baseline$totB_all <- as.numeric(baseline$totB_all)
baseline$carbon_inject <- as.numeric(baseline$carbon_inject)

# Get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Function to plot global difference
plot_difference <- function(diff_sf, value_col, title) {
  ggplot() +
    geom_sf(data = world, fill = "gray90", color = "gray50") +
    geom_sf(data = diff_sf, aes(color = !!sym(value_col)), size = 1) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    coord_sf(crs = "+proj=robin") +
    theme_minimal() +
    labs(color = title, title = title)
}

# Loop over all fishing scenarios except No Fishing
for (scenario in fishing_scenarios[-1]) {
  
  # Load scenario data
  load(paste0("data/Global_fish_biomass_", scenario, ".RData"))
  scenario_df <- out
  scenario_df$totB_all <- as.numeric(scenario_df$totB_all)
  scenario_df$carbon_inject <- as.numeric(scenario_df$carbon_inject)
  
  # Calculate differences
  diff_df <- scenario_df
  diff_df$totB_all <- scenario_df$totB_all - baseline$totB_all
  diff_df$carbon_inject <- scenario_df$carbon_inject - baseline$carbon_inject
  
  # Convert to sf
  diff_sf <- st_as_sf(diff_df, coords = c("lon", "lat"), crs = 4326)
  
  # Plot differences
  p1 <- plot_difference(diff_sf, "totB_all", paste0("Change in Total Biomass: ", scenario, " vs No Fishing"))
  p2 <- plot_difference(diff_sf, "carbon_inject", paste0("Change in Carbon Injection: ", scenario, " vs No Fishing"))
  
  print(p1)
  print(p2)
}
