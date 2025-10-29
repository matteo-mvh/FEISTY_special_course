# -----------------------------------------------------------------------------
# Description:
#   Plot biomass changes across Shelf Sea, Slope, and Open Ocean
# -----------------------------------------------------------------------------

# Set results directory
results_dir <- "C:/Users/Mmm/OneDrive/Master Studies/3. Semester/Carbon Sequesteration/Code/results_files"

# Load necessary library and data
library(FEISTY)
library(gridExtra)
library(data.table)
library(pracma)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)

load(file.path(results_dir, "FEISTY_Simulations_Results.RData"))
load(file.path(results_dir, "FEISTY_Biomass_Results.RData"))
load(file.path(results_dir, "FEISTY_Flux_Results.RData"))
load(file.path(results_dir, "FEISTY_Injection_Results.RData"))

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

# Difference (or percent if you divide by nofish and *100)
perc_change <- function(fish, nofish) {
  fish - nofish
}

biomass_names <- c("smallPel", "mesoPel", "largePel", "midwPred", "demersal")

# consistent colors by scenario
scenario_colors <- c(
  "Demersal"    = "green",
  "Forage Fish" = "orange",
  "Large Pelagic" = "red"
)

fish_colors <- c(
  "smallPel"   = "lightcoral",
  "mesoPel"    = "purple",
  "largePel"   = "lightblue",
  "midwPred"   = "blue",
  "demersal"   = "darkgreen"
)

plot_changes <- function(location) {
  location_data <- avg_Biomass[[location]]
  no_fishing <- location_data[["No Fishing"]]
  names(no_fishing) <- biomass_names
  
  # find scenarios available (besides "No Fishing")
  scenarios <- setdiff(names(location_data), "No Fishing")
  
  # compute changes
  changes <- lapply(scenarios, function(sc) {
    vals <- location_data[[sc]]
    names(vals) <- biomass_names
    perc_change(vals, no_fishing)
  })
  names(changes) <- scenarios
  
  # stack into matrix
  change_mat <- do.call(rbind, changes)
  
  # get colors for the available scenarios
  scenario_cols <- scenario_colors[scenarios]
  
  # plot
  barplot(change_mat,
          beside = TRUE,
          col = scenario_cols,
          names.arg = biomass_names,
          main = paste("Biomass Change by Fishing Scenario for", location),
          ylab = "Change",
          xlab = "Biomass Type")
  
  legend("topright",
         legend = scenarios,
         fill = scenario_cols)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Function: Compare No Fishing scenarios across all locations
# -----------------------------------------------------------------------------
plot_no_fishing <- function() {
  # extract "No Fishing" biomass from each location
  nf_data <- lapply(names(avg_Biomass), function(loc) {
    vals <- avg_Biomass[[loc]][["No Fishing"]]
    names(vals) <- biomass_names
    vals
  })
  names(nf_data) <- names(avg_Biomass)
  
  # combine into matrix (rows = locations, cols = biomass types)
  nf_mat <- do.call(rbind, nf_data)
  
  # transpose so cols = locations, rows = fish types
  nf_mat_t <- t(nf_mat)
  
  # determine y-axis max
  ymax <- ceiling(max(nf_mat_t)) * 1.1
  
  
  
  # barplot
  bp <- barplot(nf_mat_t,
                beside = TRUE,
                col = fish_colors[rownames(nf_mat_t)],
                names.arg = rownames(nf_mat),   # x-axis = locations
                ylim = c(0, ymax),
                ylab = expression("Biomass in [g" ~ m^-2 * "]"),
                border = "black")
  
  # add grid lines
  abline(h = pretty(c(0, ymax)), col = "#00000020")
  
  # calculate vertical lines between groups
  n_fish <- nrow(nf_mat_t)   # number of fish types per group
  n_loc  <- ncol(nf_mat_t)   # number of locations
  for (i in 1:(n_loc-1)) {
    # last bar of group i and first bar of group i+1
    x_sep <- (bp[i * n_fish] + bp[i * n_fish + 1]) / 2
    abline(v = x_sep, col = "black", lty = "solid")
  }
  
  # add box around plot
  box(bty = "o")
  
  # legend
  legend("topright",
         legend = names(fish_colors),
         fill = fish_colors,
         bty = "n")
}

# EXTRA PLOT: Biomass changes across all locations and scenarios

plot_biomass_changes <- function() {
  # Only keep Shelf Sea and Slope
  locations <- c("Shelf Sea", "Slope")
  
  # Collect data into long format
  change_df <- do.call(rbind, lapply(locations, function(loc) {
    loc_data <- avg_Biomass[[loc]]
    nofish <- loc_data[["No Fishing"]]
    names(nofish) <- biomass_names
    
    # scenarios = all except No Fishing
    scenarios <- setdiff(names(loc_data), "No Fishing")
    
    # compute differences
    do.call(rbind, lapply(scenarios, function(sc) {
      vals <- loc_data[[sc]]
      names(vals) <- biomass_names
      diff <- vals - nofish
      data.frame(
        Location = loc,
        Scenario = sc,
        BiomassType = biomass_names,
        Change = diff
      )
    }))
  }))
  
  # Plot with ggplot
  ggplot(change_df, aes(x = Scenario, y = Change, fill = BiomassType)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(Location ~ ., scales = "free_y") +
    scale_fill_manual(values = fish_colors) +
    geom_vline(xintercept = seq(1.5, length(unique(change_df$Scenario)) - 0.5, by = 1),
               color = "black", linetype = "solid", linewidth = 0.2) +
    theme_minimal(base_size = 14)+
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
    labs(
      y = expression("Change in Biomass in [g" ~ m^-2 * "]"),
      x = " ",
      fill = " "
    )
}




# -----------------------------------------------------------------------------

#plot_no_fishing()
#plot_biomass_changes()
#plotNetwork(sim_results[["Shelf Sea"]][["No Fishing"]])
#plotNetwork(sim_results[["Open Ocean"]][["No Fishing"]])
#plot_changes("Shelf Sea")
#plot_changes("Slope")
#plot_changes("Open Ocean")
