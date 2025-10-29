# Load data
load(file.path(results_dir, "FEISTY_Injection_Results.RData"))

# Define the fishing scenarios
scenarios <- c("No Fishing", "Demersal", "Forage Fish", "Large Pelagic")
colors <- c("black", "green", "red", "blue")  # assign colors for each scenario
locations <- c("Shelf Sea", "Slope", "Open Ocean")

for (Location in locations) {
  # Prepare plot limits
  all_depths <- carbon_inject[[Location]][[scenarios[1]]]$z
  ylim <- rev(range(all_depths))
  
  if (Location == "Shelf Sea") {
    xlim <- c(0.3, 0.5)
  } else if (Location == "Slope") {
    xlim <- c(0, 0.03)
  } else if (Location == "Open Ocean") {
    xlim <- c(0, 0.005)
  }
  
  # Plot the first scenario
  plot(carbon_inject[[Location]][[scenarios[1]]]$total,
       carbon_inject[[Location]][[scenarios[1]]]$z,
       type = "l", lwd = 2, col = colors[1],
       xlab = "Total Carbon Injection (gC/m³/yr)",
       ylab = "Depth (m)",
       main = paste(Location),
       xlim = xlim, ylim = ylim)
  
  # Add the other scenarios
  for (i in 2:length(scenarios)) {
    inject <- carbon_inject[[Location]][[scenarios[i]]]
    lines(inject$total, inject$z, col = colors[i], lwd = 2)
  }
  
  # Add legend
  legend("bottomright", legend = scenarios, col = colors, lwd = 2)
}

