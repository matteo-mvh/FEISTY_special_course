# -----------------------------------------------------------------------------
# Description:
#   This script runs FEISTY simulations for different ocean locations and 
#   fishing scenarios. The workflow is:
#     1. Define global loop settings (e.g., simulation length, life stages, etc.)
#     2. Define locations with their biological productivity and depth features
#     3. Define fishing scenarios (different exploitation strategies)
#     4. Define valid combinations of locations and fishing scenarios
#     5. Run FEISTY simulations for each valid combination
#     6. Store results in a structured list
#     7. Save results to a .RData file for later analysis
# -----------------------------------------------------------------------------

# Set results directory
results_dir <- "C:/Users/Mmm/OneDrive/Master Studies/3. Semester/Carbon Sequesteration/Code/results_files"

# Load necessary library
library(FEISTY)
setwd("C:/Users/Mmm/OneDrive/Master Studies/3. Semester/Carbon Sequesteration/Code")
source('FEISTY_carbon.R')

# Define common loop settings (shared across all runs)
loop_settings <- list(
  nStages    = 9,
  visual     = 1.5,  
  etaMature  = 0.25,  
  tEnd       = 500,   
  dfbot      = NA,    
  dfpho      = NA
)

# Define the parameter sets for 3 locations
locations <- list(
  list(name = 'Shelf Sea' , szprod = 100, lzprod = 90, bprod = 35, depth =   75, 
                            shelf =  75, Tp = NA, Tm = NA, Tb = NA, photic = 100),
  
  list(name = 'Slope'     , szprod =  60, lzprod = 50, bprod =  5, depth = 1500,
                            shelf = 250, Tp = NA, Tm = NA, Tb = NA, photic = 100),
  
  list(name = 'Open Ocean', szprod =  10, lzprod =  8, bprod =  0, depth = 3000,
                            shelf = 250, Tp = NA, Tm = NA, Tb = NA, photic = 100)
)

# Define fishing scenarios
# 1 = small pelagic, 2 = mesopelagic, 3 = large pelagic,
# 4 = midwater predators, 5 = demersal
fishing_scenarios <- list(
  list(name = 'No Fishing'   ,Fmax1 = 0.0, etaF1 = 0.05, groupidx1 = c(1), 
                              Fmax2 = 0.0, etaF2 = 0.05, groupidx2 = c(1)),
  
  list(name = 'Demersal'     ,Fmax1 = 0.3, etaF1 = 0.05, groupidx1 = c(5), 
                              Fmax2 = 0.0, etaF2 = 0.05, groupidx2 = c(1)),
  
  list(name = 'Forage Fish'  ,Fmax1 = 0.3, etaF1 = 0.05, groupidx1 = c(5), 
                              Fmax2 = 0.6, etaF2 = 0.05, groupidx2 = c(1)),
  
  list(name = 'Large Pelagic',Fmax1 = 0.3, etaF1 = 0.05, groupidx1 = c(3), 
                              Fmax2 = 0.0, etaF2 = 0.05, groupidx2 = c(5))
)

# Define valid combinations of (location, fishing scenario)
valid_combos <- list(
  "Shelf Sea"  = c("No Fishing", "Demersal", "Forage Fish"                 ),
  "Slope"      = c("No Fishing", "Demersal", "Forage Fish", "Large Pelagic"),
  "Open Ocean" = c("No Fishing",             "Forage Fish"                 )
)

# Initialize a list to store simulations
sim_results   <- list()
avg_Biomass   <- list()
total_flux    <- list()
carbon_inject <- list()

# Loop over each location
for (loc in locations) {
  # Initialize a sub-list for this location
  sim_results[[loc$name]] <- list()
  
  # Loop over each fishing scenario
  for (fish in fishing_scenarios) {
    
    # Check if this (location, fishing) combo is valid
    if (fish$name %in% valid_combos[[loc$name]]) {
      
      # Run valid simulation
      # Set up parameters for this location
      p <- setupVertical2(
        szprod     = loc$szprod,
        lzprod     = loc$lzprod,
        bprodin    = loc$bprod,
        dfbot      = loop_settings$dfbot,
        dfpho      = loop_settings$dfpho,
        nStages    = loop_settings$nStages,
        Tp         = loc$Tp,
        Tm         = loc$Tm,
        Tb         = loc$Tb,
        depth      = loc$depth,
        photic     = loc$photic,
        shelfdepth = loc$shelf,
        visual     = loop_settings$visual,
        etaMature  = loop_settings$etaMature
      )
      
      # Set fishing parameters for this scenario
      p <- setFishing(p, 
                      Fmax     = fish$Fmax1, 
                      etaF     = fish$etaF1, 
                      groupidx = fish$groupidx1
      )
      p <- setFishing(p, 
                      Fmax     = fish$Fmax2, 
                      etaF     = fish$etaF2, 
                      groupidx = fish$groupidx2
      ) 
      
      # Run simulation
      sim <- simulateFEISTY(p = p, tEnd = loop_settings$tEnd)
      
      df <- sim[["totBiomass"]]
      # take last 80 rows
      last80 <- tail(df, 80)
      # average each column
      avgBio <-colMeans(last80, na.rm = TRUE)
      
      # Calculate carbon fluxes at the position of the fish:
      sim = calcCarbonFluxes(sim) 
      totalflux = sim$fluxCarcass + sim$fluxFecal + sim$fluxRepro + sim$fluxRespiration
      
      # Calculate the injection
      inject = calcCarbonInjection(sim)
      
      #summing up the total injection inject$total
      inject$total_sum <- sum(inject$total)
  
      # Store result
      avg_Biomass  [[loc$name]][[fish$name]] <- avgBio
      sim_results  [[loc$name]][[fish$name]] <- sim
      total_flux   [[loc$name]][[fish$name]] <- totalflux
      carbon_inject[[loc$name]][[fish$name]] <- inject
      
      # Print progress update
      print(paste("Simulation completed:", loc$name, "-", fish$name))
    } else {
      
      # Report skipped simulation
      print(paste("Simulation skipped:"  , loc$name, "-", fish$name))
    }
  }
}


# Ensure the directory exists
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Save the results to that directory
save(sim_results,   file = file.path(results_dir, "FEISTY_Simulations_Results.RData"))
save(avg_Biomass,   file = file.path(results_dir, "FEISTY_Biomass_Results.RData"))
save(total_flux,    file = file.path(results_dir, "FEISTY_Flux_Results.RData"))
save(carbon_inject, file = file.path(results_dir, "FEISTY_Injection_Results.RData"))

# Final message
print(paste("All simulations completed and results saved in:", 
            file.path(results_dir, "FEISTY_Simulations_Results.RData")))

# End of script

