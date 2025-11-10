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

# --- Setup environment -------------------------------------------------------

# Set results directory
results_dir <- "C:/Users/Mmm/OneDrive/Master Studies/3. Semester/Carbon Sequesteration/FEISTY_special_course/Local/results_files"

# Load necessary library
library(FEISTY)
setwd("C:/Users/Mmm/OneDrive/Master Studies/3. Semester/Carbon Sequesteration/FEISTY_special_course/Local")
source('scripts/FEISTY_carbon.R')

# Define Fmax and etaF value ranges
<<<<<<< HEAD
Fmax_values <- seq(0, 1, by = 0.1)
etaF_values <- seq(0, 0.5, by = 0.025)
=======
Fmax_values <- seq(0, 2, by = 0.1)
etaF_values <- seq(0, 0.5, by = 0.1)
>>>>>>> de11c71 (updated the local folder)

# Create result lists before loops
avg_Biomass   <- list()
carbon_inject <- list()

# --- Define loop settings ----------------------------------------------------
loop_settings <- list(
  nStages    = 9,
  visual     = 1.5,  
  etaMature  = 0.25,  
  tEnd       = 500,   
  dfbot      = NA,    
  dfpho      = NA
)

# --- Define location parameters ----------------------------------------------
locations <- list(
  list(name = 'Shelf Sea' , szprod = 100, lzprod = 90, bprod = 35, depth =   75, 
       shelf =  75, Tp = NA, Tm = NA, Tb = NA, photic = 100),
  
  list(name = 'Slope'     , szprod =  60, lzprod = 50, bprod =  5, depth = 1500,
       shelf = 250, Tp = NA, Tm = NA, Tb = NA, photic = 100),
  
  list(name = 'Open Ocean', szprod =  10, lzprod =  8, bprod =  0, depth = 3000,
       shelf = 250, Tp = NA, Tm = NA, Tb = NA, photic = 100)
)

# --- Define fishing scenarios ------------------------------------------------
# 1 = small pelagic, 2 = mesopelagic, 3 = large pelagic,
# 4 = midwater predators, 5 = demersal
fishing_scenarios <- list(
  list(name = 'Demersal'     , groupidx1 = c(5), groupidx2 = c(5)),
  list(name = 'Forage Fish'  , groupidx1 = c(1), groupidx2 = c(5)),
<<<<<<< HEAD
  list(name = 'Large Pelagic', groupidx1 = c(3), groupidx2 = c(3))
)

# --- Define valid combinations ----------------------------------------------
valid_combos <- list(
  "Shelf Sea"  = c("Demersal", "Forage Fish"),
  "Slope"      = c("Demersal", "Forage Fish", "Large Pelagic"),
  "Open Ocean" = c("Forage Fish")
=======
  list(name = 'Large Pelagic', groupidx1 = c(3), groupidx2 = c(3)),
  list(name = 'All Fish'     , groupidx1 = c(1,2,3,4,5), groupidx2 = c(1,2,3,4,5))
  )

# --- Define valid combinations ----------------------------------------------
valid_combos <- list(
  "Shelf Sea"  = c("Demersal", "Forage Fish",                 "All Fish"),
  "Slope"      = c("Demersal", "Forage Fish", "Large Pelagic","All Fish"),
  "Open Ocean" = c(            "Forage Fish",                 "All Fish")
)


valid_combos <- list(
  "Slope"      = c("All Fish")
>>>>>>> de11c71 (updated the local folder)
)

# --- Main simulation loop ----------------------------------------------------
for (loc in locations) {
  loc_name <- loc$name
  avg_Biomass[[loc_name]]   <- list()
  carbon_inject[[loc_name]] <- list()
  
  for (fish in fishing_scenarios) {
    fish_name <- fish$name
    
    # Skip invalid combos
    if (!(fish_name %in% valid_combos[[loc_name]])) {
      cat("Skipping:", loc_name, "-", fish_name, "\n")
      next
    }
    
    cat("Running:", loc_name, "-", fish_name, "\n")
    
    # Prepare result matrices
    nF <- length(Fmax_values)
    nE <- length(etaF_values)
    
    biomass_matrix <- matrix(NA, nrow = nF, ncol = nE,
                             dimnames = list(paste0("Fmax=", Fmax_values),
                                             paste0("etaF=", etaF_values)))
    inject_matrix  <- matrix(NA, nrow = nF, ncol = nE,
                             dimnames = list(paste0("Fmax=", Fmax_values),
                                             paste0("etaF=", etaF_values)))
    
    # --- Double loop over Fmax and etaF ---
    for (i in seq_along(Fmax_values)) {
      cat("  Fmax =", Fmax_values[i], "\n")
      flush.console()
      for (j in seq_along(etaF_values)) {
        
        Fmax <- Fmax_values[i]
        etaF <- etaF_values[j]
        
        # Set up parameters for this run
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
        
        # Apply fishing pressure
        p <- setFishing(p, Fmax = Fmax, etaF = etaF, groupidx = fish$groupidx1)
        if (!identical(fish$groupidx1, fish$groupidx2)) {

          p <- setFishing(p, Fmax = 2 * Fmax, etaF = etaF, groupidx = fish$groupidx2)
        }
        
        # Run simulation safely
        sim <- try(simulateFEISTY(p = p, tEnd = loop_settings$tEnd), silent = TRUE)
        if (inherits(sim, "try-error")) next
        
        # Compute biomass average (last 80 timesteps)
        df <- sim$totBiomass
        last80 <- tail(df, 80)
        avgBio <- mean(colMeans(last80, na.rm = TRUE), na.rm = TRUE)
        
        # Compute carbon injection
        sim_flux <- calcCarbonFluxes(sim)
        inject <- calcCarbonInjection(sim_flux)
        total_inject <- sum(inject$total, na.rm = TRUE)
        
        # Store results
        biomass_matrix[i, j] <- avgBio
        inject_matrix[i, j]  <- total_inject
      }
    }
    
    # Store matrices for this scenario
    avg_Biomass[[loc_name]][[fish_name]]   <- biomass_matrix
    carbon_inject[[loc_name]][[fish_name]] <- inject_matrix
    
    cat("✅ Finished:", loc_name, "-", fish_name, "\n")
  }
}

# --- Save results ------------------------------------------------------------
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

save(avg_Biomass, carbon_inject, 
     file = file.path(results_dir, "FEISTY_FishingPressure_Matrix.RData"))

cat("\n✅ All simulations completed and saved in:",
    file.path(results_dir, "FEISTY_FishingPressure_Matrix.RData"), "\n")
