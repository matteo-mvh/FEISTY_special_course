library(FEISTY)

# Set directories
base_dir <- "C:/Users/Mmm/OneDrive/Master Studies/3. Semester/Carbon Sequesteration"
global_dir <- file.path(base_dir, "FEISTY_special_course/Global")

setwd(base_dir)
load("CTL.R")  # loads TM, matrixInject, etc.

setwd(global_dir)
source("scripts/FEISTY_carbon.R")

# List of fishing scenarios and their files
scenarios <- list(
  "No_Fishing"     = "data/Global_fish_biomass_No_Fishing.RData",
  "Large_Pelagic"  = "data/Global_fish_biomass_Large_Pelagic.RData",
  "Demersal"       = "data/Global_fish_biomass_Demersal.RData",
  "Forage_Fish"    = "data/Global_fish_biomass_Forage_Fish.RData"
)

# Loop over each scenario safely
for (scenario in names(scenarios)) {
  
  cat("Processing scenario:", scenario, "...\n")
  
  # Load scenario data
  load(scenarios[[scenario]])  # loads TM, matrixInject, etc.
  
  # Calculate sequestration
  sequestration <- calc_CarbonSequestration(TM, matrixInject)
  
  # Save the result with a descriptive name
  save_file <- file.path(global_dir, paste0("data/CarbonSequestration_", scenario, ".RData"))
  save(sequestration, file = save_file)
  
  # Remove large objects and free memory
  rm(sequestration, TM, matrixInject)
  gc()
  
  cat("Finished scenario:", scenario, "\n")
}

cat("All scenarios processed.\n")
