# ----------------------------------------
# Global run of FEISTY for all fishing scenarios
# ----------------------------------------

library(FEISTY)
source('scripts/FEISTY_carbon.R')

# Fishing Parameters
fishing_scenarios <- list(
  list(name = 'No_Fishing'   ,Fmax1 = 0.0, etaF1 = 0.05, groupidx1 = c(1), 
       Fmax2 = 0.0, etaF2 = 0.05, groupidx2 = c(1)),
  
  list(name = 'Demersal'     ,Fmax1 = 0.3, etaF1 = 0.05, groupidx1 = c(5), 
       Fmax2 = 0.0, etaF2 = 0.05, groupidx2 = c(1)),
  
  list(name = 'Forage_Fish'  ,Fmax1 = 0.3, etaF1 = 0.05, groupidx1 = c(5), 
       Fmax2 = 0.6, etaF2 = 0.05, groupidx2 = c(1)),
  
  list(name = 'Large_Pelagic',Fmax1 = 0.3, etaF1 = 0.05, groupidx1 = c(3), 
       Fmax2 = 0.0, etaF2 = 0.05, groupidx2 = c(5))
)

# Load global forcing data
glob <- read.csv(file = "data/Input_global.csv")

# Function to run FEISTY for one row
simulateFEISTY_single <- function(rowidx, glob, Fmax1, etaF1, groupidx1, Fmax2, etaF2, groupidx2) {
  p_sim = setupVertical2(
    szprod = glob[rowidx, "szprod"],
    lzprod = glob[rowidx, "lzprod"],
    dfbot  = glob[rowidx, "dfbot"],
    photic = glob[rowidx, "photic"],
    depth  = glob[rowidx, "depth"],
    Tp     = glob[rowidx, "Tp"],
    Tm     = glob[rowidx, "Tm"],
    Tb     = glob[rowidx, "Tb"])
  
  p_sim = setFishing(p_sim, Fmax = Fmax1, etaF = etaF1, groupidx = groupidx1)
  p_sim = setFishing(p_sim, Fmax = Fmax2, etaF = etaF2, groupidx = groupidx2)
  
  sim <- simulateFEISTY(p = p_sim)
  return(sim)
}

# Loop over all fishing scenarios
for (list_idx in seq_along(fishing_scenarios)) {
  
  scenario <- fishing_scenarios[[list_idx]]
  
  Fmax1     <- scenario$Fmax1
  etaF1     <- scenario$etaF1
  groupidx1 <- scenario$groupidx1
  Fmax2     <- scenario$Fmax2
  etaF2     <- scenario$etaF2
  groupidx2 <- scenario$groupidx2
  
  all_results <- list()
  
  for (rowidx in 1:nrow(glob)) {
    sim <- simulateFEISTY_single(rowidx, glob, Fmax1, etaF1, groupidx1, Fmax2, etaF2, groupidx2)
    
    sime <- calcCarbonFluxes(sim)
    inject <- calcCarbonInjection(sime)
    inject$total_sum <- sum(inject$total)
    
    res <- list(
      lon = glob[rowidx, "lon"],
      lat = glob[rowidx, "lat"],
      Biomass = colMeans(sim$totBiomass[round(0.6 * sim$nTime):sim$nTime, ]),
      totBiomass = sum(colMeans(sim$totBiomass[round(0.6 * sim$nTime):sim$nTime, ])),
      inject = inject$total_sum
    )
    
    all_results[[rowidx]] <- res
  }
  
  # Create output data frame
  out <- t(data.frame(matrix(data = unlist(all_results), nrow = 9)))
  colnames(out) <- c("lon", "lat", "totB_smpel", "totB_mesopel", 
                     "totB_largepel", "totB_midwpred", "totB_dem", "totB_all", "carbon_inject")
  out <- as.data.frame(out)
  
  # Save RData with scenario name
  file_rdata <- paste0("data/Global_fish_biomass_", scenario$name, ".RData")
  save(out, file = file_rdata)
  
  cat("✅ Finished scenario:", scenario$name, "\n")
}

