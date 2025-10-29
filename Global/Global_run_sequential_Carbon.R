# ----------------------------------------
# Global run of FEISTY (sequential version)

# Load libraries
library(FEISTY)
source('scripts/FEISTY_carbon.R')

# Fishing Parameters
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
# Select fishing scenario to run

####################
list_idx  <- 1  # Change this index to select different fishing scenarios
####################

Fmax1     <- fishing_scenarios[[list_idx]]$Fmax1
etaF1     <- fishing_scenarios[[list_idx]]$etaF1
groupidx1 <- fishing_scenarios[[list_idx]]$groupidx1
Fmax2     <- fishing_scenarios[[list_idx]]$Fmax2
etaF2     <- fishing_scenarios[[list_idx]]$etaF2
groupidx2 <- fishing_scenarios[[list_idx]]$groupidx2

# Load global forcing data
glob <- read.csv(file = "data/Input_global.csv")

# Create function to run FEISTY for one row
simulateFEISTY_single <- function(rowidx, glob) {
  p_sim = setupVertical2(
    szprod = glob[rowidx, "szprod"],
    lzprod = glob[rowidx, "lzprod"],
    dfbot  = glob[rowidx, "dfbot"],
    photic = glob[rowidx, "photic"],
    depth  = glob[rowidx, "depth"],
    Tp     = glob[rowidx, "Tp"],
    Tm     = glob[rowidx, "Tm"],
    Tb     = glob[rowidx, "Tb"])
  p_sim = setFishing(p_sim, 
                     Fmax     = Fmax1, 
                     etaF     = etaF1, 
                     groupidx = groupidx1)
  p_sim = setFishing(p_sim, 
                     Fmax     = Fmax2, 
                     etaF     = etaF2, 
                     groupidx = groupidx2)
  sim <- simulateFEISTY(p = p_sim)
  return(sim)
}

# Initialize results list
all_results <- list()

# Run model sequentially
for (rowidx in 1:nrow(glob)) {
  sim <- simulateFEISTY_single(rowidx, glob)
  
  # Calculate the injection
  # Calculate carbon fluxes at the position of the fish:
  sime = calcCarbonFluxes(sim) 
  totalflux = sime$fluxCarcass + sime$fluxFecal + sime$fluxRepro + sime$fluxRespiration
  inject = calcCarbonInjection(sime)
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

# Save results
save(out, file = "data/Global_fish_biomass.RData")
