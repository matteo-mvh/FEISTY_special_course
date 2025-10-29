# ----------------------------------------
# Global run of FEISTY (sequential version)

# Load libraries
library(FEISTY)

# Load global forcing data
glob <- read.csv(file = "data/Input_global.csv")

# Create function to run FEISTY for one row
simulateFEISTY_single <- function(rowidx, glob) {
  sim <- simulateFEISTY(
    p = setupVertical2(
      szprod = glob[rowidx, "szprod"],
      lzprod = glob[rowidx, "lzprod"],
      dfbot  = glob[rowidx, "dfbot"],
      photic = glob[rowidx, "photic"],
      depth  = glob[rowidx, "depth"],
      Tp     = glob[rowidx, "Tp"],
      Tm     = glob[rowidx, "Tm"],
      Tb     = glob[rowidx, "Tb"])
  )
  return(sim)
}

# Initialize results list
all_results <- list()

# Run model sequentially
for (rowidx in 1:nrow(glob)) {
  sim <- simulateFEISTY_single(rowidx, glob)
  res <- list(
    lon = glob[rowidx, "lon"],
    lat = glob[rowidx, "lat"],
    totBiomass = colMeans(sim$totBiomass[round(0.6 * sim$nTime):sim$nTime, ])
  )
  all_results[[rowidx]] <- res
}

# Create output data frame
out <- t(data.frame(matrix(data = unlist(all_results), nrow = 7)))
colnames(out) <- c("lon", "lat", "totB_smpel", "totB_mesopel", 
                   "totB_largepel", "totB_midwpred", "totB_dem")
out <- as.data.frame(out)

# Save results
save(out, file = "data/Global_fish_biomass.RData")
