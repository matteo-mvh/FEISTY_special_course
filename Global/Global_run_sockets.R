
# ----------------------------------------
# Global run of FEISTY 

# load libraries
library(FEISTY)
library(foreach)
library(doParallel)

# Load global forcing data
glob <- read.csv(file = "data/Input_global.csv")

# Set cores:
cl <- makeCluster(detectCores() - 2) 
registerDoParallel(cl)

# Create function:
simulateFEISTY_parallel <- function(rowidx, glob) {
  sim = simulateFEISTY(
    p = setupVertical2(
      szprod = glob[rowidx, "szprod"],
      lzprod = glob[rowidx, "lzprod"],
      dfbot  = glob[rowidx, "dfbot"],
      photic = glob[rowidx, "photic"],
      depth  = glob[rowidx, "depth"],
      Tp     = glob[rowidx, "Tp"],
      Tm     = glob[rowidx, "Tm"],
      Tb     = glob[rowidx, "Tb"]))
  return(sim)
}

# Create empty list
res         <- list()
all_results <- list()

# Run model:
all_results <- foreach(rowidx = 1:nrow(glob), .packages = c("FEISTY")) %dopar% {
  sim = simulateFEISTY_parallel(rowidx, glob)
  res$lon=glob[rowidx,"lon"]
  res$lat=glob[rowidx,"lat"]
  res$totBiomass=colMeans(sim$totBiomass[round(0.6*sim$nTime):sim$nTime,]) # last 40% simulate time
  return(res)
}

stopCluster(cl)

# Create a data frame
out <- t(data.frame(matrix(data=unlist(all_results),nrow=7)))
colnames(out) <- c("lon","lat","totB_smpel","totB_mesopel","totB_largepel","totB_midwpred", "totB_dem")
out <- as.data.frame(out)
save(out, file="data/Global_fish_biomass.RData")
