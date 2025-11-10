# ============================================================
# Load results
# ============================================================
results_dir <- "C:/Users/Mmm/OneDrive/Master Studies/3. Semester/Carbon Sequesteration/FEISTY_special_course/Local/results_files"
load(paste0(results_dir, "/FEISTY_FishingPressure_Matrix.RData"))

# Example available keys
loc_names  <- names(avg_Biomass)
fish_names <- names(avg_Biomass[[1]])

# ============================================================
# Select what to plot
# ============================================================
loc_to_plot  <- loc_names[1]        # first location
fish_to_plot <- fish_names[3]       # first fish group


biomass_matrix <- avg_Biomass[[loc_to_plot]][[fish_to_plot]]
inject_matrix  <- carbon_inject[[loc_to_plot]][[fish_to_plot]]

# ============================================================
# Define Fmax and etaF values (must match how matrices were created)
# ============================================================
# Replace these with your actual sequences
Fmax_values <- seq(0, 1, length.out = nrow(biomass_matrix))
etaF_values <- seq(0, 0.2, length.out = ncol(biomass_matrix))


# ============================================================
# Convert matrices to data frames with Fmax and etaF
# ============================================================
library(reshape2)
library(ggplot2)
library(gridExtra)

df_biomass <- melt(biomass_matrix)
df_biomass$Fmax <- Fmax_values[df_biomass$Var1]
df_biomass$etaF <- etaF_values[df_biomass$Var2]
colnames(df_biomass) <- c("Row", "Col", "Biomass", "Fmax", "etaF")

df_inject <- melt(inject_matrix)
df_inject$Fmax <- Fmax_values[df_inject$Var1]
df_inject$etaF <- etaF_values[df_inject$Var2]
colnames(df_inject) <- c("Row", "Col", "Inject", "Fmax", "etaF")

# ============================================================
# Plot biomass as a contour field
# ============================================================
p1 <- ggplot(df_biomass, aes(x = Fmax, y = etaF, z = Biomass)) +
  geom_contour_filled() +
  scale_fill_viridis_d(option = "plasma", name = "Biomass") +
  labs(x = "Fmax", y = "etaF") +
  theme_minimal(base_size = 14)

# ============================================================
# Plot carbon injection as a contour field
# ============================================================
p2 <- ggplot(df_inject, aes(x = Fmax, y = etaF, z = Inject)) +
  geom_contour_filled() +
  scale_fill_viridis_d(option = "magma", name = "Carbon Injection") +
  labs(x = "Fmax", y = "etaF") +
  theme_minimal(base_size = 14)

# ============================================================
# Arrange side by side with shared title
# ============================================================
library(gridExtra)
grid.arrange(p1, p2, ncol = 2, 
             top = paste("Location:", loc_to_plot, " | Fishing Type:", fish_to_plot))
