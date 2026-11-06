# Avg_Connectivity_Comp.R
#Spin test comparing average connectivity of HCP-D and HCP-YA

library(tidyr)
library(dplyr)
library(readxl)

# === Input data
HCP_YA_con <- read.csv("insuladev/HCP_YA/averaged_vertex_connectivity.csv", row.names = 1)
HCPD_con <- read.csv("insuladev/dropbox/HCPD_avg_con.csv")

HCPD_con$Insula <- paste(HCPD_con$Insula_Hemisphere, "_",HCPD_con$Insula_Subregion, sep = "")

HCP_YA_con_long <- HCP_YA_con %>%
  dplyr::rename(Insula = row_names) %>%
  tidyr::pivot_longer(
    cols = -Insula,
    names_to = "Glasser_ROI",
    values_to = "HCPYA_connectivity"
  )

HCP_YA_con_long$Glasser_ROI <- paste(HCP_YA_con_long$Glasser_ROI, "_ROI", sep = "")
HCP_YA_con_long <- HCP_YA_con_long %>%
  dplyr::mutate(across(
    c(Insula, Glasser_ROI),
    ~ gsub("Left", "L", gsub("Right", "R", as.character(.)))
  ))

avg_con <- merge(HCP_YA_con_long, HCPD_con, by = c("Insula", "Glasser_ROI"))
avg_con$GLASSERLABELNAME <- avg_con$Glasser_ROI

source("insuladev/code/spintests/rotate.parcellation.R")

#perm.sphere.p with slight modification
perm.sphere.p = function(x,y,perm.id,corr.type='spearman') {
  
  nroi = dim(perm.id)[1]  # number of regions
  nperm = dim(perm.id)[2] # number of permutations
  
  rho.emp = cor(x,y,method=corr.type,use="complete.obs")  # empirical correlation
  
  # permutation of measures
  x.perm = y.perm = array(NA,dim=c(nroi,nperm))
  for (r in 1:nperm) {
    for (i in 1:nroi) {
      x.perm[i,r] = x[perm.id[i,r]]
      y.perm[i,r] = y[perm.id[i,r]]
    }
  }
  
  # correlation to unpermuted measures
  rho.null.xy = rho.null.yx = vector(length=nperm)
  for (r in 1:nperm) {
    rho.null.xy[r] = cor(x.perm[,r],y,method=corr.type,use="complete.obs")
    rho.null.yx[r] = cor(y.perm[,r],x,method=corr.type,use="complete.obs")
  }
  
  # p-value definition depends on the sign of the empirical correlation
  # added +1 to numerator and denominator because otherwise p value is 0 
  # (because the observed correlation is greater than all permuted correlations)
  if (rho.emp>0) {
    p.perm.xy = (sum(rho.null.xy>rho.emp)+1)/(nperm + 1)
    p.perm.yx = (sum(rho.null.yx>rho.emp)+1)/(nperm + 1)
  } else { 
    p.perm.xy = (sum(rho.null.xy<rho.emp)+1)/(nperm + 1)
    p.perm.yx = (sum(rho.null.yx<rho.emp)+1)/(nperm + 1)
  } 
  
}

spherical.coordinates <- read.table("insuladev/dropbox/sphere_HCP.txt", header=F) #Coordinates taken from the sphere_HCP.txt file provided here https://github.com/frantisekvasa/rotate_parcellation
spherical.coordinates <- as.matrix(spherical.coordinates)

cortex_coords_l <- as.matrix(spherical.coordinates[1:180, ])
cortex_coords_r <- as.matrix(spherical.coordinates[181:360,])

set.seed(132)
perm.id <- rotate.parcellation(cortex_coords_l, cortex_coords_r, nrot=1000)

Cole_parcels = read_excel("insuladev/dropbox/Cole_Parcels_Networks.xlsx")

cole_filtered <- Cole_parcels[Cole_parcels$KEYVALUE %in% 1:360, ]

# === Spin Test ===
perform_spin_test <- function(insula_name, n_permutations = 1000) {
  # Filter the data for the given insula subregion
  PCA_Insula <- avg_con %>% filter(Insula == insula_name)
  PCA_Insula_ordered <- left_join(
    cole_filtered, 
    PCA_Insula, 
    by = "GLASSERLABELNAME", 
    keep = TRUE
  )
  observed_correlation <- cor(PCA_Insula_ordered$connectivity, 
                              PCA_Insula_ordered$HCPYA_connectivity, 
                              method = "spearman",
                              use = "complete.obs")
  
  # Perform the spin test (permutation test)
  p_value <- perm.sphere.p(
    x = PCA_Insula_ordered$connectivity,
    y = PCA_Insula_ordered$HCPYA_connectivity,
    perm.id = perm.id,
    corr.type = "spearman"
  )
  
  # Return the plot along with the observed correlation and p-value
  return(list(observed_correlation = observed_correlation, 
              p_value = p_value))
}

# Run the spin test for each insula in the desired order
insula_names <- unique(avg_con$Insula)
spin_list <- setNames(lapply(insula_names, perform_spin_test), insula_names)

results_table <- data.frame(
  Insula = names(spin_list),
  Observed_Correlation = sapply(spin_list, function(x) x$observed_correlation),
  P_Value = sapply(spin_list, function(x) x$p_value)
)

results_table$Significant <- results_table$P_Value < 0.05
results_table$FDRadjusted_p_values <- p.adjust(results_table$P_Value, method = "BH")

print(results_table)

saveRDS(results_table, "insuladev/dropbox/HCP_YA_D_results_table.RDS")