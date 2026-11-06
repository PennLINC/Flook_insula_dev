###Parcellated Connectivity and GAM
#Run after ROI_GAM_Setup.R (via ex.sh SLURM script)
#Purpose: Fit GAMs to insula-ROI connectivity, extract age-related effects,
# perform PVD (two PCAs), and compare to SA axis and Cole Networks

# === Load Required Libraries ===
library(mgcv)
library(gratia)
library(factoextra)
library(ciftiTools)
library(cifti)
library(readxl)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')

# === Load Data ===
ROI_connectivity <- read_csv("insuladev/ROI_GAM/ROI_single_line.csv")
phenotypes <- read.csv("insuladev/HCPDbold_1/cor_results/model_array_demographics.csv")

S.A.axis.cifti <- cifti::read_cifti("insuladev/dropbox/SensorimotorAssociation.Axis.Glasser360.pscalar.nii") #S-A_ArchetypalAxis repo
#file can be found at https://github.com/PennLINC/S-A_ArchetypalAxis/
S.A.axis <- data.frame(SA.axis = rank(S.A.axis.cifti$data), orig_parcelname = names(S.A.axis.cifti$Parcel))


# === Preprocess Data ===
for_PCA_cleaned <- ROI_connectivity %>%
  select(-subject) %>%
  select(-contains("???"))  # remove medial wall

# === GAM Wrapper Function ===
run_ROI_gam <- function(column_name) {
  current_data <- data.frame(
    connectivity = for_PCA_cleaned[[column_name]],
    age = phenotypes$age,
    sex = as.factor(phenotypes$sex),
    mean_fd = phenotypes$mean_fd
  )
  
  #full model
  gam_model <- gam(connectivity ~ s(age, k = 3) + sex + mean_fd, 
                   data = current_data, method = "REML")
  gam_results <- summary(gam_model)
  
  #null model (without age)
  nullmodel <- gam(connectivity ~ sex + mean_fd, method = "REML", 
                   data = current_data)
  
  #compare null and full model to get partial R2
  RSS_model <- sum((gam_model$y - gam_model$fitted.values)^2)
  RSS_nullmodel <- sum((nullmodel$y - nullmodel$fitted.values)^2)
  partialRsq <- (RSS_nullmodel - RSS_model) / RSS_nullmodel
  
  #calculate effect direction
  derv <- derivatives(gam_model, interval = "simultaneous", unconditional = FALSE)
  mean_derivative <- mean(derv$.derivative)
  if (mean_derivative < 0) partialRsq <- -partialRsq
  
  anova_smooth_pvalue <- anova.gam(nullmodel, gam_model, test = 'F')$`Pr(>F)`[2]
  
  # Prediction function 
    new_data <- data.frame(
      age = seq(min(current_data$age), max(current_data$age), length.out = 100),
      sex = rep("Male", 100),
      mean_fd = mean(current_data$mean_fd)
    )
    preds <- predict(gam_model, newdata = new_data, type = "terms", se.fit = TRUE)
    smooth_values <- preds$fit[, "s(age)"]
    se_values <- preds$se.fit[, "s(age)"]
    lower_ci <- smooth_values - 1.96 * se_values
    upper_ci <- smooth_values + 1.96 * se_values
    
    return(data.frame(
      ROI = column_name,
      GAM_smooth_Fvalue = gam_results$s.table[3],
      GAM_smooth_pvalue = gam_results$s.table[4],
      GAM_smooth_partialR2 = partialRsq,
      Anova_smooth_pvalue = anova_smooth_pvalue,
      predicted_values = smooth_values,
      Lower_CI = lower_ci,
      Upper_CI = upper_ci,
      predicted_age = new_data$age
    ))
}

# === Run GAMs ===
column_names <- colnames(for_PCA_cleaned)
results_list <- lapply(column_names, run_ROI_gam)
final_results <- do.call(rbind, results_list)

# === Remove Insula-Insula Connections ===
insula_regions <- c("L_Ig", "L_PoI1", "L_PoI2", "L_MI", "L_AVI", "L_AAIC",
                    "R_Ig", "R_PoI1", "R_PoI2", "R_MI", "R_AVI", "R_AAIC")
pattern <- paste0("^(" , paste(insula_regions, collapse = "|"), ")_(.*)_ROI$")
final_results_filtered <- final_results %>% 
  mutate(second_region = str_match(ROI, pattern)[, 3]) %>%
  filter(!(second_region %in% insula_regions))

# === FDR Correction ===
# Apply FDR correction to one p-value per ROI
fdr_df <- final_results_filtered %>%
  group_by(ROI) %>%
  summarise(raw_p = first(Anova_smooth_pvalue)) %>%
  mutate(FDRadjusted_p = p.adjust(raw_p, method = "BH"))

# Join back to full predictions df
final_results_filtered <- final_results_filtered %>%
  left_join(fdr_df, by = "ROI") %>%
  mutate(significant_count = ifelse(FDRadjusted_p < 0.05, FDRadjusted_p, NA))

# === Prepare for Visualization ===
partialR2df <- final_results_filtered %>%
  select(ROI, GAM_smooth_partialR2) %>%
  rename(ROI_name = ROI, partialR2 = GAM_smooth_partialR2) %>%
  extract(ROI_name, into = c("Insula", "Brain"), regex = "^(.*?_.*?)_(.*)$") %>%
  group_by(Insula, Brain) %>%
  slice(1) %>%
  ungroup()

partialR2_wide <- partialR2df %>%
  select(Insula, Brain, partialR2) %>%
  pivot_wider(names_from = Brain, values_from = partialR2)

FDRdf <- final_results_filtered %>%
  select(ROI, significant_count) %>%
  rename(ROI_name = ROI, FDR_sig = significant_count) %>%
  extract(ROI_name, into = c("Insula", "Brain"), regex = "^(.*?_.*?)_(.*)$") %>%
  group_by(Insula, Brain) %>%
  slice(1) %>%
  ungroup()

FDRdf_wide <- FDRdf %>%
  select(Insula, Brain, FDR_sig) %>%
  pivot_wider(names_from = Brain, values_from = FDR_sig)

partialR2_masked <- partialR2_wide
partialR2_masked[is.na(FDRdf_wide)] <- NA

#output for how many values are significant
mean(!is.na(partialR2_masked[,2:349])) * 100

#save final_results with additional columns for visualization
final_results_filtered <- final_results_filtered %>% extract(ROI, into = c("Insula", "Brain"), regex = "^(.*?_.*?)_(.*)$")
final_results_filtered <- merge(final_results_filtered, S.A.axis, by.x = "Brain", by.y = "orig_parcelname", all.x = TRUE)

# === Save Outputs ===
saveRDS(partialR2_masked, "insuladev/dropbox/partialR2_masked.RDS")
saveRDS(partialR2_wide, "insuladev/dropbox/partialR2_wide.RDS")
saveRDS(partialR2df, "insuladev/dropbox/partialR2_sep.RDS")
saveRDS(final_results_filtered, "insuladev/dropbox/all_predictions.RDS")


# === PCA Analysis ===
pred_value_ROI <- final_results_filtered %>%
  select(Insula, Brain, predicted_values, predicted_age) %>%
  pivot_wider(
    names_from = c(Insula, Brain),
    values_from = predicted_values
  )

preds_ROI <- pred_value_ROI %>% select(-predicted_age)
rownames(preds_ROI) <- pred_value_ROI$predicted_age
transposed_preds_ROI <- t(preds_ROI)

Glasser_PCA <- prcomp(transposed_preds_ROI, scale. = FALSE) 
#scale is false as there are likely importance differences in the 
#variability across the sample, and the ages are 
#inherently on the same scale already

factoextra::fviz_eig(Glasser_PCA)
saveRDS(Glasser_PCA, "insuladev/dropbox/GlasserPVD_PCA1.RDS")

Glasser_pca_scores <- as.data.frame(Glasser_PCA$x)
Glasser_pca_scores$regions <- rownames(Glasser_pca_scores)
Glasser_pca_scores <- Glasser_pca_scores %>%
  extract(regions, into = c("Insula", "Brain"), regex = "^(.*?_.*?)_(.*)$")


# === SA Axis Comparison ===
PCA_scores_SA <- merge(S.A.axis, Glasser_pca_scores, by.x = "orig_parcelname", by.y = "Brain", all.x = TRUE, all.y = TRUE)
PCA_scores_SA_forviz <- PCA_scores_SA[complete.cases(PCA_scores_SA),]

saveRDS(PCA_scores_SA_forviz, "insuladev/dropbox/PCA_scores_SA_for_viz.RDS")

ordered_insula_names <- c("L_AAIC", "L_AVI", "L_MI", "L_PoI1", "L_PoI2", "L_Ig",
                          "R_AAIC", "R_AVI", "R_MI", "R_PoI1", "R_PoI2", "R_Ig")

insula_labels <- c("L_Ig"   = "Left Granular Insula", 
                   "L_PoI1" = "Left Posterior Insula 1", 
                   "L_PoI2" = "Left Posterior Insula 2", 
                   "L_MI"   = "Left Middle Insula", 
                   "L_AVI"  = "Left Anteroventral Insula", 
                   "L_AAIC" = "Left Agranular Insula",
                   "R_Ig"   = "Right Granular Insula",
                   "R_PoI1" = "Right Posterior Insula 1", 
                   "R_PoI2" = "Right Posterior Insula 2", 
                   "R_MI"   = "Right Middle Insula", 
                   "R_AVI"  = "Right Anteroventral", 
                   "R_AAIC" = "Right Agranular Insula")

#files for rotate parcellation, perm sphere and sphere HCP can be found 
#at https://github.com/frantisekvasa
source("insuladev/code/spintests/rotate.parcellation.R")
source("insuladev/code/spintests/perm.sphere.p.R")

spherical.coordinates <- read.table("insuladev/dropbox/sphere_HCP.txt", header=F) #Coordinates taken from the sphere_HCP.txt file provided here https://github.com/frantisekvasa/rotate_parcellation
spherical.coordinates <- as.matrix(spherical.coordinates)

cortex_coords_l <- as.matrix(spherical.coordinates[1:180, ])
cortex_coords_r <- as.matrix(spherical.coordinates[181:360,])

PCA_scores_SA_forviz$GLASSERLABELNAME <- PCA_scores_SA_forviz$orig_parcelname

Cole_parcels = read_excel("insuladev/dropbox/Cole_Parcels_Networks.xlsx")

cole_filtered <- Cole_parcels[Cole_parcels$KEYVALUE %in% 1:360, ]

PCA_scores_SA_forviz_ordered <- left_join(
  cole_filtered, 
  PCA_scores_SA_forviz, 
  by = "GLASSERLABELNAME", 
  keep = TRUE
)

set.seed(132)

perm.id <- rotate.parcellation(cortex_coords_l, cortex_coords_r, nrot=1000)

combined_insula_regions <- c("AVI", "AAIC", "MI", "PoI1", "PoI2","Ig")

# Initialize an empty results dataframe
results_df <- data.frame(
  Insula = character(),
  Observed_Correlation = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (region in combined_insula_regions) {
  left_insula  <- paste0("L_", region)
  right_insula <- paste0("R_", region)
  
  # Make a copy of the data for modification
  temp_df_L <- PCA_scores_SA_forviz_ordered
  temp_df_R <- PCA_scores_SA_forviz_ordered
  
  # Fill in missing Insula labels based on parcel names
  temp_df_L <- temp_df_L %>% mutate(Insula=replace_na(Insula, left_insula))
  temp_df_R <- temp_df_R %>% mutate(Insula=replace_na(Insula, right_insula))
  
  # Now index by these temporary labels
  left_data  <- temp_df_L[temp_df_L$Insula == left_insula, ]
  right_data <- temp_df_R[temp_df_R$Insula == right_insula, ]
  
  # Average PC1 across hemispheres
  mean_PC1 <- rowMeans(cbind(left_data$PC1, right_data$PC1), na.rm = TRUE)
  
  # Calculate the empirical correlation (left vs right data should have same axis)
  observed_correlation <- cor(mean_PC1, left_data$SA.axis, method = "spearman", use = "complete.obs")
  
  # Run the spin test using perm.sphere.p
  p_value <- perm.sphere.p(
    x = mean_PC1,
    y = left_data$SA.axis,
    perm.id = perm.id,
    corr.type = "spearman"
  )
  
  # Store the result in the results dataframe
  results_df <- rbind(
    results_df,
    data.frame(
      Insula = region,
      Observed_Correlation = observed_correlation,
      P_Value = p_value
    )
  )
}

# Display the final results
print(results_df)
results_df$Significant <- results_df$P_Value < 0.05
results_df$FDRp <- p.adjust(results_df$P_Value, method = "BH")

saveRDS(results_df,"insuladev/dropbox/SAaxis_spin_df_PC1.RDS")

# === repeat for 2nd PC (quadratic)
results_df2 <- data.frame(
  Insula = character(),
  Observed_Correlation = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (region in combined_insula_regions) {
  left_insula  <- paste0("L_", region)
  right_insula <- paste0("R_", region)
  
  # Make a copy of the data for modification
  temp_df_L <- PCA_scores_SA_forviz_ordered
  temp_df_R <- PCA_scores_SA_forviz_ordered
  
  # Fill in missing Insula labels based on parcel names
  temp_df_L <- temp_df_L %>% mutate(Insula=replace_na(Insula, left_insula))
  temp_df_R <- temp_df_R %>% mutate(Insula=replace_na(Insula, right_insula))
  
  # Now index by these temporary labels
  left_data  <- temp_df_L[temp_df_L$Insula == left_insula, ]
  right_data <- temp_df_R[temp_df_R$Insula == right_insula, ]
  
  # Average PC2 across hemispheres
  mean_PC2 <- rowMeans(cbind(left_data$PC2, right_data$PC2), na.rm = TRUE)
  
  # Calculate the empirical correlation (left vs right data should have same axis)
  observed_correlation <- cor(mean_PC2, left_data$SA.axis, method = "spearman", use = "complete.obs")
  
  # Run the spin test using perm.sphere.p
  p_value <- perm.sphere.p(
    x = mean_PC2,
    y = left_data$SA.axis,
    perm.id = perm.id,
    corr.type = "spearman"
  )
  
  # Store the result in the results dataframe
  results_df2 <- rbind(
    results_df2,
    data.frame(
      Insula = region,
      Observed_Correlation = observed_correlation,
      P_Value = p_value
    )
  )
}

results_df2$FDRp <- p.adjust(results_df2$P_Value, method = "BH")

saveRDS(results_df2,"insuladev/dropbox/SAaxis_spin_df_PC2.RDS")
# === altered perm.sphere.p function to use mean comparison (rank test) vs cor ===

perm.sphere.rank = function(x,group_labels,perm.id) {
  nroi = dim(perm.id)[1]  # number of regions
  nperm = dim(perm.id)[2] # number of permutations
  groups <- unique(group_labels)
  groups <- groups[!is.na(groups)]
  
  observed_stats <- tapply(x, group_labels, mean) #empiric correlation
  
  permuted_stats <- matrix(NA, nrow = length(groups), ncol = nperm)
  rownames(permuted_stats) <- groups
  
  # for each permutation, rearrange network assignment and recalculate mean
  for (r in 1:nperm) {
    x_perm <- x[perm.id[, r]]
    for (g_idx in seq_along(groups)) {
      g <- groups[g_idx]
      group_inds <- which(group_labels == g)
        permuted_stats[g_idx, r] <- mean(x_perm[group_inds], na.rm = TRUE)
    }
  }
  
  p_vals <- numeric(length(groups))
  names(p_vals) <- groups
  for (g in groups) {
    obs <- observed_stats[g]
    null_dist <- permuted_stats[g, ]
    obs_demean <- obs - mean(null_dist,na.rm = TRUE)
    null_demean <- null_dist - mean(null_dist,na.rm = TRUE)
    p_vals[g] <- mean(abs(null_demean) >= abs(obs_demean), na.rm = TRUE)
    
  }
  return(list(
    observed = observed_stats,
    null_distributions = permuted_stats,
    p_values = p_vals
  ))
}

# === networks vs PC1 for only insula poles (AAIC, AVI, Ig) ===

#Create bilateral merged dataframe
PCA_scores_combined <- PCA_scores_SA_forviz_ordered %>%
  mutate(Combined_Insula = str_remove(Insula, "^L_|^R_"))

#Group by the combined insula and GLASSER parcel and average
PCA_scores_averaged <- PCA_scores_combined %>%
  group_by(Combined_Insula, GLASSERLABELNAME.x) %>%
  summarise(Average_PC1 = mean(PC1, na.rm = TRUE), 
            NETWORK = first(NETWORK),
            KEYVALUE = first(KEYVALUE),
            SA.Axis = first(SA.axis),
            .groups = "drop")

PCA_scores_averaged <- do.call(rbind, lapply(
  split(PCA_scores_combined,
        list(PCA_scores_combined$Combined_Insula, PCA_scores_combined$GLASSERLABELNAME.x),
        drop = TRUE),
  function(x) {
    data.frame(
      Combined_Insula = x$Combined_Insula[1],
      GLASSER_ROI     = x$GLASSERLABELNAME.x[1],
      NETWORK         = x$NETWORK[1],
      KEYVALUE        = x$KEYVALUE[1],
      PC1_averaged    = mean(x$PC1, na.rm = TRUE),
      SA.Axis         = x$SA.axis[1]
    )
  }
))
row.names(PCA_scores_averaged) <- NULL

network_names <- unique(PCA_scores_combined$NETWORK) 

#with combined hemispheres
PCA_scores_averaged$NETWORK <- as.factor(PCA_scores_averaged$NETWORK)
PCA_scores_averaged$Combined_Insula <- as.factor(PCA_scores_averaged$Combined_Insula)
anova_result_interaction_combined <- aov(PC1_averaged ~ NETWORK * Combined_Insula, data = PCA_scores_averaged)
summary(anova_result_interaction_combined)

bilat_insula_regions_networks <- c("AAIC", "Ig", "AVI")  # Adjust as needed

network_spin_test_results_combined <- data.frame(
  Insula = character(),
  Network = character(),
  Observed_Avg_PC1 = numeric(),
  Spin_P_Value = numeric(),
  SD = numeric(),
  Mean = numeric(),
  stringsAsFactors = FALSE
)

network_spin_null_long <- list()

# Loop through each insula region
for (insula in bilat_insula_regions_networks) {
  # Filter data for the specific insula
  insula_data <- PCA_scores_averaged %>%
    filter(Combined_Insula == insula) %>%
    select(Combined_Insula, SA.Axis, KEYVALUE, NETWORK, PC1_averaged) %>%
    complete(KEYVALUE = 1:360) %>%         # adds missing KEYVALUE rows (insula rows are missing and will cause error)
    arrange(KEYVALUE)
  
  if (nrow(insula_data) != nrow(perm.id)) {
    stop(paste0("Row mismatch for ", insula, 
        ": ", nrow(insula_data), " parcels vs ", 
        nrow(perm.id), " in perm.id"))
  }
  
  result <- perm.sphere.rank(
    x = insula_data$PC1_averaged,
    group_labels = insula_data$NETWORK,
    perm.id = perm.id
  )
  
  # per-network null mean/sd, named by network
  null_mean <- rowMeans(result$null_distributions, na.rm = TRUE)
  null_sd   <- apply(result$null_distributions, 1, sd, na.rm = TRUE)

  networks <- names(result$observed)

  region_df <- data.frame(
    Insula = insula,
    Network = networks,
    Observed_Avg_PC1 = as.numeric(result$observed[networks]),
    Spin_P_Value     = as.numeric(result$p_values[networks]),
    SD               = as.numeric(null_sd[networks]),
    Mean             = as.numeric(null_mean[networks]),
    row.names = NULL
  )
  network_spin_test_results_combined <- rbind(network_spin_test_results_combined, region_df)

  # Long-format null distribution for each Insula × Network (one row per permutation)
  mat_null <- result$null_distributions
  netw_rows <- rownames(mat_null)
  null_chunk <- dplyr::bind_rows(lapply(seq_len(nrow(mat_null)), function(k) {
    tibble::tibble(
      Insula = insula,
      Network = netw_rows[k],
      null_value = as.numeric(mat_null[k, ])
    )
  }))
  network_spin_null_long[[length(network_spin_null_long) + 1L]] <- null_chunk
}

print(network_spin_test_results_combined)

network_spin_test_results_combined <- network_spin_test_results_combined %>%
  mutate(FDR_all = p.adjust(Spin_P_Value, method = "BH")) 

saveRDS(network_spin_test_results_combined, "insuladev/dropbox/network_spin_results.RDS")

network_spin_null_long <- dplyr::bind_rows(network_spin_null_long)
saveRDS(network_spin_null_long, "insuladev/dropbox/network_spin_null_long.RDS")
# === PVD Second PCA === 
# == PC1

PC1 <- data.frame(PC1 = Glasser_pca_scores$PC1, Brain = Glasser_pca_scores$Brain, Insula = Glasser_pca_scores$Insula)
PC1_wide <- pivot_wider(PC1, names_from = Insula, values_from = PC1)
PC1_wide <- column_to_rownames(PC1_wide, var = "Brain")

Glasser_pca_squared <- prcomp(PC1_wide, scale. = FALSE)
factoextra::fviz_eig(Glasser_pca_squared)
Glasser_PC1_loadings <- as.data.frame(Glasser_pca_squared$rotation)

# == PC2 

PC2 <- data.frame(PC2 = Glasser_pca_scores$PC2, Brain = Glasser_pca_scores$Brain, Insula = Glasser_pca_scores$Insula)
PC2_wide <- pivot_wider(PC2, names_from = Insula, values_from = PC2)
PC2_wide <- column_to_rownames(PC2_wide, var = "Brain")

Glasser_pca2_squared <- prcomp(PC2_wide, scale. = FALSE)
factoextra::fviz_eig(Glasser_pca2_squared)
Glasser_PC2_loadings <- as.data.frame(Glasser_pca2_squared$rotation)

Glasser_PCAofPC1_scores <-Glasser_pca_squared$x
Glasser_PCAofPC2_scores <-Glasser_pca2_squared$x

# === Save Outputs ===
saveRDS(Glasser_PCAofPC1_scores, "insuladev/dropbox/PCAofPC1scores.RDS")
saveRDS(Glasser_PCAofPC2_scores, "insuladev/dropbox/PCAofPC2scores.RDS")
saveRDS(Glasser_PC1_loadings, "insuladev/dropbox/PCAofPC1_loadings.RDS")
saveRDS(Glasser_PC2_loadings, "insuladev/dropbox/PCAofPC2_loadings.RDS")

# == spin test for 2nd PC loadings == 
Glasser_PC1_loadings$Brain <- paste0(row.names(Glasser_PC1_loadings), "_ROI")
Glasser_PC2_loadings$Brain <- paste0(row.names(Glasser_PC2_loadings), "_ROI")

PCA_SA_1 <- merge(S.A.axis, Glasser_PC1_loadings, by.x = "orig_parcelname", by.y = "Brain", all.x = TRUE, all.y = TRUE)
PCA_SA_2 <- merge(S.A.axis, Glasser_PC2_loadings, by.x = "orig_parcelname", by.y = "Brain", all.x = TRUE, all.y = TRUE)

PCA_SA_1 <- PCA_SA_1[order(PCA_SA_1$SA.axis), ]
PCA_SA_2 <- PCA_SA_2[order(PCA_SA_2$SA.axis), ]

PCA_SA_1$GLASSERLABELNAME <- PCA_SA_1$orig_parcelname
PCA_SA_2$GLASSERLABELNAME <- PCA_SA_2$orig_parcelname


run_spin_test <- function(x, y, perm.id) {
  
  # observed correlation
  obs <- cor(x, y, method = "spearman", use = "complete.obs")
  
  # null distribution
  null_vals <- numeric(ncol(perm.id))
  
  for (i in seq_len(ncol(perm.id))) {
    x_perm <- x[perm.id[, i]]
    null_vals[i] <- cor(x_perm, y, method = "spearman", use = "complete.obs")
  }
  
  # corrected p-value
  k <- sum(abs(null_vals) >= abs(obs), na.rm = TRUE)
  p_val <- (k + 1) / (length(null_vals) + 1)
  
  return(list(
    observed = obs,
    p_value = p_val,
    null = null_vals
  ))
}

pc1_result <- run_spin_test(
  PCA_SA_1$PC1,
  PCA_SA_1$SA.axis,
  perm.id
)

pc2_result <- run_spin_test(
  PCA_SA_2$PC1,
  PCA_SA_2$SA.axis,
  perm.id
)

saveRDS(pc1_result, "insuladev/dropbox/SAaxis_spin_2ndPCA.RDS")
