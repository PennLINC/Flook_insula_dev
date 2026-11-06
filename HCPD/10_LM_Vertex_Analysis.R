#LM Analysis, Vertex Level
#run locally
# Purpose: Extract, normalize, and analyze vertex-wise linear model coefficients from voxel-wise .RDS outputs

# === Load Required Libraries ===
library(ciftiTools)
library(readxl)
library(tibble)
library(factoextra)

ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')

# === Load Model Outputs ===
Polynomial_demeans <- list.files(path = "insuladev/HCPDbold_1/Polynomial_Results",  full.names = TRUE, 
                                 pattern = "LM_demeaned_results_voxel\\d+\\.RDS", 
                                 recursive = TRUE)
Age_est_de <- list() #linear age term
Age_sq_est_de <- list() #quadratic age term

# === Extract coefficients from each voxel ===
# Function to extract linear age estimates
voxel_age_PCAs_de <- function(file) {
  data <- readRDS(file)
  age_estimate <- data$age_demeaned.estimate 
  voxel_id <- gsub(".*(voxel\\d+).*", "\\1", basename(file))
  Age_est_de[[voxel_id]] <- age_estimate
}

# Function to extract quadratic age^2 estimates
voxel_agesq_PCAs_de <- function(file) {
  data <- readRDS(file)
  agesq_estimate <- data$`I(age_demeaned^2).estimate`
  voxel_id <- gsub(".*(voxel\\d+).*", "\\1", basename(file))
  Age_sq_est_de[[voxel_id]] <- agesq_estimate
}

# === Apply extraction ===
age_PCA_df_de <- lapply(Polynomial_demeans, voxel_age_PCAs_de)
age_sq_PCA_df_de <- lapply(Polynomial_demeans, voxel_agesq_PCAs_de)

# === Assemble DataFrames ===
Age_est_de <- as.data.frame(age_PCA_df_de)
Age_sq_est_de <- as.data.frame(age_sq_PCA_df_de)

voxel_ids_de <- gsub(".*(voxel\\d+).*", "\\1", basename(Polynomial_demeans)) # Extract vertex IDs from file names

# Assign voxel IDs as column names
colnames(Age_est_de) <- paste(voxel_ids_de)
colnames(Age_sq_est_de) <- paste(voxel_ids_de)

# Create new row names
new_rownames_de <- paste0(1:nrow(Age_est_de), "_age")
new_rownames_sq_de <- paste0(1:nrow(Age_sq_est_de), "_agesq")
rownames(Age_est_de) <- new_rownames_de
rownames(Age_sq_est_de) <- new_rownames_sq_de

# === Save raw coefficient matrices ===
write.csv(as.data.frame(Age_est_de), "insuladev/dropbox/demeaned_age_coefficients.csv")
write.csv(as.data.frame(Age_sq_est_de), "insuladev/dropbox/demeaned_age_sq_coefficients.csv")

# === Create combined (linear and quadratic terms) matrix and remove NA ===
combined_Age_coef_de <- cbind(t(Age_est_de),t(Age_sq_est_de)) #create dataframe with insula as rows and linear + quad coefficients as columns

na_cols_de <- apply(combined_Age_coef_de, 2, function(x) all(is.na(x))) #remove columns that are all NA
combined_Age_coef_de <- combined_Age_coef_de[, !na_cols_de]

#check to ensure that all remaining rows only have two NA values
all(rowSums(is.na(combined_Age_coef_de)) == 2) #if TRUE, then all remaining NAs refer to insula columns (NA because self-referential to rows)
#remove all remaining columns with NAs
combined_Age_coef_de_final <- combined_Age_coef_de[, colSums(is.na(combined_Age_coef_de)) == 0]

# === Normalize columns ===
#separate back into linear vs quad terms to normalize, then recombine
age_cols <- length(grep("_age$", colnames(combined_Age_coef_de_final), value = TRUE)) #determine number of cols
sd_age <- sd(as.matrix(combined_Age_coef_de_final[, 1:age_cols])) # Compute standard deviation for age columns
sd_agesq <- sd(as.matrix(combined_Age_coef_de_final[, (age_cols + 1):(2 * age_cols)])) # Compute standard deviation for agesq column

combined_Age_demeaned <- combined_Age_coef_de_final # Create a copy of the original dataframe
combined_Age_demeaned [, 1:age_cols] <- combined_Age_demeaned [, 1:age_cols] / sd_age 
combined_Age_demeaned[, (age_cols+1):(2*age_cols)] <- combined_Age_demeaned[, (age_cols+1):(2*age_cols)] / sd_agesq

# === Inspect age and age^2 collinearity ===
diag_corr <- numeric(age_cols) # create a vector for storing the diagonal correlations

for (i in 1:age_cols) {
  diag_corr[i] <- cor(combined_Age_demeaned[, i],
                      combined_Age_demeaned[, i + age_cols],
                      use = "pairwise.complete.obs")
}

mean_diag_corr <- mean(diag_corr)
print(mean_diag_corr) #-0.19, aka not very correlated which is reassuring

# === PCA on normalized coefficients ===
PCA_comb_age_normalized <- prcomp(combined_Age_demeaned, center = TRUE, scale = FALSE)
factoextra::fviz_eig(PCA_comb_age_normalized)

#scores
PC1_scores_norm <- as.data.frame(PCA_comb_age_normalized$x[,1])
PC1_scores_norm$Insula <- as.numeric(gsub("voxel", "", rownames(PC1_scores_norm)))
PC1_scores_norm <- PC1_scores_norm[order(PC1_scores_norm$Insula),]

# === Save output ===
saveRDS(PC1_scores_norm, "insuladev/dropbox/PC1_scores_vertex.RDS")

#loadings
PC1_loadings_norm <- as.data.frame(PCA_comb_age_normalized$rotation[,1])

PC1_split_norm <- do.call(rbind, strsplit(rownames(PC1_loadings_norm), "_"))
PC1_split_norm <- as.data.frame(PC1_split_norm) # Convert to a dataframe and assign column names
colnames(PC1_split_norm) <- c("Voxel", "Term")  # Rename columns appropriately

PC1_loadings_norm <- cbind(PC1_split_norm, PC1_loadings_norm)
PC1_loadings_norm$Voxel <- as.numeric(PC1_loadings_norm$Voxel)

PC1_loadings_age_norm <- subset(PC1_loadings_norm, Term == "age")
PC1_loadings_age_sq_norm <- subset(PC1_loadings_norm, Term == "agesq")

# === Save output ===
saveRDS(PC1_loadings_age_norm, "insuladev/dropbox/PC1_loadings_age_vertex.RDS")
saveRDS(PC1_loadings_age_sq_norm, "insuladev/dropbox/PC1_loadings_age_sq_vertex.RDS")

