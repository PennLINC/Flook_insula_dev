#Compare Vertex findings to 

library(gifti)
library(ciftiTools)
library(readxl)
ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')

PC1_loadings_age <- readRDS("insuladev/dropbox/PC1_loadings_age_vertex.RDS")
PC1_loadings_age_sq <- readRDS("insuladev/dropbox/PC1_loadings_age_sq_vertex.RDS")
PC1_scores_age <- readRDS('insuladev/dropbox/PC1_scores_vertex.RDS')

Insula_Labels <- read_excel("insuladev/dropbox/Labels.xlsx")
Vertex_Labels <- read_excel("insuladev/dropbox/Vertex_Labels.xlsx")
Vertex_Assignments <- read_excel("insuladev/dropbox/Vertex_Assignments.xlsx")
Cifti_series <- ciftiTools::read_cifti("insuladev/dropbox/fordata.dtseries.nii")

#hist_data is in space with medial wall, so need to add it back in to convert
mw_mask_L <- Cifti_series$meta$cortex$medial_wall_mask$left
mw_mask_R <- Cifti_series$meta$cortex$medial_wall_mask$right

PC1_vis_scores_vertex <- merge(Insula_Labels, PC1_scores_age, by.x = "ColumnNumber", by.y = "Insula", all.y = TRUE)
PC1_visualization_scores_vertex <- merge(PC1_vis_scores_vertex, Vertex_Assignments, by = "Vertex", all.y = TRUE)
PC1_brain_scores_vertex <- PC1_visualization_scores_vertex[order(PC1_visualization_scores_vertex$Vertex), ]

brain_data_scores_vertex <- newdata_xifti(Cifti_series, PC1_brain_scores_vertex$`PCA_comb_age_normalized$x[, 1]` , newnames = NULL) #replace cifti_series data with new values

full_scores_L <- rep(NA, 32492)
full_scores_R <- rep(NA, 32492)

Left_vertex_labels <- Vertex_Labels[1:29696,] 
Right_vertex_labels <- Vertex_Labels[29697:59412,] 
Labels_with_mw_L <- rep(NA, 32492)
Labels_with_mw_R <- rep(NA, 32492)
Labels_with_mw_L[mw_mask_L] <- Left_vertex_labels$Label
Labels_with_mw_R[mw_mask_R] <- Right_vertex_labels$Label

# Define the label names of interest
target_labels <- c("L_AAIC_ROI", "L_AVI_ROI", "L_MI_ROI", 
                   "L_PoI2_ROI", "L_PoI1_ROI", "L_Ig_ROI")
target_labels_R <- c("R_AAIC_ROI", "R_AVI_ROI", "R_MI_ROI", 
                   "R_PoI2_ROI", "R_PoI1_ROI", "R_Ig_ROI")

# Get matching row indices
target_indices_L <- which(Labels_with_mw_L %in% target_labels)

# Create data frame
label_df_L <- data.frame(
  Vertex = target_indices_L,
  Label = Labels_with_mw_L[target_indices_L]
)

# Get matching row indices
target_indices_R <- which(Labels_with_mw_R %in% target_labels_R)

# Create data frame
label_df_R <- data.frame(
  Vertex = target_indices_R,
  Label = Labels_with_mw_R[target_indices_R]
)

# Insert values at non-medial wall positions
full_scores_L[mw_mask_L] <- brain_data_scores_vertex$data$cortex_left
full_scores_R[mw_mask_R] <- brain_data_scores_vertex$data$cortex_right


#compare to micro
MicroG1_L <- read_gifti("insuladev/HCPDbold_1/BigBrain/tpl-fs_LR_hemi-L_den-32k_desc-Micro_G1.shape.gii")
MicroG1_R <- read_gifti("insuladev/HCPDbold_1/BigBrain/tpl-fs_LR_hemi-R_den-32k_desc-Micro_G1.shape.gii")
MicroG2_L <- read_gifti("insuladev/HCPDbold_1/BigBrain/tpl-fs_LR_hemi-L_den-32k_desc-Micro_G2.shape.gii")
MicroG2_R <- read_gifti("insuladev/HCPDbold_1/BigBrain/tpl-fs_LR_hemi-R_den-32k_desc-Micro_G2.shape.gii")

cor(full_scores_L, MicroG2_L$data$shape, use="complete.obs",method="spearman") 
cor(full_scores_R, MicroG2_R$data$shape, use="complete.obs",method="spearman") 
cor(full_scores_L, MicroG1_L$data$shape, use="complete.obs",method="spearman") 
cor(full_scores_R, MicroG1_R$data$shape, use="complete.obs",method="spearman") 

#statistical test results (Moran run in python)
#Micro G1 and PC1 — left: rho = 0.604, p_moran = 0.015; right: rho = 0.655, p_moran = 0.003
#Micro G2 and PC1 — left: rho = 0.476, p_moran = 0.014; right: rho = 0.543, p_moran = 0.035

cc_L_M1 <- complete.cases(full_scores_L, MicroG1_L$data$shape)
df_MicroG1_L <- data.frame(
  MicroG1 = MicroG1_L$data$shape[cc_L_M1],
  PC1 = full_scores_L[cc_L_M1]
)

cc_R_M1 <- complete.cases(full_scores_R, MicroG1_R$data$shape)
df_MicroG1_R <- data.frame(
  MicroG1 = MicroG1_R$data$shape[cc_R_M1],
  PC1 = full_scores_R[cc_R_M1]
)

cc_L_M2 <- complete.cases(full_scores_L, MicroG2_L$data$shape)
df_MicroG2_L <- data.frame(
  MicroG2 = MicroG2_L$data$shape[cc_L_M2],
  PC1 = full_scores_L[cc_L_M2]
)

cc_R_M2 <- complete.cases(full_scores_R, MicroG2_R$data$shape)
df_MicroG2_R <- data.frame(
  MicroG2 = MicroG2_R$data$shape[cc_R_M2],
  PC1 = full_scores_R[cc_R_M2]
)

saveRDS(df_MicroG1_R, "insuladev/dropbox/microG1_R.RDS")
saveRDS(df_MicroG2_R, "insuladev/dropbox/microG2_R.RDS")
saveRDS(df_MicroG1_L, "insuladev/dropbox/microG1_L.RDS")
saveRDS(df_MicroG2_L, "insuladev/dropbox/microG2_L.RDS")

MicroG1_L_masked <- MicroG1_L$data$shape[mw_mask_L]
MicroG1_R_masked <- MicroG1_R$data$shape[mw_mask_R]
MicroG2_L_masked <- MicroG2_L$data$shape[mw_mask_L]
MicroG2_R_masked <- MicroG2_R$data$shape[mw_mask_R]


Viewing_Micro <- data.frame("Vertex" = 1:59412, "MicroG1" = c(MicroG1_L_masked,MicroG1_R_masked), "MicroG2" = c(MicroG2_L_masked,MicroG2_R_masked))
Micro_vertex <- merge(Insula_Labels, Viewing_Micro, by.x = "Vertex", all.x = TRUE) #only keep insula-vertices
Micro_vis_vertex <- merge(Micro_vertex, Vertex_Assignments, by = "Vertex", all.y = TRUE)
Micro_all_vertex <- merge(Viewing_Micro, Vertex_Assignments, by = "Vertex", all.y = TRUE)
Micro_all_vertex_order <- Micro_all_vertex[order(Micro_all_vertex$Vertex),]
Micro_brain_vertex <- Micro_vis_vertex[order(Micro_vis_vertex$Vertex), ]

## ------- AFTER RUNNING PYTHON SCRIPT (big_brain_moran) WITH MORAN TO GET NULLS ------------
#convert python objects back into R

#download python results to get the null distributions
 Micro_lh <- reticulate::py_load_object("insuladev/dropbox/Micro_G1_lh_moran.pkl")
 Micro_rh <- reticulate::py_load_object("insuladev/dropbox/Micro_G1_rh_moran.pkl")
        
 saveRDS(list(observed = Micro_lh$observed, p_value = Micro_lh$p_value, null = unlist(Micro_lh$null)), 
                "insuladev/dropbox/microG1_L_moran.RDS")
 saveRDS(list(observed = Micro_rh$observed, p_value = Micro_rh$p_value, null = unlist(Micro_rh$null)), 
                "insuladev/dropbox/microG1_R_moran.RDS")