#ROI visualization 
#Run after ROI_GAM_Setup.R (via ex.sh SLURM script)

library(plyr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(readr)
library(oro.nifti)
library(colorspace)
library(ciftiTools)

ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')

# === input data === 
ROI_connectivity <- read.csv("insuladev/ROI_GAM/ROI_single_line.csv")

# === HCP-D Averaged connectivity === 
average_connectivity <- as.data.frame(colMeans(ROI_connectivity[,!names(ROI_connectivity) %in% c("subject")]))

average_connectivity$Names_to_filter <- rownames(average_connectivity)
average_connectivity_renamed <- average_connectivity %>% 
  separate_wider_delim(
    Names_to_filter, 
    delim = "_", 
    names = c("Insula_Hemisphere", "Insula_Subregion", "Glasser_ROI"), 
    too_many = "merge")
colnames(average_connectivity_renamed) <- c("connectivity","Insula_Hemisphere", "Insula_Subregion", "Glasser_ROI")
average_connectivity_renamed$Glasser_ROI <- gsub("L_a9.46v_ROI","L_a9-46v_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)
average_connectivity_renamed$Glasser_ROI <- gsub("L_OP2.3_ROI","L_OP2-3_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)
average_connectivity_renamed$Glasser_ROI <- gsub("L_9.46d_ROI","L_9-46d_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)
average_connectivity_renamed$Glasser_ROI <- gsub("L_s6.8_ROI","L_s6-8_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)
average_connectivity_renamed$Glasser_ROI <- gsub("L_i6.8_ROI","L_i6-8_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)
average_connectivity_renamed$Glasser_ROI <- gsub("L_p9.46v_ROI","L_p9-46v_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)
average_connectivity_renamed$Glasser_ROI <- gsub("R_a9.46v_ROI","R_a9-46v_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)
average_connectivity_renamed$Glasser_ROI <- gsub("R_OP2.3_ROI","R_OP2-3_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)
average_connectivity_renamed$Glasser_ROI <- gsub("R_9.46d_ROI","R_9-46d_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)
average_connectivity_renamed$Glasser_ROI <- gsub("R_s6.8_ROI","R_s6-8_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)
average_connectivity_renamed$Glasser_ROI <- gsub("R_i6.8_ROI","R_i6-8_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)
average_connectivity_renamed$Glasser_ROI <- gsub("R_p9.46v_ROI","R_p9-46v_ROI", average_connectivity_renamed$Glasser_ROI, fixed = TRUE)

rois <- c("L_PoI1", "R_PoI1", "L_PoI2", "R_PoI2", 
          "L_Ig", "R_Ig", "L_AAIC", "R_AAIC", 
          "L_AVI", "R_AVI", "L_MI", "R_MI")
rois_with_suffix <- paste0(rois, "_ROI")

final_average <- average_connectivity_renamed[!(average_connectivity_renamed$Glasser_ROI %in% rois_with_suffix | average_connectivity_renamed$Glasser_ROI == "???"| average_connectivity_renamed$Glasser_ROI == "..."), ]
write_csv(final_average, "insuladev/dropbox/HCPD_avg_con.csv")