### Creating files for model array to run

library(plyr)
library(tidyverse)
library(dplyr)
library(readr)

# Load subject data from a single file
subject_data <- readRDS("insuladev/dropbox/HCPD_demographics.RDS")

#Load FD_Means of scans into object qc_FDmeans
FD_Means_List <- readRDS("insuladev/dropbox/HCPD_rs_file_names.RDS")
FD_Means_Files <- data.frame(matrix(nrow = 0, ncol = 1)) 
colnames(FD_Means_Files) <- "Files"
for(i in 1:nrow(FD_Means_List)) {
  FD_Means_Files[i,] <- paste0("insuladev/HCPDbold_1/qc_files/sub-",FD_Means_List$Subs[i],"/sub-",FD_Means_List$Subs[i],"_ses-V1_task-rest_dir-", FD_Means_List$Dir[i],"_run-", FD_Means_List$Run[i],"_space-fsLR_den-91k_desc-linc_qc.csv")
}
FD_Means_Files$Files <- gsub("sub-1529756_ses-V1_task-rest_","sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_", FD_Means_Files$Files)
qc_FDmeans <- lapply(FD_Means_Files$Files, function(file) read_csv(file, show_col_types = FALSE)) %>%
  rbind.fill()
print("Created qc_FDmeans")

#average FD means for each participant
qc_average <- aggregate(qc_FDmeans$mean_fd, list(qc_FDmeans$sub), mean)
print("Created qc_average")

#create single file with demographic data and fd data 
subject_data_qc <- merge(subject_data, qc_average, by.x = "participant_id", by.y = "Group.1")
colnames(subject_data_qc)[colnames(subject_data_qc) == 'x'] <- 'mean_fd'
print("Created combined demo + qc file")

#create conn matrix filename list
###need to update depending on where running from#####
conmat_files <- list.files(
  path = "insuladev/HCPDbold_1/cor_results",
  pattern = "^correlation_\\d+\\.h5$"
)

file_df <- data.frame(
  participant_id = sub("^0+", "", sub("^correlation_(\\d+)\\.h5$", "\\1", conmat_files)),
  source_file = file.path("cor_results", conmat_files),
  stringsAsFactors = FALSE
)

subject_data_qc$participant_id <- sub("^0+", "", as.character(subject_data_qc$participant_id))
model_array_demo <- merge(subject_data_qc, file_df, by = "participant_id", sort = FALSE)
model_array_demo <- model_array_demo[order(model_array_demo$source_file), ]
model_array_demo$scalar_name <- "Pearson_R"

write.csv(model_array_demo, "insuladev/HCPDbold_1/cor_results/model_array_demographics.csv", row.names = FALSE)