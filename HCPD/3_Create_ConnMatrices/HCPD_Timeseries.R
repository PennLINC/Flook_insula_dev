### Creating connectivity matrices for insula vertices - whole brain

library(cifti) 
library(dplyr)
library(gifti)
library(stringr)
library(data.table)

#1) concatenate all usable rs files for each participant (one for insula and one for whole brain)
#2) create connectivity matrices between insula vertices and whole brain vertices

############################
# get subjects and filepaths
############################
rs_subs_final_files <- readRDS("insuladev/dropbox/HCPD_rs_file_names.RDS")
#add acqvariant to filename for subject 1529756
rs_subs_final_files$filenames[grep("1529756",rs_subs_final_files$filenames)] <- sub("_dir", "_acq-VARIANTNumVolumes_dir", rs_subs_final_files$filenames[grep("1529756", rs_subs_final_files$filenames)])

HCPD_filenames  <- as.list(rs_subs_final_files$filenames)
participants <- unique(rs_subs_final_files$Subs)

############################
# concatenate
############################

#function for whole brain files
timeseries_files <- function(subject_id) {
 
  # Filter files for the specific subject
  subject_files <- rs_subs_final_files %>% 
    filter(Subs == subject_id)
  
  # Read and cbind files for the subject
  combined_data <- lapply(subject_files$filenames, function(file) {
    # Read each file
    cifti_data <- read_cifti(paste0("insuladev/HCPDbold_1/rs_files/sub-", subject_id, "/", file))
    timeseries <- cifti_data$data
    return(timeseries)
  })
  combined_subject_data <- do.call(cbind, combined_data)
  
  return(combined_subject_data)
}

#function for insula-specific files
timeseries_insula_files <- function(subject_id) {
  
  # Filter files for the specific subject
  subject_files <- rs_subs_final_files %>% 
    filter(Subs == subject_id)
  
  # Read and cbind files for the subject
  combined_insula_data <- lapply(subject_files$filenames, function(file) {
    # Read each file
    insula_right_filename <- paste0("insuladev/HCPDbold_1/rs_files/sub-", subject_id, "/", sub('\\.dtseries.nii$', '', file), "_insula_Combined_Right.dtseries.nii")
    insula_cifti_right <- read_cifti(insula_right_filename)
    insula_left_filename <- paste0("insuladev/HCPDbold_1/rs_files/sub-", subject_id, "/", sub('\\.dtseries.nii$', '', file), "_insula_Combined_Left.dtseries.nii")
    insula_cifti_left <- read_cifti(insula_left_filename)
    #extract timeseries for left and right insula and combine
    timeseries_insula_left <- insula_cifti_left$data
    timeseries_insula_right <- insula_cifti_right$data
    bilateral_insula_ts <- rbind(timeseries_insula_left, timeseries_insula_right)
    return(bilateral_insula_ts)
  })
  combined_subject_insula_data <- do.call(cbind, combined_insula_data)
  
  return(combined_subject_insula_data)
}

#correlate insula and whole brain (transpose, otherwise will error due to different shape of files)
dir.create("insuladev/HCPDbold_1/cor_results", recursive = TRUE, showWarnings = FALSE)
lapply(participants, function(subject_id){
  print(paste0(subject_id, "running conn matrices"))
  if (! file.exists(paste0("insuladev/HCPDbold_1/cor_results/correlation_", subject_id, ".csv"))){
  combined_data_temp <- (timeseries_files(subject_id)) #run whole brain function
  combined_data_t <- t(combined_data_temp)
  combined_insula_data_temp <- timeseries_insula_files(subject_id) #run insula function
  combined_insula_data_t <- t(combined_insula_data_temp)
  cor_result <- round(cor(combined_data_t,combined_insula_data_t,use="complete.obs"),5)
  write.csv(cor_result, file = paste0("insuladev/HCPDbold_1/cor_results/correlation_", subject_id, ".csv"), row.names = FALSE) }
  print(paste0(subject_id, " correlation matrix DONE"))
})