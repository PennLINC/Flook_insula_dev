### Sample Selections including demographics and QC for HCPD
#demographs, rs_files, and qc_files downloaded from HCP-D

################################################
############ STEPS INCLUDED ####################
################################################
# 1) import HCPD demographics
# 2) import RS scans and filter demographics by subs with rs scans
# 3) import medical info and filter out subs with med exclusion
# 4) import qc files and exclude subs with fd means >/=2 #
# 5) exclude participants with less than 6 minutes of rs scans #
# 6) create final files for correlation matrix creation

library(readr)
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)


# 1) Import HCPD demographics #

hcpd_participants <- read_delim("insuladev/dropbox/hcpd_participants.tsv", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)

## total participants 652
sum(is.na(hcpd_participants$age)) # Number of missing values, 0

####summary table
mean(hcpd_participants$age) #mean=14.4
sd(hcpd_participants$age)# SD=4.1
length(which(hcpd_participants$sex=="Female")) #Female n = 351


# 2) import rs scans and filter by subs with rs scans #

rs_subjects_files <- list.files(path = "insuladev/HCPDbold_1/rs_files", 
            full.names = TRUE, pattern = "denoised_bold.dtseries.nii$", recursive = TRUE) %>%
  rbind.fill()

rs_subjects <- data.frame(do.call(rbind, strsplit(rs_subjects_files, 
                    "insuladev/HCPDbold_1/rs_files/", fixed=TRUE)))
rs_subs <- separate_wider_delim(rs_subjects, cols = "X2", delim = c("/"), 
                    names = c("sub", "session"))

HCPD_rs_participants <- hcpd_participants %>% filter(participant_id %in% unique(rs_subs$sub))

#### participants with rs scans: 631 ####
mean(HCPD_rs_participants$age) #mean=14.6
sd(HCPD_rs_participants$age)# SD=3.9
length(which(HCPD_rs_participants$sex=="Female")) #Female n = 338

# 3) import medical info and filter out subs with med exclusion #

### filter out medical exclusions ###
medh01 <- read_delim("insuladev/dropbox/medh01.txt", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

# cancer/leukemia
cancer <- medh01$src_subject_id[c(which(medh01$medhis_2e == 1))] # HCD0641341

# lead poisoning
lead<- medh01$src_subject_id[c(which(medh01$medhis_2k==1))] # "HCD0392144" "HCD0641341"

# sickle cell anemia
sickle <- medh01$src_subject_id[c(which(medh01$medhis_2p==1))] # HCD0641341

# accidental poisoning
poisoning <- medh01$src_subject_id[c(which(medh01$medhis_6q ==1))] # "HCD0529751" "HCD2855875"

# multiple sclerosis
multiple_sclerosis <- medh01$src_subject_id[c(which(medh01$ms ==1))]  # "HCD0641341" "HCD1617450"

# seizure
seizure_ph9 <- medh01$src_subject_id[c(which(medh01$ph_9 ==1))]  # "HCD0360838" "HCD0478356" - each participant had one seizure
epilepsy_cfmh <- medh01$src_subject_id[c(which(medh01$cfmh_chd_seizure ==1))] # "HCD0641341" "HCD1277149"  

# brain injury
brain_injury <- medh01$src_subject_id[c(which(medh01$seq1c_2 ==1))] # 14 participants 

medical_exclusion <- unique(c(cancer, lead, sickle, poisoning, 
                              multiple_sclerosis, seizure_ph9, epilepsy_cfmh, brain_injury))
medical_exclusion_sub <- sub("HCD", "sub-", 
                             medical_exclusion, fixed = TRUE) ##21 participants

#demographics following med exclusion
HCPD_rs_medexclusion <- HCPD_rs_participants[-c(which(HCPD_rs_participants$participant_id %in% medical_exclusion_sub)),]
### n with resting state and after filtering out med exclusions: 610 participants
mean(HCPD_rs_medexclusion$age) #mean=14.7
sd(HCPD_rs_medexclusion$age)# SD=3.9
length(which(HCPD_rs_medexclusion$sex=="Female")) #Female n = 327


# 4) import qc files and exclude subs with fd means >/=.2 #

###filter by QC measures
qc_FDmeans <- list.files(path = "insuladev/HCPDbold_1/qc_files/", 
                         full.names = TRUE, pattern = "qc.csv$", recursive = TRUE) %>%
  lapply(read_csv) %>%
  rbind.fill()

qc_FDmean_2 <- filter(qc_FDmeans, mean_fd < 0.2)

summarise(qc_FDmean_2, length(unique(unlist(sub))))
###any subjects with 1+ run with FD mean <0.2 = 583

###combine demographics with subjects included by qc filter
HCPD_rs_medexclusion$participant_id <- sub("sub-", "", HCPD_rs_medexclusion$participant_id)
qc_participants <- qc_FDmean_2 %>% 
  group_by(sub) %>% 
  slice (1)
HCPD_rs_med_qc <- merge(HCPD_rs_medexclusion, qc_participants, by.x = "participant_id", by.y = "sub")

###participants after filter for head motion >0.2, rs scan, med exclusion: 563


# 5) exclude participants with less than 6 minutes of rs scans #

#filter subjects by total scan time after QC
#except for 0577358 and 1529756, all remaining subjects have 478 volumes in each rs scan with TR 0.8 (6.37min per scan)
#0577358 only has 3.37 minutes of scanning after QC (only 1 shortened scan) so will be excluded
#1529756 has 4 shortened scans after QC totaling 13.49 minutes, so will be included

rs_subs_keep_files <- data.frame(qc_FDmean_2$sub, qc_FDmean_2$run, qc_FDmean_2$dir)
names(rs_subs_keep_files)[1:3] <- c("Subs", "Run", "Dir")
rs_subs_keep_files$filenames <- paste0("sub-", 
              rs_subs_keep_files$Subs, "_ses-V1_task-rest_dir-", 
              rs_subs_keep_files$Dir,"_run-", rs_subs_keep_files$Run, 
              "_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii")
rs_subs_rm_shortscans <- rs_subs_keep_files[!(rs_subs_keep_files$Subs == "0577358"),]


# 6) final file list for connectivity matrices and sub list for demographics #

rs_subs_final_files <- rs_subs_rm_shortscans %>% 
  filter(Subs %in% unique(HCPD_rs_med_qc$participant_id))

#remaining scans: 1915 from 562 participants
table(table(rs_subs_final_files$Subs))
#377 subjects have 4 rs scans, 86 subjects have 3 rs scans, 50 subjects have 2 rs scans, and 49 subjects have 1 rs scan

###final participants (demographics)
HCPD_final_demo <- HCPD_rs_med_qc[!(HCPD_rs_med_qc$participant_id == "0577358"),1:14]
###total subjects: 562
mean(HCPD_final_demo$age) #mean=14.8, 
sd(HCPD_final_demo$age) #SD=3.9
length(which(HCPD_final_demo$sex=="Female")) #Female n = 300


saveRDS(HCPD_final_demo, file = "insuladev/dropbox/HCPD_demographics.RDS")
saveRDS(rs_subs_final_files, file = "insuladev/dropbox/HCPD_rs_file_names.RDS")

#make file list for wb_command command line loop
rs_file_list <- paste0("sub-", rs_subs_final_files$Subs, "/", rs_subs_final_files$filenames)
write.table(rs_file_list, file = "insuladev/dropbox/HCPD_rs_file_list.csv", sep = "\t", row.names=FALSE, col.names=FALSE)
