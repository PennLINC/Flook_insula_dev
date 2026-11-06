###Subcortical_GAMs
#combine all matrices into one file and run GAM on insula-subcortex connectivity
library(dplyr)
library(readr)
library(mgcv)
library(readxl)
library(ggplot2)

average_by_label <- function(mat) {
  X <- mat
  lbl <- X[, 1]
  
  # drop label col 
  X <- X[, -1, drop = FALSE]
  X <- apply(X, 2, as.numeric) #ensure that reading connectivity values as numeric
  X <- as.matrix(X)
  
  #sum values for each label (ignore NA)
  summed <- rowsum(X, group = lbl, reorder = TRUE, na.rm = TRUE) 
  #count number of non-na cells
  counts <- rowsum(1L * !is.na(X), group = lbl, reorder = TRUE) 
  
  means <- summed / counts #compute mean
  means[is.nan(means)] <- NA #convert NANs (created if there is a row of NA) to NA
  
  out <- data.frame(Label = rownames(means), means, check.names = FALSE)
  rownames(out) <- NULL
  out
}
filtered_matrices1 <- readRDS("insuladev/dropbox/subcortex_matrices.RDS")
averaged_matrices1 <- lapply(filtered_matrices1, average_by_label)
rm(filtered_matrices1)

filtered_matrices2 <- readRDS("insuladev/dropbox/subcortex_matrices2.RDS")
averaged_matrices2 <- lapply(filtered_matrices2, average_by_label)
rm(filtered_matrices2)

filtered_matrices3 <- readRDS("insuladev/dropbox/subcortex_matrices3.RDS")
averaged_matrices3 <- lapply(filtered_matrices3, average_by_label)
rm(filtered_matrices3)

filtered_matrices4 <- readRDS("insuladev/dropbox/subcortex_matrices4.RDS")
averaged_matrices4 <- lapply(filtered_matrices4, average_by_label)
rm(filtered_matrices4)

filtered_matrices5 <- readRDS("insuladev/dropbox/subcortex_matrices5.RDS")
averaged_matrices5 <- lapply(filtered_matrices5, average_by_label)
rm(filtered_matrices5)

###reorganize matrices such that Tien regions are each in their own dataframe

subcortical_dfs <- vector("list", 50)

for (i in 1:50) {
  combined_data <- NULL
  
  # Loop through each of the 5 lists and extract the i-th row from each matrix
  for (j in 1:5) {
    current_list <- get(paste0("averaged_matrices", j))  # Get list of matrices from each matrix list

    for (k in 1:length(current_list)) {
      current_matrix <- current_list[[k]]  # run through each matrix
      row_data <- current_matrix[current_matrix$Label == i, -1, drop = FALSE]  # Extract specified row (not first column which just lists the Tien Label)
      combined_data <- rbind(combined_data, row_data)  # Combine rows
    }
  }
  
  rownames(combined_data) <- paste0("Subject_", 1:562)
  
  subcortical_dfs[[i]] <- combined_data
}

subject_data <- readRDS("insuladev/dropbox/HCPD_demographics.RDS")
subject_data$participant_id <- sub("^0+", "", subject_data$participant_id) #remove leading 0s to combined with modelarray

# Load FD_Means of scans
FD_Means_List <- read_csv("insuladev/HCPDbold_1/cor_results/model_array_demographics.csv")

subject_data_qc <- merge(subject_data, FD_Means_List[, c("participant_id", "mean_fd")], by = "participant_id")

Insula_Labels <- read_excel("insuladev/dropbox/Labels.xlsx")
insula_vertex_ids <- as.character(Insula_Labels$Label)

#combine across subcortical regions 
# attach family labels to each of the 50 region slots, morphing into the 16 region atlas
subcortex_labels <- data.frame(
  Row = 1:50,
  SubcortexFamily = c(
    rep("HIP", 4),
    "THAL_post", "THAL_post", "THAL_ant", "THAL_ant", "THAL_ant", "THAL_ant",
    rep("PUT", 4), rep("CAU", 4),
    "AMY", "AMY", "THAL_post", "NAc", "NAc", "GP", "GP",  
    rep("HIP", 4),
    "THAL_post", "THAL_post", "THAL_ant", "THAL_ant", "THAL_ant", "THAL_ant",
    rep("PUT", 4), rep("CAU", 4),
    "AMY", "AMY", "THAL_post", "NAc", "NAc", "GP", "GP")
)

# get the unique family names to iterate over
family_names <- unique(subcortex_labels$SubcortexFamily)

# for each family, average the connectivity matrices across its member regions
subcortical_dfs_family <- lapply(family_names, function(fam) {
  idx <- subcortex_labels$Row[subcortex_labels$SubcortexFamily == fam]
  Reduce("+", subcortical_dfs[idx]) / length(idx)
})
names(subcortical_dfs_family) <- family_names

#combined across bilateral insula subregions such that each matrix submitted to GAM is 
# 6 insula regions x 562 times points, one for each subcortical regions
insula_averages <- function(mat, Insula_Labels) {
  df <- as.data.frame(mat)
  
  # map each column (1..1478) to an insula ROI code
  g <- Insula_Labels$Label[match(seq_len(ncol(df)), Insula_Labels$ColumnNumber)]
  
  out <- sapply(split(seq_len(ncol(df)), g), function(idx) {
    rowMeans(df[, idx, drop = FALSE], na.rm = TRUE)
  })
  
  out <- as.data.frame(out)
  colnames(out) <- as.character(sort(unique(g)))
  
  # rename columns using ROI names
  insula_roi_names <- c(
    "106"="R_PoI2", "109"="R_MI", "111"="R_AVI", "112"="R_AAIC",
    "167"="R_PoI1", "168"="R_Ig",
    "286"="L_PoI2", "289"="L_MI", "291"="L_AVI", "292"="L_AAIC",
    "347"="L_PoI1", "348"="L_Ig"
  )
  colnames(out) <- insula_roi_names[colnames(out)]
  
  # average left and right homologous regions
  bilateral_regions <- c("PoI2", "MI", "AVI", "AAIC", "PoI1", "Ig")
  
  out_bilateral <- sapply(bilateral_regions, function(roi) {
    rowMeans(out[, c(paste0("L_", roi), paste0("R_", roi)), drop = FALSE], na.rm = TRUE)
  })
  
  as.data.frame(out_bilateral)
}
subcortex_averaged <- lapply(subcortical_dfs_family, insula_averages, Insula_Labels = Insula_Labels)
                  
process_subcort_gam <- function(subject_data_qc, subcortical_avgs, i) {
  data <- subcortical_avgs[[i]]                
  gam_results <- list() 
  for(col in 1:ncol(data)){
   subcortical_GAM_test <- data.frame(
     con_value = data[[col]],
     age = subject_data_qc$age,
     sex = as.factor(subject_data_qc$sex),
     scan_quality = subject_data_qc$mean_fd
   )
   if (all(is.na(subcortical_GAM_test$con_value))) {
        gam_results[[col]] <- NA # Skip if all values are NA
        } else {
          gam_model <- gam(con_value ~ s(age, k= 3) + sex + scan_quality, data = subcortical_GAM_test, method = "REML", na.action = na.exclude)
          gam_results[[col]] <- gam_model
                      }
            }
      # Return the list of GAM results
     return(gam_results)
}

process_subcort_null <- function(subject_data_qc, subcortical_avgs, i) {
  data <- subcortical_avgs[[i]] 
  null_results <- list()
  
  for(col in 1:ncol(data)){
    subcortical_GAM_test <- data.frame(
      con_value = data[[col]],
      age = subject_data_qc$age,
      sex = as.factor(subject_data_qc$sex),
      scan_quality = subject_data_qc$mean_fd
    )
    if (all(is.na(subcortical_GAM_test$con_value))) {
      null_results[[col]] <- NA # Skip if all values are NA
    } else {
      null_model <- gam(con_value ~ sex + scan_quality, data = subcortical_GAM_test, method = "REML", na.action = na.exclude)
      null_results[[col]] <- null_model
    }
  }

return(null_results)
} 

subcortical_results <- function(subject_data_qc, subcortical_avgs, i) {
    gam_models <- process_subcort_gam(subject_data_qc, subcortical_avgs, i)
    null_models <- process_subcort_null(subject_data_qc, subcortical_avgs, i)
  
    gam_smooth_results_list <- list()
    
    for(col in 1:length(gam_models)) {
      gam_model <- gam_models[[col]]
      null_model <- null_models[[col]]
      gam_result <- summary(gam_model)
      if (all(is.na(gam_model$model$con_value)) || all(is.na(null_model$model$con_value))) {
        gam_smooth_results_list[[col]] <- NA
      } else{
        new_data <- data.frame(age = seq(min(gam_model$model$age), max(gam_model$model$age), length.out = 100))
        new_data$sex <- rep("Male", length.out = nrow(new_data))
        new_data$scan_quality <- mean(gam_model$model$scan_quality) 

        pred <- predict.gam(gam_model, newdata = new_data, type = "terms", se.fit = TRUE)
        smooth_values <- pred$fit[, "s(age)"]
        se_values     <- pred$se.fit[, "s(age)"]
      
        derivative <- diff(smooth_values) / diff(new_data$age)
        se_values <- pred$se.fit    # Compute standard errors for predictions
      
        lower_ci <- smooth_values - 1.96 * se_values
        upper_ci <- smooth_values + 1.96 * se_values
     
        gam_smooth_F <- gam_result$s.table[3]
        gam_smooth_pvalue <- gam_result$s.table[4]
      
        anova_smooth_pvalue <- anova.gam(null_model,gam_model,test='F')$`Pr(>F)`[2]
      
        RSS_model <- sum((gam_model$y - gam_model$fitted.values)^2) #residual sum of squares for age model by comparing actual values (y) and predicted values (fitted.values)
        RSS_nullmodel <- sum((null_model$y - null_model$fitted.values)^2) #residual sum of squares of null model (no age)
        partialRsq <- abs((RSS_nullmodel - RSS_model) / RSS_nullmodel)
      
        mean_derivative <- mean(derivative)
        if(mean_derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
        partialRsq <- partialRsq*-1}
      
        insula_name <- colnames(subcortical_avgs[[i]])[col]
        
        gam_smooth_results <- data.frame(Row = i, Col = col, Insula = insula_name, GAM_smooth_Fvalue = as.numeric(gam_smooth_F), GAM_smooth_pvalue = as.numeric(gam_smooth_pvalue), 
                                       GAM_smooth_partialR2 = as.numeric(partialRsq), Anova_smooth_pvalue = as.numeric(anova_smooth_pvalue), Predicted_Age = as.numeric(new_data$age), Predicted_smooth = as.numeric(smooth_values), 
                                       Upper_CI_smooth = as.numeric(upper_ci), Lower_CI_smooth = as.numeric(lower_ci))
        gam_smooth_results_list[[col]] <- gam_smooth_results
    }
    } 
    return(gam_smooth_results_list) 
    }

subcortex_findings <- list()

for (i in family_names) {
  region_results <- subcortical_results(subject_data_qc, subcortex_averaged, i)
  
  subcortex_findings[[paste0("gam_results_region_", i)]] <- region_results
}

all_edges <- Filter(is.data.frame, do.call(c, subcortex_findings))
subcortex_df <- bind_rows(all_edges)

#apply FDR correction for each subcortical region
fdr_table <- subcortex_df %>%
  distinct(Insula, Row, Anova_smooth_pvalue, GAM_smooth_pvalue)
fdr_table$FDRp <- p.adjust(fdr_table$Anova_smooth_pvalue, method = "BH")


write_rds(subcortex_df, "dropbox/subcortical_GAM_results_new.RDS")
write_rds(fdr_table, "dropbox/subcortical_GAM_fdr.RDS")

#---------- if keeping full 50 subregions --------------
subcortex_labels <- data.frame(
  Row = 1:50,
  Subcortex = c(
    "HIP-head-m-rh","HIP-head-l-rh","HIP-body-rh","HIP-tail-rh",
    "THA-VPm-rh","THA-VPl-rh","THA-VAi-rh","THA-VAs-rh","THA-DAm-rh","THA-DAl-rh",
    "PUT-VA-rh","PUT-DA-rh","PUT-VP-rh","PUT-DP-rh",
    "CAU-VA-rh","CAU-DA-rh","CAU-body-rh","CAU-tail-rh",
    "lAMY-rh","mAMY-rh","THA-DP-rh","NAc-shell-rh","NAc-core-rh","pGP-rh","aGP-rh",
    "HIP-head-m-lh","HIP-head-l-l","HIP-body-l","HIP-tail-lh",
    "THA-VPm-lh","THA-VPl-lh","THA-VAi-lh","THA-VAs-lh","THA-DAm-lh","THA-DAl-lh",
    "PUT-VA-lh","PUT-DA-l","PUT-VP-lh","PUT-DP-lh",
    "CAU-VA-lh","CAU-DA-lh","CAU-body-lh","CAU-tail-lh",
    "lAMY-lh","mAMY-lh","THA-DP-lh","NAc-shell-lh","NAc-core-lh","pGP-lh","aGP-lh"
  )
)
