#HCP_YA
#create insula connectivity matrix averaged across all participants

library(dplyr)
library(ciftiTools)
library(data.table)

conmat_files <- list.files(path = "insuladev/HCP_YA/matrices", full.names = TRUE, "^.*\\.tsv$", recursive = TRUE)

# Read the first CSV to get dimensions
first_conmat <- fread(conmat_files[1], header = TRUE)  
nrow_conmat <- nrow(first_conmat[,-1])
ncol_conmat <- ncol(first_conmat) - 1

row_names <- first_conmat[[1]]

#create empty average matrix
sum_conmat <- matrix(0, nrow = nrow_conmat, ncol = ncol_conmat)
valid_entries_count <- matrix(0, nrow = nrow_conmat, ncol = ncol_conmat)

for (i in 1:length(conmat_files)) {
  conmat <- fread(conmat_files[i], header = TRUE) # Read current file as a data.table
  conmat_matrix <- as.matrix(conmat[,-1, with=FALSE]) #remove first column (rownames)
  valid_entries_count[!is.na(conmat_matrix)] <- valid_entries_count[!is.na(conmat_matrix)] + 1 # Count valid (non-NA) entries and update valid_entries_count
  
  sum_conmat <- sum_conmat + conmat_matrix # Add the current matrix to the sum matrix
  print(paste("Finished:", conmat_files[i]))
}

# Compute the element-wise average
averaged_conmat <- sum_conmat / valid_entries_count # Divide the sum by the count of valid (non-NA) entries
averaged_conmat_df <- data.table(row_names, averaged_conmat)

insula_names <- c("Left_Ig", "Left_PoI1", "Left_PoI2", "Left_MI", "Left_AVI", "Left_AAIC","Right_Ig", "Right_PoI1", "Right_PoI2", "Right_MI", "Right_AVI", "Right_AAIC")

filtered_conmat <- averaged_conmat_df[which(row_names %in% insula_names),]
filtered_column_indices <- which(colnames(filtered_conmat) %in% insula_names)
final_conmat <- filtered_conmat[,-c(filtered_column_indices), with = FALSE]

# Write the result to a CSV file
write.csv(as.data.frame(final_conmat), "insuladev/HCP_YA/averaged_vertex_connectivity.csv", row.names = TRUE)
