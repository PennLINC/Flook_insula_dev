#Insula ROI to Subcortex connectivity
#create connectivity matrices for the last rows of the matrices

library(dplyr)
library(readxl)
library(data.table)
library(rhdf5)

conmat_files <- list.files(path = "insuladev/HCPDbold_1/cor_results", full.names = TRUE, "^correlation_\\d+\\.h5$", recursive = TRUE)

Insula_Labels <- read_excel("insuladev/dropbox/Labels.xlsx")
Insula_roi_names <- c("106" = 'R_PoI2',"109" = 'R_MI',"111" =  'R_AVI',"112" = 'R_AAIC', 
                      "167" = "R_PoI1", "168" = 'R_Ig', "286" = 'L_PoI2', "289" = 'L_MI',	"291" = 'L_AVI',
                      "292" = 'L_AAIC', "347" = 'L_PoI1', "348" = 'L_Ig')
Insula_Labels <- Insula_Labels %>%
  mutate(Label = recode(as.character(Label), !!!Insula_roi_names))

Vertex_Labels <- read_excel("insuladev/dropbox/Vertex_Assignments.xlsx", sheet = "Tien_vertex_code")

Tien_labels = merge(Vertex_Labels, Insula_Labels[,1:2], by.x = "Vertex", by.y = "Vertex", all.x = TRUE)

rm(Vertex_Labels, Insula_Labels)
gc()

non_zero_labels <- Tien_labels %>%
  filter(Label.x > 0)

row_numbers <- non_zero_labels$Vertex # Get the row numbers from non-zero labels
labels <- non_zero_labels$Label.x 

subset_connectivity_matrix <- function(conmat_file, non_zero_labels) {
  conmat <- h5read(conmat_file, "insula_pearsonr")

  filtered_conmat <- conmat[row_numbers, ] # Subset the connectivity matrix to only keep rows with non-zero labels
  filtered_conmat <- cbind(Label = labels, filtered_conmat)
  
  return(filtered_conmat)
  print0("Finished file")
}

# Use lapply to apply the function to all connectivity matrices
filtered_matrices <- lapply(conmat_files[401:length(conmat_files)], subset_connectivity_matrix, non_zero_labels = non_zero_labels)

saveRDS(filtered_matrices, file = "insuladev/dropbox/subcortex_matrices5.RDS")