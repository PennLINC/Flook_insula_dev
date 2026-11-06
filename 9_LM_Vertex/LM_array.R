#running LM on vertex data via array

library(ModelArray)
library(rhdf5)

#get participant ID from array shell command 
args <- commandArgs(trailingOnly = TRUE)
voxel_id <- args[1] #will loop through all vertices called in exec script array argument

# filenames
h5_path <- paste0("insuladev/HCPDbold_1/H5_files/voxel",voxel_id, ".h5")
csv_path <- "insuladev/HCPDbold_1/cor_results/model_array_demographics.csv"
vertex_labels_path <- "insuladev/dropbox/Vertex_Labels.xlsx"

# create a ModelArray object
modelarray <- ModelArray(h5_path, scalar_types = "Pearson_R", analysis_names = "results_gam")

# load CSV demographic file:
phenotypes <- read.csv(csv_path)
phenotypes$sex <- factor(phenotypes$sex, levels = c("Male", "Female"))
phenotypes$age_demeaned <- phenotypes$age - mean(phenotypes$age) #calculate demeaned age

# LM formulas:
demeaned.lm <- Pearson_R ~ age_demeaned + I(age_demeaned^2) + sex + mean_fd
formula.lm <- Pearson_R ~ poly(age,2) + sex + mean_fd
null.lm <- Pearson_R ~ sex + mean_fd
linear.lm <- Pearson_R ~ age + sex + mean_fd
quad.lm <- Pearson_R ~ I(age^2) + sex + mean_fd

# run LM fitting with ModelArray.lm(), running only cortex (first 59412 vertices)
print(paste("running models for voxel", voxel_id))
myfull <- ModelArray.lm(formula.lm, modelarray, phenotypes, "Pearson_R", verbose = F, element.subset = 1:59412)

mynull <- ModelArray.lm(null.lm, modelarray, phenotypes, "Pearson_R", verbose = F, element.subset = 1:59412)

myquad <- ModelArray.lm(quad.lm, modelarray, phenotypes, "Pearson_R", verbose = F, element.subset = 1:59412)

mylinear <- ModelArray.lm(linear.lm, modelarray, phenotypes, "Pearson_R", verbose = F, element.subset = 1:59412)

mydemean <- ModelArray.lm(demeaned.lm, modelarray, phenotypes, "Pearson_R", verbose = F, element.subset = 1:59412)

#save outputs both to H5 file and as separate RDS files
writeResults(h5_path, df.output = myfull, analysis_name = "results_full_lm")
writeResults(h5_path, df.output = mynull, analysis_name = "results_null")
writeResults(h5_path, df.output = mylinear, analysis_name = "results_linear")
writeResults(h5_path, df.output = myquad, analysis_name = "results_quad")

outpath <- paste0("insuladev/HCPDbold_1/Polynomial_Results/LM_full_results_voxel",voxel_id,".RDS")
saveRDS(myfull, file = outpath)
outpath_null <- paste0("insuladev/HCPDbold_1/Polynomial_Results/LM_null_results_voxel",voxel_id,".RDS")
saveRDS(mynull, file = outpath_null)
outpath_linear <- paste0("insuladev/HCPDbold_1/Polynomial_Results/LM_linear_results_voxel",voxel_id,".RDS")
saveRDS(mylinear, file = outpath_linear)
outpath_quad <- paste0("insuladev/HCPDbold_1/Polynomial_Results/LM_quad_results_voxel",voxel_id,".RDS")
saveRDS(myquad, file = outpath_quad)
outpath <- paste0("insuladev/HCPDbold_1/Polynomial_Results/LM_demeaned_results_voxel",voxel_id,".RDS")
saveRDS(mydemean, file = outpath)