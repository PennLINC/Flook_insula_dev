#Change row of 1s to NANs before running GAM or will error

library(rhdf5)
library(stringr)
library(fs)

# Define the list of voxel file paths
file_paths <- list.files("insuladev/HCPDbold_1/H5_files", pattern = "voxel\\d+\\.h5", full.names = TRUE)

# Function to update each voxel file
update_voxel_file <- function(file_path) {
  # Open the H5 file with read-write access
  h5file <- H5Fopen(file_path, "H5F_ACC_RDWR")
  
  # Read the values matrix from the specified group
  values <- h5read(h5file, "/scalars/Pearson_R/values")
  
  # Identify rows that are all 1s
  rows_to_modify <- which(apply(values, 1, function(row) all(row == 1)))
  print(rows_to_modify)
  
  # Replace those rows with NaN
  values[rows_to_modify, ] <- NaN
  
  # Check if any rows are still all 1s (for debugging)
  remaining_ones <- which(apply(values, 1, function(row) all(row == 1)))
  print(paste("Remaining all-1 rows in file:", file_path))
  print(remaining_ones)
  
  # Write the modified values back to the H5 file
  h5write(values, h5file, "/scalars/Pearson_R/values")
  
  # Close the H5 file
  H5Fclose(h5file)
  
  # Clean up
  h5closeAll()
}

# Apply the update function to each file in the list
lapply(file_paths, update_voxel_file)