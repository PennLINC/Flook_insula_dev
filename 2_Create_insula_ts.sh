##create insula-restricted time series from wholebrain using workbench
#files in insuladev/HCPDbold_1/rs_files/

List_Filenames="insuladev/dropbox/HCPD_rs_file_list.csv"
Insula_Folder="insuladev/HCPDbold_1/Insula_Masks/"

while IFS=, read -r filepath; do 
filepath=$(echo "$filepath" | tr -d '"')
subject_folder=$(dirname "$filepath")
subject_folder=$(basename "$subject_folder")
    
    if [ -d "$subject_folder" ]; then
    filename=$(basename "$filepath" .dtseries.nii)
    echo "filename is $filename"    
        for cifti_roi_file in "${Insula_Folder}Combined_Right.dscalar.nii" "${Insula_Folder}Combined_Left.dscalar.nii"; do
            cifti_roi_filename=$(basename "$cifti_roi_file" .dscalar.nii)
            insula_filename="${subject_folder}/${filename}_insula_${cifti_roi_filename}.dtseries.nii"
            if [ ! -f "$insula_filename" ]; then
                echo "Processing $filepath - $cifti_roi_filename"
                wb_command -cifti-restrict-dense-map "$filepath" COLUMN "$insula_filename" -cifti-roi "$cifti_roi_file" 
            else   
                echo "The file $insula_filename already exists"
            fi
        done
    else
        echo "Subject folder $subject_folder does not exist"
    fi
done < "$List_Filenames"

##for subjects with atypical name due to acquisition variant
wb_command -cifti-restrict-dense-map sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-AP_run-01_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii COLUMN sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-AP_run-01_space-fsLR_den-91k_desc-denoised_bold_insula_Combined_Right.dtseries.nii -cifti-roi insuladev/HCPDbold_1/Insula_Masks/Combined_Right.dscalar.nii

wb_command -cifti-restrict-dense-map sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-02_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii COLUMN sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-02_space-fsLR_den-91k_desc-denoised_bold_insula_Combined_Right.dtseries.nii -cifti-roi insuladev/HCPDbold_1/Insula_Masks/Combined_Right.dscalar.nii

wb_command -cifti-restrict-dense-map sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-03_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii COLUMN sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-03_space-fsLR_den-91k_desc-denoised_bold_insula_Combined_Right.dtseries.nii -cifti-roi insuladev/HCPDbold_1/Insula_Masks/Combined_Right.dscalar.nii

wb_command -cifti-restrict-dense-map sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-05_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii COLUMN sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-05_space-fsLR_den-91k_desc-denoised_bold_insula_Combined_Right.dtseries.nii -cifti-roi insuladev/HCPDbold_1/Insula_Masks/Combined_Right.dscalar.nii

wb_command -cifti-restrict-dense-map sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-AP_run-01_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii COLUMN sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-AP_run-01_space-fsLR_den-91k_desc-denoised_bold_insula_Combined_Left.dtseries.nii -cifti-roi insuladev/HCPDbold_1/Insula_Masks/Combined_Left.dscalar.nii

wb_command -cifti-restrict-dense-map sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-02_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii COLUMN sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-02_space-fsLR_den-91k_desc-denoised_bold_insula_Combined_Left.dtseries.nii -cifti-roi insuladev/HCPDbold_1/Insula_Masks/Combined_Left.dscalar.nii

wb_command -cifti-restrict-dense-map sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-03_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii COLUMN sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-03_space-fsLR_den-91k_desc-denoised_bold_insula_Combined_Left.dtseries.nii -cifti-roi insuladev/HCPDbold_1/Insula_Masks/Combined_Left.dscalar.nii

wb_command -cifti-restrict-dense-map sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-05_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii COLUMN sub-1529756/sub-1529756_ses-V1_task-rest_acq-VARIANTNumVolumes_dir-PA_run-05_space-fsLR_den-91k_desc-denoised_bold_insula_Combined_Left.dtseries.nii -cifti-roi insuladev/HCPDbold_1/Insula_Masks/Combined_Left.dscalar.nii
