# ExploreDTI_scripts
This repository holds several explanations of ExploreDTI scripts.
* Note that several scripts need editing depending on your data!!
* With help from the ExploreDTI forum; https://groups.google.com/forum/#!forum/e_dti
* Also credits to the Cardiff; http://sites.cardiff.ac.uk/cubric/cubric-users/user-documentation/mri-resources/mri-how-to/using-exploredti-via-the-command-line/

## For batch analysing tracts:
```
E_DTI_batch_tracts_analysis(A,B,C,D,E,F) 
```
with these arguments:
* "A" is path of folder of *.nii masks (AND/NOT/OR ROIs = value: 1/2/3 with 0 for background) - only one ROI per *.nii mask.
* "B" is path of tracts *.mat file name.
* "C" is path of DTI *.mat file name.
* "D" is path of output (analyzed) tracts *.mat file name.
* "E" is 0 (no) or 1 (yes) for "segment only" option (works on the "outer" AND ROIs).
* "F" is length range (in units mm) option (e.g., [0 inf] or [50 500]).

## Freesurfer labels-derived connectivity matrix approach
```
E_DTI_Network_analysis_exe(dti_files, tract_files, vol_files, A_T, A_L, Lab, VDims_T, out_folder, vol_ext, selected_labels, ACh, Temp_fol)
```
with these arguments
* dti_files = name of DTI *.mat file
* tract_files = name of tract *.mat file
* vol_files = [] (set empty)
* A_T = 3D matrix of atlas template. To get this variable, do:
* [A_T, VDims_T] = E_DTI_load_nii(file_name_of_template_nii_file); % see templates folder in ExploreDTI
* A_T = double(A_T);
* A_T = A_T/max(A_T(:));
* A_L = 3D matrix of atlas labels (zero = background, regions with values: 1, 2, 3, …, N)
* Lab = cell array with label names, row number corresponds to label value in A_L (e.g. Lab{1,1} = 'Some_structure'; Lab{2,1} = 'Some_other_structure'; etc)
* VDims_T  = voxel size of atlas template (see above)
* out_folder = output folder
* vol_ext = ‘’ (set empty)
* selected_labels = [] (keep empty)
* ACh = 1 (PASS) or 2 (END)
* Temp_fol = temporary folder

## Convert folder of DICOMs to initial ExploreDTI .mat file
```
dti_filename = E_DTI_Script_Get_DTI_folders(source_dir,target_dir)
```
* source_dir = Source directory containing folders of DICOMs
* target_dir = Directory to save DTI .mat file
* dti_filename = Filename of the output .mat file.

## Flip/permute nifti files
```
E_DTI_flip_permute_nii_file_exe(t1_filename_new,param);
```
* t1_filename_new= 'full_path_to_nii_file'.
* param=[];
* param.suff= '_FP';
* param.permute= [1 2 3];
* param.flip= [0 0 0];

## Mask the background of the T1_FP files to zero to reduce the computation time of the EPI correction
```
E_DTI_mask_3D_nii_file_exe(t1_filename_new,paramaters);
```
* t1_filename_new= 'full_path_to_nii_file'.
* paramaters.mfs=3; % Kernel size: 3     
* paramaters.Threshold=0.02; % Threshold: 0.02
You can play around with the threshold based on the SNR of your data, but I found that this suggested value works ok- it’s a 2% threshold

## SM/EC/EPI correction via .txt file
```
E_DTI_SMECEPI_Main(parameter_filename);
```
* parameter_filename = name of the text file containing the parameters for the SM/EC/EPI correction. The file is as saved in the GUI using Settings > SM/EC/EPI correction > Export parameter file. 
* par.RESTORE_option = 1;
* par.TE.NS = 4;
* par.TE.TS = 4;
* par.R2D.type = 3;
* par.R2D.FN = '_T1_FP_masked.nii';
* par.R2D.contrast = 2;
* par.EPI.Num_Resol = 4;
* par.EPI.Deriv_Scales = [1  0  0];

## Whole-brain tractography (dRL)
```
WholeBrainTracking_dRL_SH_FOD(dti_filename, tracts_filename, paramaters);
```
* dti_filename = Name of the DTI .mat file
* tracts_filename = Name of the output tracts .mat file
* paramaters = a structure containing the parameter for the tractography algorithm;

Below are the parameters with suggested values for dRL
* paramaters.SeedPointRes = [2 2 2]
* paramaters.StepSize = 1
* paramaters.FODThresh = 0.0500
* paramaters.AngleThresh = 45
* paramaters.FiberLengthRange = [20 500]
* paramaters.lambda = 0.0019
* paramaters.beta = 1.777405329610160e-04
* paramaters.iter = 200
* paramaters.rr = 8
* paramaters.nn = 0.0400


## Whole-brain tractography (DTI)
```
WholeBrainTrackingDTI(filename_in, filename_out, parameters);
```
* paramaters.SeedPointRes = [3 3 3]
* paramaters.StepSize = 1
* paramaters.FAThresh = 0.2000
* paramaters.AngleThresh = 45
* paramaters.FiberLengthRange = [50 500]

## Whole-brain tractography (CSD)
```
WholeBrainTrackingCSD_fast(filename_in, filename_out, parameters);
```
* paramaters.SeedPointRes = [2 2 2]
* paramaters.StepSize = 1
* paramaters.AngleThresh = 30
* paramaters.FiberLengthRange = [50 500]

## Sampling a volume along tracts
```
E_DTI_Export_Tract_Volume_Info(nii_filename, tracts_filename,dti_filename, output_directory);
```
* nii_filename = NIFTI image of volume to sample (must be in same space as DTI data)
* tracts_filename = Tracts .mat filename
* dti_filename = DTI .mat filename
* output_directory = Directory to save results to

## Create Connectivity Matrices
```
E_DTI_Network_analysis_exe(dti_filename,tracts_filename,img_filename,A_T,A_L,L,VDims_A,mat_dir,img_suffix,selected_labels,ACh);
```
* dti_filename = name of the DTI .mat file
* tracts_filename = name of the DTI tracts .mat file
* img_filename = name of another image modality (e..g MWM; in the same space as the DTI data) to estimate weights for. Use [] to disable.
* A_T = 3D volume of the atlas (e.g. AAL) reference template. Note this volume must be read into MATLAB using the command: A_T=E_DTI_load_nii(‘template_filename’)
* A_L = 3D volume of the atlas (e.g. AAL) labels. As above, the volume need to be loaded into MATLAB * A_L=E_DTI_load_nii(‘labels_filename’)
* L = cell array of label names for each region in A_L.
* VDims_A = 1×3 vector of atlas label dimension.
* mat_dir = Directory to save the connectivity matrices
* img_suffix = Suffix for other imaging filename (e..g "_MWM.nii"). Use [] to disable.
* selected_labels = a n_labelsx2 matrix of pairs of label indices to create an individual tract .mat file for that pair of labels. Use [] to disable. Use nchoosek(1:n_atlas_labels,2) to output for every pair of labels.
* ACh = Which criteria to use determine connectivity: 1 = PASS, 2=END, 3=BOTH

