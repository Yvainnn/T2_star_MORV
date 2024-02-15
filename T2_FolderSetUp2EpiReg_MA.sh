#!/bin/bash

# - T2* MultiTE images pre-processing, 14/02/2024

# -------------------------------
# (1) Folder set up
#
#	|_ Subj directory
#		|_Session(n)
#		|_Session(n)
#			|_anat
#			|_fieldmap_2t2
#			|_fieldmap_diffphase
#			|_t2_multitemag
#			|_t2_singlete
#			|_multite_interp
#		|_Mask
#		|_....
#
# -------------------------------
#
# (2) DICOM to NII 
#	All files has been converted by using DICOM2NIIx (MRIcron)
#
# -------------------------------
#
# (3) Preparing folders usefull to preprocess all the images
#
# -------------------------------
#
# (4) MOTION CORRECT
#	performed by FSL mcflirt
#
# -------------------------------
#
# (5) DISTORT CORRECTION
#	(5.1) fsl_prepare_fieldmap use the fieldmaps images (magnitude and phase) to create a fieldmap image {This step could be very time consuming if runned by using not extracted mask and brain}
#	(5.2) Run of Epi_reg line
#
#  -----------------------------
#
# (6) APPLY WARP
#	The images Have been resampled in structural space and the distort correction has been applyed by using epi_reg t22str matrix 
#
# -----------------------------
#
# (7) Bias Correction
#
#-----------------------------




# (1)------- FOLDER set up ----------

for ((i=0;i<=45;i++));
	do 
		mkdir VP00"$i"
		cd /home/malberti/Unix_Folders/REMORV_processing/VP00"$i"

	for ((j=0;j<=6;j++));
		do
			mkdir /home/malberti/Unix_Folders/REMORV_processing/VP00"$i"/sess0"$j"
			cd /home/malberti/Unix_Folders/REMORV_processing/VP00"$i"/sess0"$j"
			sess_path=/home/malberti/Unix_Folders/REMORV_processing/VP00"$i"/sess0"$j"
			mkdir "$sess_path"/anat
			mkdir "$sess_path"/fieldmap_2t2
			mkdir "$sess_path"/fieldmap_diffphase
			mkdir "$sess_path"/t2_multitemag
			mkdir "$sess_path"/t2_singlete
			mkdir  "$sess_path"/multite_interp
			cd /home/malberti/Unix_Folders/REMORV_processing
		done
	echo "$i"
done


# (2)------ DICOM2NIIx -------- [The script is setted on REMORV prj - T2*] 

REMORV_dcm=/home/malberti/Unix_Folders/MORV/dcm  # Dicom main directory
main_path=/home/malberti/Unix_Folders/MORV_processing/ # Nifti and processing main directory
fmap_2t2=GRE_FIELD_MAPPING_0012 # Name of fieldmap images (fmap and phase) as dicom 
fmap_phase=GRE_FIELD_MAPPING_0013
mag_name=T2_SWI3D_MULTITE_LOWRES_0010 # Name of T2-MultiTE dicom


for ((z=1;z<=3;z++)); #session
	do

		for ((i=6;i<=9;i++)); #subj ID <VP0010. These loop convert DICOM to NIFTI, moves/rename files in $main_path bids directory
			do 
				./dcm2niix "$REMORV_dcm"/REMORV_VP00"$i"_REMORV_VP00"$i"/SS"$z"/"$fmap_2t2"
				mv "$REMORV_dcm"/REMORV_VP00"$i"_REMORV_VP00"$i"/SS"$z"/"$fmap_2t2"/*.nii /"$main_path"/VP00"$i"/sess0"$s"/fieldmap_2t2
		
				./dcm2niix "$REMORV_dcm"/REMORV_VP00"$i"_REMORV_VP00"$i"/SS"$z"/"$mag_name"
				mv "$REMORV_dcm"/REMORV_VP00"$i"_REMORV_VP00"$i"/SS"$z"/"$mag_name"/*.nii /"$main_path"/VP00"$i"/sess0"$s"/t2_multitemag
		
				./dcm2niix "$REMORV_dcm"/REMORV_VP00"$i"_REMORV_VP00"$i"/SS"$z"/"$fmap_phase"
				mv "$REMORV_dcm"/REMORV_VP00"$i"_REMORV_VP00"$i"/SS"$z"/"$fmap_phase"/*.nii /"$main_path"/VP00"$i"/sess0"$s"/fieldmap_diffphase
			done
	
		for ((i=10;i<=45;i++)); #subj ID >VP0010. These loop convert DICOM to NIFTI, moves/rename files in $main_path bids directory
			do 

				./dcm2niix "$REMORV_dcm"/REMORV_VP0"$i"_REMORV_VP0"$i"/SS"$z"/"$fmap_2t2"
				mv "$REMORV_dcm"/REMORV_VP0"$i"_REMORV_VP0"$i"/SS"$z"/"$fmap_2t2"/*.nii /"$main_path"/VP00"$i"/sess0"$s"/fieldmap_2t2
		
				./dcm2niix "$REMORV_dcm"/REMORV_VP0"$i"_REMORV_VP0"$i"/SS"$z"/"$mag_name"
				mv "$REMORV_dcm"/REMORV_VP0"$i"_REMORV_VP0"$i"/SS"$z"/"$mag_name"/*.nii /"$main_path"/VP00"$i"/sess0"$s"/t2_multitemag
		
				./dcm2niix "$REMORV_dcm"/REMORV_VP0"$i"_REMORV_VP0"$i"/SS"$z"/"$fmap_phase"
				mv "$REMORV_dcm"/REMORV_VP0"$i"_REMORV_VP0"$i"/SS"$z"/"$fmap_phase"/*.nii /"$main_path"/VP00"$i"/sess0"$s"/fieldmap_diffphase
			done
	done



# (3)-----Set up of folders and some files

cd /home/malberti/Unix_Folders/MORV_processing/VP00"$i"

for ((i=3;i<=34;i++)); #subj ID
	do
		echospacing_multite=0.0000164 #set up of fieldmap and epireg parameters
		echospacing_singlete=0.00017
		pedir_multite=x
		pedir_singlete=-y
		fsl_mni=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz
		fsl_mni_mask=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz
		sbj_path=/home/malberti/Unix_Folders/MORV_processing/VP00"$i" 
	
		for ((j=4;j<=6;j++)); #sessions ID
			do	
				sess_path=/home/malberti/Unix_Folders/MORV_processing/VP00"$i"/sess0"$j" # Path and folder set up
		
				mkdir "$sess_path"/motioncorr_multite 
				mkdir "$sess_path"/distortcor
				mkdir "$sess_path"/fmap_setup
				mkdir "$sess_path"/multite_interp	
				moco_out="$sess_path"/motioncorr_multite/moco_qt2.nii.gz #motion correct output
				er_out="$sess_path"/distortcor #epi_reg output folder 
				fmap_out="$sess_path"/fmap_setup #fieldmap_prepare input/output folder
				fmap_2t2_sb="$sess_path"/fieldmap_2t2_sb #fieldmap for svenja's B. niftiD folder
				mat_epi2std="$sess_path"/distortcor/epi2std_init.mat
				t1_2mni_warp="$sbj_path"/sess02/anat/mnispace/xfms
		
				fslmerge -t "$sess_path"/fieldmap_2t2/fmap_mag.nii.gz "$sess_path"/fieldmap_2t2/GRE_FIELD_MAPPING*.nii #merge of fieldmap files
				mv "$sess_path"/fieldmap_diffphase/GRE_FIELD_MAPPING_0*ph.nii "$sess_path"/fieldmap_diffphase/fmap_phase.nii.gz

				for ((mag=1;mag<=9;mag++));		#rename of multite magnitude files and merge them in a timeseries
					do 
						mv T2_SWI3D_MULTITE_LOWRES_0*e"$mag".nii "$sess_path"/t2_multitemag/qt2_te"$mag".nii
					done
		
				fslmerge -t "$sess_path"/t2_multitemag/qt2_multite_mag.nii.gz "$sess_path"/t2_multitemag/qt2_te*.nii
		
				echo "$sess_path" "    " "subj"$i""	"setted - NEXT motion correct..."	# Just to taking notes about how the script is working
	
# (4)------- MOTION CORRECT-----------


				mv "$sess_path"/t2_multitemag/qt2_multite_mag.nii.gz "$sess_path"/pre_moco_mte.nii.gz  #rename time series, obtained by merging magnitude(9 diff te)
				mcflirt -in "$sess_path"/pre_moco_mte.nii.gz -refvol 4 -mats -plots -o "$sess_path"/motioncorr_multite/moco_qt2.nii.gz  -report #Motion correct by using FSL, as refvol u should use an images in the middle of the 'time series', this image has 9 volume, you have to select the 5th volume. fsl start counting from 0, it means that the 5th volume is the image number 4
		
				#	----- Check the motion correct results
					#	By using grep i select the translations/rotation angles values, then i create a file that contain them,the script do it for each MAT0000. 
					#	before running this script pls be sure that the "$moco_path" contains the matrix list as txt file and check if it's right {I wrote this script BUT i'm not sure if it still works}
		
				for matrix in `cat "$moco_path"/matrix_list.txt`; # Quality check of motion correction
					do 
						avscale --allparams "$moco_path"/output.mat/"$matrix" | grep "Translations" > "$moco_path"/avscale_trans_"$filename".txt  # These files (avscale*) are temporaty, each loop overwrite them, if don't wanna that, just switch ">" into ">>" 
						avscale --allparams "$moco_path"/output.mat/"$matrix" | grep "Rotation Angles" > "$moco_path"/avscale_trans_"$filename".txt
						awk '{print $5, $6, $7}' "$moco_path"/avscale_trans_"$filename".txt >> "$moco_path"/Translations_"$filename".txt #select only the column that contain the translations and the rotation angles 
						awk '{print $6, $7, $8}' "$moco_path"/avscale_rot_angles_"$filename".txt >> "$moco_path"/Rotation_Angles_"$filename".txt
					done
				I=1 
					#these 2 loop check if these files contain values > 1, if values are >1 it prints the values and the subj ID
				for values in `cat "$moco_path"/Translations_"$filename".txt`;
					do
						if [ "$values" > 1 ] ; then
							echo "Check Motion-correct, Translations:"$values" sbj "$i", session "$j" "
						fi
					done

				for values in `cat "$moco_path"/Rotation_Angles_"$filename".txt`;
					do
						if [ "$values" > 1 ]; then
							echo "Check Motion-correct, Rotation Angles:"$values"sbj "$i", session "$j" "
						fi
					done
				echo "subj VP00"$i" sess0"$j" NEXT EPI_REG"
				fslroi "$moco_out" "$er_out"/sbref.nii.gz 4 1 # multi TE reference image
				bet "$moco_out" "$er_out"/sbref_brain -R -f 0.2 -g 0 -m # create a mask to interpolatation on matlab
	
				fslmaths "$sess_path"/fieldmap_2t2/fmap_mag.nii.gz  -Tmean "$fmap_out"/fmap_mean_mag.nii.gz 
				echo "$sess_path" "    " "subj"$i""
				bet "$fmap_out"/fmap_mean_mag.nii.gz "$fmap_out"/fmap_mean_mag_brain  -R -f 0.6 -g 0 -m #create extracted brain image, usefull for fsl_prepare_fieldmap
				fslmaths "$sess_path"/fieldmap_diffphase/fmap_phase.nii.gz -mas  "$fmap_out"/fmap_mean_mag_brain_mask.nii.gz  "$fmap_out"/fmap_phase_masked.nii.gz  #these should be the fieldmap_diffphase image, afeter bet	
				fsl_prepare_fieldmap SIEMENS  "$fmap_out"/fmap_phase_masked.nii.gz "$fmap_out"/fmap_mean_mag_brain.nii.gz  "$fmap_out"/fieldmap.nii.gz 2.46   # prepare fieldmap from 

			epi_reg --epi="$er_out"/sbref.nii.gz --t1=/home/malberti/Unix_Folders/MORV_processing/VP00"$i"/sess"$j"/anat/mprage/t1w_acpc.nii.gz \
					--t1brain=/home/malberti/Unix_Folders/MORV_processing/VP00"$i"/sess0"$j"/anat/mprage/t1w_acpc_brain.nii.gz \
					--out="$er_out"/epi2std \
					--fmap="$fmap_out"/fieldmap.nii.gz  \
					--fmapmag="$fmap_out"/fmap_mean_mag.nii.gz \
					--fmapmagbrain="$fmap_out"/fmap_mean_mag_brain_mask.nii.gz \
					--echospacing=0.0000164 --pedir=x #> epi_reg_log_VP00"$i"_sess0"$j".txt

#------ APPLY distort correction 		

			applywarp --ref=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz \
					--in="$m_qt2"/VP00"$i"_sess0"$j"_m_qt2_tc.nii \
					--out="$m_qt2"/VP00"$i"_sess0"$j"_m_qt2_tc_2mni.nii.gz \
					--warp="$t1_2mni_warp"/acpc_dc2standard.nii.gz  \  		#	Resample into standard space
					--premat="$mat_epi2std" # Apply distort correction 
		
			echo "subj VP00"$i" sess0"$j" done"


#------ Other correction:  ( i didn't apply it/them to multite T2*)
			mkdir "$m_qt2"/biascorr 
			N4BiasFieldCorrection -d 3 -i "$m_qt2"/VP00"$i"_sess0"$j"_m_qt2_tc_2mni.nii.gz -o "$m_qt2"/biascorr/VP00"$i"_sess0"$j"_m_qt2_tc_2mni_unbiased.nii.gz -v
	done
done 

# ----- RUN MATLAB fitting_t2.m ------- Multite interp	
		

 
	
##----------------RANDOM INFO and MOTION COR check script--------------------------

#COMMAND: 
#     N4BiasFieldCorrection
#          N4 is a variant of the popular N3 (nonparameteric nonuniform normalization) 
#          retrospective bias correction algorithm. Based on the assumption that the 
#          corruption of the low frequency bias field can be modeled as a convolution of 
#          the intensity histogram by a Gaussian, the basic algorithmic protocol is to 
#          iterate between deconvolving the intensity histogram by a Gaussian, remapping 
#          the intensities, and then spatially smoothing this result by a B-spline modeling 
#          of the bias field itself. The modifications from and improvements obtained over 
#          the original N3 algorithm are described in the following paper: N. Tustison et 
#          al., N4ITK: Improved N3 Bias Correction, IEEE Transactions on Medical Imaging, 
#          29(6):1310-1320, June 2010. 
#
#OPTIONS: 
#     -d, --image-dimensionality 2/3/4
#          This option forces the image to be treated as a specified-dimensional image. If 
#          not specified, N4 tries to infer the dimensionality from the input image. 

#     -i, --input-image inputImageFilename
#          A scalar image is expected as input for bias correction. Since N4 log transforms 
#          the intensities, negative values or values close to zero should be processed 
#          prior to correction. 

#     -x, --mask-image maskImageFilename
#          If a mask image is specified, the final bias correction is only performed in the 
#          mask region. If a weight image is not specified, only intensity values inside 
#          the masked region are used during the execution of the algorithm. If a weight 
#          image is specified, only the non-zero weights are used in the execution of the 
#          algorithm although the mask region defines where bias correction is performed in 
#          the final output. Otherwise bias correction occurs over the entire image domain. 
#          See also the option description for the weight image. If a mask image is *not* 
#          specified then the entire image region will be used as the mask region. Note 
#          that this is different than the N3 implementation which uses the results of Otsu 
#          thresholding to define a mask. However, this leads to unknown anatomical regions 
#          being included and excluded during the bias correction. 

#     -r, --rescale-intensities 0/(1)
#          At each iteration, a new intensity mapping is calculated and applied but there 
#          is nothing which constrains the new intensity range to be within certain values. 
#          The result is that the range can "drift" from the original at each iteration. 
#          This option rescales to the [min,max] range of the original image intensities 
#          within the user-specified mask. A mask is required to perform rescaling. 

#     -w, --weight-image weightImageFilename
#          The weight image allows the user to perform a relative weighting of specific 
#          voxels during the B-spline fitting. For example, some studies have shown that N3 
#          performed on white matter segmentations improves performance. If one has a 
#          spatial probability map of the white matter, one can use this map to weight the 
#          b-spline fitting towards those voxels which are more probabilistically 
#          classified as white matter. See also the option description for the mask image. 

#     -s, --shrink-factor 1/2/3/(4)/...
#          Running N4 on large images can be time consuming. To lessen computation time, 
#          the input image can be resampled. The shrink factor, specified as a single 
#          integer, describes this resampling. Shrink factors <= 4 are commonly used.Note 
#          that the shrink factor is only applied to the first two or three dimensions 
#          which we assume are spatial. 

#     -c, --convergence [<numberOfIterations=50x50x50x50>,<convergenceThreshold=0.0>]
#          Convergence is determined by calculating the coefficient of variation between 
#          subsequent iterations. When this value is less than the specified threshold from 
#          the previous iteration or the maximum number of iterations is exceeded the 
#          program terminates. Multiple resolutions can be specified by using 'x' between 
#          the number of iterations at each resolution, e.g. 100x50x50. 

#     -b, --bspline-fitting [splineDistance,<splineOrder=3>]
#                           [initialMeshResolution,<splineOrder=3>]
#          These options describe the b-spline fitting parameters. The initial b-spline 
#          mesh at the coarsest resolution is specified either as the number of elements in 
#          each dimension, e.g. 2x2x3 for 3-D images, or it can be specified as a single 
#          scalar parameter which describes the isotropic sizing of the mesh elements. The 
#          latter option is typically preferred. For each subsequent level, the spline 
#          distance decreases in half, or equivalently, the number of mesh elements doubles 
#          Cubic splines (order = 3) are typically used. The default setting is to employ a 
#          single mesh element over the entire domain, i.e., -b [1x1x1,3]. 

#     -t, --histogram-sharpening [<FWHM=0.15>,<wienerNoise=0.01>,<numberOfHistogramBins=200>]
#          These options describe the histogram sharpening parameters, i.e. the 
#          deconvolution step parameters described in the original N3 algorithm. The 
#          default values have been shown to work fairly well. 

#     -o, --output correctedImage
#                  [correctedImage,<biasField>]
#          The output consists of the bias corrected version of the input image. 
#          Optionally, one can also output the estimated bias field. 

#     --version 
#          Get Version Information. 

#     -v, --verbose (0)/1
#          Verbose output. 

#     -h 
#          Print the help menu (short version). 

#     --help 
#          Print the help menu

#		ls "$moco_path"/output.mat/ > "$moco_path"/matrix_list.txt
#
#		rm "$moco_path"/avscale_trans_"$filename".txt
#		rm "$moco_path"/avscale_trans_"$filename".txt

		
#	DISTORTION CORRECTION with regular gradient echo field map
#	HCP uses correction with fm output jacobians and bias field
#      from T1*T2 imagesN$ SiemensFieldMapPreprocessing muss gelaufen se

#	echo-spacing/dwell time and phase encoding direction (from protocol) for EPI in seconds
# 	dwell= 0.00017 #for singlete, for multi no idea; effective echospacing = echospacing/acceleration factor=ipat/GRAPPA (multiband has no influence %hier 0.00072/2
# 	pe_dir = {'-y', 'x'}; % %normal sollte A>>P = -y; P>>A = y; R>>L = x; L>>R = -x;
# 	tr = [4.1 0.05];  %2.3;vp003r2 2.0; TR in s %or get from image with fslval $img pixdim4 (to see what keyword need to use: fslhd)
#	fmri_finalres = 2; %final resolution in mm -- if not 2 or 1mm need to resample T1w image accordingly so have a ref file. see below

#	How to create a fieldmap magnitude image: 
#	Usage: fsl_prepare_fieldmap <scanner> <phase_image> <magnitude_image> <out_image> <deltaTE (in ms)> [--nocheck]
# 	Prepares a fieldmap suitable for FEAT from SIEMENS or GEHC data - saves output in rad/s format
# 	<scanner> must be SIEMENS or GEHC_FIELDMAPHZ
# 	<phase_image> should be the phase difference for SIEMENS and the fieldmap in HERTZ for GEHC_FIELDMAPHZ
# 	<magnitude image> should be Brain Extracted (with BET or otherwise)
#  	<deltaTE> is the echo time difference of the fieldmap sequence - find this out form the operator
#             (defaults are *usually* 2.46ms on SIEMENS)
#            (defaults are *usually* 2.304ms for GEHC 2D-B0MAP at 3.0T and 2.272 ms GEHC 3D B0MAP at 3.0T)
#  	--nocheck supresses automatic sanity checking of image size/range/dimensions
#   	e.g. fsl_prepare_fieldmap SIEMENS images_3_gre_field_mapping images_4_gre_field_mapping fmap_rads 2.65
#   	e.g. fsl_prepare_fieldmap GEHC_FIELDMAPHZ 3dB0map_fieldmaphz mag_3dB0map fmap_rads 2.272

#	as fieldmap you should use mag and phase fieldmap ottenute dallo scanner, then you can create your combined fieldmap by using fsl_prepare_fieldmap 
#	t1_path=/home/malberti/Unix_Folders/REMORV_processing/VP003/sess01/anat/mprage\

		
#output created by using old Svenja's fiedlmap
		#epi_reg --epi="$er_out"/sbref.nii.gz \
		#--t1=/home/malberti/Unix_Folders/MORV_processing/VP00"$i"/sess0"$j"/anat/mprage/t1w_acpc.nii.gz \
		#--t1brain=/home/malberti/Unix_Folders/MORV_processing/VP00"$i"/sess0"$j"/anat/mprage/t1w_acpc_brain.nii.gz \
		#--out="$er_out"/epi2std_sb \
		#--fmap="$fmap_2t2_sb"/fieldmap.nii.gz \
		#--fmapmag="$fmap_2t2_sb"/meanmagnitude.nii.gz \
		#--fmapmagbrain="$fmap_2t2_sb"/meanmagnitude_brain.nii.gz \
		#--echospacing=0.0000164 --pedir=x



