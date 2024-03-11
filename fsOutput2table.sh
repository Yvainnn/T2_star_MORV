#!/bin/bash

%sb can you make the file names more transparent please.

# This script extracts the volume values from freesurfer|recon-all output and prepare table to load in matlab. 
    # I created 4 differnt tab, containing:
        #  brain_reg_  rh/lh as different hemi = {'rh_bankssts_volume', 'rh_caudalanteriorcingulate_volume', 'rh_caudalmiddlefrontal_volume', \
        #                                            'rh_cuneus_volume', 'rh_entorhinal_volume', 'rh_fusiform_volume', 'rh_inferiorparietal_volume', \
        #                                            'rh_inferiortemporal_volume', 'rh_isthmuscingulate_volume', 'rh_lateraloccipital_volume', \
        #                                            'rh_lateralorbitofrontal_volume', 'rh_lingual_volume', 'rh_medialorbitofrontal_volume', \
        #                                            'rh_middletemporal_volume', 'rh_parahippocampal_volume', 'rh_paracentral_volume', \
        #                                            'rh_parsopercularis_volume', 'rh_parsorbitalis_volume', 'rh_parstriangularis_volume', \
        #                                            'rh_pericalcarine_volume', 'rh_postcentral_volume', 'rh_posteriorcingulate_volume', \
        #                                            'rh_precentral_volume', 'rh_precuneus_volume', 'rh_rostralanteriorcingulate_volume', \
        #                                            'rh_rostralmiddlefrontal_volume', 'rh_superiorfrontal_volume', 'rh_superiorparietal_volume', \
        #                                            'rh_superiortemporal_volume', 'rh_supramarginal_volume', 'rh_frontalpole_volume', \
        #                                            'rh_temporalpole_volume', 'rh_transversetemporal_volume', 'rh_insula_volume'}; 
       
        # brain_aseg (rh+lh) =  {'Left-Lateral-Ventricle', 'Left-Inf-Lat-Vent','Left-Cerebellum-White-Matter', 'Left-Cerebellum-Cortex', \
        #                           'Left-Thalamus-Proper', 'Left-Caudate', 'Left-Putamen', 'Left-Pallidum', '3rd-Ventricle', '4th-Ventricle', \
        #                           'Brain-Stem', 'Left-Hippocampus', 'Left-Amygdala', 'CSF', 'Left-Accumbens-area', 'Left-VentralDC', 'Left-vessel', 'Left-choroid-plexus', \
        #                           'Right-Lateral-Ventricle', 'Right-Inf-Lat-Vent', 'Right-Cerebellum-White-Matter', 'Right-Cerebellum-Cortex', 'Right-Thalamus-Proper', \
        #                           'Right-Caudate', 'Right-Putamen', 'Right-Pallidum', 'Right-Hippocampus', 'Right-Amygdala', 'Right-Accumbens-area', 'Right-VentralDC', \
        #                           'Right-vessel', 'Right-choroid-plexus', '5th-Ventricle', 'WM-hypointensities', 'Left-WM-hypointensities	Right-WM-hypointensities', \
        #                           'non-WM-hypointensities', 'Left-non-WM-hypointensities', 'Right-non-WM-hypointensities', 'Optic-Chiasm', 'CC_Posterior', \
        #                           'CC_Mid_Posterior', 'CC_Central', 'CC_Mid_Anterior', 'CC_Anterior'};

        # gm_wm_morv_aseg = { "Brain Segmentation Volume Without Ventricles,  "Total cortical white matter volume," \
        #                     "Total cortical gray matter volume",  "Subcortical gray matter volume", "Estimated Total Intracranial Volume", "Total gray matter volume"} 


# Create a file that contain a colon with all subj name, WITHOOUT INFO ABOUT SESSION: 
    #     ls > subj.txt 
    #     gedit subj.txt # the order of subj has been setted manually, following the order (VPS): 
    #
    #                    vps = [6 7 14:21 23 28 31 37 43 44 4 5 8:11 13 22 24 26 32:35 38 40]; 
    #                    sleep = [6 7 14:21 23 28:31 36 37 43:45];
    #                    wake = [3:5 8:13 22 24 26 27 32:35 38:41];
    #
    # Moreover, aslo the session "ID" (_SS0*) has been removed manually  %sb why?
    # Here is a script that should do it (haven't tried yet) %sb then do it?
    # 
    # vps=(6 7 14 15 16 17 18 19 20 21 23 28 31 37 43 44 45 4 5 8 9 10 11 13 22 24 26 32 33 34 35 38 40)  # Bash array don't use ':' to say "values from x to y" %sb but there must be an equivalent. This way will easily produce errors.
    # for ((f=0;f<=n;f++)); #set the number fo subject %sb what value is n???
    #   do 
    #       if [ "$f" < 10 ] ; then %sb this is not the way to deal with the leading zeros. find the solution. i can only tell you in matlab.
    #            folder_filename=VP00  #(folder ID <10)
	#		else		
    #            folder_filename=VP0  #(folder ID >10)
    #        fi
    #       
    #        subj_name="$folder_filename"${vps[f]}
    #        echo "$subj_name" >> subj.txt #this file will be used only for gm_wm_morv_aseg file. 
    #
    #   done
    #
    # cat subj.txt #check IT pls %sb: you check it. did you check that for a random subject that the correct values go in there? %ma check it was a reminder to me ahaahah, yes i did
    
#!/bin/bash 



vps=(6 7 {14..21} 23 {28..31} 36 37 {43..45} {3..5} {8..13} 22 24 26 27 {32..35} {38..41}) #the numberof sbj should be different between t2* and t1 output. es. no t2* for subj 45
sleep=(6 7 {14..21} 23 {28:31} 36 37 {43..45})
wake=({3..5} {8..13} 22 24 26 27 {32..35} {38..41})

main_path=/home/malberti/Desktop #this is the output path for txt files
fs_path=/home/malberti/Desktop/dicom
SUBJECTS_DIR=/home/malberti/Unix_Folders/niftiE/freesurfer_mprage/subjects

filename=VP0 
sessions=2

for ((i=0; i<${#vps[@]}; i++)); 
do 
	echo "${vps[i]}"


done


for ((j=1;j<="$sessions";j++)); 
do 
    echo "BS_volumewithoutvent" > "$main_path"/BS_volumewithoutvent.txt  # these lines create/clean temporari files  %sb for what do we use this value? %ma these files will contain the fs output (volume/thickness) i added this line in order to be sure to create a new one if i need to run it several times
    echo "TT_cortical_wm" > "$main_path"/TT_cortical_wm.txt
    echo "TT_cortical_gm" > "$main_path"/TT_cortical_gm.txt
    echo "TT_subcortical_gm" >  "$main_path"/TT_subcortical_gm.txt
    echo "TT_ETIV" > "$main_path"/TT_ETIV.txt #%sb I don't get your logic of naming. this is ETIV right? why is it called _gm as if it was gray matter????
    echo "TT_gm" > TT_gm.txt
    echo "subj order" > "$main_path"/subj_order.txt

   
	for ((i=0; i<${#vps[@]}; i++)); 
	do 
		if [ ${vps[i]} -gt 9 ]; then
			stats_dir=$SUBJECTS_DIR/"$filename""${vps[i]}"_SS0"$j"/stats/aseg.stats 
		else 
			stats_dir=$SUBJECTS_DIR/"$filename""0""${vps[i]}"_SS0"$j"/stats/aseg.stats 
		fi 
		
		#echo "$stats_dir" " right?" "Press 'y/n' to continue/close:"
		#read -r choice


		#if [ "$choice" = "y" ]; then
		#	echo "Running folder setup."
		#elif [ "$choice" = "n" ]; then
		#	echo "Exiting the script..."
		#	exit 0
		#else
		#	echo "Invalid choice. Please enter 'y/n' or 'all'."
		#fi
	
		grep "Brain Segmentation Volume Without Ventricles," "$stats_dir"  | awk '{print $10}'  >> "$main_path"/BS_volumewithoutvent.txt # grep select the line, by looking for the "" string, and then, awk print the column 9 or 10 that contain the volume 
		grep "Total cortical white matter volume," "$stats_dir" | awk '{print $10}'  >> "$main_path"/TT_cortical_wm.txt
		grep "Total cortical gray matter volume," "$stats_dir"  | awk '{print $10}'  >> "$main_path"/TT_cortical_gm.txt
		grep "Subcortical gray matter volume," "$stats_dir" | awk '{print $9}'  >> "$main_path"/TT_subcortical_gm.txt
		grep "Estimated Total Intracranial Volume," "$stats_dir" | awk '{print $9}'  >> "$main_path"/TT_ETIV.txt
		grep "Total gray matter volume," "$stats_dir" | awk '{print $9}'  >> "$main_path"/TT_gm.txt
		echo "subj:" "${vps[i]}" >> "$main_path"/subj_order.txt
		 
		if [ ${vps[i]} -gt 9 ]; then
			FS_sub="$filename""${vps[i]}"
			fs_array+="$filename""${vps[i]}"" "
		else 
			FS_sub="$filename""0""${vps[i]}"
			fs_array+="$filename""0""${vps[i]}"" "
		fi  
		
		echo "$FS_sub"_SS0"$j" >> "$main_path"/aparc_list_sess"$j".txt
        done
   	
  
   	paste "$main_path"/BS_volumewithoutvent.txt \
   		"$main_path"/TT_cortical_wm.txt \
   		"$main_path"/TT_cortical_gm.txt \
   		"$main_path"/TT_subcortical_gm.txt \
   		"$main_path"/TT_ETIV.txt TT_gm.txt > "$main_path"/gm_wm_morv_aseg_ss"$j".csv # this line merge all the txt files in table.csv file that could be imported on matlab or excel, the following line create a 2th file that doesn't contain the first row and the subj ID
    	sed -i '1d' gm_wm_morv_aseg_ss"$j".csv 
 	#read -r
	aparcstats2table --subjectsfile="$main_path"/aparc_list_sess"$j".txt -m volume --hemi lh --common-parc -v --skip -t "$main_path"/lh_morv_volume_ss0"$j".csv #print lh.parc.stats file 
	aparcstats2table --subjectsfile="$main_path"/aparc_list_sess"$j".txt -m volume --hemi rh --common-parc -v --skip -t "$main_path"/rh_morv_volume_ss0"$j".csv

    	cut -f 2- "$main_path"/lh_morv_volume_ss0"$j".csv > "$main_path"/NOID_lh_morv_volume_ss0"$j".csv #remove the first colon
    	cut -f 2- "$main_path"/rh_morv_volume_ss0"$j".csv > "$main_path"/NOID_rh_morv_volume_ss0"$j".csv #remove the first colon
    
   	asegstats2table  --subjectsfile="$main_path"/aparc_list_sess"$j".txt -m volume --skip --common-segs -t "$main_path"/aseg_morv_volume_ss0"$j".csv -v
    	cut -f 2- "$main_path"/aseg_morv_volume_ss0"$j".csv > "$main_path"/NOID_aseg_morv_volume_ss0"$j".csv
	
	rm "$main_path"/aparc_list_sess"$j".txt
done


%sb where is CSF? can you put this also into the same textfiles as the other values we are currently using?
%ma yes ofc, csf is in aseg files

