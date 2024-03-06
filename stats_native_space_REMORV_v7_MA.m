%stats in native space

%% setup
clear all
vps = [6 7 14:21 23 28 31 37 43:45 4 5 8:11 13 22 24 26 32:35 38 40];
sleep = [6 7 14:21 23 28:31 36 37 43:45];
wake = [3:5 8:13 22 24 26 27 32:35 38:41];
nsub=length(vps);
cond = {'L0','L1','L2','C0','C1','C2'; 'L1-L0','L2-L0','L2-L1','C1-C0','C2-C0','C2-C1'};
whichval = {'raw','zallm','zgm'};
%mask_name = {'allm', 'total_gm', 'wm_aseg','precuneus','l_inf_parietal','l_fusiform','lingual','l_lat_occ'};
mask_name = {'allm', 'total_gm', 'wm_aseg','precuneus','l_inf_parietal','l_fusiform','lingual','l_lat_occ','hc','IPL','SPL','PCC','RSC','cortical_gm','subcortical_gm','CSF'};
 fs_path= 'H:\Unix_Folders\niftiE\freesurfer_mprage\subjects';
 nii_path = 'H:\Unix_Folders\niftiE\niftiE';
 d_path = 'H:\Unix_Folders\MORV\niftiD\niftiD\'; %mask path is different from t2_singlete nii files
%H:\Unix_Folders\MORV\niftiD\diffusion_analysis\nativespace_2mm T2star
% NiftiD H:\Unix_Folders\MORV\niftiD\niftiD

% fs output brain region, use this insted of mask_name
brain_reg_rh = {'rh_bankssts_volume', 'rh_caudalanteriorcingulate_volume', 'rh_caudalmiddlefrontal_volume', 'rh_cuneus_volume', 'rh_entorhinal_volume', 'rh_fusiform_volume', 'rh_inferiorparietal_volume', 'rh_inferiortemporal_volume', 'rh_isthmuscingulate_volume', 'rh_lateraloccipital_volume', 'rh_lateralorbitofrontal_volume', 'rh_lingual_volume', 'rh_medialorbitofrontal_volume', 'rh_middletemporal_volume', 'rh_parahippocampal_volume', 'rh_paracentral_volume', 'rh_parsopercularis_volume', 'rh_parsorbitalis_volume', 'rh_parstriangularis_volume', 'rh_pericalcarine_volume', 'rh_postcentral_volume', 'rh_posteriorcingulate_volume', 'rh_precentral_volume', 'rh_precuneus_volume', 'rh_rostralanteriorcingulate_volume', 'rh_rostralmiddlefrontal_volume', 'rh_superiorfrontal_volume', 'rh_superiorparietal_volume', 'rh_superiortemporal_volume', 'rh_supramarginal_volume', 'rh_frontalpole_volume', 'rh_temporalpole_volume', 'rh_transversetemporal_volume', 'rh_insula_volume'}; 
brain_reg_lh = {'lh_bankssts_volume', 'lh_caudalanteriorcingulate_volume', 'lh_caudalmiddlefrontal_volume', 'lh_cuneus_volume', 'lh_entorhinal_volume', 'lh_fusiform_volume', 'lh_inferiorparietal_volume', 'lh_inferiortemporal_volume', 'lh_isthmuscingulate_volume', 'lh_lateraloccipital_volume', 'lh_lateralorbitofrontal_volume', 'lh_lingual_volume', 'lh_medialorbitofrontal_volume', 'lh_middletemporal_volume', 'lh_parahippocampal_volume', 'lh_paracentral_volume', 'lh_parsopercularis_volume', 'lh_parsorbitalis_volume', 'lh_parstriangularis_volume', 'lh_pericalcarine_volume', 'lh_postcentral_volume', 'lh_posteriorcingulate_volume', 'lh_precentral_volume', 'lh_precuneus_volume', 'lh_rostralanteriorcingulate_volume', 'lh_rostralmiddlefrontal_volume', 'lh_superiorfrontal_volume', 'lh_superiorparietal_volume', 'lh_superiortemporal_volume', 'lh_supramarginal_volume', 'lh_frontalpole_volume', 'lh_temporalpole_volume', 'lh_transversetemporal_volume', 'lh_insula_volume'}; 
brain_aseg= {'Left-Lateral-Ventricle', 'Left-Inf-Lat-Vent','Left-Cerebellum-White-Matter', 'Left-Cerebellum-Cortex', 'Left-Thalamus-Proper', 'Left-Caudate', 'Left-Putamen', 'Left-Pallidum', '3rd-Ventricle', '4th-Ventricle', 'Brain-Stem', 'Left-Hippocampus', 'Left-Amygdala', 'CSF', 'Left-Accumbens-area', 'Left-VentralDC', 'Left-vessel', 'Left-choroid-plexus', 'Right-Lateral-Ventricle', 'Right-Inf-Lat-Vent', 'Right-Cerebellum-White-Matter', 'Right-Cerebellum-Cortex', 'Right-Thalamus-Proper', 'Right-Caudate', 'Right-Putamen', 'Right-Pallidum', 'Right-Hippocampus', 'Right-Amygdala', 'Right-Accumbens-area', 'Right-VentralDC', 'Right-vessel', 'Right-choroid-plexus', '5th-Ventricle', 'WM-hypointensities', 'Left-WM-hypointensities	Right-WM-hypointensities', 'non-WM-hypointensities', 'Left-non-WM-hypointensities', 'Right-non-WM-hypointensities', 'Optic-Chiasm', 'CC_Posterior', 'CC_Mid_Posterior', 'CC_Central', 'CC_Mid_Anterior', 'CC_Anterior'};
gm_wm_m= {'BS_volumewithoutvent', 'TT_cortical_wm', 'TT_cortical_gm', 'TT_subcortical_gm', 'TT_intracranial_gm', 'TT_gm'};
%This loops should create the rawvals matrix
for i=1:length(gm_wm_m)   % [6 7 14:21 23 28 31 37 43 44 4 5 8:11 13 22 24 26 32:35 38 40 45]
        for t=1:6
            aseg_morv=readtable(fullfile(fs_path,sprintf('gm_wm_morv_aseg_ss0%d.csv', t)));
            orig_data=readtable(fullfile(fs_path,sprintf('subj_ss0%d.txt', t)));
            %for j=length(vps)
            subj_col=aseg_morv(:,i);
            orig_dat_sbj(:,t)=orig_data(:,1);
            subj_array=table2array(subj_col);
            region(:,t)=subj_array;       
        %end
        end
        brain_gm_wm(i)={region};
end 

% gm_volume=cell2mat(brain_gm_wm(6));
% intracranial_volume=cell2mat(brain_gm_wm(5));
% tt_div=gm_volume./intracranial_volume;
% brain_gm_wm(7)={tt_div};

save(fullfile(nii_path,'sleepana_rawvals_morv_VBM_gm_wm.mat'),'vps','sleep','wake',...
  'gm_wm_m', 'brain_gm_wm','whichval','cond', 'orig_dat_sbj');


% load(fullfile(d_path,'rawvals_remorv.mat'))
% load(fullfile(d_path,'rawvals_cleanvals_remorv.mat'))
% load(fullfile(d_path,'sleepana_rawvals_cleanvals_remorv.mat'))

%% get rawvalues in masks
%inputnames = spm_select('FPList',d_path,'^*_0mm\.nii.gz');
mask_size = NaN(nsub,length(mask_name));
for i = 1:length(vps)
    vp_path = fullfile(d_path,sprintf('VP%03d',vps(i)));
    if ismember(vps(i),wake) %wake
        sw = 'w';
    else %sleep
        sw = 's';
    end
    for m = 1:length(mask_name) %masks         
        %load mask
        mask_path = fullfile(vp_path,'masks','native','fs','2mm',sprintf('%s.nii.gz', mask_name{m}));
        mas = niftiread(mask_path);
        mask_size(i,m) = sum(mas(:));        
        for t = 1:6 %timepoints
            %load diffimg
            %if t <= 3 %SS1-3 /L1-3
                t2_path = fullfile(nii_path,sprintf('VP%03d',vps(i)),sprintf('SS%d',t),'t2_singlete','rnative_m_qt2_tc.nii.gz');
                %%%%% LOAD SESS1-3
            %else
                %%%% load 4-6
               % diff_path = fullfile(d_path,sprintf('MD_C%d_%s_VP%03d_0mm.nii.gz',t-4,sw,vps(i)));
           % end
            img = niftiread(t2_path);
            rawvals{i,m}(:,t) = reshape(img(mas==1),[mask_size(i,m),1]); 
        end        
    end                             
end

% save(fullfile(d_path,'rawvals_remorv.mat'),'vps','sleep','wake',...
%     'mask_name', 'mask_size', 'rawvals');
save(fullfile(nii_path,'sleepana_rawvals_morv_MA.mat'),'vps','sleep','wake',...
   'mask_name', 'mask_size', 'rawvals','whichval','cond');

% %save

%% clean vals

%perc_incl  = NaN(length(vps),length(brain_reg_rh),2);%percent of mask voxels included
cvals = cell(length(rawvals),size(rawvals,2),2,3); %cleaned values;
%dim3: single conditions L0;L1.../ differences L1-0, L2-0...;
%dim4: rawvalues / zallm / zgm
params_ind2 = NaN(length(vps),6,length(mask_name),4,3); %to test single cond with ind2 exclusions
%dim4:mean/median/sd/sem; dim5:raw/zallm/zgm
% 
% for i = 3:length(vps)
%     for m = 3:length(mask_name)
%         v = rawvals{i,m};
%         vdiff = [v(:,2)-v(:,1) v(:,3)-v(:,1) v(:,3)-v(:,2)...
%             v(:,5)-v(:,4) v(:,6)-v(:,4) v(:,6)-v(:,5)];
%         %exclude implausible values
%        if m < length(mask_name)
%            imp1 = [30 90]; %raw val
%            imp2 = 50;%difference
%        else
%            imp1 = [30 500]; %raw val
%            imp2 = 200;%difference - so no
%        end
%         ind1 = find(v(:,1)>imp1(1) & v(:,2)>imp1(1) & v(:,3)>imp1(1) & v(:,4)>imp1(1) & v(:,5)>imp1(1) & v(:,6)>imp1(1) & ...
%             v(:,1)<imp1(2) & v(:,2)<imp1(2) & v(:,3)<imp1(2) & v(:,4)<imp1(2) & v(:,5)<imp1(2) & v(:,6)<imp1(2));
%         ind2 = find(v(:,1)>imp1(1) & v(:,2)>imp1(1) & v(:,3)>imp1(1) & v(:,4)>imp1(1) & v(:,5)>imp1(1) & v(:,6)>imp1(1) & ...
%             v(:,1)<imp1(2) & v(:,2)<imp1(2) & v(:,3)<imp1(2) & v(:,4)<imp1(2) & v(:,5)<imp1(2) & v(:,6)<imp1(2) & ...
%             abs(v(:,2)-v(:,1))<imp2 & abs(v(:,3)-v(:,1))<imp2 & abs(v(:,3)-v(:,2))<imp2 & abs(v(:,5)-v(:,4))<imp2 & abs(v(:,6)-v(:,4))<imp2 & abs(v(:,6)-v(:,5))<imp2);
%      %  perc_incl(i,m) = [length(ind1)/length(v) length(ind2)/length(v)]; 
% 
% 
%         cvals{i,m,1,1}= v(ind1,:);
%         cvals{i,m,2,1}= vdiff(ind2,:);
%         params_ind2(i,:,m,1,1) = mean(v(ind2,:));
%         params_ind2(i,:,m,2,1) = median(v(ind2,:));
%         params_ind2(i,:,m,3,1) = std(v(ind2,:));
%         params_ind2(i,:,m,4,1) = std(v(ind2,:))/sqrt(length(ind2));

for i = 1:length(vps)
    for m = 1:length(mask_name)
        v = rawvals{i,m};
        vdiff = [v(:,2)-v(:,1) v(:,3)-v(:,1) v(:,3)-v(:,2)...
            v(:,5)-v(:,4) v(:,6)-v(:,4) v(:,6)-v(:,5)];
        %exclude implausible values
        if m < length(mask_name)
            imp1 = [30 90]; %raw val
            imp2 = 50;%difference
        else
            imp1 = [300 500]; %raw val
            imp2 = 200;%difference - so no
        end
        ind1 = find(v(:,1)>imp1(1) & v(:,2)>imp1(1) & v(:,3)>imp1(1) & v(:,4)>imp1(1) & v(:,5)>imp1(1) & v(:,6)>imp1(1) & ...
            v(:,1)<imp1(2) & v(:,2)<imp1(2) & v(:,3)<imp1(2) & v(:,4)<imp1(2) & v(:,5)<imp1(2) & v(:,6)<imp1(2));
        ind2 = find(v(:,1)>imp1(1) & v(:,2)>imp1(1) & v(:,3)>imp1(1) & v(:,4)>imp1(1) & v(:,5)>imp1(1) & v(:,6)>imp1(1) & ...
            v(:,1)<imp1(2) & v(:,2)<imp1(2) & v(:,3)<imp1(2) & v(:,4)<imp1(2) & v(:,5)<imp1(2) & v(:,6)<imp1(2) & ...
            abs(v(:,2)-v(:,1))<imp2 & abs(v(:,3)-v(:,1))<imp2 & abs(v(:,3)-v(:,2))<imp2 & abs(v(:,5)-v(:,4))<imp2 & abs(v(:,6)-v(:,4))<imp2 & abs(v(:,6)-v(:,5))<imp2);
        perc_incl(i,m,:) = [length(ind1)/length(v) length(ind2)/length(v)]; 

        cvals{i,m,1,1}= v(ind1,:);
        cvals{i,m,2,1}= vdiff(ind2,:);
        params_ind2(i,:,m,1,1) = mean(v(ind2,:));
        params_ind2(i,:,m,2,1) = median(v(ind2,:));
        params_ind2(i,:,m,3,1) = std(v(ind2,:));
        params_ind2(i,:,m,4,1) = std(v(ind2,:))/sqrt(length(ind2));

    end
end
        % if m < length(mask_name)
        %     %zstand: zallm
        %     if m == 1
        %         mallm = mean(v(ind2,:));
        %         stdallm = sqrt(var(v(ind2,:)));
        %         v1 = v; % to calculate zgm also for allm
        %     end 
        %     vallm = NaN(length(v),6);
        %     for c1 = 1:6 %for all columns
        %     vallm(:,c1) = (v(:,c1)-mallm(c1))/stdallm(c1);
        %     end
        %     vdiffallm = [vallm(:,2)-vallm(:,1) vallm(:,3)-vallm(:,1) vallm(:,3)-vallm(:,2)...
        %         vallm(:,5)-vallm(:,4) vallm(:,6)-vallm(:,4) vallm(:,6)-vallm(:,5)];
        %     cvals{i,m,1,2}= vallm(ind2,:);
        %     cvals{i,m,2,2}= vdiffallm(ind2,:);
        %     params_ind2(i,:,m,1,2) = mean(vallm(ind2,:));
        %     params_ind2(i,:,m,2,2) = median(vallm(ind2,:));
        %     params_ind2(i,:,m,3,2) = std(vallm(ind2,:));
        %     params_ind2(i,:,m,4,2) = std(vallm(ind2,:))/sqrt(length(ind2));
        %     %zstand: zgm        
        %     if m == 2
        %         mgm = mean(v(ind2,:));
        %         stdgm = sqrt(var(v(ind2,:)));
        %         %calculate for allm (m==1)
        %         vgm = NaN(length(v1),6);
        %         for c1 = 1:6 %for all columns
        %             vgm(:,c1) = (v1(:,c1)-mgm(c1))/stdgm(c1);
        %         end
        %         vdiffgm = [vgm(:,2)-vgm(:,1) vgm(:,3)-vgm(:,1) vgm(:,3)-vgm(:,2)...
        %             vgm(:,5)-vgm(:,4) vgm(:,6)-vgm(:,4) vgm(:,6)-vgm(:,5)];
        %         cvals{i,m-1,1,3}= vgm(ind2,:);
        %         cvals{i,m-1,2,3}= vdiffgm(ind2,:);
        %         params_ind2(i,:,m-1,1,3) = mean(vgm(ind2,:));
        %         params_ind2(i,:,m-1,2,3) = median(vgm(ind2,:));
        %         params_ind2(i,:,m-1,3,3) = std(vgm(ind2,:));
        %         params_ind2(i,:,m-1,4,3) = std(vgm(ind2,:))/sqrt(length(ind2));
        %     end
        %     vgm = NaN(length(v),6);
        %     if m >= 2
        %         for c1 = 1:6 %for all columns
        %         vgm(:,c1) = (v(:,c1)-mgm(c1))/stdgm(c1);
        %         end
        %         vdiffgm = [vgm(:,2)-vgm(:,1) vgm(:,3)-vgm(:,1) vgm(:,3)-vgm(:,2)...
        %             vgm(:,5)-vgm(:,4) vgm(:,6)-vgm(:,4) vgm(:,6)-vgm(:,5)];
        %         cvals{i,m,1,3}= vgm(ind1,:);
        %         cvals{i,m,2,3}= vdiffgm(ind2,:);
        %         params_ind2(i,:,m,1,3) = mean(vgm(ind2,:));
        %         params_ind2(i,:,m,2,3) = median(vgm(ind2,:));
        %          params_ind2(i,:,m,3,3) = std(vgm(ind2,:));
        %         params_ind2(i,:,m,4,3) = std(vgm(ind2,:))/sqrt(length(ind2));
        %     end
        % end                       

% 
% % %save
% % % save(fullfile(d_path,'rawvals_cleanvals_remorv.mat'),'vps','sleep','wake',...
% % %     'mask_name', 'mask_size', 'perc_incl','rawvals','cvals','params_ind2');
 save(fullfile(nii_path,'sleepana_rawvals_cleanvals_remorv_MA.mat'),'vps','sleep','wake',...
     'mask_name', 'mask_size', 'perc_incl','rawvals','cvals','params_ind2','-v7.3');

%% additionally exclude 'border voxels' to be sure to remove partial volume effects
%inputnames = spm_select('FPList',d_path,'^*_0mm\.nii.gz');
mask_size = NaN(nsub,length(mask_name));
cbvals = cell(length(rawvals),size(rawvals,2),2,3); %cleaned and no border values;
%dim3: single conditions L0;L1.../ differences L1-0, L2-0...;
%dim4: rawvalues / zallm / zgm
params_ind3 = NaN(length(vps),6,length(mask_name),4,3); %to test single cond with ind3 exclusions
params_ind4 = NaN(length(vps),6,length(mask_name),4,3); %to test single cond with ind2+3 exclusions
%dim4:mean/median/sd/sem; dim5:raw/zallm/zgm
 
for i = 1:length(vps)
    vp_path = fullfile(d_path,sprintf('VP%03d',vps(i)));
    if ismember(vps(i),wake) %wake
        sw = 'w';
    else %sleep
        sw = 's';
    end
    %to exclude border voxels, get fsbrainmask, cut out CSF mask, erode once
    mask_path = fullfile(vp_path,'masks','native','fs','2mm','fs_brainmask.nii.gz');
    brainmask = niftiread(mask_path);
    mask_path = fullfile(vp_path,'masks','native','fs','2mm','CSF.nii.gz');
    csfmask = niftiread(mask_path);
    brainmask(csfmask==1) = 0;
%     SE = strel('square',3); %always >90% included
    SE = strel('cube',3);
    brainmaskErode = imerode(brainmask,SE);
    indErode = find(brainmaskErode==1);
    for m = 1:length(mask_name)
        v = rawvals{i,m};
        vdiff = [v(:,2)-v(:,1) v(:,3)-v(:,1) v(:,3)-v(:,2)...
            v(:,5)-v(:,4) v(:,6)-v(:,4) v(:,6)-v(:,5)];  
        %exclude implausible values
        if m < length(mask_name)
            imp1 = [30 100]; %raw val
            imp2 = 50;%difference
        else
            imp1 = [30 500]; %raw val
            imp2 = 200;%difference - so no
        end
%         ind1 = find(v(:,1)>imp1(1) & v(:,2)>imp1(1) & v(:,3)>imp1(1) & v(:,4)>imp1(1) & v(:,5)>imp1(1) & v(:,6)>imp1(1) & ...
%             v(:,1)<imp1(2) & v(:,2)<imp1(2) & v(:,3)<imp1(2) & v(:,4)<imp1(2) & v(:,5)<imp1(2) & v(:,6)<imp1(2));
        ind2 = find(v(:,1)>imp1(1) & v(:,2)>imp1(1) & v(:,3)>imp1(1) & v(:,4)>imp1(1) & v(:,5)>imp1(1) & v(:,6)>imp1(1) & ...
            v(:,1)<imp1(2) & v(:,2)<imp1(2) & v(:,3)<imp1(2) & v(:,4)<imp1(2) & v(:,5)<imp1(2) & v(:,6)<imp1(2) & ...
            abs(v(:,2)-v(:,1))<imp2 & abs(v(:,3)-v(:,1))<imp2 & abs(v(:,3)-v(:,2))<imp2 & abs(v(:,5)-v(:,4))<imp2 & abs(v(:,6)-v(:,4))<imp2 & abs(v(:,6)-v(:,5))<imp2);
        %exclude border voxels
        %load mask
        mask_path = fullfile(vp_path,'masks','native','fs','2mm',sprintf('%s.nii.gz', mask_name{m}));
        mas = niftiread(mask_path);
        indRaw = find(mas==1);
        
%         %if erode ROI mask itself - but that leaves for most masks below 10%
% %         SE = strel('cube',3);
%         SE = strel('square',3);
%         imgErode = imerode(mas,SE);
%         ind3 = find(imgErode==1);
%         indRaw(~ismember(indRaw,ind3))= NaN;
        %if just take away brain border voxels
        if m < length(mask_name)
            indRaw(~ismember(indRaw,indErode))= NaN;
        else %CSF - erode CSF mask
          %SE = strel('cube',3);
            SE = strel('square',3);
            imgErode = imerode(mas,SE);
            ind3 = find(imgErode==1);
            indRaw(~ismember(indRaw,ind3))= NaN;
        end
               
        params_ind3(i,:,m,1,1) = mean(v(~isnan(indRaw),:));
        params_ind3(i,:,m,2,1) = median(v(~isnan(indRaw),:));
        params_ind3(i,:,m,3,1) = std(v(~isnan(indRaw),:));
        params_ind3(i,:,m,4,1) = std(v(~isnan(indRaw),:))/sqrt(sum(~isnan(indRaw)));
        
        %combine ind2 and ind3
        tt = zeros(size(v,1),1);
        tt(ind2)=1;
        ind4 = ~isnan(indRaw) & tt==1;
                
        cbvals{i,m,1,1}= v(ind4,:);
        cbvals{i,m,2,1}= vdiff(ind4,:);
        params_ind4(i,:,m,1,1) = mean(v(ind4,:));
        params_ind4(i,:,m,2,1) = median(v(ind4,:));
        params_ind4(i,:,m,3,1) = std(v(ind4,:));
        params_ind4(i,:,m,4,1) = std(v(ind4,:))/sqrt(length(ind4));
        
        perc_incl(i,m,3:4) = [sum(~isnan(indRaw))/length(v) sum(ind4)/length(v)]; 
       
%         % for normalised data - needs to be adapted
%         if m < length(mask_name)
%             %zstand: zallm
%             if m == 1
%                 mallm = mean(v(ind2,:));
%                 stdallm = sqrt(var(v(ind2,:)));
%                 v1 = v; % to calculate zgm also for allm
%             end
%             vallm = NaN(length(v),6);
%             for c1 = 1:6 %for all columns
%             vallm(:,c1) = (v(:,c1)-mallm(c1))/stdallm(c1);
%             end
%             vdiffallm = [vallm(:,2)-vallm(:,1) vallm(:,3)-vallm(:,1) vallm(:,3)-vallm(:,2)...
%                 vallm(:,5)-vallm(:,4) vallm(:,6)-vallm(:,4) vallm(:,6)-vallm(:,5)];
%             cvals{i,m,1,2}= vallm(ind2,:);
%             cvals{i,m,2,2}= vdiffallm(ind2,:);
%             params_ind2(i,:,m,1,2) = mean(vallm(ind2,:));
%             params_ind2(i,:,m,2,2) = median(vallm(ind2,:));
%             params_ind2(i,:,m,3,2) = std(vallm(ind2,:));
%             params_ind2(i,:,m,4,2) = std(vallm(ind2,:))/sqrt(length(ind2));
%             %zstand: zgm        
%             if m == 2
%                 mgm = mean(v(ind2,:));
%                 stdgm = sqrt(var(v(ind2,:)));
%                 %calculate for allm (m==1)
%                 vgm = NaN(length(v1),6);
%                 for c1 = 1:6 %for all columns
%                     vgm(:,c1) = (v1(:,c1)-mgm(c1))/stdgm(c1);
%                 end
%                 vdiffgm = [vgm(:,2)-vgm(:,1) vgm(:,3)-vgm(:,1) vgm(:,3)-vgm(:,2)...
%                     vgm(:,5)-vgm(:,4) vgm(:,6)-vgm(:,4) vgm(:,6)-vgm(:,5)];
%                 cvals{i,m-1,1,3}= vgm(ind2,:);
%                 cvals{i,m-1,2,3}= vdiffgm(ind2,:);
%                 params_ind2(i,:,m-1,1,3) = mean(vgm(ind2,:));
%                 params_ind2(i,:,m-1,2,3) = median(vgm(ind2,:));
%                 params_ind2(i,:,m-1,3,3) = std(vgm(ind2,:));
%                 params_ind2(i,:,m-1,4,3) = std(vgm(ind2,:))/sqrt(length(ind2));
%             end
%             vgm = NaN(length(v),6);
%             if m >= 2
%                 for c1 = 1:6 %for all columns
%                 vgm(:,c1) = (v(:,c1)-mgm(c1))/stdgm(c1);
%                 end
%                 vdiffgm = [vgm(:,2)-vgm(:,1) vgm(:,3)-vgm(:,1) vgm(:,3)-vgm(:,2)...
%                     vgm(:,5)-vgm(:,4) vgm(:,6)-vgm(:,4) vgm(:,6)-vgm(:,5)];
%                 cvals{i,m,1,3}= vgm(ind1,:);
%                 cvals{i,m,2,3}= vdiffgm(ind2,:);
%                 params_ind2(i,:,m,1,3) = mean(vgm(ind2,:));
%                 params_ind2(i,:,m,2,3) = median(vgm(ind2,:));
%                  params_ind2(i,:,m,3,3) = std(vgm(ind2,:));
%                 params_ind2(i,:,m,4,3) = std(vgm(ind2,:))/sqrt(length(ind2));
%             end
%         end                       
    end
end

% %save
% save(fullfile(d_path,'sleepana_rawvals_cleanvals_bordervals_remorv.mat'),'vps','sleep','wake',...
%     'mask_name', 'mask_size', 'perc_incl','rawvals','cvals','cbvals','params_ind2','params_ind3','params_ind4','-v7.3');


%'analysis','paper2_sleep','res','diffusion',
%%%%%%%%%%%%%%%%CHANGE
save(fullfile(nii_path,'sleepana_rawvals_cleanvals_bordervals_remorv_MA.mat'),'vps','sleep','wake',...
    'mask_name', 'mask_size', 'perc_incl','rawvals','cvals','cbvals','params_ind2','params_ind3','params_ind4','-v7.3');




% %% compute distribution parameters and test for normal distribution
% warning('off','stats:lillietest:OutOfRangePLow');
% warning('off','stats:lillietest:OutOfRangePHigh');
% warning('off','stats:adtest:OutOfRangePLow');
% warning('off','stats:adtest:OutOfRangePHigh');
% params = NaN(length(vps),6,length(mask_name),6,2,3); % last dim 1=cols = raw values L0-2, C0-2; 2= cols = difference values L1-L0, L2-L0, L2-L1, C1-0, C2-0, C2-1
% %dim4: 1 mean
% % 2 med 
% % 3 mod
% % 4 iqr %interquartile range, mittlere 50%,  difference between the 75th and the 25th percentiles 
% % 5 skew %pos re, neg li
% % 6 kurt %wölbung, steilheit, outliers, 3 = normalvert, <3 - weniger outliers, treut gleichm, >3 mehr outliers, streuung wg wenigen extremen werten
% %dim5: single conditions L0;L1.../ differences L1-0, L2-0...;
% %dim6: rawvalues / zallm / zgm
% meas = {'mean', 'med','mod','iqr','skew','kurt'};
% nonorm_lillie = NaN(length(vps),6,length(mask_name),2,3); %dim4 = single/difference vals; dim5: raw, zallm, zgm
% nonorm_ad = NaN(length(vps),6,length(mask_name),2,3);
% 
% for i = 1:length(vps)
%     for m = 1:length(mask_name)
%         for i1 = 1:3 %for whichvals: raw/zallm/zgm        
%             for i2 = 1:2 %for single vals or difference vals
%                 %calculate parameters
%                 params(i,:,m,1,i2,i1) = mean(cvals{i,m,i2,i1});
%                 params(i,:,m,2,i2,i1) = median(cvals{i,m,i2,i1});
%                 params(i,:,m,3,i2,i1) = mode(cvals{i,m,i2,i1});
%                 params(i,:,m,4,i2,i1) = iqr(cvals{i,m,i2,i1});
%                 params(i,:,m,5,i2,i1) = skewness(cvals{i,m,i2,i1});
%                 params(i,:,m,6,i2,i1) = kurtosis(cvals{i,m,i2,i1});
%         
%                 %test for normal distribution           
%                 for c  = 1:length(cond)
%                     %     %Kolmogorov-Smirnov kstest(2) by default tests for standard normal distribution; is not accurate when cdf is estimated from data as in our case.
%                     %     %Jarque-Bera-test jbtest
%                     %Anderson-Darling-test
%                     [h,p,adstat,cv] = adtest(cvals{i,m,i2,i1}(:,c));
%                     if p <= 0.05
%                         fprintf('vp%02d %s %s %s Anderson-Darling p = %0.3f\n',...
%                             vps(i), mask_name{m}, whichval{i1}, cond{i2,c}, p);
%                         nonorm_ad(i,c,m,i2,i1) = 1;
%                     else
%                         nonorm_ad(i,c,m,i2,i1) = 0;
%                     end
%                     %Lilliefors test (mini strenger als adtest 98% übereinstimmung)
%                     [h,p,kstat,critval] = lillietest(cvals{i,m,i2,i1}(:,c)); 
%                     if p <= 0.05
%                         fprintf('vp%02d %s %s %s Lilliefors p = %0.3f\n',...
%                             vps(i), mask_name{m}, whichval{i1}, cond{i2,c}, p);
%                         nonorm_lillie(i,c,m,i2,i1) = 1;
%                     else
%                         nonorm_lillie(i,c,m,i2,i1) = 0;
%                     end
%                 end
%             end
%         end         
%     end
% end
% warning('on','stats:lillietest:OutOfRangePLow');
% warning('on','stats:lillietest:OutOfRangePHigh');
% warning('on','stats:adtest:OutOfRangePLow');
% warning('on','stats:adtest:OutOfRangePHigh');
% % 
% % save(fullfile(d_path,'rawvals_cleanvals_full_remorv_allmasks.mat'), ...
% %     'cond', 'd_path', 'mask_name', 'meas', 'nonorm_ad', 'nonorm_lillie',...
% %     'params', 'perc_incl', 'rawvals', 'sleep', 'vps', 'wake','whichval','cvals','params_ind2',...
% %     '-v7.3');
% 
% 
% 
% %% test parameters
% warning('off','stats:lillietest:OutOfRangePLow');
% warning('off','stats:lillietest:OutOfRangePHigh');
% singlediff = 1; %1 single cond values, 2 differences between conditions
% rawz = 3; %whichval: raw/zallm/zgm
% c1 = 4; %which cols 1=L0/L1-0; 2 = L1/L2-0 ...
% c2 = 5;
% for m=1:length(mask_name)
%     %in how many cases = vps is underlying distribution not normal?
% %     fprintf('\n%s %s %s %02d / %s %02d cases/vps not lilliefors normal\n',...
% %         mask_name{m}, whichval{rawz}, cond{singlediff,c1}, sum(nonorm_lillie(:,c1,m,singlediff,rawz)), cond{singlediff,c2}, sum(nonorm_lillie(:,c2,m,singlediff,rawz)));
% %     fprintf('%s %s %s %02d / %s %02d cases/vps not anderson-darling normal\n',...
% %         mask_name{m}, whichval{rawz}, cond{singlediff,c1}, sum(nonorm_ad(:,c1,m,singlediff,rawz)), cond{singlediff,c2}, sum(nonorm_ad(:,c2,m,singlediff,rawz)));
%     for ms = 1:2%1:length(meas)
% %         % test for normal distribution of distr parameters themselves
% %         [h,p,kstat,critval] = lillietest(params(:,c1,m,ms,singlediff,rawz));
% %             if p <= 0.05
% %                 fprintf('%s %s %s %s Lilliefors p = %0.3f\n',...
% %                     mask_name{m}, whichval{rawz}, meas{ms}, cond{singlediff,c1}, p);
% %             end
% %         [h,p,kstat,critval] = lillietest(params(:,c2,m,ms,singlediff));
% %         if p <= 0.05
% %             fprintf('%s %s %s %s Lilliefors p = %0.3f\n',...
% %                 mask_name{m}, whichval{rawz}, meas{ms}, cond{singlediff,c2}, p);
% %         end
%         %ttests    
% %         [h,p,ci,stats] = ttest(params(:,c1,m,ms,singlediff,rawz),params(:,c2,m,ms,singlediff,rawz));
% % %         if p < 0.1
% %             fprintf('%s %s %s ttest p = %0.3f %s %0.6f %s %0.6f\n',...
% %                 mask_name{m}, whichval{rawz}, meas{ms}, p, cond{singlediff,c1}, mean(params(:,c1,m,ms,singlediff,rawz)), cond{singlediff,c2}, mean(params(:,c2,m,ms,singlediff,rawz)));
% % %         end
%             [h,p,ci,stats] = ttest(params_ind2(:,c1,m,ms,rawz),params_ind2(:,c2,m,ms,rawz));
%             fprintf('%s %s %s ttest p = %0.3f %s %0.6f %s %0.6f\n',...
%                 mask_name{m}, whichval{rawz}, meas{ms}, p, cond{1,c1}, mean(params_ind2(:,c1,m,ms,rawz)), cond{1,c2}, mean(params_ind2(:,c2,m,ms,rawz)));
%     end
% end
% warning('on','stats:lillietest:OutOfRangePLow');
% warning('on','stats:lillietest:OutOfRangePHigh');


%% 
% %% PLOT L1 L0
% 
% for i = 1:4%length(vps)
%     figure
%     for m = 1:length(mask_name)
%         v = rawvals{i,m};
%         ax{m} =  subplot(3,4,m); %,'XLim',[0 6]);
%         hold on
%         if m < length(mask_name)
%             ind = find(v(:,1)>=0.0005 & v(:,1)<=0.001 & v(:,2)>=0.0005 & v(:,2)<=0.001);
%             xlim([0.0005 0.001]);
%             %    ind2 = find(v(:,1)<0.0005 | v(:,1)>0.001 | v(:,2)<0.0005 | v(:,2)>0.001);
%         else
%             ind = find(v(:,1)>=0.001 & v(:,1)<=0.0035 & v(:,2)>=0.001 & v(:,2)<=0.0035);
%             xlim([0.001 0.0035]);
%         end
%         
%         h1 = histogram(v(ind,1), 'FaceColor', [0 0 1]);
%         h2 = histogram(v(ind,2), 'FaceColor', [1 0 0]);
%         % h1 = histogram(v(ind,5)-v(ind,4), 'FaceColor', [0.5 0.5 0.5]); %C1-C0
%         % h2 = histogram(v(ind,2)-v(ind,1), 'FaceColor', [1 0 0]); %L1-L0
%         h1.BinWidth = min(h1.BinWidth, h2.BinWidth);
%         h2.BinWidth = min(h1.BinWidth, h2.BinWidth);
%         
%         set(get(gca,'YLabel'),'String','freq')
%         set(get(gca,'XLabel'),'String','MD raw')
%         title([mask_name{m} sprintf(' %.1f%% of %d vox', length(ind)/length(rawvals{i,m})*100, length(rawvals{i,m}))]);
%         
%     end    
%     hl = legend('L0', 'L1');
%     hl.Location = 'east';
%     suptitle(sprintf('VP %02d', vps(i)));
% end
% 
% 
% %% PLOT L1-L0, C1-C0
% 
% for i = 1:4%length(vps)
%     figure
%     for m = 1:length(mask_name)
%         v = rawvals{i,m};
%         ax{m} =  subplot(3,4,m); %,'XLim',[0 6]);
%         hold on
%         if m < length(mask_name)
%             ind = find(v(:,1)>=0.0005 & v(:,1)<=0.001 & v(:,2)>=0.0005 & v(:,2)<=0.001 & ...
%                 v(:,4)>=0.0005 & v(:,4)<=0.001 & v(:,5)>=0.0005 & v(:,5)<=0.001 & ...
%                 abs(v(:,2)-v(:,1))<=0.0002 & abs(v(:,5)-v(:,4))<=0.0002);
%             xlim([-0.0002 0.0002]);
%             %    ind2 = find(v(:,1)<0.0005 | v(:,1)>0.001 | v(:,2)<0.0005 | v(:,2)>0.001);
%         else
%             ind = find(v(:,1)>=0.001 & v(:,1)<=0.0035 & v(:,2)>=0.001 & v(:,2)<=0.0035 & ...
%                 v(:,4)>=0.001 & v(:,4)<=0.0035 & v(:,5)>=0.001 & v(:,5)<=0.0035);
% %             xlim([0.001 0.0035]);
%         end
%         
%         h1 = histogram(v(ind,5)-v(ind,4), 'FaceColor', [0.1 0.1 0.1]); %C1-C0
%         h2 = histogram(v(ind,2)-v(ind,1), 'FaceColor', [1 0 0]); %L1-L0
%         h1.BinWidth = min(h1.BinWidth, h2.BinWidth);
%         h2.BinWidth = min(h1.BinWidth, h2.BinWidth);
%         
%         set(get(gca,'YLabel'),'String','freq')
%         set(get(gca,'XLabel'),'String','MD voxelwise diff')
%         title([mask_name{m} sprintf(' %.1f%% of %d vox', length(ind)/length(rawvals{i,m})*100, length(rawvals{i,m}))]);
%         
%     end       
%     hl = legend('C1-C0', 'L1-L0');
%     hl.Location = 'east';
%     suptitle(sprintf('VP %02d', vps(i)));
% end
% 
% %% PLOT L1 L0 all vp together
%     figure
%     for m = 1:length(mask_name)
%         if mod(m,4)== 1
%             figure            
%         end
%         if mod(m,4)== 0
%             n=4;
%         else
%             n = mod(m,4);
%         end
%         
%         ax{m} =  subplot(2,2,n); %,'XLim',[0 6]);
%         hold on
%         for i = 1:length(vps)
%             if i == 1
%                 v = rawvals{i,m};
%             else
%                 v = [v; rawvals{i,m}];
%             end                          
%         end
%         if m < length(mask_name)
%             ind = find(v(:,1)>=0.0005 & v(:,1)<=0.001 & v(:,2)>=0.0005 & v(:,2)<=0.001);
%             xlim([0.0005 0.001]);
%             %    ind2 = find(v(:,1)<0.0005 | v(:,1)>0.001 | v(:,2)<0.0005 | v(:,2)>0.001);
%         else
%             ind = find(v(:,1)>=0.001 & v(:,1)<=0.0035 & v(:,2)>=0.001 & v(:,2)<=0.0035);
%             xlim([0.001 0.0035]);
%         end
%         h1 = histogram(v(ind,1), 'FaceColor', [0 0 1]);
%         h2 = histogram(v(ind,2), 'FaceColor', [1 0 0]);
%         % h1 = histogram(v(ind,5)-v(ind,4), 'FaceColor', [0.5 0.5 0.5]); %C1-C0
%         % h2 = histogram(v(ind,2)-v(ind,1), 'FaceColor', [1 0 0]); %L1-L0
%         h1.BinWidth = min(h1.BinWidth, h2.BinWidth);
%         h2.BinWidth = min(h1.BinWidth, h2.BinWidth);
%         
%         set(get(gca,'YLabel'),'String','freq')
%         set(get(gca,'XLabel'),'String','MD raw')
%         title([mask_name{m} sprintf(' %.1f%% of %d vox', length(ind)/length(v)*100, length(v))]);
%         if mod(m,4)==0
%         hl = legend('L0', 'L1');
%         hl.Location = 'east';
%         suptitle('all REMORV');
%         end
%     end    
%     
% 
% %% PLOT L1-L0, C1-C0 all vp
% figure
% for m = 1:length(mask_name)
%     if mod(m,4)== 1
%         figure            
%     end
%     if mod(m,4)== 0
%         n=4;
%     else
%         n = mod(m,4);
%     end        
%     ax{m} =  subplot(2,2,n); %,'XLim',[0 6]);
%     hold on
%     for i = 1:length(vps)
%         if i == 1
%             v = rawvals{i,m};
%         else
%             v = [v; rawvals{i,m}];
%         end                          
%     end
%     if m < length(mask_name)
%         ind = find(v(:,1)>=0.0005 & v(:,1)<=0.001 & v(:,2)>=0.0005 & v(:,2)<=0.001 & ...
%             v(:,4)>=0.0005 & v(:,4)<=0.001 & v(:,5)>=0.0005 & v(:,5)<=0.001 & ...
%             abs(v(:,2)-v(:,1))<=0.0002 & abs(v(:,5)-v(:,4))<=0.0002);
%         xlim([-0.0002 0.0002]);
%         %    ind2 = find(v(:,1)<0.0005 | v(:,1)>0.001 | v(:,2)<0.0005 | v(:,2)>0.001);
%     else
%         ind = find(v(:,1)>=0.001 & v(:,1)<=0.0035 & v(:,2)>=0.001 & v(:,2)<=0.0035 & ...
%                 v(:,4)>=0.001 & v(:,4)<=0.0035 & v(:,5)>=0.001 & v(:,5)<=0.0035);
% %             xlim([0.001 0.0035]);
%     end
% 
%     h1 = histogram(v(ind,5)-v(ind,4), 'FaceColor', [0.1 0.1 0.1]); %C1-C0
%     h2 = histogram(v(ind,2)-v(ind,1), 'FaceColor', [1 0 0]); %L1-L0
%     h1.BinWidth = min(h1.BinWidth, h2.BinWidth);
%     h2.BinWidth = min(h1.BinWidth, h2.BinWidth);
% 
%     set(get(gca,'YLabel'),'String','freq')
%     set(get(gca,'XLabel'),'String','MD voxelwise diff')
%     title([mask_name{m} sprintf(' %.1f%% of %d vox', length(ind)/length(v)*100, length(v))]);
%     if mod(m,4)==0
%     hl = legend('C1-C0', 'L1-L0');
%     hl.Location = 'east';
%     suptitle('all REMORV');
%     end
% end       
    