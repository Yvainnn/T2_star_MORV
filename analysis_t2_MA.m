%% load data
clear all
fs_path= 'H:\Unix_Folders\niftiE\freesurfer_mprage\subjects';
nii_path = 'H:\Unix_Folders\niftiE\niftiE';
d_path = 'H:\Unix_Folders\MORV\niftiD\niftiD\'; %mask path is different from t2_singlete nii files
%H:\Unix_Folders\MORV\niftiD\diffusion_analysis\nativespace_2mm T2star
% NiftiD H:\Unix_Folders\MORV\niftiD\niftiD

% fs output brain region, use this insted of mask_name
%brain_reg_rh = {'rh_bankssts_volume', 'rh_caudalanteriorcingulate_volume', 'rh_caudalmiddlefrontal_volume', 'rh_cuneus_volume', 'rh_entorhinal_volume', 'rh_fusiform_volume', 'rh_inferiorparietal_volume', 'rh_inferiortemporal_volume', 'rh_isthmuscingulate_volume', 'rh_lateraloccipital_volume', 'rh_lateralorbitofrontal_volume', 'rh_lingual_volume', 'rh_medialorbitofrontal_volume', 'rh_middletemporal_volume', 'rh_parahippocampal_volume', 'rh_paracentral_volume', 'rh_parsopercularis_volume', 'rh_parsorbitalis_volume', 'rh_parstriangularis_volume', 'rh_pericalcarine_volume', 'rh_postcentral_volume', 'rh_posteriorcingulate_volume', 'rh_precentral_volume', 'rh_precuneus_volume', 'rh_rostralanteriorcingulate_volume', 'rh_rostralmiddlefrontal_volume', 'rh_superiorfrontal_volume', 'rh_superiorparietal_volume', 'rh_superiortemporal_volume', 'rh_supramarginal_volume', 'rh_frontalpole_volume', 'rh_temporalpole_volume', 'rh_transversetemporal_volume', 'rh_insula_volume'}; 
%brain_reg_lh = {'lh_bankssts_volume', 'lh_caudalanteriorcingulate_volume', 'lh_caudalmiddlefrontal_volume', 'lh_cuneus_volume', 'lh_entorhinal_volume', 'lh_fusiform_volume', 'lh_inferiorparietal_volume', 'lh_inferiortemporal_volume', 'lh_isthmuscingulate_volume', 'lh_lateraloccipital_volume', 'lh_lateralorbitofrontal_volume', 'lh_lingual_volume', 'lh_medialorbitofrontal_volume', 'lh_middletemporal_volume', 'lh_parahippocampal_volume', 'lh_paracentral_volume', 'lh_parsopercularis_volume', 'lh_parsorbitalis_volume', 'lh_parstriangularis_volume', 'lh_pericalcarine_volume', 'lh_postcentral_volume', 'lh_posteriorcingulate_volume', 'lh_precentral_volume', 'lh_precuneus_volume', 'lh_rostralanteriorcingulate_volume', 'lh_rostralmiddlefrontal_volume', 'lh_superiorfrontal_volume', 'lh_superiorparietal_volume', 'lh_superiortemporal_volume', 'lh_supramarginal_volume', 'lh_frontalpole_volume', 'lh_temporalpole_volume', 'lh_transversetemporal_volume', 'lh_insula_volume'}; 
%brain_aseg= {'Left-Lateral-Ventricle', 'Left-Inf-Lat-Vent','Left-Cerebellum-White-Matter', 'Left-Cerebellum-Cortex', 'Left-Thalamus-Proper', 'Left-Caudate', 'Left-Putamen', 'Left-Pallidum', '3rd-Ventricle', '4th-Ventricle', 'Brain-Stem', 'Left-Hippocampus', 'Left-Amygdala', 'CSF', 'Left-Accumbens-area', 'Left-VentralDC', 'Left-vessel', 'Left-choroid-plexus', 'Right-Lateral-Ventricle', 'Right-Inf-Lat-Vent', 'Right-Cerebellum-White-Matter', 'Right-Cerebellum-Cortex', 'Right-Thalamus-Proper', 'Right-Caudate', 'Right-Putamen', 'Right-Pallidum', 'Right-Hippocampus', 'Right-Amygdala', 'Right-Accumbens-area', 'Right-VentralDC', 'Right-vessel', 'Right-choroid-plexus', '5th-Ventricle', 'WM-hypointensities', 'Left-WM-hypointensities	Right-WM-hypointensities', 'non-WM-hypointensities', 'Left-non-WM-hypointensities', 'Right-non-WM-hypointensities', 'Optic-Chiasm', 'CC_Posterior', 'CC_Mid_Posterior', 'CC_Central', 'CC_Mid_Anterior', 'CC_Anterior'};

% load(fullfile(d_path,'rawvals_cleanvals_full_remorv_allmasks.mat'))
% load(fullfile(d_path,'sleepana_rawvals_cleanvals_remorv.mat'))
load(fullfile(nii_path,'sleepana_rawvals_cleanvals_bordervals_remorv_MA.mat'))
%cvals: 33 vps x 15 masks x 2 cond singlevsdiff x 3 raw/zallm/zgm
%params_ind2: 33vps x 6cond x 16 masks x 4 mean/median/sd/sem x 3 raw/zallm/zgm -
%with implausible vals excluded
%params_ind3: 33vps x 6cond x 16 masks x 4 mean/median/sd/sem x 3 raw[/zallm/zgm] -
%with border voxels excluded
%params_ind4: with both excluded
% vps n=33, sleep then wake

load(fullfile(anapath,'behavres','perf.mat'),'performance');
%right vps and right order
performance = performance(ismember(performance(:,9),vps),:);
performance = [performance(17:end,:);performance(1:16,:)];

% Change brain regions mask, fs output brain region, use this insted of mask_name
%brain_reg_rh = {'rh_bankssts_volume', 'rh_caudalanteriorcingulate_volume', 'rh_caudalmiddlefrontal_volume', 'rh_cuneus_volume', 'rh_entorhinal_volume', 'rh_fusiform_volume', 'rh_inferiorparietal_volume', 'rh_inferiortemporal_volume', 'rh_isthmuscingulate_volume', 'rh_lateraloccipital_volume', 'rh_lateralorbitofrontal_volume', 'rh_lingual_volume', 'rh_medialorbitofrontal_volume', 'rh_middletemporal_volume', 'rh_parahippocampal_volume', 'rh_paracentral_volume', 'rh_parsopercularis_volume', 'rh_parsorbitalis_volume', 'rh_parstriangularis_volume', 'rh_pericalcarine_volume', 'rh_postcentral_volume', 'rh_posteriorcingulate_volume', 'rh_precentral_volume', 'rh_precuneus_volume', 'rh_rostralanteriorcingulate_volume', 'rh_rostralmiddlefrontal_volume', 'rh_superiorfrontal_volume', 'rh_superiorparietal_volume', 'rh_superiortemporal_volume', 'rh_supramarginal_volume', 'rh_frontalpole_volume', 'rh_temporalpole_volume', 'rh_transversetemporal_volume', 'rh_insula_volume'}; 
%brain_reg_lh = {'lh_bankssts_volume', 'lh_caudalanteriorcingulate_volume', 'lh_caudalmiddlefrontal_volume', 'lh_cuneus_volume', 'lh_entorhinal_volume', 'lh_fusiform_volume', 'lh_inferiorparietal_volume', 'lh_inferiortemporal_volume', 'lh_isthmuscingulate_volume', 'lh_lateraloccipital_volume', 'lh_lateralorbitofrontal_volume', 'lh_lingual_volume', 'lh_medialorbitofrontal_volume', 'lh_middletemporal_volume', 'lh_parahippocampal_volume', 'lh_paracentral_volume', 'lh_parsopercularis_volume', 'lh_parsorbitalis_volume', 'lh_parstriangularis_volume', 'lh_pericalcarine_volume', 'lh_postcentral_volume', 'lh_posteriorcingulate_volume', 'lh_precentral_volume', 'lh_precuneus_volume', 'lh_rostralanteriorcingulate_volume', 'lh_rostralmiddlefrontal_volume', 'lh_superiorfrontal_volume', 'lh_superiorparietal_volume', 'lh_superiortemporal_volume', 'lh_supramarginal_volume', 'lh_frontalpole_volume', 'lh_temporalpole_volume', 'lh_transversetemporal_volume', 'lh_insula_volume'}; 
%brain_aseg= {'Left-Lateral-Ventricle', 'Left-Inf-Lat-Vent','Left-Cerebellum-White-Matter', 'Left-Cerebellum-Cortex', 'Left-Thalamus-Proper', 'Left-Caudate', 'Left-Putamen', 'Left-Pallidum', '3rd-Ventricle', '4th-Ventricle', 'Brain-Stem', 'Left-Hippocampus', 'Left-Amygdala', 'CSF', 'Left-Accumbens-area', 'Left-VentralDC', 'Left-vessel', 'Left-choroid-plexus', 'Right-Lateral-Ventricle', 'Right-Inf-Lat-Vent', 'Right-Cerebellum-White-Matter', 'Right-Cerebellum-Cortex', 'Right-Thalamus-Proper', 'Right-Caudate', 'Right-Putamen', 'Right-Pallidum', 'Right-Hippocampus', 'Right-Amygdala', 'Right-Accumbens-area', 'Right-VentralDC', 'Right-vessel', 'Right-choroid-plexus', '5th-Ventricle', 'WM-hypointensities', 'Left-WM-hypointensities	Right-WM-hypointensities', 'non-WM-hypointensities', 'Left-non-WM-hypointensities', 'Right-non-WM-hypointensities', 'Optic-Chiasm', 'CC_Posterior', 'CC_Mid_Posterior', 'CC_Central', 'CC_Mid_Anterior', 'CC_Anterior'};
gm_wm_m= {'BS_volumewithoutvent', 'TT_cortical_wm', 'TT_cortical_gm', 'TT_subcortical_gm', 'TT_intracranial_gm', 'TT_gm'};

%prepare plots
c1 = uint8([128 128 128]);
% c1 = uint8([50 65 75]);
% c2 = uint8([165 30 55]);
c2 = uint8([200 0 0]);
% c3 = uint8([180 77 80]);
% c3 = uint8([0 105 170]);
% c3 = uint8([65 90 140]);
c3 = uint8([0 0 128]);

group = ones(33,1); group(18:end)=2; %1 sleep 2 wake

%% morning vs evening
for k=14:length(gm_wm_m) %--> run with VBM
%setup mask
mask = 1;
val = 1; %mean/median/std %%%%we look at 1 and 3!
norm = 1; %raw
vdata = cell2mat(gm_wm_m(1,k));
%vdata = params_ind4(:,:,mask,val,norm);

% vtimes = time_from_wake;



%control morning, control evening, learn morning, learn evening
v = [vdata(1:17,6) vdata(1:17,5) vdata(1:17,3) vdata(1:17,2); vdata(18:end,4) vdata(18:end,6) vdata(18:end,1) vdata(18:end,3)];

%plot
figure
hold on
means = [nanmean(v)]; 
sems = [nanstd(v)/sqrt(length(v))];
model_series = [means(1:2); means(3:4)];
model_error = [sems(1:2);sems(3:4)];
h = bar(model_series);
set(h,'BarWidth',0.7);
h(1).FaceColor = c1-40;
h(2).FaceColor = c1+40;
% h(3).FaceColor = c3;
h(1).EdgeColor = 'none';
h(2).EdgeColor = 'none';
% h(3).EdgeColor = 'none';
set(gca, 'Box', 'off');
set(gca,'XTick',[1:2]);
set(gca,'XTicklabel',{'control', 'learn'})
set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
set(get(gca,'YLabel'),'String','t2*')
% region_name= strjoin(gm_wm_m(k));
% set(get(gca,'YLabel'),'String',sprintf('T2* - %s',region_name))

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    e(i) =errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','linewidth',3);
    %for significance bars
    bars(:,i) = x;   
end
% ylim([0.00076 0.00078])
% ylim([0.000057 0.000063])

% stats: 2 mor-eve x 2 l-c x 2 s-w
% main effect within and interaction
% create rm model
groupstr = cellstr(num2str(group));
between = table(groupstr,v(:,1),v(:,2),v(:,3),v(:,4), 'VariableNames',{'group','cm','ce','lm','le'}); % input data table
within = table([1 1 -1 -1]',[1 -1 1 -1]','VariableNames',{'condition','tod'});
modelspec = 'cm-le~group'; %formula for model spec response vars~predictor; %x1,x3 - variable x1 and x3, x1-x3 vars x1 to x3
rm = fitrm(between,modelspec,'WithinDesign',within,'WithinModel','condition*tod'); %repeated measures model object ;separate means default - computes separate mean for each group
%ranova
[withintbl,~,~,~] = ranova(rm,'WithinModel','condition*tod'); %A = between-subjects model; C = within-subjects model, D = hypothesis = 0;
pTod = table2array(withintbl('(Intercept):tod','pValue'));
pTodCond = table2array(withintbl('(Intercept):condition:tod','pValue'));
p3IA = table2array(withintbl('group:condition:tod','pValue'));
%anova
betweentbl = anova(rm);
pGroup = betweentbl.pValue(betweentbl.Between=='group');
%ttests for tod per cond
%paired sample across tps
[~,pTc,~,statsTc] = ttest(v(:,1),v(:,2));
[~,pTl,~,statsTl] = ttest(v(:,3),v(:,4));

%significance bars
xtext = h.XData;
factorLines = 0.0000005;
factorStars = 0.0000006;
factorIA = 0.0000005;
signfontsize = 40;
signfontweight = 'bold';
%x-axis groups: l c across tps
pvals = [pTc pTl];
for igroup = 1:2
    sign = sb_signsign(pvals(igroup));
    if ~isempty(sign)
        % up or down
        tt = abs(model_series(igroup,:))+model_error(igroup,:);
        maxidx = tt==max(tt);
        maxval= model_series(igroup,maxidx);
        if maxval > 0            
            ty = maxval+model_error(igroup,maxidx);
            tfLines = factorLines;
            tfStars = factorStars;
        else
            ty = maxval-model_error(igroup,maxidx);
            tfLines = -factorLines;
            tfStars = -factorStars;
        end
        text(xtext(igroup),ty+tfStars,sign,'FontSize', signfontsize,'FontWeight',signfontweight,'HorizontalAlignment','center');
        line(bars(igroup,:),[ty ty]+tfLines,'LineWidth',3,'Color','k');       
    end
end
%main effect of tod
sign = sb_signsign(pTod);
if ~isempty(sign)
    % up or down
    tt = abs(model_series(:))+model_error(:);
    maxidx = tt==max(tt);
    maxval= model_series(maxidx);
    if maxval > 0
        ty = maxval+model_error(maxidx);
        tfIA = factorIA;
        tfLines = factorLines;
        tfStars = factorStars;
    else
        ty = maxval-model_error(maxidx);
        tfIA = -factorIA;
        tfLines = -factorLines;
        tfStars = -factorStars;
    end
        text(mean(xtext),ty+tfIA+tfStars,['tod ' sign],'FontSize', signfontsize,'FontWeight',signfontweight,'HorizontalAlignment','center');
end

% 
% hl = legend('morning','evening');
% hl.Location = 'northeast';
% hl.Box = 'off';
% hl.FontSize = 40;
% hl.FontWeight = 'bold';

fprintf('\nmain effect tod p=%.3f\nIA tod*cond p=%.3f\nIA tod*cond*group p=%.3f\n',pTod,pTodCond,p3IA);

end 
%% wake circadian linear decrease
% setup mask
mask = 1;
val = 1; %mean/median/std %%%%we look at 1 and 3!
norm = 1; %raw
% vdata = params_ind2(:,:,mask,val,norm);
vdata = cell2mat(gm_wm_m(1,k));

% vtimes = time_from_wake;

cmap = colormap(jet(12));
figure
hold on
v = [vdata(18:end,1:3); vdata(18:end,4:6)];
%plot
figure
hold on
h1 = bar(1,mean(v(:,1)),'FaceColor',cmap(end,:),'EdgeColor','none','BarWidth',0.7);
h2 = bar(2,mean(v(:,2)),'FaceColor',cmap(end-1,:),'EdgeColor','none','BarWidth',0.7);
h3 = bar(4,mean(v(:,3)),'FaceColor',cmap(end-2,:),'EdgeColor','none','BarWidth',0.7);
means = nanmean(v); 
sems = nanstd(v)/sqrt(length(v));
he = errorbar([1 2 4],means,sems, 'k', 'linestyle', 'none','linewidth',3);
set(gca, 'Box', 'off');
set(gca,'XTick',[1 2 4]);
set(gca,'XTicklabel',{'9:00','12:00','22:00'})
set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
set(get(gca,'YLabel'),'String','T2*')
set(get(gca,'XLabel'),'String','time of day')
ylim([2500 3000])
% ylim([0.000057 0.000062])


%stats
% main effect within and interaction
%create rm model
dummy = ones(length(v)/2,1);
dummystr = cellstr(num2str(dummy));
between = table(v(1:15,1),v(1:15,2),v(1:15,3),v(16:end,1),v(16:end,2),v(16:end,3),... %16 if u are working with VBM and 15 if u wanna process T2* maps (sbj 45 miss)
'VariableNames',{'lt0','lt1','lt2','ct0','ct1','ct2'}); % input data table
% within = table([1 1 1 -1 -1 -1]', [1 0 -1 1 0 -1]','VariableNames',{'condition','timepoints'}); %Design for within-subject factors, rows = nr of response vars, cols = within factors, cells = labels
within = table({'L';'L';'L';'C';'C';'C'}, {'T0';'T1';'T2';'T0';'T1';'T2'},'VariableNames',{'condition','timepoints'}); %Design for within-subject factors, rows = nr of response vars, cols = within factors, cells = labels
modelspec = 'lt0-ct2~1'; %formula for model spec response vars~predictor; %x1,x3 - variable x1 and x3, x1-x3 vars x1 to x3
rm = fitrm(between,modelspec,'WithinDesign',within,'WithinModel','condition*timepoints'); %repeated measures model object ;separate means default - computes separate mean for each group
%ranova
[withintbl,~,~,~] = ranova(rm,'WithinModel','condition*timepoints'); %A = between-subjects model; C = within-subjects model, D = hypothesis = 0;
pTp = table2array(withintbl('(Intercept):timepoints','pValue'));
% pIA = table2array(withintbl('cond:timepoints','pValue'));
pIA = table2array(withintbl('(Intercept):condition:timepoints','pValue'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lincontrast,~,~,~] = ranova(rm,'WithinModel',[1 0 -1 1 0 -1]'); %A = between-subjects model; C = within-subjects model, D = hypothesis = 0;
pLin = table2array(lincontrast('(Intercept):Time','pValue'));
% contrasts = anova(rm,'WithinModel','orthogonalcontrasts'); %A = between-subjects model; C = within-subjects model, D = hypothesis = 0;

% main effect between
betweentbl = anova(rm);
pCond = betweentbl.pValue(betweentbl.Between=='cond');

% significance bars
factorLines = 1.0022;
factorStars = 1.0024;
signfontsize = 40;
signfontweight = 'bold';
%linear contrast
    sign = sb_signsign(pLin);
    line([1 4],[means(1)+sems(1) means(3)+sems(3)]*1.0022,'LineWidth',3,'Color','k');
    text(2.5,(means(1)+sems(1)+means(3)+sems(3))/2*1.0024,sign,'FontSize', signfontsize,'FontWeight',signfontweight,'HorizontalAlignment','center');

    fprintf('\n\nmain effect tp p=%.3f\nIA tp*cond p=%.3f\nIA pCond p=%.3f\n',pTp,pIA,pCond);
%% sleep quadratic increase

cmap = colormap(jet(12));
v = [vdata(1:17,1:3); vdata(1:17,4:6)]; % 32sbj,instead of 33
%plot
figure
hold on
bar(1,mean(v(:,1)),'FaceColor',cmap(1,:),'EdgeColor','none','BarWidth',0.7);
bar(2,mean(v(:,2)),'FaceColor',cmap(2,:),'EdgeColor','none','BarWidth',0.7);
bar(4,mean(v(:,3)),'FaceColor',cmap(3,:),'EdgeColor','none','BarWidth',0.7);
means = nanmean(v); 
sems = nanstd(v)/sqrt(length(v));
errorbar([1 2 4],means,sems, 'k', 'linestyle', 'none','linewidth',3);
set(gca, 'Box', 'off');
set(gca,'XTick',[1 2 4]);
set(gca,'XTicklabel',{'20:00','23:00','9:00'})
set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
set(get(gca,'YLabel'),'String','MD raw')
set(get(gca,'XLabel'),'String','time of day')
ylim([2500 3000])
line([1 2 4],[means(1)+sems(1) means(2)+sems(2) means(3)+sems(3)]*1.0022,'LineWidth',3,'Color','k');
text(2.5,(means(1)+sems(1)+means(3)+sems(3))/2*1.0005,'***','FontSize', 40,'FontWeight','bold','HorizontalAlignment','center');

% ylim([0.000059 0.000063])
% line([1 2 4],[means(1)+sems(1) means(2)+sems(2) means(3)+sems(3)]*1.0022,'LineWidth',3,'Color','k');
% text(2.5, [means(2)+sems(2)]*1.0005,'*','FontSize', 40,'FontWeight','bold','HorizontalAlignment','center');

%% significance bars
xtext = h.XData;
factorLines = 5;
factorStars = 6;
factorIA = 5;
signfontsize = 40;
signfontweight = 'bold';
%x-axis groups: wake sleep across tps
%group 1: wake s1-s2
pvals = [pTw pTs];
for igroup = 1:2
    sign = sb_signsign(pvals(igroup));
    if ~isempty(sign)
        % up or down
        tt = abs(model_series(igroup,:))+model_error(igroup,:);
        maxidx = tt==max(tt);
        maxval= model_series(igroup,maxidx);
        if maxval > 0            
            ty = maxval+model_error(igroup,maxidx);
            tfLines = factorLines;
            tfStars = factorStars;
        else
            ty = maxval-model_error(igroup,maxidx);
            tfLines = -factorLines;
            tfStars = -factorStars;
        end
        text(xtext(igroup),ty+tfStars,sign,'FontSize', signfontsize,'FontWeight',signfontweight,'HorizontalAlignment','center');
        line(bars(igroup,:),[ty ty]+tfLines,'LineWidth',3,'Color','k');       
    end
end
% similar bars: s1/s2 across groups
pvals = [pTsess1 pTsess2];
for igroup = 1:2
    sign = sb_signsign(pvals(igroup));
    if ~isempty(sign)
        % up or down
        tt = abs(model_series(:,igroup))+model_error(:,igroup);
        maxidx = tt==max(tt);
        maxval= model_series(maxidx,igroup);
        if maxval > 0
            ty = maxval+model_error(maxidx,igroup);
            tfLines = factorLines;
            tfStars = factorStars;
        else
            ty = maxval-model_error(maxidx,igroup);
            tfLines = -factorLines;
            tfStars = -factorStars;
        end
            text(mean(bars(:,igroup)),ty+tfStars,sign,'FontSize', signfontsize,'FontWeight',signfontweight,'HorizontalAlignment','center');
            line(bars(:,igroup),[ty ty]+tfLines,'LineWidth',3,'Color','k');        
    end
end
% IA
sign = sb_signsign(pIA);
if ~isempty(sign)
    % up or down
    tt = abs(model_series(:))+model_error(:);
    maxidx = tt==max(tt);
    maxval= model_series(maxidx);
    if maxval > 0
        ty = maxval+model_error(maxidx);
        tfIA = factorIA;
        tfLines = factorLines;
        tfStars = factorStars;
    else
        ty = maxval-model_error(maxidx);
        tfIA = -factorIA;
        tfLines = -factorLines;
        tfStars = -factorStars;
    end
        for b= 1:2        
        line(bars(b,:),[ty ty]+tfLines,'LineWidth',3,'Color','k');
        line([xtext(b) xtext(b)],[ty ty+tfIA]+tfLines,'LineWidth',3,'Color','k');
        end
        text(mean(xtext),ty+tfIA+tfStars,sign,'FontSize', signfontsize,'FontWeight',signfontweight,'HorizontalAlignment','center');
        line(xtext,[ty ty]+tfLines+tfIA,'LineWidth',3,'Color','k');
end


%% rm ANOVA
% main effect within and interaction
%create rm model
group = ones(32,1); group(18:end)=2; %1 sleep 2 wake

groupstr = cellstr(num2str(group));
between = table(groupstr,v(:,1),v(:,2),...
'VariableNames',{'group','sess1','sess2'}); % input data table
within = table([1 2]','VariableNames',{'sessions'}); %Design for within-subject factors, rows = nr of response vars, cols = within factors, cells = labels
modelspec = 'sess1-sess2~group'; %formula for model spec response vars~predictor; %x1,x3 - variable x1 and x3, x1-x3 vars x1 to x3
rm = fitrm(between,modelspec,'WithinDesign',within,'WithinModel','separatemeans'); %repeated measures model object ;separate means default - computes separate mean for each group
%mauchly's test: only rm, sphericity is met if equal variances of differences between all possible pairs of groups are equal --> only if more than 2 levels
mauchlytbl = mauchly(rm);
%ranova
[withintbl,~,~,~] = ranova(rm,'WithinModel','separatemeans'); %A = between-subjects model; C = within-subjects model, D = hypothesis = 0;
pSessions = table2array(withintbl('(Intercept):sessions','pValue'));
pIA = table2array(withintbl('group:sessions','pValue'));

% main effect between
%Levene test, separately for all within variables, needs also to be done separately for each var;
[pLev1,statsLev1] = vartestn(v(:,1),groupstr,'TestType','LeveneAbsolute','Display','off'); 
[pLev2,statsLev2] = vartestn(v(:,2),groupstr,'TestType','LeveneAbsolute','Display','off'); 
%anova
betweentbl = anova(rm);
pGroup = betweentbl.pValue(betweentbl.Between=='group');
%ttests
%paired sample across tps
[~,pTw,~,statsTw] = ttest(v(group==1,1),v(group==1,2));
[~,pTs,~,statsTs] = ttest(v(group==2,1),v(group==2,2));
%independent samples
[~,pTsess1,~,statsTsess1] = ttest2(v(group==1,1),v(group==2,1));
[~,pTsess2,~,statsTsess2] = ttest2(v(group==1,2),v(group==2,2));
% %descriptives per group
% descrtbl = grpstats(rm,'group');
% %test normality of Residuals with Shapiro Wilk - per tp and group?
% [ypred,~] = predict(rm,between,'WithinDesign', within, 'WithinModel','orthogonalcontrasts'); %default model 'separatemeans,'orthogonalcontrasts' if within design only contains one single numeric fact
% resid = v-ypred;
% [~,pSW1,statSW1] = swtest(resid(group==1,1));
% [~,pSW2,statSW2] = swtest(resid(group==2,1));
% [~,pSW3,statSW3] = swtest(resid(group==1,2));
% [~,pSW4,statSW4] = swtest(resid(group==2,2));



%% plots OHBM 2019 poster

clear all
d_path = 'D:\Data\MORV\nifti\diffusion_analysis\nativespace_2mm';
anapath = 'D:\Data\MORV\analysis';
% load(fullfile(d_path,'rawvals_cleanvals_full_remorv_allmasks.mat'))
load(fullfile(d_path,'sleepana_rawvals_cleanvals_remorv.mat'))
load(fullfile(anapath,'totalsleeptime_remorv.mat'))
%cvals: 33 vps x 15 masks x 2 cond singlevsdiff x 3 raw/zallm/zgm
%params_ind2: 33vps x 6cond x 16 masks x 2 mean/median x 3 raw/zallm/zgm
% vps n=33, sleep then wake
load('D:\Data\MORV\analysis\scripts\times.mat')

%IDS n=41, sleep then wake
%changes times IDs to diff vps
vpind = find(ismember(IDs,vps));
times = times(vpind,:); times24 = times24(vpind,:); time_from_up = time_from_up(vpind,:); time_from_up24 = time_from_up24(vpind,:);time_from_wake = time_from_wake(vpind,:);time_from_wake24 = time_from_wake24(vpind,:);
time_from_wake(14,5)=17+55/60; %vp37 got up 7:15, t1 scan at 00:10, +24h messed up the table.
times24(times24==0.166666667000000)=24.166666667000000;
time_from_wake(time_from_wake==0)=NaN;
time_from_up(time_from_up==0)=NaN;
clear IDs

load(fullfile(anapath,'behavres','perf.mat'),'performance');
%right vps and right order
performance = performance(ismember(performance(:,9),vps),:);
performance = [performance(17:end,:);performance(1:16,:)];

%prepare plots
c1 = uint8([128 128 128]);
% c1 = uint8([50 65 75]);
% c2 = uint8([165 30 55]);
c2 = uint8([200 0 0]);
% c3 = uint8([180 77 80]);
% c3 = uint8([0 105 170]);
% c3 = uint8([65 90 140]);
c3 = uint8([0 0 128]);

group = ones(32,1); group(18:end)=2; %1 sleep 2 wake


%% measurement timepoints
cmap = colormap(jet(12));
figure
hold on
% histogram(times24(18:end,1),5,'EdgeColor','none','FaceAlpha',0.6,'FaceColor',c2)

set(gca, 'Box', 'off');
set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
set(get(gca,'YLabel'),'String','hours awake')
set(get(gca,'XLabel'),'String','time of day (hours)')
scatter(reshape(times24(18:end,[1 4]),[],1),reshape(time_from_wake(18:end,[1 4]),[],1),40,cmap(end,:),'filled')
scatter(reshape(times24(18:end,[2 5]),[],1),reshape(time_from_wake(18:end,[2 5]),[],1),40,cmap(end-1,:),'filled')
scatter(reshape(times24(18:end,[3 6]),[],1),reshape(time_from_wake(18:end,[3 6]),[],1),40,cmap(end-2,:),'filled')
scatter(reshape(times24(1:17,[1 4]),[],1),reshape(time_from_wake(1:17,[1 4]),[],1),40,cmap(1,:),'filled')
scatter(reshape(times24(1:17,[2 5]),[],1),reshape(time_from_wake(1:17,[2 5]),[],1),40,cmap(2,:),'filled')
scatter(reshape(times24(1:17,[3 6]),[],1),reshape(time_from_wake(1:17,[3 6]),[],1),40,cmap(3,:),'filled')
x = xlim;
y = ylim;
line([x(1) nanmean(reshape(times24(18:end,[1 4]),[],1))], [nanmean(reshape(time_from_wake(18:end,[1 4]),[],1)) nanmean(reshape(time_from_wake(18:end,[1 4]),[],1))],'Color',cmap(end,:),'LineWidth',3,'LineStyle',':')
line([nanmean(reshape(times24(18:end,[1 4]),[],1)) nanmean(reshape(times24(18:end,[1 4]),[],1))], [y(1) nanmean(reshape(time_from_wake(18:end,[1 4]),[],1))],'Color',cmap(end,:),'LineWidth',3,'LineStyle',':')
line([x(1) nanmean(reshape(times24(18:end,[2 5]),[],1))], [nanmean(reshape(time_from_wake(18:end,[2 5]),[],1)) nanmean(reshape(time_from_wake(18:end,[2 5]),[],1))],'Color',cmap(end-1,:),'LineWidth',3,'LineStyle',':')
line([nanmean(reshape(times24(18:end,[2 5]),[],1)) nanmean(reshape(times24(18:end,[2 5]),[],1))], [y(1) nanmean(reshape(time_from_wake(18:end,[2 5]),[],1))],'Color',cmap(end-1,:),'LineWidth',3,'LineStyle',':')
line([x(1) nanmean(reshape(times24(18:end,[3 6]),[],1))], [nanmean(reshape(time_from_wake(18:end,[3 6]),[],1)) nanmean(reshape(time_from_wake(18:end,[3 6]),[],1))],'Color',cmap(end-2,:),'LineWidth',3,'LineStyle',':')
line([nanmean(reshape(times24(18:end,[3 6]),[],1)) nanmean(reshape(times24(18:end,[3 6]),[],1))], [y(1) nanmean(reshape(time_from_wake(18:end,[3 6]),[],1))],'Color',cmap(end-2,:),'LineWidth',3,'LineStyle',':')
line([x(1) nanmean(reshape(times24(1:17,[1 4]),[],1))], [nanmean(reshape(time_from_wake(1:17,[1 4]),[],1)) nanmean(reshape(time_from_wake(1:17,[1 4]),[],1))],'Color',cmap(1,:),'LineWidth',3,'LineStyle',':')
line([nanmean(reshape(times24(1:17,[1 4]),[],1)) nanmean(reshape(times24(1:17,[1 4]),[],1))], [y(1) nanmean(reshape(time_from_wake(1:17,[1 4]),[],1))],'Color',cmap(1,:),'LineWidth',3,'LineStyle',':')
line([x(1) nanmean(reshape(times24(1:17,[2 5]),[],1))], [nanmean(reshape(time_from_wake(1:17,[2 5]),[],1)) nanmean(reshape(time_from_wake(1:17,[2 5]),[],1))],'Color',cmap(2,:),'LineWidth',3,'LineStyle',':')
line([nanmean(reshape(times24(1:17,[2 5]),[],1)) nanmean(reshape(times24(1:17,[2 5]),[],1))], [y(1) nanmean(reshape(time_from_wake(1:17,[2 5]),[],1))],'Color',cmap(2,:),'LineWidth',3,'LineStyle',':')
line([x(1) nanmean(reshape(times24(1:17,[3 6]),[],1))], [nanmean(reshape(time_from_wake(1:17,[3 6]),[],1)) nanmean(reshape(time_from_wake(1:17,[3 6]),[],1))],'Color',cmap(3,:),'LineWidth',3,'LineStyle',':')
line([nanmean(reshape(times24(1:17,[3 6]),[],1)) nanmean(reshape(times24(1:17,[3 6]),[],1))], [y(1) nanmean(reshape(time_from_wake(1:17,[3 6]),[],1))],'Color',cmap(3,:),'LineWidth',3,'LineStyle',':')
scatter(reshape(times24(18:end,[1 4]),[],1),reshape(time_from_wake(18:end,[1 4]),[],1),40,cmap(end,:),'filled')
scatter(reshape(times24(18:end,[2 5]),[],1),reshape(time_from_wake(18:end,[2 5]),[],1),40,cmap(end-1,:),'filled')
scatter(reshape(times24(18:end,[3 6]),[],1),reshape(time_from_wake(18:end,[3 6]),[],1),40,cmap(end-2,:),'filled')
scatter(reshape(times24(1:17,[1 4]),[],1),reshape(time_from_wake(1:17,[1 4]),[],1),40,cmap(1,:),'filled')
scatter(reshape(times24(1:17,[2 5]),[],1),reshape(time_from_wake(1:17,[2 5]),[],1),40,cmap(2,:),'filled')
scatter(reshape(times24(1:17,[3 6]),[],1),reshape(time_from_wake(1:17,[3 6]),[],1),40,cmap(3,:),'filled')
xlim([7 35])
ylim([0 19])

xticks =[nanmean(reshape(times24(18:end,[1 4]),[],1)) nanmean(reshape(times24(18:end,[2 5]),[],1)) nanmean(reshape(times24(1:17,[1 4]),[],1)) nanmean(reshape(times24(18:end,[3 6]),[],1)) nanmean(reshape(times24(1:17,[2 5]),[],1)) nanmean(reshape(times24(1:17,[3 6]),[],1))];
rxticks = round(xticks);

set(gca,'XTick',rxticks); 
set(gca,'XTicklabel',{'9','12','20','','23','9'})

yticks =[nanmean(reshape(time_from_wake(18:end,[1 4]),[],1)) nanmean(reshape(time_from_wake(18:end,[2 5]),[],1)) nanmean(reshape(time_from_wake(1:17,[1 4]),[],1)) nanmean(reshape(time_from_wake(18:end,[3 6]),[],1)) nanmean(reshape(time_from_wake(1:17,[2 5]),[],1)) nanmean(reshape(time_from_wake(1:17,[3 6]),[],1))];
ryticks = round(yticks);
set(gca,'YTick',ryticks([1:3 5 4])); 
set(gca,'YTicklabel',{'2','5','12','15','16'})

text(29,18,'wake t0','FontSize', 35,'FontWeight','bold','Color',cmap(end,:));
text(29,16.5,'wake t1','FontSize', 35,'FontWeight','bold','Color',cmap(end-1,:));
text(29,15,'wake t2','FontSize', 35,'FontWeight','bold','Color',cmap(end-2,:));
text(29,13.5,'sleep t0','FontSize', 35,'FontWeight','bold','Color',cmap(1,:));
text(29,12,'sleep t1','FontSize', 35,'FontWeight','bold','Color',cmap(2,:));
text(29,10.5,'sleep t2','FontSize', 35,'FontWeight','bold','Color',cmap(3,:));


%% 
mask = 3;
val = 3; %mean/median/std %%%%we look at 1 and 3!
norm = 1; %raw
vdata = params_ind2(:,:,mask,val,norm);
vtimes = time_from_wake;

% datamatrix = [vps' group vdata];
% save for spss
% save(fullfile(anapath,'spss_allm_MD_std.mat'), 'datamatrix');


%% morning higher than evening
%control morning, control evening, learn morning, learn evening
v = [vdata(1:17,6) vdata(1:17,5) vdata(1:17,3) vdata(1:17,2); vdata(18:end,4) vdata(18:end,6) vdata(18:end,1) vdata(18:end,3)];

%plot
figure
hold on
means = [nanmean(v)]; 
sems = [nanstd(v)/sqrt(length(v))];
model_series = [means(1:2); means(3:4)];
model_error = [sems(1:2);sems(3:4)];
h = bar(model_series);
set(h,'BarWidth',0.7);
h(1).FaceColor = c1-40;
h(2).FaceColor = c1+40;
% h(3).FaceColor = c3;
h(1).EdgeColor = 'none';
h(2).EdgeColor = 'none';
% h(3).EdgeColor = 'none';
set(gca, 'Box', 'off');
set(gca,'XTick',[1:2]);
set(gca,'XTicklabel',{'control', 'learn'})
set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
set(get(gca,'YLabel'),'String','MD raw')

numgroups = size(model_series, 1); 
numbars = size(model_series, 2); 
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    e(i) =errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','linewidth',3);
    %for significance bars
    bars(:,i) = x;   
end
% ylim([0.00076 0.00078])
% ylim([0.000057 0.000063])
%lines
xtext = h.XData;
for b = 1:numgroups
% % line(bars(b,:),[max(model_series(b,:))+max(model_error(b,:)) max(model_series(b,:))+max(model_error(b,:))]*1.0005,'LineWidth',3,'Color','k');
% % text(xtext(b),(max(model_series(b,:))+max(model_error(b,:)))*1.0007,'***','FontSize', 40,'FontWeight','bold','HorizontalAlignment','center');
% if b==1
% % text(xtext(b),(max(model_series(b,:))+max(model_error(b,:)))*1.0015,'*','FontSize', 40,'FontWeight','bold','HorizontalAlignment','center');
% % line(bars(b,:),[max(model_series(b,:))+max(model_error(b,:)) max(model_series(b,:))+max(model_error(b,:))]*1.001,'LineWidth',3,'Color','k');
% end
end
hl = legend('morning','evening');
hl.Location = 'northeast';
hl.Box = 'off';
hl.FontSize = 40;
hl.FontWeight = 'bold';

%test
fprintf('\nmorning vs evening\nCONTROL\n')
%morning higher than evening over all vps
evmor = [vdata(1:17,5)-vdata(1:17,6);vdata(18:end,6)-vdata(18:end,4)];
[~, pval, ci, stats]= ttest(evmor);
fprintf('ev lower morning %.6f pval %.3f\n',mean(evmor),pval)
fprintf('LEARN\n')
evmor = [vdata(1:17,2)-vdata(1:17,3);vdata(18:end,3)-vdata(18:end,1)];
[~, pval, ci, stats]= ttest(evmor);
fprintf('ev lower morning %.6f pval %.3f\n',mean(evmor),pval)

datamatrix = [vps' group v];
%save for spss
% save(fullfile(anapath,'spss_evening_morning_.mat'), 'datamatrix');

%% correlation MD and time spent awake - single regression for every vp and condition

%colormap = parula
cmap = colormap(parula(33));
clear r b m
figure
hold on
for i = 1:size(vdata,1)*2
    if i <= size(vdata,1)
        col = 1:3;
        vpind = i;
    else
        col = 4:6;
        vpind = i-size(vdata,1);
    end
    if i==15
        continue %no times since wake for that one
    end
[r(i),m(i),b(i)] = regression(vtimes(vpind,col),vdata(vpind,col),'one');
li = [min(vtimes(vpind,col)) max(vtimes(vpind,col))];
h2 = line(li, [0 m(i)*li(2)],'Color','k','LineWidth', 3, 'Color',cmap(vpind,:)); %we plot from zero, otherwise
scatter(vtimes(vpind,col),vdata(vpind,col)-vdata(vpind,find(vtimes(vpind,:)==li(1))),20,cmap(vpind,:),'filled')
end
rz = 0.5*[log(1+r)-log(1-r)];
[~,p,ci,stats]=ttest(rz);
rzmean = mean(rz);
rmean = tanh(rzmean);
% text(2,-1.25*0.00001,sprintf('r=%.3f***',rmean),'FontSize', 40,'FontWeight','bold'),
text(2,-0.75*0.00001,sprintf('r=%.3f, p=%.2f',rmean,p),'FontSize', 40,'FontWeight','bold'),
set(gca, 'Box', 'off');
set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
set(get(gca,'YLabel'),'String','delta MD')
set(get(gca,'XLabel'),'String','hours awake')
xlim([0 19])
set(gca,'XTick',[1:3:19]);

%% wake linear decrease over day
cmap = colormap(jet(12));
figure
hold on
v = [vdata(18:end,1:3); vdata(18:end,4:6)];
%plot
figure
hold on
bar(1,mean(v(:,1)),'FaceColor',cmap(end,:),'EdgeColor','none','BarWidth',0.7);
bar(2,mean(v(:,2)),'FaceColor',cmap(end-1,:),'EdgeColor','none','BarWidth',0.7);
bar(4,mean(v(:,3)),'FaceColor',cmap(end-2,:),'EdgeColor','none','BarWidth',0.7);
means = nanmean(v); 
sems = nanstd(v)/sqrt(length(v));
errorbar([1 2 4],means,sems, 'k', 'linestyle', 'none','linewidth',3);
set(gca, 'Box', 'off');
set(gca,'XTick',[1 2 4]);
set(gca,'XTicklabel',{'9:00','12:00','22:00'})
set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
set(get(gca,'YLabel'),'String','MD raw')
set(get(gca,'XLabel'),'String','time of day')
ylim([0.00077 0.000781])
% ylim([0.000057 0.000062])
line([1 4],[means(1)+sems(1) means(3)+sems(3)]*1.0022,'LineWidth',3,'Color','k');
text(2.5,(means(1)+sems(1)+means(3)+sems(3))/2*1.0024,'***','FontSize', 40,'FontWeight','bold','HorizontalAlignment','center');

%% sleep quadratic increase

cmap = colormap(jet(12));
v = [vdata(1:17,1:3); vdata(1:17,4:6)];
%plot
figure
hold on
bar(1,mean(v(:,1)),'FaceColor',cmap(1,:),'EdgeColor','none','BarWidth',0.7);
bar(2,mean(v(:,2)),'FaceColor',cmap(2,:),'EdgeColor','none','BarWidth',0.7);
bar(4,mean(v(:,3)),'FaceColor',cmap(3,:),'EdgeColor','none','BarWidth',0.7);
means = nanmean(v); 
sems = nanstd(v)/sqrt(length(v));
errorbar([1 2 4],means,sems, 'k', 'linestyle', 'none','linewidth',3);
set(gca, 'Box', 'off');
set(gca,'XTick',[1 2 4]);
set(gca,'XTicklabel',{'20:00','23:00','9:00'})
set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
set(get(gca,'YLabel'),'String','MD raw')
set(get(gca,'XLabel'),'String','time of day')
ylim([0.00076 0.00078])
line([1 2 4],[means(1)+sems(1) means(2)+sems(2) means(3)+sems(3)]*1.0022,'LineWidth',3,'Color','k');
text(2.5,(means(1)+sems(1)+means(3)+sems(3))/2*1.0005,'***','FontSize', 40,'FontWeight','bold','HorizontalAlignment','center');

% ylim([0.000059 0.000063])
% line([1 2 4],[means(1)+sems(1) means(2)+sems(2) means(3)+sems(3)]*1.0022,'LineWidth',3,'Color','k');
% text(2.5, [means(2)+sems(2)]*1.0005,'*','FontSize', 40,'FontWeight','bold','HorizontalAlignment','center');


%% wb: correlation voxels with higher dec before sleep higher increase over sleep

fpathname = 'D:\Data\MORV\nifti\diffusion_analysis\stats\new_sleep';
dpathname = 'D:\Data\MORV\nifti\diffusion_analysis\stats\new_sleep';
diffnames = {
     '\glm_C2vsC1_s_MD_0mm_6mm\t_glm_C2vsC1_s_MD_0mm_6mm.nii.gz',...
       '\glm_C2vsC1_w_MD_0mm_6mm\t_glm_C2vsC1_w_MD_0mm_6mm.nii.gz',...
     '\glm_L2vsL1_s_MD_0mm_6mm\t_glm_L2vsL1_s_MD_0mm_6mm.nii.gz',...   
    '\glm_L2vsL1_w_MD_0mm_6mm\t_glm_L2vsL1_w_MD_0mm_6mm.nii.gz'}; %.gz
funcnames = {'\glm_C1vsC0_s_MD_0mm_6mm\t_glm_C1vsC0_s_MD_0mm_6mm.nii.gz',...
 '\glm_C1vsC0_w_MD_0mm_6mm\t_glm_C1vsC0_w_MD_0mm_6mm.nii.gz',...    
'\glm_L1vsL0_s_MD_0mm_6mm\t_glm_L1vsL0_s_MD_0mm_6mm.nii.gz',...
    '\glm_L1vsL0_w_MD_0mm_6mm\t_glm_L1vsL0_w_MD_0mm_6mm.nii.gz'}; %.gz

funclabel = 'decrease t0-t1';
struclabel = 'increase t1-t2';

% figure
% hold on
for i = 1:4
 
if i == 1
    reso = 55; %49
elseif i == 2
    reso = 47;
elseif i == 3
    reso = 50;
elseif i == 4
    reso = 55;
end    

    
mnimask = niftiread('D:\Data\MORV\analysis\masks\mni_templates_fsl\MNI152_T1_2mm_brain_mask.nii.gz');
allmmask = niftiread('D:\Data\MORV\nifti\overlap_masks\MORV39\30overlap_allm.nii.gz'); %_ero1vox?

diff = niftiread([dpathname diffnames{i}]); %cols = L1/L2
diff = diff(:,:,:,1);
diff(mnimask == 0 | allmmask == 0)= NaN;
func = niftiread([fpathname funcnames{i}]); %rows = mem/act
func = func(:,:,:,2);
func(mnimask == 0 | allmmask == 0)= NaN;
func(func==0) = NaN;

[a,sortind] = sort(diff(:));
pairs = [diff(sortind), func(sortind)];
pairs2 = pairs(~isnan(pairs(:,1))& ~isnan(pairs(:,2)),:);
% [r, m ,b] = regression(pairs2(:,1)', pairs2(:,2)');
[r, m ,b] = regression(pairs2(:,2)', pairs2(:,1)');

% ax{i} = subplot(2,2,i);
figure
hold on

% hexscatter(pairs2(:,1), pairs2(:,2),'res',reso);
hexscatter(pairs2(:,2), pairs2(:,1),'res',reso);

% if ismember(i,[1,3])
%  ylabel(funclabel)
% end
% if ismember(i,[3,4])
% xlabel(struclabel);
% end

set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
% if ismember(i,[3 4])
%     ylim([-30 30])
% end
% ylim([-7 7]);
% xlim([-7 7]);
ylimit = ylim;
xlimit = xlim;

tx = xlimit(2) - 0.5*(abs(xlimit(1))+xlimit(2));
ty = ylimit(1) + 0.2*(abs(ylimit(1))+ylimit(2));
text(tx,ty,sprintf('r=%.3f',r),'FontSize',40, 'FontWeight','bold')

x = xlim;
y = x*m+b;
line(x',y','Color','k','LineWidth',5)
if i == 1
%     title('sleep control');
    ax{i}.XTick = [];
elseif i ==2
%     title('wake control');
ax{i}.XTick = [];
ax{i}.YTick = [];
elseif i ==3
    
%     title('sleep learn');
elseif i == 4
    ax{i}.YTick = [];
%     title('wake learn');
end
end
% end
% suptitle(sup);


%% correlation task-related activity and inc over sleep

% fpathname = 'D:\Data\MORV\func_spm_models\RFx69\';
% dpathname = 'D:\Data\MORV\nifti\diffusion_analysis\MORV39\stats';
% diffnames = {'\glm_L1vsL0__MD_6mm_8mm\t_glm_L1vsL0__MD_6mm_8mm.nii.gz',...
%     '\glm_L2vsL0__MD_6mm_8mm\t_glm_L2vsL0__MD_6mm_8mm.nii.gz'}; %.gz
% funcnames = {'spmT_0011.nii', 'spmT_0035.nii'};
% % sup = 'RFx69 recall';
fpathname = 'E:\Data\MORV\func_spm_models\RFx97';
dpathname = 'D:\Data\MORV\nifti\diffusion_analysis\stats\new_sleep';
diffnames = {'\glm_L2-L1vsC2-C1_s_MD_6mm_8mm\t_glm_L2-L1vsC2-C1_s_MD_6mm_8mm.nii.gz'}; %.gz
funcnames = {'\spmT_0001.nii'}; %.gz

funclabel = 't-value BOLD';
struclabel = 't-value MD (t1-t2)x(L-C)';

% figure
% hold on
for i = 1%1:4
 
res=50;

    
mnimask = niftiread('D:\Data\MORV\analysis\masks\mni_templates_fsl\MNI152_T1_2mm_brain_mask.nii.gz');
allmmask = niftiread('D:\Data\MORV\nifti\overlap_masks\MORV39\30overlap_allm.nii.gz'); %_ero1vox?

diff = niftiread([dpathname diffnames{i}]); %cols = L1/L2
diff = diff(:,:,:,1);
diff(mnimask == 0 | allmmask == 0)= NaN;
func = niftiread([fpathname funcnames{i}]); %rows = mem/act
func(mnimask == 0 | allmmask == 0)= NaN;
func(func==0) = NaN;

[a,sortind] = sort(diff(:));
pairs = [diff(sortind), func(sortind)];
pairs2 = pairs(~isnan(pairs(:,1))& ~isnan(pairs(:,2)),:);
% [rlearn, m ,b] = regression(pairs2(:,1)', pairs2(:,2)');
[rlearn, m ,b] = regression(pairs2(:,2)', pairs2(:,1)');

% ax{i} = subplot(2,2,i);
figure
hold on

% hexscatter(pairs2(:,1), pairs2(:,2),'res',reso);
hexscatter(pairs2(:,2), pairs2(:,1),'res',reso);


%  ylabel(funclabel)
  xlabel(funclabel)
% xlabel(struclabel);
ylabel(struclabel);

set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
% if ismember(i,[3 4])
%     ylim([-30 30])
% end
% ylim([-7 7]);
% xlim([-7 7]);
ylim([-6 6])
ylimit = ylim;
xlimit = xlim;

% tx = xlimit(1) + 0.1*(abs(xlimit(1))+xlimit(2));
% ty = ylimit(2) - 0.2*(abs(ylimit(1))+ylimit(2));
% text(tx,ty,sprintf('r=%.3f',rlearn),'FontSize',40, 'FontWeight','bold')
tx = xlimit(2) - 0.3*(abs(xlimit(1))+xlimit(2));
ty = ylimit(1) + 0.2*(abs(ylimit(1))+ylimit(2));
text(tx,ty,sprintf('r=%.3f***',rlearn),'FontSize',40, 'FontWeight','bold')

x = xlim;
y = x*m+b;
line(x',y','Color','k','LineWidth',5)
if i == 1
%     title('sleep control');
%     ax{i}.XTick = [];
elseif i ==2
%     title('wake control');
ax{i}.XTick = [];
ax{i}.YTick = [];
elseif i ==3
    
%     title('sleep learn');
elseif i == 4
    ax{i}.YTick = [];
%     title('wake learn');
end
end
% end
% suptitle(sup);


%% correlation performance und increase over sleep
mask = 9;
val = 1; %mean/median/std
norm = 1; %raw
vdata = params_ind2(:,:,mask,val,norm);
vvp = 1:17;%18:size(vdata,1);
sleepmin = (performance(vvp,5)-performance(vvp,4))*100;
%learn
v = vdata(vvp,3)-vdata(vvp,2); %control 6-5, learn 3-2
[rlearn, plearn] = corr(sleepmin,v,'rows','pairwise','type','Spearman');
fprintf('correlation perf and diff increase over sleep r=%.3f p=%.3f n=%d\n',rlearn,plearn,size(v,1));

%control
vcontrol = vdata(vvp,6)-vdata(vvp,5); %control 6-5, learn 3-2
[rcontrol, pcontrol] = corr(sleepmin,vcontrol,'rows','pairwise','type','Spearman');


   %steigers z test
    [rnoise, pnoise] =  corr(v,vcontrol,'rows','pairwise','type','Spearman');
    n = length(vvp);
    rm2 = (rlearn*rlearn + rcontrol*rcontrol)/2;
    f = (1 - rnoise)/(2*(1-rm2));
    h = (1-(f*rm2))/(1-rm2);
    z = (atanh(rlearn)-atanh(rcontrol))*sqrt(n-3)/(sqrt(2*(1-rnoise)*h));
    pz = normcdf(abs(z),'upper')*2;
 
fprintf('rlearn=%.3f, rcontrol=%.3f, rnoise=%.3f, n=%d, steigers z=%.3f p=%.3f\n',rlearn,rcontrol, rnoise, n, z, pz);



figure
hold on
scatter(v,sleepmin,70,[0 0 0],'filled');
[r1, m ,b] = regression(v, sleepmin,'one');
x = xlim;
x = [x(1)+0.000002 x(2)-0.000002];
y = x*m+b;
line(x',y','Color','k','LineWidth',5)
scatter(v,sleepmin,70,[0 0 0],'filled');
set(gca, 'Box', 'off');
% set(gca,'XTick',[1 2 4]);
% set(gca,'XTicklabel',{'20h','23h','9h'})
set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
set(get(gca,'XLabel'),'String','MD increase over sleep')
set(get(gca,'YLabel'),'String','performance gain %')
% ylim([0.00076 0.00078])
text(0.000015,5,sprintf('r(17)=%.3f*',rlearn),'FontSize', 40,'FontWeight','bold','HorizontalAlignment','center');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% identify precu and hippocampus voxels
% fpathname = 'D:\Data\MORV\func_spm_models\RFx69\';
% dpathname = 'D:\Data\MORV\nifti\diffusion_analysis\MORV39\stats';
% diffnames = {'\glm_L1vsL0__MD_6mm_8mm\t_glm_L1vsL0__MD_6mm_8mm.nii.gz',...
%     '\glm_L2vsL0__MD_6mm_8mm\t_glm_L2vsL0__MD_6mm_8mm.nii.gz'}; %.gz
% funcnames = {'spmT_0011.nii', 'spmT_0035.nii'};
% % sup = 'RFx69 recall';
fpathname = 'D:\Data\MORV\func_spm_models\RFx66\';
dpathname = 'D:\Data\MORV\nifti\diffusion_analysis\MORV39\stats';
diffnames = {'\glm_L1vsL0__MD_6mm_8mm\t_glm_L1vsL0__MD_6mm_8mm.nii.gz',...
    '\glm_L2vsL0__MD_6mm_8mm\t_glm_L2vsL0__MD_6mm_8mm.nii.gz'}; %.gz
funcnames = {'spmT_0009.nii', 'spmT_0033.nii'};
% sup = 'RFx66 learn';

funclabel = 't BOLD';
struclabel = 't MD raw';

%%%%%%%%%%%%%identify precu voxels
% pvoxels = 'E:\Data\MORV\1paper\images\diff\erodedmask\glm_MORV_L0_L1_dec_zstand_uncp_t_erodedmask_clindex.nii.gz';
% clindex = [16 14];
% clindex = [16 14 1 15 7 5 10 4 2]; %precu 3, moccgy 1, fusi 2, lingual 3
pvoxels = 'E:\Data\MORV\1paper\images\diff\erodedmask\glm_MORV_L0_L1_dec_raw_uncp_t_erodedmask_clindex.nii.gz';
% clindex = [16 14];
clindex = [16 14 13 3 2 1]; %precu 2, moccgy 1, fusi 3
pvox = niftiread(pvoxels);
pvox = ismember(pvox,clindex);

for i = 1:4
 
if i == 1
    reso = 50; %49
elseif i == 2
    reso = 57;
elseif i == 3
    reso = 59;
elseif i == 4
    reso = 70;
end    
% if i == 1
%     reso = 49; %50
% elseif i == 2
%     reso = 57;
% elseif i == 3
%     reso = 51;
% elseif i == 4
%     reso = 59;
% end
    
mnimask = niftiread('D:\Data\MORV\analysis\masks\mni_templates_fsl\MNI152_T1_2mm_brain_mask.nii.gz');
allmmask = niftiread('D:\Data\MORV\nifti\overlap_masks\MORV39\30overlap_allm_ero1vox.nii.gz');
epath = 'E:\Data\MORV\analysis\functional\multimages_t-values';
cd(epath)
diff = niftiread([dpathname diffnames{mod(i,2)*-1+2}]); %cols = L1/L2
diff = diff(:,:,:,2);
diff(mnimask == 0 | allmmask == 0)= NaN;
%%%%%%%%%%%%%%%%%
precudiff = diff;
precudiff(pvox==0)=NaN;

func = niftiread([fpathname funcnames{ceil(i/2)}]); %rows = mem/act
func(mnimask == 0 | allmmask == 0)= NaN;
func(func==0) = NaN;
%%%%%%%%%%
precufunc = func;
precufunc(pvox==0) = NaN;

[a,sortind] = sort(diff(:));
pairs = [diff(sortind), func(sortind)];
pairs2 = pairs(~isnan(pairs(:,1))& ~isnan(pairs(:,2)),:);
[rlearn, m ,b] = regression(pairs2(:,1)', pairs2(:,2)');

%%%%%%%%%%%%%%%%%%%
pairs3 = [precudiff(sortind), precufunc(sortind)];
pairs4 = pairs3(~isnan(pairs3(:,1))& ~isnan(pairs3(:,2)),:);

% ax{i} = subplot(2,2,i);
figure
hold on
hexscatter(pairs2(:,1), pairs2(:,2),'res',reso);
%%%%%%%%%%%%%%%%%%%%
scatter(pairs4(:,1), pairs4(:,2),15,'filled','r');

if ismember(i,[1 2])
    yy = [-10 10];
else
    yy = [-30 30];
end
line([0 0],yy,'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',3);
line([-5 5],[0 0],'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',3);


% if ismember(i,[1,3])
 ylabel(funclabel)
% end
% if ismember(i,[3,4])
xlabel(struclabel);
% end

set(gca,'FontSize', 40);
set(gca,'FontWeight','bold');
if ismember(i,[3 4])
    ylim([-30 30])
end
ylimit = ylim;
xlimit = xlim;
tx = xlimit(2) - 0.25*(abs(xlimit(1))+xlimit(2));
ty = ylimit(1) + 0.4*(abs(ylimit(1))+ylimit(2));
text(tx,ty,sprintf('r=%.3f',rlearn),'FontSize',40, 'FontWeight','bold')

x = xlim;
y = x*m+b;
line(x',y','Color','k','LineWidth',5)
% if i == 1
%     title('t1 & mem');
%     ax{i}.XTick = 
% elseif i ==2
%     title('t2 & mem');
%     ax{i}.XTick = [];
%     ax{i}.YTick = [];
% elseif i ==3
%     title('t1 & act');
% elseif i == 4
%     title('t2 & act');
%     ax{i}.YTick = [];
end
% end
% suptitle(sup);




