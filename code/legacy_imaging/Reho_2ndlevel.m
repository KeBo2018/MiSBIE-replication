% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
% file names into fmri_data
clear

% pain_t = dir('D:\CANlab_Working\Data\PIP_GLM\*\spmT_0002.nii');    
Emotion_t = dir('F:\Mito_Rest\Result_B_Filter_Detrend_Noprefilter\Reho\Reho*.nii');    
% Emotion_t = dir('F:\Mito_Rest\Result_B_Filter_Detrend_Noprefilter\ALFF\ALFF*.nii');    

% Emotion_t = dir('F:\Mito_Rest\Result_Final\ALFF\fALFF*.nii');    
load('E:\Mito_DICOM\SecondLevelSave\MitoRest_FD.mat')

Emotion_fldr = {Emotion_t.folder}; 
fname = {Emotion_t.name};
Emotion_scan_files = strcat(Emotion_fldr,'\', fname)';
Reho = fmri_data(Emotion_scan_files);

T = readtable('F:\Mito_Rest\Meta_data\Misbie_Meta_Data.xlsx')
[a b c]=xlsread('F:\Mito_Rest\Meta_data\Misbie_Meta_Data.xlsx')

PatientIndex = zeros(91,1);  % Preallocate
for i = 1:91
    if strcmp(b{i+1,2},'''Control''') ==1
        PatientIndex(i) = 1;
    else
        PatientIndex(i) = 0;
    end
end
Control=find(PatientIndex==1);
Patient=find(PatientIndex==0); 

% Reho_T=ttest(Reho);
% figure
% montage(Reho_T)

Reho.metadata_table=T;
NMDS=table2array(Reho.metadata_table(:,"NMDAS"));
[group_metrics, individual_metrics, global_gm_wm_csf_values] = qc_metrics_second_level(Reho)
Reho.X(:,2)=[ global_gm_wm_csf_values(:,3)];
Reho.X(:,1)=NMDS;
Contrast_NanF=get_wh_image(Reho,find(isnan(NMDS)==0));
out=regress(Contrast_NanF)
GroupCon_T=get_wh_image(out.t,1);
% GroupCon_T_Thre=threshold(GroupCon_T,[-2 2], 'raw-outside');
    GroupCon_T_Thre=threshold(GroupCon_T, 0.001, 'unc','k', 20);
figure

montage(GroupCon_T_Thre,'compact2')


Reho_Selected=get_wh_image(Reho,Patient)
nonNanPatientIndex=find(isnan(Reho_Selected.X(:,1))==0)
Reho_Selected=get_wh_image(Reho_Selected,nonNanPatientIndex)

out=regress(Reho_Selected)
GroupCon_T=get_wh_image(out.t,1);
% GroupCon_T_Thre=threshold(GroupCon_T,[-2 2], 'raw-outside');
GroupCon_T_Thre=threshold(GroupCon_T, 0.05, 'unc','k', 20);
figure
montage(GroupCon_T_Thre,'compact2')

% figure
% montage(GroupCon_T_Thre,'full')


Reho.X(:,2)=[ global_gm_wm_csf_values(:,3)];
Reho.X(:,1)=PatientIndex;
% Reho.X(:,1)=NMDS;
Contrast_NanF=get_wh_image(Reho,find(isnan(PatientIndex)==0));
out=regress(Reho)
GroupCon_T=get_wh_image(out.t,1);
% GroupCon_T_Thre=threshold(GroupCon_T,[-2 2], 'raw-outside');
GroupCon_T_Thre=threshold(GroupCon_T, 0.001, 'unc','k', 20);
figure
montage(GroupCon_T_Thre,'compact2')



figure
montage(GroupCon_T_Thre,'compact2')

[bucknermaps, networknames] = load_image_set('bucknerlab');
bucknermaps.image_names = char(networknames{:});
bucknermaps_R=resample_space(bucknermaps, Reho)



%% Outlier %%
% Restingstate_Activation(:,4)=[];
NMDS(4)=[];
PatientIndex(4)=[];
Control=find(PatientIndex==1);
Patient=find(PatientIndex==0);
for i=1:7
Restingstateindex{i}=find(bucknermaps_R.dat(:,i)==1)
end
for i=1:7
temp=Reho.dat(find(bucknermaps_R.dat(:,i)==1),:);
temp(:,4)=[];
Restingstate_Activation(i,:)=mean(temp,1);

[h p ci stat]=ttest2(Restingstate_Activation(i,Control),Restingstate_Activation(i,Patient));
pval_RN(i)=p;
tval_RN(i)=stat.tstat;

[R p]=corr(NMDS(Patient),Restingstate_Activation(i,Patient)','rows','complete','type','spearman');
Rvalue_RN(i)=R;
PCorr_RN(i)=p;
end
figure
scatter(NMDS(Patient),Restingstate_Activation(4,Patient))
[R P]=corr(NMDS(Patient),Restingstate_Activation(2,Patient)','rows','complete')

[R P]=corr(NMDS,Restingstate_Activation(2,:)','rows','complete')
figure
scatter(NMDS(:),Restingstate_Activation(4,:))

for i=1:size(Reho.dat,1)
    [R P]=corr(Reho.dat(i,:)',meanFD,'type','spearman');
    Rvalue(i)=R;
    Pvalue(i)=P;
end
figure
% Histogram of R values
subplot(1,2,1);
histogram(Rvalue, 50, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
xline(0, 'k--', 'LineWidth', 1.5);
xlabel('Spearman R');
ylabel('Voxel Count');
title('Distribution of R values (ReHo vs meanFD)');
box off;

% Add summary stats as text
xL = xlim; yL = ylim;
text(xL(2)*0.6, yL(2)*0.85, ...
    sprintf('Mean = %.3f\nMedian = %.3f\nSD = %.3f', ...
    mean(Rvalue), median(Rvalue), std(Rvalue)), ...
    'FontSize', 9);

% Q-Q plot to check normality
subplot(1,2,2); 
qqplot(Rvalue);
title('Q-Q Plot of R values');
box off;

sgtitle('ReHo ~ meanFD Correlation Distribution');

[h, p, ci, stats] = ttest(Rvalue(:));
fprintf('One-sample t-test vs zero:\n');
fprintf('  t(%d) = %.3f, p = %.4f, 95%% CI [%.3f, %.3f]\n', ...
    stats.df, stats.tstat, p, ci(1), ci(2));
fprintf('  Mean R = %.4f\n', mean(Rvalue(:)));

% FDR correction across voxels
p_fdr = mafdr(Pvalue(:), 'BHFDR', true);
n_sig = sum(p_fdr < 0.05);
fprintf('Voxels significant after FDR correction: %d / %d (%.1f%%)\n', ...
    n_sig, numel(Pvalue), 100*n_sig/numel(Pvalue));
