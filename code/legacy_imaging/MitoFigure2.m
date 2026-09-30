% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
clear 
% %%%%%%  Load 1st level contrast map and meta data %%%%%%%%%%
 load('E:\Mito_DICOM\SecondLevelSave\Nback_TwoRun_Right_order_500Filter_PercentageChange_SPM_Robust.mat')
%  load('E:\Mito_DICOM\SecondLevelSave\Cold_Onemodel_Robust_Transition_SPM.mat')
%  load('E:\Mito_DICOM\SecondLevelSave\Multisensory_OneModel_Filter_Robust_SPM.mat')
%  load('E:\Mito_DICOM\SecondLevelSave\Stress_onemodel_Robust_SPM.mat')

% load('E:\Mito_DICOM\SecondLevelSave\Nback_TwoRun_Right_order_Robust.mat')
% load('F:\SecondLevelSave\Multisensory_Onemodel_GroupID.mat')
% Contrast=fmri_data(image_math(Cold,Recover,'minus'))
% Contrast=fmri_data(Cold);
Contrast=fmri_data(Contrast);

Condition1=fmri_data(Twoback);Condition2=fmri_data(Zeroback);
% Condition1=fmri_data(Cold);Condition2=fmri_data(Control);
% Condition1=fmri_data(Task);Condition2=fmri_data(Rest);
% Condition1=fmri_data(Stress);Condition2=fmri_data(Control);


Contrast=fmri_data(Contrast);
T=ttest(Contrast)
Contrast.removed_images=[];
Contrast.removed_images=zeros(size(Contrast.dat,2),1);
% Condition1.removed_images=[];
% Condition2.removed_images=[];
T_Threshold=threshold(T,0.05, 'fdr','k',10);
figure
montage(T_Threshold,'full')

figure
surface(T_Threshold)
lightRestoreSingle

figure
[group_metrics individual_metrics values gwcsf gwcsfmean gwcsf_l2norm]=qc_metrics_second_level(Contrast)