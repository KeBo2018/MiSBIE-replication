% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
clear
mainpath= 'F:\Mito_FmriPrep';
Dir=dir(mainpath)
% subfolder= 'sub-WITHHELD_SUBJECT'
% funcfile= strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-coldpressor_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii')
% Covariatefile= 'ses-PicardMiSBIE\func\*task-coldpressor_run-1_desc-confounds_timeseries.tsv'
% tsvData = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');

condition1 = [0, 90];
condition2 = [120, 120];
condition3 = [270, 90];
Transition = [90 30; 240 30]; %% Model transition period

onsets = {condition1; condition2; condition3; Transition};

% Define TR and total length of the run
TR = 0.46;
totalRunLength = 782 * TR;  % Total duration of the fMRI run

% Call the onsets2fmridesign function with the total run length
[X, delta, delta_hires, hrf] = onsets2fmridesign(onsets, TR, totalRunLength);


[MetaA MetaB MetaC]=xlsread('LOCAL_SOURCE\Dropbox (Dartmouth College)\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\MiSBIE MRI Meta Data 1-31-24.xlsm')
metaID= MetaB(2:end,1);

% Initialize an empty array to store the extracted numbers
IDnumbers = zeros(1, length(metaID));

% Loop through each cell in the cell array
for i = 1:length(metaID)
    % Extract numbers using regular expression
    numStr = regexp(metaID{i}, '\d+', 'match');
    
    % Convert the string to a number and store it in the numbers array
    IDnumbers(i) = str2double(numStr{1});
end
k=1;

savefile_Path='F:\Mito_Preprocessed\ColdPressor\';
 for i=1:100

    subfolder=Dir(i*2+6).name;
funcfile= strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-coldpressor_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz')

    if ~exist(fullfile(mainpath, subfolder, funcfile), 'file')
        fprintf('File does not exist: %s\n', fullfile(mainpath, subfolder, funcfile));
        continue;  % Skip to the next iteration
    end
        subname=str2num(subfolder(7:9));
    metaid(i)=find(IDnumbers==subname);

Covariatefile=strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-coldpressor_run-1_desc-confounds_timeseries.tsv')
CovariateTable=readtable(fullfile(mainpath,subfolder,Covariatefile),'FileType', 'text', 'Delimiter', '\t')


colsToExtract = contains(CovariateTable.Properties.VariableNames, 'rot') | ...
                contains(CovariateTable.Properties.VariableNames, 'trans') | ...
                contains(CovariateTable.Properties.VariableNames, 'csf');

% Extract these columns
Covariate = table2array(CovariateTable(:, colsToExtract));
Covariate(:,5)=[];


Cold_Data=fmri_data(fullfile(mainpath,subfolder,funcfile))
Cold_Data=preprocess(Cold_Data,'smooth',8)
% Cold_Data=preprocess(Cold_Data,'smooth',8,'hpfilter',480,'outliers')

[ds, expectedds, p, wh_outlier_uncorr, wh_outlier_corr] = mahal(Cold_Data, 'noplot');
        SpikeIndex=find(wh_outlier_corr==1);
        if length(SpikeIndex)>0;
            A=zeros(size(Cold_Data.dat,2),length(SpikeIndex));
            for Covi=1:length(SpikeIndex)


                A(SpikeIndex(Covi),Covi)=1;
            end
        end
%%
Covariate_Spike=[Covariate A];
Covariate_Spike_Grand{i}=Covariate_Spike;

Cold_Data.covariates=Covariate_Spike(11:end,:);
Cold_Data_DummyExclude=get_wh_image(Cold_Data,11:size(Cold_Data.dat,2));

Cold_Data_DummyExclude.images_per_session=size(Cold_Data_DummyExclude.dat,2);
Cold_Data_DummyExclude=rescale(Cold_Data_DummyExclude,'session_grand_mean_scaling_spm_style');

% [preprocessed_dat, roi_val, maskdat, beta_dat, beta_roi_val]=canlab_connectivity_preproc(Cold_Data_DummyExclude, 'linear_trend','regressors',X(11:end,:));
Cold_Data_DummyExclude_Detrend=canlab_connectivity_preproc(Cold_Data_DummyExclude, 'linear_trend');
WholeVariable=[X Covariate_Spike];
Cold_Data_DummyExclude_Detrend.X=WholeVariable(11:end,:);
% out=regress(Cold_Data_DummyExclude_Detrend);
out=regress(Cold_Data_DummyExclude_Detrend,'robust');
% Cold_Temp=get_wh_image(beta_dat,2)
% Control_Temp=get_wh_image(beta_dat,1)
% Recover_Temp=get_wh_image(beta_dat,3)
Cold_Temp=get_wh_image(out.b,2);
Control_Temp=get_wh_image(out.b,1);
Recover_Temp=get_wh_image(out.b,3);
Transition_Temp=get_wh_image(out.b,4);

Contrast_Temp=image_math(Cold_Temp,Control_Temp,'minus');
    if i==1
        Recover=Recover_Temp;
        Control=Control_Temp;
        Cold=Cold_Temp;
        Contrast=Contrast_Temp;
        Transition=Transition_Temp;
    else
        Recover=image_math(Recover,Recover_Temp,'cat');
        Control=image_math(Control,Control_Temp,'cat');
        Cold=image_math(Cold,Cold_Temp,'cat');
        Contrast=image_math(Contrast,Contrast_Temp,'cat');
         Transition=image_math(Transition,Transition_Temp,'cat');
    end

    if strcmp(MetaC(metaid(i)+1,2),'Control')==1
        GroupID(k)=1;
        k=k+1;
    else
        GroupID(k)=0;
        k=k+1;
    end

    metaGrand(i,:)=MetaC(metaid(i)+1,:)
 end
        Recover=fmri_data(Recover);
        Control=fmri_data(Control);
        Cold=fmri_data(Cold);
        Transition=fmri_data(Transition);
        Contrast=fmri_data(Contrast);
 T=ttest(fmri_data(contrast));
 figure
 orthviews(T)

Contrast=image_math(fmri_data(Cold),fmri_data(Control),'minus');
Contrast=Transition;
T=ttest(fmri_data(Contrast))
figure
orthviews(T);
T_Threshold=threshold(T,0.05, 'fdr');
orthviews(T_Threshold);
figure
montage(T_Threshold,'full')
figure
surface(T_Threshold)

NPS_Average=apply_nps(Contrast)
mean(NPS_Average{1})