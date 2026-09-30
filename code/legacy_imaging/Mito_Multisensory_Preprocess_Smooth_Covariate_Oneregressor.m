% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
clear
mainpath= 'F:\Mito_FmriPrep';
Dir=dir(mainpath)

rand
% funcfile= strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-multisensory_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii')
% Covariatefile= 'ses-PicardMiSBIE\func\*task-multisensory_run-1_desc-confounds_timeseries.tsv'
% tsvData = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');


TR = 0.46; % Repetition time in seconds
total_duration = 652 * TR; % Total duration of the experiment in seconds, based on number of scans and TR

% Onsets for rest and task conditions with durations. Each row in a cell array: [onset, duration]
% Here, durations are explicitly included as 30 seconds for each block
ons = {
    [0 30; 60 30; 120 30; 180 30; 240 30], % Rest periods: onset times with their durations
    [30 30; 90 30; 150 30; 210 30; 270 30]  % Task periods: onset times with their durations
};

[MetaA MetaB MetaC]=xlsread('LOCAL_SOURCE\Dropbox (Dartmouth College)\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\MiSBIE MRI Meta Data 1-31-24.xlsm')
metaID= MetaB(2:end,1);

% Assuming you're using SPM's canonical HRF without modifications
[X, delta, delta_hires, hrf] = onsets2fmridesign(ons, TR, total_duration, 'hrf');

X(1:10,:)=[];
savefile_Path='F:\Mito_Preprocessed_Onemodel\Multisensory\';

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

 for i=1:100

    subfolder=Dir(i*2+6).name;
    subname=str2num(subfolder(7:9));
%     subid=find(subID==subname);
    metaid=find(IDnumbers==subname);
    
funcfile= strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-multisensory_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz')

    if ~exist(fullfile(mainpath, subfolder, funcfile), 'file')
        fprintf('File does not exist: %s\n', fullfile(mainpath, subfolder, funcfile));
        continue;  % Skip to the next iteration
    end

Covariatefile=strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-multisensory_run-1_desc-confounds_timeseries.tsv')
CovariateTable=readtable(fullfile(mainpath,subfolder,Covariatefile),'FileType', 'text', 'Delimiter', '\t')


colsToExtract = contains(CovariateTable.Properties.VariableNames, 'rot') | ...
                contains(CovariateTable.Properties.VariableNames, 'trans') | ...
                contains(CovariateTable.Properties.VariableNames, 'csf');

% Extract these columns
Covariate = table2array(CovariateTable(:, colsToExtract));
Covariate(:,5)=[];

MS_Data=fmri_data(fullfile(mainpath,subfolder,funcfile))
% MS_Data=preprocess(MS_Data,'smooth',8,'hpfilter',480)
MS_Data=preprocess(MS_Data,'smooth',8)

%%% Spike detect
[ds, expectedds, p, wh_outlier_uncorr, wh_outlier_corr] = mahal(MS_Data, 'noplot');
        SpikeIndex=find(wh_outlier_corr==1);
        if length(SpikeIndex)>0;
            A=zeros(size(MS_Data.dat,2),length(SpikeIndex));
            for Covi=1:length(SpikeIndex)


                A(SpikeIndex(Covi),Covi)=1;
            end
        end
%%
Covariate_Spike=[Covariate A];
% MS_Data.covariates=Covariate_Spike(11:end,:);
MS_Data_DummyExclude=get_wh_image(MS_Data,11:size(MS_Data.dat,2));
MS_Data_DummyExclude.images_per_session=size(MS_Data_DummyExclude.dat,2);
MS_Data_DummyExclude=rescale(MS_Data_DummyExclude,'session_grand_mean_scaling_spm_style'); % Tor asked to do


% [preprocessed_dat, roi_val, maskdat, beta_dat, beta_roi_val]=canlab_connectivity_preproc(MS_Data_DummyExclude,'hpf', 0.004, TR,'regressors',X);
preprocessed_dat=canlab_connectivity_preproc(MS_Data_DummyExclude,'hpf', 0.002, TR,'detrend');


WholeVariable=[X(:,1:2) Covariate_Spike(11:end,:)]; %Put the design matrix and coviarte together


[X_Filtered, ~, ~, ~] = hpfilter(WholeVariable, TR, 500, size(WholeVariable,1));
MS_Data_DummyExclude.X=X_Filtered;
% out=regress(MS_Data_DummyExclude);
out=regress(MS_Data_DummyExclude,'robust');  %%%option: robust, or OLS (Default)

Task_Temp=get_wh_image(out.b,2);
Rest_Temp=get_wh_image(out.b,1);

Contrast_Temp=image_math(Task_Temp,Rest_Temp,'minus');

    if i==1

        Contrast=Contrast_Temp;
        Task=Task_Temp;
        Rest=Rest_Temp;
    else
        Contrast=image_math(Contrast,Contrast_Temp,'cat');
        Task=image_math(Task,Task_Temp,'cat');
        Rest=image_math(Rest,Rest_Temp,'cat');
    end

    metaGrand(i,:)=MetaC(metaid+1,:)
 end


 
T=ttest(fmri_data(Contrast));
figure
surface(T);

Contrast.dat=-1*Contrast.dat
T=ttest(Contrast);
figure
T_Threshold=threshold(T,0.05, 'fdr');
surface(T_Threshold);