% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
clear
mainpath= 'F:\Mito_FmriPrep';
Dir=dir(mainpath)
% subfolder= 'sub-WITHHELD_SUBJECT'
% funcfile= strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-coldpressor_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii')
% Covariatefile= 'ses-PicardMiSBIE\func\*task-coldpressor_run-1_desc-confounds_timeseries.tsv'
% tsvData = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');
Cue = [0, 1; 26, 1; 67, 1; 93, 1; 134, 1; 160, 1; 201, 1; 227, 1];
TwoBack = [1, 25; 68, 25; 161, 25; 228, 25];
ZeroBack = [27, 25; 94, 25; 135, 25; 202, 25];
onsets1 = {Cue;TwoBack; ZeroBack};

Cue = [0, 1; 26, 1; 67, 1; 93, 1; 134, 1; 160, 1; 201, 1; 227, 1];
TwoBack = [1, 25; 68, 25; 135, 25; 161, 25];
ZeroBack = [27, 25; 94, 25; 202, 25; 228, 25];
onsets2 = {Cue;TwoBack; ZeroBack};

% Define TR and total length of the run
TR = 0.46;
totalRunLength = 583 * TR;  % Total duration of the fMRI run

% Call the onsets2fmridesign function with the total run length
[X1, delta, delta_hires, hrf] = onsets2fmridesign(onsets1, TR, totalRunLength);
[X2, delta, delta_hires, hrf] = onsets2fmridesign(onsets2, TR, totalRunLength);
X1(1:10,:)=[];X2(1:10,:)=[];
X=[X1;X2];
X(:,4)=[];
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


 for i=1:100

    subfolder=Dir(i*2+6).name;
funcfile1= strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-nback_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz')

    if ~exist(fullfile(mainpath, subfolder, funcfile1), 'file')
        fprintf('File does not exist: %s\n', fullfile(mainpath, subfolder, funcfile1));
        continue;  % Skip to the next iteration
    end

funcfile2= strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-nback_run-2_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz')

    if ~exist(fullfile(mainpath, subfolder, funcfile2), 'file')
        fprintf('File does not exist: %s\n', fullfile(mainpath, subfolder, funcfile2));
        continue;  % Skip to the next iteration
    end

        subname=str2num(subfolder(7:9));
%         subid=find(subID==subname);
    metaid=find(IDnumbers==subname);

Covariatefile1=strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-nback_run-1_desc-confounds_timeseries.tsv')
CovariateTable1=readtable(fullfile(mainpath,subfolder,Covariatefile1),'FileType', 'text', 'Delimiter', '\t')

Covariatefile2=strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-nback_run-2_desc-confounds_timeseries.tsv')
CovariateTable2=readtable(fullfile(mainpath,subfolder,Covariatefile2),'FileType', 'text', 'Delimiter', '\t')

colsToExtract1 = contains(CovariateTable1.Properties.VariableNames, 'rot') | ...
                contains(CovariateTable1.Properties.VariableNames, 'trans') | ...
                contains(CovariateTable1.Properties.VariableNames, 'csf');
colsToExtract2 = contains(CovariateTable2.Properties.VariableNames, 'rot') | ...
                contains(CovariateTable2.Properties.VariableNames, 'trans') | ...
                contains(CovariateTable2.Properties.VariableNames, 'csf');
% Extract these columns
Covariate1 = table2array(CovariateTable1(:, colsToExtract1));
Covariate1(:,5)=[];
Covariate2 = table2array(CovariateTable2(:, colsToExtract2));
Covariate2(:,5)=[];

    if size(Covariate2,1)<500
        fprintf('Run2 does not have intact data does not exist: %s\n', fullfile(mainpath, subfolder, funcfile2));
        continue;  % Skip to the next iteration
    end

N_Back_Data_Run1=fmri_data(fullfile(mainpath,subfolder,funcfile1))
N_Back_Data_Run1=preprocess(N_Back_Data_Run1,'smooth',8)
% Cold_Data=preprocess(Cold_Data,'smooth',8,'hpfilter',480,'outliers')

[ds, expectedds, p, wh_outlier_uncorr, wh_outlier_corr] = mahal(N_Back_Data_Run1, 'noplot');
        SpikeIndex=find(wh_outlier_corr==1);
        if length(SpikeIndex)>0;
            A1=zeros(size(N_Back_Data_Run1.dat,2),length(SpikeIndex));
            for Covi=1:length(SpikeIndex)


                A1(SpikeIndex(Covi),Covi)=1;
            end
        end
%%
Covariate_Spike1=[Covariate1 A1];

N_Back_Data_Run2=fmri_data(fullfile(mainpath,subfolder,funcfile2))
N_Back_Data_Run2=preprocess(N_Back_Data_Run2,'smooth',8)
% Cold_Data=preprocess(Cold_Data,'smooth',8,'hpfilter',480,'outliers')

[ds, expectedds, p, wh_outlier_uncorr, wh_outlier_corr] = mahal(N_Back_Data_Run2, 'noplot');
        SpikeIndex=find(wh_outlier_corr==1);
        if length(SpikeIndex)>0;
            A2=zeros(size(N_Back_Data_Run2.dat,2),length(SpikeIndex));
            for Covi=1:length(SpikeIndex)
                A2(SpikeIndex(Covi),Covi)=1;
            end
        end
%%
Covariate_Spike2=[Covariate2 A2];
Covariate_Spike2(1:10,:)=[];
Covariate_Spike1(1:10,:)=[];
Covariate_All=blkdiag(Covariate_Spike1,Covariate_Spike2)
N_Back_DummyExclude_Run1=get_wh_image(N_Back_Data_Run1,11:size(N_Back_Data_Run1.dat,2));
N_Back_DummyExclude_Run2=get_wh_image(N_Back_Data_Run2,11:size(N_Back_Data_Run2.dat,2));
% N_Back_Data_DummyExclude_Detrend_Run1=canlab_connectivity_preproc(N_Back_DummyExclude_Run1, 'linear_trend');
% N_Back_Data_DummyExclude_Detrend_Run2=canlab_connectivity_preproc(N_Back_DummyExclude_Run2, 'linear_trend');
N_Back_DummyExclude_Run1.images_per_session=573;
N_Back_DummyExclude_Run2.images_per_session=573;
% N_Back_DummyExclude_Run1=rescale(N_Back_DummyExclude_Run1,'percentchange');
% N_Back_DummyExclude_Run2=rescale(N_Back_DummyExclude_Run2,'percentchange');
N_Back_DummyExclude_Run1=rescale(N_Back_DummyExclude_Run1,'session_grand_mean_scaling_spm_style');
N_Back_DummyExclude_Run2=rescale(N_Back_DummyExclude_Run2,'session_grand_mean_scaling_spm_style');

N_Back_Data_DummyExclude_Detrend_Run1=canlab_connectivity_preproc(N_Back_DummyExclude_Run1, 'hpf',0.002,0.46,'linear_trend');
N_Back_Data_DummyExclude_Detrend_Run2=canlab_connectivity_preproc(N_Back_DummyExclude_Run2, 'hpf',0.002,0.46,'linear_trend');

N_Back_Data_DummyExclude_Detrend_Run1=canlab_connectivity_preproc(N_Back_DummyExclude_Run1, 'hpf',0.005,0.46,'linear_trend');
N_Back_Data_DummyExclude_Detrend_Run2=canlab_connectivity_preproc(N_Back_DummyExclude_Run2, 'hpf',0.005,0.46,'linear_trend');



N_Back_Data_DummyExclude_Detrend=image_math(N_Back_Data_DummyExclude_Detrend_Run1,N_Back_Data_DummyExclude_Detrend_Run2,'concatenate');
N_Back_Data_DummyExclude_Detrend.images_per_session=573;
% N_Back_Data_DummyExclude_Detrend_Rescaled=rescale(N_Back_Data_DummyExclude_Detrend,'session_grand_mean_scaling_spm_style');
runIndicator=zeros(size(N_Back_Data_DummyExclude_Detrend.dat,2),2);
runIndicator(1:size(N_Back_Data_DummyExclude_Detrend_Run1.dat,2),1)=1;
runIndicator(size(N_Back_Data_DummyExclude_Detrend_Run1.dat,2)+1:end,2)=1;

WholeVariable=[X Covariate_All runIndicator];
N_Back_Data_DummyExclude_Detrend.X=WholeVariable;
% out=regress(N_Back_Data_DummyExclude_Detrend);
out=regress(N_Back_Data_DummyExclude_Detrend,'robust');
RunIndic=get_wh_image(out.b,size(out.b.dat,2)); 

% Cold_Temp=get_wh_image(beta_dat,2)
% Control_Temp=get_wh_image(beta_dat,1)
% Recover_Temp=get_wh_image(beta_dat,3)
Twoback_Temp=get_wh_image(out.b,2);
Cue_Temp=get_wh_image(out.b,1);
Zeroback_Temp=get_wh_image(out.b,3);
Contrast_Temp=image_math(Twoback_Temp,Zeroback_Temp,'minus');
    if i==1
        Twoback=Twoback_Temp;
        Cue=Cue_Temp;
        Zeroback=Zeroback_Temp;
        Contrast=Contrast_Temp;
    else
        Twoback=image_math(Twoback,Twoback_Temp,'cat');
        Cue=image_math(Cue,Cue_Temp,'cat');
        Zeroback=image_math(Zeroback,Zeroback_Temp,'cat');
        Contrast=image_math(Contrast,Contrast_Temp,'cat');
    end

    if strcmp(MetaC(metaid+1,2),'Control')==1
        GroupID(k)=1;
        k=k+1;
    else
        GroupID(k)=0;
        k=k+1;
    end

    metaGrand(i,:)=MetaC(metaid+1,:)
    
 end


        Contrast=fmri_data(Contrast);
 T=ttest(fmri_data(Contrast));
 figure
 orthviews(T)

Contrast2=image_math(fmri_data(Twoback),fmri_data(Zeroback),'minus');
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