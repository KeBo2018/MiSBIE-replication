% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
clear
mainpath= 'F:\Mito_FmriPrep';
Dir=dir(mainpath)

% funcfile= strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-multisensory_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii')
% Covariatefile= 'ses-PicardMiSBIE\func\*task-multisensory_run-1_desc-confounds_timeseries.tsv'
% tsvData = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');


TR = 0.46; % Repetition time in seconds
total_duration = 652 * TR; % Total duration of the experiment in seconds, based on number of scans and TR



[MetaA MetaB MetaC]=xlsread('LOCAL_SOURCE\Dropbox (Dartmouth College)\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\MiSBIE MRI Meta Data 1-31-24.xlsm')
metaID= MetaB(2:end,1);
% Initialize an empty array to store the extracted numbers
IDnumbers = zeros(1, length(metaID));

k=1;
% Loop through each cell in the cell array
for i = 1:length(metaID)
    % Extract numbers using regular expression
    numStr = regexp(metaID{i}, '\d+', 'match');
    
    % Convert the string to a number and store it in the numbers array
    IDnumbers(i) = str2double(numStr{1});
end

%%%% Resting state analysis parameter setting%%%%
ASamplePeriod=0.46;
HighCutoff=0.1; LowCutoff=0.01;
AMaskFilename='E:\Mito_DICOM\Code\Resampled_CanlabMask.nii'; 
AResultFilename_Reho='F:\Mito_Rest\Result\Reho';
AResultFilename_ALFF='F:\Mito_Rest\Result\ALFF';
NVoxel=27;
TemporalMask=''; 
 Band=[0.01 0.08];
 IsNeedDetrend=1;
 TR=0.46;
ScrubbingMethod=1;
ScrubbingTiming=[];
 %%%%%%%%%%%%%%%
SampleHeaderInfo=niftiinfo('E:\Mito_DICOM\Code\Resampled_CanlabMask.nii');

Entropy_header_info=SampleHeaderInfo;
Entropy_header_info.Filename = 'entropy_matrix.nii';
Entropy_header_info.ImageSize = SampleHeaderInfo.ImageSize(1:3);
Entropy_header_info.PixelDimensions = SampleHeaderInfo.PixelDimensions(1:3); % Use only the first three dimensions
Entropy_header_info.Datatype = 'double'; 

Max_header_info=SampleHeaderInfo;
Max_header_info.Filename = 'Max_matrix.nii';
Max_header_info.ImageSize = SampleHeaderInfo.ImageSize(1:3);
Max_header_info.PixelDimensions = SampleHeaderInfo.PixelDimensions(1:3); % Use only the first three dimensions
Max_header_info.Datatype = 'double'; 
 for i=1:100

    subfolder=Dir(i*2+6).name;
    subname=str2num(subfolder(7:9));

    metaid=find(IDnumbers==subname);
funcfile= strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz')

    if ~exist(fullfile(mainpath, subfolder, funcfile), 'file')
        fprintf('File does not exist: %s\n', fullfile(mainpath, subfolder, funcfile));
        continue;  % Skip to the next iteration
    end

Covariatefile=strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-rest_run-1_desc-confounds_timeseries.tsv')
CovariateTable=readtable(fullfile(mainpath,subfolder,Covariatefile),'FileType', 'text', 'Delimiter', '\t')


colsToExtract = contains(CovariateTable.Properties.VariableNames, 'rot') | ...
                contains(CovariateTable.Properties.VariableNames, 'trans') | ...
                contains(CovariateTable.Properties.VariableNames, 'csf');

% Extract these columns
Covariate = table2array(CovariateTable(:, colsToExtract));
Covariate(:,5)=[];

MS_Data=fmri_data(fullfile(mainpath,subfolder,funcfile))


%%% Spike detect
[ds, expectedds, p, wh_outlier_uncorr, wh_outlier_corr] = mahal(MS_Data, 'noplot');
        SpikeIndex=find(wh_outlier_corr==1);
        A=[];
        if length(SpikeIndex)>0;
            A=zeros(size(MS_Data.dat,2),length(SpikeIndex));
            for Covi=1:length(SpikeIndex)


                A(SpikeIndex(Covi),Covi)=1;
            end
        end
%%
Covariate_Spike=[Covariate A];
MS_Data.covariates=Covariate_Spike(11:end,:);
MS_Data_DummyExclude=get_wh_image(MS_Data,11:size(MS_Data.dat,2));
% preprocessed_dat=canlab_connectivity_preproc(MS_Data_DummyExclude, TR,'additional_nuisance',Covariate_Spike(11:end,:),'bpf',[0.01 0.1],0.46,'detrend');
preprocessed_dat=canlab_connectivity_preproc(MS_Data_DummyExclude, TR,'additional_nuisance',Covariate_Spike(11:end,:),'detrend');

preprocessed_dat_smooth=preprocess(preprocessed_dat,'smooth',4)

% fname=strcat(savefile_Path,funcfile(23:59),'NuissanceRegress.nii');
% fname_Smooth=strcat(savefile_Path,funcfile(23:59),'NuissanceRegress_Smooth.nii');
% save(filename,variables) 

%%%%%

% ScrubbingMethod, 
% Header;
% CUTNUMBER;
AResultFilename_Reho=strcat('F:\Mito_Rest\Result_B_Filter_Detrend_Noprefilter\Reho\','Reho_',funcfile(23:59));
AResultFilename_ALFF=strcat('F:\Mito_Rest\Result_B_Filter_Detrend_Noprefilter\ALFF\','ALFF_',funcfile(23:59));

Header=preprocessed_dat.volInfo;
preprocessed_dat_smooth_mat=reconstruct_image(preprocessed_dat_smooth);
preprocessed_dat_mat=reconstruct_image(preprocessed_dat);

%%%% ALFF and ReHo
[ALFFBrain, fALFFBrain, Header] = y_alff_falff(preprocessed_dat_smooth_mat,ASamplePeriod, HighCutoff, LowCutoff, AMaskFilename, AResultFilename_ALFF, TemporalMask, ScrubbingMethod, Header)
[ReHoBrain, Header] = y_reho(preprocessed_dat_mat, NVoxel, AMaskFilename, AResultFilename_Reho, IsNeedDetrend, Band, TR, TemporalMask, ScrubbingMethod, ScrubbingTiming, Header)
%%%%% max value and entropy %%
% Assuming you have the 4-D matrix called 'data'
% data is of size 65x77x65x1390

% Reshape the 4-D matrix into a 2-D matrix where each column is a voxel time series
reshaped_data = reshape(preprocessed_dat_smooth_mat, [], size(preprocessed_dat_smooth_mat,4));
for vnum=1:size(reshaped_data,1)
    max_matrix(vnum)=max(reshaped_data(vnum,:));
    entropy_matrix(vnum)=compute_entropy(reshaped_data(vnum,:));
    
end
max_matrix_3D = reshape(max_matrix, size(preprocessed_dat_smooth_mat,1), size(preprocessed_dat_smooth_mat,2), size(preprocessed_dat_smooth_mat,3));
entropy_matrix_3D = reshape(entropy_matrix, size(preprocessed_dat_smooth_mat,1), size(preprocessed_dat_smooth_mat,2), size(preprocessed_dat_smooth_mat,3));

niftiwrite(entropy_matrix_3D, strcat('F:\Mito_Rest\Result_B_Filter_Detrend\Entropy\','Entropy_',funcfile(23:59),'.nii'),Entropy_header_info);
niftiwrite(max_matrix_3D, strcat('F:\Mito_Rest\Result_B_Filter_Detrend\Max\','Max_',funcfile(23:59),'.nii'),Max_header_info);

  if strcmp(MetaC(metaid+1,2),'Control')==1
        GroupID(k)=1;
        k=k+1;
  else 
        GroupID(k)=0;
        k=k+1;
    end 

    metaGrand(i,:)=MetaC(metaid+1,:)

 end