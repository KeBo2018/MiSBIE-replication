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

k=1;
 for i=1:100

    subfolder=Dir(i*2+6).name;
    subname=str2num(subfolder(7:9));

    metaid=find(IDnumbers==subname);
% funcfile= strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz')
% 
%     if ~exist(fullfile(mainpath, subfolder, funcfile), 'file')
%         fprintf('File does not exist: %s\n', fullfile(mainpath, subfolder, funcfile));
%         continue;  % Skip to the next iteration
%     end

Covariatefile=strcat('ses-PicardMiSBIE\func\',subfolder,'_ses-PicardMiSBIE_task-rest_run-1_desc-confounds_timeseries.tsv')

    if ~isfile(fullfile(mainpath,subfolder,Covariatefile))
        warning('Covariate file not found, skipping subject %s', subfolder);
        continue
    end

CovariateTable=readtable(fullfile(mainpath,subfolder,Covariatefile),'FileType', 'text', 'Delimiter', '\t')


colsToExtract = contains(CovariateTable.Properties.VariableNames, 'framewise_displacement')

% Extract these columns
FD = table2array(CovariateTable(:, colsToExtract));
% Covariate(:,5)=[];
meanFD(i,:)=nanmean(FD);

  if strcmp(MetaC(metaid+1,2),'Control')==1
        GroupID(k)=1;
        k=k+1;
  else 
        GroupID(k)=0;
        k=k+1;
    end 

    metaGrand(i,:)=MetaC(metaid+1,:)

 end