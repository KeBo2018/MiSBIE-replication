% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
% Load the Buckner resting-state network atlas
RestingstateNetwork = load_atlas('canlab2024');

% Directory where time series data are saved
DataDir = dir('F:\Mito_Rest\Restingstate_Timeseries\sub-*_RestingState.mat');

% Preallocate (optional: determine number of ROIs first)
num_subjects = length(DataDir);
example_data = load(fullfile(DataDir(1).folder, DataDir(1).name));
example_ts = apply_parcellation(example_data.preprocessed_dat_smooth, RestingstateNetwork);
num_ROIs = size(example_ts, 2);
Rvalue = zeros(num_ROIs, num_ROIs, num_subjects);

% Loop through each subject
for sub = 1:num_subjects

    % Load subject resting-state data
    filename = fullfile(DataDir(sub).folder, DataDir(sub).name);
    load(filename, 'preprocessed_dat_smooth');

    % Extract average time series per ROI/network
    Network_Timeseries = apply_parcellation(preprocessed_dat_smooth, RestingstateNetwork);  % [time x ROIs]

    % Compute ROI-to-ROI Pearson correlation matrix
    R = corr(Network_Timeseries);  % fast vectorized version

    % Store
    Rvalue(:, :, sub) = R;
  
end