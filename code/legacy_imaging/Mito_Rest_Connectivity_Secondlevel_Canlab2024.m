% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
% Load connectivity and metadata
clear
RestingstateNetwork = load_atlas('canlab2024');
load('E:\\Mito_DICOM\\SecondLevelSave\\Canlab2024_NetFC.mat');
load('E:\Mito_DICOM\SecondLevelSave\MitoRest_FD.mat')

T = readtable('F:\\Mito_Rest\\Meta_data\\Misbie_Meta_Data.xlsx');
[a, b, c] = xlsread('F:\\Mito_Rest\\Meta_data\\Misbie_Meta_Data.xlsx');

%% Define groups
PatientIndex = strcmp(b(2:end,2), '''Control''');
Control = find(PatientIndex == 1);
Patient = find(PatientIndex == 0);
GroupID      = double(~PatientIndex)';  % 0=Control, 1=Patient

%% Perform connectivity comparisons between groups

% Preallocate results
nRegions = size(Rvalue, 1);
bval = zeros(nRegions, nRegions);    % Group coefficient
pval = zeros(nRegions, nRegions);    % p-value of group effect
tval = zeros(nRegions, nRegions);    % t-statistic of group effect

% Loop through all ROI pairs
for i = 1:nRegions
    for j = 1:nRegions
        if i == j
            continue;  % Skip diagonal
        end

        % Extract connectivity values for all subjects for ROI pair (i,j)
        conn = squeeze(Rvalue(i, j, :));  % [91 x 1]

        % Design matrix: [group, meanFD]
        X = [GroupID', meanFD];

        % Use robust regression (adds intercept by default)
        [b, stats] = robustfit(X, conn);  % b(2)=group, b(3)=FD

        % Store group-level statistics
        bval(i, j) = b(2);                  % Coefficient for group
        tval(i, j) = stats.t(2);            % t-statistic for group
        pval(i, j) = stats.p(2);            % p-value for group
    end
end

% Exclude diagonal
pval(logical(eye(nRegions))) = NaN;
tval(logical(eye(nRegions))) = 0;

%% FDR Correction: Only upper triangle
upperIdx = find(triu(ones(nRegions), 1));

pval_vector = pval(upperIdx);
adj_pval_vector = mafdr(pval_vector, 'BHFDR', true);

adj_pval = nan(nRegions);
adj_pval(upperIdx) = adj_pval_vector;
adj_pval = adj_pval + adj_pval'; % mirror upper to lower triangle

significant_mask = adj_pval < 0.05;
if ~any(significant_mask(:))
    significant_mask = pval < 0.001;
end

%% Plot t-values with significant pairs circled
figure;
imagesc(tval*-1, [-4, 4]);
colorbar;
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm);
% 
% set(gca, 'XTick', 1:nRegions, 'XTickLabel', RestingstateNetwork.labels, ...
%          'YTick', 1:nRegions, 'YTickLabel', RestingstateNetwork.labels, ...
%          'XTickLabelRotation', 90, 'FontSize', 6);

set(gca, 'XTick', [], ...
         'YTick', [], ...
         'XTickLabelRotation', 90, 'FontSize', 6);

axis square;
title('Group Differences (Controls vs. Patients )');

hold on;
[ySig, xSig] = find(significant_mask);
scatter(xSig, ySig, 20, 'ks');
hold off;

%% Correlate connectivity with NMDAS within Patients
NMDAS = a(:,4);
% Preallocate
bval_corr = zeros(nRegions);    % beta for NMDAS
pval_corr = ones(nRegions);     % p-value for NMDAS
tval_corr = zeros(nRegions);    % t-statistic for NMDAS

% Extract patient-only data
NMDAS_pat = NMDAS(Patient);
FD_pat = meanFD(Patient);

% Loop over ROI pairs
for i = 1:nRegions
    for j = 1:nRegions
        if i == j
            continue;
        end

        conn_pat = squeeze(Rvalue(i, j, Patient));  % [nPatients x 1]

        % Design matrix: [NMDAS, FD]
        X = [NMDAS_pat, FD_pat];

        % Robust regression
        [b, stats] = robustfit(X, conn_pat);  % b(2)=NMDAS, b(3)=FD

        bval_corr(i, j) = b(2);        % beta for NMDAS
        tval_corr(i, j) = stats.t(2);  % t for NMDAS
        pval_corr(i, j) = stats.p(2);  % p for NMDAS
    end
end

% Clean up diagonal
bval_corr(logical(eye(nRegions))) = 0;
tval_corr(logical(eye(nRegions))) = 0;
pval_corr(logical(eye(nRegions))) = NaN;

% FDR correction
upperIdx_corr = find(triu(ones(nRegions), 1));
p_corr_vector = pval_corr(upperIdx_corr);
adj_p_corr_vector = mafdr(p_corr_vector, 'BHFDR', true);

adj_p_corr = nan(nRegions);
adj_p_corr(upperIdx_corr) = adj_p_corr_vector;
adj_p_corr = adj_p_corr + adj_p_corr';

significant_corr_mask = adj_p_corr < 0.05;
if ~any(significant_corr_mask(:))
    significant_corr_mask = pval_corr < 0.001;
end

% Plot t-value matrix with significant connections
figure;
imagesc(tval_corr, [-4, 4]);
colorbar;
cm = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
colormap(cm);

% set(gca, 'XTick', 1:nRegions, 'XTickLabel', RestingstateNetwork.labels, ...
%          'YTick', 1:nRegions, 'YTickLabel', RestingstateNetwork.labels, ...
%          'XTickLabelRotation', 90, 'FontSize', 6);
set(gca, 'XTick', [], ...
         'YTick', [], ...
         'XTickLabelRotation', 90, 'FontSize', 6);
axis square;
title('NMDAS–Connectivity (controlling for FD)');

hold on;
[yCorrSig, xCorrSig] = find(significant_corr_mask);
scatter(xCorrSig, yCorrSig, 20, 'ks');
hold off;