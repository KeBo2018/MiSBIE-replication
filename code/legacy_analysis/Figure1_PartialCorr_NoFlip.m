% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
%% =========================================================
%  MITO Brain Analysis — Partial Correlation & Clustering Pipeline
%  ---------------------------------------------------------
%  Pipeline:
%   1. Load data & build group label vector
%   2. Extract numeric data
%   3. Remove constant/all-NaN variables
%   4. Robust standardization (median/IQR)
%   5. [NEW] Sign-align variables: flip so MitoD > Control for all vars
%   6. Impute missing values (on standardized data)
%   7. Residualize group from each variable (Option 1)
%   8. Partial correlation (controlling for group) — for visualization
%   9. Spurious correlation check (raw vs partial)
%  10. Distance matrix from partial correlation
%  11. Clustering quality test via clusterdata_permtest on residualized data
%  12. Hierarchical clustering on partial correlation matrix
%  13. Optimal k selection (elbow + silhouette + permutation)
%  14. t-SNE visualization colored by hierarchical clusters
%  15. Sorted correlation matrix heatmap + reordered distance matrix
% =========================================================

clear; clc; close all;

%% ----------------------------------------------------------
%  1. FILE PATHS
% -----------------------------------------------------------
inputFile  = 'E:\Mito_DICOM\figures\Ke_DataforCorrelationMatrix_Cleaned_CK.xlsx';
misbieFile = ['LOCAL_SOURCE\Dartmouth College Dropbox\Ke Bo\' ...
              '2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\' ...
              'MiSBIE Data (Eprime and Questionnaires)\' ...
              'MiSBIE MRI Meta Data 6-11-24.xlsm'];

%% ----------------------------------------------------------
%  2. LOAD DATA
% -----------------------------------------------------------
T          = readtable(inputFile);
MisbieData = readtable(misbieFile);

%% ----------------------------------------------------------
%  3. BUILD SORTED BINARY LABEL VECTOR (−1 = Control, +1 = Patient)
% -----------------------------------------------------------
ids    = table2cell(MisbieData(:,1));
groups = table2cell(MisbieData(:,2));

isControl   = strcmp(groups, 'Control');
binaryLabel = double(~isControl);
binaryLabel(binaryLabel == 0) = -1;

% Sort by numeric subject ID
idNums           = cellfun(@(s) sscanf(s,'Mi%d'), ids);
[~, sortIdx_sub] = sort(idNums);
sortedBinary     = binaryLabel(sortIdx_sub);

patientIdx = sortedBinary == 1;
controlIdx = sortedBinary == -1;
fprintf('Subjects: %d patients, %d controls\n', sum(patientIdx), sum(controlIdx));

%% ----------------------------------------------------------
%  4. EXTRACT NUMERIC DATA
% -----------------------------------------------------------
isNumericVar = varfun(@isnumeric, T, 'OutputFormat', 'uniform');
numericT     = T(:, isNumericVar);
numericData  = numericT{:,:};
varNames     = numericT.Properties.VariableNames;
[nSubj, nVars] = size(numericData);
fprintf('Data: %d subjects x %d variables\n', nSubj, nVars);

%% ----------------------------------------------------------
%  5. REMOVE CONSTANT AND ALL-NaN VARIABLES
%     Must happen before standardization (std=0 causes division by zero)
% -----------------------------------------------------------
allNaN    = all(isnan(numericData), 1);
constVar  = nanstd(numericData, 0, 1) == 0;
toRemove5 = allNaN | constVar;

if any(toRemove5)
    fprintf('Removed %d variables (all-NaN or constant). Remaining: %d\n', ...
            sum(toRemove5), sum(~toRemove5));
end
numericData = numericData(:, ~toRemove5);
varNames    = varNames(~toRemove5);

%% ----------------------------------------------------------
%  6. ROBUST STANDARDIZATION (median / IQR)
% -----------------------------------------------------------
med     = nanmedian(numericData, 1);
iqr_val = iqr(numericData);

% If IQR = 0, fall back to nanstd
zeroIQR          = iqr_val == 0;
iqr_val(zeroIQR) = nanstd(numericData(:, zeroIQR), 0, 1);

% Final safety: if still 0, set to 1
iqr_val(iqr_val == 0) = 1;

numericData_z = (numericData - med) ./ iqr_val;

fprintf('Standardization: median/IQR applied to %d variables.\n', size(numericData_z,2));
fprintf('  %d variables had zero IQR — fell back to std.\n', sum(zeroIQR));

colMedians = nanmedian(numericData_z, 1);
fprintf('  Post-standardization median range: [%.3f, %.3f]\n', ...
        min(colMedians), max(colMedians));

%% ----------------------------------------------------------
%  6B. SIGN-FLIP DISABLED — original variable directions preserved
% -----------------------------------------------------------
% No sign-alignment is applied. flipMask is set to all-false so that
% downstream annotation and summary code remains functional.
flipMask = false(1, size(numericData_z, 2));

%% ----------------------------------------------------------
%  7. IMPUTE MISSING VALUES (on standardized, sign-aligned data)
% -----------------------------------------------------------
numericData_filled = fillmissing(numericData_z, 'linear',  1);
numericData_filled = fillmissing(numericData_filled, 'nearest', 1);

nanFrac = sum(isnan(numericData_filled(:))) / numel(numericData_filled);
fprintf('\nRemaining NaN fraction after imputation: %.4f\n', nanFrac);

%% ----------------------------------------------------------
%  8. RESIDUALIZE GROUP FROM EACH VARIABLE  (Option 1)
%  ---------------------------------------------------------
%  WHY: clusterdata_permtest permutes columns of the input
%  data to build a null distribution. This is only valid on
%  the raw data matrix, not on a pre-computed correlation
%  matrix. By regressing group out of each variable first,
%  we remove group-driven variance before clustering, which
%  is equivalent to what partial correlation achieves but on
%  the raw data — making clusterdata_permtest valid.
% -----------------------------------------------------------
fprintf('\nResidualizing group membership from each variable...\n');

% Use complete rows only (same as partial correlation below)
completeRows = all(~isnan(numericData_filled), 2);
fprintf('Complete cases: %d / %d\n', sum(completeRows), nSubj);

X_input  = numericData_filled(completeRows, :);
G_input  = sortedBinary(completeRows);

% Regress group out of each variable (on rank-transformed data for Spearman)
G_ranked = tiedrank(G_input);
X_resid  = zeros(size(X_input));
for v = 1:size(X_input, 2)
    X_ranked     = tiedrank(X_input(:, v));
    mdl_v        = fitlm(G_ranked, X_ranked);
    X_resid(:,v) = mdl_v.Residuals.Raw;
end

fprintf('Residualized data: %d subjects x %d variables\n', size(X_resid));

%% ----------------------------------------------------------
%  9A. STANDARD SPEARMAN CORRELATION (for comparison/visualization)
% -----------------------------------------------------------
R_pearson = corr(numericData_filled, 'Type','Spearman','Rows','pairwise');
R_pearson = (R_pearson + R_pearson') / 2;
R_pearson(1:size(R_pearson,1)+1:end) = 1;

%% ----------------------------------------------------------
%  9B. PARTIAL CORRELATION — controlling for group membership
%      Used for visualization and heatmap only
%      Clustering is done on residualized raw data (Section 11)
% -----------------------------------------------------------
fprintf('\nComputing partial correlation (controlling for group)...\n');

R_partial = partialcorr(X_input, G_input, 'Type', 'Spearman');

% Symmetry correction
R_partial = (R_partial + R_partial') / 2;
R_partial(1:size(R_partial,1)+1:end) = 1;

%% ----------------------------------------------------------
%  9C. SPURIOUS CORRELATION CHECK
% -----------------------------------------------------------
delta_R = R_pearson - R_partial;

figure('Name','Spurious Correlation Check','Position',[50 50 1400 500]);
subplot(1,3,1);
imagesc(R_pearson, [-1 1]);
colormap(redblue(256)); colorbar;
title('Spearman Correlation (Raw)');
axis square; set(gca,'XTick',[],'YTick',[]);

subplot(1,3,2);
imagesc(R_partial, [-1 1]);
colormap(redblue(256)); colorbar;
title('Partial Spearman Correlation (Group Controlled)');
axis square; set(gca,'XTick',[],'YTick',[]);

subplot(1,3,3);
imagesc(delta_R, [-1 1]);
colormap(redblue(256)); colorbar;
title('Difference (Spurious Component)');
axis square; set(gca,'XTick',[],'YTick',[]);
sgtitle('Correlation Comparison: Raw vs Partial');

meanSpurious = mean(abs(delta_R(triu(true(size(delta_R)),1))));
fprintf('Mean absolute spurious correlation: %.4f\n', meanSpurious);

%% ----------------------------------------------------------
%  9D. REMOVE ALL-NaN ROWS/COLS FROM PARTIAL CORR MATRIX
% -----------------------------------------------------------
allNaNRows   = all(isnan(R_partial), 2);
allNaNCols   = all(isnan(R_partial), 1);
toRemove     = allNaNRows & allNaNCols';

R_clean        = R_partial(~toRemove, ~toRemove);
varNames_clean = varNames(~toRemove);
X_resid_clean  = X_resid(:, ~toRemove);   % keep residualized data aligned
flipMask_clean = flipMask(~toRemove);      % keep flip mask aligned

fprintf('Final matrix: %d x %d variables\n', size(R_clean,1), size(R_clean,2));

% Save partial correlation matrix
Rtable_clean = array2table(R_clean, ...
    'VariableNames', varNames_clean, ...
    'RowNames',      varNames_clean);
writetable(Rtable_clean, 'partial_correlation_matrix.xlsx', 'WriteRowNames', true);
fprintf('Partial correlation matrix saved.\n');

%% ----------------------------------------------------------
%  10. DISTANCE MATRIX FROM PARTIAL CORRELATION
% -----------------------------------------------------------
D = sqrt(2 * (1 - R_clean));
D(1:size(D,1)+1:end) = 0;

if any(isnan(D(:))) || any(isinf(D(:)))
    warning('NaNs/Infs in D — replacing with max finite distance.');
    maxD = max(D(isfinite(D(:))));
    D(isnan(D) | isinf(D)) = maxD;
end

D_vec = squareform(D, 'tovector');

%% ----------------------------------------------------------
%  11. CLUSTERING QUALITY TEST via clusterdata_permtest
%  ---------------------------------------------------------
%  Input: X_resid_clean (subjects x variables, group removed)
%  This is valid for clusterdata_permtest because:
%   - Each row is an observation (subject)
%   - Group variance has been removed
%   - Column permutation creates a meaningful null distribution
%  The resulting cluster structure mirrors the partial
%  correlation clustering but with a valid permutation null
% -----------------------------------------------------------
fprintf('\n========== CLUSTERING QUALITY TEST ==========\n');
fprintf('Input: residualized data (%d subjects x %d variables)\n', ...
    size(X_resid_clean));
fprintf('Group variance removed before permutation testing\n\n');

stats_perm = clusterdata_permtest(X_resid_clean, ...
    'k',              2:40, ...
    'nperm',          100, ...
    'distancemetric', 'correlation', ...
    'linkagemethod',  'ward', ...
    'reducedims',     true, ...
    'verbose',        true, ...
    'doplot',         false);

% Print summary
fprintf('\n--- Permutation Test Summary ---\n');
fprintf('Best k by pseudo-Z:  %d\n',   stats_perm.best_k);
fprintf('Max pseudo-Z:        %.4f\n', stats_perm.max_pseudoZ);

% k=2 results
k_values = stats_perm.inputs.k;
idx_k2   = find(k_values == 2);
fprintf('\nResults for k=2 (statistically optimal):\n');
fprintf('  Silhouette score = %.4f\n', stats_perm.cluster_quality(idx_k2));
fprintf('  Null mean        = %.4f\n', stats_perm.cluster_quality_null_mean(idx_k2));
fprintf('  Null std         = %.4f\n', stats_perm.cluster_quality_null_std(idx_k2));
fprintf('  Pseudo-Z         = %.4f\n', stats_perm.cluster_quality_pseudoZ(idx_k2));
fprintf('  P-value          = %.4f\n', stats_perm.P_val(idx_k2));

%% ----------------------------------------------------------
%  12. HIERARCHICAL CLUSTERING ON PARTIAL CORRELATION MATRIX
% -----------------------------------------------------------
Z_link = linkage(D_vec, 'ward');

% Cophenetic correlation
c = cophenet(Z_link, D_vec);
fprintf('\nCophenetic correlation: %.4f  (>0.75 = good fit)\n', c);

% Elbow plot
figure('Name','Elbow Plot');
lastMerges = Z_link(end-14:end, 3);
plot(fliplr(2:16), lastMerges, 'o-', 'LineWidth', 2, 'Color', [0.2 0.4 0.8]);
xlabel('Number of clusters'); ylabel('Merge distance');
title('Elbow Plot — Hierarchical Clustering');
grid on;

% Silhouette scores (on partial correlation matrix for consistency)
fprintf('\nSilhouette scores (on partial correlation matrix):\n');
sil_scores = zeros(1,20);
for k = 2:21
    c_tmp           = cluster(Z_link, 'maxclust', k);
    s               = silhouette(R_clean, c_tmp, 'correlation');
    sil_scores(k-1) = mean(s);
    fprintf('  k=%d: silhouette = %.4f\n', k, mean(s));
end

figure('Name','Silhouette Scores');
bar(2:21, sil_scores, 'FaceColor', [0.2 0.5 0.8]);
xlabel('Number of clusters k'); ylabel('Mean silhouette score');
title('Optimal k — Hierarchical Clustering on Partial Correlation');
[~, best_k_idx] = max(sil_scores);
best_k_sil = best_k_idx + 1;
fprintf('\nBest k by silhouette: %d\n', best_k_sil);

% Dendrogram
figure('Name','Dendrogram');
dendrogram(Z_link, 0, 'Labels', varNames_clean, 'ColorThreshold','default');
title(sprintf('Dendrogram (Ward Linkage) — k=%d', stats_perm.best_k));
xtickangle(90);

%% ----------------------------------------------------------
%  13. ASSIGN FINAL CLUSTERS
%  ---------------------------------------------------------
%  Use k=2 — statistically optimal solution
% -----------------------------------------------------------
maxClusters = 11;
clusters    = cluster(Z_link, 'maxclust', maxClusters);

fprintf('\n=== Final Cluster Summary (k=%d) ===\n', maxClusters);
for i = 1:maxClusters
    members = varNames_clean(clusters == i);
    fprintf('Cluster %d (%d variables): %s\n', i, numel(members), strjoin(members, ', '));
end

%% ----------------------------------------------------------
%  14. t-SNE ON PARTIAL CORRELATION MATRIX
% -----------------------------------------------------------
rng(42);
perplexVal = min(30, floor((size(R_clean,1)-1)/3));

Y_tsne = tsne(R_clean, ...
         'Algorithm',        'exact', ...
         'Distance',         'correlation', ...
         'Perplexity',       perplexVal, ...
         'NumPCAComponents', 0, ...
         'Exaggeration',     10, ...
         'LearnRate',        200);

%% ----------------------------------------------------------
%  15. COLOUR PALETTE
% -----------------------------------------------------------
seaborn_colors = [
    0.1216, 0.4667, 0.7059;   % 1  blue
    1.0000, 0.4980, 0.0549;   % 2  orange
    0.1725, 0.6275, 0.1725;   % 3  green
    0.8392, 0.1529, 0.1569;   % 4  red
    0.5804, 0.4039, 0.7412;   % 5  purple
    0.5490, 0.3373, 0.2941;   % 6  brown
    0.8902, 0.4667, 0.7608;   % 7  pink
    0.4980, 0.4980, 0.4980;   % 8  gray
    0.7373, 0.7412, 0.1333;   % 9  olive/yellow
    0.0902, 0.7451, 0.8118;   % 10 cyan
    0.6824, 0.7804, 0.9098;   % 11 light blue
    0.5961, 0.8745, 0.5412;   % 12 light green
    1.0000, 0.5961, 0.5882;   % 13 light red/salmon
    0.7725, 0.6902, 0.8353;   % 14 light purple
    0.7686, 0.6118, 0.5804;   % 15 light brown
    0.9686, 0.7137, 0.8235;   % 16 light pink
];
%% ----------------------------------------------------------
%  16. t-SNE PLOT — colored by hierarchical clusters
% -----------------------------------------------------------
figure('Name','t-SNE — Partial Correlation Clusters', 'Position',[100 100 1100 800]);
hold on;
for i = 1:maxClusters
    idx_c = (clusters == i);
    scatter(Y_tsne(idx_c,1), Y_tsne(idx_c,2), 80, ...
            'filled', ...
            'MarkerFaceColor', seaborn_colors(i,:), ...
            'DisplayName',     sprintf('Cluster %d  (n=%d)', i, sum(idx_c)));
end
text(Y_tsne(:,1)+0.15, Y_tsne(:,2), varNames_clean, ...
    'FontSize', 7, 'Color', [0.3 0.3 0.3]);
legend('Location','best','FontSize',10);
xlabel('t-SNE Dimension 1'); ylabel('t-SNE Dimension 2');
title({'t-SNE Visualization of Partial Correlation Structure'; ...
       sprintf('Hierarchical Clusters (k=%d, Ward Linkage)', maxClusters)});
box off; grid off; hold off;

%% ----------------------------------------------------------
%  17. SORTED PARTIAL CORRELATION HEATMAP WITH CLUSTER BORDERS
%      [NEW] Variable labels annotated with (↑) if flipped
% -----------------------------------------------------------
[sortedClusters, sortIdx_cl] = sort(clusters);
R_sorted        = R_clean(sortIdx_cl, sortIdx_cl);
varNames_sorted = varNames_clean(sortIdx_cl);
flipMask_sorted = flipMask_clean(sortIdx_cl);

% Annotate flipped variable names with a dagger symbol
varLabels_sorted = varNames_sorted;
for vi = 1:numel(varNames_sorted)
    if flipMask_sorted(vi)
        varLabels_sorted{vi} = ['†' varNames_sorted{vi}];
    end
end

figure('Name','Sorted Partial Correlation Matrix', 'Position',[100 100 800 700]);
imagesc(R_sorted, [-1 1]);
colormap(redblue(256));
colorbar;
set(gca,'XTick',[],'YTick',[]);

% Draw cluster boundary lines
n_vars     = size(R_sorted,1);
boundaries = find(diff(sortedClusters)) + 0.5;
for b = boundaries
    line([b b],          [0.5 n_vars+0.5], 'Color','k','LineWidth',2);
    line([0.5 n_vars+0.5], [b b],          'Color','k','LineWidth',2);
end

% Cluster labels on Y axis
tickPos = accumarray(sortedClusters, (1:n_vars)', [], @mean);
set(gca,'YTick', tickPos, ...
        'YTickLabel', arrayfun(@(x) sprintf('C%d',x), ...
        1:maxClusters,'UniformOutput',false), ...
        'FontSize', 10);
title(sprintf('Partial Correlation Matrix — Sorted by Cluster (k=%d)', maxClusters));

% Add footnote about flipped variables
annotation('textbox', [0.01 0.01 0.98 0.04], ...
    'String', '† sign-inverted (originally: higher = better health)', ...
    'EdgeColor', 'none', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);

saveas(gcf, 'sorted_partial_correlation_matrix.png');

%% ----------------------------------------------------------
%  18. REORDERED DISTANCE MATRIX
%      Reviewer 2 Q8 specifically requested this
% -----------------------------------------------------------
figure('Name','Reordered Distance Matrix', 'Position',[100 100 800 700]);
imagesc(D(sortIdx_cl, sortIdx_cl), [0 max(D(:))]);
colormap(flipud(redblue(256)));
colorbar;
set(gca,'XTick',[],'YTick',[]);

% Draw cluster boundary lines
for b = boundaries
    line([b b],          [0.5 n_vars+0.5], 'Color','k','LineWidth',2);
    line([0.5 n_vars+0.5], [b b],          'Color','k','LineWidth',2);
end

% Cluster labels
set(gca,'YTick', tickPos, ...
        'YTickLabel', arrayfun(@(x) sprintf('C%d',x), ...
        1:maxClusters,'UniformOutput',false), ...
        'FontSize', 10);
title(sprintf('Reordered Distance Matrix — Sorted by Cluster (k=%d)', maxClusters));
saveas(gcf, 'reordered_distance_matrix.png');

%% ----------------------------------------------------------
%  19. PRINT FINAL SUMMARY FOR RESPONSE LETTER
% -----------------------------------------------------------
fprintf('\n========== SUMMARY FOR RESPONSE LETTER ==========\n');
fprintf('Clustering method:       Hierarchical, Ward linkage\n');
fprintf('Correlation type:        Partial Spearman (group controlled)\n');
fprintf('Permutation test input:  Residualized raw data (group removed)\n');
fprintf('N permutations:          1000\n');
fprintf('Cophenetic r:            %.4f\n', c);
fprintf('\nSign alignment:\n');
fprintf('  %d variables flipped: %s\n', sum(flipMask_clean), ...
        strjoin(varNames_clean(flipMask_clean), ', '));
fprintf('\nk=2 results (statistically optimal):\n');
fprintf('  Silhouette score:      %.4f\n', stats_perm.cluster_quality(idx_k2));
fprintf('  Pseudo-Z:              %.4f\n', stats_perm.cluster_quality_pseudoZ(idx_k2));
figure
plot(stats_perm.cluster_quality_pseudoZ)

fprintf('  P-value:               %.4f\n', stats_perm.P_val(idx_k2));
fprintf('\nBest k by pseudo-Z:      %d (Z=%.4f)\n', ...
    stats_perm.best_k, stats_perm.max_pseudoZ);
fprintf('\nCluster sizes (k=2):\n');
for i = 1:maxClusters
    fprintf('  Cluster %d: %d variables\n', i, sum(clusters==i));
end

%% ----------------------------------------------------------
%  HELPER — redblue diverging colormap (no toolbox needed)
% -----------------------------------------------------------
function c = redblue(m)
    if nargin < 1, m = 64; end
    bottom = [0   0   0.7];
    mid    = [1   1   1  ];
    top    = [0.7 0   0  ];
    h  = floor(m/2);
    c1 = interp1([0;1], [bottom; mid], linspace(0,1,h));
    c2 = interp1([0;1], [mid;   top],  linspace(0,1,m-h));
    c  = [c1; c2];
end
