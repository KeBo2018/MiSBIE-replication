%% =========================================================
%  MITO Brain Analysis — Partial Correlation & Clustering Pipeline
%  ---------------------------------------------------------
%  Pipeline:
%   1. Load data & build group label vector
%   2. Extract numeric data
%   3. Remove constant/all-NaN variables
%   4. Robust standardization (median/IQR)
%   5. Sign-align variables: flip so MitoD > Control for all vars
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

clearvars; clc; close all;
% User-confirmed Figure 1 source, adapted only for portable inputs/outputs,
% runtime fixes and explicitly reported optional-toolbox diagnostics.
repoRoot=misbieFindPackageRoot(mfilename('fullpath'));
addpath(fullfile(repoRoot,'code'),'-begin');
out=fullfile(repoRoot,'outputs','Figure1');if ~isfolder(out),mkdir(out);end
assert(exist('partialcorr','file')==2,'Statistics and Machine Learning Toolbox is required.');
runPermutationTest=exist('clusterdata_permtest','file')==2;
runTargetCorrelations=true; % Extra original NMDAS diagnostic; not Tables S2/S3.


%% ----------------------------------------------------------
%  1. FILE PATHS
% -----------------------------------------------------------
inputFile=fullfile(repoRoot,'data','phenotype','phenotype_scores.csv');
groupFile=fullfile(repoRoot,'data','participants','cohort_groups.csv');
T=readtable(inputFile,'VariableNamingRule','preserve');
Groups=readtable(groupFile,'VariableNamingRule','preserve');
% IDs are pseudonymized consistently across both CSVs. Preserve original
% participant row order: the source interpolates missing values along rows.
[found,groupOrder]=ismember(string(T.ReleaseID),string(Groups.ReleaseID));
assert(all(found) && numel(unique(string(T.ReleaseID)))==height(T), ...
    'Participant IDs in the phenotype and group files do not align.');
sortedBinary=1-2*Groups.GroupID(groupOrder); % -1 Control, +1 MitoD
patientIdx=sortedBinary==1;controlIdx=sortedBinary==-1;
assert(height(T)==110 && sum(patientIdx)==40 && sum(controlIdx)==70, ...
    'Expected the original 110-participant phenotyping cohort.');
fprintf('Subjects: %d patients, %d controls\n',sum(patientIdx),sum(controlIdx));
writetable(table(string(T.ReleaseID),sortedBinary,'VariableNames', ...
    {'ReleaseID','GroupCode'}),fullfile(out,'participant_order.csv'));

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

zeroIQR          = iqr_val == 0;
iqr_val(zeroIQR) = nanstd(numericData(:, zeroIQR), 0, 1);
iqr_val(iqr_val == 0) = 1;

numericData_z = (numericData - med) ./ iqr_val;

fprintf('Standardization: median/IQR applied to %d variables.\n', size(numericData_z,2));
fprintf('  %d variables had zero IQR — fell back to std.\n', sum(zeroIQR));

colMedians = nanmedian(numericData_z, 1);
fprintf('  Post-standardization median range: [%.3f, %.3f]\n', ...
        min(colMedians), max(colMedians));

%% ----------------------------------------------------------
%  6B. SIGN-ALIGN VARIABLES SO ALL POINT IN "WORSE = HIGHER" DIRECTION
%  ---------------------------------------------------------
%  Rationale: Some variables are naturally inversely scaled relative to
%  mitochondrial health (e.g., GDF15 is high when health is poor, while
%  neurological health scores are high when health is good). To make the
%  correlation matrix and heatmap easier to interpret — so that all
%  variables share the convention "higher = worse mitochondrial health" —
%  we flip the sign of any variable where the MitoD group mean is lower
%  than the Control group mean on the standardized data.
%
%  This step is retained because it is part of the user-confirmed analysis.
%  It changes signed correlations and may change Ward clustering based on
%  sqrt(2*(1-r)); it is not merely a cosmetic change to the heatmap.
%  The direction is determined by the observed group mean difference, not
%  an independent clinical definition of whether a measure is worse/better.
% -----------------------------------------------------------
fprintf('\n--- Sign-Aligning Variables (Step 6B) ---\n');

meanPatient = nanmean(numericData_z(patientIdx, :), 1);
meanControl = nanmean(numericData_z(controlIdx, :), 1);
groupDiff   = meanPatient - meanControl;

flipMask = groupDiff < 0;
numericData_z(:, flipMask) = -numericData_z(:, flipMask);

fprintf('  Variables aligned:  %d (no flip needed)\n', sum(~flipMask));
fprintf('  Variables flipped:  %d (sign inverted)\n',  sum(flipMask));
if any(flipMask)
    fprintf('  Flipped variables:\n');
    flippedNames = varNames(flipMask);
    flippedDifferences=groupDiff(flipMask);
    for fi = 1:numel(flippedNames)
        fprintf('    - %s  (MitoD-Control = %.3f -> flipped to +%.3f)\n', ...
                flippedNames{fi}, flippedDifferences(fi), -flippedDifferences(fi));
    end
end
fprintf('  Post-flip difference range: [%.3f, %.3f]\n', ...
        min(groupDiff .* (1 - 2*flipMask)), ...
        max(groupDiff .* (1 - 2*flipMask)));

%% ----------------------------------------------------------
%  7. IMPUTE MISSING VALUES
% -----------------------------------------------------------
numericData_filled = fillmissing(numericData_z, 'linear', 1, 'EndValues','extrap');
numericData_filled = fillmissing(numericData_filled, 'nearest', 1);

nanFrac = sum(isnan(numericData_filled(:))) / numel(numericData_filled);
fprintf('\nRemaining NaN fraction after imputation: %.4f\n', nanFrac);

%% ----------------------------------------------------------
%  8. RESIDUALIZE GROUP FROM EACH VARIABLE
% -----------------------------------------------------------
fprintf('\nResidualizing group membership from each variable...\n');

completeRows = all(~isnan(numericData_filled), 2);
fprintf('Complete cases: %d / %d\n', sum(completeRows), nSubj);

X_input  = numericData_filled(completeRows, :);
G_input  = sortedBinary(completeRows);

G_ranked = tiedrank(G_input);
X_resid  = zeros(size(X_input));
for v = 1:size(X_input, 2)
    X_ranked     = tiedrank(X_input(:, v));
    mdl_v        = fitlm(G_ranked, X_ranked);
    X_resid(:,v) = mdl_v.Residuals.Raw;
end

fprintf('Residualized data: %d subjects x %d variables\n', size(X_resid));

%% ----------------------------------------------------------
%  9A. STANDARD SPEARMAN CORRELATION (for comparison)
% -----------------------------------------------------------
R_pearson = corr(numericData_filled, 'Type','Spearman','Rows','pairwise');
R_pearson = (R_pearson + R_pearson') / 2;
R_pearson(1:size(R_pearson,1)+1:end) = 1;

%% ----------------------------------------------------------
%  9B. PARTIAL CORRELATION — controlling for group membership
% -----------------------------------------------------------
fprintf('\nComputing partial correlation (controlling for group)...\n');

R_partial = partialcorr(X_input, G_input, 'Type', 'Spearman');
R_partial = (R_partial + R_partial') / 2;
R_partial(1:size(R_partial,1)+1:end) = 1;

%% ----------------------------------------------------------
%  9C. SPURIOUS CORRELATION CHECK
% -----------------------------------------------------------
delta_R = R_pearson - R_partial;

figure('Name','Spurious Correlation Check','Position',[50 50 1400 500]);
subplot(1,3,1);
imagesc(R_pearson, [-1 1]); colormap(redblue(256)); colorbar;
title('Spearman Correlation (Raw)');
axis square; set(gca,'XTick',[],'YTick',[]);

subplot(1,3,2);
imagesc(R_partial, [-1 1]); colormap(redblue(256)); colorbar;
title('Partial Spearman Correlation (Group Controlled)');
axis square; set(gca,'XTick',[],'YTick',[]);

subplot(1,3,3);
imagesc(delta_R, [-1 1]); colormap(redblue(256)); colorbar;
title('Difference (Spurious Component)');
axis square; set(gca,'XTick',[],'YTick',[]);
sgtitle('Correlation Comparison: Raw vs Partial');

meanSpurious = mean(abs(delta_R(triu(true(size(delta_R)),1))));
fprintf('Mean absolute spurious correlation: %.4f\n', meanSpurious);

%% ----------------------------------------------------------
%  9D. REMOVE ALL-NaN ROWS/COLS FROM PARTIAL CORR MATRIX
% -----------------------------------------------------------
allNaNRows = all(isnan(R_partial), 2);
allNaNCols = all(isnan(R_partial), 1);
toRemove   = allNaNRows & allNaNCols';

R_clean        = R_partial(~toRemove, ~toRemove);
varNames_clean = varNames(~toRemove);
X_resid_clean  = X_resid(:, ~toRemove);
flipMask_clean = flipMask(~toRemove);

fprintf('Final matrix: %d x %d variables\n', size(R_clean,1), size(R_clean,2));

Rtable_clean = array2table(R_clean, ...
    'VariableNames', varNames_clean, ...
    'RowNames',      varNames_clean);
writetable(Rtable_clean, fullfile(out,'partial_correlation_matrix.csv'), 'WriteRowNames', true);
fprintf('Partial correlation matrix saved.\n');

%% ----------------------------------------------------------
%  10. DISTANCE MATRIX FROM PARTIAL CORRELATION
% -----------------------------------------------------------
D = sqrt(max(0,2 * (1 - R_clean))); % Clamp floating-point roundoff at r=1.
D(1:size(D,1)+1:end) = 0;

if any(isnan(D(:))) || any(isinf(D(:)))
    warning('NaNs/Infs in D — replacing with max finite distance.');
    maxD = max(D(isfinite(D(:))));
    D(isnan(D) | isinf(D)) = maxD;
end

D_vec = squareform(D, 'tovector');

%% ----------------------------------------------------------
%  11. CLUSTERING QUALITY TEST via clusterdata_permtest
% -----------------------------------------------------------
fprintf('\n========== CLUSTERING QUALITY TEST ==========\n');
fprintf('Input: residualized data (%d subjects x %d variables)\n', size(X_resid_clean));

stats_perm=[];idx_k2=[];
if runPermutationTest
    rng(42,'twister'); % Pin this otherwise unseeded diagnostic for new runs.
stats_perm = clusterdata_permtest(X_resid_clean, ...
    'k',              2:40, ...
    'nperm',          100, ...
    'distancemetric', 'correlation', ...
    'linkagemethod',  'ward', ...
    'reducedims',     true, ...
    'verbose',        true, ...
    'doplot',         false);

fprintf('\n--- Permutation Test Summary ---\n');
fprintf('Best k by pseudo-Z:  %d\n',   stats_perm.best_k);
fprintf('Max pseudo-Z:        %.4f\n', stats_perm.max_pseudoZ);

k_values = stats_perm.inputs.k;
idx_k2   = find(k_values == 2);
fprintf('\nResults for k=2 (statistically optimal):\n');
fprintf('  Silhouette score = %.4f\n', stats_perm.cluster_quality(idx_k2));
fprintf('  Null mean        = %.4f\n', stats_perm.cluster_quality_null_mean(idx_k2));
fprintf('  Null std         = %.4f\n', stats_perm.cluster_quality_null_std(idx_k2));
fprintf('  Pseudo-Z         = %.4f\n', stats_perm.cluster_quality_pseudoZ(idx_k2));
fprintf('  P-value          = %.4f\n', stats_perm.P_val(idx_k2));

    save(fullfile(out,'cluster_permutation.mat'),'stats_perm');
else
    fprintf('CANlab clusterdata_permtest unavailable: permutation diagnostic skipped.\n');
end

%% ----------------------------------------------------------
%  12. HIERARCHICAL CLUSTERING ON PARTIAL CORRELATION MATRIX
% -----------------------------------------------------------
Z_link = linkage(D_vec, 'ward');

c = cophenet(Z_link, D_vec);
fprintf('\nCophenetic correlation: %.4f  (>0.75 = good fit)\n', c);

figure('Name','Elbow Plot');
lastMerges = Z_link(end-14:end, 3);
plot(fliplr(2:16), lastMerges, 'o-', 'LineWidth', 2, 'Color', [0.2 0.4 0.8]);
xlabel('Number of clusters'); ylabel('Merge distance');
title('Elbow Plot — Hierarchical Clustering');
grid on;

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

figure('Name','Dendrogram');
dendrogram(Z_link, 0, 'Labels', varNames_clean, 'ColorThreshold','default');
title('Dendrogram (Ward Linkage)');
xtickangle(90);

%% ----------------------------------------------------------
%  13. ASSIGN FINAL CLUSTERS  (k = 11)
% -----------------------------------------------------------
maxClusters = 11;
clusters    = cluster(Z_link, 'maxclust', maxClusters);

fprintf('\n=== Final Cluster Summary (k=%d) ===\n', maxClusters);
for i = 1:maxClusters
    fprintf('  Cluster %d: %d variables\n', i, sum(clusters==i));
end

%% ----------------------------------------------------------
%  CLUSTER LABELS — confirmed against k=11 membership output
%  ---------------------------------------------------------
%  C1  (n=8)  : ln_b2m_copies_ul_serum, all timepoints
%  C2  (n=11) : ln_nd1_copies_ul_plasma + met_feo2/tidal/vo2_kg
%  C3  (n=18) : NAB, RBANS, Trail Making, VF, WASI, CWI, TOPF
%  C4  (n=54) : Body comp, labs, CBC, pulmonary, lifestyle (heterogeneous)
%  C5  (n=16) : ln_nd1 & ln_b2m in saliva, all timepoints
%  C6  (n=6)  : hr_avg across all task windows
%  C7  (n=14) : STAI-state, DSM-5 traits, CTQ, neuroticism, SSQ
%  C8  (n=14) : REE, VO2, ventilation, plasma B2M, early lactate
%  C9  (n=28) : MFIS, BDI, PCL, COMPASS, NMDAS, PFS, SCR
%  C10 (n=15) : ln_nd1_serum, ln_b2m_plasma, late lactate
%  C11 (n=18) : FGF-21 + GDF-15, all timepoints
% -----------------------------------------------------------
clusterLabels = {
    'Serum nDNA';       % C1
    'Plasma mtDNA & Aerobic Capacity';% C2
    'Cognitive Function';                      % C3
    'General Clinical & Metabolic Status';           % C4
    'Saliva nDNA and mtDNA ';            % C5
    'Heart Rate';                     % C6
    'Stress & Psychopathology';            % C7
    'Metabolism & Plasma nDNA & Immediate Lactate';  % C8
    'Multisystem Symptom Burden';              % C9
    'Serum mtDNA & Sustained Lactate';         % C10
    'FGF21/GDF15 Mitokine';  % C11
};

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
%  15. COLOUR PALETTE — 11 maximally distinct colors
% -----------------------------------------------------------
seaborn_colors = [
    0.1216, 0.4667, 0.7059;   %  1  blue
    0.8392, 0.1529, 0.1569;   %  2  red
    0.1725, 0.6275, 0.1725;   %  3  green
    1.0000, 0.4980, 0.0549;   %  4  orange
    0.5804, 0.4039, 0.7412;   %  5  purple
    0.0902, 0.7451, 0.8118;   %  6  cyan
    0.8902, 0.4667, 0.7608;   %  7  pink
    0.5490, 0.3373, 0.2941;   %  8  brown
    0.7373, 0.7412, 0.1333;   %  9  olive
    0.4980, 0.4980, 0.4980;   % 10  gray
    0.9490, 0.7020, 0.0000;   % 11  gold
];

%% ----------------------------------------------------------
%  16. t-SNE PLOT — colored by clusters, with centroid labels
% -----------------------------------------------------------
figure('Name','t-SNE — Partial Correlation Clusters', 'Position',[100 100 1300 900]);
hold on;

for i = 1:maxClusters
    idx_c = (clusters == i);
    scatter(Y_tsne(idx_c,1), Y_tsne(idx_c,2), 80, ...
            'filled', ...
            'MarkerFaceColor', seaborn_colors(i,:), ...
            'DisplayName',     sprintf('C%d: %s  (n=%d)', i, clusterLabels{i}, sum(idx_c)));
end

% Variable name labels
% text(Y_tsne(:,1)+0.15, Y_tsne(:,2), varNames_clean, ...
%     'FontSize', 6, 'Color', [0.35 0.35 0.35]);

% Cluster centroid labels
for i = 1:maxClusters
    idx_c    = (clusters == i);
    cx       = mean(Y_tsne(idx_c, 1));
    cy       = mean(Y_tsne(idx_c, 2));
    labelStr = sprintf('C%d\n%s', i, clusterLabels{i});
    text(cx, cy, labelStr, ...
        'FontSize',            15, ...
        'FontWeight',          'bold', ...
        'Color',               seaborn_colors(i,:) * 0.65, ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor',     [1 1 1], ...
        'EdgeColor',           'none');
end

% legend('Location','eastoutside', 'FontSize', 9);
% xlabel('t-SNE Dimension 1');
% ylabel('t-SNE Dimension 2');
% title({'t-SNE Visualization of Partial Correlation Structure'; ...
%        sprintf('Hierarchical Clusters (k=%d, Ward Linkage)', maxClusters)});
% box off; grid off; hold off;
axis off
saveas(gcf, fullfile(out,'tsne_partial_correlation_clusters.png'));

%% ----------------------------------------------------------
%  17. SORTED PARTIAL CORRELATION HEATMAP WITH CLUSTER BORDERS
% -----------------------------------------------------------
[sortedClusters, sortIdx_cl] = sort(clusters);
R_sorted        = R_clean(sortIdx_cl, sortIdx_cl);
varNames_sorted = varNames_clean(sortIdx_cl);
flipMask_sorted = flipMask_clean(sortIdx_cl);

varLabels_sorted = varNames_sorted;
for vi = 1:numel(varNames_sorted)
    if flipMask_sorted(vi)
        varLabels_sorted{vi} = [char(8224) varNames_sorted{vi}];
    end
end

figure('Name','Sorted Partial Correlation Matrix', 'Position',[100 100 900 800]);
imagesc(R_sorted, [-1 1]);
colormap(redblue(256));
colorbar;
set(gca,'XTick',[],'YTick',[]);

n_vars     = size(R_sorted, 1);
boundaries = find(diff(sortedClusters)) + 0.5;
for b = boundaries(:)'
    line([b b],            [0.5 n_vars+0.5], 'Color','k','LineWidth',2);
    line([0.5 n_vars+0.5], [b b],            'Color','k','LineWidth',2);
end

tickPos  = accumarray(sortedClusters, (1:n_vars)', [], @mean);
tickLbls = arrayfun(@(x) sprintf('C%d: %s', x, clusterLabels{x}), ...
                    1:maxClusters, 'UniformOutput', false);
set(gca, 'YTick',      tickPos, ...
         'YTickLabel', tickLbls, ...
         'FontSize',   15);

title(sprintf('Partial Correlation Matrix — Sorted by Cluster (k=%d)', maxClusters));
annotation('textbox', [0.01 0.01 0.98 0.04], ...
    'String',    [char(8224) ' sign-inverted where the original MitoD mean was lower than Control'], ...
    'EdgeColor', 'none', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
saveas(gcf, fullfile(out,'sorted_partial_correlation_matrix.png'));

%% ----------------------------------------------------------
%  18. REORDERED DISTANCE MATRIX
% -----------------------------------------------------------
figure('Name','Reordered Distance Matrix', 'Position',[100 100 900 800]);
imagesc(D(sortIdx_cl, sortIdx_cl), [0 max(D(:))]);
colormap(flipud(redblue(256)));
colorbar;
set(gca,'XTick',[],'YTick',[]);

for b = boundaries(:)'
    line([b b],            [0.5 n_vars+0.5], 'Color','k','LineWidth',2);
    line([0.5 n_vars+0.5], [b b],            'Color','k','LineWidth',2);
end

set(gca, 'YTick',      tickPos, ...
         'YTickLabel', tickLbls, ...
         'FontSize',   8.5);
title(sprintf('Reordered Distance Matrix — Sorted by Cluster (k=%d)', maxClusters));
saveas(gcf, fullfile(out,'reordered_distance_matrix.png'));

%% ----------------------------------------------------------
%  19. PRINT FINAL SUMMARY FOR RESPONSE LETTER
% -----------------------------------------------------------
fprintf('\n========== SUMMARY FOR RESPONSE LETTER ==========\n');
fprintf('Clustering method:       Hierarchical, Ward linkage\n');
fprintf('Correlation type:        Partial Spearman (group controlled)\n');
fprintf('Permutation test input:  Residualized raw data (group removed)\n');
fprintf('N permutations:          100\n');
fprintf('Cophenetic r:            %.4f\n', c);
fprintf('\nSign alignment:\n');
fprintf('  %d variables flipped: %s\n', sum(flipMask_clean), ...
        strjoin(varNames_clean(flipMask_clean), ', '));
if runPermutationTest
fprintf('\nk=2 results (statistically optimal):\n');
fprintf('  Silhouette score:      %.4f\n', stats_perm.cluster_quality(idx_k2));
fprintf('  Pseudo-Z:              %.4f\n', stats_perm.cluster_quality_pseudoZ(idx_k2));
fprintf('  P-value:               %.4f\n', stats_perm.P_val(idx_k2));
fprintf('\nBest k by pseudo-Z:      %d (Z=%.4f)\n', ...
    stats_perm.best_k, stats_perm.max_pseudoZ);
end
fprintf('\nCluster sizes (k=%d):\n', maxClusters);
for i = 1:maxClusters
    fprintf('  C%d: %s — %d variables\n', i, clusterLabels{i}, sum(clusters==i));
end

if runPermutationTest
figure;
plot(stats_perm.cluster_quality_pseudoZ);
xlabel('k index'); ylabel('Pseudo-Z');
title('Pseudo-Z by k');


end

if runTargetCorrelations
targetVar = 'nmdas_i_ii_iii_score';

targetIdx = find(strcmpi(varNames_clean, targetVar));

if isempty(targetIdx)
    fprintf('Variable "%s" not found. Available names containing "NMDAS":\n', targetVar);
    hits = varNames_clean(contains(varNames_clean, 'NMDAS', 'IgnoreCase', true));
    fprintf('  %s\n', hits{:});
else
    r_target = R_clean(targetIdx, :);

    % Compute p-values from partial correlation
    % df = n - 2 - num_covariates (1 covariate = group membership)
    n_obs = size(X_input, 1);
    df    = n_obs - 2 - 1;
    t_stat = r_target .* sqrt(df) ./ sqrt(1 - r_target.^2);
    p_vals = 2 * (1 - tcdf(abs(t_stat), df));   % two-tailed

    % Preserve the original mafdr default (Storey method, not BH).
    mask      = true(1, numel(p_vals));
    mask(targetIdx) = false;                     % exclude self
    p_sub     = p_vals(mask);
    if exist('mafdr','file')==2
        [p_fdr_sub, ~] = mafdr(p_sub);
        fdrMethod='Original mafdr default';
    else
        p_fdr_sub=nan(size(p_sub));
        fdrMethod='Not calculated: Bioinformatics mafdr unavailable';
        fprintf('%s. Raw correlations and p values are exported.\n',fdrMethod);
    end
    p_fdr     = nan(size(p_vals));
    p_fdr(mask) = p_fdr_sub;

    targetTable=table(string(varNames_clean(:)),r_target(:),p_vals(:),p_fdr(:), ...
        repmat(string(fdrMethod),numel(p_vals),1), ...
        'VariableNames',{'Variable','PartialRho','P','OriginalMafdrOutput','FDRMethod'});
    writetable(targetTable,fullfile(out,'NMDAS_phenotype_correlations.csv'));
    % Sort by absolute correlation (descending)
    [~, sortOrd] = sort(abs(r_target), 'descend');

    fprintf('\n=== Partial correlations with "%s" (n=%d, df=%d) ===\n', ...
            varNames_clean{targetIdx}, n_obs, df);
    fprintf('%-50s  %+6s   %8s   %8s\n', 'Variable', 'r', 'p', 'p_FDR');
    fprintf('%s\n', repmat('-', 1, 80));
    for i = 1:numel(sortOrd)
        j = sortOrd(i);
        if j == targetIdx, continue; end
        fprintf('%-50s  %+.4f   %8.4f   %8.4f\n', ...
                varNames_clean{j}, r_target(j), p_vals(j), p_fdr(j));
    end
end

end

% Export the original panels and explicit data/membership provenance.
signTable=table(string(varNames(:)),groupDiff(:),flipMask(:), ...
    'VariableNames',{'Variable','MitoDMinusControlBeforeFlip','SignInverted'});
writetable(signTable,fullfile(out,'sign_alignment_audit.csv'));
clusterText=string(clusterLabels(clusters));
membership=table(string(varNames_clean(:)),clusters(:),clusterText(:), ...
    flipMask_clean(:),Y_tsne(:,1),Y_tsne(:,2),'VariableNames', ...
    {'Variable','Cluster','OriginalClusterLabel','SignInverted','TSNE1','TSNE2'});
writetable(membership,fullfile(out,'phenotype_clusters.csv'));
save(fullfile(out,'Figure1_analysis.mat'),'R_clean','D','Z_link','clusters', ...
    'Y_tsne','varNames_clean','flipMask_clean','clusterLabels','X_resid_clean', ...
    'stats_perm','numericData_filled','sortedBinary');
figures=findall(groot,'Type','figure');
for figureIndex=1:numel(figures)
    h=figures(figureIndex);label=get(h,'Name');
    if isempty(label),label=sprintf('Figure1_diagnostic_%d',figureIndex);end
    label=regexprep(label,'[^A-Za-z0-9]+','_');
    exportgraphics(h,fullfile(out,[label '.png']),'Resolution',300);
    exportgraphics(h,fullfile(out,[label '.pdf']),'ContentType','vector');
    savefig(h,fullfile(out,[label '.fig']));
end
fprintf('Figure 1 source-aligned outputs: %s\n',out);
fprintf('This source generates the clustering panels; the archived Figure 1A selection is separate.\n');

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

% BEGIN SHARED PACKAGE LOCATOR -- embedded so each analysis remains standalone
function repoRoot = misbieFindPackageRoot(scriptFile, allowDialog, useSavedLocation)
% Find a complete data tree without relying on the MATLAB current directory.
% Select the extracted package, its data folder, or its code folder if asked.
if nargin<2, allowDialog=true; end
if nargin<3, useSavedLocation=true; end
if isempty(scriptFile), scriptDir=pwd; else, scriptDir=fileparts(scriptFile); end
seeds={scriptDir,pwd};
prefGroup='MiSBIE_Reproducibility';prefName='PackageRoot';
if useSavedLocation && ispref(prefGroup,prefName)
    saved=getpref(prefGroup,prefName);
    if ischar(saved) && ~isempty(saved), seeds{end+1}=saved; end
end
candidates=misbieRootCandidates(seeds);
for candidateIndex=1:numel(candidates)
    [ok,~]=misbieHasData(candidates{candidateIndex});
    if ok
        repoRoot=candidates{candidateIndex};
        misbieRememberRoot(repoRoot,prefGroup,prefName,useSavedLocation);
        return
    end
end
if ~allowDialog || ~usejava('awt')
    error('MiSBIE:PackageNotFound', ...
        ['Complete MiSBIE data folder not found. Extract the entire sharing ZIP, ' ...
        'keep code/ and data/ together, and run START_HERE.m in the extracted ' ...
        'MiSBIE_data_code folder. Script location: %s'],scriptDir);
end
fprintf(['\nThe complete MiSBIE data folder was not found beside this script.\n' ...
    'Extract the entire sharing ZIP first. Select MiSBIE_data_code, data, or code.\n' ...
    'The selected location will be remembered for the other MiSBIE scripts.\n']);
selected=uigetdir(scriptDir,'Select the EXTRACTED MiSBIE_data_code folder (or its data/code folder)');
if isequal(selected,0)
    error('MiSBIE:FolderSelectionCancelled', ...
        'Folder selection cancelled. Extract the whole ZIP and run its START_HERE.m.');
end
candidates=misbieRootCandidates({selected});
bestRoot=selected;bestMissing=misbieRequiredDataFiles();
for candidateIndex=1:numel(candidates)
    [ok,missing]=misbieHasData(candidates{candidateIndex});
    if ok
        repoRoot=candidates{candidateIndex};
        misbieRememberRoot(repoRoot,prefGroup,prefName,useSavedLocation);
        return
    end
    if numel(missing)<numel(bestMissing)
        bestRoot=candidates{candidateIndex};bestMissing=missing;
    end
end
error('MiSBIE:IncompletePackage', ...
    ['The selected location does not contain all required data.\n' ...
    'Checked package location: %s\nMissing files:\n  %s\n' ...
    'Extract the entire latest ZIP, including data/. Copying only .m files is insufficient.'], ...
    bestRoot,strjoin(bestMissing,sprintf('\n  ')));
end

function candidates = misbieRootCandidates(seeds)
candidates={};
% Covers code/, package root, data/, and the extra folder Windows extraction creates.
wrappers={'MiSBIE_data_code_GitHub_package','MiSBIE_data_code_GitHub_package_PathFix','MiSBIE_data_code_Figure1_KmeansFix'};
for seedIndex=1:numel(seeds)
    base=char(seeds{seedIndex});
    for level=1:4
        if isempty(base), break; end
        candidates{end+1}=base; %#ok<AGROW>
        candidates{end+1}=fullfile(base,'MiSBIE_data_code'); %#ok<AGROW>
        for wrapperIndex=1:numel(wrappers)
            candidates{end+1}=fullfile(base,wrappers{wrapperIndex},'MiSBIE_data_code'); %#ok<AGROW>
        end
        parent=fileparts(base);
        if strcmp(parent,base), break; end
        base=parent;
    end
end
candidates=unique(candidates,'stable');
end

function [ok,missing] = misbieHasData(repoRoot)
required=misbieRequiredDataFiles();missing={};
for fileIndex=1:numel(required)
    candidate=fullfile(repoRoot,required{fileIndex});
    if ~isfile(candidate)
        missing{end+1}=required{fileIndex}; %#ok<AGROW>
    else
        info=dir(candidate);
        if isempty(info) || info(1).bytes==0
            missing{end+1}=required{fileIndex}; %#ok<AGROW>
        end
    end
end
ok=isempty(missing);
end

function required = misbieRequiredDataFiles()
required={ ...
    'data/participants/analysis_data.csv', ...
    'data/participants/cohort_groups.csv', ...
    'data/participants/structural_summaries.csv', ...
    'data/phenotype/phenotype_scores.csv', ...
    'data/phenotype/clinical_items.csv', ...
    'data/figure2/Paired_condition_scores.mat', ...
    'data/figure_assets/Figure2_background.png', ...
    'data/summary/table_s1_reconciled.csv', ...
    'data/summary/table_s1_original_summary.csv', ...
    'data/summary/archival_Hedgeg_Ke_Preprocessed_WedgePlot.csv', ...
    'data/summary/reported_Table_S2.csv', ...
    'data/summary/reported_Table_S3.csv'};
end

function misbieRememberRoot(repoRoot,prefGroup,prefName,remember)
fprintf('MiSBIE package: %s\n',repoRoot);
if remember
    try
        setpref(prefGroup,prefName,repoRoot);
    catch
        warning('MiSBIE:PreferenceNotSaved', ...
            'Data were found, but MATLAB could not remember the folder for later runs.');
    end
end
end
% END SHARED PACKAGE LOCATOR
