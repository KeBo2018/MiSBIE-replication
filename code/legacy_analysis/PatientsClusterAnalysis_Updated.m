% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
%% ========================================================================
%  PATIENT CLUSTER ANALYSIS - COMPLETE CHARACTERIZATION
%  ------------------------------------------------------------------------
%  Subgroups: resilient (mild) vs severe, defined by k-means (k=2) on NMDAS
%  for the patients included in the N-back analysis (N-back BehOK sample).
%
%  Outputs (all on the same n=25 clustering used for Figure 4D):
%    1. Cluster solution + silhouette coefficient
%    2. Bootstrap stability of cluster assignments (10,000 iterations)
%    3. Age:      mild vs severe (two-sample t-test)
%    4. Sex:      mild vs severe (chi-square + Fisher's exact)
%    5. Genotype: mild vs severe (chi-square)
%    6. GDF15:    mild vs severe (two-sample t-test)  [merged by patient ID]
%    7. Violin plots (numerical variables) + pie charts (categorical)
%
%  Addresses Reviewer 1 Q7 and Reviewer 3 Q3/Q8.
%% ========================================================================

clear; clc;

%% ---- File paths --------------------------------------------------------
clusterFile  = 'Clustered_Patients_Datamatrix.xlsx';   % 25 N-back patients + Patient_Group
behaviorFile = 'DataMatrix_Behavior.xlsx';             % full sample, contains GDF15

T = readtable(clusterFile);

% Group variable: 1 = Resilient (mild), 2 = Severe
groupVar    = 'Patient_Group';
groupLabels = {'Resilient','Severe'};

% Patient ID column (first column, MATLAB auto-names it Var1 if header is blank)
idVar = T.Properties.VariableNames{1};
T.PatientID = erase(string(T.(idVar)), "'");   % strip stray quotes -> 'WITHHELD_SUBJECT' etc.

resilient = T.(groupVar) == 1;
severe    = T.(groupVar) == 2;

fprintf('============================================================\n');
fprintf('  PATIENT CLUSTER CHARACTERIZATION (N-back BehOK sample)\n');
fprintf('============================================================\n');
fprintf('Total patients: %d  (Resilient n=%d, Severe n=%d)\n', ...
        height(T), sum(resilient), sum(severe));


%% ========================================================================
% 1. K-MEANS SOLUTION + SILHOUETTE COEFFICIENT
% ========================================================================
fprintf('\n========== 1. CLUSTER SOLUTION + SILHOUETTE ==========\n');

NMDAS_vec = T.NMDAS;
n         = numel(NMDAS_vec);

% Re-fit k-means on NMDAS to obtain the reference solution
rng(42);
[origIdx, origCent] = kmeans_Matlab(NMDAS_vec, 2, 'Replicates', 50);

% Align labels: cluster 1 = Resilient (lower NMDAS), cluster 2 = Severe
if origCent(1) > origCent(2)
    origIdx  = 3 - origIdx;
    origCent = flip(origCent);
end

% Silhouette coefficient
sil_vals = silhouette(NMDAS_vec, origIdx);
mean_sil = mean(sil_vals);

fprintf('Cluster centers (NMDAS): Resilient=%.2f, Severe=%.2f\n', ...
        origCent(1), origCent(2));
fprintf('Cluster sizes:           Resilient n=%d, Severe n=%d\n', ...
        sum(origIdx==1), sum(origIdx==2));
fprintf('Mean silhouette coefficient = %.3f  (>0.7 = strong separation)\n', mean_sil);

for c = 1:2
    nmd = NMDAS_vec(origIdx==c);
    fprintf('  %s: n=%d, NMDAS mean=%.2f +/- %.2f, range [%.1f, %.1f]\n', ...
        groupLabels{c}, numel(nmd), mean(nmd), std(nmd), min(nmd), max(nmd));
end


%% ========================================================================
% 2. BOOTSTRAP STABILITY OF CLUSTER ASSIGNMENTS
% ========================================================================
fprintf('\n========== 2. BOOTSTRAP STABILITY (10,000 iterations) ==========\n');

rng(42);
nIter = 10000;

stabilityCount = zeros(n, 1);   % times each subject matches original assignment
totalCount     = zeros(n, 1);   % times each subject appears in a bootstrap sample

for i = 1:nIter
    bootIdx    = randsample(n, n, true);
    NMDAS_boot = NMDAS_vec(bootIdx);

    [bootLabels, bootCent] = kmeans_Matlab(NMDAS_boot, 2, 'Replicates', 10);

    % Align labels: cluster 1 = lower NMDAS (Resilient)
    if bootCent(1) > bootCent(2)
        bootLabels = 3 - bootLabels;
    end

    for j = 1:n
        subj                 = bootIdx(j);
        totalCount(subj)     = totalCount(subj) + 1;
        if bootLabels(j) == origIdx(subj)
            stabilityCount(subj) = stabilityCount(subj) + 1;
        end
    end
end

subjectStability = stabilityCount ./ totalCount;

fprintf('Per-subject assignment stability:\n');
fprintf('  Mean = %.1f%%,  Min = %.1f%%,  Max = %.1f%%\n', ...
        100*mean(subjectStability), ...
        100*min(subjectStability), ...
        100*max(subjectStability));


%% ========================================================================
% 3. AGE: mild vs severe
% ========================================================================
fprintf('\n========== 3. AGE: Resilient vs Severe ==========\n');

age_res = T.Age(resilient);  age_res = age_res(~isnan(age_res));
age_sev = T.Age(severe);     age_sev = age_sev(~isnan(age_sev));

[~, p_age, ci_age, st_age] = ttest2(age_res, age_sev);

fprintf('Resilient: %.2f +/- %.2f years (n=%d)\n', ...
        mean(age_res), std(age_res), numel(age_res));
fprintf('Severe:    %.2f +/- %.2f years (n=%d)\n', ...
        mean(age_sev), std(age_sev), numel(age_sev));
fprintf('Two-sample t-test: t(%d)=%.3f, p=%.4f, 95%% CI [%.3f, %.3f]\n', ...
        st_age.df, st_age.tstat, p_age, ci_age(1), ci_age(2));


%% ========================================================================
% 4. SEX: mild vs severe  (chi-square + Fisher's exact)
% ========================================================================
fprintf('\n========== 4. SEX: Resilient vs Severe ==========\n');

sex_all = categorical(string(T.sex));
sexCats = categories(sex_all);

% 2 x nSexCats contingency table
sexTable = zeros(2, numel(sexCats));
for c = 1:numel(sexCats)
    sexTable(1,c) = sum(sex_all(resilient) == sexCats{c});
    sexTable(2,c) = sum(sex_all(severe)    == sexCats{c});
end

disp(array2table(sexTable, 'RowNames', groupLabels, 'VariableNames', sexCats));

% Chi-square (manual)
[chi2_sex, df_sex, p_sex_chi, expMin_sex] = chi2_contingency(sexTable);
fprintf('Chi-square: chi2(%d)=%.3f, p=%.4f  (min expected cell = %.2f)\n', ...
        df_sex, chi2_sex, p_sex_chi, expMin_sex);

% Fisher's exact test (valid for 2x2 with small expected counts)
if isequal(size(sexTable), [2 2])
    [~, p_sex_fisher] = fishertest(sexTable);
    fprintf('Fisher''s exact test: p=%.4f  (RECOMMENDED - small expected counts)\n', ...
            p_sex_fisher);
end


%% ========================================================================
% 5. GENOTYPE: mild vs severe
% ========================================================================
fprintf('\n========== 5. GENOTYPE: Resilient vs Severe ==========\n');

% Variable may be named 'geneticDiagnosis' or 'genetic diagnosis' -> auto-detect
genCol = T.Properties.VariableNames( ...
            contains(lower(T.Properties.VariableNames), 'genetic'));
gen_all = categorical(string(T.(genCol{1})));
genCats = categories(gen_all);

genTable = zeros(2, numel(genCats));
for c = 1:numel(genCats)
    genTable(1,c) = sum(gen_all(resilient) == genCats{c});
    genTable(2,c) = sum(gen_all(severe)    == genCats{c});
end

disp(array2table(genTable, 'RowNames', groupLabels, 'VariableNames', genCats));

[chi2_gen, df_gen, p_gen, expMin_gen] = chi2_contingency(genTable);
fprintf('Chi-square: chi2(%d)=%.3f, p=%.4f  (min expected cell = %.2f)\n', ...
        df_gen, chi2_gen, p_gen, expMin_gen);
if expMin_gen < 5
    fprintf('  NOTE: min expected count < 5 - interpret chi-square with caution.\n');
end


%% ========================================================================
% 6. GDF15: mild vs severe   (merged from DataMatrix_Behavior.xlsx by ID)
% ========================================================================
fprintf('\n========== 6. GDF15: Resilient vs Severe ==========\n');

B = readtable(behaviorFile);

% Match patient IDs between the two files
B.PatientID = erase(string(B.Var1), "'");

% GDF15 column in the behavior file
GDF15_raw = B.t_fasting_plasma;
GDF15_log_B = log10(GDF15_raw);

% Map GDF15 onto the 25 clustered patients by ID
T.GDF15_log = nan(height(T), 1);
for k = 1:height(T)
    hit = find(B.PatientID == T.PatientID(k), 1);
    if ~isempty(hit)
        T.GDF15_log(k) = GDF15_log_B(hit);
    end
end

gdf_res = T.GDF15_log(resilient);  gdf_res = gdf_res(~isnan(gdf_res));
gdf_sev = T.GDF15_log(severe);     gdf_sev = gdf_sev(~isnan(gdf_sev));

if numel(gdf_res) > 1 && numel(gdf_sev) > 1
    [~, p_gdf, ci_gdf, st_gdf] = ttest2(gdf_sev, gdf_res);  % severe vs mild -> positive t
    fprintf('Resilient: log GDF15 = %.3f +/- %.3f (n=%d)\n', ...
            mean(gdf_res), std(gdf_res), numel(gdf_res));
    fprintf('Severe:    log GDF15 = %.3f +/- %.3f (n=%d)\n', ...
            mean(gdf_sev), std(gdf_sev), numel(gdf_sev));
    fprintf('Two-sample t-test (Severe vs Resilient): t(%d)=%.3f, p=%.4f, 95%% CI [%.3f, %.3f]\n', ...
            st_gdf.df, st_gdf.tstat, p_gdf, ci_gdf(1), ci_gdf(2));
else
    fprintf('Not enough GDF15 data to run the comparison.\n');
    fprintf('  Resilient with GDF15: n=%d, Severe with GDF15: n=%d\n', ...
            numel(gdf_res), numel(gdf_sev));
end


%% ========================================================================
% 7. FIGURES: violin plots (numerical) + pie charts (categorical)
% ========================================================================
LineWidth = 2;  Fontsize = 12;  pointsize = 4;
dotcolor1 = [255, 212, 159] / 255;
dotcolor2 = [203, 104, 104] / 255;
colorcoding = [dotcolor1; dotcolor2];

% ---- Violin plots: numerical variables ----
numVars = {'Age','Nback_Behavioral_Acc','GDF15_log','cns_eyes','pfs_phys','pfs_mental'};
numVars = numVars(ismember(numVars, T.Properties.VariableNames));

nVars = numel(numVars);
nCols = 3;
nRows = ceil(nVars / nCols);

figure('Color','w','Position',[200 100 350*nCols 320*nRows]);
for i = 1:nVars
    subplot(nRows, nCols, i);

    y_res = T{resilient, numVars{i}};  y_res = y_res(~isnan(y_res));
    y_sev = T{severe,    numVars{i}};  y_sev = y_sev(~isnan(y_sev));

    violinplot({y_res, y_sev}, groupLabels, ...
        'mc','k', 'facecolor', colorcoding, ...
        'plotlegend', 0, 'pointsize', pointsize);
    title(strrep(numVars{i},'_',' '), 'FontSize', 14);

    if numel(y_res) > 1 && numel(y_sev) > 1
        [~, p] = ttest2(y_res, y_sev);
    else
        p = NaN;
    end

    if ~isnan(p) && p < 0.05
        yl   = ylim;
        ySig = max([y_res; y_sev]) + 0.08 * range(yl);
        hold on; plot([1 2], [ySig ySig], 'k-', 'LineWidth', 1);
        if     p < 0.001, sigLabel = '***';
        elseif p < 0.01,  sigLabel = '**';
        else              sigLabel = '*';
        end
        text(1.5, ySig + 0.02*range(yl), sigLabel, ...
            'HorizontalAlignment','center', 'FontSize', 12);
        hold off;
    end
end
sgtitle('Numerical variables: Resilient vs Severe (two-sample t-tests)', 'FontSize', 14);
saveas(gcf, 'cluster_characterization_violins.png');

% ---- Pie charts: categorical variables ----
catVars = {'sex'};
genColName = T.Properties.VariableNames( ...
                contains(lower(T.Properties.VariableNames),'genetic'));
catVars = [catVars, genColName];
catVars = catVars(ismember(catVars, T.Properties.VariableNames));

for v = 1:numel(catVars)
    figure('Color','w','Position',[200 200 900 400]);
    for g = 1:2
        subplot(1,2,g);
        grpData = categorical(string(T{T.(groupVar)==g, catVars{v}}));
        grpData = grpData(~ismissing(grpData));
        h = pie(countcats(grpData), categories(grpData));
        set(findobj(h,'Type','text'),'FontSize',14);
        title(groupLabels{g}, 'FontSize', 16);
    end
    sgtitle(sprintf('%s distribution by subgroup', strrep(catVars{v},'_',' ')), ...
            'FontSize', 18);
    saveas(gcf, sprintf('cluster_characterization_pie_%s.png', catVars{v}));
end

fprintf('\n============================================================\n');
fprintf('  END OF CLUSTER CHARACTERIZATION\n');
fprintf('============================================================\n');


%% ========================================================================
% LOCAL FUNCTION: manual chi-square test for a contingency table
% ========================================================================
function [chi2, df, p, expMin] = chi2_contingency(observed)
    rowSums = sum(observed, 2);
    colSums = sum(observed, 1);
    nTot    = sum(observed(:));
    expected = (rowSums * colSums) / nTot;
    chi2 = sum((observed(:) - expected(:)).^2 ./ expected(:));
    df   = (size(observed,1)-1) * (size(observed,2)-1);
    p    = 1 - chi2cdf(chi2, df);
    expMin = min(expected(:));
end
