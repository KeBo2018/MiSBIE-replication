% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
%% ========================================================================
%  MASTER SCRIPT: Reproduce all major results from MiSBIE paper
%  Author: Ke Bo | 2025
%  This script:
%   1. Loads unified data table
%   2. Cleans variables
%   3. Runs all major analyses:
%       - Group comparisons
%       - Disease severity correlations (patients only)
%       - N-back behavior correlations
%       - GDF15 correlations + regressions
%       - K-means subgroup analysis
%   4. Produces all key plots automatically
% =========================================================================

clear; clc; close all;

%% =======================
% LOAD DATA
% ========================

T = readtable('LOCAL_SOURCE\Dartmouth College Dropbox\Ke Bo\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\Manuscript\Submission2_Science\Manuscript revision\Analysis\DataMatrix_Behavior.xlsx');

% Clean MiID format
% T.MiID = erase(string(T.MiID), "'");

% Extract variables
Multi  = T.Multisensory;
Cold   = T.Cold;
Nback  = T.Nback;
Age    = T.Age;
Group  = T.GroupID;                  % 0 = patient, 1 = control
NMDAS  = T.NMDAS;
BehAcc = T.Nback_ACC_;
BehOK  = T.x_Behavioral_Available_ == 1;
GDF15  = T.t_fasting_plasma;
GDF15_log = log10(GDF15);            % matches paper


%% ========================================================================
% 1. GROUP COMPARISONS (Multisensory, Cold, N-back)
% ========================================================================

fprintf('\n================ GROUP COMPARISONS ================\n');

print_ttest(Multi(Group==0), Multi(Group==1), 'Multisensory');
print_ttest(Cold(Group==0), Cold(Group==1), 'Cold');
print_ttest(Nback(Group==0 & BehOK), Nback(Group==1 & BehOK), 'N-back');

%% =======================
% GROUP COMPARISONS
% =======================

fprintf('\n================ GROUP COMPARISONS ================\n');

print_ttest(Multi(Group==0), Multi(Group==1), 'Multisensory');
print_ttest(Cold(Group==0), Cold(Group==1), 'Cold');
print_ttest(Nback(Group==0 & BehOK), Nback(Group==1 & BehOK), 'N-back');

figure; 
subplot(1,3,1); boxplot(Multi, Group); title('Multisensory'); ylabel('Brain Activation');
subplot(1,3,2); boxplot(Cold, Group); title('Cold Pain');
subplot(1,3,3); boxplot(Nback(BehOK), Group(BehOK)); title('N-back');
sgtitle('Group Comparisons: Patients (0) vs Controls (1)');
saveas(gcf, 'group_comparisons.png');

%% ========================================================================
% 2. CORRELATIONS WITH DISEASE SEVERITY (patients only)
% ========================================================================

pat = Group == 0;

fprintf('\n================ NMDAS CORRELATIONS (Age-corrected) ================\n');

% Make sure to align vectors and use complete rows
valid1 = pat & ~isnan(Multi) & ~isnan(NMDAS) & ~isnan(Age);
[r1,p1] = partialcorr(Multi(valid1), NMDAS(valid1), Age(valid1), ...
                      'Type','Spearman','Rows','complete');

valid2 = pat & ~isnan(Cold) & ~isnan(NMDAS) & ~isnan(Age);
[r2,p2] = partialcorr(Cold(valid2), NMDAS(valid2), Age(valid2), ...
                      'Type','Spearman','Rows','complete');

valid3 = pat & BehOK & ~isnan(Nback) & ~isnan(NMDAS) & ~isnan(Age);
[r3,p3] = partialcorr(Nback(valid3), NMDAS(valid3), Age(valid3), ...
                      'Type','Spearman','Rows','complete');

fprintf('NMDAS vs Multisensory (age-corrected): r=%.3f, p=%.4f, n=%d\n', r1, p1, sum(valid1));
fprintf('NMDAS vs Cold (age-corrected):         r=%.3f, p=%.4f, n=%d\n', r2, p2, sum(valid2));
fprintf('NMDAS vs N-back (age-corrected):       r=%.3f, p=%.4f, n=%d\n', r3, p3, sum(valid3));

% Plot
figure;
subplot(1,3,1); scatter(NMDAS(pat), Multi(pat),'filled'); lsline; xlabel('NMDAS'); ylabel('Multisensory');
subplot(1,3,2); scatter(NMDAS(pat), Cold(pat),'filled'); lsline; xlabel('NMDAS'); ylabel('Cold');
subplot(1,3,3); scatter(NMDAS(pat & BehOK), Nback(pat & BehOK),'filled'); lsline; xlabel('NMDAS'); ylabel('N-back');
sgtitle('Disease Severity Correlations');
saveas(gcf, 'severity_correlations.png');


%% ========================================================================
% 3. N-BACK BEHAVIOR CORRELATION
% ========================================================================

[rB,pB] = partialcorr(NMDAS(valid3), BehAcc(valid3),Age(valid3) ,'Type','spearman','Rows','complete');
fprintf('\nN-back brain vs behavior: r=%.3f, p=%.4f\n', rB, pB);

figure;
scatter(Nback(BehOK), BehAcc(BehOK),'filled');
lsline; xlabel('N-back Brain Activation'); ylabel('Behavior Accuracy');
title(sprintf('Brain–Behavior Correlation (r=%.2f, p=%.4f)', rB,pB));
saveas(gcf, 'brain_behavior_correlation.png');


%% ========================================================================
% 4. GDF15 CORRELATIONS
% ========================================================================

fprintf('\n================ GDF15 CORRELATIONS ================\n');

% Whole sample
[rG,pG] = corr(GDF15_log(BehOK), Nback(BehOK), 'Type','Pearson','Rows','complete');
fprintf('Whole-sample GDF15-log vs N-back: r=%.3f, p=%.4f\n', rG,pG);

% Patients only
[rGp,pGp] = corr(GDF15_log(pat & BehOK), Nback(pat & BehOK), 'Type','Pearson','Rows','complete');
fprintf('Patients: r=%.3f, p=%.4f\n', rGp, pGp);

% Controls only
con = Group == 1;
[rGc,pGc] = corr(GDF15_log(con & BehOK), Nback(con & BehOK), 'Type','Pearson','Rows','complete');
fprintf('Controls: r=%.3f, p=%.4f\n', rGc, pGc);

% Plot
figure;
scatter(GDF15_log(BehOK), Nback(BehOK), 60, Group(BehOK), 'filled');
lsline; xlabel('log10(GDF15)'); ylabel('N-back Brain Activation');
title('GDF15 vs N-back Brain Activation');
colormap([1 0.5 0; 0 0.5 1]); colorbar('Ticks',[0,1],'TickLabels',{'Patient','Control'});
saveas(gcf, 'GDF15_brain_correlation.png');


%% ========================================================================
% 4B. BH-FDR CORRECTION ACROSS PRIMARY CORRELATIONS
% ========================================================================
% Addresses Reviewer Question 7: multiple comparisons correction.
%
% Two test families are corrected separately, matching the manuscript's
% pre-specified hypotheses:
%
%   Family 1 — NMDAS correlations (patients only, age-corrected Spearman):
%     (1) NMDAS vs. N-back brain
%     (2) NMDAS vs. N-back behavior
%     (3) NMDAS vs. Multisensory brain
%     (4) NMDAS vs. Cold brain
%
%   Family 2 — GDF15 ~ N-back brain (Pearson):
%     (5) Whole sample
%     (6) Patients only
%     (7) Controls only
%
%   Family 3 — GDF15 ~ N-back behavior (Pearson):
%     (8) Whole sample
%     (9) Patients only
%    (10) Controls only
%
% Additionally, all 10 tests are corrected together as the most
% conservative option (Option B).
%
% NOTE: r3, p3 = NMDAS vs N-back brain    (from Section 2)
%       r1, p1 = NMDAS vs Multisensory     (from Section 2)
%       r2, p2 = NMDAS vs Cold             (from Section 2)
%       rB, pB = N-back brain vs behavior  (from Section 3)
%       rG, pG = GDF15 vs N-back (whole)   (from Section 4)
%       rGp,pGp= GDF15 vs N-back (patients)(from Section 4)
%       rGc,pGc= GDF15 vs N-back (controls)(from Section 4)
%
% GDF15 ~ behavior p-values need to be computed here:
% ========================================================================

fprintf('\n========== BH-FDR CORRECTION: PRIMARY CORRELATIONS ==========\n');

% --- GDF15 ~ N-back BEHAVIOR correlations (not in Section 4, computed here) ---
valid_beh_all = BehOK & ~isnan(GDF15_log) & ~isnan(BehAcc);
[rGB,  pGB ] = corr(GDF15_log(valid_beh_all), BehAcc(valid_beh_all), ...
                    'Type','Pearson','Rows','complete');

valid_beh_pat = pat & BehOK & ~isnan(GDF15_log) & ~isnan(BehAcc);
[rGBp, pGBp] = corr(GDF15_log(valid_beh_pat), BehAcc(valid_beh_pat), ...
                    'Type','Pearson','Rows','complete');

valid_beh_con = con & BehOK & ~isnan(GDF15_log) & ~isnan(BehAcc);
[rGBc, pGBc] = corr(GDF15_log(valid_beh_con), BehAcc(valid_beh_con), ...
                    'Type','Pearson','Rows','complete');

fprintf('GDF15 vs N-back Behavior (whole):    r=%.3f, p=%.4f\n', rGB,  pGB);
fprintf('GDF15 vs N-back Behavior (patients): r=%.3f, p=%.4f\n', rGBp, pGBp);
fprintf('GDF15 vs N-back Behavior (controls): r=%.3f, p=%.4f\n', rGBc, pGBc);

% --- Assemble test families ---
% Labels, r-values, raw p-values for each family
fam1_labels = {'NMDAS vs N-back Brain (age-corr)', ...
               'NMDAS vs N-back Behavior (age-corr)', ...
               'NMDAS vs Multisensory Brain (age-corr)', ...
               'NMDAS vs Cold Brain (age-corr)'};
fam1_r = [r3;  rB;  r1;  r2];
fam1_p = [p3;  pB;  p1;  p2];

fam2_labels = {'GDF15 vs N-back Brain (whole)', ...
               'GDF15 vs N-back Brain (patients)', ...
               'GDF15 vs N-back Brain (controls)'};
fam2_r = [rG;  rGp; rGc];
fam2_p = [pG;  pGp; pGc];

fam3_labels = {'GDF15 vs N-back Behavior (whole)', ...
               'GDF15 vs N-back Behavior (patients)', ...
               'GDF15 vs N-back Behavior (controls)'};
fam3_r = [rGB;  rGBp; rGBc];
fam3_p = [pGB;  pGBp; pGBc];

% --- BH-FDR function (inline) ---
% Returns q-values using Benjamini-Hochberg procedure
bh_fdr = @(p) bh_fdr_fn(p);

% --- Apply FDR within each family ---
fam1_q = bh_fdr(fam1_p);
fam2_q = bh_fdr(fam2_p);
fam3_q = bh_fdr(fam3_p);

% --- Option B: all 10 tests together ---
all_labels = [fam1_labels, fam2_labels, fam3_labels];
all_r      = [fam1_r; fam2_r; fam3_r];
all_p      = [fam1_p; fam2_p; fam3_p];
all_q_B    = bh_fdr(all_p);

% --- Print results ---
fprintf('\n--- FAMILY 1: NMDAS correlations (n=%d tests) ---\n', length(fam1_p));
fprintf('%-45s  %7s  %7s  %7s  %s\n', 'Test', 'r', 'p', 'q(BH)', 'Sig(q<0.05)');
for i = 1:length(fam1_p)
    sig = '';
    if fam1_q(i) < 0.05; sig = '*** FDR-sig'; end
    fprintf('  %-43s  %+6.3f  %7.4f  %7.4f  %s\n', ...
        fam1_labels{i}, fam1_r(i), fam1_p(i), fam1_q(i), sig);
end

fprintf('\n--- FAMILY 2: GDF15 ~ N-back Brain (n=%d tests) ---\n', length(fam2_p));
fprintf('%-45s  %7s  %7s  %7s  %s\n', 'Test', 'r', 'p', 'q(BH)', 'Sig(q<0.05)');
for i = 1:length(fam2_p)
    sig = '';
    if fam2_q(i) < 0.05; sig = '*** FDR-sig'; end
    fprintf('  %-43s  %+6.3f  %7.4f  %7.4f  %s\n', ...
        fam2_labels{i}, fam2_r(i), fam2_p(i), fam2_q(i), sig);
end

fprintf('\n--- FAMILY 3: GDF15 ~ N-back Behavior (n=%d tests) ---\n', length(fam3_p));
fprintf('%-45s  %7s  %7s  %7s  %s\n', 'Test', 'r', 'p', 'q(BH)', 'Sig(q<0.05)');
for i = 1:length(fam3_p)
    sig = '';
    if fam3_q(i) < 0.05; sig = '*** FDR-sig'; end
    fprintf('  %-43s  %+6.3f  %7.4f  %7.4f  %s\n', ...
        fam3_labels{i}, fam3_r(i), fam3_p(i), fam3_q(i), sig);
end

fprintf('\n--- OPTION B: All %d tests corrected together (most conservative) ---\n', length(all_p));
fprintf('%-45s  %7s  %7s  %7s  %s\n', 'Test', 'r', 'p', 'q(BH)', 'Sig(q<0.05)');
for i = 1:length(all_p)
    sig = '';
    if all_q_B(i) < 0.05; sig = '*** FDR-sig'; end
    fprintf('  %-43s  %+6.3f  %7.4f  %7.4f  %s\n', ...
        all_labels{i}, all_r(i), all_p(i), all_q_B(i), sig);
end

fprintf('\n=============================================================================\n');

%% ========================================================================
% 5. MULTIPLE REGRESSION: Group + Age + GDF15
% ========================================================================

valid = BehOK & ~isnan(GDF15_log) & ~isnan(Age) & ~isnan(Nback);

Y = zscore(Nback(valid));
X = [zscore(GDF15_log(valid)), zscore(Age(valid)), zscore(Group(valid)), ones(sum(valid),1)];

[b,~,~,~,stats] = regress(Y, X);

fprintf('\nMultiple regression predicting N-back brain activation:\n');
fprintf('GDF15 p = %.4f\n', stats(3));  % p-value of whole model
disp('Coefficient order: [GDF15, Age, Group, Constant]');
disp(b);


% %% ========================================================================
% % 6. K-MEANS SUBGROUP (Mild vs Severe patients)
% % ========================================================================
% 
% % K-means clustering on NMDAS (patients)
[idxC,cent] = kmeans_Matlab(NMDAS(pat), 2);

% Map back to full size 
idx_full = nan(size(Group));
idx_full(pat) = idxC;

mild   = idx_full == 1 & pat & BehOK;
severe = idx_full == 2 & pat & BehOK;

[H,P,CI,STATS] = ttest2(Nback(mild), Nback(severe));

mild   = idx_full == 1 & pat ;
severe = idx_full == 2 & pat ;

[H,P,CI,STATS] = ttest2(Cold(mild), Cold(severe));

mild   = idx_full == 1 & pat ;
severe = idx_full == 2 & pat ;

[H,P,CI,STATS] = ttest2(Multi(mild), Multi(severe));

fprintf('\nK-Means Subgroups (N-back): t(%d)=%.3f, p=%.4f\n', ...
    STATS.df, STATS.tstat, P);

figure;
boxplot([Nback(mild); Nback(severe)], [zeros(sum(mild),1); ones(sum(severe),1)]);
title('Mild vs Severe MitoD (N-back brain)');
set(gca,'XTickLabel',{'Mild','Severe'});
ylabel('N-back Brain Activation');
saveas(gcf, 'kmeans_mild_vs_severe.png');
%% ========================================================================
% PERMUTATION TEST: NMDAS vs N-back (patients only, age-controlled)
% ========================================================================
% Addresses Reviewer 2 Q5 - permutation-based validation of Spearman correlation

rng(42);
nPerm = 10000;

% Use same valid index as Section 2
brain_data  = Nback(valid3);
nmdas_data  = NMDAS(valid3);
age_data    = Age(valid3);
n_valid     = sum(valid3);

% Observed partial correlation (already computed above as r3)
r_obs = r3;

% Permutation loop - permute NMDAS labels only
r_perm = nan(nPerm, 1);
for i = 1:nPerm
    nmdas_perm = nmdas_data(randperm(n_valid));
    r_perm(i)  = partialcorr(brain_data, nmdas_perm, age_data, ...
                             'Type','Spearman','Rows','complete');
end

% Permutation p-value (two-tailed)
p_perm = mean(abs(r_perm) >= abs(r_obs));

fprintf('\n========== PERMUTATION TEST RESULTS ==========\n');
fprintf('Observed partial r = %.3f\n', r_obs);
fprintf('Permutation p      = %.4f (based on %d permutations)\n', p_perm, nPerm);
fprintf('(Parametric p      = %.4f)\n', p3);

%% ========================================================================
% 2C. SPLIT-HALF ROBUSTNESS: NMDAS correlations (patients only)
% ========================================================================
% Uses same sample and variables as Section 2 above
% Addresses Reviewer 2 Q6

rng(42);
nIter = 1000;

% Use same valid indices as Section 2
valid_brain = valid3;  % pat & BehOK & ~isnan(Nback) & ~isnan(NMDAS) & ~isnan(Age)
valid_beh   = pat & BehOK & ~isnan(BehAcc) & ~isnan(NMDAS) & ~isnan(Age);

% Convert to indices
brain_idx = find(valid_brain);
beh_idx   = find(valid_beh);

n_brain = length(brain_idx);
n_beh   = length(beh_idx);

rSplit_brain = nan(nIter,1);
rSplit_beh   = nan(nIter,1);

% --- Brain loop ---
for i = 1:nIter
    idx_shuf = brain_idx(randperm(n_brain));
    half1 = idx_shuf(1:floor(n_brain/2));
    half2 = idx_shuf(floor(n_brain/2)+1:end);

    r1 = partialcorr(Nback(half1), NMDAS(half1), Age(half1), ...
                     'Type','Spearman','Rows','complete');
    r2 = partialcorr(Nback(half2), NMDAS(half2), Age(half2), ...
                     'Type','Spearman','Rows','complete');
    rSplit_brain(i) = mean([r1, r2]);
end

% --- Behavior loop ---
for i = 1:nIter
    idx_shuf = beh_idx(randperm(n_beh));
    half1 = idx_shuf(1:floor(n_beh/2));
    half2 = idx_shuf(floor(n_beh/2)+1:end);

    r1 = partialcorr(BehAcc(half1), NMDAS(half1), Age(half1), ...
                     'Type','Spearman','Rows','complete');
    r2 = partialcorr(BehAcc(half2), NMDAS(half2), Age(half2), ...
                     'Type','Spearman','Rows','complete');
    rSplit_beh(i) = mean([r1, r2]);
end

% --- Print ---
fprintf('\n========== SPLIT-HALF ROBUSTNESS (age-controlled) ==========\n');

fprintf('\nNMDAS vs Brain Activation (N-back):\n');
fprintf('  Full-sample partial r = %.3f\n', r3);   % already computed above
fprintf('  Median split-half r   = %.3f\n', median(rSplit_brain));
fprintf('  95%% CI               = [%.3f, %.3f]\n', ...
    prctile(rSplit_brain,2.5), prctile(rSplit_brain,97.5));
fprintf('  %% same direction     = %.1f%%\n', ...
    100 * mean(rSplit_brain < 0));

fprintf('\nNMDAS vs Behavioral Performance:\n');
fprintf('  Full-sample r         = %.3f (p=%.4f)\n', ...
    corr(BehAcc(valid_beh), NMDAS(valid_beh), 'Type','Spearman','Rows','complete'), ...
    corr(BehAcc(valid_beh), NMDAS(valid_beh), 'Type','Spearman','Rows','complete'));
fprintf('  Median split-half r   = %.3f\n', median(rSplit_beh));
fprintf('  95%% CI               = [%.3f, %.3f]\n', ...
    prctile(rSplit_beh,2.5), prctile(rSplit_beh,97.5));
fprintf('  %% same direction     = %.1f%%\n', ...
    100 * mean(rSplit_beh < 0));

%%%%%%%%%%% Motion regression %%%%

%% ========================================================================
% 1B. GROUP COMPARISONS (Correcting for Motion)
% ========================================================================

fprintf('\n================ GROUP COMPARISONS (Motion-corrected) ================\n');

%% ========================================================================
% 1B. GROUP COMPARISONS Using fitlm (Motion + Age controlled)
% ========================================================================

fprintf('\n================ GROUP COMPARISONS (fitlm, motion-controlled) ================\n');

%% -------- MULTISENSORY --------------------------------------------------
valid = ~isnan(Multi) & ~isnan(Group) & ~isnan(T.multisensory_meanFD) & ~isnan(Age);

tbl = table( ...
    Multi(valid), ...
    Group(valid), ...
    T.multisensory_meanFD(valid), ...
    Age(valid), ...
    'VariableNames', {'Y','Group','FD','Age'} ...
);

mdl_multi = fitlm(tbl, 'Y ~ Group + FD + Age');
disp(mdl_multi)

fprintf('Multisensory group effect: t=%.3f, p=%.4f\n', ...
    mdl_multi.Coefficients.tStat('Group'), ...
    mdl_multi.Coefficients.pValue('Group') );


%% -------- COLD ----------------------------------------------------------
valid = ~isnan(Cold) & ~isnan(Group) & ~isnan(T.coldpressor_meanFD) & ~isnan(Age);

tbl = table( ...
    Cold(valid), ...
    Group(valid), ...
    T.coldpressor_meanFD(valid), ...
    Age(valid), ...
    'VariableNames', {'Y','Group','FD','Age'} ...
);

mdl_cold = fitlm(tbl, 'Y ~ Group + FD + Age');
disp(mdl_cold)

fprintf('Cold group effect: t=%.3f, p=%.4f\n', ...
    mdl_cold.Coefficients.tStat('Group'), ...
    mdl_cold.Coefficients.pValue('Group') );


%% -------- N-BACK --------------------------------------------------------
valid = BehOK & ~isnan(Nback) & ~isnan(Group) & ~isnan(T.nback_meanFD) & ~isnan(Age);

tbl = table( ...
    Nback(valid), ...
    Group(valid), ...
    T.nback_meanFD(valid), ...
    Age(valid), ...
    'VariableNames', {'Y','Group','FD','Age'} ...
);

mdl_nback = fitlm(tbl, 'Y ~ Group + FD + Age');
disp(mdl_nback)

fprintf('N-back group effect: t=%.3f, p=%.4f\n', ...
    mdl_nback.Coefficients.tStat('Group'), ...
    mdl_nback.Coefficients.pValue('Group') );


%% ========================================================================
% 2B. NMDAS CORRELATIONS Using fitlm (Age + Motion controlled)
% ========================================================================

fprintf('\n================ NMDAS CORRELATIONS (fitlm, age + motion) ================\n');

%% -------- MULTISENSORY --------------------------------------------------
valid = pat & ~isnan(Multi) & ~isnan(NMDAS) & ...
         ~isnan(Age) & ~isnan(T.multisensory_meanFD);

tbl = table( ...
    Multi(valid), ...
    NMDAS(valid), ...
    Age(valid), ...
    T.multisensory_meanFD(valid), ...
    'VariableNames', {'Y','NMDAS','Age','FD'} ...
);

mdl_multi_sev = fitlm(tbl, 'Y ~ NMDAS + Age + FD');
disp(mdl_multi_sev)

fprintf('NMDAS → Multisensory: t=%.3f, p=%.4f\n', ...
    mdl_multi_sev.Coefficients.tStat('NMDAS'), ...
    mdl_multi_sev.Coefficients.pValue('NMDAS') );


%% -------- COLD ----------------------------------------------------------
valid = pat & ~isnan(Cold) & ~isnan(NMDAS) & ...
         ~isnan(Age) & ~isnan(T.coldpressor_meanFD);

tbl = table( ...
    Cold(valid), ...
    NMDAS(valid), ...
    Age(valid), ...
    T.coldpressor_meanFD(valid), ...
    'VariableNames', {'Y','NMDAS','Age','FD'} ...
);

mdl_cold_sev = fitlm(tbl, 'Y ~ NMDAS + Age + FD');
disp(mdl_cold_sev)

fprintf('NMDAS → Cold: t=%.3f, p=%.4f\n', ...
    mdl_cold_sev.Coefficients.tStat('NMDAS'), ...
    mdl_cold_sev.Coefficients.pValue('NMDAS') );


%% -------- N-BACK --------------------------------------------------------
valid = pat & BehOK & ~isnan(Nback) & ~isnan(NMDAS) & ...
         ~isnan(Age) & ~isnan(T.nback_meanFD);

tbl = table( ...
    Nback(valid), ...
    NMDAS(valid), ...
    Age(valid), ...
    T.nback_meanFD(valid), ...
    'VariableNames', {'Y','NMDAS','Age','FD'} ...
);

mdl_nback_sev = fitlm(tbl, 'Y ~ NMDAS + Age + FD');
disp(mdl_nback_sev)

fprintf('NMDAS → N-back: t=%.3f, p=%.4f\n', ...
    mdl_nback_sev.Coefficients.tStat('NMDAS'), ...
    mdl_nback_sev.Coefficients.pValue('NMDAS') );


% Load your table
data = T;  % replace with your table variable


% Group coding
patients = data.GroupID == 0;
controls = data.GroupID == 1;

% Task-specific inclusion rules
include_multisensory = ~isnan(data.multisensory_meanFD);
include_cold         = ~isnan(data.coldpressor_meanFD);
include_nback        = data.x_Behavioral_Available_ == 1 & ~isnan(data.nback_meanFD);

tasks = {
    'Multisensory', 'multisensory_meanFD', 'multisensory_meanDVARS', include_multisensory;
    'Cold Pressor', 'coldpressor_meanFD',   'coldpressor_meanDVARS', include_cold;
    'N-back',       'nback_meanFD',         'nback_meanDVARS',        include_nback
};

fprintf('\n=============== QC GROUP COMPARISONS (TASK-SPECIFIC SAMPLES) ===============\n');
LineWidth=2;
Fontsize=12;
dotcolor2 = [0, 114, 189] / 255;    % Control - Blue
dotcolor1 = [217, 83, 25] / 255;    % Patient - Orange/Red
colorcoding=[dotcolor1; dotcolor2];
pointsize=4;


% --- Prepare figure ---
nTasks = size(tasks,1);
figure('Position',[100 100 1200 300*nTasks]);

for i = 1:nTasks

    taskName   = tasks{i,1};
    FD_var     = tasks{i,2};
    DVARS_var  = tasks{i,3};
    include    = tasks{i,4};

    fprintf('\n--- %s ---\n', taskName);

    % Apply inclusion
    pat_idx = patients & include;
    con_idx = controls & include;

    % Extract QC measures
    FD_pat = data.(FD_var)(pat_idx);
    FD_con = data.(FD_var)(con_idx);

    DV_pat = data.(DVARS_var)(pat_idx);
    DV_con = data.(DVARS_var)(con_idx);

    % Compare FD
    [~,p_fd,~,stats_fd] = ttest2(FD_pat, FD_con);
    fprintf('FD: t(%d) = %.3f, p = %.4f\n', stats_fd.df, stats_fd.tstat, p_fd);

    % Compare DVARS
    [~,p_dv,~,stats_dv] = ttest2(DV_pat, DV_con);
    fprintf('DVARS: t(%d) = %.3f, p = %.4f\n', stats_dv.df, stats_dv.tstat, p_dv);


    %% ============ FD subplot (column 1) ============
    subplot(nTasks, 2, (i-1)*2 + 1)
    vp = violinplot({FD_pat, FD_con}, {'Patients','Controls'}, ...
                    'facecolor', colorcoding, ...
                    'mc', 'k', ...
                    'plotlegend', 0, ...
                    'pointsize', pointsize);
    ylabel('Mean FD');
    title(sprintf('%s – FD', taskName));


    %% ============ DVARS subplot (column 2) ============
    subplot(nTasks, 2, (i-1)*2 + 2)
    vp = violinplot({DV_pat, DV_con}, {'Patients','Controls'}, ...
                    'facecolor', colorcoding, ...
                    'mc', 'k', ...
                    'plotlegend', 0, ...
                    'pointsize', pointsize);
    ylabel('Mean DVARS');
    title(sprintf('%s – DVARS', taskName));

end


%% ========================================================================
% GROUP x PERFORMANCE INTERACTION (Reviewer 1 Q5)
% ========================================================================
% Tests whether group differences in brain activation are driven by
% differences in task performance or altered brain-behavior coupling

fprintf('\n========== GROUP x PERFORMANCE INTERACTION ==========\n');

% Valid subjects for this analysis
valid_gxp = BehOK & ~isnan(Nback) & ~isnan(BehAcc) & ~isnan(Age);

% Zscore predictors for interpretability
Y     = Nback(valid_gxp);
G     = Group(valid_gxp);
Acc   = BehAcc(valid_gxp);
A     = Age(valid_gxp);

% Build table
tbl_gxp = table(Y, G, Acc, A, ...
    'VariableNames', {'Brain','Group','Accuracy','Age'});

% Model 1: Main effects only (Group + Accuracy + Age)
mdl1 = fitlm(tbl_gxp, 'Brain ~ Group + Accuracy + Age');
fprintf('\nModel 1: Main effects (Group + Accuracy + Age)\n');
fprintf('Group effect:    t(%.0f)=%.3f, p=%.4f\n', ...
    mdl1.Coefficients.tStat('Group'), ...
    mdl1.Coefficients.tStat('Group'), ...
    mdl1.Coefficients.pValue('Group'));
fprintf('Accuracy effect: t(%.0f)=%.3f, p=%.4f\n', ...
    mdl1.Coefficients.tStat('Accuracy'), ...
    mdl1.Coefficients.tStat('Accuracy'), ...
    mdl1.Coefficients.pValue('Accuracy'));

% Model 2: With Group x Accuracy interaction
mdl2 = fitlm(tbl_gxp, 'Brain ~ Group + Accuracy + Age + Group:Accuracy');
fprintf('\nModel 2: With Group x Accuracy interaction\n');
fprintf('Group effect:             t(%.0f)=%.3f, p=%.4f\n', ...
    mdl2.Coefficients.tStat('Group'), ...
    mdl2.Coefficients.tStat('Group'), ...
    mdl2.Coefficients.pValue('Group'));
fprintf('Accuracy effect:          t(%.0f)=%.3f, p=%.4f\n', ...
    mdl2.Coefficients.tStat('Accuracy'), ...
    mdl2.Coefficients.tStat('Accuracy'), ...
    mdl2.Coefficients.pValue('Accuracy'));
fprintf('Group x Accuracy interaction: t(%.0f)=%.3f, p=%.4f\n', ...
    mdl2.Coefficients.tStat('Group:Accuracy'), ...
    mdl2.Coefficients.tStat('Group:Accuracy'), ...
    mdl2.Coefficients.pValue('Group:Accuracy'));

% Model comparison
fprintf('\nModel comparison (main effects vs interaction):\n');
% With this:
fprintf('\nModel comparison (main effects vs interaction):\n');
fprintf('Model 1 R-squared: %.4f\n', mdl1.Rsquared.Ordinary);
fprintf('Model 2 R-squared: %.4f\n', mdl2.Rsquared.Ordinary);
fprintf('Model 1 AIC: %.4f\n', mdl1.ModelCriterion.AIC);
fprintf('Model 2 AIC: %.4f\n', mdl2.ModelCriterion.AIC);

% F-test for model comparison
rss1 = sum(mdl1.Residuals.Raw.^2);
rss2 = sum(mdl2.Residuals.Raw.^2);
df1  = mdl1.DFE;
df2  = mdl2.DFE;
F    = ((rss1 - rss2) / (df1 - df2)) / (rss2 / df2);
p_ftest = 1 - fcdf(F, df1 - df2, df2);
fprintf('\nF-test comparing Model 1 vs Model 2:\n');
fprintf('F(%.0f, %.0f) = %.3f, p = %.4f\n', df1-df2, df2, F, p_ftest);
fprintf('(Non-significant p confirms interaction term does not improve model fit)\n');

fprintf('\nN subjects in this analysis: %d\n', sum(valid_gxp));
fprintf('N patients: %d, N controls: %d\n', ...
    sum(valid_gxp & Group==0), sum(valid_gxp & Group==1));
fprintf('\n=============================================================================\n');


function print_ttest(x, y, name)
    [~, p, ~, stats] = ttest2(x, y);
    [~,p,ci,stats] = ttest2(x, y);
    fprintf('%s: t(%d)=%.3f, p=%.4f\n', name, stats.df, stats.tstat, p);
end

function q = bh_fdr_fn(p)
% Benjamini-Hochberg FDR correction.
% Input:  p — column vector of raw p-values
% Output: q — BH-adjusted q-values (same order as input)
%
% Procedure:
%   1. Rank p-values from smallest to largest
%   2. q(i) = p(i) * n / rank(i)
%   3. Enforce monotonicity from the largest p downward
%   4. Cap at 1
    p = p(:);
    n = numel(p);
    [p_sorted, sort_idx] = sort(p);
    q_sorted = min(1, p_sorted .* n ./ (1:n)');
    % Enforce monotonicity: q at rank i <= q at rank i+1
    for i = n-1:-1:1
        q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
    end
    % Map back to original order
    q = nan(n, 1);
    q(sort_idx) = q_sorted;
end
