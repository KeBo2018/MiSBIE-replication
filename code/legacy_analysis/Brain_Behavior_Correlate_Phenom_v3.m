% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
%% ========================================================================
%  CORRECTED: Phenome-wide correlation of Brain & Behavior with all variables
%  Fully aligned with SimpleCodeForMajorAnlysis.m
%
%  NO .mat files needed for main analysis — everything derived from:
%    1. DataMatrix_Behavior.xlsx  (brain, behavior, group, age, NMDAS, etc.)
%    2. Ke_DataforCorrelationMatrix_Cleaned_CK_Combined_SingleItem.xlsx (phenome)
%
%  KEY FIXES vs original:
%   1. Uses partialcorr with Age as covariate (matching SimpleCode)
%   2. All variables sourced from DataMatrix_Behavior.xlsx (same as SimpleCode)
%   3. No dependency on BrainAndBehavior_74.mat or MitoID_Nback.mat
%   4. Explicit NaN handling before correlation
% =========================================================================

clear; clc; close all;

%% =======================
% SET PATHS (update for your machine)
% ========================
mainpath = 'LOCAL_SOURCE\Dartmouth College Dropbox\Ke Bo\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\Manuscript\Submission2_Science\Manuscript revision\Analysis';

%% =======================
% LOAD DATA — identical to SimpleCodeForMajorAnlysis.m
% ========================

T_main = readtable(fullfile(mainpath, 'DataMatrix_Behavior.xlsx'));

% Extract variables (identical variable names as SimpleCode)
Multi  = T_main.Multisensory;
Cold   = T_main.Cold;
Nback  = T_main.Nback;
Age    = T_main.Age;
Group  = T_main.GroupID;                  % 0 = patient, 1 = control
NMDAS  = T_main.NMDAS;
BehAcc = T_main.Nback_ACC_;
BehOK  = T_main.x_Behavioral_Available_ == 1;

%% =======================
% DERIVE THE 74 N-BACK SUBJECTS (replaces .mat loading)
% ========================

% BehOK == 1 identifies the 74 subjects with N-back data
% (verified: identical to subjects in MitoID_Nback.mat)
idx_74 = find(BehOK);          % row indices into the 110-subject table
n_74   = length(idx_74);       % should be 74

% Brain activation & behavioral accuracy for the 74
Brain_74    = Nback(idx_74);       % same as Brain from BrainAndBehavior_74.mat
Behavior_74 = BehAcc(idx_74);     % same as Behavior from BrainAndBehavior_74.mat
Age_74      = Age(idx_74);
Group_74    = Group(idx_74);       % same as X1_S from MitoID_Nback.mat

% Patient mask within the 74
pat_mask = Group_74 == 0;

fprintf('N-back subjects: %d (patients: %d, controls: %d)\n', ...
    n_74, sum(pat_mask), sum(~pat_mask));

%% =======================
% LOAD PHENOME TABLE & SUBSET
% ========================

inputFile = 'LOCAL_SOURCE\Dartmouth College Dropbox\Ke Bo\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\Ke_DataforCorrelationMatrix_Cleaned_CK_Combined_SingleItem.xlsx';
T_phenom = readtable(inputFile);

% Subset phenome table to the 74 N-back subjects
% idx_74 are row indices into the 110-row table; phenome table also has
% 110 rows (WITHHELD_SUBJECT-WITHHELD_SUBJECT), so indexing is direct.
T_S = T_phenom(idx_74, :);

% Further subset to patients only
T_S_Patient      = T_S(pat_mask, :);
Brain_Patient    = Brain_74(pat_mask);
Behavior_Patient = Behavior_74(pat_mask);
Age_Patient      = Age_74(pat_mask);
NMDAS_pat        = NMDAS(idx_74(pat_mask));   % from DataMatrix (same source as SimpleCode)

fprintf('Patient subsample: n=%d\n', sum(pat_mask));


%% ========================================================================
% SANITY CHECK: Reproduce SimpleCode NMDAS vs N-back Brain result
% ========================================================================

fprintf('\n========== SANITY CHECK ==========\n');

valid_check = ~isnan(Brain_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient);
[r_check, p_check] = partialcorr(Brain_Patient(valid_check), NMDAS_pat(valid_check), ...
                                   Age_Patient(valid_check), ...
                                   'Type', 'Spearman', 'Rows', 'complete');
fprintf('partialcorr(Brain, NMDAS | Age): r=%.3f, p=%.4f, n=%d\n', ...
    r_check, p_check, sum(valid_check));
fprintf('(Should match SimpleCode Section 2 valid3 result exactly)\n');


%% ========================================================================
% PHENOME-WIDE CORRELATIONS: Brain ~ Phenome (patients, age-corrected)
% ========================================================================

fprintf('\n========== BRAIN ~ PHENOME (age-corrected Spearman partial corr) ==========\n');

nVars = size(T_S_Patient, 2);
Rvalue_Brain = nan(nVars, 1);
Pvalue_Brain = nan(nVars, 1);
Nobs_Brain   = nan(nVars, 1);

for CNum = 2:nVars
    currentData = T_S_Patient.(CNum);
    
    if isnumeric(currentData) || islogical(currentData)
        currentData = double(currentData);
        
        validRows = ~isnan(currentData) & ~isnan(Brain_Patient) & ~isnan(Age_Patient);
        
        if sum(validRows) > 3
            [R, P] = partialcorr(currentData(validRows), Brain_Patient(validRows), ...
                                 Age_Patient(validRows), ...
                                 'Type', 'Spearman', 'Rows', 'complete');
            Rvalue_Brain(CNum) = R;
            Pvalue_Brain(CNum) = P;
            Nobs_Brain(CNum)   = sum(validRows);
        end
    end
end

% FDR correction (only on non-NaN p-values)
FDR_Brain = ones(nVars, 1);
valid_p = ~isnan(Pvalue_Brain);
FDR_Brain(valid_p) = fdr(Pvalue_Brain(valid_p));

varNames = T_S_Patient.Properties.VariableNames';

outFile_Brain = fullfile(mainpath, 'MISBIE_Brain_Correlation_FDR_AgeCorrected.xlsx');
save_corr_results_excel(varNames, Rvalue_Brain, Pvalue_Brain, FDR_Brain, Nobs_Brain, outFile_Brain);


%% ========================================================================
% PHENOME-WIDE CORRELATIONS: Behavior ~ Phenome (patients, age-corrected)
% ========================================================================

fprintf('\n========== BEHAVIOR ~ PHENOME (age-corrected Spearman partial corr) ==========\n');

Rvalue_Beh = nan(nVars, 1);
Pvalue_Beh = nan(nVars, 1);
Nobs_Beh   = nan(nVars, 1);

for CNum = 2:nVars
    currentData = T_S_Patient.(CNum);
    
    if isnumeric(currentData) || islogical(currentData)
        currentData = double(currentData);
        
        validRows = ~isnan(currentData) & ~isnan(Behavior_Patient) & ~isnan(Age_Patient);
        
        if sum(validRows) > 3
            [R, P] = partialcorr(currentData(validRows), Behavior_Patient(validRows), ...
                                 Age_Patient(validRows), ...
                                 'Type', 'Spearman', 'Rows', 'complete');
            Rvalue_Beh(CNum) = R;
            Pvalue_Beh(CNum) = P;
            Nobs_Beh(CNum)   = sum(validRows);
        end
    end
end

% FDR correction
FDR_Beh = ones(nVars, 1);
valid_p = ~isnan(Pvalue_Beh);
FDR_Beh(valid_p) = fdr(Pvalue_Beh(valid_p));

outFile_Behavior = fullfile(mainpath, 'MISBIE_Behavior_Correlation_FDR_AgeCorrected.xlsx');
save_corr_results_excel(varNames, Rvalue_Beh, Pvalue_Beh, FDR_Beh, Nobs_Beh, outFile_Behavior);


%% ========================================================================
% ADDITIONAL CHECKS: Vision and Duration
% ========================================================================

Vision = table2array(T_S_Patient(:, 175));

[R, P] = partialcorr(Brain_Patient, Vision, Age_Patient, ...
                      'Type', 'Spearman', 'Rows', 'complete');
fprintf('\nBrain vs Vision (age-corrected): R=%.4f, P=%.4f\n', R, P);

[R, P] = partialcorr(NMDAS_pat, Brain_Patient, Age_Patient, ...
                      'Type', 'Spearman', 'Rows', 'complete');
fprintf('NMDAS vs Brain (age-corrected):  R=%.4f, P=%.4f\n', R, P);

[R, P] = partialcorr(NMDAS_pat, Behavior_Patient, [Vision, Age_Patient], ...
                      'Type', 'Spearman', 'Rows', 'complete');
fprintf('NMDAS vs Behavior | Vision+Age:  R=%.4f, P=%.4f\n', R, P);


%% ========================================================================
% SPLIT-HALF ROBUSTNESS: NMDAS vs Brain & Behavior (patients only)
% ========================================================================

nIter = 1000;
rng(42);

valid_brain = find(~isnan(Brain_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient));
valid_beh   = find(~isnan(Behavior_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient));

n_brain = length(valid_brain);
n_beh   = length(valid_beh);

rSplit_brain = nan(nIter, 1);
rSplit_beh   = nan(nIter, 1);

for i = 1:nIter
    idx_shuf = valid_brain(randperm(n_brain));
    half1 = idx_shuf(1:floor(n_brain/2));
    half2 = idx_shuf(floor(n_brain/2)+1:end);
    
    r1 = partialcorr(Brain_Patient(half1), NMDAS_pat(half1), ...
                     Age_Patient(half1), 'Type', 'Spearman', 'Rows', 'complete');
    r2 = partialcorr(Brain_Patient(half2), NMDAS_pat(half2), ...
                     Age_Patient(half2), 'Type', 'Spearman', 'Rows', 'complete');
    rSplit_brain(i) = mean([r1, r2]);
end

for i = 1:nIter
    idx_shuf = valid_beh(randperm(n_beh));
    half1 = idx_shuf(1:floor(n_beh/2));
    half2 = idx_shuf(floor(n_beh/2)+1:end);
    
    r1 = partialcorr(Behavior_Patient(half1), NMDAS_pat(half1), ...
                     Age_Patient(half1), 'Type', 'Spearman', 'Rows', 'complete');
    r2 = partialcorr(Behavior_Patient(half2), NMDAS_pat(half2), ...
                     Age_Patient(half2), 'Type', 'Spearman', 'Rows', 'complete');
    rSplit_beh(i) = mean([r1, r2]);
end

fprintf('\n========== SPLIT-HALF ROBUSTNESS (age-controlled) ==========\n');

fprintf('\nNMDAS vs Brain Activation (N-back):\n');
fprintf('  Full-sample partial r = %.3f\n', r_check);
fprintf('  Median split-half r   = %.3f\n', median(rSplit_brain));
fprintf('  95%% CI               = [%.3f, %.3f]\n', ...
    prctile(rSplit_brain, 2.5), prctile(rSplit_brain, 97.5));
fprintf('  %% same direction     = %.1f%%\n', 100 * mean(rSplit_brain < 0));

[r_beh_full, p_beh_full] = partialcorr(Behavior_Patient(valid_beh), ...
    NMDAS_pat(valid_beh), Age_Patient(valid_beh), ...
    'Type', 'Spearman', 'Rows', 'complete');

fprintf('\nNMDAS vs Behavioral Performance:\n');
fprintf('  Full-sample partial r = %.3f (p=%.4f)\n', r_beh_full, p_beh_full);
fprintf('  Median split-half r   = %.3f\n', median(rSplit_beh));
fprintf('  95%% CI               = [%.3f, %.3f]\n', ...
    prctile(rSplit_beh, 2.5), prctile(rSplit_beh, 97.5));
fprintf('  %% same direction     = %.1f%%\n', 100 * mean(rSplit_beh < 0));

%% --- Plot split-half distributions ---
figure('Color', 'w', 'Position', [200 200 800 350]);

subplot(1,2,1)
histogram(rSplit_brain, 30, 'FaceColor', [217 83 25]/255, 'EdgeColor', 'w')
hold on
xline(median(rSplit_brain), 'k--', 'LineWidth', 2)
xline(prctile(rSplit_brain, 2.5),  'k:', 'LineWidth', 1.5)
xline(prctile(rSplit_brain, 97.5), 'k:', 'LineWidth', 1.5)
xline(0, 'r-', 'LineWidth', 1)
xlabel('Mean split-half Spearman r')
ylabel('Count')
title('NMDAS vs Brain Activation')
legend({'Distribution','Median','95% CI','Null (r=0)'}, 'Location', 'northwest')
set(gca, 'FontSize', 11)

subplot(1,2,2)
histogram(rSplit_beh, 30, 'FaceColor', [0 114 189]/255, 'EdgeColor', 'w')
hold on
xline(median(rSplit_beh), 'k--', 'LineWidth', 2)
xline(prctile(rSplit_beh, 2.5),  'k:', 'LineWidth', 1.5)
xline(prctile(rSplit_beh, 97.5), 'k:', 'LineWidth', 1.5)
xline(0, 'r-', 'LineWidth', 1)
xlabel('Mean split-half Spearman r')
ylabel('Count')
title('NMDAS vs Behavioral Performance')
set(gca, 'FontSize', 11)

sgtitle('Split-Half Robustness: NMDAS Correlations (1000 iterations)', 'FontSize', 13)
saveas(gcf, 'SplitHalf_Robustness_NMDAS.png')


%% ========================================================================
% WHITE MATTER AND NODDI STRUCTURAL ANALYSES
% ========================================================================
% NOTE: These sections still need .mat files for structural data
% (WhiteMatter_86Sub.mat, NODDI_86Sub.mat) — those contain DTI/NODDI
% metrics not available in DataMatrix_Behavior.xlsx.
% Subject matching uses MiID strings derived from idx_74.

% Build MitoID strings for the 25 patients (replaces Nback_S_MitoID_Patients)
patient_row_ids = idx_74(pat_mask);   % row numbers in the 110-subject table
Nback_S_MitoID_Patients = cell(sum(pat_mask), 1);
for i = 1:sum(pat_mask)
    Nback_S_MitoID_Patients{i} = sprintf('Mi%03d', patient_row_ids(i));
end

last_three_chars = cell(size(Nback_S_MitoID_Patients));
for i = 1:length(Nback_S_MitoID_Patients)
    last_three_chars{i} = Nback_S_MitoID_Patients{i}(end-2:end);
end

%% --- WHITE MATTER ---
load(fullfile(mainpath, 'WhiteMatter_86Sub.mat'))

[common_ids, idx_subject_ids, idx_MitoID_Patients] = intersect(subject_ids, last_three_chars);

NMDAS_Common  = NMDAS_pat(idx_MitoID_Patients);
Brain_common  = Brain_Patient(idx_MitoID_Patients);
White_common1 = ValsCT(idx_subject_ids, 1);
White_common2 = ValsFA(idx_subject_ids, 2);

fprintf('\n========== WHITE MATTER STRUCTURAL ANALYSES ==========\n');
fprintf('N subjects: %d\n', length(NMDAS_Common));

[R,P] = corr(NMDAS_Common, White_common1, 'type', 'Spearman');
fprintf('NMDAS vs Cortical Thickness:            R=%.4f, P=%.4f\n', R, P);

[R,P] = corr(NMDAS_Common, White_common2, 'type', 'Spearman');
fprintf('NMDAS vs FA:                            R=%.4f, P=%.4f\n', R, P);

[R,P] = corr(Brain_common, White_common1, 'type', 'Spearman');
fprintf('Brain Activation vs Cortical Thickness: R=%.4f, P=%.4f\n', R, P);

[R,P] = corr(Brain_common, White_common2, 'type', 'Spearman');
fprintf('Brain Activation vs FA:                 R=%.4f, P=%.4f\n', R, P);

[R,P] = partialcorr(Brain_common, NMDAS_Common, ...
                    [White_common1 White_common2], 'type', 'Spearman');
fprintf('Brain vs NMDAS controlling for CT+FA:   R=%.4f, P=%.4f\n', R, P);

figure; scatter(Brain_common, White_common1, 200, '.'); h=lsline;
xlabel('Working memory brain pattern response'); ylabel('Average cortical thickness');
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'LineWidth', 1); set(h(1), 'color', '#FF594C', 'linewidth', 2)

figure; scatter(Brain_common, White_common2, 200, '.'); h=lsline;
xlabel('Working memory brain pattern response'); ylabel('Average FA');
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'LineWidth', 1); set(h(1), 'color', '#FF594C', 'linewidth', 2)

figure; scatter(NMDAS_Common, White_common1, 200, '.'); h=lsline;
xlabel('NMDAS'); ylabel('Average cortical thickness');
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'LineWidth', 1); set(h(1), 'color', '#FF594C', 'linewidth', 2)

figure; scatter(NMDAS_Common, White_common2, 200, '.'); h=lsline;
xlabel('NMDAS'); ylabel('Average FA');
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'LineWidth', 1); set(h(1), 'color', '#FF594C', 'linewidth', 2)


%% --- NODDI ---
load(fullfile(mainpath, 'NODDI_86Sub.mat'))

[common_ids, idx_subject_ids, idx_MitoID_Patients] = intersect(subject_ids, last_three_chars);

NMDAS_Common  = NMDAS_pat(idx_MitoID_Patients);
Brain_common  = Brain_Patient(idx_MitoID_Patients);
White_common1 = ValsICVF(idx_subject_ids, 2);
White_common2 = ValsMax(idx_subject_ids, 3);

fprintf('\n========== NODDI STRUCTURAL ANALYSES ==========\n');
fprintf('N subjects: %d\n', length(NMDAS_Common));

[R,P] = corr(NMDAS_Common, White_common1, 'type', 'Spearman');
fprintf('NMDAS vs ICVF:                          R=%.4f, P=%.4f\n', R, P);

[R,P] = corr(NMDAS_Common, White_common2, 'type', 'Spearman');
fprintf('NMDAS vs Max:                           R=%.4f, P=%.4f\n', R, P);

[R,P] = corr(Brain_common, White_common1, 'type', 'Spearman');
fprintf('Brain Activation vs ICVF:               R=%.4f, P=%.4f\n', R, P);

[R,P] = corr(Brain_common, White_common2, 'type', 'Spearman');
fprintf('Brain Activation vs Max:                R=%.4f, P=%.4f\n', R, P);

[R,P] = partialcorr(Brain_common, NMDAS_Common, ...
                    [White_common1 White_common2], 'type', 'Spearman');
fprintf('Brain vs NMDAS controlling for ICVF+Max: R=%.4f, P=%.4f\n', R, P);

figure; scatter(Brain_common, White_common1, 200, '.'); h=lsline;
xlabel('Working memory brain pattern response'); ylabel('ICVF');
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'LineWidth', 1); set(h(1), 'color', '#FF594C', 'linewidth', 2)

figure; scatter(Brain_common, White_common2, 200, '.'); h=lsline;
xlabel('Working memory brain pattern response'); ylabel('Max (NODDI)');
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'LineWidth', 1); set(h(1), 'color', '#FF594C', 'linewidth', 2)

figure; scatter(NMDAS_Common, White_common1, 200, '.'); h=lsline;
xlabel('NMDAS'); ylabel('ICVF');
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'LineWidth', 1); set(h(1), 'color', '#FF594C', 'linewidth', 2)

figure; scatter(NMDAS_Common, White_common2, 200, '.'); h=lsline;
xlabel('NMDAS'); ylabel('Max (NODDI)');
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'LineWidth', 1); set(h(1), 'color', '#FF594C', 'linewidth', 2)


%% ========================================================================
% SENSITIVITY ANALYSES: Partial out potential confounders
% ========================================================================
% Addresses Reviewer 2 Q2 & Reviewer 3 Q3:
%   "How can we be sure that the imaging or cognitive findings do not
%    simply result from differences in visual ability, basic neurological
%    function, or levels of fatigue?"
%
% Strategy: For each potential confound, repeat the two key analyses:
%   (A) NMDAS vs Brain  (partialcorr, Spearman, controlling Age + confound)
%   (B) NMDAS vs Behavior (partialcorr, Spearman, controlling Age + confound)
% Compare R and p with the baseline (Age-only) model.
%
% Also test direct association of each confound with Brain and Behavior
% to see which confounds are themselves related to the outcomes.
%
% Sources:
%   - 5 variables available in DataMatrix_Behavior.xlsx
%   - 13 additional variables from the phenome table (T_S_Patient)
% ========================================================================

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  SENSITIVITY ANALYSES: CONFOUND-BY-CONFOUND PARTIAL CORR\n');
fprintf('================================================================\n');

%% --- Define confound variables ---
% Each row: { 'Label', source, column_spec }
%   source = 'main'   -> from DataMatrix_Behavior.xlsx (T_main, indexed by idx_74 & pat_mask)
%   source = 'phenom' -> from phenome table (T_S_Patient, already patient-subset)
%   column_spec: variable name (for 'main') or column index (for 'phenom')

confounds = {
    % ---- DEMOGRAPHIC ----
    % Sex: binary variable. Column name 'sex' from DataMatrix_Behavior.xlsx.
    % safe_to_double() handles string ('M'/'F') or numeric coding automatically.
    'Sex',                                           'main', 'sex';
    % Genotype: nominal 3-class — split into dummy variables (one vs rest).
    % Avoids imposing false ordinal structure. Dummies computed below from
    % 'geneticDiagnosis' column and stored in a Dummy_Vars struct.
    'Genotype: Point Mutation (vs others)',          'dummy', 'is_PointMutation';
    'Genotype: Single Deletion (vs others)',         'dummy', 'is_Deletion';
    'Genotype: MELAS (vs others)',                   'dummy', 'is_MELAS';
    % ---- VISION ----
    'nmdas_cf_1 (Vision self-report)',          'phenom', 377;
    'nmdas_cca_1 (Visual acuity clinical)',     'phenom', 396;
    'nmdas_cca_2 (Ptosis)',                     'phenom', 397;
    'nmdas_cca_3 (CPEO)',                       'phenom', 398;
    'cns_eyes (CNS Eyes exam)',                 'phenom', 210;
    'cns_crnl_nrv_ii (Cranial nerve II)',       'phenom', 248;
    'cns_crnl_nrv_iii (Cranial nerve III)',     'phenom', 250;
    'cns_crnl_nrv_iv (Cranial nerve IV)',       'phenom', 252;
    'cns_crnl_nrv_vi (Cranial nerve VI)',       'phenom', 256;
    % ---- FATIGUE ----
    'pfs_phys (Physical fatigue)',              'phenom',  48;
    'pfs_mental (Mental fatigue)',              'phenom',  49;
    % ---- NEUROLOGICAL / MOTOR ----
    'nmdas_cf_9 (Exercise tolerance)',          'phenom', 385;
    'nmdas_cf_10 (Gait stability)',             'phenom', 386;
    'nmdas_cca_5 (Myopathy)',                   'phenom', 400;
    'nmdas_cca_6 (Cerebellar ataxia)',          'phenom', 401;
    'nmdas_cca_7 (Neuropathy)',                 'phenom', 402;
    'cns_gait_statn (Gait and station)',        'phenom', 270;
    % ---- PSYCHIATRIC ----
    'nmdas_ssi_1 (Psychiatric symptoms)',       'phenom', 387;
    % ---- VISION ----
    'nmdas_cf_1 (Vision self-report)',          'phenom', 377;
    'nmdas_cca_1 (Visual acuity clinical)',     'phenom', 396;
    'nmdas_cca_2 (Ptosis)',                     'phenom', 397;
    'nmdas_cca_3 (CPEO)',                       'phenom', 398;
    'cns_eyes (CNS Eyes exam)',                 'phenom', 210;
    'cns_crnl_nrv_ii (Cranial nerve II)',       'phenom', 248;
    'cns_crnl_nrv_iii (Cranial nerve III)',     'phenom', 250;
    'cns_crnl_nrv_iv (Cranial nerve IV)',       'phenom', 252;
    'cns_crnl_nrv_vi (Cranial nerve VI)',       'phenom', 256;
    % ---- FATIGUE ----
    'pfs_phys (Physical fatigue)',              'phenom',  48;
    'pfs_mental (Mental fatigue)',              'phenom',  49;
    % ---- NEUROLOGICAL / MOTOR ----
    'nmdas_cf_9 (Exercise tolerance)',          'phenom', 385;
    'nmdas_cf_10 (Gait stability)',             'phenom', 386;
    'nmdas_cca_5 (Myopathy)',                   'phenom', 400;
    'nmdas_cca_6 (Cerebellar ataxia)',          'phenom', 401;
    'nmdas_cca_7 (Neuropathy)',                 'phenom', 402;
    'cns_gait_statn (Gait and station)',        'phenom', 270;
    % ---- PSYCHIATRIC ----
    'nmdas_ssi_1 (Psychiatric symptoms)',       'phenom', 387;
};

nConf = size(confounds, 1);

%% --- Build genotype dummy variables (patients only) ---
% Extract geneticDiagnosis for the 25 patients from DataMatrix_Behavior.xlsx
% Adjust the column name below if it differs in your spreadsheet.
genoDiag_pat = T_main.geneticDiagnosis(idx_74(pat_mask));
genoDiag_pat = safe_to_string(genoDiag_pat);   % ensure cell array of strings

fprintf('\n[Genotype coding] Unique values in geneticDiagnosis (patients):\n');
disp(unique(genoDiag_pat));

% Binary dummy variables (1 = that class, 0 = other patient classes)
% *** Update the string literals below to exactly match your data values ***
Dummy_Vars.is_PointMutation = double(strcmpi(genoDiag_pat, 'PointMutation') | ...
                                     strcmpi(genoDiag_pat, 'Point Mutation') | ...
                                     strcmpi(genoDiag_pat, 'point_mutation'));
Dummy_Vars.is_Deletion      = double(strcmpi(genoDiag_pat, 'Deletion') | ...
                                     strcmpi(genoDiag_pat, 'Single Deletion') | ...
                                     strcmpi(genoDiag_pat, 'single_deletion'));
Dummy_Vars.is_MELAS         = double(strcmpi(genoDiag_pat, 'MELAS'));

fprintf('  Point Mutation n=%d, Deletion n=%d, MELAS n=%d (of %d patients)\n', ...
    sum(Dummy_Vars.is_PointMutation), sum(Dummy_Vars.is_Deletion), ...
    sum(Dummy_Vars.is_MELAS), sum(pat_mask));

% Also build all-74 versions for group comparison and GDF15 sections
genoDiag_74 = T_main.geneticDiagnosis(idx_74);
genoDiag_74 = safe_to_string(genoDiag_74);
Dummy_Vars_74.is_PointMutation = double(strcmpi(genoDiag_74, 'PointMutation') | ...
                                        strcmpi(genoDiag_74, 'Point Mutation') | ...
                                        strcmpi(genoDiag_74, 'point_mutation'));
Dummy_Vars_74.is_Deletion      = double(strcmpi(genoDiag_74, 'Deletion') | ...
                                        strcmpi(genoDiag_74, 'Single Deletion') | ...
                                        strcmpi(genoDiag_74, 'single_deletion'));
Dummy_Vars_74.is_MELAS         = double(strcmpi(genoDiag_74, 'MELAS'));

%% --- Baseline results (Age only) for comparison ---
valid_base_brain = ~isnan(Brain_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient);
[r_base_brain, p_base_brain] = partialcorr( ...
    Brain_Patient(valid_base_brain), NMDAS_pat(valid_base_brain), ...
    Age_Patient(valid_base_brain), 'Type', 'Spearman', 'Rows', 'complete');

valid_base_beh = ~isnan(Behavior_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient);
[r_base_beh, p_base_beh] = partialcorr( ...
    Behavior_Patient(valid_base_beh), NMDAS_pat(valid_base_beh), ...
    Age_Patient(valid_base_beh), 'Type', 'Spearman', 'Rows', 'complete');

fprintf('\nBASELINE (controlling Age only):\n');
fprintf('  NMDAS vs Brain:    r=%.3f, p=%.4f, n=%d\n', r_base_brain, p_base_brain, sum(valid_base_brain));
fprintf('  NMDAS vs Behavior: r=%.3f, p=%.4f, n=%d\n', r_base_beh, p_base_beh, sum(valid_base_beh));

%% --- Preallocate results table ---
Conf_Label     = cell(nConf, 1);
% Direct association: confound vs Brain / Behavior (controlling Age)
R_conf_brain   = nan(nConf, 1);
P_conf_brain   = nan(nConf, 1);
R_conf_beh     = nan(nConf, 1);
P_conf_beh     = nan(nConf, 1);
% Main result after controlling confound: NMDAS vs Brain / Behavior (controlling Age + confound)
R_adj_brain    = nan(nConf, 1);
P_adj_brain    = nan(nConf, 1);
N_adj_brain    = nan(nConf, 1);
R_adj_beh      = nan(nConf, 1);
P_adj_beh      = nan(nConf, 1);
N_adj_beh      = nan(nConf, 1);

%% --- Loop through confounds ---
fprintf('\n%-45s | %22s | %22s | %22s | %22s\n', ...
    'Confound', 'Conf~Brain (r, p)', 'NMDAS~Brain|Conf (r,p)', ...
    'Conf~Behav (r, p)', 'NMDAS~Behav|Conf (r,p)');
fprintf('%s\n', repmat('-', 1, 145));

for c = 1:nConf
    label  = confounds{c, 1};
    source = confounds{c, 2};
    colspec = confounds{c, 3};
    
    Conf_Label{c} = label;
    
    % --- Extract confound variable ---
    if strcmp(source, 'phenom')
        conf_var = double(table2array(T_S_Patient(:, colspec)));
    elseif strcmp(source, 'dummy')
        conf_var = Dummy_Vars.(colspec);
    else
        % 'main' source — extract from T_main using idx_74 and pat_mask
        raw = T_main.(colspec)(idx_74(pat_mask));
        conf_var = safe_to_double(raw);
    end
    
    % --- (A) Direct association: confound vs Brain (controlling Age) ---
    v = ~isnan(conf_var) & ~isnan(Brain_Patient) & ~isnan(Age_Patient);
    if sum(v) > 3
        [R_conf_brain(c), P_conf_brain(c)] = partialcorr( ...
            conf_var(v), Brain_Patient(v), Age_Patient(v), ...
            'Type', 'Spearman', 'Rows', 'complete');
    end
    
    % --- (B) Direct association: confound vs Behavior (controlling Age) ---
    v = ~isnan(conf_var) & ~isnan(Behavior_Patient) & ~isnan(Age_Patient);
    if sum(v) > 3
        [R_conf_beh(c), P_conf_beh(c)] = partialcorr( ...
            conf_var(v), Behavior_Patient(v), Age_Patient(v), ...
            'Type', 'Spearman', 'Rows', 'complete');
    end
    
    % --- (C) NMDAS vs Brain, controlling Age + confound ---
    v = ~isnan(Brain_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient) & ~isnan(conf_var);
    if sum(v) > 4   % need df for 2 covariates
        [R_adj_brain(c), P_adj_brain(c)] = partialcorr( ...
            Brain_Patient(v), NMDAS_pat(v), [Age_Patient(v), conf_var(v)], ...
            'Type', 'Spearman', 'Rows', 'complete');
        N_adj_brain(c) = sum(v);
    end
    
    % --- (D) NMDAS vs Behavior, controlling Age + confound ---
    v = ~isnan(Behavior_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient) & ~isnan(conf_var);
    if sum(v) > 4
        [R_adj_beh(c), P_adj_beh(c)] = partialcorr( ...
            Behavior_Patient(v), NMDAS_pat(v), [Age_Patient(v), conf_var(v)], ...
            'Type', 'Spearman', 'Rows', 'complete');
        N_adj_beh(c) = sum(v);
    end
    
    % --- Print row ---
    fprintf('%-45s | r=%6.3f, p=%6.4f | r=%6.3f, p=%6.4f | r=%6.3f, p=%6.4f | r=%6.3f, p=%6.4f\n', ...
        label, ...
        R_conf_brain(c), P_conf_brain(c), ...
        R_adj_brain(c), P_adj_brain(c), ...
        R_conf_beh(c), P_conf_beh(c), ...
        R_adj_beh(c), P_adj_beh(c));
end

%% ========================================================================
% BOOTSTRAP TEST: Does controlling each confound significantly attenuate
%                 the NMDAS~Brain and NMDAS~Behavior associations?
%
% Strategy (5,000 resamples per confound):
%   For each bootstrap iteration, resample patients with replacement.
%   Compute:
%     r_base_boot  = partialcorr(Brain, NMDAS | Age)            [baseline]
%     r_adj_boot   = partialcorr(Brain, NMDAS | Age + Confound) [adjusted]
%     delta_boot   = |r_base_boot| - |r_adj_boot|
%   95% CI of delta across 5000 iterations.
%   Bootstrap p = proportion of iterations where delta <= 0.
%   Sig_Decrease = TRUE when 95% CI lower bound > 0.
%
% This matches the bootstrap columns in Supplementary Table S2.
% ========================================================================

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  BOOTSTRAP: Does each confound significantly attenuate NMDAS correlations?\n');
fprintf('  (5000 resamples, seed=42)\n');
fprintf('================================================================\n');

nBoot = 5000;
rng(42);   % reproducibility

% Preallocate bootstrap output vectors
DeltaR_brain_obs = nan(nConf, 1);
DeltaR_brain_lo  = nan(nConf, 1);
DeltaR_brain_hi  = nan(nConf, 1);
DeltaR_brain_p   = nan(nConf, 1);
DeltaR_brain_sig = false(nConf, 1);

DeltaR_beh_obs   = nan(nConf, 1);
DeltaR_beh_lo    = nan(nConf, 1);
DeltaR_beh_hi    = nan(nConf, 1);
DeltaR_beh_p     = nan(nConf, 1);
DeltaR_beh_sig   = false(nConf, 1);

for c = 1:nConf
    label   = confounds{c, 1};
    source  = confounds{c, 2};
    colspec = confounds{c, 3};

    % Extract confound (patient subset only)
    if strcmp(source, 'phenom')
        conf_var = double(table2array(T_S_Patient(:, colspec)));
    elseif strcmp(source, 'dummy')
        conf_var = Dummy_Vars.(colspec);
    else
        raw = T_main.(colspec)(idx_74(pat_mask));
        conf_var = safe_to_double(raw);
    end

    % ---- BRAIN bootstrap ----
    v_brain = ~isnan(Brain_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient) & ~isnan(conf_var);
    if sum(v_brain) > 5
        Br = Brain_Patient(v_brain);
        Nm = NMDAS_pat(v_brain);
        Ag = Age_Patient(v_brain);
        Cv = conf_var(v_brain);
        n  = sum(v_brain);

        % Observed delta
        r_base_obs = partialcorr(Br, Nm, Ag,       'Type','Spearman','Rows','complete');
        r_adj_obs  = partialcorr(Br, Nm, [Ag, Cv], 'Type','Spearman','Rows','complete');
        DeltaR_brain_obs(c) = abs(r_base_obs) - abs(r_adj_obs);

        % Bootstrap loop
        delta_boot = nan(nBoot, 1);
        for b = 1:nBoot
            idx_b  = randsample(n, n, true);
            r_b    = partialcorr(Br(idx_b), Nm(idx_b), Ag(idx_b),                'Type','Spearman','Rows','complete');
            r_ba   = partialcorr(Br(idx_b), Nm(idx_b), [Ag(idx_b), Cv(idx_b)],   'Type','Spearman','Rows','complete');
            delta_boot(b) = abs(r_b) - abs(r_ba);
        end
        DeltaR_brain_lo(c)  = prctile(delta_boot,  2.5);
        DeltaR_brain_hi(c)  = prctile(delta_boot, 97.5);
        DeltaR_brain_p(c)   = mean(delta_boot <= 0);
        DeltaR_brain_sig(c) = DeltaR_brain_lo(c) > 0;
    end

    % ---- BEHAVIOR bootstrap ----
    v_beh = ~isnan(Behavior_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient) & ~isnan(conf_var);
    if sum(v_beh) > 5
        Bh = Behavior_Patient(v_beh);
        Nm = NMDAS_pat(v_beh);
        Ag = Age_Patient(v_beh);
        Cv = conf_var(v_beh);
        n  = sum(v_beh);

        r_base_obs = partialcorr(Bh, Nm, Ag,       'Type','Spearman','Rows','complete');
        r_adj_obs  = partialcorr(Bh, Nm, [Ag, Cv], 'Type','Spearman','Rows','complete');
        DeltaR_beh_obs(c) = abs(r_base_obs) - abs(r_adj_obs);

        delta_boot = nan(nBoot, 1);
        for b = 1:nBoot
            idx_b  = randsample(n, n, true);
            r_b    = partialcorr(Bh(idx_b), Nm(idx_b), Ag(idx_b),                'Type','Spearman','Rows','complete');
            r_ba   = partialcorr(Bh(idx_b), Nm(idx_b), [Ag(idx_b), Cv(idx_b)],   'Type','Spearman','Rows','complete');
            delta_boot(b) = abs(r_b) - abs(r_ba);
        end
        DeltaR_beh_lo(c)  = prctile(delta_boot,  2.5);
        DeltaR_beh_hi(c)  = prctile(delta_boot, 97.5);
        DeltaR_beh_p(c)   = mean(delta_boot <= 0);
        DeltaR_beh_sig(c) = DeltaR_beh_lo(c) > 0;
    end

    fprintf('%-45s | Brain: dR=%.3f [%.3f,%.3f] p=%.3f Sig=%d | Beh: dR=%.3f [%.3f,%.3f] p=%.3f Sig=%d\n', ...
        label, ...
        DeltaR_brain_obs(c), DeltaR_brain_lo(c), DeltaR_brain_hi(c), DeltaR_brain_p(c), DeltaR_brain_sig(c), ...
        DeltaR_beh_obs(c),   DeltaR_beh_lo(c),   DeltaR_beh_hi(c),   DeltaR_beh_p(c),   DeltaR_beh_sig(c));
end

fprintf('\nBootstrap complete.\n');

%% --- Save intermediate sensitivity results (without bootstrap) ---
SensTable = table( ...
    Conf_Label, ...
    R_conf_brain, P_conf_brain, ...
    R_adj_brain, P_adj_brain, N_adj_brain, ...
    R_conf_beh, P_conf_beh, ...
    R_adj_beh, P_adj_beh, N_adj_beh, ...
    'VariableNames', { ...
        'Confound', ...
        'R_Confound_vs_Brain', 'P_Confound_vs_Brain', ...
        'R_NMDAS_Brain_AdjConf', 'P_NMDAS_Brain_AdjConf', 'N_Brain', ...
        'R_Confound_vs_Behavior', 'P_Confound_vs_Behavior', ...
        'R_NMDAS_Behav_AdjConf', 'P_NMDAS_Behav_AdjConf', 'N_Behav'} ...
);

outFile_Sens = fullfile(mainpath, 'MISBIE_Sensitivity_Confounds.xlsx');
writetable(SensTable, outFile_Sens);
fprintf('\nSensitivity results saved to:\n%s\n', outFile_Sens);

%% --- Summary printout ---
fprintf('\n========== SENSITIVITY SUMMARY ==========\n');
fprintf('Baseline NMDAS vs Brain (Age only):    r=%.3f, p=%.4f\n', r_base_brain, p_base_brain);
fprintf('Baseline NMDAS vs Behavior (Age only): r=%.3f, p=%.4f\n\n', r_base_beh, p_base_beh);

sig_brain = P_adj_brain < 0.05;
sig_beh   = P_adj_beh < 0.05;

fprintf('NMDAS~Brain remains significant after controlling for:\n');
for c = 1:nConf
    if sig_brain(c)
        fprintf('  [PASS] %-45s r=%.3f, p=%.4f (n=%d)\n', ...
            Conf_Label{c}, R_adj_brain(c), P_adj_brain(c), N_adj_brain(c));
    else
        fprintf('  [FAIL] %-45s r=%.3f, p=%.4f (n=%d)\n', ...
            Conf_Label{c}, R_adj_brain(c), P_adj_brain(c), N_adj_brain(c));
    end
end

fprintf('\nNMDAS~Behavior after controlling for each confound:\n');
for c = 1:nConf
    if sig_beh(c)
        fprintf('  [PASS] %-45s r=%.3f, p=%.4f (n=%d)\n', ...
            Conf_Label{c}, R_adj_beh(c), P_adj_beh(c), N_adj_beh(c));
    else
        fprintf('  [FAIL] %-45s r=%.3f, p=%.4f (n=%d)\n', ...
            Conf_Label{c}, R_adj_beh(c), P_adj_beh(c), N_adj_beh(c));
    end
end


%% ========================================================================
% SENSITIVITY: GROUP COMPARISONS controlling for each confound
% ========================================================================
% Repeats the key group comparison (N-back Brain: patients vs controls)
% from SimpleCode Section 1, adding each confound as covariate.
% Uses fitlm: Brain ~ Group + Age + Confound
% This uses ALL 74 N-back subjects (not just patients).
% ========================================================================

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  SENSITIVITY: GROUP COMPARISON (N-back) WITH CONFOUND CONTROL\n');
fprintf('================================================================\n');

% Full 74-subject data
Brain_all_74 = Brain_74;
Age_all_74   = Age_74;
Group_all_74 = Group_74;

% Baseline: Group effect without confound
v_base = ~isnan(Brain_all_74) & ~isnan(Age_all_74) & ~isnan(Group_all_74);
tbl_base = table(Brain_all_74(v_base), Group_all_74(v_base), Age_all_74(v_base), ...
    'VariableNames', {'Brain','Group','Age'});
mdl_base = fitlm(tbl_base, 'Brain ~ Group + Age');
fprintf('\nBASELINE (Group + Age):\n');
fprintf('  Group effect: t=%.3f, p=%.4f\n', ...
    mdl_base.Coefficients.tStat('Group'), mdl_base.Coefficients.pValue('Group'));

% For group comparison, confound variables come from the 74-subject subset
T_S_all74 = T_phenom(idx_74, :);    % phenome for all 74

fprintf('\n%-45s | %22s | %22s\n', 'Confound', 'Group t-stat', 'Group p-value');
fprintf('%s\n', repmat('-', 1, 95));

Grp_t = nan(nConf, 1);
Grp_p = nan(nConf, 1);

for c = 1:nConf
    label  = confounds{c, 1};
    source = confounds{c, 2};
    colspec = confounds{c, 3};
    
    % Extract confound for all 74
    if strcmp(source, 'phenom')
        conf_all74 = double(table2array(T_S_all74(:, colspec)));
    elseif strcmp(source, 'dummy')
        conf_all74 = Dummy_Vars_74.(colspec);
    else
        raw = T_main.(colspec)(idx_74);
        conf_all74 = safe_to_double(raw);
    end
    
    v = ~isnan(Brain_all_74) & ~isnan(Age_all_74) & ~isnan(Group_all_74) & ~isnan(conf_all74);
    
    if sum(v) > 5
        tbl_c = table(Brain_all_74(v), Group_all_74(v), Age_all_74(v), conf_all74(v), ...
            'VariableNames', {'Brain','Group','Age','Confound'});
        mdl_c = fitlm(tbl_c, 'Brain ~ Group + Age + Confound');
        Grp_t(c) = mdl_c.Coefficients.tStat('Group');
        Grp_p(c) = mdl_c.Coefficients.pValue('Group');
    end
    
    fprintf('%-45s | t = %7.3f            | p = %.4f\n', label, Grp_t(c), Grp_p(c));
end

% Save group sensitivity
GrpSensTable = table(Conf_Label, Grp_t, Grp_p, ...
    'VariableNames', {'Confound', 'Group_tstat', 'Group_pvalue'});
outFile_GrpSens = fullfile(mainpath, 'MISBIE_Sensitivity_GroupComparison.xlsx');
writetable(GrpSensTable, outFile_GrpSens);
fprintf('\nGroup sensitivity results saved to:\n%s\n', outFile_GrpSens);


%% ========================================================================
% SENSITIVITY: GDF15 ~ Brain controlling for each confound
% ========================================================================
% Repeats the GDF15 correlation from SimpleCode Section 4
% Runs THREE versions: Whole sample, Patients only, Controls only
% ========================================================================

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  SENSITIVITY: GDF15 ~ N-back Brain WITH CONFOUND CONTROL\n');
fprintf('================================================================\n');

GDF15_all_74 = T_main.t_fasting_plasma(idx_74);
GDF15_log_74 = log10(GDF15_all_74);

% --- Baselines ---
% Whole sample
v_gdf = ~isnan(GDF15_log_74) & ~isnan(Brain_all_74);
[r_gdf_base, p_gdf_base] = corr(GDF15_log_74(v_gdf), Brain_all_74(v_gdf), ...
    'Type', 'Pearson', 'Rows', 'complete');
fprintf('\nBASELINE GDF15-log vs N-back Brain (whole):    r=%.3f, p=%.4f, n=%d\n', ...
    r_gdf_base, p_gdf_base, sum(v_gdf));

% Patients only
v_gdf_pat = ~isnan(GDF15_log_74) & ~isnan(Brain_all_74) & pat_mask;
[r_gdf_pat_base, p_gdf_pat_base] = corr(GDF15_log_74(v_gdf_pat), Brain_all_74(v_gdf_pat), ...
    'Type', 'Pearson', 'Rows', 'complete');
fprintf('BASELINE GDF15-log vs N-back Brain (patients): r=%.3f, p=%.4f, n=%d\n', ...
    r_gdf_pat_base, p_gdf_pat_base, sum(v_gdf_pat));

% Controls only
con_mask = ~pat_mask;
v_gdf_con = ~isnan(GDF15_log_74) & ~isnan(Brain_all_74) & con_mask;
[r_gdf_con_base, p_gdf_con_base] = corr(GDF15_log_74(v_gdf_con), Brain_all_74(v_gdf_con), ...
    'Type', 'Pearson', 'Rows', 'complete');
fprintf('BASELINE GDF15-log vs N-back Brain (controls): r=%.3f, p=%.4f, n=%d\n', ...
    r_gdf_con_base, p_gdf_con_base, sum(v_gdf_con));

% --- Preallocate: whole, patient, control ---
R_gdf_adj     = nan(nConf, 1);  P_gdf_adj     = nan(nConf, 1);  N_gdf_adj     = nan(nConf, 1);
R_gdf_pat_adj = nan(nConf, 1);  P_gdf_pat_adj = nan(nConf, 1);  N_gdf_pat_adj = nan(nConf, 1);
R_gdf_con_adj = nan(nConf, 1);  P_gdf_con_adj = nan(nConf, 1);  N_gdf_con_adj = nan(nConf, 1);

fprintf('\n%-35s | %20s | %20s | %20s\n', ...
    'Confound', 'Whole (r, p, n)', 'Patients (r, p, n)', 'Controls (r, p, n)');
fprintf('%s\n', repmat('-', 1, 105));

for c = 1:nConf
    label  = confounds{c, 1};
    source = confounds{c, 2};
    colspec = confounds{c, 3};
    
    if strcmp(source, 'phenom')
        conf_all74 = double(table2array(T_S_all74(:, colspec)));
    elseif strcmp(source, 'dummy')
        conf_all74 = Dummy_Vars_74.(colspec);
    else
        raw = T_main.(colspec)(idx_74);
        conf_all74 = safe_to_double(raw);
    end
    
    % Whole sample
    v = ~isnan(GDF15_log_74) & ~isnan(Brain_all_74) & ~isnan(conf_all74);
    if sum(v) > 4
        [R_gdf_adj(c), P_gdf_adj(c)] = partialcorr( ...
            GDF15_log_74(v), Brain_all_74(v), conf_all74(v), ...
            'Type', 'Pearson', 'Rows', 'complete');
        N_gdf_adj(c) = sum(v);
    end
    
    % Patients only
    v_p = ~isnan(GDF15_log_74) & ~isnan(Brain_all_74) & ~isnan(conf_all74) & pat_mask;
    if sum(v_p) > 4
        [R_gdf_pat_adj(c), P_gdf_pat_adj(c)] = partialcorr( ...
            GDF15_log_74(v_p), Brain_all_74(v_p), conf_all74(v_p), ...
            'Type', 'Pearson', 'Rows', 'complete');
        N_gdf_pat_adj(c) = sum(v_p);
    end
    
    % Controls only
    v_c = ~isnan(GDF15_log_74) & ~isnan(Brain_all_74) & ~isnan(conf_all74) & con_mask;
    if sum(v_c) > 4
        [R_gdf_con_adj(c), P_gdf_con_adj(c)] = partialcorr( ...
            GDF15_log_74(v_c), Brain_all_74(v_c), conf_all74(v_c), ...
            'Type', 'Pearson', 'Rows', 'complete');
        N_gdf_con_adj(c) = sum(v_c);
    end
    
    fprintf('%-35s | r=%5.3f p=%5.3f n=%2d | r=%5.3f p=%5.3f n=%2d | r=%5.3f p=%5.3f n=%2d\n', ...
        label, ...
        R_gdf_adj(c), P_gdf_adj(c), N_gdf_adj(c), ...
        R_gdf_pat_adj(c), P_gdf_pat_adj(c), N_gdf_pat_adj(c), ...
        R_gdf_con_adj(c), P_gdf_con_adj(c), N_gdf_con_adj(c));
end


%% ========================================================================
% OUTPUT: 5 MAJOR TABLES
% ========================================================================

fprintf('\n\n');
fprintf('================================================================\n');
fprintf('  SAVING 5 OUTPUT TABLES\n');
fprintf('================================================================\n');

% --- TABLE 1: Brain ~ Phenome (age-corrected, patients only) ---
% Already saved above
fprintf('\nTable 1 (Brain~Phenome): %s\n', outFile_Brain);

% --- TABLE 2: Behavior ~ Phenome (age-corrected, patients only) ---
% Already saved above
fprintf('Table 2 (Behavior~Phenome): %s\n', outFile_Behavior);

% --- TABLE 3: Confound sensitivity for NMDAS~Brain (patients only) ---
% Convert logical Sig_Decrease to double (1/0) for clean Excel output
Table3 = table( ...
    Conf_Label, ...
    R_conf_brain, P_conf_brain, ...
    R_adj_brain, P_adj_brain, N_adj_brain, ...
    DeltaR_brain_obs, DeltaR_brain_lo, DeltaR_brain_hi, DeltaR_brain_p, double(DeltaR_brain_sig), ...
    'VariableNames', { ...
        'Confound', ...
        'R_Confound_vs_Brain', 'P_Confound_vs_Brain', ...
        'R_NMDAS_Brain_Adjusted', 'P_NMDAS_Brain_Adjusted', 'N', ...
        'Bootstrap_DeltaR', 'Bootstrap_CI_Lo', 'Bootstrap_CI_Hi', 'Bootstrap_p', 'Sig_Decrease'} ...
);
outFile_T3 = fullfile(mainpath, 'Table3_Sensitivity_Brain_NMDAS.xlsx');
writetable(Table3, outFile_T3);
fprintf('Table 3 (Brain-NMDAS sensitivity): %s\n', outFile_T3);

% --- TABLE 4: Confound sensitivity for NMDAS~Behavior (patients only) ---
Table4 = table( ...
    Conf_Label, ...
    R_conf_beh, P_conf_beh, ...
    R_adj_beh, P_adj_beh, N_adj_beh, ...
    DeltaR_beh_obs, DeltaR_beh_lo, DeltaR_beh_hi, DeltaR_beh_p, double(DeltaR_beh_sig), ...
    'VariableNames', { ...
        'Confound', ...
        'R_Confound_vs_Behavior', 'P_Confound_vs_Behavior', ...
        'R_NMDAS_Behav_Adjusted', 'P_NMDAS_Behav_Adjusted', 'N', ...
        'Bootstrap_DeltaR', 'Bootstrap_CI_Lo', 'Bootstrap_CI_Hi', 'Bootstrap_p', 'Sig_Decrease'} ...
);
outFile_T4 = fullfile(mainpath, 'Table4_Sensitivity_Behavior_NMDAS.xlsx');
writetable(Table4, outFile_T4);
fprintf('Table 4 (Behavior-NMDAS sensitivity): %s\n', outFile_T4);

% --- TABLE 5: GDF15~Brain confound sensitivity (whole, patients, controls) ---
Table5 = table( ...
    Conf_Label, ...
    R_gdf_adj, P_gdf_adj, N_gdf_adj, ...
    R_gdf_pat_adj, P_gdf_pat_adj, N_gdf_pat_adj, ...
    R_gdf_con_adj, P_gdf_con_adj, N_gdf_con_adj, ...
    'VariableNames', { ...
        'Confound', ...
        'R_Whole', 'P_Whole', 'N_Whole', ...
        'R_Patients', 'P_Patients', 'N_Patients', ...
        'R_Controls', 'P_Controls', 'N_Controls'} ...
);
outFile_T5 = fullfile(mainpath, 'Table5_Sensitivity_GDF15_Brain.xlsx');
writetable(Table5, outFile_T5);
fprintf('Table 5 (GDF15-Brain sensitivity): %s\n', outFile_T5);

fprintf('\n========== ALL 5 TABLES SAVED ==========\n');


%% ========================================================================
% HELPER FUNCTION
% ========================================================================

function save_corr_results_excel(varNames, R, P, Pfdr, Nobs, outFile)
    Sig_FDR = double(Pfdr < 0.05);
    ResultTable = table( ...
        varNames(:), R(:), P(:), Pfdr(:), Sig_FDR(:), Nobs(:), ...
        'VariableNames', {'Variable', 'R', 'p', 'p_FDR', 'Sig_FDR', 'N'} );
    writetable(ResultTable, outFile);
    fprintf('Results saved to:\n%s\n', outFile);
end

function out = safe_to_double(x)
% Convert a variable that may be numeric or a cell array of strings to double.
% String categories are converted to integer codes (1, 2, 3, ...) via unique().
% This handles Sex ('M'/'F') and GenotypeClass ('PointMutation','Deletion','MELAS').
    if isnumeric(x) || islogical(x)
        out = double(x);
    elseif iscell(x)
        % Map unique string labels to integers; NaN for empty/missing strings
        [cats, ~, idx] = unique(x);
        out = double(idx);
        % Mark empty strings as NaN
        empty_mask = cellfun(@(s) isempty(s) || (ischar(s) && all(isspace(s))), x);
        out(empty_mask) = NaN;
        % Print the coding so the user can verify
        fprintf('  [safe_to_double] String->numeric coding:\n');
        for k = 1:numel(cats)
            fprintf('    %d = %s\n', k, cats{k});
        end
    else
        out = double(x);
    end
end

function out = safe_to_string(x)
% Ensure x is a cell array of strings regardless of input type.
% Handles: cell of strings, categorical, numeric, char arrays.
    if iscell(x)
        out = x;
    elseif iscategorical(x)
        out = cellstr(x);
    elseif isnumeric(x)
        out = arrayfun(@num2str, x, 'UniformOutput', false);
    elseif ischar(x)
        out = cellstr(x);
    else
        out = cellstr(string(x));
    end
end
