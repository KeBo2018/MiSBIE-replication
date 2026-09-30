% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
%% ========================================================================
%  COMBINED SENSITIVITY ANALYSIS: Table S3
%  Regressing out potential confounders from NMDAS~Brain and NMDAS~Behavior
%
%  This script extends Brain_Behavior_Correlate_Phenom_v3.m by adding
%  structural MRI confounders (cortical thickness, FA, NODDI ICVF, NODDI Max)
%  to the existing confound list (sex, genotype, vision, fatigue, neurological,
%  psychiatric), and outputs a single unified Table S3 matching the format
%  of Supplementary Table S2.
%
%  DATA SOURCES:
%    1. DataMatrix_Behavior.xlsx          — brain, behavior, NMDAS, age, sex, genotype
%    2. Ke_DataforCorrelationMatrix_...   — phenome variables (vision, fatigue, etc.)
%    3. WhiteMatter_86Sub.mat             — ValsCT (cortical thickness), ValsFA (FA)
%    4. NODDI_86Sub.mat                   — ValsICVF, ValsMax
%
%  OUTPUT:
%    TableS3_Sensitivity_Brain_NMDAS.xlsx     — NMDAS~Brain sensitivity
%    TableS3_Sensitivity_Behavior_NMDAS.xlsx  — NMDAS~Behavior sensitivity
% ========================================================================

clear; clc; close all;

%% ========================================================================
% 1. LOAD DATA (identical to Brain_Behavior_Correlate_Phenom_v3.m)
% ========================================================================

mainpath = 'LOCAL_SOURCE\Dartmouth College Dropbox\Ke Bo\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\Manuscript\Submission2_Science\Manuscript revision\Analysis';

T_main = readtable(fullfile(mainpath, 'DataMatrix_Behavior.xlsx'));

Multi  = T_main.Multisensory;
Cold   = T_main.Cold;
Nback  = T_main.Nback;
Age    = T_main.Age;
Group  = T_main.GroupID;       % 0 = patient, 1 = control
NMDAS  = T_main.NMDAS;
BehAcc = T_main.Nback_ACC_;
BehOK  = T_main.x_Behavioral_Available_ == 1;

% Subset to 74 N-back subjects
idx_74   = find(BehOK);
Brain_74    = Nback(idx_74);
Behavior_74 = BehAcc(idx_74);
Age_74      = Age(idx_74);
Group_74    = Group(idx_74);
pat_mask    = Group_74 == 0;

% Patient-only subsets
Brain_Patient    = Brain_74(pat_mask);
Behavior_Patient = Behavior_74(pat_mask);
Age_Patient      = Age_74(pat_mask);
NMDAS_pat        = NMDAS(idx_74(pat_mask));

fprintf('N-back subjects: %d (patients: %d, controls: %d)\n', ...
    length(idx_74), sum(pat_mask), sum(~pat_mask));

% Load phenome table
inputFile = 'LOCAL_SOURCE\Dartmouth College Dropbox\Ke Bo\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\Ke_DataforCorrelationMatrix_Cleaned_CK_Combined_SingleItem.xlsx';
T_phenom    = readtable(inputFile);
T_S         = T_phenom(idx_74, :);
T_S_Patient = T_S(pat_mask, :);

%% ========================================================================
% 2. LOAD STRUCTURAL MRI DATA & BUILD PATIENT-LEVEL STRUCTURAL VARIABLES
% ========================================================================
% WhiteMatter_86Sub.mat contains:
%   subject_ids  — cell array of 3-digit ID strings (e.g. '003')
%   ValsCT       — cortical thickness, col 1 = average
%   ValsFA       — fractional anisotropy, col 2 = average
%
% NODDI_86Sub.mat contains:
%   subject_ids  — same format
%   ValsICVF     — intracellular volume fraction, col 2 = average
%   ValsMax      — max NODDI metric, col 3 = average

% Build 3-digit MitoID strings for patients (25 subjects)
patient_row_ids = idx_74(pat_mask);
last_three_ids  = arrayfun(@(x) sprintf('%03d', x), patient_row_ids, ...
                           'UniformOutput', false);

% --- White Matter (CT + FA) ---
load(fullfile(mainpath, 'WhiteMatter_86Sub.mat'));  % loads subject_ids, ValsCT, ValsFA
[~, idx_WM_sub, idx_WM_pat] = intersect(subject_ids, last_three_ids);

CT_Patient  = nan(sum(pat_mask), 1);
FA_Patient  = nan(sum(pat_mask), 1);
CT_Patient(idx_WM_pat)  = ValsCT(idx_WM_sub, 1);
FA_Patient(idx_WM_pat)  = ValsFA(idx_WM_sub, 2);

fprintf('White matter matched: %d/%d patients\n', length(idx_WM_pat), sum(pat_mask));

% --- NODDI (ICVF + Max) ---
load(fullfile(mainpath, 'NODDI_86Sub.mat'));  % loads subject_ids, ValsICVF, ValsMax
[~, idx_ND_sub, idx_ND_pat] = intersect(subject_ids, last_three_ids);

ICVF_Patient = nan(sum(pat_mask), 1);
Max_Patient  = nan(sum(pat_mask), 1);
ICVF_Patient(idx_ND_pat) = ValsICVF(idx_ND_sub, 2);
Max_Patient(idx_ND_pat)  = ValsMax(idx_ND_sub, 3);

fprintf('NODDI matched: %d/%d patients\n', length(idx_ND_pat), sum(pat_mask));

%% ========================================================================
% 3. BUILD GENOTYPE DUMMY VARIABLES
% ========================================================================

genoDiag_pat = T_main.geneticDiagnosis(idx_74(pat_mask));
genoDiag_pat = safe_to_string(genoDiag_pat);

fprintf('\n[Genotype] Unique values: '); disp(unique(genoDiag_pat)');

Dummy_Vars.is_PointMutation = double(strcmpi(genoDiag_pat, 'PointMutation')   | ...
                                     strcmpi(genoDiag_pat, 'Point Mutation')  | ...
                                     strcmpi(genoDiag_pat, 'point_mutation'));
Dummy_Vars.is_Deletion      = double(strcmpi(genoDiag_pat, 'Deletion')        | ...
                                     strcmpi(genoDiag_pat, 'Single Deletion') | ...
                                     strcmpi(genoDiag_pat, 'single_deletion'));
Dummy_Vars.is_MELAS         = double(strcmpi(genoDiag_pat, 'MELAS'));

%% ========================================================================
% 4. DEFINE CONFOUND LIST
%    Each row: { 'Label', source, column_spec }
%      'main'     -> T_main column name (string)
%      'phenom'   -> T_S_Patient column index (integer)
%      'dummy'    -> Dummy_Vars field name (string)
%      'struct'   -> pre-built patient-level vector variable name (string)
% ========================================================================

confounds = {
    % ---- DEMOGRAPHIC ----
    'Sex',                                           'main',   'sex';
    'Genotype: Point Mutation (vs others)',          'dummy',  'is_PointMutation';
    'Genotype: Single Deletion (vs others)',         'dummy',  'is_Deletion';
    'Genotype: MELAS (vs others)',                   'dummy',  'is_MELAS';
    % ---- VISION ----
    'nmdas_cf_1 (Vision self-report)',               'phenom', 377;
    'nmdas_cca_1 (Visual acuity clinical)',          'phenom', 396;
    'nmdas_cca_2 (Ptosis)',                          'phenom', 397;
    'nmdas_cca_3 (CPEO)',                            'phenom', 398;
    'cns_eyes (CNS Eyes exam)',                      'phenom', 210;
    'cns_crnl_nrv_ii (Cranial nerve II)',            'phenom', 248;
    'cns_crnl_nrv_iii (Cranial nerve III)',          'phenom', 250;
    'cns_crnl_nrv_iv (Cranial nerve IV)',            'phenom', 252;
    'cns_crnl_nrv_vi (Cranial nerve VI)',            'phenom', 256;
    % ---- FATIGUE ----
    'pfs_phys (Physical fatigue)',                   'phenom',  48;
    'pfs_mental (Mental fatigue)',                   'phenom',  49;
    % ---- NEUROLOGICAL / MOTOR ----
    'nmdas_cf_9 (Exercise tolerance)',               'phenom', 385;
    'nmdas_cf_10 (Gait stability)',                  'phenom', 386;
    'nmdas_cca_5 (Myopathy)',                        'phenom', 400;
    'nmdas_cca_6 (Cerebellar ataxia)',               'phenom', 401;
    'nmdas_cca_7 (Neuropathy)',                      'phenom', 402;
    'cns_gait_statn (Gait and station)',             'phenom', 270;
    % ---- PSYCHIATRIC ----
    'nmdas_ssi_1 (Psychiatric symptoms)',            'phenom', 387;
    % ---- STRUCTURAL MRI (White Matter) ----
    'Cortical Thickness (average, DTI)',             'struct', 'CT_Patient';
    'Fractional Anisotropy (average, DTI)',          'struct', 'FA_Patient';
    % ---- STRUCTURAL MRI (NODDI) ----
    'ICVF (intracellular volume fraction, NODDI)',   'struct', 'ICVF_Patient';
    'Max NODDI metric',                              'struct', 'Max_Patient';
};

nConf = size(confounds, 1);

%% ========================================================================
% 5. BASELINE (Age only)
% ========================================================================

valid_base_brain = ~isnan(Brain_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient);
[r_base_brain, p_base_brain] = partialcorr( ...
    Brain_Patient(valid_base_brain), NMDAS_pat(valid_base_brain), ...
    Age_Patient(valid_base_brain), 'Type','Spearman','Rows','complete');

valid_base_beh = ~isnan(Behavior_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient);
[r_base_beh, p_base_beh] = partialcorr( ...
    Behavior_Patient(valid_base_beh), NMDAS_pat(valid_base_beh), ...
    Age_Patient(valid_base_beh), 'Type','Spearman','Rows','complete');

fprintf('\nBASELINE (Age only):\n');
fprintf('  NMDAS vs Brain:    r=%.3f, p=%.4f, n=%d\n', r_base_brain, p_base_brain, sum(valid_base_brain));
fprintf('  NMDAS vs Behavior: r=%.3f, p=%.4f, n=%d\n', r_base_beh,   p_base_beh,   sum(valid_base_beh));

%% ========================================================================
% 6. SENSITIVITY LOOP: partial correlations + bootstrap
% ========================================================================

% Preallocate
Conf_Label       = cell(nConf, 1);
R_conf_brain     = nan(nConf,1); P_conf_brain  = nan(nConf,1);
R_conf_beh       = nan(nConf,1); P_conf_beh    = nan(nConf,1);
R_adj_brain      = nan(nConf,1); P_adj_brain   = nan(nConf,1); N_adj_brain = nan(nConf,1);
R_adj_beh        = nan(nConf,1); P_adj_beh     = nan(nConf,1); N_adj_beh   = nan(nConf,1);

DeltaR_brain_obs = nan(nConf,1); DeltaR_brain_lo = nan(nConf,1);
DeltaR_brain_hi  = nan(nConf,1); DeltaR_brain_p  = nan(nConf,1);
DeltaR_brain_sig = false(nConf,1);

DeltaR_beh_obs   = nan(nConf,1); DeltaR_beh_lo   = nan(nConf,1);
DeltaR_beh_hi    = nan(nConf,1); DeltaR_beh_p    = nan(nConf,1);
DeltaR_beh_sig   = false(nConf,1);

nBoot = 5000;
rng(42);

fprintf('\n========== SENSITIVITY LOOP ==========\n');
fprintf('%-48s | %12s | %12s | %12s | %12s\n', ...
    'Confound', 'r_Conf~Brain', 'r_Adj~Brain', 'r_Conf~Beh', 'r_Adj~Beh');
fprintf('%s\n', repmat('-', 1, 105));

for c = 1:nConf
    label   = confounds{c,1};
    source  = confounds{c,2};
    colspec = confounds{c,3};
    Conf_Label{c} = label;

    % --- Extract confound variable ---
    conf_var = extract_confound(source, colspec, T_S_Patient, T_main, ...
                                idx_74, pat_mask, Dummy_Vars);

    % --- (A) Confound vs Brain (age-corrected) ---
    [R_conf_brain(c), P_conf_brain(c)] = safe_partial(conf_var, Brain_Patient,    Age_Patient);

    % --- (B) Confound vs Behavior (age-corrected) ---
    [R_conf_beh(c),   P_conf_beh(c)  ] = safe_partial(conf_var, Behavior_Patient, Age_Patient);

    % --- (C) NMDAS vs Brain | Age + confound ---
    v = ~isnan(Brain_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient) & ~isnan(conf_var);
    if sum(v) > 4
        [R_adj_brain(c), P_adj_brain(c)] = partialcorr( ...
            Brain_Patient(v), NMDAS_pat(v), [Age_Patient(v), conf_var(v)], ...
            'Type','Spearman','Rows','complete');
        N_adj_brain(c) = sum(v);
    end

    % --- (D) NMDAS vs Behavior | Age + confound ---
    v = ~isnan(Behavior_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient) & ~isnan(conf_var);
    if sum(v) > 4
        [R_adj_beh(c), P_adj_beh(c)] = partialcorr( ...
            Behavior_Patient(v), NMDAS_pat(v), [Age_Patient(v), conf_var(v)], ...
            'Type','Spearman','Rows','complete');
        N_adj_beh(c) = sum(v);
    end

    % --- (E) Bootstrap: Brain ---
    v = ~isnan(Brain_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient) & ~isnan(conf_var);
    if sum(v) > 5
        [DeltaR_brain_obs(c), DeltaR_brain_lo(c), DeltaR_brain_hi(c), ...
         DeltaR_brain_p(c),   DeltaR_brain_sig(c)] = ...
            bootstrap_delta(Brain_Patient(v), NMDAS_pat(v), Age_Patient(v), conf_var(v), nBoot);
    end

    % --- (F) Bootstrap: Behavior ---
    v = ~isnan(Behavior_Patient) & ~isnan(NMDAS_pat) & ~isnan(Age_Patient) & ~isnan(conf_var);
    if sum(v) > 5
        [DeltaR_beh_obs(c), DeltaR_beh_lo(c), DeltaR_beh_hi(c), ...
         DeltaR_beh_p(c),   DeltaR_beh_sig(c)] = ...
            bootstrap_delta(Behavior_Patient(v), NMDAS_pat(v), Age_Patient(v), conf_var(v), nBoot);
    end

    fprintf('  %-46s | %+9.3f     | %+9.3f     | %+9.3f     | %+9.3f\n', ...
        label, R_conf_brain(c), R_adj_brain(c), R_conf_beh(c), R_adj_beh(c));
end

%% ========================================================================
% 7. SAVE OUTPUT TABLES (matching Table S2 format)
% ========================================================================

% Table S3a: NMDAS ~ Brain
TableS3a = table( ...
    Conf_Label, ...
    R_conf_brain, P_conf_brain, ...
    R_adj_brain, P_adj_brain, N_adj_brain, ...
    DeltaR_brain_obs, DeltaR_brain_lo, DeltaR_brain_hi, ...
    DeltaR_brain_p, double(DeltaR_brain_sig), ...
    'VariableNames', { ...
        'Confound', ...
        'R_Confound_vs_Brain', 'P_Confound_vs_Brain', ...
        'R_NMDAS_Brain_Adjusted', 'P_NMDAS_Brain_Adjusted', 'N', ...
        'Bootstrap_DeltaR', 'Bootstrap_CI_Lo', 'Bootstrap_CI_Hi', ...
        'Bootstrap_p', 'Sig_Decrease'});

outFile_S3a = fullfile(mainpath, 'TableS3_Sensitivity_Brain_NMDAS.xlsx');
writetable(TableS3a, outFile_S3a);
fprintf('\nTable S3a saved: %s\n', outFile_S3a);

% Table S3b: NMDAS ~ Behavior
TableS3b = table( ...
    Conf_Label, ...
    R_conf_beh, P_conf_beh, ...
    R_adj_beh, P_adj_beh, N_adj_beh, ...
    DeltaR_beh_obs, DeltaR_beh_lo, DeltaR_beh_hi, ...
    DeltaR_beh_p, double(DeltaR_beh_sig), ...
    'VariableNames', { ...
        'Confound', ...
        'R_Confound_vs_Behavior', 'P_Confound_vs_Behavior', ...
        'R_NMDAS_Behav_Adjusted', 'P_NMDAS_Behav_Adjusted', 'N', ...
        'Bootstrap_DeltaR', 'Bootstrap_CI_Lo', 'Bootstrap_CI_Hi', ...
        'Bootstrap_p', 'Sig_Decrease'});

outFile_S3b = fullfile(mainpath, 'TableS3_Sensitivity_Behavior_NMDAS.xlsx');
writetable(TableS3b, outFile_S3b);
fprintf('Table S3b saved: %s\n', outFile_S3b);

fprintf('\nBaseline (Brain):    r=%.3f, p=%.4f, n=%d\n', r_base_brain, p_base_brain, sum(valid_base_brain));
fprintf('Baseline (Behavior): r=%.3f, p=%.4f, n=%d\n', r_base_beh,   p_base_beh,   sum(valid_base_beh));
fprintf('\n========== DONE ==========\n');


%% ========================================================================
% HELPER FUNCTIONS
% ========================================================================

function conf_var = extract_confound(source, colspec, T_S_Patient, T_main, ...
                                     idx_74, pat_mask, Dummy_Vars)
% Extract a confound variable for the patient subset based on source type.
    switch source
        case 'phenom'
            conf_var = double(table2array(T_S_Patient(:, colspec)));
        case 'dummy'
            conf_var = Dummy_Vars.(colspec);
        case 'struct'
            % Variable already built as a patient-level vector in the workspace
            % Use evalin to access caller workspace variable by name
            conf_var = evalin('caller', colspec);
        case 'main'
            raw = T_main.(colspec)(idx_74(pat_mask));
            conf_var = safe_to_double(raw);
        otherwise
            error('Unknown source type: %s', source);
    end
end

function [R, P] = safe_partial(conf_var, outcome, age)
% Partial correlation of conf_var ~ outcome | age.
% Falls back to simple corr for binary/near-zero variance variables.
    v = ~isnan(conf_var) & ~isnan(outcome) & ~isnan(age);
    R = nan; P = nan;
    if sum(v) < 4; return; end
    cv = conf_var(v);
    if numel(unique(cv)) < 3 || var(cv) < 1e-10
        [R, P] = corr(cv, outcome(v), 'Type','Spearman','Rows','complete');
    else
        [R, P] = partialcorr(cv, outcome(v), age(v), ...
                             'Type','Spearman','Rows','complete');
    end
end

function [delta_obs, ci_lo, ci_hi, boot_p, sig] = ...
        bootstrap_delta(Brain, NMDAS, Age, Conf, nBoot)
% Bootstrap test of whether controlling Conf significantly attenuates NMDAS~Brain.
% Returns observed delta, 95% CI, bootstrap p, and Sig_Decrease flag.
    n = length(Brain);
    r_base = partialcorr(Brain, NMDAS, Age,         'Type','Spearman','Rows','complete');
    r_adj  = partialcorr(Brain, NMDAS, [Age, Conf], 'Type','Spearman','Rows','complete');
    delta_obs = abs(r_base) - abs(r_adj);

    delta_boot = nan(nBoot, 1);
    for b = 1:nBoot
        idx_b  = randsample(n, n, true);
        r_b    = partialcorr(Brain(idx_b), NMDAS(idx_b), Age(idx_b), ...
                             'Type','Spearman','Rows','complete');
        r_ba   = partialcorr(Brain(idx_b), NMDAS(idx_b), [Age(idx_b), Conf(idx_b)], ...
                             'Type','Spearman','Rows','complete');
        delta_boot(b) = abs(r_b) - abs(r_ba);
    end
    ci_lo  = prctile(delta_boot,  2.5);
    ci_hi  = prctile(delta_boot, 97.5);
    boot_p = mean(delta_boot <= 0);
    sig    = ci_lo > 0;
end

function out = safe_to_double(x)
% Convert numeric or string cell array to double. Strings mapped to integer codes.
    if isnumeric(x) || islogical(x)
        out = double(x);
    elseif iscell(x)
        [cats, ~, idx] = unique(x);
        out = double(idx);
        empty_mask = cellfun(@(s) isempty(s) || (ischar(s) && all(isspace(s))), x);
        out(empty_mask) = NaN;
        fprintf('  [safe_to_double] String coding: ');
        for k = 1:numel(cats); fprintf('%d=%s  ', k, cats{k}); end; fprintf('\n');
    else
        out = double(x);
    end
end

function out = safe_to_string(x)
% Ensure output is a cell array of strings regardless of input type.
    if iscell(x);         out = x;
    elseif iscategorical(x); out = cellstr(x);
    elseif isnumeric(x);  out = arrayfun(@num2str, x, 'UniformOutput', false);
    elseif ischar(x);     out = cellstr(x);
    else;                  out = cellstr(string(x));
    end
end
