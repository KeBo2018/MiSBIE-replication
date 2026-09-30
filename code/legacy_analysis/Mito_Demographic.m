% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
% ========================================================================
% compare_groups.m
% 
% Compares "group 0" vs. "groups 1–3 combined" on age, sex, income, BMI,
% adjusted weight, and l_cat, using your Excel sheet:
%   MiSBIEProject-DemographicsForKe_DATA_2025-06-06_17.xlsx
%
% - Continuous variables: two‐sample t‐test (ttest2)
% - Categorical variables (sex, l_cat): chi‐square test (crosstab → chi2, p)
%
% Before running, be sure the Excel file is either in the current folder
% or that you give its full path below.
% ========================================================================

%% 1) READ THE EXCEL FILE
% If needed, replace 'MiSBIEProject-DemographicsForKe_DATA_2025-06-06_17.xlsx'
% with the full path.
clear
fname = 'LOCAL_SOURCE\Dartmouth College Dropbox\Ke Bo\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\MiSBIEProject-DemographicsForKe_DATA_2025-06-06_17.xlsx';

% ========================================================================
% compare_groups.m
% 
% Compares "group 0" vs. "groups 1–3 combined" on age, sex, income, BMI,
% adjusted weight, and l_cat, using your Excel sheet:
%   MiSBIEProject-DemographicsForKe_DATA_2025-06-06_17.xlsx
%
% - Continuous variables: two‐sample t‐test (ttest2)
% - Categorical variables (sex, l_cat): chi‐square test (crosstab → chi2, p)
% - Plots group comparisons using violin plots for continuous variables
%
% Before running, be sure the Excel file is either in the current folder
% or that you give its full path below. Also make sure violinplot.m is on
% your MATLAB path.
% ========================================================================

%% 1) READ THE EXCEL FILE
% fname = 'MiSBIEProject-DemographicsForKe_DATA_2025-06-06_17.xlsx';
T = readtable(fname);

%% 2) EXTRACT GROUP INDICATOR (0 vs. 1–3)
grpAll       = T.pi_geneticdiagnosistype; 
isGroup0     = (grpAll == 0); 
isGroupOther = (grpAll >= 1 & grpAll <= 3);

if any(~(isGroup0 | isGroupOther))
    warning('Some rows have pi_geneticdiagnosistype outside 0–3.');
end

%% 3) PULL OUT EACH VARIABLE
age_all    = T.dcf_age;    % continuous
sex_all    = T.pi_sex;     % categorical (e.g. 1=M, 2=F, etc.)
income_all = T.pi_income;  % ordinal/categorical or numeric coding
bmi_all    = T.dcf_bmi;    % continuous
wgt_all    = T.dcg_adj_wgt;% continuous (adjusted weight)
lcat_all   = T.l_cat;      % categorical (integer codes)

% Split into group0 / groupOther
age_0    = age_all(isGroup0);   age_O    = age_all(isGroupOther);
income_0 = income_all(isGroup0);income_O = income_all(isGroupOther);
bmi_0    = bmi_all(isGroup0);   bmi_O    = bmi_all(isGroupOther);
wgt_0    = wgt_all(isGroup0);   wgt_O    = wgt_all(isGroupOther);
sex_0    = sex_all(isGroup0);   sex_O    = sex_all(isGroupOther);
lcat_0   = lcat_all(isGroup0);  lcat_O   = lcat_all(isGroupOther);

%% 4) TWO‐SAMPLE T‐TESTS FOR CONTINUOUS VARIABLES
fprintf('\n===== TWO‐SAMPLE T‐TESTS (group0 vs. groupOther) =====\n\n');

% 4.1 Age
[h_age, p_age, ci_age, stats_age] = ttest2(age_0, age_O);
fprintf('Age:         t(%d)=%.2f, p=%.3f   (group0 mean=%.2f, groupOther mean=%.2f)\n', ...
    stats_age.df, stats_age.tstat, p_age, mean(age_0), mean(age_O));

% 4.2 Income
[h_inc, p_inc, ci_inc, stats_inc] = ttest2(income_0, income_O);
fprintf('Income:      t(%d)=%.2f, p=%.3f   (group0 mean=%.2f, groupOther mean=%.2f)\n', ...
    stats_inc.df, stats_inc.tstat, p_inc, mean(income_0,'omitnan'), mean(income_O,'omitnan'));

% 4.3 BMI
[h_bmi, p_bmi, ci_bmi, stats_bmi] = ttest2(bmi_0, bmi_O);
fprintf('BMI:         t(%d)=%.2f, p=%.3f   (group0 mean=%.2f, groupOther mean=%.2f)\n', ...
    stats_bmi.df, stats_bmi.tstat, p_bmi, mean(bmi_0,'omitnan'), mean(bmi_O,'omitnan'));

% 4.4 Adjusted weight
[h_wgt, p_wgt, ci_wgt, stats_wgt] = ttest2(wgt_0, wgt_O);
fprintf('AdjWeight:   t(%d)=%.2f, p=%.3f   (group0 mean=%.2f, groupOther mean=%.2f)\n', ...
    stats_wgt.df, stats_wgt.tstat, p_wgt, mean(wgt_0,'omitnan'), mean(wgt_O,'omitnan'));

%% 5) CHI‐SQUARE TESTS FOR CATEGORICAL VARIABLES
fprintf('\n===== CHI‐SQUARE TESTS (group0 vs. groupOther) =====\n\n');

% 5.1 Sex (e.g. 1=M, 2=F, etc.)
[tab_sex, chi2_sex, p_sex] = crosstab(sex_all, isGroupOther);
fprintf('Sex (1 vs. 2):      chi2=%.2f, p=%.3f\n', chi2_sex, p_sex);
disp('    Contingency table for Sex (rows=sex codes, cols=[group0, groupOther]):');
disp(tab_sex);

% 5.2 l_cat (assuming integer codes)
[tab_lcat, chi2_lcat, p_lcat] = crosstab(lcat_all, isGroupOther);
fprintf('l_cat:       chi2=%.2f, p=%.3f\n', chi2_lcat, p_lcat);
disp('    Contingency table for l_cat (rows=l_cat categories, cols=[group0, groupOther]):');
disp(tab_lcat);

%% 6) OPTIONAL: SAVE RESULTS TO A TXT OR MAT FILE
% You could uncomment and adjust these lines if you want to save summary results:
% res.age      = struct('tstat',stats_age.tstat,'df',stats_age.df,'p',p_age, 'mean0',mean(age_0),'meanO',mean(age_O));
% res.income   = struct('tstat',stats_inc.tstat,'df',stats_inc.df,'p',p_inc, 'mean0',mean(income_0),'meanO',mean(income_O));
% res.bmi      = struct('tstat',stats_bmi.tstat,'df',stats_bmi.df,'p',p_bmi, 'mean0',mean(bmi_0),'meanO',mean(bmi_O));
% res.weight   = struct('tstat',stats_wgt.tstat,'df',stats_wgt.df,'p',p_wgt, 'mean0',mean(wgt_0),'meanO',mean(wgt_O));
% res.sex      = struct('chi2',chi2_sex, 'p', p_sex, 'table', tab_sex);
% res.l_cat    = struct('chi2',chi2_lcat,'p', p_lcat,'table', tab_lcat);
% save('group0_vs_other_results.mat','res');

%% 7) PLOT COMPARISONS WITH VIOLIN PLOTS
% For continuous variables (age, income, BMI, adjusted weight), use violinplot()
% For each, package data into a cell array {group0, groupOther}, label x-axis,
% and annotate title with p-value from the t-test above.

figure('Name','Continuous Variable Comparisons (Violin Plots)','NumberTitle','off');
pointsize=5;
dotcolor2 = [0, 114, 189] / 255;    % Control - Blue
dotcolor1 = [217, 83, 25] / 255;    % Patient - Orange/Red
colorcoding=[dotcolor2; dotcolor1];

% 7.1 Age
subplot(1,3,1);
% Combine into cell array
Y_age = { age_0,age_O};
violinplot(Y_age, 'xlabel', {'Group0','GroupOther'}, ...
           'facecolor', colorcoding, ...
           'edgecolor', 'k', 'linewidth', 1.5, ...
           'mc', 'k', 'pointsize',pointsize, 'plotlegend',0);
title(sprintf('Age (p=%.3f)', p_age));
ylabel('Age');
set(gca,'xtick',[])
names={'Control';'MitoD'}
set(gca,'xtick',[1:2],'xticklabel',names)

% 7.2 Income
% subplot(2,2,2);
% Y_inc = { income_0,income_O};
% violinplot(Y_inc, 'xlabel', {'Group0','GroupOther'}, ...
%            'facecolor', colorcoding, ...
%            'edgecolor', 'k', 'linewidth', 1.5, ...
%            'mc', 'k', 'pointsize', pointsize, 'plotlegend',0);
% title(sprintf('Income (p=%.3f)', p_inc));
% ylabel('Income');
% set(gca,'xtick',[])
% names={'Control';'MitoD'}
% set(gca,'xtick',[1:2],'xticklabel',names)
% 7.3 BMI
subplot(1,3,2);
Y_bmi = {bmi_0,bmi_O};
violinplot(Y_bmi, 'xlabel', {'Group0','GroupOther'}, ...
           'facecolor', colorcoding, ...
           'edgecolor', 'k', 'linewidth', 1.5, ...
           'mc', 'k', 'pointsize', pointsize, 'plotlegend',0);
title(sprintf('BMI (p=%.3f)', p_bmi));
ylabel('BMI');
set(gca,'xtick',[])
names={'Control';'MitoD'}
set(gca,'xtick',[1:2],'xticklabel',names)
% 7.4 Adjusted Weight
subplot(1,3,3);
Y_wgt = {wgt_0, wgt_O};
violinplot(Y_wgt, 'xlabel', {'Group0','GroupOther'}, ...
           'facecolor', colorcoding, ...
           'edgecolor', 'k', 'linewidth', 1.5, ...
           'mc', 'k', 'pointsize', pointsize, 'plotlegend',0);
title(sprintf('Adjusted Weight (p=%.3f)', p_wgt));
ylabel('Adjusted Weight');
set(gca,'xtick',[])
names={'Control';'MitoD'}
set(gca,'xtick',[1:2],'xticklabel',names)
fprintf('\nViolin plots generated for continuous variables.\n\n');



%% 9) PIE CHARTS FOR SEX DISTRIBUTION BY GROUP
% tab_sex is sized [nSexCodes × 2] where col-1 = group0, col-2 = groupOther
sexLabels = {'Male','Female'};  % adjust if your coding is different

figure('Name','Sex Distribution by Group','NumberTitle','off');
% Group 0
subplot(1,2,1);
pie(tab_sex(:,1), sexLabels);
title('Group (Control)');

% Group Other
subplot(1,2,2);
pie(tab_sex(:,2), sexLabels);
title('Groups(MitoD)');


fprintf('\nDone.\n');
% ========================================================================
% End of compare_groups.m
% ========================================================================
