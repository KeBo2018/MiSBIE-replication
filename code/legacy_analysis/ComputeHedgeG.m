% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
%=====================================================================
% MATLAB Script to compute Hedge's g and two‐sample t‐test p‐values
% (with a check of the actual column names)
%=====================================================================

% 1. Read in the original data table
inputFilename = 'hedges_g_summary_with_descriptives_and_p_values_CK_Final.xlsx';
T = readtable(inputFilename);

% 2. Display what column names MATLAB actually imported
disp('Imported column names:')
disp(T.Properties.VariableNames')

% Pause so you can copy exactly what each name is (e.g. including spaces)
pause

% 3. Now define the exact names as you saw them above. For example, if
%    you saw a leading space, write it inside the quotes here:
colMean1  = 'Mean_Group1';   % <<— replace with exactly what you saw
colMean2  = 'Mean_Group2';   % <<— replace with exactly what you saw
colStd1   = 'Std_Group1';    % <<— replace with exactly what you saw
colStd2   = 'Std_Group2';    % <<— replace with exactly what you saw
colN1     = 'N_Group1';      % <<— replace with exactly what you saw
colN2     = 'N_Group2';      % <<— replace with exactly what you saw

% 4. Preallocate result vectors
nRows   = height(T);
HedgesG = nan(nRows,1);
P_Value = nan(nRows,1);

% 5. Loop through each row to compute Hedge’s g and p‐value
for i = 1:nRows
    % Extract summary stats from row i using dynamic field‐referencing:
    m1  = T.(colMean1)(i);
    m2  = T.(colMean2)(i);
    sd1 = T.(colStd1)(i);
    sd2 = T.(colStd2)(i);
    n1  = T.(colN1)(i);
    n2  = T.(colN2)(i);
    
    % (a) Compute pooled standard deviation
    df_pool = n1 + n2 - 2;
    if df_pool > 0
        pooledSD = sqrt( ((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2) / df_pool );
    else
        pooledSD = NaN;
    end
    
    % (b) Compute Cohen’s d
    if pooledSD > 0
        d = (m1 - m2) / pooledSD;
    else
        d = NaN;
    end
    
    % (c) Small‐sample correction factor J for Hedge’s g
    denomJ = 4*(n1 + n2) - 9;
    if denomJ ~= 0
        J = 1 - (3 / denomJ);
    else
        J = NaN;
    end
    
    % (d) Hedge's g
    HedgesG(i) = d * J;
    
    % (e) Comp ute two‐sample t statistic (pooled‐variance)
    if pooledSD > 0 && n1 > 0 && n2 > 0
        seMeanDiff = pooledSD * sqrt(1/n1 + 1/n2);
        tStat      = (m1 - m2) / seMeanDiff;
    else
        tStat = NaN;
    end
    
    % (f) Degrees of freedom
    df_t = df_pool;   % = n1 + n2 − 2
    
    % (g) Two‐tailed p‐value
    if ~isnan(tStat) && df_t > 0
        P_Value(i) = 2 * (1 - tcdf(abs(tStat), df_t));
    else
        P_Value(i) = NaN;
    end
end

% 6. Attach new columns
T.HedgesG = HedgesG;
T.P_Value = P_Value;

% 7. Write the augmented table back out
outputFilename = 'hedges_g_with_pvals_results.xlsx';
writetable(T, outputFilename);

fprintf('Done! Saved as "%s".\n', outputFilename);
