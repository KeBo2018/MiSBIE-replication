% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
%=====================================================================
% Combined Script:
%   1) Load original Excel (“Hedgeg_Ke_Preprocessed_ItemLevel.xlsx”)
%   2) Compute FDR‐adjusted p‐values and a PassFDR flag
%   3) Write a new Excel with all original columns plus FDR_PValue and PassFDR
%   4) Plot Hedges’ g in sorted order by Subcategory, coloring bars by
%      subcategory unless they fail FDR correction (adjusted p ≥ 0.05),
%      in which case they are colored grey.
%
% Assumes:
%   • The original Excel has columns named exactly:
%       Variable (string or cell), Hedges_g (numeric), Subcategory (string),
%       P_Value (numeric), plus whatever other columns you already have.
%   • MATLAB’s Statistics and Machine Learning Toolbox is installed (for mafdr).
%=====================================================================

clear; clc;

%% --------------------------------------------------------------------
% STEP 1: Load the original Excel file, which already contains
%         Variable | Hedges_g | Subcategory | P_Value | … (other columns)
%% --------------------------------------------------------------------
inputXLS = 'E:\Mito_DICOM\Figure1_Phenome\Hedgeg_Ke_Preprocessed_ItemLevel.xlsx';
T = readtable(inputXLS);

% Verify that P_Value exists
if ~ismember('P_Value', T.Properties.VariableNames)
    error('Table must contain a column named ''P_Value''.');
end

%% --------------------------------------------------------------------
% STEP 2: Compute FDR‐adjusted p‐values (Benjamini–Hochberg) and a PassFDR flag
%% --------------------------------------------------------------------
rawP   = T.P_Value;                          % raw p‐values
adjP   = mafdr(rawP, 'BHFDR', true);          % adjusted p‐values
passF  = double(adjP < 0.05);                 % 1 if adjP < 0.05, else 0

% Append new columns to the table T
T.FDR_PValue = adjP;
T.PassFDR    = passF;

%% --------------------------------------------------------------------
% STEP 3: Write out a new Excel file that contains all original columns
%         plus the two new columns (FDR_PValue, PassFDR)
%% --------------------------------------------------------------------
outputXLS = 'E:\Mito_DICOM\Figure1_Phenome\Hedgeg_Ke_Preprocessed_WithFDR.xlsx';
writetable(T, outputXLS);
fprintf('Wrote augmented Excel to:\n  %s\n', outputXLS);

%% --------------------------------------------------------------------
% STEP 4: Plotting code (as before), now using T.P_Value and adjP
%% --------------------------------------------------------------------

% Extract relevant columns for plotting
variables   = T.Variable;         % cell array or string array of variable names
values      = T.Hedges_g;         % numeric vector of Hedges’ g
subfields   = T.Subcategory;      % cell array or string array of subcategories
pvals       = T.P_Value;          % raw p‐values (numeric vector)

% Compute significance mask using the already‐computed adjP
isSignificant = adjP < 0.05;      % logical vector, true if passes FDR

% Prepare for sorting within each subcategory
uniqueSubfields = unique(subfields, 'stable');
numSubfields    = numel(uniqueSubfields);

sortedAbsValues = [];
sortedSubfields = string.empty(0,1);
sortedSigFlags  = [];   % will hold the significance flags (true/false) after sorting

subfieldMidpoints = zeros(numSubfields, 1);
runningIndex      = 0;

for i = 1:numSubfields
    thisField = uniqueSubfields{i};
    
    % Find all rows in this subcategory
    maskThis       = strcmp(subfields, thisField);
    fieldValues    = values(maskThis);
    fieldSigFlags  = isSignificant(maskThis);
    
    % Sort by absolute effect size (ascending)
    [sortedVals, sortOrder] = sort(abs(fieldValues));
    
    % Reorder the significance flags accordingly
    sortedFieldSig = fieldSigFlags(sortOrder);
    
    % Append sorted data to the overall lists
    countThis         = numel(sortedVals);
    sortedAbsValues   = [sortedAbsValues; sortedVals];
    sortedSubfields   = [sortedSubfields; repmat(string(thisField), countThis, 1)];
    sortedSigFlags    = [sortedSigFlags; sortedFieldSig];
    
    % Compute midpoint for xtick labeling
    runningIndex                = runningIndex + countThis;
    subfieldMidpoints(i)        = runningIndex - floor(countThis/2);
end

% Generate a colormap for the subcategories
colors = lines(numSubfields);

% Plot bars one by one
figure;
hold on;

totalBars = numel(sortedAbsValues);
for j = 1:totalBars
    if ~sortedSigFlags(j)
        % Grey if not significant
        bar(j, sortedAbsValues(j), 'FaceColor', [0.5, 0.5, 0.5]);
    else
        % Find which subcategory this bar belongs to
        idxField = find(strcmp(uniqueSubfields, sortedSubfields(j)));
        bar(j, sortedAbsValues(j), 'FaceColor', colors(idxField, :));
    end
end

% Customize axes and labels
ylabel('|\itHedges''\it\itg\rm|', 'FontSize', 14);
xticks(subfieldMidpoints);
xticklabels(uniqueSubfields);
xtickangle(30);
set(gca, 'FontSize', 12);

hold off;
