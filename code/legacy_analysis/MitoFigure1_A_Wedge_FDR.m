% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
%==========================================================================
% Script: Filter by FDR, then plot up to 8 highest‐|Hedges’ g| per subcategory
%==========================================================================
clear; clc;

%% STEP 1: Load the table with FDR information
% Replace with the actual path if needed
inputXLS = 'E:\Mito_DICOM\Figure1_Phenome\Hedgeg_Ke_Preprocessed_WithFDR.xlsx';
T = readtable(inputXLS);

% Verify necessary columns exist
requiredCols = {'Variable','Hedges_g','Subcategory','PassFDR'};
if ~all(ismember(requiredCols, T.Properties.VariableNames))
    error('Table must contain columns: Variable, Hedges_g, Subcategory, PassFDR');
end

%% STEP 2: Keep only rows that passed FDR (q < 0.05)
maskPassed = T.PassFDR == 1;
Tf         = T(maskPassed, :);

% If no rows survived FDR, exit early
if isempty(Tf)
    error('No variables passed FDR < 0.05.');
end

%% STEP 3: For each subcategory, select up to 8 variables with largest |Hedges_g|
allSubfields = unique(Tf.Subcategory, 'stable');
selectedVars = table();  % will accumulate filtered rows

for i = 1:numel(allSubfields)
    sf    = allSubfields{i};
    mask  = strcmp(Tf.Subcategory, sf);
    subset = Tf(mask, :);
    
    % Compute absolute effect sizes and sort descending
    [~, sortIdxDesc] = sort(abs(subset.Hedges_g), 'descend');
    subsetSortedDesc = subset(sortIdxDesc, :);
    
    % Take up to top 8
    nTake = min(11, height(subsetSortedDesc));
    toTake = subsetSortedDesc(1:nTake, :);
    
    % For plotting, we will want them sorted ascending by |Hedges_g|
    [~, sortIdxAsc] = sort(abs(toTake.Hedges_g), 'ascend');
    toTakePlot = toTake(sortIdxAsc, :);
    
    % Append to the overall list
    selectedVars = [selectedVars; toTakePlot];
end

% Now 'selectedVars' contains only PassFDR==1 rows, at most 8 per subcategory,
% sorted ascending by |Hedges_g| within each subcategory.

%% STEP 4: Prepare sorted vectors for plotting
values     = selectedVars.Hedges_g;
subfields  = selectedVars.Subcategory;
absValues  = abs(values);

% Re‐compute unique subfields in the order they appear in 'selectedVars'
uniqueSubfields = unique(subfields, 'stable');
numSubfields    = numel(uniqueSubfields);

% Build sortedAbsValues & sortedSubfields (they’re already grouped & sorted)
sortedAbsValues  = absValues;       % Already ascending within each block
sortedSubfields  = subfields;       % Already in corresponding order

% If desired, verify grouping by subfield:
%   Each block of rows in 'selectedVars' corresponds to one subfield, in the order of uniqueSubfields.

%% STEP 5: Generate a color for each subfield via seaborn_colors
% seaborn_colors should return an N_subfields-by-3 matrix of RGB triplets
subfield_colors = seaborn_colors(numSubfields);

% Build a color cell array matching each entry in 'sortedSubfields'
colors = cell(height(selectedVars),1);
for i = 1:numSubfields
    sf   = uniqueSubfields{i};
    idx  = strcmp(sortedSubfields, sf);
    colors(idx) = repmat(subfield_colors(i,:), sum(idx), 1);
end

%% STEP 6: Create one label per subcategory (at median index in each block)
% Initialize all labels empty
labels = repmat({''}, height(selectedVars), 1);

for i = 1:numSubfields
    sf   = uniqueSubfields{i};
    idx  = find(strcmp(sortedSubfields, sf));
    % pick median index within that block
    mid = idx(round(numel(idx)/2));
    labels{mid} = sf;
end

%% STEP 7: Plot using tor_wedge_plot
figure;
[handles, key_points] = tor_wedge_plot(sortedAbsValues, labels, ...
    'outer_circle_radius', 2, 'colors', colors);
% 
% %% STEP 8: Darken and enlarge labels
% for i = 1:length(handles)
%     if isfield(handles(i), 'texth') && ~isempty(handles(i).texth)
%         lbl = labels{i};
%         if ~isempty(lbl)
%             idxSF = find(strcmp(uniqueSubfields, lbl), 1);
%             if ~isempty(idxSF)
%                 baseColor    = subfield_colors(idxSF, :);
%                 darkerColor  = min(max(baseColor * 0.7, 0), 1);
%                 set(handles(i).texth, 'Color', darkerColor, ...
%                                        'FontSize', 16, ...
%                                        'FontWeight', 'Bold');
%             end
%         end
%     end
% end
% 
% %% (Optional) Display which variables were plotted:
% fprintf('Plotted variables (up to 8 per subfield, all passed FDR):\n');
% for i = 1:height(selectedVars)
%     fprintf('  %s (Subcategory: %s, Hedges_g = %.3f, FDR q = %.3e)\n', ...
%         selectedVars.Variable{i}, selectedVars.Subcategory{i}, ...
%         selectedVars.Hedges_g(i), selectedVars.FDR_PValue(i));
% end
