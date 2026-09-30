%% Check the MiSBIE package inputs before running any analysis
% Base MATLAB only. This check does not calculate new scientific results.
repoRoot=misbieFindPackageRoot(mfilename('fullpath'));
T=readtable(fullfile(repoRoot,'data','participants','analysis_data.csv'));
P=readtable(fullfile(repoRoot,'data','phenotype','phenotype_scores.csv'));
G=readtable(fullfile(repoRoot,'data','participants','cohort_groups.csv'));
assert(height(T)==110 && height(P)==110 && height(G)==110, ...
    'MiSBIE:InputMismatch','Expected the 110-participant release tables.');
assert(numel(unique(string(T.ReleaseID)))==110 && ...
    isequal(sort(string(T.ReleaseID)),sort(string(P.ReleaseID))) && ...
    isequal(sort(string(T.ReleaseID)),sort(string(G.ReleaseID))), ...
    'MiSBIE:InputMismatch','Release participant identifiers do not align.');
S=load(fullfile(repoRoot,'data','figure2','Paired_condition_scores.mat'));
keys={'Nback','Multisensory','Cold'};counts=[88 90 91];
for task=1:3
    assert(isfield(S,[keys{task} '_Control']) && isfield(S,[keys{task} '_Task']), ...
        'MiSBIE:InputMismatch','Missing paired fields for %s.',keys{task});
    assert(numel(S.([keys{task} '_Control']))==counts(task) && ...
        numel(S.([keys{task} '_Task']))==counts(task), ...
        'MiSBIE:InputMismatch','Unexpected paired sample size for %s.',keys{task});
end
fprintf('PASS: all %d required data files found; participant tables and paired scores loaded.\n', ...
    numel(misbieRequiredDataFiles()));
addpath(fullfile(repoRoot,'code'),'-begin');
assert(isfile(fullfile(repoRoot,'code','kmeans_Matlab.m')), ...
    'MiSBIE:MissingKmeansHelper','Extract code/kmeans_Matlab.m from the complete new ZIP.');
fprintf('kmeans_Matlab helper: %s\n',which('kmeans_Matlab'));
if exist('partialcorr','file')==2
    [kmLabels,kmCenters]=kmeans_Matlab([0;1;2;100;101;102],2, ...
        'Start',[1;101],'Replicates',1);
    assert(max(abs(sort(kmCenters)-[1;101]))<1e-10 && ...
        all(kmLabels(1:3)==kmLabels(1)) && all(kmLabels(4:6)==kmLabels(4)) && ...
        kmLabels(1)~=kmLabels(4),'MiSBIE:KmeansCheckFailed', ...
        'The supplied kmeans_Matlab function failed the fixed-start check.');
    fprintf('PASS: kmeans_Matlab fixed-start clustering check.\n');
end
fprintf('MATLAB: %s\n',version);
fprintf('partialcorr: %s\n',which('partialcorr'));
fprintf('mediation: %s\n',which('mediation'));
fprintf('Output folder: %s\n',fullfile(repoRoot,'outputs'));
if exist('partialcorr','file')~=2
    warning('MiSBIE:MissingStatisticsToolbox', ...
        'Data setup passed. Analyses other than Figure 2 require Statistics and Machine Learning Toolbox.');
end
if exist('mediation','file')~=2
    warning('MiSBIE:MissingCANlab', ...
        'CANlab mediation.m was not found. Figure 3 will label its alternative bootstrap method.');
end
fprintf('Open scripts from: %s\n',fullfile(repoRoot,'code'));

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
