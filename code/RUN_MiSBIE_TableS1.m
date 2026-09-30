%% Reproduce the reconciled Table S1 with its original 211-test FDR family
clearvars;
repoRoot=misbieFindPackageRoot(mfilename('fullpath'));
out=fullfile(repoRoot,'outputs','TableS1');if ~isfolder(out),mkdir(out);end
P=readtable(fullfile(repoRoot,'data','phenotype','phenotype_scores.csv'),'VariableNamingRule','preserve');
G=readtable(fullfile(repoRoot,'data','participants','cohort_groups.csv'));
[found,ix]=ismember(string(P.ReleaseID),string(G.ReleaseID));assert(all(found));g=G.GroupID(ix);
R=readtable(fullfile(repoRoot,'data','summary','table_s1_reconciled.csv'),'TextType','string','VariableNamingRule','preserve');
O=readtable(fullfile(repoRoot,'data','summary','table_s1_original_summary.csv'),'TextType','string','VariableNamingRule','preserve');
assert(height(R)==214 && height(O)==214);
S=R;S.validation=repmat("Original summary retained; exact raw variable absent",214,1);
matched=0;
for i=1:height(R)
    v=R.variable(i);
    orig=[O.N_Control(i),O.N_Patient(i),O.Mean_Control(i),O.Mean_Patient(i),O.SD_Control(i),O.SD_Patient(i)];
    z=orig;gg=O.Hedges_g(i);
    if ismember(v,string(P.Properties.VariableNames))
        matched=matched+1;x=P.(v);a=x(g==1 & isfinite(x));b=x(g==0 & isfinite(x));
        observed=[numel(a),numel(b),mean(a),mean(b),std(a,0),std(b,0)];
        if any(abs(observed-orig)>1e-4)
            z=observed;pool=sqrt(((z(1)-1)*z(5)^2+(z(2)-1)*z(6)^2)/(z(1)+z(2)-2));
            gg=(z(3)-z(4))/pool*(1-3/(4*(z(1)+z(2))-9));
            S.validation(i)="Recalculated from participant records";
        else
            S.validation(i)="Participant records reproduce original descriptives";
        end
    end
    S.n_control(i)=z(1);S.n_patient(i)=z(2);S.mean_control(i)=z(3);S.mean_patient(i)=z(4);
    S.sd_control(i)=z(5);S.sd_patient(i)=z(6);S.g(i)=gg;
    if z(1)>1 && z(2)>1 && all(isfinite(z(3:6)))
        se2=z(5)^2/z(1)+z(6)^2/z(2);
        df=se2^2/((z(5)^2/z(1))^2/(z(1)-1)+(z(6)^2/z(2))^2/(z(2)-1));
        S.p_welch(i)=2*tcdf(-abs((z(3)-z(4))/sqrt(se2)),df);
    else
        S.p_welch(i)=NaN;
    end
end
assert(matched==194,'Unexpected count of available raw variables.');
valid=isfinite(S.p_welch);assert(sum(valid)==211);
[p,order]=sort(S.p_welch(valid));q=flipud(cummin(flipud(p.*numel(p)./(1:numel(p))')));
idx=find(valid);S.q(:)=NaN;S.q(idx(order))=min(1,q);
for name={'g','p_welch','q'}
    x=S.(name{1});y=R.(name{1});m=isfinite(x)&isfinite(y);
    assert(isequal(isfinite(x),isfinite(y)) && max(abs(x(m)-y(m)))<1e-8, ...
        'Recalculation differs from reference: %s',name{1});
end
writetable(S,fullfile(out,'TableS1_reproduced.csv'));
fprintf('Reproduced 214 rows: 194 raw-variable matches; 20 retained summaries; 211 tests.\n');

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
