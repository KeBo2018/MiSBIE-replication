%% Supplementary numerical analyses using the released data
% Reproduces participant-level tests where source variables are available.
% All originally reported Tables S1-S5 are also supplied separately.
% Missing raw variables, including unavailable NODDI summaries, are explicitly
% marked rather than replaced by zeros or another measure.
clearvars;
repoRoot=misbieFindPackageRoot(mfilename('fullpath'));
addpath(fullfile(repoRoot,'code'),'-begin');
out=fullfile(repoRoot,'outputs','Supplementary');if ~isfolder(out),mkdir(out);end
T=readtable(fullfile(repoRoot,'data','participants','analysis_data.csv'),'VariableNamingRule','preserve');
P=readtable(fullfile(repoRoot,'data','phenotype','phenotype_scores.csv'),'VariableNamingRule','preserve');
I=readtable(fullfile(repoRoot,'data','phenotype','clinical_items.csv'),'VariableNamingRule','preserve');
[found,ix]=ismember(string(T.ReleaseID),string(P.ReleaseID));assert(all(found));P=P(ix,:);
[found,ix]=ismember(string(T.ReleaseID),string(I.ReleaseID));assert(all(found));I=I(ix,:);
pat=T.GroupID==0;control=T.GroupID==1;eligible=T.Behavioral_Available==1&isfinite(T.Nback);
% Recompute each of the 304 reported phenome associations, keeping its name.
for outcome=1:2
    old=readtable(fullfile(repoRoot,'data','summary',sprintf('reported_Table_S%d.csv',outcome+1)), ...
        'VariableNamingRule','preserve','TextType','string');
    y=T.Nback;if outcome==2,y=T.Nback_ACC;end
    n=height(old);R=nan(n,1);p=R;N=zeros(n,1);status=repmat("Raw variable unavailable",n,1);
    variables=string(old{:,1});
    for j=1:n
        v=variables(j);v=regexprep(v,'^.*\.','');
        if ismember(v,string(P.Properties.VariableNames)),x=P.(v);
        elseif ismember(v,string(I.Properties.VariableNames)),x=I.(v);
        else,continue;end
        m=pat&eligible&isfinite(x)&isfinite(y)&isfinite(T.Age);N(j)=sum(m);
        if sum(m)<5 || numel(unique(x(m)))<2,status(j)="Insufficient variation or observations";continue;end
        [R(j),p(j)]=partialcorr(x(m),y(m),T.Age(m),'Type','Spearman');status(j)="Recomputed from released observations";
    end
    % The archived fdr() output differs numerically from standard BH q values.
    % Retain it explicitly alongside BH, not as a silently replaced column.
    qBH=bh(p);reportedQ=old.p_FDR;
    assert(all(isfinite(p)) && numel(p)==304,'The supplied release supports all 304 tests.');
    assert(all(abs(R-old.R)<.00051) && all(N==old.N) && all(abs(p-old.p)<.000051), ...
        'Correlations differ from the reported table.');
    assert(isequal(qBH<.05,reportedQ<.05),'FDR significance classifications changed.');
    S=table(variables,N,R,p,qBH,reportedQ,status,'VariableNames', ...
        {'Variable','N','RhoAge','PAge','QBH304Tests','ReportedFDRValue','Status'});
    writetable(S,fullfile(out,sprintf('TableS%d_available_raw_variables.csv',outcome+1)));
end
% Motion: same imaging samples as the manuscript, before behavior exclusion.
rows={};names={'Nback','Multisensory','Cold'};
for task=1:3
    for metric={'meanFD','meanDVARS'}
        key=[names{task} '_' metric{1}];x=T.(key);m=isfinite(T.(names{task}))&isfinite(x);
        [~,p,~,st]=ttest2(x(m&control),x(m&pat));
        rows(end+1,:)={names{task},metric{1},sum(m&control),sum(m&pat),st.tstat,st.df,p}; %#ok<SAGROW>
    end
end
writetable(cell2table(rows,'VariableNames',{'Task','Metric','NControl','NPatient','T','DF','P'}),fullfile(out,'head_motion_tests.csv'));
% Inclusion sensitivity: compare the analysis samples without changing flags.
rows={};
for useEligibility=[false true]
    m=isfinite(T.Nback)&pat&isfinite(T.NMDAS)&isfinite(T.Age);
    if useEligibility,m=m&eligible;end
    for outcome=1:2
        y=T.Nback;if outcome==2,y=T.Nback_ACC;end
        v=m&isfinite(y);[r,p]=partialcorr(T.NMDAS(v),y(v),T.Age(v),'Type','Spearman');
        rows(end+1,:)={useEligibility,outcome,sum(v),r,p}; %#ok<SAGROW>
    end
    % The source refits k=2 within each available sample for sensitivity.
    rng(20260923,'twister');[id,c]=kmeans_Matlab(T.NMDAS(m),2,'Replicates',100);
    [~,hi]=max(c);idx=find(m);severe=false(height(T),1);severe(idx(id==hi))=true;
    con=control&isfinite(T.Nback);if useEligibility,con=con&eligible;end
    [~,p,~,st]=ttest2(T.Nback(severe),T.Nback(con));
    fprintf('Eligibility applied=%d: severe/control brain t=%.4f, p=%.6g.\n',useEligibility,st.tstat,p);
end
writetable(cell2table(rows,'VariableNames',{'BehaviorEligibilityApplied','Outcome1Brain2Behavior','N','RhoAge','PAge'}),fullfile(out,'Nback_inclusion_sensitivity.csv'));
% Age plus individual covariate sensitivity on patients.
S=readtable(fullfile(repoRoot,'data','participants','structural_summaries.csv'));
[found,ix]=ismember(string(T.ReleaseID),string(S.ReleaseID));ct=nan(height(T),1);fa=ct;
ct(found)=S.ValsCT_1(ix(found));fa(found)=S.ValsFA_2(ix(found));
covs=[T.NMDAS_Vision,T.NMDAS_Psychiatric,T.PhysicalFatigability,T.MentalFatigability,T.CNS_Eyes,ct,fa];
covNames={'NMDAS vision','NMDAS psychiatric','Physical fatigability','Mental fatigability','CNS eyes','Cortical thickness','Fractional anisotropy'};
rows={};
for j=1:size(covs,2)
    for outcome=1:2
        y=T.Nback;if outcome==2,y=T.Nback_ACC;end
        m=pat&eligible&all(isfinite([y,T.NMDAS,T.Age,covs(:,j)]),2);
        [r,p]=partialcorr(T.NMDAS(m),y(m),[T.Age(m),covs(m,j)],'Type','Spearman');
        rows(end+1,:)={covNames{j},outcome,sum(m),r,p}; %#ok<SAGROW>
    end
end
writetable(cell2table(rows,'VariableNames',{'AdditionalCovariate','Outcome1Brain2Behavior','N','Rho','P'}),fullfile(out,'available_covariate_sensitivity.csv'));
fprintf('NODDI ICVF/Max source summaries were not located; those tests were not rerun.\n');
function q=bh(p)
    valid=isfinite(p);q=nan(size(p));ix=find(valid);[s,o]=sort(p(valid));
    qq=flipud(cummin(flipud(s.*numel(s)./(1:numel(s))')));q(ix(o))=min(1,qq);
end

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
