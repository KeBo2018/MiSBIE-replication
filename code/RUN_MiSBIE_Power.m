%% MiSBIE: independent MitoD-versus-control power analysis
% Click Run. Requires Statistics and Machine Learning Toolbox.
% Uses the bundled analysis CSV and the same task-specific eligibility as Figure 3.
% RequiredNPerGroup is for a balanced TWO-GROUP study; RequiredTotalN is twice it.
% Actual-sample power uses the observed unequal patient/control sample sizes.
% Task dz values are fixed sensitivity benchmarks, NOT observed between-group d.
% The two standardizations differ: dz uses the SD of paired condition differences;
% independent-groups d uses the pooled within-group SD of task-contrast scores.
clearvars; clc;
scriptFolder=fileparts(mfilename('fullpath'));
repoRoot=misbieFindPackageRoot(mfilename('fullpath'));
dataFile=fullfile(repoRoot,'data','participants','analysis_data.csv');
assert(isfile(dataFile),'Bundled data file is missing: %s',dataFile);
T=readtable(dataFile,'VariableNamingRule','preserve');
Group=column(T,'GroupID'); valid=ismember(Group,[0 1]);
Brain=[column(T,'Nback'),column(T,'Multisensory'),column(T,'Cold')];
eligible=isfinite(Brain)&valid;
eligible(:,1)=eligible(:,1)&column(T,'Behavioral_Available')==1;
taskNames={'Working memory';'Multisensory';'Cold pain'};
contrastNames={'2-back minus 0-back';'Multisensory minus rest';'Cold minus room temperature'};
taskDz=[1.67;2.18;.92]; % Original manuscript: mean paired effects across CV repetitions.
fractions=[.25 .5]; alpha=.05; targetPower=.80;
rows=cell(6,14);k=0;
for task=1:3
    nPatient=sum(eligible(:,task)&Group==0);
    nControl=sum(eligible(:,task)&Group==1);
    assert(nPatient>=2&&nControl>=2,'Insufficient valid observations.');
    detectableD=fzero(@(d) twoGroupPower(d,nPatient,nControl,alpha)-targetPower,[.001 5]);
    for fraction=fractions
        k=k+1; d=fraction*taskDz(task);
        upperN=2;
        while twoGroupPower(d,upperN,upperN,alpha)<targetPower,upperN=upperN*2;end
        lowN=2;
        while lowN<upperN
            mid=floor((lowN+upperN)/2);
            if twoGroupPower(d,mid,mid,alpha)>=targetPower,upperN=mid;else,lowN=mid+1;end
        end
        requiredN=lowN;
        actualPower=twoGroupPower(d,nPatient,nControl,alpha);
        powerCheck=sampsizepwr('t2',[0 1],d,[],min(nPatient,nControl), ...
            'Ratio',max(nPatient,nControl)/min(nPatient,nControl),'Alpha',alpha,'Tail','both');
        assert(abs(actualPower-powerCheck)<1e-8,'Noncentral-t and sampsizepwr disagree.');
        rows(k,:)={taskNames{task},contrastNames{task},taskDz(task),fraction,d, ...
            requiredN,2*requiredN,twoGroupPower(d,requiredN,requiredN,alpha), ...
            twoGroupPower(d,requiredN-1,requiredN-1,alpha),nPatient,nControl, ...
            actualPower,detectableD,detectableD/taskDz(task)};
    end
end
PowerResults=cell2table(rows,'VariableNames',{'Task','Contrast','TaskDz', ...
    'Fraction','AssumedBetweenGroupD','RequiredNPerGroup','RequiredTotalN', ...
    'PowerAtRequiredN','PowerAtNMinus1','AvailablePatients','AvailableControls', ...
    'PowerAtAvailableN','DFor80PercentAtAvailableN','FractionFor80Percent'});
outputFolder=fullfile(repoRoot,'outputs','Power');
if ~isfolder(outputFolder),mkdir(outputFolder);end
writetable(PowerResults,fullfile(outputFolder,'Power_between_groups.csv'));
save(fullfile(outputFolder,'Power_between_groups.mat'),'PowerResults','alpha','targetPower','dataFile');
disp(PowerResults);
fprintf('\nN is PER GROUP for a balanced independent-groups t test, not total N.\n');
fprintf('PowerAtAvailableN uses the actual unequal groups. Benchmark fractions are assumptions.\n');
fprintf('Output: %s\n',outputFolder);

%% Local helpers; this file remains a directly runnable script
function p=twoGroupPower(d,n1,n2,alpha)
df=n1+n2-2; critical=tinv(1-alpha/2,df);
noncentrality=abs(d)/sqrt(1/n1+1/n2);
p=nctcdf(-critical,df,noncentrality)+nctcdf(critical,df,noncentrality,'upper');
assert(isfinite(p)&&p>=0&&p<=1,'Invalid noncentral-t power.');
end
function x=column(T,key)
names=regexprep(lower(string(T.Properties.VariableNames)),'[^a-z0-9]','');
idx=find(names==regexprep(lower(string(key)),'[^a-z0-9]',''));
assert(numel(idx)==1,'Missing or ambiguous column: %s',key);
v=T{:,idx};if isnumeric(v)||islogical(v),x=double(v);else,x=str2double(string(v));end
x=x(:);
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
