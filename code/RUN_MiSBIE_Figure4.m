%% GDF15 results -- manuscript sample definitions and covariates
% Panel A: original N-back imaging cohort, before behavioral eligibility.
% GDF15/NMDAS in panel A is Spearman in the original plotting source.
% Panel B: eligible N-back cohort; Pearson within groups; partial Pearson
% controlling Group for the combined cohort. Group coding is 0 patient/1 control.
% Main-effect and interaction regressions additionally include Age.
clearvars;
repoRoot=misbieFindPackageRoot(mfilename('fullpath'));
out=fullfile(repoRoot,'outputs','Figure4');if ~isfolder(out),mkdir(out);end
T=readtable(fullfile(repoRoot,'data','participants','analysis_data.csv'));
g=T.GroupID;x=log10(T.GDF15);brain=T.Nback;beh=T.Nback_ACC;age=T.Age;nm=T.NMDAS;
base=isfinite(brain)&ismember(g,[0 1])&isfinite(x);
eligible=base&T.Behavioral_Available==1&isfinite(beh);
a=base&g==0;b=base&g==1;
[~,p,ci,st]=ttest2(x(a),x(b));
groupStats=table(sum(a),sum(b),mean(x(a)),mean(x(b)),st.tstat,st.df,p,ci(1),ci(2), ...
    'VariableNames',{'NPatient','NControl','MeanPatient','MeanControl','T','DF','P','CILower','CIUpper'});
writetable(groupStats,fullfile(out,'GDF15_group_test.csv'));
rows=cell(8,7);j=0;
for method={'Spearman','Pearson'}
    m=a&isfinite(nm);[r,p]=corr(x(m),nm(m),'Type',method{1});j=j+1;
    rows(j,:)={'GDF15 and NMDAS','Patients before behavioral exclusions',method{1},'None',sum(m),r,p};
end
ys=[brain,100*beh];ynames={'N-back brain','N-back behavior'};
for yidx=1:2
    for groupIndex=1:3
        m=eligible&isfinite(ys(:,yidx));pop='All eligible';adjust='Group';
        if groupIndex==2,m=m&g==0;pop='Patients';adjust='None';end
        if groupIndex==3,m=m&g==1;pop='Controls';adjust='None';end
        if groupIndex==1,[r,p]=partialcorr(x(m),ys(m,yidx),g(m),'Type','Pearson');
        else,[r,p]=corr(x(m),ys(m,yidx),'Type','Pearson');end
        j=j+1;rows(j,:)={ynames{yidx},pop,'Pearson',adjust,sum(m),r,p};
    end
end
C=cell2table(rows,'VariableNames',{'Outcome','Population','Correlation','Covariates','N','R','P'});
writetable(C,fullfile(out,'GDF15_correlations.csv'));disp(groupStats);disp(C);
% Match mean-centering in the archived code; coefficients' tests are invariant
% to centering in Model 1. Model 2 uses the archived full factorial expansion.
for yidx=1:2
    m=eligible&isfinite(age)&isfinite(ys(:,yidx));
    A=table(x(m)-mean(x(eligible)),age(m)-mean(age(eligible)),2*g(m)-1,ys(m,yidx)-mean(ys(eligible,yidx)), ...
        'VariableNames',{'GDF15','Age','Group','Outcome'});
    mdl1=fitlm(A,'Outcome ~ GDF15 + Age + Group');mdl2=fitlm(A,'Outcome ~ GDF15*Age*Group');
    writetable(mdl1.Coefficients,fullfile(out,sprintf('Outcome%d_main_effects.csv',yidx)),'WriteRowNames',true);
    writetable(mdl2.Coefficients,fullfile(out,sprintf('Outcome%d_interactions.csv',yidx)),'WriteRowNames',true);
    save(fullfile(out,sprintf('Outcome%d_models.mat',yidx)),'mdl1','mdl2');
end
f=figure('Color','w','Position',[100 100 950 750]);blue=[0 114 189]/255;orange=[217 83 25]/255;
subplot(2,2,1);hold on;rng(42);
scatter(1+.1*randn(sum(b),1),x(b),22,blue,'filled');scatter(2+.1*randn(sum(a),1),x(a),22,orange,'filled');
plot([.8 1.2],mean(x(b))*[1 1],'k-','LineWidth',2);plot([1.8 2.2],mean(x(a))*[1 1],'k-','LineWidth',2);
set(gca,'XTick',[1 2],'XTickLabel',{'Controls','Patients'});xlim([.5 2.5]);ylabel('Plasma GDF15 (log_{10})');title('A  Plasma GDF15');box off;
subplot(2,2,2);m=a&isfinite(nm);scatter(x(m),nm(m),25,orange,'filled');hold on;fitline(x(m),nm(m),orange);
xlabel('Plasma GDF15 (log_{10})');ylabel('NMDAS');title(sprintf('Disease severity: rho = %.2f',C.R(1)));box off;
for yidx=1:2
    subplot(2,2,2+yidx);hold on;
    for group=0:1
        m=eligible&g==group;color=orange;if group==1,color=blue;end
        scatter(x(m),ys(m,yidx),25,color,'filled');fitline(x(m),ys(m,yidx),color);
    end
    xlabel('Plasma GDF15 (log_{10})');ylabel(ynames{yidx});title(ynames{yidx});box off;
end
exportgraphics(f,fullfile(out,'Figure4_recomputed.png'),'Resolution',300);
exportgraphics(f,fullfile(out,'Figure4_recomputed.pdf'),'ContentType','vector');savefig(f,fullfile(out,'Figure4_recomputed.fig'));
function fitline(x,y,color)
    z=polyfit(x,y,1);xx=linspace(min(x),max(x),100);plot(xx,polyval(z,xx),'Color',color,'LineWidth',1.5);
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
