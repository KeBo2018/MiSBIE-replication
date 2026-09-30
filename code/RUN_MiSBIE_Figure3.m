%% MiSBIE: complete revised Figure 3 (formerly Figure 4)
% Open this file and click Run. This is a SCRIPT, not a function.
% MATLAB R2020b+; Statistics and Machine Learning Toolbox.
% All plotting helpers are local to this file. No previous workspace needed.
% Loads ../data/participants/analysis_data.csv automatically.
% Keep the repository directory structure intact.
%
% Revision 2026-09-29: original subgroup terminology; Bayes factors in text outputs only.
% Layout follows the original Figure 4 styling.
% Row 1: A, N-back NMDAS distribution; B, patient and control brain-behavior.
% Row 2: C, NMDAS versus N-back brain and behavior, patients only.
% D: Age-adjusted group -> NMDAS -> outcome path models (brain and behavior).
% E: Original resilient and severe subgroups, N-back brain and behavior.
% F: Multisensory and cold-pain severity correlations and subgroups.
% Supplement: other tasks' NMDAS distributions, overall group contrasts,
% and the pooled brain-behavior association. Within-group plots are in row 1.
%
% Provenance: Mito_Figure4_Forrevision.mlx / _CiRSquare.m for layout, colors,
% regressions, mediation models and task-specific k=2 severity clustering;
% SimpleCodeForMajorAnlysis_WithFDR.m for unified data and final WM eligibility.
% Known legacy errors are NOT propagated: mislabeled brain-behavior statistic,
% cluster labels based on row position, and stale multisensory severity values
% reused for cold pain. All displayed numbers are calculated from selected data.

%% Settings: normally no edits needed
dataFile = '';                 % Blank: use the release table; select package folder only if needed.
bootstrapSamples = 5000;       % Participant bootstrap for mediation CI.
randomSeed = 20260923;
bayesPriorScale = 1;           % JZS g ~ InvGamma(1/2, n*r^2/2).
useCANlabWhenAvailable = true; % Original inference via mediation.m on MATLAB path.
showControlContext = true;    % Matches original severity plots; tests are patients only.
fontName = 'Arial';

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder), scriptFolder = pwd; end
repoRoot=misbieFindPackageRoot(mfilename('fullpath'));
addpath(fullfile(repoRoot,'code'),'-begin');
outputFolder = fullfile(repoRoot,'outputs','Figure3');
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
assert(exist('partialcorr','file')==2 && exist('ksdensity','file')==2, ...
    'Statistics and Machine Learning Toolbox is required.');
if isempty(dataFile), dataFile=fullfile(repoRoot,'data','participants','analysis_data.csv'); end
assert(isfile(dataFile),'Release analysis CSV is missing.');
fprintf('Input table: %s\n',dataFile);
T=readtable(dataFile,'VariableNamingRule','preserve');
Group=column(T,'GroupID'); Age=column(T,'Age'); NMDAS=column(T,'NMDAS');
Accuracy=column(T,'Nback_ACC'); BehaviorFlag=column(T,'Behavioral_Available');
Brain=[column(T,'Nback'),column(T,'Multisensory'),column(T,'Cold')];
assert(all(ismember(Group(isfinite(Group)),[0 1])),'Expected GroupID 0=MitoD, 1=control.');
assert(all(ismember(BehaviorFlag(isfinite(BehaviorFlag)),[0 1])),'Expected behavioral flag 0/1.');
assert(all(Accuracy(isfinite(Accuracy))>=0 & Accuracy(isfinite(Accuracy))<=1), ...
    'Expected accuracy proportions in [0,1], not percentages.');
BehOK=BehaviorFlag==1; Patient=Group==0; Control=Group==1;
thresholdDisagreement=isfinite(Brain(:,1)) & (BehOK~=(isfinite(Accuracy)&Accuracy>.40));
if any(thresholdDisagreement)
    warning('%d recorded behavioral flags differ from accuracy > .40. Flags retained.',sum(thresholdDisagreement));
end
rng(randomSeed,'twister');
taskNames={'Working memory','Multisensory','Cold pain'};
outcomeNames={'N-back brain','N-back behavior','Multisensory brain','Cold-pain brain'};
blue=[0 114 189]/255; orange=[217 83 25]/255;
mildColor=[255 153 51]/255; severeColor=[204 0 0]/255;
groupColors=[blue;mildColor;severeColor];
imaging=isfinite(Brain)&repmat(Patient|Control,1,3);
imaging(:,1)=imaging(:,1)&BehOK;
severity=imaging&repmat(Patient&isfinite(NMDAS)&isfinite(Age),1,3);
behaviorSeverity=severity(:,1)&isfinite(Accuracy);
Y=[Brain(:,1),100*Accuracy,Brain(:,2),Brain(:,3)];
outcomeMasks=[severity(:,1),behaviorSeverity,severity(:,2),severity(:,3)];
expected=[74 25;90 28;91 29];
for ii=1:3
    found=[sum(imaging(:,ii)),sum(severity(:,ii))];
    fprintf('%s: imaging n=%d, severity patients n=%d\n',taskNames{ii},found);
    if any(found~=expected(ii,:)), warning('%s sample differs from manuscript; inspect input version.',taskNames{ii}); end
end

%% Continuous associations (same Spearman/age adjustment as final analyses)
N=zeros(4,1); rho=zeros(4,1); pRaw=rho; rhoAge=rho; pAge=rho;
for ii=1:4
    m=outcomeMasks(:,ii); N(ii)=sum(m);
    assert(N(ii)>=5,'Too few observations for a severity analysis.');
    [rho(ii),pRaw(ii)]=corr(NMDAS(m),Y(m,ii),'Type','Spearman');
    [rhoAge(ii),pAge(ii)]=partialcorr(NMDAS(m),Y(m,ii),Age(m),'Type','Spearman');
end
assert(all(isfinite([rho;pRaw;rhoAge;pAge])),'Undefined severity correlation.');
qRaw=bh(pRaw); qAge=bh(pAge); Outcome=outcomeNames';
severityStatistics=table(Outcome,N,rho,pRaw,qRaw,rhoAge,pAge,qAge);
writetable(severityStatistics,fullfile(outputFolder,'severity_correlations.csv'));

% The older master script plotted Brain vs Accuracy but annotated NMDAS vs
% Accuracy. Here both plotted variables enter the actual correlation.
bbMask=imaging(:,1)&isfinite(Accuracy)&isfinite(Age);
bbNames={'All participants';'MitoD only';'Controls only'};
bbMasks=[bbMask,bbMask&Patient,bbMask&Control];
bbN=zeros(3,1); bbR=zeros(3,1); bbP=bbR; bbAdjustedR=bbR; bbAdjustedP=bbR;
for ii=1:3
    m=bbMasks(:,ii); bbN(ii)=sum(m);
    assert(bbN(ii)>=5,'Too few observations for brain-behavior analysis.');
    [bbR(ii),bbP(ii)]=corr(Brain(m,1),Accuracy(m),'Type','Spearman');
    covariates=Age(m);
    if ii==1, covariates=[covariates Group(m)]; end
    [bbAdjustedR(ii),bbAdjustedP(ii)]=partialcorr(Brain(m,1),Accuracy(m),covariates,'Type','Spearman');
end
Adjustment={'Age and group';'Age';'Age'};
brainBehaviorStatistics=table(bbNames,bbN,bbR,bbP,bbAdjustedR,bbAdjustedP,Adjustment);
writetable(brainBehaviorStatistics,fullfile(outputFolder,'brain_behavior_correlations.csv'));

%% Age-adjusted rank-regression Bayes factors within each task
% Age is removed from BOTH ranked NMDAS and ranked brain response.
% Common intercept/age nuisance coefficients have matching flat priors.
% Tested slopes use a JZS mixture of g-priors; no BIC approximation is used.
% These are Gaussian models of ranks supporting the partial-Spearman tests.
BF10=nan(4,1);BF01=BF10;bfSensitivity=nan(3,3);priorScales=[.5 1 2];
for ii=[1 3 4]
    m=outcomeMasks(:,ii);[rr,nu]=ageResidualRanks([NMDAS(m),Y(m,ii)],Age(m));
    assert(abs(corr(rr(:,1),rr(:,2))-rhoAge(ii))<1e-10, ...
        'Bayes residualization does not match the age-adjusted correlation.');
    BF10(ii)=rankModelBF(rr(:,1),rr(:,2),nu,sum(m),bayesPriorScale);
    BF01(ii)=1/BF10(ii);
end
brainRows=[1 3 4];
for task=1:3
    m=severity(:,task);[rr,nu]=ageResidualRanks([NMDAS(m),Brain(m,task)],Age(m));
    for ss=1:3,bfSensitivity(task,ss)=rankModelBF(rr(:,1),rr(:,2),nu,sum(m),priorScales(ss));end
end
bayesTaskStatistics=table(taskNames',N(brainRows),rhoAge(brainRows),pAge(brainRows), ...
    BF10(brainRows),BF01(brainRows),bfSensitivity(:,1),bfSensitivity(:,2),bfSensitivity(:,3), ...
    'VariableNames',{'Task','N','RhoAge','PAge','BF10','BF01','BF10ScaleHalf','BF10ScaleOne','BF10ScaleTwo'});
writetable(bayesTaskStatistics,fullfile(outputFolder,'bayes_age_adjusted_tasks.csv'));
fprintf('Age-adjusted BF10 (WM, multisensory, cold): %.4g, %.4g, %.4g\n',BF10(brainRows));


%% Resilient and severe k=2 NMDAS subgroups: re-fit within each task as in Figure 4 code
% Labels depend on ordered centers, never participant row positions.
clusterID=nan(height(T),3); centers=nan(2,3);
for ii=1:3
    m=severity(:,ii);
    [ids,c]=kmeans_Matlab(NMDAS(m),2,'Distance','sqeuclidean','Replicates',100,'MaxIter',1000);
    [c,order]=sort(c); map=zeros(2,1); map(order)=[1;2];
    clusterID(m,ii)=map(ids); centers(:,ii)=c;
end
% Audit exact memberships and sample counts. Cluster 1 = resilient; cluster 2 = severe.
ReleaseID=string(T.ReleaseID);
rowAudit=table(ReleaseID,Group,Age,NMDAS,Accuracy,BehOK, ...
    imaging(:,1),imaging(:,2),imaging(:,3),severity(:,1),severity(:,2),severity(:,3), ...
    bbMask,clusterID(:,1),clusterID(:,2),clusterID(:,3), ...
    'VariableNames',{'ReleaseID','GroupID','Age','NMDAS','Accuracy','BehaviorIncluded', ...
    'WMImaging','MultiImaging','ColdImaging','WMSeverity','MultiSeverity','ColdSeverity', ...
    'BrainBehavior','WMCluster','MultiCluster','ColdCluster'});
writetable(rowAudit,fullfile(outputFolder,'analysis_row_audit.csv'));
Task=taskNames'; NImaging=sum(imaging,1)'; NSeverity=sum(severity,1)';
ResilientCenter=centers(1,:)'; SevereCenter=centers(2,:)';
sampleCounts=table(Task,NImaging,NSeverity,ResilientCenter,SevereCenter);
writetable(sampleCounts,fullfile(outputFolder,'task_sample_counts.csv'));
subgroupSets=cell(4,3); taskForOutcome=[1 1 2 3];
for ii=1:4
    k=taskForOutcome(ii);
    subgroupSets{ii,1}=Y(imaging(:,k)&Control&isfinite(Y(:,ii)),ii);
    subgroupSets{ii,2}=Y(clusterID(:,k)==1&isfinite(Y(:,ii)),ii);
    subgroupSets{ii,3}=Y(clusterID(:,k)==2&isfinite(Y(:,ii)),ii);
end
subgroupStatistics=groupTests(subgroupSets,outcomeNames);
writetable(subgroupStatistics,fullfile(outputFolder,'subgroup_tests.csv'));

%% Original mediation model, with explicit coding and complete-case sample
% Control=1, MitoD=0 reproduces the negative a-path orientation in the deck.
% Both outcomes use the SAME complete cases. Accuracy remains a proportion
% in path models, as in the old script; displayed scatter/violins use percent.
medMask=imaging(:,1)&isfinite(Accuracy)&isfinite(NMDAS)&isfinite(Age);
fprintf('Mediation complete cases: n=%d (%d patients, %d controls).\n', ...
    sum(medMask),sum(medMask&Patient),sum(medMask&Control));
medBrain=pathModel(Group(medMask),Brain(medMask,1),NMDAS(medMask),Age(medMask),bootstrapSamples);
medBehavior=pathModel(Group(medMask),Accuracy(medMask),NMDAS(medMask),Age(medMask),bootstrapSamples);
medBrain.Outcome='N-back brain'; medBehavior.Outcome='N-back behavior (proportion)';
% Use the original CANlab inference in the figure when the toolbox is present.
% Otherwise show explicitly labeled percentile CIs, with no invented p/stars.
medBrain.canlabP=nan(1,5);medBehavior.canlabP=nan(1,5);
medBrain.canlabSE=nan(1,5);medBehavior.canlabSE=nan(1,5);
CANlab=struct('available',exist('mediation','file')==2,'completed',false);
if useCANlabWhenAvailable && CANlab.available
    [CANlab.brainPaths,CANlab.brainStats]=mediation(Group(medMask),Brain(medMask,1), ...
        NMDAS(medMask),'boottop','covs',Age(medMask),'bootsamples',bootstrapSamples);
    [CANlab.behaviorPaths,CANlab.behaviorStats]=mediation(Group(medMask),Accuracy(medMask), ...
        NMDAS(medMask),'boottop','covs',Age(medMask),'bootsamples',bootstrapSamples);
    assert(max(abs(CANlab.brainPaths(:)-medBrain.paths(:)))<1e-7, ...
        'CANlab and OLS path estimates differ; inspect model implementation.');
    assert(max(abs(CANlab.behaviorPaths(:)-medBehavior.paths(:)))<1e-7, ...
        'CANlab and OLS behavior path estimates differ; inspect model implementation.');
    medBrain.canlabP=CANlab.brainStats.p(1,1:5);
    medBehavior.canlabP=CANlab.behaviorStats.p(1,1:5);
    medBrain.canlabSE=CANlab.brainStats.ste(1,1:5);
    medBehavior.canlabSE=CANlab.behaviorStats.ste(1,1:5);
    CANlab.completed=true;
elseif useCANlabWhenAvailable
    fprintf('CANlab mediation.m not on path: using self-contained OLS + percentile bootstrap; no CANlab p values claimed.\n');
end
mediationStatistics=[medTable(medBrain);medTable(medBehavior)];
writetable(mediationStatistics,fullfile(outputFolder,'mediation_paths.csv'));
if CANlab.completed
    mediationDescription=sprintf(['Path significance uses the original CANlab mediation routine, ' ...
        'with age covariates and bias-corrected bootstrap tests (initial resamples = %d). ' ...
        'Stars on paths indicate unadjusted bootstrap p values.'],bootstrapSamples);
else
    mediationDescription=sprintf(['Indirect-effect brackets are percentile 95%% confidence intervals ' ...
        'from %d participant bootstrap resamples. CANlab bias-corrected p values are unavailable; ' ...
        'no path significance stars are shown.'],bootstrapSamples);
end

%% Assemble the main figure: fixed centimeter positions keep rows separated
% A and B occupy the first row; C gets two wider axes on the second row.
% The original mediation, subgroup, and other-task results remain below.
fig=figure('Color','w','Name','MiSBIE revised Figure 3','Units','centimeters', ...
    'Position',[1 1 18 24],'Renderer','painters');
cols=[1.25 6.92 12.59]; colWidth=4.60;
panelHeading(fig,'A','NMDAS distribution',[.25 23.05 5.9 .5],fontName);
panelHeading(fig,'B','N-back brain and behavior',[6.02 23.05 11.3 .5],fontName);
ax=figureAxes(fig,[cols(1) 19.20 colWidth 3.35]);
distribution(ax,NMDAS(severity(:,1)),clusterID(severity(:,1),1), ...
    'N-back patients',mildColor,severeColor,fontName);
bbTitles={'Patients','Controls'};
bbXLim=[floor(min(Brain(bbMask,1)))-.2 ceil(max(Brain(bbMask,1)))+.2];
% Reserve an empty band above 100% for the two-line statistical inset.
bbYLim=[min(40,floor(min(100*Accuracy(bbMask))/10)*10) 116];
for ii=2:3
    m=bbMasks(:,ii); ax=figureAxes(fig,[cols(ii) 19.20 colWidth 3.35]);
    brainBehaviorPanel(ax,Brain(m,1),100*Accuracy(m),Group(m), ...
        bbR(ii),bbP(ii),bbAdjustedR(ii),bbAdjustedP(ii),blue,orange,fontName, ...
        bbTitles{ii-1});
    xlim(ax,bbXLim);ylim(ax,bbYLim);set(ax,'YTick',40:20:100);
end

panelHeading(fig,'C','Disease severity and N-back responses',[.25 17.55 17.1 .5],fontName);
severityX=[1.25 10.15]; severityTitles={'Brain pattern expression','Behavioral performance'};
for ii=1:2
    ax=figureAxes(fig,[severityX(ii) 13.40 7.1 3.55]);
    severityPanel(ax,NMDAS,Y(:,ii),outcomeMasks(:,ii), ...
        imaging(:,1)&Control&isfinite(NMDAS)&isfinite(Y(:,ii)),showControlContext, ...
        severityTitles{ii},rho(ii),pRaw(ii),rhoAge(ii),pAge(ii),orange,blue,fontName,BF10(ii));
    if ii==2,ylim(ax,[40 102]);set(ax,'YTick',40:20:100);end
end

panelHeading(fig,'D','Severity as a mediator',[.25 12.00 6.1 .5],fontName);
panelHeading(fig,'E','Resilient patients and Severe Patients',[6.40 12.00 11.0 .5],fontName);
ax=figureAxes(fig,[.35 7.12 5.90 4.58]);
mediationPanel(ax,medBrain,medBehavior,fontName);
subgroupX=[7.25 12.90]; subgroupTitles={'Brain pattern expression','Behavioral performance'};
for ii=1:2
    ax=figureAxes(fig,[subgroupX(ii) 7.65 4.40 3.75]);
    subgroupPanel(ax,subgroupSets(ii,:),subgroupTitles{ii},groupColors, ...
        subgroupStatistics.p((ii-1)*3+(1:3)),fontName);
end

panelHeading(fig,'F','Multisensory',[.25 5.98 8.5 .5],fontName);
panelHeading(fig,'','Cold pain',[9.05 5.98 8.5 .5],fontName);
bottomX=[1.10 5.60 10.10 14.45]; bottomWidth=[3.25 3.25 3.10 3.10];
for ii=1:2
    k=ii+2; task=ii+1;
    ax=figureAxes(fig,[bottomX(2*ii-1) 2.00 bottomWidth(2*ii-1) 3.45]);
    severityPanel(ax,NMDAS,Y(:,k),outcomeMasks(:,k), ...
        imaging(:,task)&Control&isfinite(NMDAS)&isfinite(Y(:,k)),showControlContext, ...
        taskNames{task},rho(k),pRaw(k),rhoAge(k),pAge(k),orange,blue,fontName,BF10(k));
    title(ax,'');style(ax,fontName,7);xlabel(ax,{'NMDAS','(Disease severity)'});
    ax=figureAxes(fig,[bottomX(2*ii) 2.00 bottomWidth(2*ii) 3.45]);
    subgroupPanel(ax,subgroupSets(k,:),taskNames{task},groupColors, ...
        subgroupStatistics.p((k-1)*3+(1:3)),fontName);
    title(ax,'');style(ax,fontName,7);
end
saveFigure(fig,outputFolder,'Figure3_complete');

%% Supplement: NMDAS distributions for the other two task samples
distFig=figure('Color','w','Name','Supplementary NMDAS distributions','Units','centimeters', ...
    'Position',[2 2 18 8],'Renderer','painters');
distLetters='AB';
for ii=1:2
    k=ii+1; ax=figureAxes(distFig,[1.35+(ii-1)*8.6 1.45 6.7 5.1]);
    distribution(ax,NMDAS(severity(:,k)),clusterID(severity(:,k),k), ...
        taskNames{k},mildColor,severeColor,fontName);
    panelHeading(distFig,distLetters(ii),'',[.35+(ii-1)*8.6 7.15 .7 .5],fontName);
end
saveFigure(distFig,outputFolder,'FigureS4_NMDAS');

%% Supplement: original overall group comparisons remain available
groupFig=figure('Color','w','Name','MiSBIE Figure S8','Units','centimeters', ...
    'Position',[2 2 18 15],'Renderer','painters');
overallRows=cell(4,7); order=[1 3 4 2];
panelLetters='ABCD';
for ii=1:4
    o=order(ii); k=taskForOutcome(o); m=imaging(:,k)&isfinite(Y(:,o));
    a=Y(m&Patient,o); b=Y(m&Control,o);
    [~,pv,ci,st]=ttest2(a,b,'Vartype','equal');
    overallRows(ii,:)={outcomeNames{o},numel(a),numel(b),st.tstat,st.df,pv,ci(:)'};
    ax=subplot(2,2,ii,'Parent',groupFig); hold(ax,'on');
    violin(ax,b,1,blue,[]); violin(ax,a,2,orange,[]);
    yline(ax,0,':','Color',[.6 .6 .6]);
    set(ax,'XTick',[1 2],'XTickLabel',{'Control','MitoD'});xlim(ax,[.5 2.5]);
    ylabel(ax,outcomeLabel(o));
    title(ax,{sprintf('%s  %s',panelLetters(ii),outcomeNames{o}), ...
        sprintf('t(%d) = %.2f, p %s',st.df,st.tstat,pPrefix(pv))},'FontWeight','normal');
    style(ax,fontName);
end
overallStatistics=cell2table(overallRows,'VariableNames', ...
    {'Outcome','NPatient','NControl','tPatientMinusControl','df','p','CI95'});
writetable(overallStatistics,fullfile(outputFolder,'overall_group_tests.csv'));
saveFigure(groupFig,outputFolder,'FigureS7_group_comparisons');

% Keep separate exports of row 1's within-group panels for PowerPoint reuse.
bbFig=figure('Color','w','Name','Brain-behavior within groups','Units','centimeters', ...
    'Position',[2 2 18 8],'Renderer','painters');
for ii=2:3
    m=bbMasks(:,ii); ax=subplot(1,2,ii-1,'Parent',bbFig);
    brainBehaviorPanel(ax,Brain(m,1),100*Accuracy(m),Group(m),bbR(ii),bbP(ii), ...
        bbAdjustedR(ii),bbAdjustedP(ii),blue,orange,fontName, ...
        bbTitles{ii-1});
    xlim(ax,bbXLim);ylim(ax,bbYLim);set(ax,'YTick',40:20:100);
end
saveFigure(bbFig,outputFolder,'Brain_behavior_within_groups');

pooledFig=figure('Color','w','Name','Supplementary pooled brain-behavior','Units','centimeters', ...
    'Position',[2 2 10 8],'Renderer','painters');
ax=figureAxes(pooledFig,[1.5 1.65 7.8 5.1]);
brainBehaviorPanel(ax,Brain(bbMask,1),100*Accuracy(bbMask),Group(bbMask),bbR(1),bbP(1), ...
    bbAdjustedR(1),bbAdjustedP(1),blue,orange,fontName,'All participants');
xlim(ax,bbXLim);ylim(ax,bbYLim);set(ax,'YTick',40:20:100);
saveFigure(pooledFig,outputFolder,'Supplement_brain_behavior_pooled');

%% Save all estimates and an automatically generated caption
results=struct('inputFile',dataFile,'randomSeed',randomSeed,'bootstrapSamples',bootstrapSamples, ...
    'severityStatistics',severityStatistics,'brainBehaviorStatistics',brainBehaviorStatistics, ...
    'subgroupStatistics',subgroupStatistics,'mediationStatistics',mediationStatistics, ...
    'overallStatistics',overallStatistics,'sampleCounts',sampleCounts,'rowAudit',rowAudit, ...
    'medBrain',medBrain,'medBehavior',medBehavior,'CANlab',CANlab, ...
    'bayesTaskStatistics',bayesTaskStatistics,'bayesPriorScale',bayesPriorScale);
save(fullfile(outputFolder,'Figure3_results.mat'),'results');
fid=fopen(fullfile(outputFolder,'Figure3_caption.txt'),'w');
assert(fid>=0,'Cannot write caption file.');
fprintf(fid,['Figure 3. The influence of mitochondrial disease severity on brain activation and behavioral performance.\n' ...
    '(A) NMDAS distribution in N-back patients. Amber and red indicate the resilient and severe patient subgroups identified using k-means clustering. ' ...
    '(B) Correlations between N-back brain pattern expression and task performance in patients and controls. ' ...
    '(C) Correlations between NMDAS and N-back brain pattern expression (left) and task performance (right) within patients. ' ...
    'In B, C, and F, age-adjusted partial Spearman correlations and p values are shown first, with unadjusted values in parentheses. ' ...
    'Points show individual observations; lines show linear fits. Controls on NMDAS axes are included for visual context only. ' ...
    '(D) Mediation analyses of disease severity, MitoD status, and N-back brain pattern expression or task performance. ' ...
    'Models control for age; control=1 and MitoD=0. Path coefficients are unstandardized; accuracy is expressed as a proportion. ' ...
    'Parentheses indicate bootstrap standard errors. %s ' ...
    '(E) N-back brain pattern expression and performance in control, resilient, and severe groups. ' ...
    '(F) Corresponding multisensory and cold-pain analyses. Dots represent participants and horizontal bars indicate means. ' ...
    'Blue, controls; orange, MitoD; amber, resilient patients; red, severe patients. ' ...
    '* p<.05, ** p<.01, *** p<.001 (two-sided t tests for group comparisons; bootstrap inference for mediation).\n'],mediationDescription);
fclose(fid);
% Statistical checks and Bayes factors are kept separate from the figure caption.
fid=fopen(fullfile(outputFolder,'Figure3_analysis_checks.txt'),'w');
assert(fid>=0,'Cannot write analysis checks.');
fprintf(fid,'Panel B samples: patients n=%d, controls n=%d. Panel D complete cases: n=%d.\n',bbN(2),bbN(3),sum(medMask));
fprintf(fid,'JZS prior scale r=%.3g.\n',bayesPriorScale);
for ii=1:3
    fprintf(fid,'%s: n=%d; age-adjusted rho=%.6f; BF10=%.6g.\n',taskNames{ii},N(brainRows(ii)),rhoAge(brainRows(ii)),BF10(brainRows(ii)));
end
fprintf(fid,'Supplementary pooled brain-behavior: n=%d, raw rho=%.4f, age/group-adjusted rho=%.4f, p=%.6g.\n', ...
    bbN(1),bbR(1),bbAdjustedR(1),bbAdjustedP(1));
fprintf(fid,'Patient-only brain-behavior: n=%d, rho=%.4f, age-adjusted rho=%.4f, p=%.6g.\n', ...
    bbN(2),bbR(2),bbAdjustedR(2),bbAdjustedP(2));
fprintf(fid,'Control-only brain-behavior: n=%d, rho=%.4f, age-adjusted rho=%.4f, p=%.6g.\n', ...
    bbN(3),bbR(3),bbAdjustedR(3),bbAdjustedP(3));
fprintf(fid,['Mediation n=%d. Missing control NMDAS values are excluded, never replaced with zero. ' ...
    'The existing manuscript/deck coefficients may use another input version; reconcile before replacing text.\n'],sum(medMask));
fprintf(fid,['Raw WM severity rho=%.4f (the earlier manuscript reports -.68). ' ...
    'CANlab inference used in panel D=%d. Percentile bootstrap CIs in the CSV are not CANlab bias-corrected p values.\n'],rho(1),CANlab.completed);
fclose(fid);
disp(brainBehaviorStatistics); disp(severityStatistics); disp(mediationStatistics);
fprintf('\nSaved the complete Figure 3 and supplementary figures (FIG, PNG, PDF), CSV statistics, and caption to:\n%s\n',outputFolder);

%% Local helpers (Run executes the script above; no function call is needed)
function x=column(T,key)
names=regexprep(lower(string(T.Properties.VariableNames)),'[^a-z0-9]','');
key=regexprep(lower(string(key)),'[^a-z0-9]',''); idx=find(names==key);
assert(numel(idx)==1,'Missing or ambiguous required column: %s',char(key));
raw=T{:,idx}; if isnumeric(raw)||islogical(raw),x=double(raw);else,x=str2double(string(raw));end
x=x(:);
end
function q=bh(p)
[s,idx]=sort(p(:)); a=flipud(cummin(flipud(s.*numel(s)./(1:numel(s))')));
q=nan(size(p));q(idx)=min(1,a);
end
function [rr,nu]=ageResidualRanks(values,age)
Z=[ones(numel(age),1),tiedrank(age(:))];
assert(rank(Z)==2,'Age nuisance design is rank deficient.');
R=tiedrank(values);rr=R-Z*(Z\R);nu=size(values,1)-rank(Z);
end
function bf=rankModelBF(x,y,nu,n,priorScale)
% Univariate Gaussian rank regression; g ~ InvGamma(1/2,n*r^2/2).
% Intercept/age nuisance priors match under the two hypotheses.
x=x(:);y=y(:);
assert(numel(x)==n && numel(y)==n && nu>1,'Invalid Bayes analysis sample.');
R2=(x'*y)^2/((x'*x)*(y'*y));
assert(isfinite(R2)&&R2>=0&&R2<1,'Invalid rank-regression R squared.');
a=n*priorScale^2;
% u=sqrt(n*r^2/g) has a half-normal density; tail above 12 is negligible.
fun=@(u) sqrt(2/pi).*exp(-u.^2/2 + .5.*(log(u.^2)-log(u.^2+a)) ...
    + nu/2.*(log(u.^2+a)-log(u.^2+a*(1-R2))));
bf=integral(fun,0,12,'RelTol',1e-9,'AbsTol',1e-10);
assert(isfinite(bf)&&bf>0,'Bayes-factor quadrature failed.');
end
function s=pText(p)
if p<.001,s='<.001';else,s=sprintf('%.3f',p);end
end
function ax=figureAxes(fig,box)
ax=axes('Parent',fig,'Units','centimeters','Position',box,'PositionConstraint','innerposition');
end
function style(ax,fontName,fontSize)
if nargin<3,fontSize=7.5;end
set(ax,'FontName',fontName,'FontSize',fontSize,'LineWidth',1,'TickDir','out', ...
    'TickLength',[.018 .018],'Box','off','XColor','k','YColor','k', ...
    'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1.12,'TitleFontWeight','bold');
end
function panelHeading(fig,letter,heading,box,fontName)
% Coordinates are centimeters; panel letters are separate from the heading.
if ~isempty(letter)
    annotation(fig,'textbox','Units','centimeters','Position',[box(1) box(2)-.04 .6 box(4)+.12], ...
        'String',letter,'EdgeColor','none','FontName',fontName,'FontSize',13, ...
        'FontWeight','bold','Margin',0);
end
annotation(fig,'textbox','Units','centimeters','Position',[box(1)+.68 box(2) box(3)-.68 box(4)], ...
    'String',heading,'EdgeColor','none','FontName',fontName,'FontSize',9, ...
    'FontWeight','bold','Margin',0,'FitBoxToText','off','HorizontalAlignment','center');
end
function distribution(ax,x,cluster,name,lowColor,highColor,fontName)
hold(ax,'on'); edges=-.5:5:max(44.5,ceil(max(x)/5)*5+4.5);
counts=[histcounts(x(cluster==1),edges);histcounts(x(cluster==2),edges)]';
h=bar(ax,edges(1:end-1)+2.5,counts,.95,'stacked','EdgeColor','w','LineWidth',.5);
h(1).FaceColor=lowColor;h(2).FaceColor=highColor;
[v,idx]=sort(x); c=cluster(idx); rgb=repmat(lowColor,numel(x),1);rgb(c==2,:)=repmat(highColor,sum(c==2),1);
scatter(ax,v,-.22-.15*mod((1:numel(x))',3),13,rgb,'filled');
ylim(ax,[-.8 max(1,max(sum(counts,2))*1.12)]);xlim(ax,[-1 max(41,max(x)+1)]);
set(ax,'XTick',0:10:max(40,max(x)));xlabel(ax,{'NMDAS','(Disease severity)'});ylabel(ax,'Number of patients');
title(ax,name,'FontWeight','bold','FontSize',8.5);
style(ax,fontName);
end
function fittedLine(ax,x,y,color)
if numel(unique(x))<2,return;end
c=polyfit(x,y,1);xx=[min(x) max(x)];plot(ax,xx,polyval(c,xx),'-','Color',color,'LineWidth',1.3);
end
function brainBehaviorPanel(ax,x,y,g,r,p,ra,pa,blue,orange,fontName,name)
hold(ax,'on');scatter(ax,x(g==1),y(g==1),15,blue,'filled');scatter(ax,x(g==0),y(g==0),15,orange,'filled');
lineColor=[1 89/255 76/255];
if all(g==1),lineColor=blue;elseif any(g==1)&&any(g==0),lineColor=[.25 .25 .25];end
fittedLine(ax,x,y,lineColor);
xlabel(ax,{'Brain pattern expression','(N-back)'});ylabel(ax,'N-back performance (%)');
title(ax,name,'FontWeight','bold','FontSize',8.5);style(ax,fontName);
adjustment='age';if any(g==0)&&any(g==1),adjustment='adj';end
correlationNote(ax,r,p,ra,pa,'left',fontName,nan,adjustment);
end
function correlationNote(ax,r,p,ra,pa,side,fontName,bf,adjustment)
% Age-adjusted rho/p first; unadjusted rho/p in parentheses.
if nargin<8,bf=nan;end
if nargin<9,adjustment='age';end
pos=.97;if strcmp(side,'left'),pos=.04;end
oldUnits=ax.Units;ax.Units='centimeters';axisWidth=ax.Position(3);ax.Units=oldUnits;
fontSize=7; if axisWidth<3.5,fontSize=6;end
label={sprintf('\\rho = %.2f (%.2f)',ra,r),sprintf('p %s (%s)',pPrefix(pa),pText(p))};
% Bayes factors are exported for the manuscript text, not displayed in Figure 3.
text(ax,pos,.97,label,'Units','normalized', ...
    'HorizontalAlignment',side,'VerticalAlignment','top','FontName',fontName, ...
    'FontSize',fontSize,'Interpreter','tex','BackgroundColor','w','Margin',1);
end
function s=pPrefix(p)
if p<.001,s='< .001';else,s=['= ' sprintf('%.3f',p)];end
end
function severityPanel(ax,x,y,m,controls,showControls,name,r,p,ra,pa,orange,blue,fontName,bf)
hold(ax,'on');
if showControls,scatter(ax,x(controls),y(controls),13,blue,'filled');end
scatter(ax,x(m),y(m),15,orange,'filled');fittedLine(ax,x(m),y(m),[1 89/255 76/255]);
xlim(ax,[-1 max(41,max(x(m))+1)]);set(ax,'XTick',0:10:max(40,max(x(m))));
xlabel(ax,'NMDAS (Disease severity)');
if contains(lower(name),'behavior'),ylabel(ax,'N-back performance (%)');else,ylabel(ax,'Brain pattern expression');end
title(ax,name,'FontWeight','bold','FontSize',8.5);style(ax,fontName);
oldUnits=ax.Units;ax.Units='centimeters';axisWidth=ax.Position(3);ax.Units=oldUnits;
if axisWidth<3.5
    visible=y(m);if showControls,visible=[visible;y(controls)];end
    span=max(eps,max(visible)-min(visible));
    ylim(ax,[min(visible)-.06*span max(visible)+.5*span]);
end
correlationNote(ax,r,p,ra,pa,'right',fontName,bf);
end
function violin(ax,v,x,color,bw)
v=v(isfinite(v));hold(ax,'on');
if numel(unique(v))>1
    if isempty(bw),[d,yy]=ksdensity(v);else,[d,yy]=ksdensity(v,'Bandwidth',bw);end
    width=.29*d/max(d);patch(ax,[x-width fliplr(x+width)],[yy fliplr(yy)],color, ...
        'FaceAlpha',.60,'EdgeColor',[.1 .1 .1],'LineWidth',.8);
end
scatter(ax,x+.15*sin((1:numel(v))'*2.3999632297),v,12,'k','filled');
plot(ax,[x-.21 x+.21],[mean(v) mean(v)],'k-','LineWidth',1.6);
end
function subgroupPanel(ax,sets,name,colors,p,fontName)
bw=.4;if contains(lower(name),'behavior'),bw=5;elseif contains(lower(name),'multi'),bw=.8;end
hold(ax,'on');for j=1:3,violin(ax,sets{j},j,colors(j,:),bw);end
% Place brackets using observed values, not automatic density-tail limits.
% This avoids the large empty accuracy range in the earlier assembled figure.
v=vertcat(sets{:});v=v(isfinite(v));span=max(eps,max(v)-min(v));
pairs=[1 2;1 3;2 3];y=max(v)+.12*span;n=0;
for j=1:3
    if p(j)<.05
        n=n+1; yy=y+(n-1)*.10*span;xx=pairs(j,:);
        plot(ax,[xx(1) xx(1) xx(2) xx(2)],[yy-.025*span yy yy yy-.025*span],'k-','LineWidth',1.1);
        stars='*';if p(j)<.001,stars='***';elseif p(j)<.01,stars='**';end
        text(ax,mean(xx),yy+.012*span,stars,'HorizontalAlignment','center','FontSize',8,'FontName',fontName);
    end
end
ylim(ax,[min(v)-.18*span y+(max(0,n-1)*.10+.12)*span]);xlim(ax,[.5 3.5]);
set(ax,'XTick',1:3,'XTickLabel',{'Control',sprintf('Resilient\npatients'),sprintf('Severe\nPatients')},'XTickLabelRotation',0);
xlabel(ax,'Group');
if contains(lower(name),'behavior'),ylabel(ax,'N-back performance (%)');else,ylabel(ax,'Brain pattern expression');end
title(ax,name,'FontWeight','bold','FontSize',8.5);style(ax,fontName);
end
function out=groupTests(sets,names)
pairs=[1 2;1 3;2 3]; groupNames={'Control','Resilient patients','Severe Patients'};r=cell(12,9);i=0;
for o=1:4
    for j=1:3
        i=i+1;a=sets{o,pairs(j,1)};b=sets{o,pairs(j,2)};
        [~,p,ci,s]=ttest2(a,b,'Vartype','equal');
        r(i,:)={names{o},[groupNames{pairs(j,1)} ' minus ' groupNames{pairs(j,2)}],numel(a),numel(b),s.tstat,s.df,p,ci(1),ci(2)};
    end
end
out=cell2table(r,'VariableNames',{'Outcome','Contrast','N1','N2','t','df','p','CILow','CIHigh'});
out.qBH12=bh(out.p);
end
function paths=olsPaths(x,y,m,age)
D=[ones(numel(x),1),x,age];F=[D,m];
assert(rank(D)==3&&rank(F)==4,'Rank-deficient mediation design.');
ba=D\m;bc=D\y;bb=F\y;
paths=[ba(2),bb(4),bb(2),bc(2),ba(2)*bb(4)];
end
function s=pathModel(x,y,m,age,B)
s.paths=olsPaths(x,y,m,age);s.n=numel(x);draws=nan(B,5);b=0;attempts=0;
while b<B
    attempts=attempts+1;assert(attempts<5*B,'Too many singular bootstrap resamples.');
    id=randi(numel(x),numel(x),1);D=[ones(numel(id),1),x(id),age(id),m(id)];
    if rank(D)<4,continue;end
    b=b+1;draws(b,:)=olsPaths(x(id),y(id),m(id),age(id));
end
s.bootstrapSE=std(draws,0,1);s.ci=prctile(draws,[2.5 97.5],1);s.bootstrapSamples=B;
end
function t=medTable(s)
Outcome=repmat({s.Outcome},5,1);Path={'a';'b';'c_prime';'c';'a_times_b'};
N=repmat(s.n,5,1);Estimate=s.paths';BootstrapSE=s.bootstrapSE';
CILow=s.ci(1,:)';CIHigh=s.ci(2,:)';
CANlabP=s.canlabP';CANlabSE=s.canlabSE';
t=table(Outcome,Path,N,Estimate,BootstrapSE,CILow,CIHigh,CANlabP,CANlabSE);
end
function mediationPanel(ax,a,b,fontName)
axis(ax,[0 1 0 1]);axis(ax,'off');hold(ax,'on');
models={a,b};labels={{'Brain pattern','expression'},{'Behavioral','performance'}};
nodeColor=[8 44 59]/255;
for k=1:2
    z=.49*(2-k)+.055;s=models{k};
    quiver(ax,.265,z+.164,.143,.139,0,'k','MaxHeadSize',.4,'LineWidth',1.7);
    quiver(ax,.617,z+.303,.128,-.139,0,'k','MaxHeadSize',.4,'LineWidth',1.7);
    quiver(ax,.31,z+.11,.38,0,0,'k','MaxHeadSize',.25,'LineWidth',1.7);
    ovalNode(ax,[.16 z+.11],{'Diagnostic','group'},nodeColor,fontName);
    ovalNode(ax,[.51 z+.36],{'Disease','severity'},nodeColor,fontName);
    ovalNode(ax,[.85 z+.11],labels{k},nodeColor,fontName);
    text(ax,.01,z+.27,pathLabel(s,1,'a'),'FontSize',6,'FontName',fontName);
    text(ax,.99,z+.27,pathLabel(s,2,'b'),'HorizontalAlignment','right','FontSize',6,'FontName',fontName);
    text(ax,.51,z+.01,pathLabel(s,3,'c'''),'HorizontalAlignment','center','FontSize',6,'FontName',fontName);
    if all(isfinite(s.canlabP))
        indirectLabel=pathLabel(s,5,'ab');
    else
        indirectLabel=sprintf('ab %.3g [%.3g, %.3g]',s.paths(5),s.ci(1,5),s.ci(2,5));
    end
    text(ax,.51,z+.205,indirectLabel, ...
        'HorizontalAlignment','center','FontSize',6,'FontName',fontName);
end
end
function ovalNode(ax,center,label,color,fontName)
rectangle(ax,'Position',[center(1)-.15 center(2)-.08 .30 .16], ...
    'Curvature',[1 1],'EdgeColor',color,'FaceColor','w','LineWidth',1.7);
text(ax,center(1),center(2),label,'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontName',fontName,'FontSize',6);
end
function label=pathLabel(s,index,name)
if isfinite(s.canlabSE(index))
    label=sprintf('%s %.3g (%.2g)%s',name,s.paths(index),s.canlabSE(index),sigStars(s.canlabP(index)));
else
    label=sprintf('%s %.3g',name,s.paths(index));
end
end
function s=sigStars(p)
s='';if p<.001,s='***';elseif p<.01,s='**';elseif p<.05,s='*';end
end
function label=outcomeLabel(i)
if i==2,label='Accuracy (%)';else,label='Task contrast pattern expression';end
end
function saveFigure(fig,folder,name)
drawnow;savefig(fig,fullfile(folder,[name '.fig']));
exportgraphics(fig,fullfile(folder,[name '.png']),'Resolution',600,'BackgroundColor','white');
exportgraphics(fig,fullfile(folder,[name '.pdf']),'ContentType','vector','BackgroundColor','white');
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
