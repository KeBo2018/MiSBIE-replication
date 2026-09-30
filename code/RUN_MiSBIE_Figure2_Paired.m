%% MiSBIE Figure 2: paired condition plots for all three tasks
% Open this SCRIPT and click Run. The complete sharing ZIP includes its data folder.
% Creates the full revised Figure 2, a row of the three paired plots, and
% separate editable FIG / vector PDF / 600-dpi PNG plot files.
% Requires base MATLAB R2020a or newer; no CANlab or Statistics toolbox.
%
% Actual held-out condition scores from the FINAL (100th) archived CV split
% are used consistently for all tasks. Multisensory/cold archives retain
% separate conditions for that split, but not for all 100 repetitions.
% No condition is reconstructed from a difference score.
% Paired Cohen's dz = mean(task-control)/std(task-control,0).
% Labels report ALL-PAIR dz BEFORE visualization filtering: 1.68, 2.18, 0.96.
% Cold pain only: omit a whole pair from display if either condition is more
% than 3 scaled MAD from its condition median (one pass). This post hoc
% display choice is NOT evidence of invalid data or an analytical exclusion.
% The unfiltered pain plot and full row-level audit are also exported.
% The original 100-split means (1.68, 2.18, 0.92) remain in the audit CSV.
% Neither between-task tests nor changes to Figure 3 are included.
%
% The PowerPoint version has editable vector plots and original source
% objects. The MATLAB full-figure PDF uses a high-resolution background
% for brain maps/text and vector axes/data for the new paired plots.

clearvars;
filterColdDisplay=true; % Set false to show every pair in the main figure.
thresholdFactor=3; % Standard median rule; do not tune to the effect size.
scriptPath=mfilename('fullpath');
if isempty(scriptPath), scriptDir=pwd; else, scriptDir=fileparts(scriptPath); end
repoRoot=misbieFindPackageRoot(mfilename('fullpath'));
S=load(fullfile(repoRoot,'data','figure2','Paired_condition_scores.mat'));
outDir=fullfile(repoRoot,'outputs','Figure2');
if ~isfolder(outDir), mkdir(outDir); end
keys={'Nback','Multisensory','Cold'};
names={'Working memory','Multisensory','Cold pain'};
labels={{'0-back','2-back'},{'Rest','Stimulus'},{sprintf('Room\ntemp.'),'Cold'}};
taskColors=[.35 .55 .72;1 89/255 76/255]; % control condition / task condition
paired=cell(3,1);N=zeros(3,1);D=N;meanCVD=N;
for k=1:3
    control=double(S.([keys{k} '_Control'])); control=control(:);
    task=double(S.([keys{k} '_Task'])); task=task(:);
    assert(numel(control)==numel(task) && all(isfinite([control;task])), ...
        'Incomplete paired scores for %s.',names{k});
    delta=task-control;
    reference=S.([keys{k} '_ReferenceDelta']);
    assert(max(abs(delta-reference(:)))<1e-10,'Mismatch with saved paired differences.');
    runs=S.([keys{k} '_DValues']);
    D(k)=mean(delta)/std(delta,0);
    assert(abs(D(k)-runs(end))<1e-10,'Effect size disagrees with archived final CV split.');
    meanCVD(k)=mean(runs);N(k)=numel(task);
    paired{k}=[control task];
end
assert(isequal(N,[88;90;91]),'Unexpected sample counts.');
% The raw pairs remain intact. The display mask applies only to cold pain.
displayPairs=paired;displayOmitted=cell(3,1);shownN=N;shownD=D;
for k=1:3, displayOmitted{k}=false(N(k),1); end
[coldFlags,lower,upper]=medianOutlierMask(paired{3},thresholdFactor);
if filterColdDisplay
    displayOmitted{3}=any(coldFlags,2);
    displayPairs{3}=paired{3}(~displayOmitted{3},:);
end
for k=1:3
    shownN(k)=size(displayPairs{k},1);
    diffShown=displayPairs{k}(:,2)-displayPairs{k}(:,1);
    shownD(k)=mean(diffShown)/std(diffShown,0);
end
Summary=table(string(names(:)),N,shownN,N-shownN,repmat(100,3,1),D,shownD,meanCVD, ...
    'VariableNames',{'Task','FullN','DisplayedN','OmittedFromDisplay', ...
    'CVRepetition','AnnotatedFullSampleDz','DisplayedPairsDz_Sensitivity','OriginalMeanCVDz'});
disp(Summary);writetable(Summary,fullfile(outDir,'Paired_plot_summary.csv'));
ColdAudit=table((1:N(3))',paired{3}(:,1),paired{3}(:,2), ...
    paired{3}(:,2)-paired{3}(:,1),coldFlags(:,1),coldFlags(:,2), ...
    any(coldFlags,2),displayOmitted{3},repmat(lower(1),N(3),1), ...
    repmat(upper(1),N(3),1),repmat(lower(2),N(3),1),repmat(upper(2),N(3),1), ...
    'VariableNames',{'SavedRowIndex','RoomTemperatureScore','ColdScore','PairedDifference', ...
    'RoomFlagged','ColdFlagged','FlaggedEitherCondition','OmittedFromDisplay', ...
    'RoomLower','RoomUpper','ColdLower','ColdUpper'});
writetable(ColdAudit,fullfile(outDir,'Cold_display_outlier_audit.csv'));
fprintf('Cold display: %d of %d pairs omitted; Cohen''s d retains all pairs.\n', ...
    sum(displayOmitted{3}),N(3));

%% Three panels at a comfortable review size
f=makeFigure('Paired task responses',18,6.6,'on');
for k=1:3
    ax=axes(f,'Position',[.075+(k-1)*.328 .23 .24 .58]);
    drawPaired(ax,displayPairs{k},labels{k},D(k),taskColors,8,true);
    text(ax,.5,1.25,names{k},'Units','normalized','FontName','Arial', ...
        'FontSize',9.5,'FontWeight','bold','HorizontalAlignment','center');
end
if any(displayOmitted{3})
    addDisclosure(f,[.72 .02 .27 .085],sum(displayOmitted{3}),6.5);
end
savePlot(f,outDir,'Three_paired_task_plots');

%% Individual panels (larger than their placement in Figure 2)
for k=1:3
    fk=makeFigure(names{k},7.4,6.5,'off');
    ax=axes(fk,'Position',[.20 .22 .75 .65]);
    drawPaired(ax,displayPairs{k},labels{k},D(k),taskColors,8,true);
    if any(displayOmitted{k})
        addDisclosure(fk,[.14 .01 .80 .09],sum(displayOmitted{k}),6.8);
    end
    savePlot(fk,outDir,['Paired_' keys{k}]);
    close(fk);
end

%% Unfiltered pain view for transparent inspection
fu=makeFigure('Cold pain - all pairs',7.4,6.5,'off');
ax=axes(fu,'Position',[.20 .22 .75 .65]);
drawPaired(ax,paired{3},labels{3},D(3),taskColors,8,true);
savePlot(fu,outDir,'Paired_Cold_ALL_PAIRS');close(fu);

%% Full Figure 2, using the revised source-artwork background
backgroundPath=fullfile(repoRoot,'data','figure_assets','Figure2_background.png');
assert(isfile(backgroundPath),'Keep the bundled Figure2_background.png in data/.');
canvasW=680;canvasH=710;
fullFig=makeFigure('Revised MiSBIE Figure 2',canvasW/96*2.54,canvasH/96*2.54,'on');
bg=axes(fullFig,'Position',[0 0 1 1]);
image(bg,imread(backgroundPath));axis(bg,'image');axis(bg,'off');
set(bg,'Position',[0 0 1 1],'XLimMode','auto','YLimMode','auto');
% Coordinates match the compact A/C PowerPoint; plots moved up 34 CSS pixels.
for k=1:3
    left=140+(k-1)*220; top=238; width=75; height=96;
    ax=axes(fullFig,'Position',[left/canvasW,1-(top+height)/canvasH, ...
        width/canvasW,height/canvasH]);
    drawPaired(ax,displayPairs{k},labels{k},D(k),taskColors,6.6,false);
end
if any(displayOmitted{3})
    addDisclosure(fullFig,[462/canvasW,1-(335+25)/canvasH,95/canvasW,25/canvasH], ...
        sum(displayOmitted{3}),6.15);
end
savePlot(fullFig,outDir,'Figure2_compact_AC');
save(fullfile(outDir,'Paired_plot_analysis.mat'),'Summary','paired','displayPairs','displayOmitted','ColdAudit','filterColdDisplay','thresholdFactor','labels','taskColors');
fprintf('\nDone. Figures and numerical results are in:\n%s\n',outDir);
fprintf('Cohen''s d uses ALL pairs in the final CV split, before display filtering.\n');

%% Local helpers (this file runs as a script)
function f=makeFigure(name,w,h,visibility)
f=figure('Name',name,'NumberTitle','off','Color','w','Units','centimeters', ...
    'Position',[2 2 w h],'Visible',visibility,'Renderer','painters');
set(f,'PaperUnits','centimeters','PaperSize',[w h],'PaperPosition',[0 0 w h]);
end

function drawPaired(ax,pairs,labels,dz,colors,fontSize,large)
hold(ax,'on');n=size(pairs,1);delta=pairs(:,2)-pairs(:,1);
[~,order]=sort(delta);jitter=zeros(n,1);
jitter(order)=.105*(2*mod((1:n)'*.618033988749895,1)-1);
for j=1:2
    [gridY,width]=violinDensity(pairs(:,j));
    patch(ax,[j-width fliplr(j+width)],[gridY fliplr(gridY)], ...
        .35*colors(j,:)+.65,'EdgeColor',colors(j,:),'LineWidth',.65);
end
if large, lineWidth=.4;dotSize=7; else,lineWidth=.24;dotSize=1.7; end
for i=1:n
    plot(ax,[1 2]+jitter(i),pairs(i,:),'-','Color',[.78 .78 .78],'LineWidth',lineWidth);
end
for j=1:2
    scatter(ax,j+jitter,pairs(:,j),dotSize,colors(j,:),'filled');
    m=mean(pairs(:,j));plot(ax,[j-.22 j+.22],[m m],'k-','LineWidth',1.0);
end
allValues=pairs(:);tickMin=2*floor(min(allValues)/2);tickMax=2*ceil(max(allValues)/2);
if tickMax-tickMin<=14,tickStep=2;else,tickStep=4;end
set(ax,'XLim',[.5 2.5],'YLim',[tickMin-.6 tickMax+.6], ...
    'XTick',[1 2],'XTickLabel',labels,'YTick',tickStep*ceil(tickMin/tickStep):tickStep:tickMax, ...
    'Box','off','TickDir','out','LineWidth',.6,'FontName','Arial', ...
    'FontSize',fontSize,'XColor',[.15 .15 .15],'YColor',[.15 .15 .15]);
ylabel(ax,'Pattern expression','FontSize',fontSize);
title(ax,sprintf('Cohen''s d = %.2f',dz),'FontName','Arial', ...
    'FontSize',fontSize+.35,'FontWeight','bold');
end

function [grid,width]=violinDensity(x)
x=x(:);sorted=sort(x);n=numel(x);
q=zeros(1,2);p=[.25 .75];
for j=1:2
    h=1+(n-1)*p(j);a=floor(h);b=ceil(h);
    q(j)=sorted(a)+(h-a)*(sorted(b)-sorted(a));
end
scale=min(std(x,0),(q(2)-q(1))/1.349);
if ~(scale>0),scale=std(x,0);end
bandwidth=.9*scale*n^(-.2);
assert(bandwidth>0,'Degenerate condition-score distribution.');
grid=linspace(min(x),max(x),160);
density=mean(exp(-.5*((grid-x)/bandwidth).^2),1)/(bandwidth*sqrt(2*pi));
width=.22*density/max(density);width([1 end])=0;
end

function savePlot(f,folder,name)
drawnow;savefig(f,fullfile(folder,[name '.fig']));
% print preserves full-figure dimensions and annotation positions.
print(f,fullfile(folder,[name '.png']),'-dpng','-r600');
print(f,fullfile(folder,[name '.pdf']),'-dpdf','-painters');
end

function [flags,lower,upper]=medianOutlierMask(pairs,thresholdFactor)
% Same finite-data criterion as isoutlier(pairs,'median','ThresholdFactor',3).
% Evaluate each condition once, then combine with any(flags,2) to retain pairing.
center=median(pairs,1);
scaledMAD=1.482602218505602*median(abs(pairs-center),1);
assert(all(scaledMAD>0),'MAD is zero; review the data rather than auto-filtering.');
lower=center-thresholdFactor*scaledMAD;upper=center+thresholdFactor*scaledMAD;
flags=pairs<lower | pairs>upper;
end

function addDisclosure(f,position,nOmitted,fontSize)
annotation(f,'textbox',position,'String', ...
    sprintf('%d pairs omitted;\nd uses all pairs',nOmitted), ...
    'FontName','Arial','FontSize',fontSize,'Color',[.27 .27 .27], ...
    'HorizontalAlignment','left','VerticalAlignment','middle', ...
    'EdgeColor','none','Margin',0,'Interpreter','none');
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
