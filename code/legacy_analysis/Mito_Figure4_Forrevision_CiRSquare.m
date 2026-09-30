% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
%% 
% Misbie brain project 
% 
% 
%% Ke Bo
%% Major results Figure4
% 
% We studied the relationship between task brain expression, MitoD group, Mito disease severity, and behavioral performance. Note due to lack of behavioral measurement, we'll not be possible to test behavioral performance in Cold pain and multisensory tasks.
% 
% 
% Load SVM decoded subject level pattern expression for all tasks. Note the 
% availbility of the data are different across tasks.
% 
% Multisensory: 90 total, 28 MitoD
% 
% Cold pain: 91 total, 29 MitoD
% 
% Working memory: 88 total, 28 MitoD
% 
% T = 0.54, p = 0.59, n = 90), cold pain (T = –0.11 , p = 0.59, n = 91), or 
% N-back (T = 1.38, p = 0.19, n = 88)
% 
% bf01_Multisensory=t2smpbf(0.54,28,62)
% 
% bf01_ColdPain=t2smpbf(-0.11,29,62)
% 
% bf01_WorkingMemory=t2smpbf(1.38,28,60)
% 
% power_calc
% 
% Data Loading section


clear
load('E:\Mito_DICOM\SecondLevelSave\Thresholded\Multisensory_100_All.mat')
Y1=mean(dist_Perm,1);
Score1=Allscores(:,1);
GroupID1=GroupID;
index1=find(GroupID1==0);
age1=age;
Acc_1=mean(Acc);

load('E:\Mito_DICOM\SecondLevelSave\Thresholded\Cold_100_All.mat')
Y2=mean(dist_Perm,1);
Score2=Allscores(:,1);
GroupID2=GroupID;
index2=find(GroupID2==0);
age2=age;
Acc_2=mean(Acc);

load('E:\Mito_DICOM\SecondLevelSave\Thresholded\Stress_100_All.mat')
Y3=mean(dist_Perm,1);
Score3=Allscores(:,1);
GroupID3=GroupID;
index3=find(GroupID3==0);
age3=age;

% load('E:\Mito_DICOM\SecondLevelSave\Thresholded\Nback_100_All.mat')
load('E:\Mito_DICOM\SecondLevelSave\Thresholded\Nback_100_500Filter.mat')

Y4=mean(dist_Perm,1);
Score4=Allscores(:,1);
GroupID4=GroupID;
index4=find(GroupID4==0);
age4=age;
Acc_4=mean(Acc);
%% 
% Load GDF 15 data

[a b c]=xlsread('LOCAL_SOURCE\Dartmouth College Dropbox\Ke Bo\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\GDF15_export_2024-09-24.csv');

%% 
% *Compare decoding accuracy across task*

% Sample decoding accuracies
acc = [mean(Acc_1)*100, mean(Acc_2)*100, mean(Acc_4)*100];
labels = {'Multisensory', 'Cold Pain', 'N-back'};

% Define colors: green, blue, pink
bar_colors = [0.2, 0.8, 0.2;    % green
              0.2, 0.4, 0.8;    % blue
              0.9, 0.4, 0.6];   % pink

% Create the bar plot
figure;
b = bar(acc, 'FaceColor', 'flat', 'BarWidth', 0.5);  % BarWidth < 1 makes it thinner
hold on;

% Assign custom colors to each bar
b.CData = bar_colors;

% Add dashed line for chance level
yline(50, '--k', 'Chance', ...
    'LabelHorizontalAlignment', 'left', ...
    'LabelVerticalAlignment', 'bottom', ...
    'FontSize', 11);

% Format axes
ylim([40 100]);
set(gca, 'XTick', 1:3, 'XTickLabel', labels, 'FontSize', 12);
ylabel('Decoding Accuracy (%)');

% Add numerical values on top of bars
for i = 1:length(acc)
    text(i, acc(i) + 3, sprintf('%.1f%%', acc(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
end

box off;
set(gca,'fontsize',14)
set(gca,'fontweight','bold')
%% 
% 
% 
% *Compare group difference between Patient and control*


LineWidth=2;
Fontsize=12;
dotcolor2 = [0, 114, 189] / 255;    % Control - Blue
dotcolor1 = [217, 83, 25] / 255;    % Patient - Orange/Red
colorcoding=[dotcolor2; dotcolor1];
pointsize=4;
bw=0.4;
figure
subplot(1,3,1)
violinplot({Y1(GroupID1==1),Y1(GroupID1==0)},'facecolor',colorcoding,'mc','k','bw',bw,'plotlegend',0,'pointsize',pointsize)
set(gca,'xtick',[])
names={'Control';'MitoD'}


set(gca,'xtick',[1:2],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
ylabel('Multisensory brain pattern expression')
box off
[h p ci stat]=ttest2(Y1(GroupID1==1),Y1(GroupID1==0))


subplot(1,3,2)
violinplot({Y2(GroupID2==1),Y2(GroupID2==0)},'facecolor',colorcoding,'mc','k','bw',bw,'plotlegend',0,'pointsize',pointsize)
set(gca,'xtick',[])
names={'Control';'MitoD'}

set(gca,'xtick',[1:2],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
ylabel('Cold pain brain pattern expression')
box off

[h p ci stat]=ttest2(Y2(GroupID2==1),Y2(GroupID2==0))

subplot(1,3,3)
violinplot({Y4(GroupID4==1),Y4(GroupID4==0)},'facecolor',colorcoding,'mc','k','bw',bw,'plotlegend',0,'pointsize',pointsize)
% violinplot({Brain_S(X1_S==1),Brain_S(X1_S==0)},'facecolor',colorcoding,'mc','k','bw',bw,'plotlegend',0,'pointsize',pointsize)


set(gca,'xtick',[])
names={'Control';'MitoD'}

set(gca,'xtick',[1:2],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
ylabel('N-back brain pattern expression')
box off

[h p ci stat]=ttest2(Y4(GroupID4==1),Y4(GroupID4==0))
% [h p ci stat]=ttest2(Brain_S(X1_S==1),Brain_S(X1_S==0))
%% 
% *Basic correlation: Brain - Disease severity*


%%%% Correlation across 3 tasks %%%%%
figure

dotcolor2 = [0, 114, 189] / 255;    % Control - Blue
dotcolor1 = [217, 83, 25] / 255;    % Patient - Orange/Red

dotsize=500;
subplot(1,3,1)
scatter(Score1(:),Y1(:),dotsize,dotcolor2,'.')
hold on
scatter(Score1(find(GroupID1==0)),Y1(find(GroupID1==0)),dotsize,dotcolor1,'.')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
set(h(2),'Visible', 'off')

[R1 P1]=corr(Score1(find(GroupID1==0)),Y1(find(GroupID1==0))','type','Spearman')
[R2 P2]=partialcorr(Score1(find(GroupID1==0)),Y1(find(GroupID1==0))',age1(find(GroupID1==0))','type','Spearman')

% ylabel('Mean activation in significant activated region')
ylabel('Task activation score (Multisensory)')
xlabel('NMDAS(Disease severity) score')
% ylabel('CNS score')

set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)



subplot(1,3,2)

scatter(Score2(:),Y2(:),dotsize,dotcolor2,'.')
hold on
scatter(Score2(find(GroupID2==0)),Y2(find(GroupID2==0)),dotsize,dotcolor1,'.')
[R P]=corr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))','type','Spearman')
[R P]=partialcorr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))',age2(find(GroupID2==0))','type','Spearman')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
set(h(2),'Visible', 'off')
% ylabel('Mean activation in significant activated region')
ylabel('Task activation score (Cold pain)')
xlabel('NMDAS(Disease severity) score')
% ylabel('CNS score')

set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)

subplot(1,3,3)
scatter(Score4(:),Y4(:),dotsize,dotcolor2,'.')
hold on

scatter(Score4(find(GroupID4==0)),Y4(find(GroupID4==0)),dotsize,dotcolor1,'.')
[R P]=corr(Score4(find(GroupID4==0)),Y4(find(GroupID4==0))')
[R P]=partialcorr(Score4(find(GroupID4==0)),Y4(find(GroupID4==0))',age4(find(GroupID4==0))')

% ylabel('Mean activation in significant activated region')
ylabel('Task activation score (N-back)')
xlabel('NMDAS(Disease severity) score')
% ylabel('CNS score')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
set(h(2),'Visible', 'off')
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Advanced data analysis 
% Brain, behavior, and disease severity.
% 
% For working memory task, only part of behavioral data are complete. The participants 
% with icomplete or abnormal behavioral data will be excluded (Behavioral accuracy 
% below 40%).
% 
% This leave the data be: 
% 
% Control 74 participants
% 
% MitoD 25 participants
% Data loading for working memory Behavior

%%% Behavioral evidences %%%%%%%%%

[a b c]=xlsread('LOCAL_SOURCE\Dropbox (Dartmouth College)\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\MiSBIE MRI Meta Data 6-11-24.xlsm');
cellArray1=metaGrand; 
cellArray2=b;
ids1 = cellArray1(:, 1);
ids2 = cellArray2(2:end, 1);
% Find the indices of the IDs in cellArray1 that are also in cellArray2
[~, idx] = ismember(ids1, ids2);
Behavioral_ACC=a(idx,66);
%%% Behavioral ACC larger than 40% are selected. Chance level is 50%.
Behavioral_ACC_Index_AccBase=find(Behavioral_ACC>0.4);
Behavioral_ACC_Index_Mito=intersect(Behavioral_ACC_Index_AccBase,find(GroupID4==0))

% Save participants' data only with accurate behavioral performance
% X1: MitoD ID (Mito=0; Control=1)
% 
% X2: NMDAS score (Disease severity)
% 
% Y: Behavioral performance
% 
% M: Task pattern response
% 
% 
% 
% Scatter plot and correlation analysis to first test the relationship across 
% these variables

 




[R P]=corr(Score4(GroupID4==0),Y4(GroupID4==0)','type','Spearman')
[R P]=partialcorr(Score4(GroupID4==0),Y4(GroupID4==0)',age4(GroupID4==0)','type','Spearman')

%%% Correlation between disease severity and brain %%%
dotcolor2 = [0, 114, 189] / 255;    % Control - Blue
dotcolor1 = [217, 83, 25] / 255;    % Patient - Orange/Red

X1=GroupID4(Behavioral_ACC_Index_AccBase)';
X2=Score4(Behavioral_ACC_Index_AccBase,1);

X1_All=GroupID4';
X2_All=Score4;


Behavior=Behavioral_ACC(Behavioral_ACC_Index_AccBase);
Brain=Y4(Behavioral_ACC_Index_AccBase)';
Age=age4(Behavioral_ACC_Index_AccBase)';

Brain_All=Y4';
Age_All=age4';


[R_DS_Brain P_DS_Brain]=corr(X2(X1==0),Brain(X1==0),'type','Spearman')
[R_DS_Brain_partial P_DS_Brain_partial]=partialcorr(X2(X1==0),Brain(X1==0),Age(X1==0),'type','Spearman')
figure
% scatter(X2(:),Brain,dotsize,dotcolor2,'.')
% hold on
scatter(X2(X1==0),Brain(X1==0),dotsize,dotcolor1,'.')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
% set(h(2),'Visible', 'off')
% ylabel('Mean activation in significant activated region')
ylabel('Brain task p Aattern expression (n-back)')
xlabel('NMDAS(Disease severity) score')
% ylabel('CNS score')

set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)

%%% Correlation between disease severity and behavior %%%
[R_DS_Behavior P_DS_Behavior]=corr(X2(X1==0),Behavior(X1==0),'type','Spearman');
[R_DS_Behavior_partial P_DS_Behavior_partial]=partialcorr(X2(X1==0),Behavior(X1==0),Age(X1==0),'type','Spearman')

figure
[R P]=corr(X2(X1==0),Behavior(X1==0),'type','spearman')
scatter(X2(:),Behavior*100,dotsize,dotcolor2,'.')
hold on
scatter(X2(X1==0),Behavior(X1==0)*100,dotsize,dotcolor1,'.')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
% set(h(2),'Visible', 'off')
xlabel('NMDAS(Disease severity) score')
ylabel('N-back behavioral performance (%)')
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)

% X1 = 0 -> Patients
% X1 = 1 -> Controls
  
% Advanced analysis:
% Group comparison with T test statstic
% And divide patient into two groups

% Step 1: Filter out patients (X1 == 0)
patients_severity = X2(X1 == 0);


stats = clusterdata_permtest(patients_severity, ...
    'k',        1:4, ...
    'nperm',  100, ...
    'verbose', true);


[Groupidx, centers] = kmeans_Matlab(patients_severity, 2)
severe_idx=(Groupidx==Groupidx(1));
nonsevere_idx=(Groupidx==Groupidx(4));

% Your data (1x25)
x = patients_severity;         % 1x25
res = kmeans_bootstrap_stability(x, 2000);




[h p ci stat]=ttest2(Behavior(X1 == 0),Behavior(X1 == 1))

Behavior_Patient=Behavior(X1 == 0); 
Brain_Patient=Brain(X1 == 0);
Age_Patient=Age(X1 == 0);
[h p ci stat]=ttest2(Brain(X1==1),Brain_Patient(nonsevere_idx))
[h p ci stat]=ttest2(Brain(X1==1),Brain_Patient(severe_idx))

[h p ci stat]=ttest2(Behavior(X1==1),Behavior_Patient(nonsevere_idx))
[h p ci stat]=ttest2(Behavior(X1==1),Behavior_Patient(severe_idx))

[h p ci stat]=ttest2(Behavior_Patient(severe_idx),Behavior_Patient(nonsevere_idx))
[h p ci stat]=ttest2(Brain(X1==1),Brain_Patient(nonsevere_idx))
[h p ci stat]=ttest2(Behavior(nonsevere_idx),Behavior_Patient(severe_idx))
[h p ci stat]=ttest2(Behavior_Patient(nonsevere_idx),Behavior(X1 == 1))
[h p ci stat]=ttest2(Brain_Patient(severe_idx),Brain_Patient(nonsevere_idx))
[h p ci stat]=ttest2(Age_Patient(severe_idx),Age_Patient(nonsevere_idx))


%%
patients_severity_All = X2_All(X1_All == 0);

[Groupidx, centers] = kmeans_Matlab(patients_severity_All, 2)
severe_idx_All=(Groupidx==Groupidx(1));
nonsevere_idx_All=(Groupidx==Groupidx(5));
Brain_Patient_All=Brain_All(X1_All == 0);

[h p ci stat]=ttest2(Brain_All(X1_All==1),Brain_Patient_All(nonsevere_idx_All))
[h p ci stat]=ttest2(Brain_All(X1_All==1),Brain_Patient(severe_idx_All))
[h p ci stat]=ttest2(Brain_Patient(severe_idx_All),Brain_Patient_All(nonsevere_idx_All))

figure
% violinplot({Behavior(X1==1),Behavior_Patient(nonsevere_idx),Behavior_Patient(severe_idx)},'facecolor',colorcoding,'mc','k','bw',0.05,'plotlegend',0,'pointsize',8)
violinplot({Brain_All(X1==1),Brain_Patient_All(nonsevere_idx_All),Brain_Patient_All(severe_idx_All)},'facecolor',colorcoding,'mc','k','bw',0.4,'plotlegend',0,'pointsize',8)


% Figures for comparing brain pattern response for Control, mild and severe patients.

colorcoding=[dotcolor2;[255, 153, 51] / 255;[204, 0, 0] / 255];
LineWidth=2;
Fontsize=12;
figure
% violinplot({Behavior(X1==1),Behavior_Patient(nonsevere_idx),Behavior_Patient(severe_idx)},'facecolor',colorcoding,'mc','k','bw',0.05,'plotlegend',0,'pointsize',8)
violinplot({Brain(X1==1),Brain_Patient(nonsevere_idx),Brain_Patient(severe_idx)},'facecolor',colorcoding,'mc','k','bw',0.4,'plotlegend',0,'pointsize',8)

set(gca,'xtick',[])
names={'Control';'Mild Patients';'Severe Patients'}

set(gca,'xtick',[1:3],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
box off
% ylabel('Working memory task peroformance (%)')
ylabel('Brain task pattern expression (n-back)')

% violinplot({Y_Patient(nonsevere_idx),Y_Patient(severe_idx)})
figure
violinplot({100*Behavior(X1==1),100*Behavior_Patient(nonsevere_idx),100*Behavior_Patient(severe_idx)},'facecolor',colorcoding,'mc','k','bw',5,'plotlegend',0,'pointsize',8)
set(gca,'xtick',[])
names={'Control';'Mild Patients';'Severe Patients'}

set(gca,'xtick',[1:3],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
box off
% ylabel('Working memory task peroformance (%)')
ylabel('N-back task performance (%)')




[h p ci stat]=ttest2(M(X1==1),Brain_Patient(nonsevere_idx))

colorcoding=[dotcolor2; dotcolor1];

figure
violinplot({100*Behavior(X1==1),100*Behavior(X1==0)},'facecolor',colorcoding,'mc','k','bw',6,'plotlegend',0,'pointsize',8)
set(gca,'xtick',[])
names={'Control';'MitoD'}

set(gca,'xtick',[1:2],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
ylabel('N-back task performance (%)')
box off

[h p ci stat]=ttest2(Behavior(X1==1),Behavior(X1==0))
% Step 2: Sort the patients by severity in descending order
[sorted_severity, sorted_idx] = sort(patients_severity, 'descend');



% Advanced analysis:
% Multi-regression model and mediation analysis across 
% X1: MitoD ID (Mito=0; Control=1)
% 
% X2: NMDAS score (Disease severity)
% 
% Behavior: Behavioral performance
% 
% Brain: Task pattern response
% 
% 
% 
% *Model 1: Predict behavioral performance using MitoGroup, Brain, and age*
% 
% All three factors are significant

predictors = [X1 Brain Age];
model = fitlm(predictors, Behavior)
%% 
% *Model 2: Predict behavioral performance using MitoGroup, and age.* 
% 
% Without brain as predictor, MitoD group no longer predict behavior, suggesting 
% potential brain compensation in MitoD group

predictors = [X1 Age];
model = fitlm(predictors, Behavior)

%% 
% *Mediation analysis:*
% 
% test if the brain is the mediator of disease severity to predict behavioral 
% performance
% 
% 
% 
% Test if disease severity is the mediator of using brain to predict MitoGene 
% (Partial mediation) (n=74,MitoN=24)

X1_R=X1;
X1_R(find(X1==1))=0;
X1_R(find(X1==0))=1;

[paths, stats1, stats2] = mediation( X1_R,Brain, X2, 'boottop', 'stats', 'plots','covs',Age);
%% 
% Full group mediation analysis (n=88, MitoN=28)

[paths, stats1, stats2] = mediation( GroupID4',Y4', Score4, 'boottop', 'stats', 'plots','covs',age');
%% 
% Test if disease severity is the mediator of using MitoGene to predict behavior 
% (Full mediation)

[paths, stats1, stats2] = mediation( X1,Behavior, X2, 'boottop', 'stats', 'plots','covs',Age);

[paths, stats1, stats2] = mediation( X2,Behavior, Brain, 'boottop', 'stats', 'plots','covs',Age);

%% Same analysis for multisensory 

X1=GroupID1';
X2=Score1;
Brain=Y1';
Age=age1';
patients_severity = X2(X1 == 0);

% stats = clusterdata_permtest(patients_severity, ...
%     'k',        1:4, ...
%     'nperm',  1000, ...
%     'verbose', true);

%% 
% *Mediation analysis*

[paths, stats1, stats2] = mediation( X1,Brain, X2, 'boottop', 'stats', 'plots','covs',Age);
[paths, stats1, stats2] = mediation( X1,Brain, X2, 'boottop', 'stats', 'plots');

%% 
% *Severity Group division*

Brain_Patient=Brain(X1 == 0);
Age_Patient=Age(X1 == 0);

[Groupidx, centers] = kmeans_Matlab(patients_severity, 2)

severe_idx=(Groupidx==Groupidx(1));
nonsevere_idx=(Groupidx==Groupidx(6));

%% 
% T test Compare between severe and mild patient

[h p ci stat]=ttest2(Brain_Patient(severe_idx),Brain_Patient(nonsevere_idx))

[h p ci stat]=ttest2(Brain(X1==1),Brain_Patient(nonsevere_idx))
[h p ci stat]=ttest2(Brain(X1==1),Brain_Patient(severe_idx))

figure

colorcoding=[dotcolor2;[255, 153, 51] / 255;[204, 0, 0] / 255];

LineWidth=2;
Fontsize=12;
figure
% violinplot({Y(X1==1),Y_Patient(nonsevere_idx),Y_Patient(severe_idx)},'facecolor',colorcoding,'mc','k','bw',0.05,'plotlegend',0,'pointsize',8)
violinplot({Brain(X1==1),Brain_Patient(nonsevere_idx),Brain_Patient(severe_idx)},'facecolor',colorcoding,'mc','k','bw',0.8,'plotlegend',0,'pointsize',8)
 
set(gca,'xtick',[])
names={'Control';'Mild Patients';'Severe Patients'}

set(gca,'xtick',[1:3],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
box off
% ylabel('Working memory task peroformance (%)')
% ylabel('Cold Pain brain score')
ylabel('Brain task pattern expression (Multisensory)')

% Correlation within patient population

[R_DS_Brain P_DS_Brain]=corr(X2(X1==0),Brain(X1==0),'type','Spearman')
[R_DS_Brain_partial P_DS_Brain_partial]=partialcorr(X2(X1==0),Brain(X1==0),Age(X1==0),'type','Spearman')
figure
[R P]=corr(X2(X1==0),Brain(X1==0),'type','spearman')
scatter(X2(:),Brain,dotsize,dotcolor2,'.')
hold on
scatter(X2(X1==0),Brain(X1==0),dotsize,dotcolor1,'.')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
set(h(2),'Visible', 'off')
xlabel('NMDAS(Disease severity) score')
ylabel('Brain task pattern expression (Multisensory)')
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)

%% Same analysis for Cold pain

X1=GroupID2';
X2=Score2;
Brain=Y2';
Age=age2';

%% 
% *Mediation analysis*

[paths, stats1, stats2] = mediation( X1,Brain, X2, 'boottop', 'stats', 'plots','covs',Age);
[paths, stats1, stats2] = mediation( X1,Brain, X2, 'boottop', 'stats', 'plots');

%% 
% *Severity Group division*

Brain_Patient=Brain(X1 == 0);
Age_Patient=Age(X1 == 0);

[Groupidx, centers] = kmeans_Matlab(patients_severity, 2)

severe_idx=(Groupidx==Groupidx(1));
nonsevere_idx=(Groupidx==Groupidx(6));

%% 
% T test Compare between severe and mild patient

[h p ci stat]=ttest2(Brain_Patient(severe_idx),Brain_Patient(nonsevere_idx))

[h p ci stat]=ttest2(Brain(X1==1),Brain_Patient(nonsevere_idx))
[h p ci stat]=ttest2(Brain(X1==1),Brain_Patient(severe_idx))

figure

colorcoding=[dotcolor2;[255, 153, 51] / 255;[204, 0, 0] / 255];

LineWidth=2;
Fontsize=12;
figure
% violinplot({Y(X1==1),Y_Patient(nonsevere_idx),Y_Patient(severe_idx)},'facecolor',colorcoding,'mc','k','bw',0.05,'plotlegend',0,'pointsize',8)
violinplot({Brain(X1==1),Brain_Patient(nonsevere_idx),Brain_Patient(severe_idx)},'facecolor',colorcoding,'mc','k','bw',0.4,'plotlegend',0,'pointsize',8)
 
set(gca,'xtick',[])
names={'Control';'Mild Patients';'Severe Patients'}

set(gca,'xtick',[1:3],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
box off
% ylabel('Working memory task peroformance (%)')
% ylabel('Cold Pain brain score')
ylabel('Brain task pattern expression (Cold Pain)')

% Correlation within patient population

[R_DS_Brain P_DS_Brain]=corr(X2(X1==0),Brain(X1==0),'type','Spearman')
[R_DS_Brain_partial P_DS_Brain_partial]=partialcorr(X2(X1==0),Brain(X1==0),Age(X1==0),'type','Spearman')
figure
[R P]=corr(X2(X1==0),Brain(X1==0),'type','spearman')
scatter(X2(:),Brain,dotsize,dotcolor2,'.')
hold on
scatter(X2(X1==0),Brain(X1==0),dotsize,dotcolor1,'.')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
set(h(2),'Visible', 'off')
xlabel('NMDAS(Disease severity) score')
ylabel('Brain task pattern expression (Cold Pain)')
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)
% Analysis for GDF 15


X1=GroupID4';
X2=Score4;
Brain=Y4';
Age=age4';

X1_S=GroupID4(Behavioral_ACC_Index_AccBase)';
X2_S=Score4(Behavioral_ACC_Index_AccBase);
Brain_S=Y4(Behavioral_ACC_Index_AccBase)';
Age_S=age4(Behavioral_ACC_Index_AccBase)';
%%%GDF_15%%%


[a b c]=xlsread('LOCAL_SOURCE\Dartmouth College Dropbox\Ke Bo\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\GDF15_export_2024-09-24.csv')
load('E:\Mito_DICOM\SecondLevelSave\Nback_TwoRun_Right_order_Filter_PercentageChange_SPM_robust.mat')

cellArray=metaGrand(:,1);
doubleArray = zeros(1, length(cellArray));

% Loop through each element in the cell array
for i = 1:length(cellArray)
    % Extract the numeric part from the string and convert it to a double
    doubleArray(i) = str2double(cellArray{i}(3:end));
end
GDF=a(doubleArray,5)
GDF_15=GDF(Behavioral_ACC_Index_AccBase) 
[paths, stats1, stats2] = mediation( X1,Brain, log10(GDF), 'boottop', 'stats', 'plots');
[paths, stats1, stats2] = mediation( X1_S,Brain_S, log10(GDF_15), 'boottop', 'stats', 'plots');
[paths, stats1, stats2] = mediation( log(GDF_15),Behavior, Brain_S, 'boottop', 'stats', 'plots','covs',Age_S);
[paths, stats1, stats2] = mediation( log(GDF_15),Behavior, Brain_S, 'boottop', 'stats', 'plots');
[paths, stats1, stats2] = mediation( log(GDF_15(X1_S==0)),Behavior(X1_S==0), Brain_S(X1_S==0), 'boottop', 'stats', 'plots');


[paths, stats1, stats2] = mediation( X1_S,Behavior, log10(GDF_15), 'boottop', 'stats', 'plots');
[R P]=corr(Brain_S,Behavior,'row','complete','type','spearman')



%% *Equivalence analysis using n=88 (MitoD)*
% figure 4B

[R_DS_Brain P_DS_Brain]=corr(Score4(GroupID4==0),Y4(GroupID4==0)','type','Spearman')
[R_DS_Brain_partial P_DS_Brain_partial]=partialcorr(Score4(GroupID4==0),Y4(GroupID4==0)',age(GroupID4==0)','type','Spearman')
figure
scatter(Score4,Y4',dotsize,dotcolor2,'.')
hold on
scatter(Score4(GroupID4==0),Y4(GroupID4==0)',dotsize,dotcolor1,'.')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
set(h(2),'Visible', 'off')
xlabel('NMDAS(Disease severity) score')
ylabel('Brain task pattern expression (Multisensory)')
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)
%% 
% Figure 4D

figure
% violinplot({Behavior(X1==1),Behavior_Patient(nonsevere_idx),Behavior_Patient(severe_idx)},'facecolor',colorcoding,'mc','k','bw',0.05,'plotlegend',0,'pointsize',8)
violinplot({Brain_All(X1==1),Brain_Patient_All(nonsevere_idx_All),Brain_Patient_All(severe_idx_All)},'facecolor',colorcoding,'mc','k','bw',0.4,'plotlegend',0,'pointsize',8)
set(gca,'xtick',[])
names={'Control';'Mild Patients';'Severe Patients'}

set(gca,'xtick',[1:3],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
box off
% ylabel('Working memory task peroformance (%)')
ylabel('Brain task pattern expression (n-back)')

%% 
% Figure 4C

[paths, stats1, stats2] = mediation( GroupID4',Y4', Score4, 'boottop', 'stats', 'plots','covs',age');
%% 
% Figure 3C

figure
dotcolor2 = [0, 114, 189] / 255;    % Control - Blue
dotcolor1 = [217, 83, 25] / 255;    % Patient - Orange/Red
colorcoding=[dotcolor2; dotcolor1];
violinplot({Y4(GroupID4==1),Y4(GroupID4==0)},'facecolor',colorcoding,'mc','k','bw',bw,'plotlegend',0,'pointsize',pointsize)
% violinplot({Brain_S(X1_S==1),Brain_S(X1_S==0)},'facecolor',colorcoding,'mc','k','bw',bw,'plotlegend',0,'pointsize',pointsize)


set(gca,'xtick',[])
names={'Control';'MitoD'}

set(gca,'xtick',[1:2],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
ylabel('N-back brain pattern expression')
box off

[h p ci stat]=ttest2(Y4(GroupID4==1),Y4(GroupID4==0))

%% 
% Figure 6B

figure

% Log transform GDF_15
logGDF = log10(GDF);
 
[R P]=corr(log10(GDF(X1==0)), Brain(X1==0),'rows','complete')

% Scatter plot for both groups
scatter(logGDF(GroupID4==0), Y4(GroupID4==0), dotsize, dotcolor1, '.') % Patients
hold on
scatter(logGDF(GroupID4==1), Y4(GroupID4==1), dotsize, dotcolor2, '.') % Controls

% Fit lines separately
% Patients
p_idx = GroupID4 == 0;
mdl_pat = fitlm(logGDF(p_idx), Y4(p_idx));
xfit = linspace(min(logGDF), max(logGDF), 100);
yfit_pat = predict(mdl_pat, xfit');
plot(xfit, yfit_pat, 'Color', dotcolor1, 'LineWidth', 2)

% Controls
c_idx = GroupID4 == 1;
mdl_ctrl = fitlm(logGDF(c_idx), Y4(c_idx));
yfit_ctrl = predict(mdl_ctrl, xfit');
plot(xfit, yfit_ctrl, 'Color', dotcolor2, 'LineWidth', 2)

xlabel('Plasma GDF-15 (log10)')
ylabel('Brain task pattern expression (N-back)')
set(gca, 'fontsize', 12, 'fontweight', 'bold', 'LineWidth', 1)
%% 
% 
% 
% Correlation analysis

[R P]=corr(log10(GDF(X1==0)),X2(X1==0),'row','complete','type','spearman')

[h p ci stat]=ttest2(log10(GDF(X1==0)),log10(GDF(X1==1)));

GDF_NaN_MitoD=isnan(GDF_15(X1_S==0));
GDF_NaN_Control=isnan(GDF_15(X1_S==1));
Brain_S(GDF_NaN)
mean(Brain_S)
mean(Brain_S(GDF_NaN))
[h p ci stat]=ttest2(Brain_S(X1_S==0),Brain_S(X1_S==1))
Brain_S_Control=Brain_S(X1_S==1);
Brain_S_MitoD=Brain_S(X1_S==0);
[h p ci stat]=ttest2(Brain_S_Control,Brain_S_MitoD)
[h p ci stat]=ttest2(Brain_S_Control(GDF_NaN_Control==0),Brain_S_MitoD(GDF_NaN_MitoD==0))

corr(log(GDF_15)) 

[R P]=corr(GDF(X1==0),Brain(X1==0),'row','complete','type','spearman')
[R P]=corr(log10(GDF_15(X1_S==0)),Brain_S(X1_S==0),'row','complete','type','spearman')
[R P]=corr(log(GDF_15(X1_S==1)),Brain_S(X1_S==1),'row','complete','type','spearman')
[R P]=corr(log10(GDF_15),Age_S,'row','complete')
[R P]=corr(Age_S,X1_S,'row','complete')

[R P]=corr(log(GDF),Age,'row','complete')
[R P]=corr(Age(X1==0),X2(X1==0),'row','complete','type','spearman')


[R P]=partialcorr(log10(GDF_15),Behavior,X1_S,'row','complete')
[R P]=partialcorr(log10(GDF_15),Brain_S,X1_S,'row','complete')
[R P]=corr(log10(GDF_15(X1_S==0)),Brain_S(X1_S==0),'row','complete')
[R P]=corr(log10(GDF_15(X1_S==1)),Brain_S(X1_S==1),'row','complete')


[R P]=corr(log(GDF_15(X1_S==0)),Behavior(X1_S==0),'row','complete')
[R P]=partialcorr(log(GDF_15(X1_S==0)),Behavior(X1_S==0),Age_S(X1_S==0),'row','complete')

[R P]=corr(log(GDF_15(X1_S==1)),Behavior(X1_S==1),'row','complete')
[R P]=partialcorr(log(GDF_15(X1_S==1)),Behavior(X1_S==1),Age_S(X1_S==1),'row','complete')



[R P]=corr(log(GDF_15(X1_S==1)),Brain_S(X1_S==1),'row','complete','type','spearman')
[R P]=corr(log(GDF_15(X1_S==0)),Brain_S(X1_S==0),'row','complete','type','spearman')


[R P]=corr(log(GDF_15),Behavior,'row','complete','type','spearman')
[R P]=corr(log(GDF_15),Brain_S,'row','complete','type','spearman')
[R P]=corr(Brain_S,Behavior,'row','complete')


[R P]=corr(Brain_S, Behavior,'row','complete','type','spearman')

[R P]=partialcorr(log(GDF),X2,Age,'row','complete')

[R P]=corr(log10(GDF(X1==0)),Brain(X1==0),'row','complete')
[R P]=corr(log10(GDF(X1==1)),Brain(X1==1),'row','complete')
[R P]=partialcorr(log10(GDF(X1==0)),Brain(X1==0),Age(X1==0),'row','complete')
[R P]=partialcorr(log10(GDF(X1==1)),Brain(X1==1),Age(X1==1),'row','complete')
[R P]=partialcorr(log10(GDF),Brain,[Age X1],'row','complete')


[R P]=corr(log10(GDF_15(X1_S==0)),Brain_S(X1_S==0),'row','complete')
[R P]=corr(log10(GDF_15(X1_S==1)),Brain_S(X1_S==1),'row','complete')
[R P]=partialcorr(log10(GDF_15(X1_S==0)),Brain_S(X1_S==0),Age_S(X1_S==0),'row','complete')
[R P]=partialcorr(log10(GDF_15(X1_S==1)),Brain_S(X1_S==1),Age_S(X1_S==1),'row','complete')



[R P]=partialcorr(log10(GDF_15),Behavior,[X1_S],'row','complete')

[R P]=partialcorr(log(GDF),Brain,X1,'row','complete')
[R P]=partialcorr(log(GDF_15),Behavior,X1_S,'row','complete')

[R P]=corr(Age,Brain,'row','complete')
[R P]=corr(Age(X1==1),log10(GDF(X1==1)),'row','complete')

[R P]=corr(Age_S,Behavior,'row','complete')
[R P]=corr(Age(X1==0),Brain(X1==0),'row','complete')
[R P]=corr(Age_S(X1_S==0),Behavior(X1_S==0),'row','complete')
[R P]=corr(Age(X1==1),Brain(X1==1),'row','complete')
[R P]=corr(Age_S(X1_S==1),Behavior(X1_S==1),'row','complete')

GDF_patient=GDF_15(X1_S==0)


[h p ci stat]=ttest2(GDF_15(X1_S==1),GDF_patient(severe_idx))
% [h p ci stat]=ttest2(GDF_patient(nonsevere_idx),GDF_patient(severe_idx))

Group=X1;
NMDAS=X2;
GDF_Log=log10(GDF);
T = table(Age, Group, Brain, NMDAS,GDF_Log);

% Interaction regression on brain expression
GDF_Log_M=GDF_Log-nanmean(GDF_Log);
Brain_M=Brain-nanmean(Brain);
Brain_S_M=Brain_S-nanmean(Brain_S);
GDF_15_Log=log10(GDF_15)
GDF15_Log_M=GDF_15_Log-nanmean(GDF_15_Log);
Group_M=Group;
Group_M(Group_M==0)=-1;
Age_M=Age-nanmean(Age);
Group_S=X1_S;
GDF15_Log_M=GDF_15_Log-nanmean(GDF_15_Log);
Behavior_M=Behavior-nanmean(Behavior);
Group_S_M=Group_S;
Group_S_M(Group_S_M==0)=-1;
Age_S_M=Age_S-nanmean(Age_S);
NisnanGDF=find(isnan(GDF)==0);


T = table(GDF15_Log_M, Age_S_M, Group_S_M, Brain_S_M);

% 3) Model 1: main effects only
mdl1 = fitlm(T, 'Brain_S_M ~ GDF15_Log_M + Group_S_M + Age_S_M');
disp('--- Model 1: Main effects only ---');
disp(mdl1);
fprintf('R^2 = %.3f\n\n', mdl1.Rsquared.Ordinary);

% 4) Model 2: add two interaction terms
%    • GDF15:age tests whether the GDF-15 slope differs across ages
%    • GDF15:severity tests whether the GDF-15 slope differs by disease severity
% Brain=mitogroup + GDF15 + age + GDF15*mitogroup + age * mitogroup
mdl2 = fitlm(T, 'Brain_S_M ~  GDF15_Log_M*Age_S_M*Group_S_M');
disp('--- Model 2: + interactions ---');
disp(mdl2);
fprintf('R^2 = %.3f\n\n', mdl2.Rsquared.Ordinary);


%%%%%%%%Behavior %%%
Group_S=X1_S;
GDF_15_Log=log10(GDF_15)
GDF15_Log_M=GDF_15_Log-nanmean(GDF_15_Log);
Behavior_M=Behavior-nanmean(Behavior);
Group_S_M=Group_S;
Group_S_M(Group_S_M==0)=-1;
Age_S_M=Age_S-nanmean(Age_S);
NisnanGDF=find(isnan(GDF)==0);

% 5) Compare models to see if interactions significantly improve fit
% cmp = anova(mdl1, mdl2, 'summary');
% disp('--- Model comparison (ANOVA) ---');
% disp(cmp);

T = table(GDF15_Log_M, Age_S_M, Group_S_M, Behavior_M);


% 3) Model 1: main effects only
mdl1 = fitlm(T, 'Behavior_M ~ GDF15_Log_M + Group_S_M + Age_S_M');
disp('--- Model 1: Main effects only ---');
disp(mdl1);
fprintf('R^2 = %.3f\n\n', mdl1.Rsquared.Ordinary);

% 4) Model 2: add two interaction terms
%    • GDF15:age tests whether the GDF-15 slope differs across ages
%    • GDF15:severity tests whether the GDF-15 slope differs by disease severity
% Brain=mitogroup + GDF15 + age + GDF15*mitogroup + age * mitogroup
T = table(GDF15_Log_M, Age_S_M, Group_S_M, Behavior_M);
mdl2 = fitlm(T, 'Behavior_M ~ GDF15_Log_M*Age_S_M*Group_S_M');
disp('--- Model 2: + interactions ---');
disp(mdl2);
fprintf('R^2 = %.3f\n\n', mdl2.Rsquared.Ordinary);

% 5) Compare models to see if interactions significantly improve fit
% cmp = anova(mdl1, mdl2, 'summary');
% disp('--- Model comparison (ANOVA) ---');
% disp(cmp);

fitlm([GDF15_Log_M(X1_S==0) X2_S(X1_S==0)],Brain_S_M(X1_S==0),'interactions')
fitlm([GDF_Log(X1==1) Age(X1==1)],Brain(X1==1))
fitlm([X2_S(X1_S==0) Age_S(X1_S==0)],Behavior(X1_S==0))
%% 
% %%%%%%%%%%%test

T = table(GDF_Log, Age_S, Group_S, Brain);



% 3) Model 1: main effects only
mdl1 = fitlm(T, 'Brain ~ GDF_Log + Group + Age');
disp('--- Model 1: Main effects only ---');
disp(mdl1);
fprintf('R^2 = %.3f\n\n', mdl1.Rsquared.Ordinary);

% 4) Model 2: add two interaction terms
%    • GDF15:age tests whether the GDF-15 slope differs across ages
%    • GDF15:severity tests whether the GDF-15 slope differs by disease severity
% Brain=mitogroup + GDF15 + age + GDF15*mitogroup + age * mitogroup
mdl1 = fitlm(T, 'Brain ~  GDF_Log*Group + Age*Group +Age*GDF_Log + Age + Group + GDF_Log');
disp('--- Model 2: + interactions ---');
disp(mdl2);
fprintf('R^2 = %.3f\n\n', mdl2.Rsquared.Ordinary);



%% 
% 
% 
% Figure for GDF 15 and brain

figure
dotcolor2 = [0, 114, 189] / 255;    % Control - Blue
dotcolor1 = [217, 83, 25] / 255;    % Patient - Orange/Red
colorcoding=[dotcolor2; dotcolor1];
violinplot({Age_S(X1_S==1),Age_S(X1_S==0)},'facecolor',colorcoding,'mc','k','bw',5,'plotlegend',0,'pointsize',4)
% legend('Control','MitoD')
set(gca,'xtick',[])
names={'Control';'MitoD'}

set(gca,'xtick',[1:2],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
box off
ylabel('Age (year)')

[h p ci stat]=ttest2(Age_S(X1_S==1),Age_S(X1_S==0))

figure

% Log transform GDF_15
logGDF = log10(GDF);

% [R P]=corr(log10(GDF(X1==1)), Brain(X1==1),'complete')

% Scatter plot for both groups
scatter(logGDF(X1==1), Brain(X1==1), dotsize, dotcolor2, '.') % Patients
hold on
scatter(logGDF(X1==0), Brain(X1==0), dotsize, dotcolor1, '.') % Controls

% Fit lines separately
% Patients
p_idx = X1 == 0;
mdl_pat = fitlm(logGDF(p_idx), Brain(p_idx));
xfit = linspace(min(logGDF), max(logGDF), 100);
yfit_pat = predict(mdl_pat, xfit');
plot(xfit, yfit_pat, 'Color', dotcolor2, 'LineWidth', 2)

% Controls
c_idx = X1 == 1;
mdl_ctrl = fitlm(logGDF(c_idx), Brain(c_idx));
yfit_ctrl = predict(mdl_ctrl, xfit');
plot(xfit, yfit_ctrl, 'Color', dotcolor1, 'LineWidth', 2)

xlabel('Plasma GDF-15 (log10)')
ylabel('Brain task pattern expression (N-back)')
set(gca, 'fontsize', 12, 'fontweight', 'bold', 'LineWidth', 1)

%%%% Subject exclusion %%%
figure

% Log transform GDF_15
logGDF = log10(GDF_15);

% [R P]=corr(log10(GDF(X1==1)), Brain(X1==1),'complete')

% Scatter plot for both groups
scatter(logGDF(X1_S==1), Brain_S(X1_S==1), dotsize, dotcolor2, '.') % Patients
hold on
scatter(logGDF(X1_S==0), Brain_S(X1_S==0), dotsize, dotcolor1, '.') % Controls

% Fit lines separately
% Patients
p_idx = X1_S == 0;
mdl_pat = fitlm(logGDF(p_idx), Brain_S(p_idx));
xfit = linspace(min(logGDF), max(logGDF), 100);
yfit_pat = predict(mdl_pat, xfit');
plot(xfit, yfit_pat, 'Color', dotcolor1, 'LineWidth', 2)

% Controls
c_idx = X1_S == 1;
mdl_ctrl = fitlm(logGDF(c_idx), Brain_S(c_idx));
yfit_ctrl = predict(mdl_ctrl, xfit');
plot(xfit, yfit_ctrl, 'Color', dotcolor2, 'LineWidth', 2)

xlabel('Plasma GDF-15 (log10)')
ylabel('Brain task pattern expression (N-back)')
set(gca, 'fontsize', 12, 'fontweight', 'bold', 'LineWidth', 1)



figure
scatter(logGDF(:),Behavior(:),dotsize,dotcolor2,'.')
hold on
scatter(logGDF(find(X1_S==0)),Behavior(find(X1_S==0)),dotsize,dotcolor1,'.')
% [R P]=corr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))','type','Spearman')
% [R P]=partialcorr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))',age2(find(GroupID2==0))','type','Spearman')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
set(h(2),'Visible', 'off')
% ylabel('Mean activation in significant activated region')
ylabel('Task activation score (N-back)')
xlabel('Plasma GDF-15 (log10)')
% ylabel('CNS score')

set(gca,'fontsize',12,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)

%% 
% plot showing GDF 15 is related to MitoD and diseaseverity
% 
% 


LineWidth=2;
Fontsize=12;
dotcolor2 = [0, 114, 189] / 255;    % Control - Blue
dotcolor1 = [217, 83, 25] / 255;    % Patient - Orange/Red
colorcoding=[dotcolor2; dotcolor1];
pointsize=4;
bw=0.1;
figure
violinplot({log10(GDF(GroupID4==1)),log10(GDF(GroupID4==0))},'facecolor',colorcoding,'mc','k','bw',bw,'plotlegend',0,'pointsize',pointsize)
set(gca,'xtick',[])
names={'Control';'MitoD'}

set(gca,'xtick',[1:2],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
ylabel('Plasma GDF-15 (log10)')
box off
[h p ci stat]=ttest2(log10(GDF(GroupID4==1)),log10(GDF(GroupID4==0)))

figure
[R P]=corr(log10(GDF(GroupID4==0)),Score4(GroupID4==0),'rows','complete','type','Spearman')
% scatter(log10(GDF),Score4,dotsize,dotcolor2,'.')
hold on
scatter(log10(GDF(GroupID4==0)),Score4(GroupID4==0),dotsize,dotcolor1,'.')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
% set(h(2),'Visible', 'off')
xlabel('Plasma GDF-15 (log10)')
ylabel('NMDAS (Disease severity) Score')
set(gca,'fontsize',12,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)
%% 
% Figures for GDF 15 and behavior


figure

% Log transform GDF_15
logGDF = log10(GDF_15(:));

% Scatter plot for both groups
scatter(logGDF(X1_S==0), 100*Behavior(X1_S==0), dotsize, dotcolor1, '.') % Patients
hold on
scatter(logGDF(X1_S==1), 100*Behavior(X1_S==1), dotsize, dotcolor2, '.') % Controls

% Fit lines separately
% Patients
p_idx = X1_S == 0;
mdl_pat = fitlm(logGDF(p_idx), 100*Behavior(p_idx));
xfit = linspace(min(logGDF), max(logGDF), 100);
yfit_pat = predict(mdl_pat, xfit');
plot(xfit, yfit_pat, 'Color', dotcolor1, 'LineWidth', 2)

% Controls
c_idx = X1_S == 1;
mdl_ctrl = fitlm(logGDF(c_idx), 100*Behavior(c_idx));
yfit_ctrl = predict(mdl_ctrl, xfit');
plot(xfit, yfit_ctrl, 'Color', dotcolor2, 'LineWidth', 2)

xlabel('Plasma GDF-15 (log10)')
ylabel('N-back behavior performance (%)')
set(gca, 'fontsize', 12, 'fontweight', 'bold', 'LineWidth', 1)
%% 
% 

scatter(GDF_15(:),M(:),dotsize,dotcolor2,'.')
hold on
scatter(GDF_15(find(X1==0)),M(find(X1==0)),dotsize,dotcolor1,'.')
% [R P]=corr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))','type','Spearman')
% [R P]=partialcorr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))',age2(find(GroupID2==0))','type','Spearman')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
set(h(2),'Visible', 'off')
% ylabel('Mean activation in significant activated region')
ylabel('Task activation score (N-back)')
xlabel('GDF15 (saliva PostMR)')
% ylabel('CNS score')

set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)


%%%%%

[a b c]=xlsread('LOCAL_SOURCE\Dartmouth College Dropbox\Ke Bo\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\GDF15_export_2024-09-24.csv')
% load('E:\Mito_DICOM\SecondLevelSave\Multisensory_OneModel_Filter_Robust_SPM.mat')
% load('E:\Mito_DICOM\SecondLevelSave\Nback_TwoRun_Right_order_Filter_PercentageChange_SPM_Robust.mat')
load('E:\Mito_DICOM\SecondLevelSave\Thresholded\Nback_100_500Filter.mat')


X1=GroupID4(Behavioral_ACC_Index_AccBase)';
X2=Score4(Behavioral_ACC_Index_AccBase,1);

Y=Behavioral_ACC(Behavioral_ACC_Index_AccBase);
Brain=Y4(Behavioral_ACC_Index_AccBase)';

cellArray=metaGrand(:,1);
doubleArray = zeros(1, length(cellArray));

% Loop through each element in the cell array
for i = 1:length(cellArray)
    % Extract the numeric part from the string and convert it to a double
    doubleArray(i) = str2double(cellArray{i}(3:end));
end
GDF=log10(a(doubleArray,5))

GDF_15=GDF(Behavioral_ACC_Index_AccBase);
[R P]=corr(GDF_15(X1==0),M(X1==0),'row','complete','type','spearman')
[R P]=corr(GDF_15(X1==0),Y(X1==0),'row','complete','type','spearman')
[R P]=corr(GDF_15(X1==0),X2(X1==0),'row','complete','type','spearman')
[R P]=corr(GDF_15,X2,'row','complete')



figure
scatter(GDF_15(:),M(:),dotsize,dotcolor2,'.')
hold on
scatter(GDF_15(find(X1==0)),M(find(X1==0)),dotsize,dotcolor1,'.')
% [R P]=corr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))','type','Spearman')
% [R P]=partialcorr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))',age2(find(GroupID2==0))','type','Spearman')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
set(h(2),'Visible', 'off')
% ylabel('Mean activation in significant activated region')
ylabel('Task activation score (N-back)')
xlabel('Log(GDF15 density) (Plasma)')
% ylabel('CNS score')

set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)

figure
scatter(GDF_15(:),Y(:)*100,dotsize,dotcolor2,'.')
hold on
scatter(GDF_15(find(X1==0)),Y(find(X1==0))*100,dotsize,dotcolor1,'.')
% [R P]=corr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))','type','Spearman')
% [R P]=partialcorr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))',age2(find(GroupID2==0))','type','Spearman')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
set(h(2),'Visible', 'off')
% ylabel('Mean activation in significant activated region')
ylabel('N-back behavioral performance (%)')
xlabel('Log(GDF15 density) (Plasma)')
% ylabel('CNS score')

set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)


figure
scatter(GDF_15(:),Y(:)*100,dotsize,dotcolor2,'.')
hold on
% [R P]=corr(GDF_15(find(X1==0)),Y(find(X1==0).*100),'row','complete','type','spearman')
% [R P]=corr(GDF_15,Y,'row','complete')

scatter(GDF_15(find(X1==0)),Y(find(X1==0))*100,dotsize,dotcolor1,'.')
% [R P]=corr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))','type','Spearman')
% [R P]=partialcorr(Score2(find(GroupID2==0)),Y2(find(GroupID2==0))',age2(find(GroupID2==0))','type','Spearman')
h = lsline
set(h(1),'color','#FF594C','linewidth',2)
set(h(2),'Visible', 'off')
% ylabel('Mean activation in significant activated region')
ylabel('Behavioral performance (%)')
% xlabel('GDF15 (saliva PostMR)')
% ylabel('GDF15 (saliva)')
xlabel('GDF15 (Plasma)')
% ylabel('CNS score')

set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)

%%
[a b c]=xlsread('LOCAL_SOURCE\Dartmouth College Dropbox\Ke Bo\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\Misbie_Lactate.xlsx')
load('E:\Mito_DICOM\SecondLevelSave\Nback_TwoRun_Right_order_Filter_PercentageChange_SPM_robust.mat')

cellArray=metaGrand(:,1);
doubleArray = zeros(1, length(cellArray));

% Loop through each element in the cell array
for i = 1:length(cellArray)
    % Extract the numeric part from the string and convert it to a double
    doubleArray(i) = str2double(cellArray{i}(3:end));
end
Lactate=a(doubleArray,5)