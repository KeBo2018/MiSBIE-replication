% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
clear 
% %%%%%%  Load 1st level contrast map and meta data %%%%%%%%%%
 load('E:\Mito_DICOM\SecondLevelSave\Nback_TwoRun_Right_order_500Filter_PercentageChange_SPM_Robust.mat')
%  load('E:\Mito_DICOM\SecondLevelSave\Cold_Onemodel_Robust_Transition_SPM.mat')
%  load('E:\Mito_DICOM\SecondLevelSave\Multisensory_OneModel_Filter_Robust_SPM.mat')
%  load('E:\Mito_DICOM\SecondLevelSave\Stress_onemodel_Robust_SPM.mat')

% load('E:\Mito_DICOM\SecondLevelSave\Nback_TwoRun_Right_order_Robust.mat')
% load('F:\SecondLevelSave\Multisensory_Onemodel_GroupID.mat')
% Contrast=fmri_data(image_math(Cold,Recover,'minus'))
% Contrast=fmri_data(Cold);
Contrast=fmri_data(Contrast);

Condition1=fmri_data(Twoback);Condition2=fmri_data(Zeroback);
% Condition1=fmri_data(Cold);Condition2=fmri_data(Control);
% Condition1=fmri_data(Task);Condition2=fmri_data(Rest);
% Condition1=fmri_data(Stress);Condition2=fmri_data(Control);


Contrast=fmri_data(Contrast);
T=ttest(Contrast)
Contrast.removed_images=[];
Contrast.removed_images=zeros(size(Contrast.dat,2),1);
Condition1.removed_images=[];
Condition2.removed_images=[];

T_Threshold=threshold(T,0.05, 'fdr');
orthviews(T_Threshold);
figure
montage(T_Threshold,'compact3')

%%%% Task selection
for perm=1:100
%     [paired_d, stats_SVM, optout, cat_obj]=canlab_run_paired_SVM(fmri_data(Cold),fmri_data(Control))
    [paired_d, stats_SVM, optout, cat_obj]=canlab_run_paired_SVM(fmri_data(Twoback),fmri_data(Zeroback))
%     [paired_d, stats_SVM, optout, cat_obj]=canlab_run_paired_SVM(fmri_data(Task),fmri_data(Rest))
%     [paired_d, stats_SVM, optout, cat_obj]=canlab_run_paired_SVM(fmri_data(Stress),fmri_data(Control))
    Acc(perm)=stats_SVM.paired_accuracy;
    D_Value(perm)=stats_SVM.paired_d;
    
    %
dist_Perm(perm,:)=stats_SVM.paired_hyperplane_dist_scores;
double_dist_Perm(perm,:)=stats_SVM.dist_from_hyperplane_xval;
paired_d_Perm(perm)=paired_d;
end

[paired_d, stats_SVM, optout, cat_obj]=canlab_run_paired_SVM(fmri_data(Twoback),fmri_data(Zeroback),'dobootstrap',1,'boot_n',5000)

figure
% montage(stats_SVM.weight_obj,'full')
Temp=stats_SVM.weight_obj;
Temp_T=threshold(Temp,0.05,'fdr')
figure
montage(Temp_T,'full')

Temp_T=threshold(Temp,0.005,'uncorrected')
figure
montage(Temp_T,'full')

data=mean(double_dist_Perm,1);
% Split data into two conditions (Condition 1: first 85, Condition 2: last 85)
condition1 = data(1:length(data)/2);
condition2 = data(length(data)/2+1:end);

% Create x-axis (subject ID) for each subject
subjectID = 1:length(data)/2;

% Plot
figure;
hold on;

taskColor = [255/255, 0, 0];  % Red for Condition 1 (Task)
controlColor = [0, 166/255, 68/255];  % Green for Condition 2 (Control)

% Plot condition 1 values (first condition for each subject)
plot(subjectID, condition1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Plot condition 2 values (second condition for each subject)
plot(subjectID, condition2, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', controlColor);

% Connect the points with arrows (or lines) from condition 1 to condition 2
for i = 1:length(data)/2
    % Draw lines connecting the two points for each subject
    plot([subjectID(i), subjectID(i)], [condition1(i), condition2(i)], 'k-');
end
xlabel('SubjectID')
ylabel('Distance to hyperpline')
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)

figure
violinplot(D_Value')
stats_SVM.paired_hyperplane_dist_scores=mean(dist_Perm,1)';
paired_d_Perm=mean(paired_d_Perm);
create_figure('ROC')


ROC = roc_plot(stats_SVM.yfit, logical(stats_SVM.Y > 0), 'color', 'r');
create_figure('subjects');
plot(stats_SVM.paired_hyperplane_dist_scores, 'o');
plot_horizontal_line(0);
xlabel('Participant'); ylabel('Classifier score');


% index=str2num(Cold.image_names(:,7:9))


%  [a b c]=xlsread('LOCAL_SOURCE\Dropbox (Dartmouth College)\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\MiSBIE MRI Meta Data 6-11-24.xlsm')
% c(1,:)=[];
% metaGrand=c(index,:);

[combined_severity Allscores GroupID age]=combine_SeverityScore(metaGrand);

sex=zeros(1,size(metaGrand,1));
Mito_Sub=zeros(1,size(metaGrand,1));
for i=1:size(metaGrand,1)
    sex_char=cell2mat(metaGrand(i,4));
    if strcmp(sex_char,'Male')==1
        sex(i)=1;
    else
        sex(i)=-1;
    end

    Mito_Char=cell2mat(metaGrand(i,2));
    if strcmp(Mito_Char,'Single Deletion')==1
        Mito_Sub(i)=1;
    else
        Mito_Sub(i)=-1;
    end
end
%         sex=cell2mat(metaGrand(:,4));
% 
%  %%%% ID select %%%%%
% subjectIDs1 = metaGrand(:, 1);
% subjectIDs2 = c(:, 1);
% 
% % Initialize an empty cell array to store the filtered rows
% Updated_metaGrand = {};
% 
% % Loop through the subject IDs in cellArray1
% for i = 1:length(subjectIDs1)
%     % Find the rows in cellArray2 that match the current subject ID in cellArray1
%     matchingRows = strcmp(subjectIDs2, subjectIDs1{i});
%     
%     % Append the matching rows to the filteredCellArray2
%     Updated_metaGrand = [Updated_metaGrand; c(matchingRows, :)];
% end
 
 %  load('E:\Mito_DICOM\combined_severity.mat')
% load('F:\SecondLevelSave\Cold_OneModel.mat')
 %%%% T Test for 2nd level analysis %%%%%
%  Contrast=image_math(Cold,Control,'minus')
Contrast=fmri_data(Contrast);
T=ttest(Contrast)
% Contrast.removed_images=[];
Contrast.removed_images=zeros(size(Contrast.dat,2),1);
Condition1.removed_images=[];
Condition2.removed_images=[];

T_Threshold=threshold(T,0.05, 'fdr');
orthviews(T_Threshold);
figure
montage(T_Threshold,'compact3')

%%%%%  2nd level robust regression

% Robust_Contrast=robfit_parcelwise(Contrast,'csf_wm_covs',1)
% 
% table_of_atlas_regions_covered(region(T_Threshold))

 
%%%%%%%%%%%%% Regress out white matter and csf on second level %%%%

[group_metrics, individual_metrics, global_gm_wm_csf_values] = qc_metrics_second_level(Contrast)
Contrast.X=[ global_gm_wm_csf_values(:,3:3)];
Covariate=[ age'];

GroupCon=regress(Contrast,'residual','nointercept')

GroupCon_T=ttest(GroupCon.resid);
orthviews(GroupCon_T)
figure
GroupCon_T_Thre=threshold(GroupCon_T,0.05, 'fdr');

% figure
% montage(GroupCon_T_Thre,'full')


% write(GroupCon.resid, 'fname', 'GroupCon_Residual.nii')
%%%%% Compare before and after WM/CSF regression %%%%%

% figure
% 
% corr(GroupCon_T.dat,T.dat)
% figure
% scatter(T.dat,GroupCon_T.dat)
% h=lsline
% xlabel('Before 2nd level regression')
% ylabel('After 2nd level regression ')


%%%%%%%%%%% Load group level data and group ID %%%%%%%%%%%%%%%

% NMDAS_Score=cell2mat(Updated_metaGrand(:,6));
% CNS_Score=cell2mat(Updated_metaGrand(:,7));
% NMDAS_Score2=cell2mat(metaGrand(:,6));
% CNS_Score2=cell2mat(metaGrand(:,7));
%%%%%%%% Load activated region %%%%%

ActivationIndex_Pos=find(GroupCon_T_Thre.dat>0 & GroupCon_T_Thre.sig>0);
ActivationIndex_Neg=find(GroupCon_T_Thre.dat<0 & GroupCon_T_Thre.sig>0);

PosActivation=nanmean(GroupCon.resid.dat(ActivationIndex_Pos,:),1)';
NegActivation=nanmean(GroupCon.resid.dat(ActivationIndex_Neg,:),1)';

for i=1:7
    [R1(i) P1(i)]=corr(Allscores(:,i),PosActivation,'rows','complete','type','Spearman')
    [R2(i) P2(i)]=corr(Allscores(:,i),NegActivation,'rows','complete','type','Spearman')
end
 [R1(8) P1(8)]=corr(combined_severity,PosActivation,'rows','complete','type','Spearman')
 [R2(8) P2(8)]=corr(combined_severity,NegActivation,'rows','complete','type','Spearman')
[R1(9) P1(9)]=corr(combined_severity(find(GroupID==0)),PosActivation(find(GroupID==0)),'rows','complete','type','Spearman')
[R2(9) P2(9)]=corr(combined_severity(find(GroupID==0)),NegActivation(find(GroupID==0)),'rows','complete','type','Spearman')

[coeff, score, latent, tsquared, explained, mu] =pca(Allscores(:,1:2))

[Lambda, Psi, T, stats, F]=factoran(Allscores(:,1:7),2)
[R1(10) P1(10)]=corr(Allscores(find(GroupID==0),1),PosActivation(find(GroupID==0)),'rows','complete')
[R2(10) P2(10)]=corr(Allscores(find(GroupID==0),1),NegActivation(find(GroupID==0)),'rows','complete')

[R P]=corr(score(:,1),stats_SVM.paired_hyperplane_dist_scores,'rows','complete')
for i=1:7
    [R(i) P(i)]=corr(Allscores(find(GroupID==0),i),stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)),'rows','complete','type','Spearman')
end

figure
scatter(Allscores(find(GroupID==0),1),stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)))


[R(i) P(i)]=corr(Allscores(find(GroupID==0),i),stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)),'rows','complete','type','Spearman')

[h p ci stat]=ttest2(NegActivation(find(GroupID==0)),NegActivation(find(GroupID==1)))

% 
% 
% [R(1) P(1)]=corr(CNS_Score,PosActivation,'rows','complete','type','Spearman')
% [R(2) P(2)]=corr(CNS_Score,NegActivation,'rows','complete','type','Spearman')
% 
% [R(3) P(3)]=corr(NMDAS_Score,PosActivation,'rows','complete','type','Spearman')
% [R(4) P(4)]=corr(NMDAS_Score,NegActivation,'rows','complete','type','Spearman')
% 
% [R(5) P(5)]=corr(combined_severity,PosActivation,'rows','complete','type','Spearman')
% [R(6) P(6)]=corr(combined_severity,NegActivation,'rows','complete','type','Spearman')
% 
% % [R(5) P(5)]=corr(combined_severity,PosActivation,'rows','complete')
% % [R(6) P(6)]=corr(combined_severity,NegActivation,'rows','complete')
% 
% [R(7) P(7)]=corr(combined_severity(find(GroupID==0)),PosActivation(find(GroupID==0)),'rows','complete','type','Spearman')
% [R(8) P(8)]=corr(combined_severity(find(GroupID==0)),NegActivation(find(GroupID==0)),'rows','complete','type','Spearman')

X1=Allscores(:,1);
X2=Allscores(:,2);
Y=stats_SVM.paired_hyperplane_dist_scores;

[b, ~, residuals] = regress(X1, Covariate);
X1=residuals;
[b, ~, residuals] = regress(X2, Covariate);
X2=residuals;
[b, ~, residuals] = regress(Y, Covariate);
Y=residuals;

X=[X1 GroupID']
mdl1 = fitlm([ GroupID' Y ],[X1]) 

figure
dotcolor1=[3 141 204]/256;
dotcolor2=[0.8 0.3450 0.4330];
dotsize=500;
subplot(2,2,1)
scatter(Allscores(:,1),stats_SVM.paired_hyperplane_dist_scores,dotsize,dotcolor1,'.')
% ylabel('Mean activation in significant activated region')
ylabel('Task activation score')
xlabel('NMDAS(Disease severity) score')
% ylabel('CNS score')
h=lsline
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)
hold on
scatter(Allscores(find(GroupID==0),1),stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)),dotsize,dotcolor2,'.')

subplot(2,2,2)
scatter(Allscores(find(GroupID==0),1),stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)),dotsize,dotcolor2,'.')

% ylabel('Mean activation in significant activated region')
ylabel('Task activation score')
xlabel('NMDAS(Disease severity) score')
% ylabel('CNS score')
h=lsline
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)

subplot(2,2,3)
scatter(Allscores(:,2),stats_SVM.paired_hyperplane_dist_scores,dotsize,dotcolor1,'.')
% ylabel('Mean activation in significant activated region')
ylabel('Task activation score')
xlabel('The Columbia Neurological Score')
% ylabel('CNS score')
h=lsline
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)
hold on
scatter(Allscores(find(GroupID==0),2),stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)),dotsize,dotcolor2,'.')

subplot(2,2,4)
scatter(Allscores(find(GroupID==0),2),stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)),dotsize,dotcolor2,'.')

% ylabel('Mean activation in significant activated region')
ylabel('Task activation score')
xlabel('The Columbia Neurological Score')
% ylabel('CNS score')
h=lsline
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
set(h(1),'color','#FF594C','linewidth',2)


X1=Allscores(:,1);
X2=Allscores(:,2);
Y=stats_SVM.paired_hyperplane_dist_scores;


[RR(1) PP(1)]=corr(X1,Y,'rows','complete','type','Spearman');
[RR(2) PP(2)]=corr(X2,Y,'rows','complete','type','Spearman');
[RR(3) PP(3)]=corr(X1(find(GroupID==0)),Y(find(GroupID==0)),'rows','complete','type','Spearman');
[RR(4) PP(4)]=corr(X2(find(GroupID==0)),Y(find(GroupID==0)),'rows','complete','type','Spearman');

index=find(GroupID==0);
% index=1:size(X1,1);
Datamatrix=[X1(index) X2(index) Y(index) sex(index)' age(index)' GroupID(index)' Mito_Sub(index)']
NameCell={'NMDAS','CNS','Brain score','sex','age','MitoD_Or_Not','MitoSubgroup'}
for i=1:7
    for j=1:7
        [R P]=corr(Datamatrix(:,i),Datamatrix(:,j),'type','Spearman');
        Rmatrix(i,j)=R;
        Pmatrix(i,j)=P;
    end
end


figure
imagesc(abs(Rmatrix),[0 0.8])      
colorbar
xticklabels(NameCell')
yticklabels(NameCell')
% mdl1 = fitglm([X1(find(GroupID==0)) age(find(GroupID==0))'],Y(find(GroupID==0)))
% 
% mdl1 = fitglm([X1 zscore(age)'],Y)

[RR(1) PP(1)]=partialcorr(Allscores(:,1),stats_SVM.paired_hyperplane_dist_scores,age','rows','complete','type','Spearman');
[RR(2) PP(2)]=partialcorr(Allscores(:,2),stats_SVM.paired_hyperplane_dist_scores,age','rows','complete','type','Spearman');
[RR(3) PP(3)]=partialcorr(Allscores(find(GroupID==0),1),stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)),[age(find(GroupID==0))' ],'rows','complete','type','Spearman');
[RR(4) PP(4)]=partialcorr(Allscores(find(GroupID==0),2),stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)),[age(find(GroupID==0))' ],'rows','complete','type','Spearman');



[a b c]=xlsread('LOCAL_SOURCE\Dropbox (Dartmouth College)\2023_Bo_Ke_Picard_MITO_brain_analysis_MISBIE\MiSBIE Data (Eprime and Questionnaires)\MiSBIE MRI Meta Data 6-11-24.xlsm')
cellArray1=metaGrand;
cellArray2=b;
ids1 = cellArray1(:, 1);
ids2 = cellArray2(2:end, 1);
% Find the indices of the IDs in cellArray1 that are also in cellArray2
[~, idx] = ismember(ids1, ids2);
Behavioral_ACC=a(idx,66);
Behavioral_ACC_Index=a(idx,67);
Behavioral_ACC_Index=find(Behavioral_ACC_Index==1);
[R P]=partialcorr(GroupID',stats_SVM.paired_hyperplane_dist_scores,Behavioral_ACC,'rows','complete','type','Spearman')
[R P]=corr(stats_SVM.paired_hyperplane_dist_scores(:),Behavioral_ACC(:),'rows','complete','type','Spearman')
[R P]=corr(Allscores(Behavioral_ACC_Index,1),Behavioral_ACC(Behavioral_ACC_Index),'rows','complete','type','Spearman')
[R P]=corr(stats_SVM.paired_hyperplane_dist_scores(Behavioral_ACC_Index),GroupID(Behavioral_ACC_Index)','rows','complete','type','Spearman')


[R P]=partialcorr(Behavioral_ACC,GroupID',stats_SVM.paired_hyperplane_dist_scores,'rows','complete','type','Spearman')

[R P]=corr(Allscores(:,1),Behavioral_ACC(:),'rows','complete','type','Spearman')


figure
scatter(stats_SVM.paired_hyperplane_dist_scores,Behavioral_ACC*100,dotsize,dotcolor1,'.')
h=lsline
% ylabel('Mean activation in significant activated region')
xlabel('Brain activation score')
ylabel('N-back task accuracy (%)')
hold on
scatter(stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)),Behavioral_ACC(GroupID==0)*100,dotsize,dotcolor2,'.')
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)


dotcolor1=[3 141 204]/256;
dotcolor2=[0.8 0.3450 0.4330];

figure
scatter(Allscores(:,1),Behavioral_ACC*100,dotsize,dotcolor1,'.')
h=lsline
hold on
scatter(Allscores(GroupID==0,1),Behavioral_ACC(find(GroupID==0))*100,dotsize,dotcolor2,'.')
xlabel('NMDAS')
ylabel('N-back task accuracy (%)')
set(gca,'fontsize',10,'fontweight','bold','LineWidth',1)
% [R P]=corr(Behavioral_ACC(GroupID==0),stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)),'rows','complete','type','Spearman')
% Filter out the IDs that exist in both cell arrays
selectedIDs = ids1(idx > 0);

% for i=1:5
% end
X=GroupID(Behavioral_ACC_Index)';
% X=Allscores(Behavioral_ACC_Index,1);
Y=Behavioral_ACC(Behavioral_ACC_Index);
M=stats_SVM.paired_hyperplane_dist_scores(Behavioral_ACC_Index);

X=GroupID(:);
% X=Allscores(Behavioral_ACC_Index,1);
Y=Behavioral_ACC(:);
M=stats_SVM.paired_hyperplane_dist_scores(:);

[R P]=corr(X,Y,'type','spearman')
[R P]=corr(X,Y,'type','spearman')
[R P]=partialcorr(X,Y,M,'type','spearman')

predictors = [X M];
model = fitlm(predictors, Y)

predictors = [X];
model = fitlm(predictors, Y)

predictors = [Y M];
model = fitlm(predictors, X);
disp(model);

% Behavioral_ACC,GroupID',stats_SVM.paired_hyperplane_dist_scores
Behavioral_ACC_Index;
% Y=Behavioral_ACC(Behavioral_ACC_Index);

[paths, stats1, stats2] = mediation(X, Y, M, 'boottop', 'stats', 'plots');
% The residuals represent X with the effect of Z removed
X_adjusted = residuals;

Y=Allscores(:,1);
M=stats_SVM.paired_hyperplane_dist_scores;

[paths, stats1, stats2] = mediation(X1, Y, M, 'boottop', 'stats', 'plots');

X1=GroupID';
Y=stats_SVM.paired_hyperplane_dist_scores;
M=Allscores(:,1);

[paths, stats1, stats2] = mediation(X1, Y, M, 'boottop', 'stats', 'plots');


X1=GroupID';
Y=stats_SVM.paired_hyperplane_dist_scores;
M=age';

[paths, stats1, stats2] = mediation(X1, Y, M, 'boottop', 'stats', 'plots');

figure
scatter(age',Y,dotsize,dotcolor1,'.')
[R P]=corr(age',Y)

scatter(age',Y,dotsize,dotcolor1,'.')
h=lsline
% ylabel('Mean activation in significant activated region')
ylabel('Age')
xlabel('Brain activation score')
hold on
scatter(age(GroupID==0),Y(find(GroupID==0)),dotsize,dotcolor2,'.')


figure
scatter(age',X1,dotsize,dotcolor1,'.')
[R P]=corr(age',X1)

scatter(age',X1,dotsize,dotcolor1,'.')
h=lsline
% ylabel('Task acivation score')
xlabel('Age')
ylabel('NMDAS')
hold on
scatter(age(GroupID==0),X1(find(GroupID==0)),dotsize,dotcolor2,'.')


figure
[h p ci stat]=ttest2(stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)),stats_SVM.paired_hyperplane_dist_scores(find(GroupID==1)))

% vector1 = stats_SVM.paired_hyperplane_dist_scores(find(GroupID==0)); % Example vector 1
% vector2 = stats_SVM.paired_hyperplane_dist_scores(find(GroupID==1)); % Example vector 2
vector1 = Behavioral_ACC(find(GroupID==0)); % Example vector 1
vector2 = Behavioral_ACC(find(GroupID==1)); % Example vector 2
vector1 = vector1(~isnan(vector1)); % Remove NaN values
vector2 = vector2(~isnan(vector2)); % Remove NaN values
% Calculate mean and standard error for each vector
mean1 = mean(vector1)*100;
mean2 = mean(vector2)*100;

stderr1 = std(vector1) / sqrt(length(vector1))*100;
stderr2 = std(vector2) / sqrt(length(vector2))*100;

% Means and standard errors
means = [mean1, mean2];
stderrs = [stderr1, stderr2];

% Create bar chart
figure;
bar(means);
hold on;

% Add error bars
errorbar(means, stderrs, 'k', 'linestyle', 'none');

% Customize plot
set(gca, 'XTickLabel', {'MitoD', 'Control'});
ylabel('N-back task accuracy (%)');
title('N-back task accuracy compared across group');
 