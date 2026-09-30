% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
% clear 
% %%%%%%%%%%%%%%%%
 load('E:\Mito_DICOM\SecondLevelSave\Nback_TwoRun_Right_order_Filter_PercentageChange_SPM_Robust.mat')
Contrast=fmri_data(Contrast);
% T=ttest(image_math(fmri_data(Twoback),fmri_data(Zeroback),'minus'))

%%% T-test for Contrast and FDR correction %%% 
T=ttest(fmri_data(Contrast))
T_Threshold=threshold(T,0.05, 'fdr');
montage(T_Threshold,'full')
figure
surface(T_Threshold)
Contrast.removed_images=zeros(1,88);


WM_SingleSub=image_math(fmri_data(Twoback),fmri_data(Zeroback),'minus')
WM_SingleSub.removed_images=[];

Y=[zeros(88,1)]
xval_SVM(Twoback.dat,Zeroback.dat)

A=get_wh_image(Minus,2)
figure
A.removed_images=[];
threshold(A)
orthviews(A)
figure
 surface(WM_SingleSub);

%%%%  Load disease severity score %%%%
NMDAS_Score=cell2mat(metaGrand(:,6));
CNS_Score=cell2mat(metaGrand(:,7));
%%%%%%%%%%%%%%%

%%%% Use robust regression and regress out white matter and CSF %%%%%
Robust_Contrast=robfit_parcelwise(Contrast,'csf_wm_covs',1)


%%%% Plot two groups saperetely %%%%
Control_Contrast=Contrast;
Control_Contrast.dat=Control_Contrast.dat(:,find(GroupID==1))
T1=ttest(Control_Contrast)
Control_Contrast_T_Thre=threshold(T1,0.05, 'fdr');
figure
montage(Control_Contrast_T_Thre)

GroupID(find(GroupID==0))=-1;
Mito_Contrast=Contrast;
Mito_Contrast.dat=Mito_Contrast.dat(:,find(GroupID==-1))

T2=ttest(Mito_Contrast)
orthviews(T2)
% Mito_Contrast_T_Thre=threshold(T2,0.05, 'fdr');
Mito_Contrast_T_Thre=threshold(T2,[-2 2], 'raw-outside');
figure
montage(Mito_Contrast_T_Thre)

%%%% Store group information (Or disease score) as second level regressor %%%%
Contrast.X=GroupID;

%%% Get rid of NAN %%%%%%
Contrast_NanF=get_wh_image(Contrast,find(isnan(GroupID(:,3))==0));
GroupCon=regress(Contrast)

GroupCon_T=get_wh_image(GroupCon.t,1);
orthviews(GroupCon_T)
figure
GroupCon_T_Thre=threshold(GroupCon_T,0.05, 'fdr');
GroupCon_T_Thre=threshold(GroupCon_T,[-1 1], 'raw-outside');
figure
montage(GroupCon_T_Thre,'full')
figure
surface(GroupCon_T_Thre)

figure
montage(GroupCon_T_Thre,'full')
ttest2(Control_Contrast,Mito_Contrast)


canlabAtlas=load_atlas('canlab2018');

V1 = select_atlas_subset(canlabAtlas, {'Ctx_V1'});
M_Insula = select_atlas_subset(canlabAtlas, {'MI'});
P_Insula = select_atlas_subset(canlabAtlas, {'PoI'});
Thal = select_atlas_subset(canlabAtlas, {'Thal'});
S1 = select_atlas_subset(canlabAtlas, {'Ctx_3a_', 'Ctx_3b_', 'Ctx_1_', 'Ctx_2_'});
S2 = select_atlas_subset(canlabAtlas, {'Ctx_OP1'}); 

S1Activation=apply_mask(Contrast,S1)
S1A=mean(S1Activation.dat,1)

S2Activation=apply_mask(Contrast,S2)
S2A=mean(S2Activation.dat,1)

P_Insula_Activation=apply_mask(Contrast,P_Insula)
P_Insula_A=mean(P_Insula_Activation.dat,1)

Thal_Activation=apply_mask(Contrast,Thal)
Thal_A=mean(Thal_Activation.dat,1)

Activation=[S1A; S2A; P_Insula_A; Thal_A] 
figure
violinplot(Activation','mc','k','bw',1,'plotlegend',0,'pointsize',2)

    colorcoding=[110 203 99; 0 176 240;255 89 76];
colorcoding=colorcoding/255;
G=[parcel_means_Neu(:,1) parcel_means_Neg(:,1) parcel_means_Reg(:,1)]'
SubjectMean=mean(G,1);
Grandmean=mean(mean(G));
G_within=G-ones(3,1)*SubjectMean+Grandmean*ones(3,size(G,2));

violinplot(G_within','facecolor',colorcoding(:,:),'mc','k','bw',0.08,'plotlegend',0,'pointsize',2)

%%%%
keywords = {'Ctx_3a_', 'Ctx_3b_', 'Ctx_1_', 'Ctx_2_'};
keywords = {'MI'}
% Initialize an empty vector to store indices
indices = [];

% Loop through the keywords
for i = 1:length(keywords)
    % Find indices where the keyword is found in cellsArray
    foundIndices = find(contains(canlabAtlas.labels, keywords{i}));
    
    % Append the found indices to the indices vector
    indices = [indices; foundIndices];
end
S1_Robust_L=Robust_Contrast.tscores(indices(:,1),4);
S1_Robust_R=Robust_Contrast.tscores(indices(:,2),4);

% Since the same index might be found for different keywords, remove duplicates
indices = unique(indices);

% Display the indices
disp(indices);




[group_metrics, individual_metrics, global_gm_wm_csf_values] = qc_metrics_second_level(Contrast)
Contrast.X(:,2:3)=[ global_gm_wm_csf_values(:,2:3)];
Contrast.X(:,1)=CNS_Score';
Contrast_NanF=get_wh_image(Contrast,find(isnan(CNS_Score)==0));
out=regress(Contrast_NanF)
GroupCon_T=get_wh_image(out.t,1);
GroupCon_T_Thre=threshold(GroupCon_T,[-2 2], 'raw-outside');
figure
montage(GroupCon_T_Thre,'full')



Contrast.X(:,1:2)=[ global_gm_wm_csf_values(:,2:3)];
Contrast.X=Contrast.X';
out=regress(Contrast);

GroupCon=regress(Contrast,'residual','nointercept')
GroupCon_T=GroupCon.resid;

GroupCon_T=ttest(GroupCon.resid);
orthviews(GroupCon_T)
figure
GroupCon_T_Thre=threshold(GroupCon_T,0.05, 'fdr');

GroupCon_T=get_wh_image(out.t,1);
figure
orthviews(GroupCon_T)
GroupCon_T_Thre=threshold(GroupCon_T,[-2 2], 'raw-outside');
figure
montage(GroupCon_T_Thre,'full')

figure

corr(GroupCon_T.dat,T.dat)
figure
scatter(T.dat,GroupCon_T.dat)
h=lsline
xlabel('Before 2nd level regression')
ylabel('After 2nd level regression ')
%%%%%%%%%%%% ROI analysis %%%%

S1Activation=apply_mask(GroupCon.resid,S1)
S1A=mean(S1Activation.dat,1)

S2Activation=apply_mask(GroupCon.resid,S2)
S2A=mean(S2Activation.dat,1)

P_Insula_Activation=apply_mask(GroupCon.resid,P_Insula)
P_Insula_A=mean(P_Insula_Activation.dat,1)

Thal_Activation=apply_mask(GroupCon.resid,Thal)
Thal_A=mean(Thal_Activation.dat,1)

Activation=[S1A; S2A; P_Insula_A; Thal_A] 
    colorcoding=[110 203 99; 0 176 240];
colorcoding=colorcoding/255;
figure
violinplot(Activation','mc','k','bw',1,'plotlegend',0,'pointsize',2)
ylabel('Activation beta')
ControlID=find(GroupID(:,2)==1);
MitoID=find(GroupID(:,2)==0);
figure
for i=1:4
    subplot(2,2,i)
    violinplot({Activation(i,ControlID);Activation(i,MitoID)}','facecolor',colorcoding(:,:),'mc','k','bw',1,'plotlegend',0,'pointsize',2)
    [h p ci tstat]=ttest2(Activation(i,ControlID),Activation(i,MitoID));
    pval(i)=p;
    ylabel('Activation beta')
end


Robust_Meta=get_wh_image(Robust_Contrast.t_obj,1);
figure
montage(Robust_Meta,'full')

Contrast.X=GroupID(:,2);

Contrast_NanF=get_wh_image(Contrast,find(isnan(GroupID(:,2))==0));
global_gm_wm_csf_values_NanF=global_gm_wm_csf_values(find(isnan(GroupID(:,2))==0),:);

Contrast_NanF.X=[Contrast_NanF.X global_gm_wm_csf_values_NanF(:,2:3)]
% GroupCon=regress(Contrast_NanF)
GroupCon=regress(Contrast)
GroupCon_T=get_wh_image(GroupCon.t,4);
orthviews(GroupCon_T)
% GroupCon_T_Thre=threshold(GroupCon_T,[-2 2], 'raw-outside');
GroupCon_T_Thre=threshold(GroupCon_T,0.05, 'fdr');
figure
montage(GroupCon_T_Thre,'full')

Contrast.X=GroupID(:,2);

MPA2=apply_multiaversive_mpa2_patterns(GroupCon.resid);
MPA2_Val=table2array(MPA2);
figure
violinplot(MPA2_Val,'mc','k','bw',1,'plotlegend',0,'pointsize',2)
[h p]=ttest(MPA2_Val(:,5))

Allsigniture=apply_all_signatures(GroupCon.resid,'similarity_metric','cosine_similarity');
allVectors = [];
fields = fieldnames(Allsigniture);
for i = 1:length(fields)
    currentField = Allsigniture.(fields{i});
    % Check if the field is a table
    if istable(currentField)
        % Convert the table to a numeric array
        currentArray = table2array(currentField);
        % Ensure the array is a vector (here we simply reshape it to a single row)
        currentVector = currentArray(:).';
        % Append the vector to the matrix
        allVectors = [allVectors; currentVector]; % Appends as new rows
    end
end
figure
violinplot(allVectors(1:end,:)','xlabel',Allsigniture.signaturenames(1:end),'mc','k','bw',0.05,'plotlegend',0,'pointsize',2)
title('Cold pain signiture response (Cosine Similarity)')
grid on
for i=1:14
[h p ci stat]=ttest(allVectors(i,:))
pval(i)=p;
tval(i)=stat.tstat;
d(i)=tval(i)/sqrt(91);
end

GroupCon_Control=get_wh_image(GroupCon.resid,find(GroupID(:,2)==1));
GroupCon_Mito=get_wh_image(GroupCon.resid,find(GroupID(:,2)==-1));

Control_Sign=allVectors(:,find(GroupID(:,2)==1))
Mito_Sign=allVectors(:,find(GroupID(:,2)==-1))

for i=1:14
[h p ci stat]=ttest2(Mito_Sign(i,:),Control_Sign(i,:));
pval(i)=p;
tval(i)=stat.tstat;
end

obj = load_atlas('painpathways')
[obj.labels', obj.label_descriptions]

Step 2: Generate subsets of the atlas using the the keywords from all available label fields
dpins_obj = select_atlas_subset(obj, {'aMCC_MPFC'})
% brainstem_obj = select_atlas_subset(obj, {'Brainstem'}, 'labels_2', 'flatten')
% s1_roi = select_atlas_subset(obj, {'s2'})

Pain_related_ROI=apply_mask(GroupCon.resid,dpins_obj)
Pain_mean=mean(Pain_related_ROI.dat,1)
mean(Pain_mean)
[h p]=ttest(Pain_mean)


for i=1:24
    [obj_subset, to_extract] = select_atlas_subset(obj, i);
    Pain_related_ROI=apply_mask(GroupCon.resid,obj_subset);
    Pain_mean=mean(Pain_related_ROI.dat,1)
    [h p ci stat]=ttest(Pain_mean)
    pval(i)=p;
    tval(i)=stat.tstat;
    Grand(i,:)=Pain_mean;
   
%     violinplot(Pain_mean','xlabel',obj.labels(i),'mc','k','bw',1,'plotlegend',0,'pointsize',2)
%     hold on

end

figure
for i=1:4
subplot(2,2,i)
violinplot(Grand(i*6-5:i*6,:)','xlabel',obj.labels(i*6-5:i*6),'mc','k','bw',1,'plotlegend',0,'pointsize',2,'facecolor','b')
grid on
end

sigindex=find(pval<0.05)

for i=1:24
[R p]=corr(GroupID(:,3),Grand(i,:)','rows','complete');
Rvalue(i)=R;
PCorr(i)=p;
end

GroupID(3,:)=cell2mat(metaGrand(:,6));


[bucknermaps, networknames] = load_image_set('bucknerlab');
bucknermaps.image_names = char(networknames{:});
bucknermaps_R=resample_space(bucknermaps, GroupCon.resid)
for i=1:7
Restingstateindex{i}=find(bucknermaps_R.dat(:,i)==1)
end
for i=1:7
temp=GroupCon.resid.dat(find(bucknermaps_R.dat(:,i)==1),:);
Restingstate_Activation(i,:)=mean(temp,1);

[h p ci stat]=ttest(Restingstate_Activation(i,:));
pval_RN(i)=p;
tval_RN(i)=stat.tstat;

[R p]=corr(GroupID(4,:)',Restingstate_Activation(i,:)','rows','complete');
Rvalue_RN(i)=R;
PCorr_RN(i)=p;
end

figure
[R P]=corr(GroupID(4,:)',Restingstate_Activation(6,:)','rows','complete')
scatter(GroupID(4,:)',Restingstate_Activation(6,:)',50)
h=lsline
xlabel('Columbia Neurological score')
ylabel('Frontal parietal network activation')

figure
[R P]=corr(GroupID(4,:)',Restingstate_Activation(2,:)','rows','complete')
scatter(GroupID(4,:)',Restingstate_Activation(2,:)',50)
h=lsline
xlabel('Columbia Neurological score')
ylabel('Somatosensory network activation')

figure
[R P]=corr(GroupID(1,:)',Restingstate_Activation(6,:)','rows','complete')
scatter(GroupID(1,:)',Restingstate_Activation(6,:)',50)
h=lsline
xlabel('NMDS')
ylabel('Frontal parietal network activation')


figure
violinplot(Restingstate_Activation','xlabel',networknames' ,'mc','k','bw',1,'plotlegend',0,'pointsize',2,'facecolor','r')
grid on
corr(Restingstate_Activation)

    colorcoding=[110 203 99; 0 176 240];
colorcoding=colorcoding/255;
ControlID=find(GroupID(:,2)==1);
PatientID=find(GroupID(:,2)==0);
figure
for i=1:7
    subplot(2,4,i)

violinplot({Restingstate_Activation(i,ControlID);Restingstate_Activation(i,PatientID)}','facecolor',colorcoding,'mc','k','bw',1,'plotlegend',0,'pointsize',5);
title(networknames{i})
end
[h p]=ttest2(Restingstate_Activation(4,ControlID),Restingstate_Activation(4,PatientID))
figure
for i=1:24
    subplot(4,6,i)

violinplot({Grand(i,ControlID);Grand(i,PatientID)}','facecolor',colorcoding,'mc','k','bw',1,'plotlegend',0,'pointsize',5);
title(obj.labels(i))
end

allVectors(1:end,:)','xlabel',Allsigniture.signaturenames(1:end)
figure
for i=1:4
    subplot(2,2,i)

violinplot({allVectors(i,ControlID);allVectors(i,PatientID)}','facecolor',colorcoding,'mc','k','bw',0.1,'plotlegend',0,'pointsize',5);
title(Allsigniture.signaturenames(i))
end


%%% Test topic map %%
topic_file = which('neurosynth_topics_v4.mat');
which(topic_file) % error if not found - if error, add to path

load(topic_file)
topic_obj = resample_space(topic_obj_reverseinference, Contrast);

for i=1:50
index=find(topic_obj.dat(:,i)>0);
topic_activation(i,:)=mean(Contrast.dat(index,:),1);
[R P]=corr(topic_activation(i,:)',CNS_Score,'rows','complete','type','Spearman')
Rvalue(i,1)=R;
pvalue(i,1)=P;
end
figure

[R P]=corr(topic_activation(34,:)',NMDAS_Score,'rows','complete')
scatter(topic_activation(34,:)',NMDAS_Score,'filled','k')
h=lsline

xlabel('working memory topic map activation')
% ylabel('Columbia Neurological score')
ylabel('NMDAS score')
set(gca,'fontsize',15,'fontweight','bold','LineWidth',1)
WM_Map=get_wh_image(topic_obj,34)
figure
surface(WM_Map)