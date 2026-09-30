% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
clear
load('E:\Mito_DICOM\SecondLevelSave\Thresholded\Multisensory_100.mat')
D_All(1,:)=D_Value;
Acc_All(1,:)=Acc;
k=1;
for i=1:90
    if strcmp(metaGrand(i,2),'Control')==1
        GroupID(k)=1;
        k=k+1;
    else
        GroupID(k)=0;
        k=k+1;
    end
end
Dist_Brain{1}=mean(dist_Perm,1);
GroupID_Cell{1}=GroupID;
load('E:\Mito_DICOM\SecondLevelSave\Thresholded\Cold_100.mat')
D_All(2,:)=D_Value;
Acc_All(2,:)=Acc;
GroupID_Cell{2}=GroupID;
Dist_Brain{2}=mean(dist_Perm,1);
load('E:\Mito_DICOM\SecondLevelSave\Thresholded\Stress_100.mat')
D_All(3,:)=D_Value;
Acc_All(3,:)=Acc;
GroupID_Cell{3}=GroupID;
Dist_Brain{3}=mean(dist_Perm,1);
load('E:\Mito_DICOM\SecondLevelSave\Thresholded\Nback_100.mat')
D_All(4,:)=D_Value;
Acc_All(4,:)=Acc;
GroupID_Cell{4}=GroupID;
Dist_Brain{4}=mean(dist_Perm,1);
D_All=Acc_All;
colorcoding=[110 203 99; 0 176 240;255 89 76; 1 34 23];
colorcoding=colorcoding/255;
LineWidth=2;
Fontsize=15;
colorcoding(4,:)=[];
D_All(3,:)=[];
figure
violinplot(D_All','facecolor',colorcoding,'mc','k','bw',0.001,'plotlegend',0,'pointsize',2)
set(gca,'xtick',[])
names={'Multisensory';'Cold Pain';'Social Stress';'N back'}
names(3)=[];
set(gca,'xtick',[1:4],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
box off
ylabel('Decoding accuracy (%)')

colorcoding=[110 203 99; 0 176 240;255 89 76; 1 34 23];
colorcoding=colorcoding/255;
LineWidth=2;
Fontsize=20;
figure
violinplot(Acc_All','facecolor',colorcoding,'mc','k','bw',0.0004,'plotlegend',0,'pointsize',2)
set(gca,'xtick',[])
names={'Multisensory';'Cold Pain';'Social Stress';'N back'}
set(gca,'xtick',[1:4],'xticklabel',names)
set(gca,'xticklabel',names)
set(gca,'linewidth',LineWidth,'Fontsize',Fontsize,'FontWeight','bold')
box off


colorcoding=[ 255 89 76 ;110 203 99];
colorcoding=colorcoding/255;
LineWidth=2;
Fontsize=20;
figure
subplot(2,2,1)
GroupID=GroupID_Cell{1};
violinplot({Dist_Brain{1}(find(GroupID==0)) Dist_Brain{1}(find(GroupID==1))},'facecolor',colorcoding,'mc','k','bw',0.4,'plotlegend',0,'pointsize',2)
ylabel('Task brain score')
title('Multisensory')

LineWidth=2;
Fontsize=20;
subplot(2,2,2)
GroupID=GroupID_Cell{2};
violinplot({Dist_Brain{2}(find(GroupID==0)) Dist_Brain{2}(find(GroupID==1))},'facecolor',colorcoding,'mc','k','bw',0.4,'plotlegend',0,'pointsize',2)
ylabel('Task brain score')
title('Cold pain')

LineWidth=2;
Fontsize=20;
subplot(2,2,3)
GroupID=GroupID_Cell{3};
violinplot({Dist_Brain{3}(find(GroupID==0)) Dist_Brain{3}(find(GroupID==1))},'facecolor',colorcoding,'mc','k','bw',0.4,'plotlegend',0,'pointsize',2)
ylabel('Task brain score')
title('Social Stress')

LineWidth=2;
Fontsize=20;
subplot(2,2,4)
GroupID=GroupID_Cell{4};
violinplot({Dist_Brain{4}(find(GroupID==0)) Dist_Brain{4}(find(GroupID==1))},'facecolor',colorcoding,'mc','k','bw',0.4,'plotlegend',0,'pointsize',2)
ylabel('Task brain score')
title('N-back')

%%%%%%%%% figure 2D%%%%%%%%%%
colorcoding2=[110 203 99; 0 176 240;255 89 76; 1 34 23]/255;
figure
violinplot({Acc_All(1,:)*100 Acc_All(2,:)*100 Acc_All(4,:)*100 },'facecolor',colorcoding2,'mc','k','bw',0.5,'plotlegend',0,'pointsize',2)
axis([0.5 3.5 40 110])
ylabel('Task decoding accuracy (%)')
set(gca,'fontsize',15,'fontweight','bold')
yticklabels(40:10:100)
xticklabels({'Multisensory','Cold pain','N-back'})
