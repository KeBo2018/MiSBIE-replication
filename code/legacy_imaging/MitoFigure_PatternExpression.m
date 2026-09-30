% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.

load('E:\Mito_DICOM\SecondLevelSave\Thresholded\Nback_100_500Filter.mat')


datamat=mean(double_dist_Perm,1);
datamatA(1:88,1)=datamat(1:88);
datamatA(1:88,2)=datamat(89:end);
figure;
hold on;

for i = 1:size(D, 1)
    % Draw line between the two conditions
    plot([subjectID(i), subjectID(i)], [D(i,1), D(i,2)], 'k-','LineWidth', 1.5); 

    % Control condition (e.g., blue filled circle with black border)
    plot(subjectID(i), D(i,1), 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);

    % Task condition (e.g., red filled circle with black border)
    plot(subjectID(i), D(i,2), 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
end

xlabel('SubjectID');
ylabel('N-back Pattern expression');
ylim([-6 6]);
box on;

% Optional: prettier legend using dummy plots
h1 = plot(nan, nan, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
h2 = plot(nan, nan, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
legend([h1 h2], {'Task state', 'Control state'}, 'Location', 'northeast');
set(gca,'fontsize',15,'LineWidth', 1.5)

%%%% Schematic plot for condition 1 and two
Condition1_Sub1=get_wh_image(Condition1,18);
Condition2_Sub1=get_wh_image(Condition2,18);

figure

montage(Condition1_Sub1,'full')

figure
montage(Condition2_Sub1,'full')