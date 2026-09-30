% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
% Load connectivity and metadata 
clear
load('F:\canlab2024_Labels.mat')                % atlas_obj.labels_5: 518×1 cell array of 'network1','network2',...
RestingstateNetwork = atlas_obj;
load('E:\\Mito_DICOM\\SecondLevelSave\\Canlab2024_NetFC.mat');   % Rvalue: 518×518×nSubj
load('E:\Mito_DICOM\SecondLevelSave\MitoRest_FD.mat')           % meanFD: 1×nSubj
T = readtable('F:\\Mito_Rest\\Meta_data\\Misbie_Meta_Data.xlsx');
[a, b, c] = xlsread('F:\\Mito_Rest\\Meta_data\\Misbie_Meta_Data.xlsx');

%% Define groups
PatientIndex = strcmp(b(2:end,2), '''Control''');
GroupID      = double(~PatientIndex);  % 0=Control, 1=Patient
Control      = find(PatientIndex);
Patient      = find(~PatientIndex);

nSubj = size(Rvalue,3);

%% Build network‐level connectivity (method 1)
% atlas_obj.labels_5 is 518×1 cell of 'network1','network2',...
[netLabels, ~, netIdx] = unique(atlas_obj.labels_5);
nNet = numel(netLabels);

% Preallocate: netconn(i,j,s) = mean edge between network i and j in subject s
netconn = zeros(nNet,nNet,nSubj);
for s = 1:nSubj
    for ni = 1:nNet
        idx_i = find(netIdx == ni);
        for nj = ni:nNet
            idx_j = find(netIdx == nj);
            block = Rvalue(idx_i, idx_j, s);
            if ni==nj
                % exclude the diagonal within a network
                mask = ~eye(numel(idx_i));
                block = block(mask);
            end
            netconn(ni,nj,s) = mean(block(:));
            netconn(nj,ni,s) = netconn(ni,nj,s);
        end
    end
end

%% 1) Group differences on network‐level edges
% Preallocate results
bval_net  = zeros(nNet);
tval_net  = zeros(nNet);
pval_net  = ones(nNet);

for i = 1:nNet
    for j = 1:nNet
        if i==j, continue; end
        conn = squeeze(netconn(i,j,:));       % [nSubj×1]
        % ensure both predictors are column vectors of length nSubj
        grp = GroupID(:);       % nSubj×1
        fd  = meanFD(:);        % nSubj×1
        X   = [grp, fd];        % nSubj×2
        [b,stats] = robustfit(X, conn);
        bval_net(i,j) = b(2);   % group effect
        tval_net(i,j) = stats.t(2);
        pval_net(i,j) = stats.p(2);
    end
end
pval_net(1:nNet+1:end) = NaN;  % diagonal→NaN

% FDR on upper triangle
upperIdx = find(triu(ones(nNet),1));
pv       = pval_net(upperIdx);
adjpv    = mafdr(pv,'BHFDR',true);

adj_pval_net = nan(nNet);
adj_pval_net(upperIdx) = adjpv;
adj_pval_net = adj_pval_net + adj_pval_net';

sig_net_mask = adj_pval_net < 0.05;
if ~any(sig_net_mask(:))
    sig_net_mask = pval_net < 0.001;
end

% Plot
figure;
imagesc(tval_net*-1,[-4 4]); colorbar;
cm = colormap_tor([0 0 1],[1 0 0],[1 1 1]); colormap(cm);
set(gca, ...
    'XTick',1:nNet,'XTickLabel',netLabels, ...
    'YTick',1:nNet,'YTickLabel',netLabels, ...
    'XTickLabelRotation',45,'FontSize',8);
axis square; title('Control vs Mito  Group Differences');
hold on;
[yS,xS] = find(sig_net_mask);
scatter(xS,yS,100,'ks');
hold off;

%% 2) Correlate network‐level edges with NMDAS in patients
NMDAS     = a(:,4);
NMDAS_pat = NMDAS(Patient);
FD_pat    = meanFD(Patient)';

bval_corr_net = zeros(nNet);
tval_corr_net = zeros(nNet);
pval_corr_net = ones(nNet);

for i = 1:nNet
    for j = 1:nNet
        if i==j, continue; end
        conn_pat = squeeze(netconn(i,j,Patient));  % [nPat×1]
        % ensure predictors are column vectors of length nPat
        nmd = NMDAS_pat(:);  % nPat×1
        fdp = FD_pat(:);     % nPat×1
        X   = [nmd, fdp];    % nPat×2
        [b,stats] = robustfit(X, conn_pat);
        bval_corr_net(i,j) = b(2);
        tval_corr_net(i,j) = stats.t(2);
        pval_corr_net(i,j) = stats.p(2);
    end
end
pval_corr_net(1:nNet+1:end) = NaN;

% FDR correction
upperIdx2 = find(triu(ones(nNet),1));
pv2       = pval_corr_net(upperIdx2);
adjpv2    = mafdr(pv2,'BHFDR',true);

adj_pval_corr_net = nan(nNet);
adj_pval_corr_net(upperIdx2) = adjpv2;
adj_pval_corr_net = adj_pval_corr_net + adj_pval_corr_net';

sig_corr_mask = adj_pval_corr_net < 0.05;
if ~any(sig_corr_mask(:))
    sig_corr_mask = pval_corr_net < 0.001;
end

% Plot
figure;
imagesc(tval_corr_net,[-4 4]); colorbar;
colormap(cm);
set(gca, ...
    'XTick',1:nNet,'XTickLabel',netLabels, ...
    'YTick',1:nNet,'YTickLabel',netLabels, ...
    'XTickLabelRotation',45,'FontSize',8);
axis square; title('NMDAS–Network Connectivity');
hold on;
[yC,xC] = find(sig_corr_mask);
scatter(xC,yC,100,'ks');
hold off;
