% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
function results = kmeans_bootstrap_stability(patients_severity, B)
% K-means bootstrap stability for a 1D severity vector.
% patients_severity: 1xN or Nx1 numeric vector
% B: number of bootstrap iterations (e.g., 2000)

if nargin < 2, B = 2000; end
x = patients_severity(:);
N = numel(x);
k = 2;

rng(42); % reproducible

% ==== 1) Fit the reference (original) solution ====
[idx_ref, C_ref] = kmeans_Matlab(x, k, 'Replicates', 100, 'MaxIter', 1000, ...
    'OnlinePhase','off', 'Distance','sqeuclidean');

% Helper to predict labels for any set using learned centers:
predict_labels = @(X, C) assign_to_nearest(X, C);

% ==== 2) Storage for bootstrap metrics ====
ARI = zeros(B,1);                        % Agreement with reference on ALL N
consensus = zeros(N, N);                 % Co-association counts
consensus_den = zeros(N, N);             % Pairwise “seen together” counts
assign_stability_counts = zeros(N,1);    % Out-of-bag assignment matches ref
assign_stability_den = zeros(N,1);       % Times each subject was OOB
centers_store = zeros(B,k);              % Bootstrap centers (sorted)

% ==== 3) Bootstrap loop ====
for b = 1:B
    % Sample with replacement
    idx_inbag = randsample(N, N, true);
    Xb = x(idx_inbag);

    % Fit k-means on the in-bag sample
    [idx_b_inbag, Cb] = kmeans_Matlab(Xb, k, 'Replicates', 50, 'MaxIter', 1000, ...
        'OnlinePhase','off', 'Distance','sqeuclidean');

    % Sort centers to make labels comparable over iterations
    [Cb, order] = sort(Cb(:));
    idx_b_inbag = relabel_by_center(idx_b_inbag, order);

    centers_store(b,:) = Cb(:).';

    % ---- (A) Agreement with reference on FULL SET (using prediction) ----
    idx_pred_all = predict_labels(x, Cb);
    idx_pred_all = align_to_reference(idx_pred_all, idx_ref);
    ARI(b) = adjustedRandIndex(idx_pred_all, idx_ref);

    % ---- (B) Update consensus using only in-bag subjects ----
    % Map in-bag labels back to original indices (ties handled by first occurrence)
    % For duplicates, keep their first assigned label
    first_occ = ~duplicated_flag(idx_inbag);
    ib_unique = idx_inbag(first_occ);
    lab_unique = idx_b_inbag(first_occ);

    % Build an indicator vector for each cluster among in-bag unique indices
    for c = 1:k
        members = ib_unique(lab_unique == c);
        consensus(members, members) = consensus(members, members) + 1;
    end
    % Count times pairs were seen together
    consensus_den(ib_unique, ib_unique) = consensus_den(ib_unique, ib_unique) + 1;

    % ---- (C) Per-subject assignment stability (OOB prediction) ----
    oob_mask = true(N,1);
    oob_mask(idx_inbag) = false; % those not sampled this round
    if any(oob_mask)
        idx_pred_oob = predict_labels(x(oob_mask), Cb);
        idx_pred_oob = align_to_reference_outofbag(idx_pred_oob, idx_ref(oob_mask));
        assign_stability_counts(oob_mask) = assign_stability_counts(oob_mask) + ...
            (idx_pred_oob == idx_ref(oob_mask));
        assign_stability_den(oob_mask) = assign_stability_den(oob_mask) + 1;
    end
end

% ==== 4) Aggregate metrics ====
consensus_prob = consensus ./ max(consensus_den, 1); % avoid divide-by-zero
mean_ARI = mean(ARI);
ci_ARI = quantile(ARI, [0.025 0.975]);

% Per-subject stability (only defined where subject was OOB at least once)
per_subject_stability = nan(N,1);
has_oob = assign_stability_den > 0;
per_subject_stability(has_oob) = assign_stability_counts(has_oob) ./ assign_stability_den(has_oob);

% Bootstrap CIs for centers and decision threshold
centers_sorted = sort(centers_store, 2);
c1_ci = quantile(centers_sorted(:,1), [0.025 0.5 0.975]);
c2_ci = quantile(centers_sorted(:,2), [0.025 0.5 0.975]);
thresh_boot = mean(centers_sorted, 2);  % midpoint as decision boundary in 1D
thresh_ci = quantile(thresh_boot, [0.025 0.5 0.975]);

% Optional: derive a “consensus partition” via clustering the consensus matrix
% Here we just do hierarchical clustering on 1 - consensus_prob
D = 1 - consensus_prob;
D(1:N+1:end) = 0;
Z = linkage(squareform((D + D')/2, 'tovector'), 'average');
consensus_idx = cluster(Z, 'maxclust', k);
consensus_idx = align_to_reference(consensus_idx, idx_ref);

% ==== 5) Package outputs ====
results.idx_ref = idx_ref;
results.C_ref = sort(C_ref);
results.mean_ARI = mean_ARI;
results.ci_ARI = ci_ARI;
results.per_subject_stability = per_subject_stability;   % Nx1
results.consensus_prob = consensus_prob;                 % NxN
results.consensus_partition = consensus_idx;             % Nx1
results.center_ci = struct('c1',[c1_ci(1) c1_ci(2) c1_ci(3)], ...
                           'c2',[c2_ci(1) c2_ci(2) c2_ci(3)]);
results.threshold_ci = [thresh_ci(1) thresh_ci(2) thresh_ci(3)];
results.centers_store = centers_store;                   % Bx2
results.ARI = ARI;                                       % Bx1

% ====== Nested helpers ======
function y = assign_to_nearest(X, C)
    % Assign each X to nearest center C in 1D (k x 1)
    D = abs(X - C(:)');   % N x k distances
    [~, y] = min(D, [], 2);
end

function y = relabel_by_center(labels, order)
    % If order = [2;1], swap labels 1<->2, etc., so label numbers follow
    % ascending centers.
    map = zeros(1, numel(order));
    for i = 1:numel(order)
        map(order(i)) = i;
    end
    y = map(labels).';
end

function y = align_to_reference(y, ref)
    % For k=2, flip labels if it increases agreement with ref
    agree = mean(y == ref);
    agree_flipped = mean((3 - y) == ref); % 1<->2
    if agree_flipped > agree
        y = 3 - y;
    end
end

function y = align_to_reference_outofbag(y, ref_oob)
    % Same as above but on the OOB subset
    agree = mean(y == ref_oob);
    agree_flipped = mean((3 - y) == ref_oob);
    if agree_flipped > agree
        y = 3 - y;
    end
end

function tf = duplicated_flag(v)
    % True for duplicated elements after the first occurrence
    [~, firstIdx] = unique(v, 'stable');
    tf = true(size(v));
    tf(firstIdx) = false;
end

function ari = adjustedRandIndex(u, v)
    % u, v are Nx1 integer labels (1..k)
    % Computes the Hubert & Arabie ARI.
    u = u(:); v = v(:);
    n = numel(u);
    % Contingency table
    [~,~,uc] = unique(u);
    [~,~,vc] = unique(v);
    k1 = max(uc); k2 = max(vc);
    nij = accumarray([uc vc], 1, [k1 k2]);
    ai = sum(nij,2);
    bj = sum(nij,1);
    comb2 = @(x) x.*(x-1)/2;
    sum_nij = sum(comb2(nij(:)));
    sum_ai = sum(comb2(ai));
    sum_bj = sum(comb2(bj));
    denom = 0.5*(sum_ai + sum_bj);
    if denom == 0
        ari = 0; % degenerate case
    else
        ari = (sum_nij - denom) / (0.5*(sum_ai + sum_bj) - denom);
    end
end
end
