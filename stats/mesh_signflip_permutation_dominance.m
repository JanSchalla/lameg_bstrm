function result = mesh_signflip_permutation_dominance(surf1_vals, surf2_vals, Vertices, VertConn, radius_mm, n_perm, alpha, vertex_idx)
% MESH_SIGNFLIP_PERMUTATION_DOMINANCE  Tests whether the observed pattern of
% dominance (surf1 vs surf2) is distinguishable from what arbitrary
% surf1/surf2 labeling of the same paired values would produce, using a
% spatially-coherent sign-flip permutation test.
%
% NULL HYPOTHESIS: the labeling "surf1" vs "surf2" at each vertex is
% exchangeable - i.e. surf1-surf2 and surf2-surf1 are equally likely
% labelings of the same underlying pair of values, with no preference
% for one surface over the other.
%
% INPUTS:
%   surf1_vals  - Nx1 vector, value at each vertex for surface 1
%   surf2_vals  - Nx1 vector, value at each vertex for surface 2
%   Vertices    - Nx3 vertex coordinates (mm)
%   Faces       - Mx3 face/triangle index matrix
%   radius_mm   - cluster radius in mm; clusters are flipped as whole units
%                 to preserve spatial coherence under permutation
%   n_perm      - number of permutations (e.g. 2000; use exact 2^n_clusters
%                 enumeration only if n_clusters is small, otherwise Monte Carlo)
%   alpha       - significance level
%   vertex_idx  - (optional) ROI vertex indices; if omitted, whole mesh is tested
%
% OUTPUT: struct with
%   prop_pos_obs       - observed proportion of vertices favoring surf1
%   dominance_index    - observed (n_pos - n_neg)/n_total
%   perm_prop_pos      - permutation distribution of prop_pos (n_perm x 1)
%   p_value            - two-sided p-value: P(|perm stat - 0.5| >= |obs stat - 0.5|)
%   ci_null            - [2.5, 97.5] percentile range of the null distribution
%   dominates          - true if observed falls outside ci_null

    if nargin < 6, n_perm = 2000; end
    if nargin < 7, alpha = 0.05; end
    if nargin < 8, vertex_idx = []; end

    surf1_vals = surf1_vals(:);
    surf2_vals = surf2_vals(:);
    N = numel(surf1_vals);
    diff_vals = surf1_vals - surf2_vals;

    cluster_id = grow_clusters_mm(VertConn, Vertices, N, radius_mm);

    if isempty(vertex_idx)
        test_idx = (1:N)';
    else
        test_idx = vertex_idx(:);
    end

    % Identify clusters touching the ROI, and restrict their membership to ROI vertices
    clusters_in_roi = unique(cluster_id(test_idx));
    n_clusters = numel(clusters_in_roi);
    cluster_members = cell(n_clusters, 1);
    for c = 1:n_clusters
        cid = clusters_in_roi(c);
        cluster_members{c} = intersect(find(cluster_id == cid), test_idx);
    end

    % ---- Observed statistic ----
    obs_vals = diff_vals(test_idx);
    obs_vals = obs_vals(~isnan(obs_vals) & obs_vals ~= 0);
    n_pos_obs = sum(obs_vals > 0);
    n_tot_obs = numel(obs_vals);
    prop_pos_obs = n_pos_obs / n_tot_obs;
    dom_idx_obs = (2*n_pos_obs - n_tot_obs) / n_tot_obs;
    obs_stat = abs(prop_pos_obs - 0.5);

    % ---- Permutation: randomly flip sign of each CLUSTER's contribution ----
    % Flipping a cluster = swapping which surface is "1" vs "2" for all its
    % vertices at once, preserving spatial coherence under the null.
    perm_prop_pos = nan(n_perm, 1);

    % Precompute, per cluster, the values and counts (fixed across permutations)
    cluster_vals = cell(n_clusters, 1);
    for c = 1:n_clusters
        v = diff_vals(cluster_members{c});
        v = v(~isnan(v) & v ~= 0);
        cluster_vals{c} = v;
    end

    for p = 1:n_perm
        flip = (rand(n_clusters,1) > 0.5) * (-2) + 1; % +1 or -1 per cluster
        n_pos_p = 0;
        n_tot_p = 0;
        for c = 1:n_clusters
            v = cluster_vals{c} * flip(c);
            n_pos_p = n_pos_p + sum(v > 0);
            n_tot_p = n_tot_p + numel(v);
        end
        if n_tot_p > 0
            perm_prop_pos(p) = n_pos_p / n_tot_p;
        end
    end

    perm_prop_pos = perm_prop_pos(~isnan(perm_prop_pos));

    perm_stat = abs(perm_prop_pos - 0.5);
    p_value = mean(perm_stat >= obs_stat);

    ci_null = prctile(perm_prop_pos, [100*alpha/2, 100*(1-alpha/2)]);

    result.prop_pos_obs = prop_pos_obs;
    result.dominance_index = dom_idx_obs;
    result.perm_prop_pos = perm_prop_pos;
    result.p_value = p_value;
    result.ci_null = ci_null;
    result.dominates = (prop_pos_obs < ci_null(1)) || (prop_pos_obs > ci_null(2));
    result.n_clusters_used = n_clusters;
end


function adj = build_adjacency(Faces, N)
    i = [Faces(:,1); Faces(:,2); Faces(:,3); Faces(:,2); Faces(:,3); Faces(:,1)];
    j = [Faces(:,2); Faces(:,3); Faces(:,1); Faces(:,1); Faces(:,2); Faces(:,3)];
    adj = sparse(i, j, true, N, N);
    adj = adj | adj';
end


function cluster_id = grow_clusters_mm(adj, Vertices, N, radius_mm)
    cluster_id = zeros(N, 1);
    unvisited = true(N, 1);
    next_id = 0;
    seed_order = randperm(N);

    for s = 1:N
        seed = seed_order(s);
        if ~unvisited(seed)
            continue;
        end
        next_id = next_id + 1;

        frontier = seed;
        visited_this_cluster = false(N,1);
        visited_this_cluster(seed) = true;
        seed_coord = Vertices(seed, :);

        for hop = 1:30
            neighbors = find(any(adj(frontier, :), 1));
            neighbors = neighbors(unvisited(neighbors) & ~visited_this_cluster(neighbors));
            if isempty(neighbors)
                break;
            end
            d = sqrt(sum((Vertices(neighbors,:) - seed_coord).^2, 2));
            neighbors = neighbors(d <= radius_mm);
            if isempty(neighbors)
                break;
            end
            visited_this_cluster(neighbors) = true;
            frontier = neighbors;
        end

        members = find(visited_this_cluster);
        cluster_id(members) = next_id;
        unvisited(members) = false;
    end
end