function result = mesh_block_bootstrap_dominance(diff_vals, Faces, cluster_radius, n_boot, alpha, vertex_idx)
% MESH_BLOCK_BOOTSTRAP_DOMINANCE  Directional dominance test for values defined
% on a cortical mesh, using spatially contiguous vertex clusters as the
% resampling unit (mesh analogue of a block bootstrap).
%
% INPUTS:
%   diff_vals      - Nx1 vector, surface1-surface2 difference per vertex
%                     (N = number of vertices in the FULL mesh, e.g. 15005)
%   Faces          - Mx3 face/triangle index matrix from Brainstorm surface
%   cluster_radius - approx number of graph-hops defining cluster size
%                     (analogous to block_size; tune to decorrelation length)
%   n_boot         - number of bootstrap resamples
%   alpha          - significance level
%   vertex_idx     - (optional) indices into the FULL mesh defining an ROI.
%                     If provided, the test is restricted to this ROI but
%                     clustering still respects full-mesh adjacency at the
%                     ROI's boundary (more accurate than reclustering in isolation).
%                     If omitted, the whole mesh is tested.
%
% OUTPUT: struct with prop_pos, dominance_index, boot_prop_pos, ci, p_value, dominates

    if nargin < 5, alpha = 0.05; end
    if nargin < 6, vertex_idx = []; end

    N = numel(diff_vals);

    % ---- Build adjacency list from Faces (once; cache if calling repeatedly) ----
    adj = build_adjacency(Faces, N);

    % ---- Partition the (full) mesh into spatially contiguous clusters ----
    cluster_id = grow_clusters(adj, N, cluster_radius);
    n_clusters = max(cluster_id);

    % ---- Restrict to ROI if specified ----
    if isempty(vertex_idx)
        test_idx = (1:N)';
    else
        test_idx = vertex_idx(:);
    end

    % Map: which clusters are touched by the test region, and which of
    % their member vertices fall inside it
    clusters_in_roi = unique(cluster_id(test_idx));

    % Precompute vertex lists per relevant cluster, restricted to ROI membership
    cluster_members = cell(numel(clusters_in_roi), 1);
    for c = 1:numel(clusters_in_roi)
        cid = clusters_in_roi(c);
        members_full = find(cluster_id == cid);
        cluster_members{c} = intersect(members_full, test_idx);
    end

    % ---- Observed statistic ----
    obs_vals = diff_vals(test_idx);
    obs_vals = obs_vals(~isnan(obs_vals) & obs_vals ~= 0);
    n_pos_obs = sum(obs_vals > 0);
    n_tot_obs = numel(obs_vals);
    prop_pos_obs = n_pos_obs / n_tot_obs;
    dom_idx_obs = (2*n_pos_obs - n_tot_obs) / n_tot_obs;

    % ---- Bootstrap ----
    n_avail_clusters = numel(cluster_members);
    boot_prop_pos = nan(n_boot, 1);

    for b = 1:n_boot
        % Resample clusters with replacement until we have >= n_tot_obs vertices
        sampled_vals = [];
        % Oversample cluster draws; cheaper than checking length every iteration
        draw_count = max(n_avail_clusters, ceil(1.5 * n_avail_clusters));
        draws = randi(n_avail_clusters, draw_count, 1);

        k = 1;
        total_needed = numel(test_idx); % match original ROI size
        while numel(sampled_vals) < total_needed && k <= numel(draws)
            mem = cluster_members{draws(k)};
            sampled_vals = [sampled_vals; diff_vals(mem)]; %#ok<AGROW>
            k = k + 1;
            if k > numel(draws)
                draws = randi(n_avail_clusters, draw_count, 1); % refill if needed
                k = 1;
            end
        end

        sampled_vals = sampled_vals(~isnan(sampled_vals) & sampled_vals ~= 0);
        if isempty(sampled_vals)
            continue;
        end
        boot_prop_pos(b) = sum(sampled_vals > 0) / numel(sampled_vals);
    end

    boot_prop_pos = boot_prop_pos(~isnan(boot_prop_pos));

    ci_lower = prctile(boot_prop_pos, 100*alpha/2);
    ci_upper = prctile(boot_prop_pos, 100*(1-alpha/2));

    boot_centered = boot_prop_pos - mean(boot_prop_pos) + 0.5;
    p_value = 2 * min( mean(boot_centered >= prop_pos_obs), ...
                        mean(boot_centered <= prop_pos_obs) );
    p_value = min(p_value, 1);

    result.prop_pos = prop_pos_obs;
    result.dominance_index = dom_idx_obs;
    result.boot_prop_pos = boot_prop_pos;
    result.ci = [ci_lower, ci_upper];
    result.p_value = p_value;
    result.dominates = ~(ci_lower <= 0.5 && 0.5 <= ci_upper);
    result.n_clusters_used = n_avail_clusters;
end


function adj = build_adjacency(Faces, N)
% Build sparse adjacency matrix from triangle faces
    i = [Faces(:,1); Faces(:,2); Faces(:,3); Faces(:,2); Faces(:,3); Faces(:,1)];
    j = [Faces(:,2); Faces(:,3); Faces(:,1); Faces(:,1); Faces(:,2); Faces(:,3)];
    adj = sparse(i, j, true, N, N);
    adj = adj | adj'; % ensure symmetric
end


function cluster_id = grow_clusters(adj, N, radius)
% Partition mesh vertices into contiguous clusters via BFS region growing
% from random unvisited seed vertices, each grown out to `radius` hops.
    cluster_id = zeros(N, 1);
    unvisited = true(N, 1);
    next_id = 0;

    % Randomize seed order so cluster shapes/sizes aren't biased by vertex ordering
    seed_order = randperm(N);

    for s = 1:N
        seed = seed_order(s);
        if ~unvisited(seed)
            continue;
        end
        next_id = next_id + 1;

        % BFS out to `radius` hops from seed, staying within unvisited vertices
        frontier = seed;
        visited_this_cluster = false(N,1);
        visited_this_cluster(seed) = true;

        for hop = 1:radius
            neighbors = find(any(adj(frontier, :), 1));
            neighbors = neighbors(unvisited(neighbors) & ~visited_this_cluster(neighbors));
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