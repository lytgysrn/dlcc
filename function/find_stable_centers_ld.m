function [belong,tcl] = find_stable_centers_ld(sort_a, dm_order_sub, sim_matrix, ld_value, th_pre)
% Algorithm 1 in the paper, define stable centers and group local centers

    % Default threshold range if not provided
    if nargin < 5
        th_pre = [0.4, 0.75];
    end

    % Sort points by ld_value in descending order
    save_order = sort_a;
    [~, sort_idx] = sort(ld_value, 'descend');
    sort_a = sort_a(sort_idx);

    % Initialization
    N = length(sort_a);
    stable = true(N, 1);
    ms = zeros(N, 1);
    closest_nbr = zeros(N, 1);
    
    % Precompute valid neighbors for each point
    valid_neighbors_list = cell(N, 1);
    for x = 1:N
        neighbors = dm_order_sub(:, sort_a(x));
        valid_neighbors = neighbors(ismember(neighbors, sort_a));
        valid_neighbors_list{x}= arrayfun(@(v) find(save_order == v, 1), valid_neighbors);
    end

    % Compute mean similarity (ms) and closest neighbor for each point
    for i = 1:N
        pos_idx = valid_neighbors_list{i};
        if length(pos_idx) > 1
          %  pos_idx = arrayfun(@(v) find(save_order == v, 1), valid_neighbors);
            sim_to_neighbors = sim_matrix(pos_idx(1), pos_idx(2:end));
            ms(i) = mean(sim_to_neighbors);
            [~, max_idx] = max(sim_to_neighbors);
            closest_nbr(i) = find(sort_a == save_order(pos_idx(max_idx + 1)), 1);
        else
        closest_nbr(i)=i;
        end
    end

    % Compute threshold (th_list) for each point
    th_list = arrayfun(@(x) min(ms(x), ms(closest_nbr(x))), 1:N);
    th_list = max(th_list, th_pre(1));
    th_list = min(th_list, th_pre(2));

    % Determine stable points
    for i = 1:N
        valid_neighbors = valid_neighbors_list{i};
        if length(valid_neighbors) > 1
          %  pos_idx = arrayfun(@(v) find(save_order == v, 1), valid_neighbors);
            sim_to_neighbors = sim_matrix(valid_neighbors(1), valid_neighbors(2:end));
            valid_neighbors = valid_neighbors(2:end);  % Exclude the first point itself
            valid_neighbors = valid_neighbors(sim_to_neighbors >= th_list(i));   
        end
        if i > 1
            higher_rank_neighbors = valid_neighbors(ismember(save_order(valid_neighbors), sort_a(1:i-1)));
        else
            higher_rank_neighbors = [];
        end
        if ~isempty(higher_rank_neighbors)
            stable(i) = false;
        end
    end

    % Extract stable centers
    stable_center = sort(sort_a(stable));

    % Assign each point to its closest stable center
    stable_centers = arrayfun(@(v) find(save_order == v, 1), stable_center);
 %   belong = arrayfun(@(x) stable_centers(find(sim_matrix(x, pos_simmat) == max(sim_matrix(x, pos_simmat)), 1)), 1:N);
    [belong,tcl]=simplesplit(sim_matrix,1:N,stable_centers);
    belong=struct('belong',belong,'stable_centers',stable_centers);

end
