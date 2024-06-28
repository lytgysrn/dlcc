function temp_clus = get_temp_cluster(a,dm0, temp_clus, Kclus, temp_cl, s, method,count_list)
%Obtain temporary clusters
 %Input: 
%   a: filtered centers
%   dm0: depth-based similarity matrix
%   temp_clus: Initial pools of clustering, containing all neighbors of centers in each group.
%   Kclus: the number of groups in temp_clus.
%   temp_cl: grouping results of filtered centers.
%   s: size of each subset
%   method: strategy for grouping (min or max)
%   count_list: obs frequencies in each group (only for the max strategy)

%Output:
  %temp_clus: final temporary clustering results
  
if nargin<8
count_list=[];
end

clus_ind = 1:Kclus;
cl = cellfun(@length,temp_cl);
cl = repelem(clus_ind,cl);
%remove overlapped obs from each temporary cluster
[temp_clus, ulabel]=remove_ol_obs(clus_ind,temp_clus);

if strcmp(method, 'min')
    %compute scores for obs in temporary clusters
    ascore = Assign_score(a,Kclus, dm0, temp_clus, temp_cl);
    tobs=cat(1,temp_clus{:});
    tobs=tobs';
    otherobs = setdiff(1:size(dm0,1),tobs);

    %compute positive scores for unlabelled obs
    [~, idx] = max(dm0(a, otherobs));
    group_belong = cl(idx);
    left_clus = cell(1,Kclus);
    for i = clus_ind
        left_clus{i} = otherobs(group_belong == i);
    end
    score_left = Assign_score(a,Kclus, dm0, left_clus, temp_cl);

    %median of score_left is treated as the threshold for the obs in temporary clusters
    lq = arrayfun(@(x) quantile(score_left{x},0.5), 1:Kclus);
    %move those obs with scores<lq to the left obs pool.
    l0u = [];
    for i = 1:Kclus
        index = find(ascore{i} < lq(i));
        if ~isempty(index)
            left_clus{i} = [left_clus{i}, temp_clus{i}(index).'];
            score_left{i} = [score_left{i}, ascore{i}(index)];
            temp_clus{i}(index) = [];
            ascore{i}(index) = [];
            TF = find(score_left{i} < 0);
            if ~isempty(TF)
                l0u = [l0u, left_clus{i}(TF)];
                score_left{i}(TF) = [];
                left_clus{i}(TF) = [];
            end
        end
    end
    % if l0u is not empty, relocate obs with negative scores in the left pool
    if ~isempty(l0u)
        [~, idx] = max(dm0(a, l0u));
        st2 = cl(idx);
        left_clus2 = cell(1, Kclus);
        for i = clus_ind
            left_clus2{i} = l0u(st2 == i);
        end
        score_left2 = Assign_score(a,Kclus, dm0, left_clus2, temp_cl);
        left_clus = cellfun(@(x,y) [x, y], left_clus, left_clus2, 'UniformOutput', false);
        score_left = cellfun(@(x,y) [x, y], score_left, score_left2, 'UniformOutput', false);
    end

    %check after above processes, are there any empty temporary clusters.
    lg_tmc= cellfun(@length, temp_clus);
    idx=find(lg_tmc==0, 1);
    if ~isempty(idx)
        error('0 observation in some cluster');
    end

    %set the threshold to accept some obs in the left pool
    new_border = zeros(1, Kclus);
    for i = clus_ind
        num_obs_current = length(temp_clus{i});
        sortscore = sort(ascore{i}, 'descend');
        sortscore = sortscore(floor(num_obs_current/2):num_obs_current);
        [~,idx]=min(diff(sortscore));
        cdd_1 = sortscore(idx);
        if num_obs_current < (s/2)
            temp_s = 1 - (s/2 - num_obs_current) / length(score_left{i});
            if temp_s > 0
                cdd_2 = quantile(score_left{i}, temp_s);
            else
                cdd_2 = 0;
            end
        else
            cdd_2 = 1;
        end
        new_border(i) = min(cdd_1, cdd_2);
    end
    for x = clus_ind
        index = score_left{x} > new_border(x);
        left_clus{x} = left_clus{x}(index);
    end

    temp_clus = cellfun(@transpose, temp_clus, 'UniformOutput', false);
    temp_clus = cellfun(@(x, y) [x, y], temp_clus, left_clus, 'UniformOutput', false);
elseif strcmp(method, 'max')
if ~isempty(ulabel)
    % Remove those obs in temp cluster

    % if there is a cluster without any obs
    lg_tmc= cellfun(@length, temp_clus);
    idx=find(lg_tmc==0, 1);
    if ~isempty(idx)
        error('0 observation in some cluster');
    end
% score matrix building
    score_mat = zeros(length(ulabel), length(temp_cl));
        for j = 1:Kclus
           [idx,idx2]=ismember(ulabel,temp_clus{j});
           idx2=idx2(idx);
           score_mat(idx, j) = count_list{j}(idx2);            
        end
    
%assign overlapped obs to the cluster with highest score
    clus_assign = zeros(1, size(score_mat, 1));
    for i = 1:size(score_mat, 1)
        h = find(score_mat(i, :) == max(score_mat(i, :)));
        if length(h) == 1
            clus_assign(i) = h;
        else
            % Draw case, not assign
            clus_assign(i) = 0;
        end
    end

    temp_clb = cell(1, length(clus_ind));
    for i = 1:length(clus_ind)
        temp_clb{i} = ulabel(clus_assign == clus_ind(i));
    end

    for i = 1:length(clus_ind)
        temp_clus{i} = [temp_clus{i}; temp_clb{i}'];
    end
end
    temp_clus = cellfun(@transpose, temp_clus, 'UniformOutput', false);
end

end

%delete overlapped obs
function [temp_clus, ulabel]=remove_ol_obs(clus_ind,temp_clus)
t_obs=cat(1,temp_clus{:});
[tmp_t, ~,ic] = unique(t_obs);
t_obs = accumarray(ic, 1);
ulabel=tmp_t(t_obs>1);
for i = clus_ind
    temp_clus{i}=setdiff(temp_clus{i},ulabel);
end
end

%assign scores for obs under the min strategy
function deflist = Assign_score(a,Kclus, dm0, temp_clus, temp_cl)
deflist = cell(1, Kclus);
for k = 1:Kclus
    tmp_l = length(temp_clus{k});
    if tmp_l ~= 0
        max_a = max(dm0(temp_cl{k}, temp_clus{k}), [], 1);
        max_b = max(dm0(setdiff(a, temp_cl{k}), temp_clus{k}), [], 1);
        larger = max(max_a, max_b);
        deflist{k} = (max_a - max_b) ./ larger;
    end
end
end

