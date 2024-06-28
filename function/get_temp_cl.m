function temp_cl = get_temp_cl(a,dm0,s,Th,method,dm0_order,mask)
%grouping centers after filtering

%Input:
%   a: filtered centers
%   dm0: depth-based similarity matrix
%   s: size of each subset
%   Th: threshold for grouping filtered centers
%   method: strategy for grouping (min or max)
%   dm0_order: ordering matrix corresponding to dm0
%   mask: logical parameter. If true, the output labels are based on the whole data set; otherwise, the labels are based on the filtered centers only.

%Output:

%temp_cl: cell array containing the grouping results of filtered centers.


center_length = length(a);

if nargin<6
    [~, dm0_order] = sort(dm0, 2, 'descend');
    dm0_order=dm0_order';
end

if nargin<7
    mask = true;
end

if strcmp(method,'min')
    temp_save = cell(1,center_length);
    sim_mat1 = sim_mat(center_length,dm0_order,s,a);
    %for each center finding all other centers with sim>Th
    for i = 1:center_length
        temp_save{i} = find(sim_mat1(i,:)>Th);
    end

    temp_cl = cell(1,center_length);
    for j = 1:center_length
        if length(temp_save{j})>1
            temp_cl{j} = multi_intersect(arrayfun(@(x) temp_save{x}, temp_save{j}, 'UniformOutput', false));
        else
            temp_cl{j}=temp_save{j};
        end
    end
    %drop same groups
    Tdup =find_duplicate_cells(temp_cl);
    if ~isempty(Tdup)
        temp_cl(Tdup) = [];
    end
    %delete group x if group x in y
    l_ts = length(temp_cl);
    rel_mat = zeros(l_ts,l_ts);
    for i = 1:l_ts
        for j = 1:l_ts
            rel = length(intersect(temp_cl{i},temp_cl{j}))/length(temp_cl{i});
            if(rel~=1)
                rel=0;
            end
            rel_mat(i,j) = rel;
        end
    end
    dr = find(sum(rel_mat,2)>1);
    if ~isempty(dr)
        temp_cl(dr) = [];
    end
    % build between group similarity matrix to see if any of them can be
    % merged
    temp_K = length(temp_cl);
    gap_bet_clus = zeros(temp_K,temp_K);
    for y = 1:(temp_K-1)
        for u = (y+1):temp_K
            c1 = temp_cl{y};
            c2 = temp_cl{u};
            min_sim = min(min(sim_mat1(c2,c1)));
            gap_bet_clus(u,y) = min_sim;
            gap_bet_clus(y,u) = min_sim;
        end
    end
    Thset = sort(unique(gap_bet_clus(:),'stable'),'descend');
    Thset = Thset(Thset>Th);
    lt = length(Thset);
    gap_bet_clus(1:temp_K+1:end) = 1;
    if lt>0
        m = 1;
        while m<=lt
            v = 1;
            while v <= temp_K
                if any(Thset(m)==gap_bet_clus(:,v))
                    comb_lab = find(gap_bet_clus(:,v)==Thset(m))';
                    if min(min(gap_bet_clus(comb_lab,comb_lab)))<=Th
                        comb_lab = comb_lab(1);
                    end
                    temp_cl{v} = unique([temp_cl{v}, temp_cl{comb_lab}]);
                    temp_cl(comb_lab) = [];
                    gap_bet_clus(1:temp_K,v) = arrayfun(@(x) min(gap_bet_clus(x,[v,comb_lab])), 1:temp_K);
                    gap_bet_clus(v,1:temp_K) = gap_bet_clus(1:temp_K,v);
                    gap_bet_clus = gap_bet_clus(setdiff(1:temp_K,comb_lab), setdiff(1:temp_K,comb_lab));
                    temp_K = temp_K - length(comb_lab);
                    gap_bet_clus(1:temp_K+1:end) = 1;
                else
                    v = v + 1;
                end
            end
            m = m + 1;
        end
    end


elseif strcmp(method,'max')
    sim_mat1 = sim_mat(center_length,dm0_order,s,a);
    list_group = arrayfun(@(x) find(sim_mat1(x, :) > Th), 1:center_length, 'UniformOutput', false);

    nclus = 1;
    temp_cl{nclus} = list_group{1};
    for i = 2:center_length
        for j = 1:length(temp_cl)
            %check similarity between center i and centers in temp_cl
            sim = sim_mat1(i,temp_cl{j});
            if max(sim)>Th
                temp_cl{j} = [temp_cl{j}, list_group{i}];
            end
        end
        current_a = unique([temp_cl{:}]);
        %if ith a does not includes, a new group is defined.
        if ~ismember(i, current_a)
            nclus = nclus + 1;
            temp_cl{nclus} = list_group{i};
        end
    end
end
%check if any centers are multi-assigned, if yes, drop them (min) or merge groups (max)
[ab, ~,ic] = unique([temp_cl{:}]);
counts = accumarray(ic, 1);
rep_label = ab(counts > 1);
if ~isempty(rep_label)

    if strcmp(method,'min')
        for i = 1:length(rep_label)
            cllabel = [];
            for j = 1:length(temp_cl)
                if any(rep_label(i)==temp_cl{j})
                    cllabel = [cllabel, j];
                end
            end
            for m = cllabel
                temp_cl{m} = setdiff(temp_cl{m},rep_label(i));
            end
        end
    elseif strcmp(method,'max')
        while ~isempty(rep_label)
            cllabel = [];
            for j = 1:length(temp_cl)
                if any(rep_label(1)==temp_cl{j})
                    cllabel = [cllabel, j];
                end
            end
            remain_cluster = cllabel(1);
            obs = unique([temp_cl{cllabel}]);
            temp_cl{remain_cluster} = obs;
            cllabel = cllabel(2:end);
            for w = flip(cllabel)
                temp_cl(w) = [];
            end
            [ab, ~,ic] = unique([temp_cl{:}]);
            counts = accumarray(ic, 1);
            rep_label = ab(counts > 1);
        end
    end
end
if mask
    for i = 1:length(temp_cl)
        temp_cl{i} = a(temp_cl{i});
    end
end
end



function idx = find_duplicate_cells(c)
str = cellfun(@mat2str, c, 'UniformOutput', false);
[~, idx] = unique(str);
idx = setdiff(1:length(c), idx);
end

function output = multi_intersect(cellArray)
output = cellArray{1};
for k = 2:numel(cellArray)
    output = intersect(output, cellArray{k});
end
end



