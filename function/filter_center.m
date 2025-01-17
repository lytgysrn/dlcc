function a = filter_center(dm0, dm0_order, s, ranktopobs, rank_mat, method,ld_value,th_pre)
% Local centers filtering procedures
%input:
%data:data set
%dm0: depth-based similarity matrix
%dm0_order: the corresponding ordering matrix of dm0
%ranktopobs: depth center of each subset
%rank_mat: the depth rank of the obs in its own subset.
%s: size of each subset
%method: min or max strategy
% output:
%a: filtered centers

if nargin<6
    method='min';
end
if strcmp(method,'max')
    ld_value=0;
end

if nargin<8
    th_pre=[0.4,0.75];
end

Nobs = size(dm0, 1);
%a_all: all local centers
[a_all, ~,ic] = unique(ranktopobs);
%count: frequencies of a_alln.sort.choose
count = accumarray(ic, 1)';
idx1=rank_mat(a_all) <= 2;
a_2 = a_all(idx1);
count=count(idx1);

[sort_count,sort_lab]=sort(count,"descend");
sort_a2=a_2(sort_lab);
sort_a2=sort_a2(sort_count>1);
sort_count=sort_count(sort_count>1);


if strcmp(method,'min')
    maxp=f_prop(sort_a2, s, Nobs, dm0_order);
    if maxp<0.75
        nbr_all_a = dm0_order(1:s, sort_a2);
        sim_matrix = sim_mat(length(sort_a2),nbr_all_a,s);
        [test1,~]=find_stable_centers_ld(sort_a2, dm0_order(1:floor(s/2), :), sim_matrix, ld_value(sort_a2), th_pre);
        test1=test1.belong;
        intra_avg_sim = arrayfun(@(i) (sum(sim_matrix(i, test1 == test1(i))) - 1) / (sum(test1 == test1(i)) - 1), ...
            1:size(sim_matrix, 1))';
        inter_max_sim = arrayfun(@(i) max(sim_matrix(i, test1 ~= test1(i))), 1:size(sim_matrix, 1))';
        a_sub=sort_a2((intra_avg_sim-inter_max_sim)>0);
        anbs = unique(dm0_order(1:s, a_sub));
        maxp = length(anbs) / Nobs;
        upper_p=f_prop(a_all, s, Nobs, dm0_order);
        while maxp<0.9*upper_p
            others=setdiff(a_all,a_sub);
            new_coverages = arrayfun(@(x) length(setdiff(dm0_order(1:s, x), anbs)), others);

            max_gain = max(new_coverages);
            candidates = others(new_coverages == max_gain); 

            if length(candidates) > 1
                sim_sums = arrayfun(@(x) sum(sim_matrix(x, :)), candidates); 
                [~, min_idx] = min(sim_sums);
                best_point = candidates(min_idx);
            else
                best_point = candidates(1);
            end
            a_sub = [a_sub, best_point];
            anbs = unique([anbs; dm0_order(1:s, best_point)]);
            maxp = length(anbs) / Nobs; 
        end
        ld_value=ld_value(a_sub);
        a=struct('a', a_sub, 'ld_value',ld_value);
    else
    %if two centers with same frequency, using contributed proportion to order them
    dupli_lab = find(diff(sort_count) == 0) + 1;
    if ~isempty(dupli_lab)
        for i = 1:length(dupli_lab)
            prop1 = f_prop(sort_a2(1:dupli_lab(i)-1), s, Nobs, dm0_order);
            prop2 = f_prop(sort_a2(setdiff(1:dupli_lab(i), dupli_lab(i)-1)), s, Nobs, dm0_order);
            if prop1 < prop2
                tmp = sort_a2(dupli_lab(i));
                sort_a2(dupli_lab(i)) = sort_a2(dupli_lab(i)-1);
                sort_a2(dupli_lab(i)-1) = tmp;
            end
        end
    end
    total_l = length(sort_a2);

    nbr_all_a = dm0_order(1:s, sort_a2);
    sim_matrix = sim_mat(total_l,nbr_all_a,s);
    ld_value=ld_value(sort_a2);
    [test1,tcl]=find_stable_centers_ld(sort_a2, dm0_order(1:floor(s/2), :), sim_matrix, ld_value, th_pre);
    test1=test1.belong;
    
    neighbors_list = cellfun(@(x) unique(nbr_all_a(:, x(:))), tcl, 'UniformOutput', false);
    ctn = true;
    while ctn
        % check if group should be deleted
        unique_neighbors_num = arrayfun(@(i) ...
            length(setdiff(neighbors_list{i}, unique(cell2mat(neighbors_list([1:i-1, i+1:end]))))), ...
            1:length(neighbors_list));


        to_remove = find(unique_neighbors_num < s / 10);
        [~, min_idx] = min(unique_neighbors_num(to_remove));
        to_remove = to_remove(min_idx);

        ctn = ~isempty(to_remove);
        if ctn
            save_sort = sort_a2;
            [remove_point,tcl] = max_sim_remove(tcl, sim_matrix, to_remove);

            % update
            sort_a2(remove_point) = [];
            tcl = cellfun(@(x) cell2mat(arrayfun(@(v) find(sort_a2 == save_sort(v), 1), x, 'UniformOutput', false)), tcl, 'UniformOutput', false);
            nbr_all_a = dm0_order(1:s, sort_a2);
            neighbors_list = cellfun(@(x) unique(nbr_all_a(:, x(:))), tcl, 'UniformOutput', false);
            test1(remove_point) = [];
            ld_value(remove_point) = [];
            sim_matrix(remove_point, :) = [];
            sim_matrix(:, remove_point) = [];
        end
    end
    total_l = length(sort_a2);
    props = arrayfun(@(x) f_prop(sort_a2(1:x), s, Nobs, dm0_order), 1:total_l);
    %initial start
    if (props(end)>0.75)
        [~, first_occurrence] = unique(test1, 'stable');
        begin_idx = max(find(props >= 0.75, 1, 'first')+1, first_occurrence(end) + 1);
    else
        begin_idx=total_l;
    end
    if begin_idx < total_l
        idx = begin_idx;
        notes = zeros(1, total_l);
        notes2 = zeros(1, total_l);
        stop = false;
        while ~stop
            current_group = test1(idx);
            relevant_indices = find(test1(1:idx) ~= current_group);
            max_bg_sim = round(max(sim_matrix(relevant_indices, idx)), 5);

            % Find same group points
            same_group_idx = setdiff(1:(idx-1), relevant_indices);
            %min max between group sim of other points in same group
            sim_compare_min = min(max(sim_matrix(relevant_indices, same_group_idx), [], 1));
            %same group similarity matrix
            sub_sim = sim_matrix([same_group_idx, idx], [same_group_idx, idx]);
            max_ingroup = sum(sub_sim, 1);

            notes(idx) = max_ingroup(end) == min(max_ingroup); % Class internal dissimilarity
            notes2(idx) = max_bg_sim <= sim_compare_min;

            if idx == total_l
                stop = true;
            end
            idx = idx + 1;
        end
        stop_candidates = intersect(find(notes == 1), find(notes2 == 1));

        %detect bad points in each range
        l_sc = length(stop_candidates);
        if l_sc > 0
            stop_candidates = [begin_idx-1, stop_candidates];
            bad_point = zeros(l_sc + 1, 1);
            for j = 1:length(stop_candidates)
                tcl_b = cellfun(@(x) x(x <= stop_candidates(j)), tcl, 'UniformOutput', false);
                check = score_computer(tcl_b, sim_matrix);
                bad_point(j) = sum(cellfun(@(x) sum(x < 0), check.difwb2)) + ...
                    sum(cellfun(@(x) sum(x < 0), check.difwb));
            end
            gap = diff(bad_point);
            if max(gap) > 1
                stop_candidate = stop_candidates(find(gap > 1, 1));
            else
                stop_candidate = max(stop_candidates);
            end
        else
            stop_candidate =begin_idx-1;
        end
        a = sort_a2(1:stop_candidate);
        ld_value=ld_value(1:stop_candidate);
    else
        a=sort_a2;
    end
    a=struct('a', a, 'ld_value',ld_value);
    end

elseif strcmp(method,'max')
    a=sort_a2;
    %delete isolated center
    sdm=sim_mat(length(a),dm0_order,s,a);
    No_iso = sum(sdm)~=1;
    a=a(No_iso);
    all_nbs = unique(dm0_order(1:s, a));
    prop = length(all_nbs) / Nobs;
    %if the proportion is not large enough, some regions may not be
    %covered, (happened when clusters are non-convex locally).
    if prop < 0.85
        a_2 = setdiff(a_all, a);
        nbs_a_2 = cell(1, length(a_2));
        for i = 1:length(a_2)
            nbs_a_2{i} = dm0_order(1:s, a_2(i));
        end
        l_a_2 = zeros(1, length(a_2));
        for g = 1:length(a_2)
            l_a_2(g) = sum(ismember(nbs_a_2{g}, all_nbs));
        end
        % a_3 contains centers that have overlapping neighbors with the centers in a.
        a_3 = a_2(l_a_2 ~= 0);
        % a_2 represents centers that are isolated from centers in 'a'.
        a_2= a_2(l_a_2==0);
        if ~isempty(a_2)
            sim_group_ne=sim_mat(length(a_3),dm0_order,s,a_3,a_2);
            sim_group_e=sim_mat(length(a_3),dm0_order,s,a_3,a);
            mne=max(sim_group_ne,[],2);
            %mne: max similarity between centers in a_3 and centers in a_2
            me=max(sim_group_e,[],2);
            %me: max similarity between centers in a_3 and a
            idx=mne-me>0;
            % find centers in a_3 which are more close to centers in a_2
            [~, in_groupe_idx] = find(sim_group_e(idx,:) > 0);
            in_groupe_idx=unique(in_groupe_idx);
            delete_idx=a(in_groupe_idx);
            %drop those centers in a to make sure non-convex shapes clusters
            %and convex shapes clusters will not influnece each other.
            a=setdiff(a,delete_idx);
            a_3=a_3(idx);
            a_2=union(a_2,a_3);
            a = union(a,a_2);
            sdm=sim_mat(length(a),dm0_order,s,a);
            No_iso = sum(sdm)~=1;
            a=a(No_iso);
            a_disp = [num2str(a_2(1:3)) '...'];
            warning(['Cluster(s) based on centers ' a_disp ' may be non-convex locally. The clustering result of this cluster may need further verification.'])
        end
    end
    a=sort(a);
end

end


function prop = f_prop(a, s, Nobs, dm0_order)
all_nbs = unique(dm0_order(1:s,a));
prop = length(all_nbs) / Nobs;
end

