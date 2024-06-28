function a = filter_center(dm0, dm0_order, s, ranktopobs, rank_mat, method)
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

if strcmp(method,'min')
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
    %remove local centers contribute 0 to cumulative proportion
    props = arrayfun(@(x) f_prop(sort_a2(1:x), s, Nobs, dm0_order), 1:length(sort_a2));
    d_props = diff(props);
    wd0 = find(d_props == 0);
    if ~isempty(wd0)
        sort_a2 = sort_a2(setdiff(1:length(sort_a2), wd0+1));
    end
    la=length(sort_a2);
    p2 = arrayfun(@(x) f_prop(sort_a2(1:x), s, Nobs, dm0_order), 1:la);

    % if the final cumulative proportion is less than 3/4, select all
    % local centers
    if max(p2) < 0.75
        a = sort_a2;
    else
        % give a symmetric similarity matrix
        dmsub = dm0(sort_a2, sort_a2);
        dmsub = (dmsub+dmsub.')/2;
        dmsub(1:la+1:end)=0;
        dmsub=tril(dmsub);

        mss = arrayfun(@(x) max(dmsub(x, :)), 1:la);

        dp2 = diff(p2);%cumulative proportion gap
        sort_p2 = sort(dp2, 'descend');%orderthe gap
        %if the number of remaining local centers is too small
        if length(sort_p2) <= 1
            a = sort_a2;
        else
            prop_pr = -diff(sort_p2);
            prop_pr = round(prop_pr, 8, 'significant');
            [~, p_id] = sort(prop_pr, 'descend');
            rank_pr = ones(1, la) * la;%initial proportion rank
            for i = 1:length(p_id)
                pr_id = find(dp2 == sort_p2(p_id(i))) + 1;
                %find the position
                if length(pr_id) > 1
                    labelled = find(rank_pr(pr_id) ~= la);
                    if ~isempty(labelled)
                        pr_id = setdiff(pr_id, pr_id(labelled));
                    end
                    pr_id = pr_id(1);
                end
                rank_pr(pr_id) = i;
            end

            %then define the rank based on the similarity bet new center and existed centers
            mss(1) = 1000;
            rank_dis = tiedrank(mss);
            rank_dis(1) = 1000;
            frank = rank_pr * 0.6 + 0.4 * rank_dis;
            tmp_l = find(p2 >= 0.75, 1);
            frank = frank(tmp_l:end);
            lab_pos = tmp_l - 1 + find(frank == min(frank), 1);
            a = sort_a2(1:lab_pos);
        end
    end
elseif strcmp(method,'max')
    a=sort_a2(sort_count>1);
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
end

a=sort(a);
end


function prop = f_prop(a, s, Nobs, dm0_order)
all_nbs = unique(dm0_order(1:s,a));
prop = length(all_nbs) / Nobs;
end

