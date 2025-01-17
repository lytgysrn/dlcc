function temp_cl = get_temp_cl(a,dm0,s,Th,dm0_order,mask)
%grouping centers after filtering under max strategy

%Input:
%   a: filtered centers
%   dm0: depth-based similarity matrix
%   s: size of each subset
%   Th: threshold for grouping filtered centers
%   dm0_order: ordering matrix corresponding to dm0
%   mask: logical parameter. If true, the output labels are based on the whole data set; otherwise, the labels are based on the filtered centers only.

%Output:

%temp_cl: cell array containing the grouping results of filtered centers.


center_length = length(a);

if nargin<5
    [~, dm0_order] = sort(dm0, 2, 'descend');
    dm0_order=dm0_order';
end

if nargin<6
    mask = true;
end

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


% if strcmp(method,'max')
[ab, ~,ic] = unique([temp_cl{:}]);
counts = accumarray(ic, 1);
rep_label = ab(counts > 1);
if ~isempty(rep_label)
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

    if mask
        for i = 1:length(temp_cl)
            temp_cl{i} = a(temp_cl{i});
        end
    end
end




