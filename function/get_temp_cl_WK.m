function [temp_cl]=get_temp_cl_WK(a,s,dm0_order,k, initial_n,spectral_method)
% obtain groups of filtered centers with given k (min strategy)
%Input:
%   a: filtered centers
%   s: size of each subset
%   dm0_order: ordering matrix corresponding to dm0
%   k: the given number of clusters
%   spectral_method: 1 for unnormalized cut, 2 for Normalized Cut 3 for Symmetric Normalized Cut. Default is 3

%Output:
%temp_cl: cell array containing the grouping results of filtered centers.

if nargin<6
    spectral_method=3;
end
if nargin<5
    initial_n=1;
end


lg_a=length(a);
if (lg_a>k)
    sdm=sim_mat(lg_a,dm0_order,s,a);
    groups=spec_clus_withsim(sdm,k,spectral_method);
    groups=groups{1};

    temp_cl = cell(1, k);
    for i = 1:k
        temp_cl{i} = find(groups == i);
    end

    temp_cl=group_adjust(temp_cl,sdm,spectral_method,initial_n);

    for i = 1:k
        temp_cl{i} = a(temp_cl{i});
    end
else
    temp_cl = arrayfun(@(x) a(x), 1:length(a), 'UniformOutput', false);
end
end




