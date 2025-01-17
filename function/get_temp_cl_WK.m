function [temp_cl]=get_temp_cl_WK(a_info,s,dm0_order,th_pre,K)
% obtain groups of filtered centers (min strategy)
%Input:
%   a_info: filtered centers along with local depth value
%   s: size of each subset
%   dm0_order: ordering matrix corresponding to dm0
%   K: the given number of clusters (optional)
%   th_pre: pre_defined th range when grouping local centers
%Output:
%temp_cl: cell array containing the grouping results of filtered centers.

if nargin<5
    K=[];
end
if nargin<4
    th_pre=[0.4,0.75];
end

a=a_info.a;
ld_value=a_info.ld_value;

lg_a=length(a);

sdm=sim_mat(lg_a,dm0_order,s,a);
[belong,temp_cl]=find_stable_centers_ld(a, dm0_order(1:floor(s/2), :), sdm, ld_value, th_pre);

k = length(temp_cl);
num_points=length(belong.belong);

if ~isempty(K)
    if k < K
        stable_centers=belong.stable_centers;
        while k < K
            non_stable_points = setdiff(1:num_points, stable_centers); % non-stable centers
            sub_matrix = sdm(non_stable_points, stable_centers);
            max_sims=max(sub_matrix,[],2);

            %find max sim to other centers with higher ld value
            ld_sub=ld_value(non_stable_points);
            ld_matrix = ld_value' >ld_sub; % (n x m)
            valid_sim_matrix = sdm(non_stable_points, :) .* ld_matrix'; % (m x n)
            max_sim = max(valid_sim_matrix, [], 2); 
           % define candidates
            sub_matrix2 = sdm(non_stable_points, non_stable_points);
            sub_matrix2(1:size(sub_matrix2, 1) + 1:end) = 0;
            neighbor_mask = (sub_matrix2 > max_sim);
            ld_matrix = repmat(ld_sub, size(sub_matrix2, 1), 1); 
            neighbor_ld_values = neighbor_mask .* ld_matrix; 
            nonzero_mask = neighbor_ld_values > 0; 
            ld_comparison = (neighbor_ld_values < ld_sub') & nonzero_mask; 
            less_ld_points=sum(ld_comparison, 2); 
            %ld_value_ratio = less_ld_points./ sum(nonzero_mask, 2); 
            valid_candidates=non_stable_points(less_ld_points>0);
           %compute scores, and choose the minimun one
            [~, best_idx] = min(max_sims(less_ld_points>0).*exp(-0.05 * less_ld_points(less_ld_points>0)));
            top_index = valid_candidates(best_idx); 
            %update stable center
            stable_centers= [stable_centers,top_index];
            k=k+1;
        end
        [~,temp_cl] = simplesplit(sdm,1:num_points,stable_centers);
    elseif k > K
        first_appear=cellfun(@min,temp_cl);
        [~, drop_idx] = maxk(first_appear, k-K);
        temp_cl(drop_idx) = [];
    end
end


temp_cl=group_adjust(temp_cl,sdm);
k=length(temp_cl);
for i = 1:k
    temp_cl{i} = a(temp_cl{i});
end

end




