function [groups,tcl] = simplesplit(sdm,temp_cl_sub,stable_centers)
%simple grouping when stable centers are given

sdm2=sdm(temp_cl_sub,temp_cl_sub);
num_points = size(sdm2, 1); 
groups = zeros(1, num_points);

num_stable_centers = numel(stable_centers);
for k = 1:num_stable_centers
    groups(stable_centers(k)) = k;
end


groups=arrayfun(@(x) find(sdm2(x, stable_centers) == max(sdm2(x, stable_centers)), 1), 1:num_points);

tcl = cell(num_stable_centers, 1);
for i = 1:num_stable_centers
    tcl{i} = find(groups == i);
end
