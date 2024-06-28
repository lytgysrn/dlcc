function current_label = cluster2cv(data, temp_clus)
% Convert clustering label (the labels of observations in the whole data
% set) to class vector,where each observation is assigned a label based on its cluster membership.

    current_label = zeros(size(data, 1), 1);
    for i = 1:length(temp_clus)
        current_label(temp_clus{i}) = i;
    end
end