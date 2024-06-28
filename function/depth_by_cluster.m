function DD = depth_by_cluster(X, Kclus, cluster, sub, depth)
% Compute data depth with respect to each cluster
% Inputs:
%   X: data set
%   Kclus: number of clusters
%   cluster: cell array containing observation labels in each cluster
%   sub: subset index, indicating whose data depth is calculated
%   depth: type of data depth

    subX = X(sub, :);
    DD = zeros(size(subX, 1), Kclus);
    
    for i = 1:Kclus
        if strcmp(depth, 'mahalanobis')
            DD(:, i) = Maha_d(subX, X(cluster{i}, :), 'MCD');
        elseif strcmp(depth, 'spatial')
            DD(:, i) = spatial_d(subX, X(cluster{i}, :));
        end
    end
end