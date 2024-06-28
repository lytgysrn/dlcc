function d = if_corr(data,shape)
%a simple test for determining if there are strong correlations among
%variables, indicating which type of covariance matrix will be used in
%building the Mahalanobis-depth based similarity matrix

if nargin<2
    data = zscore(data); % Standardize the data
    [~, ~, latent] = pca(data); % Perform PCA on the standardized data
    ve = latent / sum(latent); % Calculate the proportion of explained variance
else
    ve=shape/sum(shape);
end

if ve(1) < 0.6
    i = find(ve > 0.05, 1, 'last');
    if i == length(ve)
        i = i - 1;
    end
    cumve = sum(ve(1:i));
    if cumve < 0.95
        d = false;
    else
        d = true;
    end
else
    d = true;
end

end