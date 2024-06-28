function clusters = spec_clus_withsim(sdm, k, method)
%spectral clustering with given similarity mat
%sdm: similarity matrix
%k: the number of clusters
%method: spectral clustering methods

    N = size(sdm, 1);
    sdm(1:N+1:end) = 0;
    MAXiter = 500;
    rep = 50;
    cSsdm = sum(sdm,1);
    cS0_lab = find(cSsdm == 0, 1);
    if ~isempty(cS0_lab)% if there are isolated points, do this to avoid no solution
         diags = sub2ind([N N], cS0_lab, cS0_lab);
         sdm(diags) = 1;
    end
    
    clusters = {};
    % Method 1: Unnormalized Method
    if any(method == 1)
        D = diag(sum(sdm));
        L = D - sdm;
        [~, ~, V] = svd(L);
        f = size(V, 2);
        ker1 = V(:, f-k+1:f);
        group1 = kmeans(ker1, k, 'MaxIter', MAXiter,'replicates', rep);
        clusters{end+1} = group1;
    end
    
    % Method 2: Normalized Cut Method
    if any(method == 2)
        D2 = diag(1 ./ sum(sdm));
        L2 = eye(N) - D2 * sdm;
        [~, ~, V] = svd(L2);
        f = size(V, 2);
        ker2 = V(:, f-k+1:f);
        group2 = kmeans(ker2, k, 'MaxIter', MAXiter, 'replicates', rep);
        clusters{end+1} = group2;
    end
    
    % Method 3: Symmetric Normalized Cut Method
    if any(method == 3)
        D3 = diag(1 ./ sqrt(sum(sdm)));
        L3 = eye(N) - D3 * sdm * D3;
        [~, ~, V] = svd(L3);
        f = size(V, 2);
        ker3 = V(:, f-k+1:f);
        ker3 = bsxfun(@rdivide, ker3, sqrt(sum(ker3.^2, 2)));
        group3 = kmeans(ker3, k, 'MaxIter', MAXiter, 'replicates', rep);
        clusters{end+1} = group3;
    end
    

end
