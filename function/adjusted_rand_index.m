function ARI = adjusted_rand_index(S, T)
%compute Adjusted Rand index (clustering comparison metric)

    % Find number of clusters
    clustersS = unique(S);
    clustersT = unique(T);
    numClustS = length(clustersS);
    numClustT = length(clustersT);

    % Initialize contingency table
    n = length(S);
    Contingency = zeros(numClustS, numClustT);
    for i=1:numClustS
        for j=1:numClustT
            Contingency(i,j) = sum(S == clustersS(i) & T == clustersT(j));
        end
    end

    % Compute ARI
    comb2n = @(n) arrayfun(@(x) x*(x-1)/2, n);
    N = comb2n(n);
    N1 = sum(comb2n(sum(Contingency, 2)));
    N2 = sum(comb2n(sum(Contingency, 1)));
    N3 = sum(sum(comb2n(Contingency)));
    ARI = (N*N3 - N1*N2) / (N*(N1+N2)/2 - N1*N2);
end
