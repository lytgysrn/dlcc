function result = KNNdep(k, K, dm0, classes)
%KNN classification algorithm based on the depth-based similarity matrix
%Input:
 %k: the number of clusters
 %K: K for the knn algorithm (# of nbr to use)
 %dm0: n*m depth-based similarity matrix, n represents # of points for
 %classification, m represents # of labelled points.
 %classes: classes of labelled points.

%Output:Dmatrix: n*k matrix, storing # of nbrs in each cluster
%       Class_med: classification result

D = size(dm0,1);
Dmatrix = zeros(D,k);
for j = 1:D
    d = dm0(j,:);
    [~,sorted_indices] = sort(d,'descend');
    c = classes(sorted_indices(1:K));
    for i = 1:k
        Dmatrix(j,i) = sum(c==i);
    end
end

Class_med = zeros(1,D);
for x = 1:D
    result=find(Dmatrix(x,:)==max(Dmatrix(x,:))) ;
    % If there is a tie, select the nearest nbr among obs in the tied clusters.
    length_result=length(result);
    if length_result > 1
        maxsim2clus = arrayfun(@(i) max(dm0(x, (classes==result(i)))), 1:length_result);
        [~, max_index] = max(maxsim2clus);
        result = result(max_index);
    end
    Class_med(x) = result;
end
result = struct('Dmatrix',Dmatrix,'class',Class_med);
end