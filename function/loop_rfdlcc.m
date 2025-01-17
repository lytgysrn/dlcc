function c_p = loop_rfdlcc(X, temp_clus, label,ntrees)
%run 100 times random forest classification with different seeds, write it as a function for
%convenience.

seed_set = 1:100;
c_p = zeros(100,2);
for i = 1:length(seed_set)
    rng(seed_set(i)); % set the seed for reproducibility
    cluster_result=left_class(X,temp_clus,'rf','ntrees',ntrees);
    c_p(i,1) = Misclassification(cluster_result.cluster_vector, label);
    c_p(i,2) =adjusted_rand_index(cluster_result.cluster_vector, label);
end
end