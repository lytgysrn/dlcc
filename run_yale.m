%read data
clearvars, close all

load yale.mat

tic
[dm_yale, yale_L]=rspatial_dp(X,false);
toc
%202.4 seconds in  matlab R2023a, 11th Gen Intel(R) Core(TM)
%i7-11700K @ 3.60GHz   3.50 GHz RAM 16GB
yale_L=sqrt(yale_L);
tic
[yale_rm,yale_rto,yale_dmo]=getlocalcenter(X,dm_yale,200,'spatial',yale_L);


rng(2023)
[yale_dlcc_tmp,yale_dlcc_result]=DLCC(X,dm_yale,yale_dmo,yale_rto,yale_rm,200,0,'min','rf','k',10);
toc

%evaluation
Misclassification(label,yale_dlcc_result.cluster_vector)%0.0065
adjusted_rand_index(label,yale_dlcc_result.cluster_vector)%0.9856
confusionmat(label,int32(yale_dlcc_result.cluster_vector))




% Below are steps by steps example (same as applying DLCC function directly)

yale_a=filter_center(dm_yale,yale_dmo,200,yale_rto,yale_rm);
yale_temp_cl=get_temp_cl_WK(yale_a,200,yale_dmo,10);

a = cell2mat(yale_temp_cl);
s=200;
Kclus = length(yale_temp_cl);
temp_clus = cell(1,Kclus);
for i = 1:Kclus
    [C,~,~] = unique(yale_dmo(1:s,yale_temp_cl{i}));
    temp_clus{i} = C;
end

temp_clus=get_temp_cluster(a,dm_yale,temp_clus,Kclus, yale_temp_cl,s, 'min');

% we try 100 different seeds for random forest to see the randomness influence.
check=loop_rfdlcc(X,temp_clus,label,100);
check_table=array2table(check);
summary(check_table)
    % 
    %    Values:
    % 
    %         Min        0.0025 
    %         Median    0.00625 
    %         Max         0.011 
    % 
    % check2: 100×1 double
    % 
    %     Values:
    % 
    %         Min       0.97629 
    %         Median    0.98632 
    %         Max       0.99449 



