%read data
clearvars, close all

load spdata.mat


tic 
[dm_spdata, spdata_Lmatrix]=rspatial_dp(X,true);
toc 

tic
sp_info=getlocalcenter(X,dm_spdata,1300,'spatial',spdata_Lmatrix);
toc


[spdata_dlcc_tmp,spdata_dlcc_result]=DLCC(X,dm_spdata,sp_info,1300,0,'min','knn',K=6);


spdata_temp_cv=cluster2cv(X,spdata_dlcc_tmp.temp_clus);
confusionmat(spdata_temp_cv(spdata_temp_cv~=0),label2(spdata_temp_cv~=0))
Misclassification(label2(spdata_temp_cv~=0),(spdata_temp_cv(spdata_temp_cv~=0)))%  0.0881
adjusted_rand_index(label2(spdata_temp_cv~=0),int32(spdata_temp_cv(spdata_temp_cv~=0)))% 0.8598

%for the final result
Misclassification(label2,spdata_dlcc_result.cluster_vector)%0.1868
adjusted_rand_index(label2,spdata_dlcc_result.cluster_vector) %0.6564
confusionmat(label2,spdata_dlcc_result.cluster_vector)


