%read data
clearvars, close all

load pendigit.mat

tic
[dm_pendigi, pendigi_L]=rspatial_dp(X,true);
toc

%202.4 seconds in  matlab R2023a, 11th Gen Intel(R) Core(TM)
%i7-11700K @ 3.60GHz   3.50 GHz RAM 16GB
pendigi_info=getlocalcenter(X,dm_pendigi,1000,'spatial',pendigi_L);

rng(2024)
[pendigi_dlcc_tmp,pendigi_dlcc_result]=DLCC(X,dm_pendigi,pendigi_info,1000,0,'min','knn',K=10);
%evaluation
%for tmp cluster
pendigi_temp_cv=cluster2cv(X,pendigi_dlcc_tmp.temp_clus);
confusionmat(label(pendigi_temp_cv~=0)+1,(pendigi_temp_cv(pendigi_temp_cv~=0)))
Misclassification(label(pendigi_temp_cv~=0)+1,(pendigi_temp_cv(pendigi_temp_cv~=0)))% 0.0703
adjusted_rand_index(label(pendigi_temp_cv~=0)+1,(pendigi_temp_cv(pendigi_temp_cv~=0)))%0.8616

%for the final result
Misclassification(label+1,pendigi_dlcc_result.cluster_vector)%  0.1822
adjusted_rand_index(label+1,pendigi_dlcc_result.cluster_vector)%0.6661
confusionmat(label+1,(pendigi_dlcc_result.cluster_vector))


%maha version
if_corr(X)
rng(2023)
bwmodel=EM_EEV(X,6:10,500);
if_corr(X,bwmodel.shape)
dmED_pen=Maha_dmcreator(X,'Euclidean');
pendigiMD_info=getlocalcenter(X,dmED_pen,1000,'mahalanobis');

rng(2023)
[pendigi_dlcc_tmpMD,pendigi_dlccMD_result]=DLCC(X,dmED_pen,pendigiMD_info,1000,0,'min','knn','depth','mahalanobis',K_knn=6);

%for the final result
Misclassification(label+1,pendigi_dlccMD_result.cluster_vector)%0.2052
adjusted_rand_index(label+1,pendigi_dlccMD_result.cluster_vector)%0.6488
confusionmat(label+1,(pendigi_dlccMD_result.cluster_vector))

