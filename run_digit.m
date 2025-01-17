clearvars, close all

load optidigits.mat


tic %get spatial depth-based similarity matrix
[dm_digit, digit_Lmatrix]=rspatial_dp(X,true);
toc %around 250~260s in  matlab R2023a, 11th Gen Intel(R) Core(TM) i7-11700K 



digit_info=getlocalcenter(X,dm_digit,500,'spatial',digit_Lmatrix);

rng(2023)
savek=zeros(20,1);

[digit_dlcc_tmp,digit_dlcc_result]=DLCC(X,dm_digit,digit_info,500,0,'min','knn','K_knn',11,K=10);

%evaluation
%for tmp cluster
digit_temp_cv=cluster2cv(X,digit_dlcc_tmp.temp_clus);
confusionmat(int32(digit_temp_cv(digit_temp_cv~=0)),label(digit_temp_cv~=0)+1)
Misclassification(label(digit_temp_cv~=0),int32(digit_temp_cv(digit_temp_cv~=0)))% 0.0316
adjusted_rand_index(label(digit_temp_cv~=0),int32(digit_temp_cv(digit_temp_cv~=0)))% 0.9434

%for the final result
Misclassification(label+1,digit_dlcc_result.cluster_vector)%0.094
adjusted_rand_index(label,digit_dlcc_result.cluster_vector)%0.8144
                                                           
confusionmat(label+1,int32(digit_dlcc_result.cluster_vector))

