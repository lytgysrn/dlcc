clearvars, close all

load optidigits.mat


tic %get spatial depth-based similarity matrix
[dm_digit, digit_Lmatrix]=rspatial_dp(X,true);
toc %around 250~260s in  matlab R2023a, 11th Gen Intel(R) Core(TM) i7-11700K @ 3.60GHz 3.50 GHz RAM 16GB
digit_Lmatrix=sqrt(digit_Lmatrix);

tic
[digit_rm,digit_rto,digit_dmo]=getlocalcenter(X,dm_digit,500,'spatial',digit_Lmatrix);
toc

rng(2023)
[digit_dlcc_tmp,digit_dlcc_result]=DLCC(X,dm_digit,digit_dmo,digit_rto,digit_rm,500,0,'min','knn','k',10,'K_knn',5,'initial_n',7);

%evaluation
Misclassification(label+1,digit_dlcc_result.cluster_vector)%0.0870
adjusted_rand_index(label,digit_dlcc_result.cluster_vector)%0.8259
confusionmat(label+1,int32(digit_dlcc_result.cluster_vector))
