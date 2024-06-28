clearvars, close all

load data\ba.mat
load data\blend.mat
load data\smiley.mat

%smiley
rng(2023)
dmMD_smiley=Maha_dmcreator(smiley,'MCD');
[sm_rm,sm_rto,sm_dmo]=getlocalcenter(smiley,dmMD_smiley,25,'mahalanobis');
[sm_dlcc_tmp,sm_dlcc_result]=DLCC(smiley,dmMD_smiley,sm_dmo,sm_rto,sm_rm,25,0.3,'max','maxdep','depth','mahalanobis');

%evaluation
%for tmp cluster
sm_temp_cv=cluster2cv(smiley,sm_dlcc_tmp.temp_clus);
confusionmat(smiley_label(sm_temp_cv~=0),sm_temp_cv(sm_temp_cv~=0))
%for the final result
Misclassification(smiley_label,sm_dlcc_result.cluster_vector)
adjusted_rand_index(smiley_label,sm_dlcc_result.cluster_vector)
confusionmat(smiley_label,sm_dlcc_result.cluster_vector)


%blend
rng(2023)
dmMD_blend=Maha_dmcreator(blend,'MCD');
[blend_rm,blend_rto,blend_dmo]=getlocalcenter(blend,dmMD_blend,25,'mahalanobis');
[blend_dlcc_tmp,blend_dlcc_result]=DLCC(blend,dmMD_blend,blend_dmo,blend_rto,blend_rm,25,0,'max','knn','depth','mahalanobis','K_knn',3);

%evaluation
%for tmp cluster
blend_temp_cv=cluster2cv(blend,blend_dlcc_tmp.temp_clus);
confusionmat(blend_label(blend_temp_cv~=0),blend_temp_cv(blend_temp_cv~=0))
%for the final result
Misclassification(blend_label,blend_dlcc_result.cluster_vector)
adjusted_rand_index(blend_label,blend_dlcc_result.cluster_vector)
confusionmat(blend_label,blend_dlcc_result.cluster_vector)

%cassini
rng(2023)
dmMD_ba=Maha_dmcreator(ba,'MCD');
[ba_rm,ba_rto,ba_dmo]=getlocalcenter(ba,dmMD_ba,25,'mahalanobis');
[ba_dlcc_tmp,ba_dlcc_result]=DLCC(ba,dmMD_ba,ba_dmo,ba_rto,ba_rm,25,0.4,'max','maxdep','depth','mahalanobis','maxdepth',true);

%evaluation
%for tmp cluster
ba_temp_cv=cluster2cv(ba,ba_dlcc_tmp.temp_clus);
confusionmat(ba_label(ba_temp_cv~=0),ba_temp_cv(ba_temp_cv~=0))
%for the final result
Misclassification(ba_label,ba_dlcc_result.cluster_vector)
adjusted_rand_index(ba_label,ba_dlcc_result.cluster_vector)
confusionmat(ba_label,ba_dlcc_result.cluster_vector)


%spatial depth case
rng(2023)
dmSD_smiley=rspatial_dp(smiley);
[sm_rm,sm_rto,sm_dmo]=getlocalcenter(smiley,dmSD_smiley,30);
[sm_dlcc_tmp,sm_dlcc_result]=DLCC(smiley,dmSD_smiley,sm_dmo,sm_rto,sm_rm,30,0,'max','knn');

%evaluation
%for tmp cluster
sm_temp_cv=cluster2cv(smiley,sm_dlcc_tmp.temp_clus);
confusionmat(smiley_label(sm_temp_cv~=0),sm_temp_cv(sm_temp_cv~=0))
%for the final result
Misclassification(smiley_label,sm_dlcc_result.cluster_vector)
adjusted_rand_index(smiley_label,sm_dlcc_result.cluster_vector)
confusionmat(smiley_label,sm_dlcc_result.cluster_vector)


%blend
rng(2023)
[dmSD_blend,blend_L]=rspatial_dp(blend);
blend_L=sqrt(blend_L);
[blend_rm,blend_rto,blend_dmo]=getlocalcenter(blend,dmSD_blend,22,'spatial',blend_L);
[blend_dlcc_tmp,blend_dlcc_result]=DLCC(blend,dmSD_blend,blend_dmo,blend_rto,blend_rm,22,0,'max','knn','K_knn',3);

%evaluation
%for tmp cluster
blend_temp_cv=cluster2cv(blend,blend_dlcc_tmp.temp_clus);
confusionmat(blend_label(blend_temp_cv~=0),blend_temp_cv(blend_temp_cv~=0))
%for the final result
Misclassification(blend_label,blend_dlcc_result.cluster_vector)
adjusted_rand_index(blend_label,blend_dlcc_result.cluster_vector)
confusionmat(blend_label,blend_dlcc_result.cluster_vector)

%cassini
rng(2023)
[dmSD_ba,ba_L]=rspatial_dp(ba);
[ba_rm,ba_rto,ba_dmo]=getlocalcenter(ba,dmSD_ba,20,'spatial');
[ba_dlcc_tmp,ba_dlcc_result]=DLCC(ba,dmSD_ba,ba_dmo,ba_rto,ba_rm,20,0.28,'max','knn','K_knn',3);

%evaluation
%for tmp cluster
ba_temp_cv=cluster2cv(ba,ba_dlcc_tmp.temp_clus);
confusionmat(ba_label(ba_temp_cv~=0),ba_temp_cv(ba_temp_cv~=0))
%for the final result
Misclassification(ba_label,ba_dlcc_result.cluster_vector)
adjusted_rand_index(ba_label,ba_dlcc_result.cluster_vector)
confusionmat(ba_label,ba_dlcc_result.cluster_vector)


%Micro-array data
%leukemia
X=load('data\leukemia_x.txt');
Class=load('data\leukemia_y.txt');
X=X';
[dm_lc, lc_Lmatrix]=rspatial_dp(X,true);
lc_Lmatrix=sqrt(lc_Lmatrix);
[lc_rm,lc_rto,lc_dmo]=getlocalcenter(X,dm_lc,10,'spatial',lc_Lmatrix);

rng(2023)
[lc_dlcc_tmp,lc_dlcc_result]=DLCC(X,dm_lc,lc_dmo,lc_rto,lc_rm,10,0.1,'max','maxdep');
adjusted_rand_index(Class,lc_dlcc_result.cluster_vector)%0.8897
confusionmat(Class+1,lc_dlcc_result.cluster_vector)
Misclassification(Class+1,lc_dlcc_result.cluster_vector)

%sucancer
X=load('data\SuCancer.txt');
Class=load('data\SuCancer_y.txt');
X=X';

[dm_sc, sc_Lmatrix]=rspatial_dp(X,true);
sc_Lmatrix=sqrt(sc_Lmatrix);


[sc_rm,sc_rto,sc_dmo]=getlocalcenter(X,dm_sc,40,'spatial',sc_Lmatrix);

rng(2023)
[sc_dlcc_tmp,sc_dlcc_result]=DLCC(X,dm_sc,sc_dmo,sc_rto,sc_rm,40,0.54,'max','knn');
adjusted_rand_index(Class,sc_dlcc_result.cluster_vector)%0.1154
confusionmat(Class+1,sc_dlcc_result.cluster_vector)
Misclassification(Class+1,sc_dlcc_result.cluster_vector)

%
% Example about presenting graphs, though I prefer to use ggplot in R
% figure; 
% hold on; 
% colors = unique(blend_dlcc_result.cluster_vector); 
% for i = 1:length(colors)   
%     scatter(blend(blend_dlcc_result.cluster_vector == colors(i),1), blend(blend_dlcc_result.cluster_vector == colors(i),2), 'filled');
% end
% set(gca, 'Box', 'on', 'XGrid', 'off', 'YGrid', 'off');
% xlabel('x');
% ylabel('y');
% hold off;