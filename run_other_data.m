clearvars, close all

% load data

load fisheriris
species=grp2idx(species);

wine = readmatrix('wine.data', 'FileType', 'text', 'Delimiter', ',');
wine_label=wine(:,1);
wine(:,1)=[];

seed = readmatrix('seed_Data.csv');
seed_label=seed(:,end);
seed(:,end)=[];

simu_data=readmatrix('simubyPassino.csv');
simu_label=simu_data(:,end);
simu_data(:,end)=[];

%iris
dm_iris=rspatial_dp(meas);
%dm_iris_approxi=rspatial_dp_approxi(meas,false,100);

iris_info=getlocalcenter(meas,dm_iris,55);
rng(2024)
[iris_dlcc_tmp,iris_dlcc_result]=DLCC(meas,dm_iris,iris_info,55,0,'min','rf',ifloop=false);

%for the final result
Misclassification(species,iris_dlcc_result.cluster_vector)% 0.0400
adjusted_rand_index(species,iris_dlcc_result.cluster_vector)%0.8858
confusionmat(species,iris_dlcc_result.cluster_vector)


%seed
dm_seed=rspatial_dp(seed);
seed_info=getlocalcenter(seed,dm_seed,63);
[seed_dlcc_tmp,seed_dlcc_result]=DLCC(seed,dm_seed,seed_info,63,0,'min','knn');
%evaluation

%for the final result
Misclassification(seed_label+1,seed_dlcc_result.cluster_vector)%0.0857
adjusted_rand_index(seed_label+1,seed_dlcc_result.cluster_vector)%0.7626
confusionmat(seed_label+1,seed_dlcc_result.cluster_vector)



%wine
%since spatial depth is not affine invariant, for wine data set, the scales
%among variables are very different, normalization is necessary.
wine=zscore(wine);
dm_wine=rspatial_dp(wine);
wine_info=getlocalcenter(wine,dm_wine,50);
[wine_dlcc_md,wine_dlcc_mdr]=DLCC(wine,dm_wine,wine_info,50,0,'min','knn',K_knn=37);
%evaluation


Misclassification(wine_label,wine_dlcc_mdr.cluster_vector)%0.0225
adjusted_rand_index(wine_label,wine_dlcc_mdr.cluster_vector)% 0.9295
confusionmat(wine_label,wine_dlcc_mdr.cluster_vector)

%simu_data
[dm_simu_data, simu_data_Lmatrix]=rspatial_dp(simu_data);
simu_info=getlocalcenter(simu_data,dm_simu_data,165,'spatial',simu_data_Lmatrix);
rng(2024)
[simu_data_dlcc_tmp,simu_data_dlcc_result]=DLCC(simu_data,dm_simu_data,simu_info,165,0,'min','rf');


%for the final result
Misclassification(simu_label,simu_data_dlcc_result.cluster_vector)%0.1
adjusted_rand_index(simu_label,simu_data_dlcc_result.cluster_vector)%0.7521
confusionmat(simu_label,simu_data_dlcc_result.cluster_vector)

