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
[iris_rm,iris_rto,iris_dmo]=getlocalcenter(meas,dm_iris,50);
rng(2023)
[iris_dlcc_tmp,iris_dlcc_result]=DLCC(meas,dm_iris,iris_dmo,iris_rto,iris_rm,50,0.62,'min','rf');

%evaluation

Misclassification(species,iris_dlcc_result.cluster_vector)% 0.0400
adjusted_rand_index(species,iris_dlcc_result.cluster_vector)%0.8858
confusionmat(species,iris_dlcc_result.cluster_vector)


%seed
dm_seed=rspatial_dp(seed);
[seed_rm,seed_rto,seed_dmo]=getlocalcenter(seed,dm_seed,63);
rng(2023)
[seed_dlcc_tmp,seed_dlcc_result]=DLCC(seed,dm_seed,seed_dmo,seed_rto,seed_rm,63,0.3,'min','rf');
%evaluation

Misclassification(seed_label+1,seed_dlcc_result.cluster_vector)%0.0857
adjusted_rand_index(seed_label+1,seed_dlcc_result.cluster_vector)%0.7624
confusionmat(seed_label+1,seed_dlcc_result.cluster_vector)



%wine
%since spatial depth is not affine invariant, for wine data set, the scales
%among variables are very different, normalization is necessary.
wine=zscore(wine);
dm_wine=rspatial_dp(wine);
[wine_rm,wine_rto,wine_dmo]=getlocalcenter(wine,dm_wine,30);
rng(2023)
[wine_dlcc_md,wine_dlcc_mdr]=DLCC(wine,dm_wine,wine_dmo,wine_rto,wine_rm,30,0,'min','maxdep');
%evaluation


Misclassification(wine_label,wine_dlcc_mdr.cluster_vector)%0.0281
adjusted_rand_index(wine_label,wine_dlcc_mdr.cluster_vector)% 0.9125
confusionmat(wine_label,wine_dlcc_mdr.cluster_vector)

%simu_data
[dm_simu_data, simu_data_Lmatrix]=rspatial_dp(simu_data);
simu_data_Lmatrix=sqrt(simu_data_Lmatrix);
[simu_data_rm,simu_data_rto,simu_data_dmo]=getlocalcenter(simu_data,dm_simu_data,160,'spatial',simu_data_Lmatrix);
rng(2023)
[simu_data_dlcc_tmp,simu_data_dlcc_result]=DLCC(simu_data,dm_simu_data,simu_data_dmo,simu_data_rto,simu_data_rm,160,0,'min','rf','k',4);

%evaluation

Misclassification(simu_label,simu_data_dlcc_result.cluster_vector)%0.0790
adjusted_rand_index(simu_label,simu_data_dlcc_result.cluster_vector)%0.8060
confusionmat(simu_label,simu_data_dlcc_result.cluster_vector)

