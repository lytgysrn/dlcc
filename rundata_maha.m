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
if_corr(meas)
rng(2023)
dmMD_iris=Maha_dmcreator(meas,'MCD');
[iris_rm,iris_rto,iris_dmo]=getlocalcenter(meas,dmMD_iris,50,'mahalanobis');
[iris_dlcc_tmp,iris_dlcc_result]=DLCC(meas,dmMD_iris,iris_dmo,iris_rto,iris_rm,50,0.72,'min','maxdep','depth','mahalanobis','maxdepth',true);

%evaluation


Misclassification(species,iris_dlcc_result.cluster_vector)%0.0333
adjusted_rand_index(species,iris_dlcc_result.cluster_vector)%0.9039
confusionmat(species,iris_dlcc_result.cluster_vector)


%seed
% seed=zscore(seed);
if_corr(seed)
rng(2023)
dmMD_seed=Maha_dmcreator(seed,'MCD');
[seed_rm,seed_rto,seed_dmo]=getlocalcenter(seed,dmMD_seed,70,'mahalanobis');
[seed_dlcc_tmp,seed_dlcc_result]=DLCC(seed,dmMD_seed,seed_dmo,seed_rto,seed_rm,70,0.7,'min','knn','depth','mahalanobis','maxdepth',true,'ifloop',false,'K_knn',6);


Misclassification(seed_label+1,seed_dlcc_result.cluster_vector)%0.081
adjusted_rand_index(seed_label+1,seed_dlcc_result.cluster_vector)%0.7711
confusionmat(seed_label+1,seed_dlcc_result.cluster_vector)


%wine
if_corr(wine)
rng(2023)
wine=zscore(wine);
bwmodel=EM_EEV(wine,2:4,110);
if_corr(wine,bwmodel.shape)
%still not pass the test, use Euclidean distance to build the sim mat
dmED_wine=Maha_dmcreator(wine,'Euclidean');
[wine_rm,wine_rto,wine_dmo]=getlocalcenter(wine,dmED_wine,50,'mahalanobis');
[wine_dlcc_tmp,wine_dlcc_result]=DLCC(wine,dmED_wine,wine_dmo,wine_rto,wine_rm,50,0.3,'min','maxdep','ifloop',true,'maxdepth',true,'depth','mahalanobis');


%evaluation

Misclassification(wine_label,wine_dlcc_result.cluster_vector)%0.0056
adjusted_rand_index(wine_label,wine_dlcc_result.cluster_vector)%0.9817
confusionmat(wine_label,wine_dlcc_result.cluster_vector)


%simu_data
if_corr(simu_data)
rng(2023)
bmodel=EM_EEV(simu_data,3:5,100);
if_corr(simu_data,bmodel.shape)
dmMD_simu=Maha_dmcreator(simu_data,'provided',bmodel.classes,bmodel.cov);
[simu_rm,simu_rto,simu_dmo]=getlocalcenter(simu_data,dmMD_simu,250,'mahalanobis');
[simu_dlcc_tmp,simu_dlcc_result]=DLCC(simu_data,dmMD_simu,simu_dmo,simu_rto,simu_rm,250,0.33,'min','maxdep','depth','mahalanobis','ifloop',true);

%evaluation

Misclassification(simu_label,simu_dlcc_result.cluster_vector)%0.073
adjusted_rand_index(simu_label,simu_dlcc_result.cluster_vector)%0.8202
confusionmat(simu_label,simu_dlcc_result.cluster_vector)


