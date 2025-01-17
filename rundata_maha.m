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
iris_info=getlocalcenter(meas,dmMD_iris,50,'mahalanobis');
[iris_dlcc_tmp,iris_dlcc_result]=DLCC(meas,dmMD_iris,iris_info,50,0,'min','maxdep','depth','mahalanobis','maxdepth',true,K=3);


%for the final result
Misclassification(species,iris_dlcc_result.cluster_vector)%0.0333
adjusted_rand_index(species,iris_dlcc_result.cluster_vector)%0.9039
confusionmat(species,iris_dlcc_result.cluster_vector)


%seed
% seed=zscore(seed);
if_corr(seed)
rng(2023)
dmMD_seed=Maha_dmcreator(seed,'MCD');
seed_info=getlocalcenter(seed,dmMD_seed,65,'mahalanobis');

[seed_dlcc_tmp,seed_dlcc_result]=DLCC(seed,dmMD_seed,seed_info,65,0,'min','maxdep','depth','mahalanobis','maxdepth',true,'ifloop',true,K=3);


%for the final result
Misclassification(seed_label+1,seed_dlcc_result.cluster_vector)%0.0857
adjusted_rand_index(seed_label+1,seed_dlcc_result.cluster_vector)%0.7612
confusionmat(seed_label+1,seed_dlcc_result.cluster_vector)

%wine
if_corr(wine)
rng(2023)
wine=zscore(wine);
bwmodel=EM_EEV(wine,2:4,110);
if_corr(wine,bwmodel.shape)
%still not pass the test, use Euclidean distance to build the sim mat
dmED_wine=Maha_dmcreator(wine,'Euclidean');
wine_info=getlocalcenter(wine,dmED_wine,50,'mahalanobis');
[wine_dlcc_tmp,wine_dlcc_result]=DLCC(wine,dmED_wine,wine_info,50,0,'min','maxdep','ifloop',true,'maxdepth',true,'depth','mahalanobis');



%for the final result
Misclassification(wine_label,wine_dlcc_result.cluster_vector)%0.0056
adjusted_rand_index(wine_label,wine_dlcc_result.cluster_vector)%0.9817
confusionmat(wine_label,wine_dlcc_result.cluster_vector)


%simu_data
if_corr(simu_data)
rng(2023)
bmodel=EM_EEV(simu_data,3:5,100);
if_corr(simu_data,bmodel.shape)
dmMD_simu=Maha_dmcreator(simu_data,'provided',bmodel.classes,bmodel.cov);
simu_info=getlocalcenter(simu_data,dmMD_simu,250,'mahalanobis');
[simu_dlcc_tmp,simu_dlcc_result]=DLCC(simu_data,dmMD_simu,simu_info,250,0,'min','maxdep','depth','mahalanobis','ifloop',true);


%for the final result
Misclassification(simu_label,simu_dlcc_result.cluster_vector)%0.07
adjusted_rand_index(simu_label,simu_dlcc_result.cluster_vector)%0.8261
confusionmat(simu_label,simu_dlcc_result.cluster_vector)


