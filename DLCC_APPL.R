
#READ/SIMULATE DATA
#iris
#geyser
data(geyser, package = 'MASS')
#seed (https://archive.ics.uci.edu/ml/datasets/seeds)
seed<-read.csv('Seed_Data.csv')
seed_data<-seed[,-ncol(seed)]
#wine (https://archive.ics.uci.edu/ml/datasets/wine)
wine<-read.table('wine.data',sep=',')
wine_data<-wine[,-1]
wine_label<-wine[,1]


#Simulate data
set.seed(123)
p1<-mlbench.shapes(1000)
p2<-mlbench.spirals(500)
p3<-mlbench.smiley(500)

a=c(1.8,-0.3)
b=c(0,3)
C=c(-0.5,-0.5)
X<-rbind(p1$x + rep(b, each = nrow(p1$x)),p2$x+rep(a,each=nrow(p2$x)),p3$x+rep(C, each = nrow(p3$x)))
plot(X)

data_p3<-as.data.frame(p3$x)
data_X<-as.data.frame(X)

set.seed(123)
ba<-mlbench.cassini(600, relsize=c(3,2,1))
databa<-ba$x%>%as.data.frame()

#similarity matrix (global)
dm0_iris<-depth.dm(iris[,-5],method = 1)
dm0_gey<-depth.dm(geyser,method = 1)
dm_seed<-depth.dm(seed_data,method = 1)
dm_p3<-depth.dm(p3$x,method = 1)
dm_X<-depth.dm(X,method = 1)
dm_ba<-depth.dm(ba$x,method = 1)


#Parameter estimations (For max strategy)
paralist_ci<-para_select_max(databa,sizelist = seq(25,35,5),dm_ba,0.1)
bestpara_ci<-para_select_max2(databa,paralist = paralist_ci$paralist,dm_ba,method = 'maxdep',maxdepth = T)


paralist_smy<-para_select_max(data_p3,sizelist = seq(25,35,5),dm_p3,0.1)
bestpara_smy<-para_select_max2(data_p3,paralist = paralist_smy$paralist,dm_p3,method = 'maxdep',maxdepth = F)

paralist_cpx<-para_select_max(data_X,sizelist = seq(25,30,5),dm_X,0.1)
bestpara_cpx<-para_select_max2(data_X,paralist = paralist_cpx$paralist,dm_X,method = 'mclass',maxdepth = F)


#MIN STRATEGY

#IRIS
iris_real<-as.numeric(iris[,5])
#global similarity matrix DLCC
set.seed(123)
Iris_temp_clus<-DRcluster(iris[,-5],dm0_iris,0.7,50,method = 'min',class_method = 'maxdep',T)
table(Iris_temp_clus$depth_clus$cluster_vector,iris_real)
adjustedRandIndex(Iris_temp_clus$depth_clus$cluster_vector,iris_real)#0.904

#similarity matrix with pre analyses
if_corr(scale(iris[,-5]))

cov_iris<-covMcd(iris[,-5],nsamp = 'deterministic')
covm_iris<-cov_iris$cov
bic_basic<-loglik_gm(iris[,-5],cov_iris$center,covm_iris)*2

iris_de<-decomp_mat(covm_iris)
# lambda_iris*(D%*%A%*%t(D))
iris_EM<-EM_CM(iris[,-5],2:3,150,covm = covm_iris)
class<-sapply(1:nrow(iris), function(x){which.max(iris_EM$z[x,])})
cov.mat.set<-iris_EM$cov
dm_m_iris<-depth.dm_m(iris[,-5],class,cov.mat.set)
set.seed(123)
Iris_temp_clus<-DRcluster(data=iris[,-5],dm0=dm_m_iris,0.7,50,method = 'min',class_method = 'maxdep',T)
table(Iris_temp_clus$depth_clus$cluster_vector,iris_real)
adjustedRandIndex(Iris_temp_clus$depth_clus$cluster_vector,iris_real)

#comparison
#gmm
iris_mc <- Mclust(iris[,-5], G = 3)
adjustedRandIndex(iris_mc$classification,iris_real)#0.904
#pam
iris_pam<-pam(scale(iris[,-5]),3)
adjustedRandIndex(iris_pam$clustering,iris_real)#0.642


#geyser
gey_temp_clus<-DRcluster(geyser,dm0_gey,0.65,90,method = 'min',class_method = 'maxdep',T)
geyser_cl<-cbind(geyser,gey_temp_clus$depth_clus$cluster_vector)
colnames(geyser_cl)<-c('x','y','Cluster')
geyser_cl$Cluster<-factor(geyser_cl$Cluster)
#gmm 
geyser_mc <- Mclust(geyser, G = 3)
table(geyser_cl$Cluster,geyser_mc$classification)

ggplot(geyser_cl, aes(x = x, 
                      y = y, color = Cluster)) + geom_point(show.legend = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(geyser_cl, aes(x = x, 
                      y = y, color = as.factor(geyser_mc$classification))) + geom_point(show.legend = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#seed
set.seed(123)
seed_label<-seed[,ncol(seed)]
seed_temp_clus<-DRcluster(seed_data,dm_seed,Th=0.82,size =70,method = 'min',class_method = 'maxdep',maxdepth = T)
table(seed_temp_clus$depth_clus$cluster_vector,seed_label+1)

adjustedRandIndex(seed_temp_clus$depth_clus$cluster_vector,(seed_label+1))#0.763

#similarity matrix with pre analyses
if_corr(scale(seed_data))
cov_seed<-covMcd(seed_data,nsamp = 'deterministic')
covm_seed<-cov_seed$cov
bic_basic<-loglik_gm(seed_data,cov_seed$center,covm_seed)*2

seed_EM<-EM_CM(seed_data,2:3,150,covm = covm_seed)
class<-sapply(1:nrow(seed_data), function(x){which.max(seed_EM$z[x,])})
cov.mat.set<-seed_EM$cov
dm_m_seed<-depth.dm_m(seed_data,class,cov.mat.set)
set.seed(123)
seed_temp_clus<-DRcluster(seed_data,dm0=dm_m_seed,Th=0.8,size =70,method = 'min',class_method = 'maxdep',maxdepth = T)
table(seed_temp_clus$depth_clus$cluster_vector,seed_label+1)
adjustedRandIndex(seed_temp_clus$depth_clus$cluster_vector,(seed_label+1))#0.787

#compared by other algorithm
seed_mc <- Mclust(seed_data, G = 3)
adjustedRandIndex(seed_mc$classification,seed_label) #0.737
seed_pam<-pam(scale(seed_data),3)
adjustedRandIndex(seed_pam$clustering,seed_label)#0.747

#Italian wine, original data set has too large p, we use the one from UCI with 13 predictors.

if_corr(scale(wine_data))
wine_data<-scale(wine_data)%>%as.data.frame()
M<-Mclust(wine_data,G=2:5,modelNames = c('EEV','EEI'))
M$BIC
M_vari<-sort(M$parameters$variance$shape/sum(M$parameters$variance$shape),decreasing = T)
sum(M_vari[1:max(which(M_vari>0.05))])
#less than 95%, hence Euclidean distance

#normalize
dm_eu<-dist(wine_data,diag = T,upper = T)
# dm_eu<-dm_eu^2
dm_eu<-1/(1+dm_eu)%>%as.matrix()
diag(dm_eu)<-1
wine_temp_clus<-DRcluster(wine_data,dm_eu,Th=0.7,size =50,method = 'min',class_method = 'maxdep',maxdepth = T)
table(wine_temp_clus$depth_clus$cluster_vector,wine_label)
adjustedRandIndex(wine_temp_clus$depth_clus$cluster_vector,wine_label)#0.982

#compare with others
wine_mc <- Mclust(wine_data, G = 3)
adjustedRandIndex(wine_mc$classification,wine_label) #0.930
wine_pam<-pam(scale(wine_data),3)
adjustedRandIndex(wine_pam$clustering,wine_label)#0.741


# Simulated data by #' @author  Francesco Sanna Passino (Ray data)
sim_data<-read.csv('simubyPassino.csv')
zl<-sim_data$zl
sim_data<-sim_data[-ncol(sim_data)]

if_corr(scale(sim_data))#False

M_s<-Mclust(scale(sim_data), G = 3:5, modelNames = c('EEV','EEI'))
#EEV
Ms_vari<-sort(M_s$parameters$variance$shape/sum(M_s$parameters$variance$shape),decreasing = T)
sum(Ms_vari[1:max(which(Ms_vari>0.05))]) #larger than 95%
M_s<-Mclust(sim_data, G = 4, modelNames = c('EEV'))
class<-M_s$classification
cov.mat.set<-M_s$parameters$variance$sigma

cov.mat.set<-lapply(1:dim(cov.mat.set)[3],function(x){cov.mat.set[,,x]})
dm_sim<-depth.dm_m(sim_data,class,cov.mat.set)
set.seed(123)
sim_temp_clus<-DRcluster(sim_data%>%as.data.frame(),dm_sim,Th=0.7,size =250,method = 'min',class_method = 'maxdep',maxdepth = F)
adjustedRandIndex(sim_temp_clus$depth_clus$cluster_vector,zl)#0.823

#plot
V6<-sim_temp_clus$depth_clus$cluster_vector
data_new<-cbind(sim_data,V6)%>%as.data.frame()
data_new$V6<-as.factor(data_new$V6)
ggplot(data_new, aes(x = V1, 
                     y = V2, color = as.factor(zl))) + geom_point(show.legend = F)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data_new, aes(x = V1, 
                     y = V2, color = V6)) + geom_point(show.legend = F)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Comparison
sim_gmm<-Mclust(sim_data, G = 4)
adjustedRandIndex(sim_gmm$classification,zl)#0.780
sim_pam<-pam(sim_data,4)
adjustedRandIndex(sim_pam$clustering,zl)#0.484


#MAX STRATEGY

#Smiley

p3_temp_clus<-DRcluster(data_p3,dm_p3,bestpara_smy$bestpara[2],bestpara_smy$bestpara[1],method = 'max','maxdep',FALSE)

adjustedRandIndex(p3_temp_clus$depth_clus$cluster_vector,p3$classes) #1

data_p3<-cbind(p3$x,p3_temp_clus$depth_clus$cluster_vector)%>%as.data.frame()
colnames(data_p3)<-c('x','y','Cluster')
data_p3$Cluster<-as.factor(data_p3$Cluster)
ggplot(data_p3, aes(x = x, 
                    y = y, color = Cluster)) + geom_point()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#complex shape
X_temp_clus<-DRcluster(data_X,dm_X,bestpara_cpx$bestpara[2],bestpara_cpx$bestpara[1],method = 'max','mclass',FALSE)
levels(p2$classes)<-c('5','6')
levels(p3$classes)<-seq(7,10,1)%>%as.character()
X_label<-c(p1$classes,p2$classes,p3$classes)
adjustedRandIndex(X_temp_clus$depth_clus$cluster_vector,X_label) #0.858
current.label<-rep(0,nrow(data_X))
for (i in 1:length(X_temp_clus$temp.clus)) {
  current.label[X_temp_clus$temp.clus[[i]]]<-i
}

dataX<-as.data.frame(cbind(data_X,X_temp_clus$depth_clus$cluster_vector))
colnames(dataX)<-c('x','y','Cluster')
dataX$Cluster<-factor(dataX$Cluster)

#temp clus
ggplot(dataX, aes(x = x, 
                  y = y, color = factor(current.label))) + geom_point()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#final clus (distance-based classification may be better)
ggplot(dataX, aes(x = x, 
                  y = y, color = Cluster)) + geom_point()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#cassini

ba_temp_clus<-DRcluster(databa,dm_ba,bestpara_ci$bestpara[2],bestpara_ci$bestpara[1],method = 'max','maxdep',T)

data_ba<-as.data.frame(cbind(databa,ba_temp_clus$depth_clus$cluster_vector))
colnames(data_ba)<-c('x','y','Cluster')
data_ba$Cluster<-factor(data_ba$Cluster)
table(ba_temp_clus$depth_clus$cluster_vector,as.numeric(ba$classes))
adjustedRandIndex(ba_temp_clus$depth_clus$cluster_vector,as.numeric(ba$classes)) #0.993

center_label<-unlist(ba_temp_clus$temp.center)
clus.temp<-rep(NA,nrow(data_ba))
for (i in 1:length(ba_temp_clus$temp.clus)) {
  clus.temp[ba_temp_clus$temp.clus[[i]]]<-i
}


ggplot(data_ba, aes(x = x, 
                    y = y, color = Cluster)) + geom_point()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data=data_ba[center_label,],aes(x = x, 
                                             y = y,color=clus.temp[center_label]%>%as.factor(),size=0.5,alpha=0.5),show.legend = FALSE)
