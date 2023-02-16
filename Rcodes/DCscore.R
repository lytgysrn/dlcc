library(fpc)
library(tidyverse)
library(cluster)

#simulate a circle
simu.circle<-function (x,y,r,n){
  theta<-runif(n,min=0,max = 2*pi)
  x_pos<-r*cos(theta)+x
  y_pos<-r*sin(theta)+y
  return(cbind(x_pos,y_pos))
}

#cl.loop. AssignCluster are all inner functions in DBCA
cl.loop<-function(k,cluster,adj,cl){
  if (cluster[k]==-1){
    cluster[k]<-cl
    if (length(adj[[k]]!=0)){
      for (w in 1:length(adj[[k]])) {
        new.k<-adj[[k]][w]
        cluster<-cl.loop(new.k,cluster,adj,cl)
      } 
    }
  }
  return(cluster)
}
AssignCluster<-function(adj){
  cluster<-rep(-1,length(adj))
  cl<-0
  for (i in 1:length(cluster)) {
    if (cluster[i]==-1){
      if (length(adj[[i]])==0){
        cluster[i]<-0
      } else {
        cl<-cl+1
        cluster<-cl.loop(i,cluster,adj,cl)
      }
    }
  }
  return(cluster)
}
DBCA<-function(dm,theta,mn){
  #dm:similarity matrix (global)
  #theta: similarity threshold
  #mn:minimum obs in a DBCA-based cluster, which is not the same as the minPt in DBSCAN.
  if (missing(mn)){
    mn<-3
  }
  adj<-list()
  for (i in 1:nrow(dm)) {
    label<-which(dm[i,]>=theta)
    label<-setdiff(label,i)
    adj[[i]]<-label
  }
  cluster<-AssignCluster(adj)
  cls.n<-unique(cluster)%>%setdiff(c(0))
  for (j in cls.n) {
    pos<-which(cluster==cls.n[j])
    if (length(pos)<mn){
      cluster[pos]<-0
    }
  }
  return(cluster)
}

#DBCAmerge: For find the maximum theta that connect all obs in the cluster
DBCAmerge<-function(dbca,dm0.save){
  #dbca:result of dbca
  #dm0.save: similarity matrix
  if(length(dbca[dbca==0])!=0){
    dbca[dbca==0]<-100:(100+length(dbca[dbca==0])-1)
  }
  
  ab<-unique(dbca)
  temp.K<-length(ab)
  gap_bet_tk<-matrix(0,temp.K,temp.K)
  all_pos<-lapply(1:temp.K,function(x){which(dbca==ab[x])})
  
  for (y in 1:(temp.K-1)) {
    for (u in (y+1):temp.K) {
      c1<-all_pos[[y]]
      c2<-all_pos[[u]]
      max_sim<-max(dm0.save[c2,c1])
      gap_bet_tk[u,y]<-max_sim
      gap_bet_tk[y,u]<-max_sim
    }
  }
  diag(gap_bet_tk)<-1
  gapset<-gap_bet_tk%>%as.numeric()
  gapset<-gapset[gapset!=0]%>%sort(decreasing = T)%>%unique()
  gapset<-gapset[-1]
  
  m=1
  gs<-c()
  while(is_empty(gs)==T){
    v=1
    while (v <= temp.K) {
      if (gapset[m]%in%gap_bet_tk[,v]){
        comb_lab<-which(gap_bet_tk[,v]==gapset[m])
        gap_bet_tk[1:temp.K,v]<-sapply(1:temp.K,function(x){max(gap_bet_tk[x,c(v,comb_lab)])})
        gap_bet_tk[v,1:temp.K]<-gap_bet_tk[1:temp.K,v]
        gap_bet_tk<-gap_bet_tk[-comb_lab,-comb_lab]
        temp.K<-temp.K-length(comb_lab)
        if (is.matrix(gap_bet_tk)==FALSE){
          gap_bet_tk<-gap_bet_tk%>%as.matrix()
        }
      } else {
        v=v+1
      }
    }
    if (length(gap_bet_tk)==1){
      gs<-gapset[m]
    } else {
      m=m+1
    }
  }
  return(list(gs=gs,gapset=gapset))
}

#NC_score: Density based clustering score.
NC_score<-function(dm0,cl_label,ignore_out,mn){
  if (missing(mn)){
    mn<-3
  }
  
  Kest=length(unique(cl_label))
  mins<-c()
  ss<-c()
  os<-c()
  Nw<-c()
  dm0.save<-dm0
  
  dm0.sort<-sapply(1:dim(dm0.save)[1],function(i){sort(dm0.save[i,],decreasing = T)})
  dbcaot<-DBCA(dm0.save,min(dm0.sort[2,]))
  
  temp.K<-length(unique(dbcaot))
  if (temp.K==1){
    gs<-min(dm0.sort[2,])
  } else {
    gs<-DBCAmerge(dbcaot,dm0.save)$gs
  }
  
  for (w in 1:Kest){
    lc<-which(cl_label==w)
    Nw[w]<-length(lc)
    dm.c<-dm0.save[lc,lc]
    dm.c.order<-sapply(1:Nw[w],function(i){sort(dm.c[i,],decreasing = T)})
    ms<-min(dm.c.order[2,])
    DBCAclr<-DBCA(dm.c,ms,mn)
    ldbca<-length(unique(DBCAclr))
    if (ignore_out==T){
      if(ldbca==1){
        DBCAclr<-DBCA(dm.c,ms+0.0000000001,mn)
      }
      while(identical(unique(DBCAclr)%>%sort(),c(0,1))){
        ms<-min(dm.c.order[2,][dm.c.order[2,]>ms])
        DBCAclr<-DBCA(dm.c,ms,mn)
        ldbca<-length(unique(DBCAclr))
      }
      if (ldbca==1){
        ldbca=2
      }
    }
    
    if (ldbca!=1){
      if (ignore_out==T){
        save_cut<-min(dm.c.order[2,][dm.c.order[2,]>ms])
        DBCAclr<-DBCA(dm.c,save_cut)
      }
      dbcams<-DBCAmerge(DBCAclr,dm.c)
      ms<-dbcams$gs
      if (ignore_out==T){
        ldbca<-1
        if(ldbca==1){
          new_pos<-which(dbcams$gapset==dbcams$gs)-1
          DBCAclr<-DBCA(dm.c,dbcams$gapset[new_pos],mn)
          while(identical(unique(DBCAclr)%>%sort(),c(0,1))){
            new_pos<-new_pos-1
            ms<-dbcams$gapset[new_pos]
            DBCAclr<-DBCA(dm.c,ms,mn)
            ldbca<-length(unique(DBCAclr))
          }
        }
        ms<-dbcams$gapset[new_pos+1]
        DBCAclr<-DBCA(dm.c,ms,mn)
      }
    } 
    mins[w]=ms
    diag(dm0)<-0
    outlier<-which(DBCAclr==0)
    lengtho<-length(outlier)
    if (lengtho>0){
      x_length<-Nw[w]-lengtho
      lc_ro<-lc[-outlier]
      wc1<-lapply(1:x_length,function(x){dm0[lc_ro[x],lc][dm0[lc_ro[x],lc]>=ms]})
      o_lab<-lc[outlier]
      wc_o<-sapply(1:lengtho,function(x){max(dm0[o_lab[x],lc_ro])})
    } else {
      wc1<-lapply(1:Nw[w],function(x){dm0[lc[x],lc][dm0[lc[x],lc]>=ms]})
      wc_o<-c()
    }
    # mwc<-sapply(1:Nw[w],function(x){mean(wc1[[x]])})
    ss[w]<-mean(c(unlist(wc1),wc_o))#with in cluster sim
    
    
    wc2<-sapply(1:Nw[w],function(x){max(dm0[lc[x],-lc])})
    maxbc<-max(wc2)
    if (maxbc>=gs){
      overms_lab<-which(wc2>=gs)
      wc2_s<-lapply(1:length(overms_lab),function(x){dm0[lc[overms_lab[x]],-lc][dm0[lc[overms_lab[x]],-lc]>=gs]})
      oth_s<-rep(gs,Nw[w]-length(wc2_s))
      os[w]<-mean(c(oth_s,unlist(wc2_s)))
    }else{
      os[w]<-gs
    }
  }
  fs<-ss*Nw/os
  score<-sum(fs)
  return(list(score=score,fs=fs,iclus_s=ss,bclus_s=os,ms=mins,gs=gs))
}


#…………………………………………………………………………………………………………………………………………

#RUN LATER EXAMPLES AFTER GENERATING SIMULATED DATA SETS IN DLCC_APPL
#SmileyF (I did not set seed for this simulation before, so the results may be a little bit different)
set.seed(123)
circle<-simu.circle(0,0,1.7,1000)
data_graph<-rbind(circle,p3$x)%>%as.data.frame()
realgra<-c(rep(5,1000),p3$classes)

#dlcc
dm_gra<-depth.dm(data_graph,method = 1)
paralist_gra<-para_select_max(data_graph,sizelist = seq(40,50,5),dm_gra,0.1)
bestpara_gra<-para_select_max2(data_graph,paralist = paralist_gra$paralist,dm_gra,method = 'mclass',maxdepth = F)
gra<-DRcluster(data_graph,dm_gra,bestpara_gra$bestpara[2],bestpara_gra$bestpara[1],method = 'max','mclass',FALSE)

#GMM with G=5 and 16
gmm0<-Mclust(data_graph,G=5)
gmm_g<-Mclust(data_graph,G=16)
#DBCA
DBCAcl<-DBCA(dm_gra,theta = 0.98)

#real label
ggplot(data_graph, aes(x = x_pos, 
                       y = y_pos, color = realgra%>%as.factor())) + 
  geom_point(show.legend = FALSE)+   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))
#dlcc
ggplot(data_graph, aes(x = x_pos, 
                       y = y_pos, color = gra$depth_clus$cluster_vector%>%as.factor())) + 
  geom_point(show.legend = FALSE)+
  # geom_point(data=data_graph[a,],aes(x = x_pos,y=y_pos,color='',alpha=0.5,size=0.5),show.legend = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#G=5
ggplot(data_graph, aes(x = x_pos, 
                       y = y_pos, color = gmm0$classification%>%as.factor())) + 
  geom_point(show.legend = FALSE)+   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))

#G=16
ggplot(data_graph, aes(x = x_pos, 
                       y = y_pos, color = gmm_g$classification%>%as.factor())) + 
  geom_point(show.legend = FALSE)+   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))

#DBCA
ggplot(data_graph[-which(DBCAcl==0),], aes(x = x_pos, 
                                           y = y_pos))+geom_point(aes(color = DBCAcl[-which(DBCAcl==0)]%>%as.factor()),show.legend = FALSE)+geom_point(data=data_graph[which(DBCAcl==0),],aes(x=x_pos,y=y_pos,alpha=0.5),show.legend = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


dbca_0<-which(DBCAcl==0)

silhouette(DBCAcl[-dbca_0],dist(data_graph[-dbca_0,]))[,3]%>%mean()
silhouette(gmm_g$classification,dist(data_graph))[,3]%>%mean()
silhouette(gmm0$classification,dist(data_graph))[,3]%>%mean()
silhouette(realgra,dist(data_graph))[,3]%>%mean()
silhouette(gra$depth_clus$cluster_vector,dist(data_graph))[,3]%>%mean()

calinhara(data_graph[-dbca_0,],DBCAcl[-dbca_0])
calinhara(data_graph,gmm_g$classification)
calinhara(data_graph,gmm0$classification)
calinhara(data_graph,realgra)
calinhara(data_graph,gra$depth_clus$cluster_vector)


r_score<-NC_score(dm_gra,realgra,T)
dlcc_score<-NC_score(dm_gra,gra$depth_clus$cluster_vector,T)
gmm_score<-NC_score(dm_gra,gmm_g$classification,T)
gmm0_score<-NC_score(dm_gra,gmm0$classification,T)
dbca_0<-which(DBCAcl==0)
dbcacl_score<-NC_score(dm_gra[-dbca_0,-dbca_0],DBCAcl[-dbca_0],T)


#For cassini


ba_dlcc<-NC_score(dm_ba,ba_temp_clus$depth_clus$cluster_vector,T)
fake_lab<-ba_temp_clus$depth_clus$cluster_vector
fake_lab[501:600]<-2
ba_f2<-NC_score(dm_ba,fake_lab,T)
fake_lab2<-ba_temp_clus$depth_clus$cluster_vector
fake_lab2[c(561)]<-2
ba_f3<-NC_score(dm_ba,fake_lab2,T)


ggplot(data_ba, aes(x = x, 
                    y = y, color = Cluster)) + geom_point(show.legend = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggplot(data_ba, aes(x = x, 
                    y = y, color = fake_lab%>%as.factor())) + geom_point(show.legend = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data_ba, aes(x = x, 
                    y = y, color = fake_lab2%>%as.factor())) + geom_point(show.legend = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

silhouette(ba_temp_clus$depth_clus$cluster_vector,dist(databa))[,3]%>%mean()
silhouette(fake_lab,dist(databa))[,3]%>%mean()
silhouette(fake_lab2,dist(databa))[,3]%>%mean()

calinhara(databa,fake_lab)
calinhara(databa,ba_temp_clus$depth_clus$cluster_vector)
calinhara(databa,fake_lab2)

