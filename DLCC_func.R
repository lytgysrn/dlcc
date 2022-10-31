library(tidyverse)
library(ddalpha)
library(factoextra)
library(mclust)
library(pgmm)
library(ggplot2)
library(mlbench)
library(fossil)
library(cluster)
#ALL DEPTHS HERE ARE MAHALANOBIS DEPTH


#depth_by_cluster: depth with regard to each cluster
depth_by_cluster<-function(X,Kclus,cluster,sub){
  subX<-X[sub,]
  DD<- matrix(0,dim(subX)[1],Kclus) 
  for(i in (1:Kclus)){ 
    DD[,i] <- depth.Mahalanobis(subX,X[cluster[[i]],],mah.estimate = 'MCD')
  } 
  return(DD)
}

#Section: Similarity matrix building
#depth.dm: similarity matrix based on Mahalanobis depth
depth.dm<-function(data,method){
  #method: 0 for the empirical covariance, 1 for MCD covariance
  n=nrow(data)
  dm<-matrix(0,n,n)
  
  if (method == 0){
    covm<-cov(data)
    covm.i<-solve(covm)
    
    for (i in 1:(n-1)) {
      xi<-as.vector(data[i,])
      for (j in (i+1):n) {
        xj<-as.vector(data[j,])
        c<-as.numeric(xj-xi)
        dist<-1/(1+t(c)%*%as.matrix(covm.i)%*%c)
        
        dm[i,j]<-dist
        dm[j,i]<-dist
      }
    }
  } else if (method == 1){
    
    covm<-covMcd(data,nsamp = 'deterministic')$cov
    # covm<-cov.mcd(data)$cov
    covm.i<-solve(covm)
    
    for (i in 1:(n-1)) {
      xi<-as.vector(data[i,])
      for (j in (i+1):n) {
        xj<-as.vector(data[j,])
        c<-as.numeric(xj-xi)
        dist<-1/(1+t(c)%*%as.matrix(covm.i)%*%c)
        
        dm[i,j]<-dist
        dm[j,i]<-dist
      }
    }
  } 
  diag(dm)<-1
  return(dm)
}

#depth.dm: similarity matrix based on Mahalanobis depth with multiple covariance matrices
depth.dm_m<-function(data,class,cov.mat.set){
  n=nrow(data)
  dm<-matrix(0,n,n)
  K=length(cov.mat.set)
  covm.i.set<-lapply(1:K, function(x){solve(cov.mat.set[[x]])})
  
  for (i in 1:(n-1)) {
    xi<-as.vector(data[i,])
    for (j in (i+1):n) {
      xj<-as.vector(data[j,])
      c<-as.numeric(xj-xi)
      dist<-1/(1+t(c)%*%as.matrix(covm.i.set[[class[i]]])%*%c)
      dm[i,j]<-dist
      if (class[i]==class[j]){
        dm[j,i]<-dist
      } else {
        dm[j,i]<-1/(1+t(c)%*%as.matrix(covm.i.set[[class[j]]])%*%c)
      }
    }
  } 
  diag(dm)<-1
  return(dm)
}

#decomp_mat:decomp covariance matrix
decomp_mat<-function(matrix){
  eigen<-eigen(matrix)
  lambda<-abs(det(matrix))^(1/dim(matrix)[1])
  eigenvalue<-eigen$values%>%as.vector()
  ratio<-eigenvalue[-length(eigenvalue)]/eigenvalue[length(eigenvalue)]
  x<-(1/prod(ratio))^(1/dim(matrix)[1])
  eigenvalue_no<-c(ratio,1)*x
  A<-diag(eigenvalue_no)%>%as.matrix()
  D<-eigen$vectors%>%as.matrix()
  return(list(lambda=lambda,D=D,A=A))
}

loglik_gm<-function(data,mu,covm){
  d<-dmvnorm(data,mean=mu,sigma=covm)
  llh<-sum(log(d))
  return(llh)
}
#correlation test based on PCA
if_corr<-function(data){
  pr.data<-princomp(data)
  ve<-pr.data$sdev^2/sum(pr.data$sdev^2)
  if (ve[1]<0.6){
    i<-max(which(ve>0.05))
    if (i==length(ve)){
      i<-i-1
    }
    cumve<-sum(ve[1:i])
    if (cumve<0.95){
      d<-FALSE
    } else {
      d<-TRUE
    }
  } else {
    d<-TRUE
  }
  return(d)
}
#EM_CM: EM algorithm for AEEV GMM
EM_CM<-function(data,Kset,t,covm){
  #data: data of observations
  #K: # of components
  #t: maximum number of iterations
  
  X=data
  D=dim(X)[2]
  N=dim(X)[1]
  cov_de<-decomp_mat(covm)
  lambda<-cov_de$lambda
  A<-cov_de$A
  
  Qhis<-c()
  muhis<-list()
  covhis<-list()
  alphahis<-list()
  zlist<-list()
  for (l in 1:length(Kset)) {
    #Initial para
    K=Kset[l]
    alpha=rep(1/K,K)
    #use kmeans initial mu and cov
    X.kmeans=kmeans(X,K)
    X.kmeans.df<-data.frame(x=X,cluster=X.kmeans$cluster)
    
    mu=t(sapply(1:K, function(x){colMeans(X.kmeans.df[which(X.kmeans.df$cluster==x),]
                                          [,-(D+1)])}))
    cov=lapply(1:K,function(x){covm})
    
    z=matrix(rep(0,N*K),N,K)
    
    eps=1
    Q=0
    Q[1]=0
    for (h in 2:(t+1)) {
      # E step
      p=matrix(rep(0,N*K),N,K)
      for (i in 1:K) {
        p[,i]=alpha[i]*dmvnorm(X,mean=mu[i,],sigma=cov[[i]])
      }
      p_sum=rowSums(p)
      #update z
      z=p/p_sum
      #update Q and eps
      Q[h]=sum(z*log(p))
      eps=abs(Q[h]-Q[h-1])
      
      #determine stop or not
      if (eps>1e-3){
        #M step
        sumz=colSums(z)
        #update pi
        alpha=sumz/N
        
        for (j in 1:K) {
          zx=z[,j]*X
          #update mu
          mu[j,]=colSums(zx)/sumz[j]
          
          x_min_mu=sweep(X,2, mu[j,])
          z_x_min_mu=z[,j]*x_min_mu
          #update covariance matrix, get D*D diagonal matrix
          new_D=eigen(t(as.matrix(z_x_min_mu))%*%as.matrix(x_min_mu))$vectors
          cov[[j]]<-lambda*new_D%*%A%*%t(new_D)
        }
      } else {
        break
      }
    }
    Qhis[l]=Q[length(Q)]
    muhis[[l]]<-mu
    covhis[[Kset[l]%>%as.character()]]<-cov
    alphahis[[l]]<-alpha
    zlist[[l]]<-z
  }
  npara<-Kset-1+Kset*D+Kset*D*(D-1)/2
  BIC<-2*Qhis-npara*log(N)
  bml<-which.max(BIC)
  la<-Kset[bml]%>%as.character()
  best_model<-list('K'=Kset[bml],'mu'=muhis[[bml]],'cov'=covhis[[bml]],z=zlist[[bml]],loglik=Qhis[bml],BIC=BIC[bml])
  return(best_model)
}

#Section: Local center selections
#tan value
calc.k<-function(k1,k2){
  k=(k1-k2)/(1+k1*k2)
  return(k)
}  

#current proportion of a's neighbors
f.prop<-function(a,size,Nobs,dm0.order){
  all_nbs<-dm0.order[1:size,a]%>%as.numeric()%>%unique
  prop<-length(all_nbs)/Nobs
}

#Section: Clustering in DLCC
#Assign_score: scores of obs in get.temp.cluster under min strategy
Assign_score<-function(Kclus,dm0,temp.clus,temp.cl){
  a<-unlist(temp.cl)
  deflist<-list()
  for (k in 1:Kclus){
    if (length(temp.clus[[k]])!=0){
      max_a<-sapply(1:length(temp.clus[[k]]),function(x){max(dm0[temp.cl[[k]],temp.clus[[k]][x]])})
      max_b<-sapply(1:length(temp.clus[[k]]),function(x){max(dm0[setdiff(a,temp.cl[[k]]),temp.clus[[k]][x]])})
      larger<-sapply(1:length(temp.clus[[k]]),function(x){ifelse((max_a[x]-max_b[x])>0,max_a[x],max_b[x])})
      deflist[[k]]<-(max_a-max_b)/larger
    } else {
      deflist[[k]]<-NA
    }
  }
  return(deflist)
}

#get.temp.cl: obtain group results of local centers
get.temp.cl<-function(a,dm0,size,Th,method){
  
  Nobs<-nrow(dm0)
  temp.cl<-list()
  temp.save<-list()
  center_lengh<-length(a)
  dm0.order<-sapply(1:Nobs,function(i){sort.list(dm0[i,],decreasing = T)})
  
  
  if (method=='min'){
    sim.mat1<-sim.mat(N=center_lengh,dm0 = dm0,size = size,centers = a)
    for (i in 1:center_lengh) {
      temp.save[[i]]<-which(sim.mat1[,i]>Th)
    }
    temp.cl=list()
    for (j in 1:center_lengh) {
      temp.cl[[j]]<-Reduce(intersect,temp.save[temp.save[[j]]])
    }
    dupli<-duplicated(temp.cl)
    Tdup<-which(dupli==TRUE)
    if (length(Tdup)!=0){
      temp.cl<-temp.cl[-Tdup]
    }
    l.ts<-length(temp.cl)
    rel.mat<-matrix(0,l.ts,l.ts)
    for (i in 1:(l.ts)) {
      for (j in 1:l.ts) {
        rel<-length(intersect(temp.cl[[i]],temp.cl[[j]]))/length(temp.cl[[i]])
        if(rel!=1){
          rel=0
        }
        rel.mat[i,j]<-rel
      }
    }
    
    dr<-which(rowSums(rel.mat)>1)
    if (length(dr)!=0){
      temp.cl<-temp.cl[-dr]
    }
    
    temp.cl.save<-temp.cl
    for (i in 1:length(temp.cl)) {
      temp.cl[[i]]<-a[temp.cl[[i]]]
    }
    
    
  } else if (method=='max'){
    nclus<-1
    temp.cl[[nclus]]<-a[1]
    for (i in 2:center_lengh) {
      save.length<-length(unlist(temp.cl))
      for (j in 1:length(temp.cl)) {
        sim<-sapply(1:length(temp.cl[[j]]),function(x){
          length(intersect(dm0.order[1:size,temp.cl[[j]][x]],dm0.order[1:size,a[i]]))/size})
        
        if (max(sim)>Th){
          temp.cl[[j]]<-c(temp.cl[[j]],a[i])
        }
      }
      #If the length of obs in temp.cl is not change, which means a[i] is not accepted by current group, create a new group for a[i]
      new.length<-length(unlist(temp.cl))
      if (new.length==save.length){
        nclus=nclus+1
        temp.cl[[nclus]]<-a[i]
        save.length<-new.length+1
      } 
    }
  }
  
  rep.label<-as.numeric(names(which(table(unlist(temp.cl))>1)))
  if (length(rep.label)>0){
    for (i in 1:length(rep.label)) {
      cllabel<-c()
      if (method=='min'){
        for (j in 1:length(temp.cl)) {
          if (rep.label[i]%in%temp.cl[[j]]){
            cllabel<-c(cllabel,j)
            #saves the IDs of temp.cl contains the repeated center
          }
        }
        
        for (m in cllabel){
          temp.cl[[m]]<-setdiff(temp.cl[[m]],rep.label[i])
        }
        
      } else if (method=='max'){
        for (j in 1:length(temp.cl)) {
          if (rep.label[i]%in%temp.cl[[j]]){
            cllabel<-c(cllabel,j)
          }
        }
        remain.cluster<-cllabel[1]
        obs<-unique(unlist(lapply(1:length(cllabel), function(x){temp.cl[[cllabel[x]]]})))
        temp.cl[[remain.cluster]]<-obs
        cllabel<-cllabel[-1]
        for (w in rev(cllabel)) {
          temp.cl[[w]]<-NULL
        }
      }
    }
  }
  return(temp.cl)
}

#remove.ol.obs: remove overlapped nbr between groups
remove.ol.obs<-function(a,clus.ind,table.list){
  temp.clus<-lapply(clus.ind, function(x){ as.numeric(names(table.list[[x]]))})
  t.obs<-table(unlist(temp.clus))
  #Obs appear more than once in clusters 
  ulabel<-(names(t.obs[which(t.obs>1)]))%>%as.numeric()
  
  if (length(ulabel)!=0){
    for (i in clus.ind) {
      temp.clus[[i]]<-setdiff(temp.clus[[i]],ulabel)
    }
  }
  return(list(temp.clus=temp.clus,ulabel=ulabel))
}

#get.temp.cluster: obtain temporary clusters
get.temp.cluster<-function(dm0,table.list,Kclus,temp.cl,size,method){
  a<-unlist(temp.cl)
  clus.ind<-1:Kclus
  cl<-sapply(1:length(temp.cl), function(x){length(temp.cl[[x]])})
  cl<-rep(clus.ind,cl)
  
  temp<-remove.ol.obs(a,clus.ind,table.list)
  temp.clus<-temp$temp.clus
  ulabel<-temp$ulabel
  
  if (method=='max'){
    
    if (length(ulabel)!=0){
      #Remove those obs in temp cluster
      
      #if there is a cluster without any obs
      delete_label<-c()
      for (x in clus.ind) {
        if (is_empty(temp.clus[[x]])){
          delete_label<-c(delete_label,x)
          cl<-cl[-which(cl==x)]%>%as.factor()%>%as.numeric()
          Kclus<-Kclus-1
          clus.ind<-clus.ind[-length(clus.ind)]
          a<-a[-x]
        }
      }
      if (is_empty(delete_label)==FALSE){
        table.list<-list()
        temp.cl<-temp.cl[-delete_label]
        for (i in 1:Kclus) {
          table.list[[i]]<-table(dm0.order[1:size,temp.cl[[i]]])
        }
        temp.clus<-remove.ol.obs(a,1:Kclus,table.list)$temp.clus
      }
      
      ulabel<-ulabel%>%as.character()
      score.mat<-matrix(0,length(ulabel),length(table.list))
      for (i in 1:length(ulabel)) {
        for (j in 1:length(table.list)){
          if (is.na(table.list[[j]][ulabel[i]])){
            score.mat[i,j]<-0
          } else {
            score.mat[i,j]<-table.list[[j]][ulabel[i]]
          }
        }
      }
      clus.assign<-c()
      for (i in 1:nrow(score.mat)) {
        h<-which(score.mat[i,]==max(score.mat[i,]))
        if (length(h)==1){
          clus.assign[i]<-h
        } else {
          #Draw case, not assign
          clus.assign[i]<-0
        }
      }
      temp.clb<-lapply(clus.ind,function(x){as.numeric(ulabel[which(clus.assign==x)])})
      for (i in clus.ind) {
        temp.clus[[i]]<-c(temp.clus[[i]],temp.clb[[i]])
      }
    }
  }
  
  
  if (method=='min'){
    
    
    ascore<-Assign_score(Kclus,dm0,temp.clus,temp.cl)
    
    
    tobs<-unlist(temp.clus)
    otherobs<-setdiff(1:nrow(dm0),tobs)
    
    stcall<-sapply(1:length(otherobs),function(x){cl[which.max(dm0[a,otherobs[x]])]})
    #assign cluster for left obs.
    left_clus<-list()
    for (i in clus.ind){
      left_clus[[i]]<-otherobs[which(stcall==i)]
    }
    score_left<-Assign_score(Kclus,dm0,left_clus,temp.cl)
    
    #set med score_left as the threshold for obs in current clusters
    lq<-sapply(1:Kclus, function(x){ifelse((length(score_left[[x]])<2),NA,quantile(score_left[[x]],0.5))})
    for (i in 1:Kclus){
      index<-which(ascore[[i]]<lq[i])
      if (length(index)>0){
        left_clus[[i]]<-c(left_clus[[i]],temp.clus[[i]][index])
        score_left[[i]]<-c(score_left[[i]],ascore[[i]][index])
        temp.clus[[i]]<-temp.clus[[i]][-index]
        ascore[[i]]<-ascore[[i]][-index]
      }
    }
    
    
    #if there is a cluster without any obs
    delete_label<-c()
    for (x in clus.ind) {
      if (is_empty(temp.clus[[x]])){
        delete_label<-c(delete_label,x)
        cl<-cl[-which(cl==x)]%>%as.factor()%>%as.numeric()
        Kclus<-Kclus-1
        clus.ind<-clus.ind[-length(clus.ind)]
        a<-a[-x]
      }
    }
    if (is_empty(delete_label)==FALSE){
      table.list<-list()
      temp.cl<-temp.cl[-delete_label]
      for (i in 1:Kclus) {
        table.list[[i]]<-table(dm0.order[1:size,temp.cl[[i]]])
      }
      temp.clus<-remove.ol.obs(a,1:Kclus,table.list)$temp.clus
      ascore<-Assign_score(Kclus,dm0,temp.clus,temp.cl)
    }
    
    
    new.border<-c()
    for (i in clus.ind) {
      num_obs_current<-length(temp.clus[[i]])
      
      sortscore<-sort(ascore[[i]],decreasing = T)[floor(num_obs_current/2):num_obs_current]
      cdd_1<-sortscore[which.min(diff(sortscore))]
      if (num_obs_current<(size/2)){
        cdd_2<-quantile(score_left[[i]],1-(size/2-num_obs_current)/length(score_left[[i]]))
      } else {
        cdd_2=1
      }
      new.border[i]<-min(cdd_1,cdd_2)
    }
    
    for (x in clus.ind) {
      index<-which(score_left[[x]]>new.border[x])
      left_clus[[x]]<-left_clus[[x]][index]
    }
    
    for (i in clus.ind) {
      temp.clus[[i]]<-c(temp.clus[[i]],left_clus[[i]])
    }
    
  }
  return(temp.clus)
}
#get.local.center: obtain filtered local centers
get.local.center<-function(data,dm0,size,method){
  Nobs<-nrow(data)
  dm0.order<-sapply(1:Nobs,function(i){sort.list(dm0[i,],decreasing = T)})
  ls<-length(size)
  
  rank.mat<-matrix(0,Nobs,ls)
  
  ranktopobs<-c()
  for (i in size) {
    for(j in 1:Nobs){
      local.d<-depth.Mahalanobis(data[dm0.order[1:i,j],],data[dm0.order[1:i,j],],mah.estimate	='MCD')
      rank.mat[j,as.numeric(which(size==i))]<-rank(-local.d,ties.method = 'min')[1]
      ranktopobs<-c(ranktopobs,which(rank(-local.d,ties.method = 'min')==1)%>%names()%>%as.numeric())
      
    }
  }
  #local centers
  t.rank<-table(ranktopobs)
  
  
  
  if (method=='min'){
    a_2<-t.rank%>%names()%>%as.numeric()
    a_2<-a_2[rank.mat[a_2]<=2]
    t.choose<-t.rank[a_2%>%as.character()]
    
    sort.choose<-sort(t.choose,decreasing = T)
    n.sort.choose<-sort.choose%>%names()%>%as.numeric()
    
    dupli.lab<-which(duplicated(sort.choose)==T)
    if (length(dupli.lab)!=0){
      for (i in 1:length(dupli.lab)) {
        prop1<-f.prop(n.sort.choose[1:(dupli.lab[i]-1)],size,Nobs,dm0.order)
        prop2<-f.prop(n.sort.choose[-(dupli.lab[i]-1)][1:(dupli.lab[i]-1)],size,Nobs,dm0.order)
        if (prop1<prop2){
          tmp<-n.sort.choose[dupli.lab[i]]
          n.sort.choose[dupli.lab[i]]<-n.sort.choose[dupli.lab[i]-1]
          n.sort.choose[dupli.lab[i]-1]<-tmp
        }
      }
      names(sort.choose)<-n.sort.choose
    }
    
    props<-sapply(1:length(sort.choose),function(x){f.prop(n.sort.choose[1:x],size,Nobs,dm0.order)})
    d.props<-diff(props)
    k<-sapply(1:(length(d.props)-1), function(x){calc.k(d.props[x],d.props[x+1])})
    nk<-k[which(k<0)]
    nk<-nk[which(-nk>median(abs(k)))]
    k<-abs(k)
    if (length(nk)==0){
      sort.k<- sort(k,decreasing = T)
      prop.k<--diff(sort.k)
    } else {
      sort.k<-abs(sort(nk))
      prop.k<--diff(nk)
    }
    if (length(prop.k)<=1){p.id=1}else{p.id<-order(prop.k,decreasing = T)}
    lab_pos<-which(k==sort.k[p.id[1]])[1]
    if (d.props[lab_pos]>d.props[lab_pos+1]){
      a<-n.sort.choose[1:(lab_pos+1)]
    } else {
      a<-n.sort.choose[1:(lab_pos+2)]
    }
    while (props[length(a)]<0.7){
      p.id<-p.id[-1]
      if (length(p.id)==0){
        ind<-min(which(props>=0.7)[1],length(props))
        a<-n.sort.choose[1:ind]
        props[length(a)]=0.7
      } else {
        lab_pos<-which(k==sort.k[p.id[1]])[1]
        if (d.props[lab_pos]>d.props[lab_pos+1]){
          a<-n.sort.choose[1:(lab_pos+1)]
        } else {
          a<-n.sort.choose[1:(lab_pos+2)]
        }
      }
    }
  } else if (method=='max'){
    a<-t.rank%>%names()%>%as.numeric()
    a<-a[rank.mat[a]<=2]
    a_sort<-sort(t.rank[a%>%as.character()],decreasing = T)
    a<-which(a_sort>1)%>%names()%>%as.numeric()
    all_nbs<-dm0.order[1:size,a]%>%as.numeric()%>%unique
    prop<-length(all_nbs)/Nobs
    if (prop<0.85){
      a_2<-t.rank%>%names()%>%as.numeric()%>%setdiff(a)
      nbs_a_2<-lapply(1:length(a_2),function(x){dm0.order[1:size,a_2[x]]})
      l_a_2<-c()
      for (g in 1:length(a_2)) {
        l_a_2[g]<-length(intersect(nbs_a_2[[g]],all_nbs))
      }
      a_2<-a_2[which(l_a_2==0)]
      a<-union(a,a_2)
    }
  }
  a<-sort(a)
  return(a)
}

#DRcluster: The whole DLCC algorithm, include other functions.
DRcluster<-function(data,dm0,Th,size,method,class_method,maxdepth){
  #data:data
  #dm0:similarity matrix
  #Th:threshold
  #size:size
  #method:min/max
  #class_method:classification method
  #maxdepth:T/F
  Nobs<-nrow(data)
  dm0.order<-sapply(1:Nobs,function(i){sort.list(dm0[i,],decreasing = T)})
  ls<-length(size)
  
  rank.mat<-matrix(0,Nobs,ls)
  # dep.mat<-matrix(0,nrow(data),ls)
  # Ind.mat<-matrix(0,Nobs,ls) 
  ranktopobs<-c()
  for (i in size) {
    for(j in 1:Nobs){
      local.d<-depth.Mahalanobis(data[dm0.order[1:i,j],],data[dm0.order[1:i,j],],mah.estimate	='MCD')
      rank.mat[j,as.numeric(which(size==i))]<-rank(-local.d,ties.method = 'min')[1]
      #Each elements records the centrality rank of obs j within the cluster constructed by it and its neighbors
      ranktopobs<-c(ranktopobs,which(rank(-local.d,ties.method = 'min')==1)%>%names()%>%as.numeric())
      # dep.mat[j,as.numeric(which(size==i))]<-local.d[1]
      # Ind.mat[j,as.numeric(which(size==i))]<-signif(mean(local.d), 10)
      #Each element records the mean depth of the cluster constructed by obs j and its neighbors
    }
  }
  
  #local centers
  t.rank<-table(ranktopobs)
  
  
  
  if (method=='min'){
    a_2<-t.rank%>%names()%>%as.numeric()
    a_2<-a_2[rank.mat[a_2]<=2]
    t.choose<-t.rank[a_2%>%as.character()]
    
    sort.choose<-sort(t.choose,decreasing = T)
    n.sort.choose<-sort.choose%>%names()%>%as.numeric()
    
    dupli.lab<-which(duplicated(sort.choose)==T)
    if (length(dupli.lab)!=0){
      for (i in 1:length(dupli.lab)) {
        prop1<-f.prop(n.sort.choose[1:(dupli.lab[i]-1)],size,Nobs,dm0.order)
        prop2<-f.prop(n.sort.choose[-(dupli.lab[i]-1)][1:(dupli.lab[i]-1)],size,Nobs,dm0.order)
        if (prop1<prop2){
          tmp<-n.sort.choose[dupli.lab[i]]
          n.sort.choose[dupli.lab[i]]<-n.sort.choose[dupli.lab[i]-1]
          n.sort.choose[dupli.lab[i]-1]<-tmp
        }
      }
      names(sort.choose)<-n.sort.choose
    }
    
    props<-sapply(1:length(sort.choose),function(x){f.prop(n.sort.choose[1:x],size,Nobs,dm0.order)})
    d.props<-diff(props)
    k<-sapply(1:(length(d.props)-1), function(x){calc.k(d.props[x],d.props[x+1])})
    nk<-k[which(k<0)]
    nk<-nk[which(-nk>median(abs(k)))]
    k<-abs(k)
    if (length(nk)==0){
      sort.k<- sort(k,decreasing = T)
      prop.k<--diff(sort.k)
    } else {
      sort.k<-abs(sort(nk))
      prop.k<--diff(nk)
    }
    if (length(prop.k)<=1){p.id=1}else{p.id<-order(prop.k,decreasing = T)}
    lab_pos<-which(k==sort.k[p.id[1]])[1]
    if (d.props[lab_pos]>d.props[lab_pos+1]){
      a<-n.sort.choose[1:(lab_pos+1)]
    } else {
      a<-n.sort.choose[1:(lab_pos+2)]
    }
    while (props[length(a)]<0.7){
      p.id<-p.id[-1]
      if (length(p.id)==0){
        ind<-min(which(props>=0.7)[1],length(props))
        a<-n.sort.choose[1:ind]
        props[length(a)]=0.7
      } else {
        lab_pos<-which(k==sort.k[p.id[1]])[1]
        if (d.props[lab_pos]>d.props[lab_pos+1]){
          a<-n.sort.choose[1:(lab_pos+1)]
        } else {
          a<-n.sort.choose[1:(lab_pos+2)]
        }
      }
    }
    
  } else if (method=='max'){
    a<-t.rank%>%names()%>%as.numeric()
    a<-a[rank.mat[a]<=2]
    a_sort<-sort(t.rank[a%>%as.character()],decreasing = T)
    a<-which(a_sort>1)%>%names()%>%as.numeric()
    all_nbs<-dm0.order[1:size,a]%>%as.numeric()%>%unique
    prop<-length(all_nbs)/Nobs
    if (prop<0.85){
      a_2<-t.rank%>%names()%>%as.numeric()%>%setdiff(a)
      nbs_a_2<-lapply(1:length(a_2),function(x){dm0.order[1:size,a_2[x]]})
      l_a_2<-c()
      for (g in 1:length(a_2)) {
        l_a_2[g]<-length(intersect(nbs_a_2[[g]],all_nbs))
      }
      a_2<-a_2[which(l_a_2==0)]
      a<-union(a,a_2)
    }
  }
  a<-sort(a)
  
  temp.cl<-get.temp.cl(a,dm0,size,Th,method)
  if (length(temp.cl)==1){
    temp.clus=list()
    cluster=list()
    cluster_vector=list()
    cluster[[1]]<-1:Nobs
    cluster_vector[[1]]<-rep(1,Nobs)
    temp.clus[[1]]<-1:Nobs
    depth_clus<-list(cluster=cluster,cluster_vector=cluster_vector)
  } else {
    if (method=='min'){
      stop=0
      a.save<-a
      loop.i<-0
      while(stop!=1){
        Kclus=length(temp.cl)
        table.list<-list()
        for (i in 1:Kclus) {
          table.list[[i]]<-table(dm0.order[1:size,temp.cl[[i]]])
        }
        temp.clus<-get.temp.cluster(dm0,table.list = table.list,Kclus = Kclus, temp.cl,size, method = method)
        depth_clus<-DAobs(data,temp.clus,method = class_method ,maxdepth)
        Kclus<-length(temp.clus)
        d<-depth_by_cluster(data,Kclus = Kclus,depth_clus$cluster,sub = a.save)
        
        
        new.a<-sapply(1:Kclus, function(x){a.save[which.max(d[,x])]})
        if (TRUE%in%duplicated(new.a)){
          new.a[duplicated(new.a)]<-setdiff(a.save,new.a)[1]
        }
        
        if (setequal(new.a,a)!=TRUE){
          temp.cl<-lapply(1:Kclus,function(x){new.a[x]})
          a<-new.a
          loop.i<-loop.i+1
          if (loop.i>10){
            stop<-1
          }
        } else {
          stop<-1
        }
      }
    } else {
      Kclus=length(temp.cl)
      table.list<-list()
      for (i in 1:Kclus) {
        table.list[[i]]<-table(dm0.order[1:size,temp.cl[[i]]])/length(temp.cl[[i]])
      }
      temp.clus<-get.temp.cluster(dm0,table.list = table.list,Kclus = Kclus, temp.cl,size, method = method)
      depth_clus<-DAobs(data,temp.clus,method = class_method ,maxdepth)
    }
  }
  return(list('temp.center'=temp.cl,'temp.clus'=temp.clus,'depth_clus'=depth_clus))
}
#DAobs: Given temp clusters, given final clustering results
DAobs<-function(data,temp.clus,method,maxdepth){
  X<-data
  Kclus<-length(temp.clus)
  labelledobs<-unlist(temp.clus)
  left.obs<-X[-labelledobs,]
  left.obs.label<-setdiff(1:nrow(X),labelledobs)
  
  if (method=='maxdep'){
    depth.mat<-matrix(rep(0),nrow(left.obs),Kclus)
    for (i in 1:nrow(left.obs)){
      for (j in 1:Kclus){
        depth.mat[i,j]<-depth.Mahalanobis(left.obs[i,],X[temp.clus[[j]],])
      }
    }
    depth.order<-sapply(1:nrow(depth.mat),function(i){sort.list(depth.mat[i,],decreasing = T)})
    
    cluster=list()
    for (i in 1:Kclus) {
      cluster[[i]]=left.obs.label[which(depth.order[1,]==i)]
      cluster[[i]]=c(temp.clus[[i]],cluster[[i]])
    }
  } else if (method=='mclass'){
    current.label<-rep(0,nrow(X))
    for (i in 1:length(temp.clus)) {
      current.label[temp.clus[[i]]]<-i
    }
    result<-MclustDA(X[labelledobs,],class = current.label[labelledobs],modelType = "MclustDA")
    pre<-predict.MclustDA(result, newdata = X[-labelledobs,])
    cluster=list()
    for (i in 1:Kclus) {
      cluster[[i]]=left.obs.label[which(pre$classification==i)]
      cluster[[i]]=c(temp.clus[[i]],cluster[[i]])
    }
  }
  # else if (method=='ddalpha'){
  #   current.label<-rep(0,nrow(X))
  #   for (i in 1:length(temp.clus)) {
  #     current.label[temp.clus[[i]]]<-i
  #   }
  #   data.X<-cbind(X,current.label)%>%as.data.frame()
  # 
  #   dd<-ddalpha.train(current.label~.,data.X,subset = labelledobs,depth = 'zonoid',
  #                     separator = 'alpha',outsider.methods = 'kNNAff')
  #   classes<-ddalpha.classify(dd,data.X,subset = left.obs.label)
  #   ddresult<-unlist(classes)
  # 
  #   cluster=list()
  #   for (i in 1:Kclus) {
  #     cluster[[i]]=left.obs.label[which(ddresult==i)]
  #     cluster[[i]]=c(temp.clus[[i]],cluster[[i]])
  #   }
  # }
  
  cc<-rep(0,nrow(X))
  for (i in 1:length(cluster)) {
    cc[cluster[[i]]]<-i
  }
  
  if (maxdepth==TRUE){
    nelabel=0
    
    while (length(nelabel)!=0) {
      DD<- matrix(0,dim(X)[1],Kclus) 
      for(i in (1:Kclus)){ 
        DD[,i] <- depth.Mahalanobis(X,X[cluster[[i]],])
      } 
      DDmat<-sapply(1:nrow(DD),function(i){sort.list(DD[i,],decreasing = T)})
      nelabel<-which(DDmat[1,]!=cc)
      cc<-DDmat[1,]
      if (length(nelabel)!=0){
        for (i in 1:length(cluster)) {
          cluster[[i]]<-which(cc==i)
        }
      }
    }
  }
  return(list(cluster=cluster,cluster_vector=cc))
}

#sim.mat: similarity matrix of centers
sim.mat<-function(N,dm0,size,centers){
  dm0.order<-sapply(1:nrow(dm0),function(i){sort.list(dm0[i,],decreasing = T)})
  sim.mat<-matrix(0,N,N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      sim<-length(intersect(dm0.order[1:size,centers[i]],dm0.order[1:size,centers[j]]))/size
      sim.mat[i,j]<-sim
      sim.mat[j,i]<-sim
    }
    diag(sim.mat)<-1
  }
  return(sim.mat)
}

#Section: Parameter selection (max strategy). Also see DCscore.R, and read the Density based clustering score function.

#is.diag: determine if except the diag, elements of matrix in para_select_max are all zero
is.diag<-function(matrix){
  diag(matrix)<-0
  result<-sum(matrix)==0
  return(result)
}

#para_select_max: under max strategy, with given sizes guess proper thresholds 
para_select_max<-function(data,sizelist,dm0,probs){
  Nobs=nrow(data)
  paralist=list()
  Thsetlist<-list()
  d=dim(data)[2]
  
  for (q in 1:length(sizelist)) {
    centers<-get.local.center(data,dm0,sizelist[q],method = 'max')
    N=length(centers)
    
    sim.mat1<-sim.mat(N=N,dm0 = dm0,size = sizelist[q],centers = centers)
    sim.order<-sapply(1:nrow(sim.mat1),function(i){sort(sim.mat1[i,],decreasing = T)})
    # sapply(1:N,function(i){min(sim.order[,i][sim.order[,i]!=0])})
    
    
    num_nbr<-min(2*d,round(N/2))
    nbrs<-sapply(1:N, function(x){mean(sim.order[2:(num_nbr+1),x])})
    
    bc<-quantile(nbrs, probs = probs) 
    vv<-which(nbrs<bc)
    
    maxsim<-sim.order[2,][vv]
    upper<-sim.order[2,]%>%min
    
    #find the lower boundary
    tcl<-get.temp.cl(centers,dm0,sizelist[q],upper*0.99,method='max')
    temp.K<-length(tcl)
    gap_bet_clus<-matrix(0,temp.K,temp.K)
    for (y in 1:(temp.K-1)) {
      for (u in (y+1):temp.K) {
        c1<-match(tcl[[y]],centers)
        c2<-match(tcl[[u]],centers)
        max_sim<-max(sim.mat1[c2,c1])
        gap_bet_clus[u,y]<-max_sim
        gap_bet_clus[y,u]<-max_sim
      }
    }
    diag(gap_bet_clus)<-1
    
    Thset<-gap_bet_clus%>%as.numeric()%>%sort(decreasing = T)%>%unique()
    Thset<-Thset[-1]
    Thset<-Thset[Thset!=0]
    m=1
    lower=c()
    
    if (length(Thset)==1){
      Th<-Thset
      cut=1
    } else if (length(Thset)==0) {
      Th<-upper*0.99
      cut=0
    }else {
      
      while(is_empty(lower)==T){
        v=1
        while (v <= temp.K) {
          if (Thset[m]%in%gap_bet_clus[,v]){
            comb_lab<-which(gap_bet_clus[,v]==Thset[m])
            # gap_bet_clus[comb_lab,v]<-0
            
            gap_bet_clus[1:temp.K,v]<-sapply(1:temp.K,function(x){max(gap_bet_clus[x,c(v,comb_lab)])})
            gap_bet_clus[v,1:temp.K]<-gap_bet_clus[1:temp.K,v]
            gap_bet_clus<-gap_bet_clus[-comb_lab,-comb_lab]
            temp.K<-temp.K-length(comb_lab)
            if (is.matrix(gap_bet_clus)==FALSE){
              gap_bet_clus<-gap_bet_clus%>%as.matrix()
            }
          } else {
            v=v+1
          }
        }
        if (is.diag(gap_bet_clus)==TRUE){
          lower<-Thset[m]
        } else {
          m=m+1
        }
      }
      if (nrow(gap_bet_clus)>1){
        Th<-lower
        cut<-1
      } else {
        sim.mat2<-sim.mat1[vv,vv]
        diag(sim.mat2)<--1
        nrlabel<-which(colSums(sim.mat2)==-1)
        
        if (length(nrlabel)!=0){
          vv.new<-vv[-nrlabel]
          sim.mat2<-sim.mat2[-nrlabel,-nrlabel]
          maxsim<-maxsim[-nrlabel]
        } else {
          vv.new<-vv
        }
        
        if (length(vv.new)>0){
          Th=0
          #sim between edge centers
          maxsim2<-sapply(1:length(vv.new),function(x){max(sim.mat2[,x])})
          #the closest other edge centers label
          labelmaxsim<-sapply(1:length(vv.new),function(x){which.max(sim.mat2[,x])})
          gap<-sapply(1:length(vv.new),function(x){abs(maxsim2[x]-maxsim[labelmaxsim[x]])})
          select_label<-which(gap==max(gap))
          Th<-maxsim2[select_label]%>%min()
          cut<-1
          if(Th%in%Thset==F){
            Th<-0
          } 
          while (Th<lower) {
            if(max(gap)!=0){
              gap[select_label]<-0
              select_label<-which(gap==max(gap))
              Th<-maxsim2[select_label]%>%min()
              if(Th%in%Thset==F){
                Th<-0
              } 
            } else {
              Th<-1
              cut<-0
            }
          }
        } else {
          Th<-0
          cut<-0
        }
        if(Th==1){
          Th<-0
        }
        
        if (Th>=upper){
          Th<-Thset[1]
          cut<-0
        }
      }
    }
    # nbrsim<-sort(nbrs,decreasing = T)
    # Th<-nbrsim[which.max(-diff(nbrsim))]
    paralist[[q]]<-c(sizelist[q],Th,cut)
    Thsetlist[[q]]<-Thset
  }
  return(list(paralist=paralist,Thset=Thsetlist))
}
#para_select_max.K: with known k, give threshold wrt size 
para_select_max.K<-function(data,sizelist,dm0,beginTh,K,outlist){
  Nobs=nrow(data)
  Th_output<-c()
  d=dim(data)[2]
  cutpoint.list<-list()
  c.Klist<-list()
  for (q in 1:length(sizelist)) {
    centers<-get.local.center(data,dm0,sizelist[q],method = 'max')
    N=length(centers)
    sim.mat1<-sim.mat(N=N,dm0 = dm0,size = sizelist[q],centers = centers)
    tcl<-get.temp.cl(centers,dm0,sizelist[q],beginTh,method='max')
    temp.K<-length(tcl)
    gap_bet_clus<-matrix(0,temp.K,temp.K)
    for (y in 1:(temp.K-1)) {
      for (u in (y+1):temp.K) {
        c1<-match(tcl[[y]],centers)
        c2<-match(tcl[[u]],centers)
        max_sim<-max(sim.mat1[c2,c1])
        gap_bet_clus[u,y]<-max_sim
        gap_bet_clus[y,u]<-max_sim
      }
    }
    diag(gap_bet_clus)<-1
    Thset<-gap_bet_clus%>%as.numeric()%>%sort(decreasing = T)%>%unique()
    Thset<-Thset[-1]
    Thset<-Thset[Thset!=0]
    m=1
    stop=0
    Kset=temp.K
    cutpoint<-c()
    while(stop!=1){
      v=1
      while (v <= temp.K) {
        if (Thset[m]%in%gap_bet_clus[,v]){
          comb_lab<-which(gap_bet_clus[,v]==Thset[m])
          # gap_bet_clus[comb_lab,v]<-0
          
          gap_bet_clus[1:temp.K,v]<-sapply(1:temp.K,function(x){max(gap_bet_clus[x,c(v,comb_lab)])})
          gap_bet_clus[v,1:temp.K]<-gap_bet_clus[1:temp.K,v]
          gap_bet_clus<-gap_bet_clus[-comb_lab,-comb_lab]
          temp.K<-temp.K-length(comb_lab)
          if (is.matrix(gap_bet_clus)==FALSE){
            gap_bet_clus<-gap_bet_clus%>%as.matrix()
          }
        } else {
          v=v+1
        }
      }
      Kset<-c(Kset,nrow(gap_bet_clus))
      cutpoint<-c(cutpoint,Thset[m])
      if (is.diag(gap_bet_clus)==TRUE){
        stop<-1
      } else {
        m=m+1
      }
    }
    cutpoint<-c(cutpoint,0)
    Klabel<-which(Kset==K)
    if (length(Klabel)==1){
      Th_output[q]<-cutpoint[Klabel[1]]
    } else {
      Th_output[q]<-NA
    }
    if (outlist==T){
      cutpoint.list[[q]]<-cutpoint
      c.Klist[[q]]<-Kset
    }
  }
  if (outlist==F){
    return(Th_output)
  } else {
    return(list(Th=Th_output,cutpoint=cutpoint.list,K=c.Klist))
  }
}
#para_select_max2: under max strategy, find best combinations of parameters
para_select_max2<-function(data,paralist,dm0,method,maxdepth){
  fs<-c()
  npara<-length(paralist)
  d=dim(data)[2]
  for (j in 1:npara) {
    if (paralist[[j]][3]==1){
      zlist=c(1,2)
    } else {
      zlist=c(1)
    }
    avgss<-c()
    Ksave<-c()
    for (z in zlist) {
      if (z==1){
        temp_clus<-DRcluster(data,dm0,paralist[[j]][2] ,paralist[[j]][1],method = 'max',class_method=method,maxdepth = maxdepth)
      } else {
        temp_clus<-DRcluster(data,dm0,paralist[[j]][2]*0.99 ,paralist[[j]][1],method = 'max',class_method=method,maxdepth = maxdepth)
      }
      Kest<-length(temp_clus$depth_clus$cluster)
      if (Kest==1){
        avgss[z]=-1
      } else {
        Ksave[z]<-Kest
        rds<-c()
        ss<-c()
        NCscore<-NC_score(dm0,temp_clus$depth_clus$cluster_vector,T)
        avgss[z]<-NCscore$score
      }
    }
    large_label<-which.max(avgss)
    fs[j]<-avgss[large_label]
    
    if (large_label==2){
      paralist[[j]][2]<-paralist[[j]][2]*0.99
    }
  }
  bestpara<-paralist[[which.max(fs)]]
  return(list(bestpara=bestpara,fs=fs))
}
