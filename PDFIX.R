#Logic of codes here takes reference from the 2D Projection depth calculation in the matlab codes of CompPD package, author Xiaohui Liu, School of Statistics, Jiangxi University of Finance and Economics, Nanchang


#stahel-donoho outlyingness
SDO<-function(x,vec){
  u<-vec$u
  Med<-vec$Med
  MAD<-vec$MAD
  d=dim(x)[1]
  if(is.null(d)){
    temp.save<-abs(t(u)%*%x-Med)/MAD
    out<-max(temp.save)
    ulab<-which.max(temp.save)
  } else {
    out<-rep(0,d)
    ulab<-rep(0,d)
    for (i in 1:d) {
      temp.save<-abs(t(u)%*%x[i,]-Med)/MAD
      out[i]<-max(temp.save)
      ulab[i]<-which.max(temp.save)
    }
  }
  return(list(SDO=out,index=ulab))
}

#SDO with a given point
SDO_given<-function(x,vec,med){
  d=dim(x)[1]
  u<-vec$u
  
  out<-rep(0,d)
  ulab<-rep(0,d)
  Med<-t(u)%*%med
  save.mad<-c()
  for (j in 1:length(Med)) {
    utx<-t(u[,j])%*%t(x)
    save.mad[j]<-median(abs(utx-Med[j]))
  }
  
  for (i in 1:d) {
    if (i==298){
      out[i]<-0
      ulab[i]<-0
    } else {
      tx<-t(u)%*%x[i,]
      temp.save<-abs(tx-Med)/save.mad
      out[i]<-max(temp.save)
      ulab[i]<-which.max(temp.save)
    }
  }
  return(list(SDO=out,index=ulab,mad=save.mad))
}

#Projection depth based on SDO
PDV<-function(SDO){
  depth<-1/(1+SDO)
}

#Find vector set for projection depth with a given center
PD_vec_fixmed<-function(X,med){
  X<-t(X)
  p=dim(X)[1]
  size=dim(X)[2]
  u_0<-c(1,0)%>%as.matrix()
  u0x<-t(u_0)%*%X
  umed<-t(u_0)%*%med
  sort_xu<-sort(u0x)
  permu_xum<-sort.list(u0x)
  sortX<-X[,permu_xum]
  med_label<-which(sort_xu==umed[1])
  med_name<-permu_xum[med_label]
  if (med_label==1){
    sub_ind<-c(2:size)
  } else if (med_label==size) {
    sub_ind<-c(1:(med_label-1))
  } else {
    sub_ind<-c(1:(med_label-1),(med_label+1):size)
  }
  MedXM<-sortX[,sub_ind]
  MedXminMedXj<-MedXM-sortX[,med_label]
  slopemed<--MedXminMedXj[1,]/MedXminMedXj[2,]
  theta<-atan(slopemed)+pi*((slopemed<0)%>%as.numeric())
  theta0<-min(theta)
  tmpSubIndxMed<-which.min(theta)
  
  u=c()
  Med=c()
  MAD=c()
  
  #median (initial u)
  if(size%%2!=0){
    mid_ind<-floor(size/2)+1
    mid_indS<-mid_ind-1
    sub_ind2<-c(1:mid_indS,(mid_ind+1):size)
    
    MADX<-sortX-sortX[,med_label]
    if (med_label>1){
      MADX[,1:(med_label-1)]<--MADX[,1:(med_label-1)]
    }
    u0x<-t(u_0)%*%MADX
    sort_xu<-sort(u0x)
    permu_xu<-sort.list(u0x)
    MADX<-MADX[,permu_xu]
    
    #Calculate the MAD sequences lie in interval:[0, theta0) (the angle)
    
    
    MADXMat<-MADX[,sub_ind2]
    MADXMatSubMADXj<-MADXMat-MADX[,mid_ind]
    slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
    MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
    beta1<-min(MADtheta)
    SubindMAD<-which.min(MADtheta)
    while (beta1<theta0) {
      subIndx1MAD<-sub_ind2[SubindMAD]
      
      tmpVecMAD<-MADX[,subIndx1MAD]
      MADX[,subIndx1MAD]<-MADX[,mid_ind]
      MADX[,mid_ind]<-tmpVecMAD
      
      
      #Update arguments u, MedVal and MADVal
      u=cbind(u,c(cos(beta1),sin(beta1)))
      Med=cbind(Med,c(t(u[,dim(u)[2]])%*%sortX[,med_label]))
      MAD=cbind(MAD,c(t(u[,dim(u)[2]])%*%MADX[,mid_ind]))
      beta0=beta1
      
      MADXMat<-MADX[,sub_ind2]
      MADXMatSubMADXj<-MADXMat-MADX[,mid_ind]
      slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
      MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
      MADtheta[which(MADtheta<=beta0)]=2*pi
      beta1<-min(MADtheta)
      SubindMAD<-which.min(MADtheta)
    }
    
    
    # Calculate the Med sequences lie in interval:[theta0, theta1) with theta1 < pi
    
    u=cbind(u,c(cos(theta0),sin(theta0)))
    Med=cbind(Med,c(t(u[,dim(u)[2]])%*%sortX[,med_label]))
    
    subIndx1Med<-sub_ind[tmpSubIndxMed]
    
    med_label_temp<-med_label
    if (subIndx1Med<med_label){
      med_label<-med_label-1
      tmpVecMeD<-sortX[,subIndx1Med]
      tmpVecMeD2<-sortX[,med_label_temp]
      sortX[,med_label_temp]<-tmpVecMeD
      sortX[,subIndx1Med]<-sortX[,med_label]
      sortX[,med_label]<-tmpVecMeD2
    } else if (subIndx1Med>med_label){
      med_label<-med_label+1
      tmpVecMeD<-sortX[,subIndx1Med]
      tmpVecMeD2<-sortX[,med_label_temp]
      sortX[,med_label_temp]<-tmpVecMeD
      sortX[,subIndx1Med]<-sortX[,med_label]
      sortX[,med_label]<-tmpVecMeD2
    }
    
    if (med_label==1){
      sub_ind<-c(2:size)
    } else if (med_label==size) {
      sub_ind<-c(1:(med_label-1))
    } else {
      sub_ind<-c(1:(med_label-1),(med_label+1):size)
    }
    
    IsContinueMed = T
    while (IsContinueMed) {
      MADX<-sortX-sortX[,med_label]
      if (med_label>1){
        MADX[,1:(med_label-1)]<--MADX[,1:(med_label-1)]
      }
      ux<-t(u[,dim(u)[2]])%*%MADX
      sort_xu<-sort(ux)
      permu_xu<-sort.list(ux)
      MADX<-MADX[,permu_xu]
      #update MAD
      MAD <-cbind(MAD,c(t(u[,dim(u)[2]])%*%MADX[,mid_ind]))
      
      MedXM<-sortX[,sub_ind]
      MedXminMedXj<-MedXM-sortX[,med_label]
      slopemed<--MedXminMedXj[1,]/MedXminMedXj[2,]
      theta<-atan(slopemed)+pi*((slopemed<0)%>%as.numeric())
      theta[theta<=theta0]=2*pi
      theta1<-min(theta)
      tmpSubIndxMed<-which.min(theta)
      
      #Calculate the MAD sequences lie in interval: [theta0, theta1_max) with theta1_max < pi
      if (theta1<pi){
        beta0<-theta0
        MADXMat<-MADX[,sub_ind2]
        MADXMatSubMADXj<-MADXMat-MADX[,mid_ind]
        slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
        MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
        MADtheta[which(MADtheta<=beta0)]=2*pi
        beta1<-min(MADtheta)
        SubindMAD<-which.min(MADtheta)
        while (beta1<theta1) {
          subIndx1MAD<-sub_ind2[SubindMAD]
          
          tmpVecMAD<-MADX[,subIndx1MAD]
          MADX[,subIndx1MAD]<-MADX[,mid_ind]
          MADX[,mid_ind]<-tmpVecMAD
          
          #Update arguments u, MedVal and MADVal
          u=cbind(u,c(cos(beta1),sin(beta1)))
          Med=cbind(Med,c(t(u[,dim(u)[2]])%*%sortX[,med_label]))
          MAD=cbind(MAD,c(t(u[,dim(u)[2]])%*%MADX[,mid_ind]))
          beta0=beta1
          
          MADXMat<-MADX[,sub_ind2]
          MADXMatSubMADXj<-MADXMat-MADX[,mid_ind]
          slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
          MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
          MADtheta[which(MADtheta<=beta0)]=2*pi
          beta1<-min(MADtheta)
          SubindMAD<-which.min(MADtheta)
        }
        theta0<-theta1
        u=cbind(u,c(cos(theta0),sin(theta0)))
        Med=cbind(Med,c(t(u[,dim(u)[2]])%*%sortX[,med_label]))
        
        subIndx1Med<-sub_ind[tmpSubIndxMed]
        
        med_label_temp<-med_label
        if (subIndx1Med<med_label){
          med_label<-med_label-1
          tmpVecMeD<-sortX[,subIndx1Med]
          tmpVecMeD2<-sortX[,med_label_temp]
          sortX[,med_label_temp]<-tmpVecMeD
          sortX[,subIndx1Med]<-sortX[,med_label]
          sortX[,med_label]<-tmpVecMeD2
        }else if (subIndx1Med>med_label){
          med_label<-med_label+1
          tmpVecMeD<-sortX[,subIndx1Med]
          tmpVecMeD2<-sortX[,med_label_temp]
          sortX[,med_label_temp]<-tmpVecMeD
          sortX[,subIndx1Med]<-sortX[,med_label]
          sortX[,med_label]<-tmpVecMeD2
        }
        
        if (med_label==1){
          sub_ind<-c(2:size)
        } else if (med_label==size) {
          sub_ind<-c(1:(med_label-1))
        } else {
          sub_ind<-c(1:(med_label-1),(med_label+1):size)
        }
        
      } else {
        IsContinueMed<-F
      }
    }
    # Calculate the MAD sequences lie in interval:[theta1_max, pi) with theta1_max < pi
    beta0<-theta0
    MADXMat<-MADX[,sub_ind2]
    MADXMatSubMADXj<-MADXMat-MADX[,mid_ind]
    slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
    MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
    MADtheta[which(MADtheta<=beta0)]=2*pi
    beta1<-min(MADtheta)
    SubindMAD<-which.min(MADtheta)
    while (beta1<theta1) {
      subIndx1MAD<-sub_ind2[SubindMAD]
      tmpVecMAD<-MADX[,subIndx1MAD]
      MADX[,subIndx1MAD]<-MADX[,mid_ind]
      MADX[,mid_ind]<-tmpVecMAD
      
      #Update arguments u, MedVal and MADVal
      u=cbind(u,c(cos(beta1),sin(beta1)))
      Med=cbind(Med,c(t(u[,dim(u)[2]])%*%sortX[,med_label]))
      MAD=cbind(MAD,c(t(u[,dim(u)[2]])%*%MADX[,mid_ind]))
      beta0=beta1
      
      MADXMat<-MADX[,sub_ind2]
      MADXMatSubMADXj<-MADXMat-MADX[,mid_ind]
      slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
      MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
      MADtheta[which(MADtheta<=beta0)]=2*pi
      beta1<-min(MADtheta)
      SubindMAD<-which.min(MADtheta)
    }
  } else {
    mid_ind<-size/2
    mid_indS<-mid_ind-1
    mid_indP1<-mid_ind+1
    mid_indP2<-mid_indP1+1
    sizeS2<-size-2
    sub_ind2<-c(1:mid_indS,mid_indP2:size,1:mid_indS,mid_indP2:size)
    
    MADX<-sortX-sortX[,med_label]
    if (med_label>1){
      MADX[,1:(med_label-1)]<--MADX[,1:(med_label-1)]
    }
    u0x<-t(u_0)%*%MADX
    sort_xu<-sort(u0x)
    permu_xu<-sort.list(u0x)
    MADX<-MADX[,permu_xu]
    #Calculate the MAD sequences lie in interval:[0, theta0) (the angle)
    
    MADXMat<-MADX[,sub_ind2[1:sizeS2]]
    MADXMatSubMADXj<-cbind(MADXMat-MADX[,mid_ind],MADXMat-MADX[,mid_indP1])
    slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
    MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
    beta1<-min(MADtheta)
    SubindMAD<-which.min(MADtheta)
    while (beta1<theta0) {
      subIndx1MAD<-sub_ind2[SubindMAD]
      subIndx2MAD<-mid_ind+(SubindMAD>sizeS2)%>%as.numeric()
      
      tmpVecMAD<-MADX[,subIndx1MAD]
      MADX[,subIndx1MAD]<-MADX[,subIndx2MAD]
      MADX[,subIndx2MAD]<-tmpVecMAD
      
      MADXj = (MADX[,mid_ind] + MADX[,mid_indP1]) / 2
      #Update arguments u, MedVal and MADVal
      u=cbind(u,c(cos(beta1),sin(beta1)))
      Med=cbind(Med,c(t(u[,dim(u)[2]])%*%sortX[,med_label]))
      MAD=cbind(MAD,c(t(u[,dim(u)[2]])%*%MADXj))
      beta0=beta1
      
      MADXMat<-MADX[,sub_ind2[1:sizeS2]]
      MADXMatSubMADXj<-cbind(MADXMat-MADX[,mid_ind],MADXMat-MADX[,mid_indP1])
      slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
      MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
      MADtheta[which(MADtheta<=beta0)]=2*pi
      beta1<-min(MADtheta)
      SubindMAD<-which.min(MADtheta)
    }
    #update permu
    
    u=cbind(u,c(cos(theta0),sin(theta0)))
    Med=cbind(Med,c(t(u[,dim(u)[2]])%*%sortX[,med_label]))
    
    subIndx1Med<-sub_ind[tmpSubIndxMed]
    
    
    med_label_temp<-med_label
    if (subIndx1Med<med_label){
      med_label<-med_label-1
      tmpVecMeD<-sortX[,subIndx1Med]
      tmpVecMeD2<-sortX[,med_label_temp]
      sortX[,med_label_temp]<-tmpVecMeD
      sortX[,subIndx1Med]<-sortX[,med_label]
      sortX[,med_label]<-tmpVecMeD2
    } else if (subIndx1Med>med_label){
      med_label<-med_label+1
      tmpVecMeD<-sortX[,subIndx1Med]
      tmpVecMeD2<-sortX[,med_label_temp]
      sortX[,med_label_temp]<-tmpVecMeD
      sortX[,subIndx1Med]<-sortX[,med_label]
      sortX[,med_label]<-tmpVecMeD2
    }
    
    if (med_label==1){
      sub_ind<-c(2:size)
    } else if (med_label==size) {
      sub_ind<-c(1:(med_label-1))
    } else {
      sub_ind<-c(1:(med_label-1),(med_label+1):size)
    }
    
    
    IsContinueMed = T
    while (IsContinueMed) {
      MADX<-sortX-sortX[,med_label]
      if (med_label>1){
        MADX[,1:(med_label-1)]<--MADX[,1:(med_label-1)]
      }
      ux<-t(u[,dim(u)[2]])%*%MADX
      sort_xu<-sort(round(ux,10))
      permu_xu<-sort.list(ux)
      MADX<-MADX[,permu_xu]
      #update MAD
      MADXj = (MADX[,mid_ind] + MADX[,mid_indP1]) / 2
      MAD <-cbind(MAD,c(t(u[,dim(u)[2]])%*%MADXj))
      
      MedXM<-sortX[,sub_ind]
      MedXminMedXj<-MedXM-sortX[,med_label]
      slopemed<--MedXminMedXj[1,]/MedXminMedXj[2,]
      theta<-atan(slopemed)+pi*((slopemed<0)%>%as.numeric())
      theta[theta<=theta0]=2*pi
      theta1<-min(theta)
      tmpSubIndxMed<-which.min(theta)
      #Calculate the MAD sequences lie in interval: [theta0, theta1_max) with theta1_max < pi
      if (theta1<pi){
        beta0<-theta0
        MADXMat<-MADX[,sub_ind2[1:sizeS2]]
        MADXMatSubMADXj<-cbind(MADXMat-MADX[,mid_ind],MADXMat-MADX[,mid_indP1])
        slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
        MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
        MADtheta[which(MADtheta<=beta0)]=2*pi
        beta1<-min(MADtheta)
        SubindMAD<-which.min(MADtheta)
        while (beta1<theta1) {
          subIndx1MAD<-sub_ind2[SubindMAD]
          subIndx2MAD<-mid_ind+(SubindMAD>sizeS2)%>%as.numeric()
          
          tmpVecMAD<-MADX[,subIndx1MAD]
          MADX[,subIndx1MAD]<-MADX[,subIndx2MAD]
          MADX[,subIndx2MAD]<-tmpVecMAD
          
          MADXj = (MADX[,mid_ind] + MADX[,mid_indP1]) / 2
          #Update arguments u, MedVal and MADVal
          u=cbind(u,c(cos(beta1),sin(beta1)))
          Med=cbind(Med,c(t(u[,dim(u)[2]])%*%sortX[,med_label]))
          MAD=cbind(MAD,c(t(u[,dim(u)[2]])%*%MADXj))
          beta0=beta1
          
          MADXMat<-MADX[,sub_ind2[1:sizeS2]]
          MADXMatSubMADXj<-cbind(MADXMat-MADX[,mid_ind],MADXMat-MADX[,mid_indP1])
          slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
          MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
          MADtheta[which(MADtheta<=beta0)]=2*pi
          beta1<-min(MADtheta)
          SubindMAD<-which.min(MADtheta)
        }
        theta0<-theta1
        u=cbind(u,c(cos(theta0),sin(theta0)))
        Med=cbind(Med,c(t(u[,dim(u)[2]])%*%sortX[,med_label]))
        
        subIndx1Med<-sub_ind[tmpSubIndxMed]
        
        med_label_temp<-med_label
        if (subIndx1Med<med_label){
          med_label<-med_label-1
          tmpVecMeD<-sortX[,subIndx1Med]
          tmpVecMeD2<-sortX[,med_label_temp]
          sortX[,med_label_temp]<-tmpVecMeD
          sortX[,subIndx1Med]<-sortX[,med_label]
          sortX[,med_label]<-tmpVecMeD2
        }else if (subIndx1Med>med_label){
          med_label<-med_label+1
          tmpVecMeD<-sortX[,subIndx1Med]
          tmpVecMeD2<-sortX[,med_label_temp]
          sortX[,med_label_temp]<-tmpVecMeD
          sortX[,subIndx1Med]<-sortX[,med_label]
          sortX[,med_label]<-tmpVecMeD2
        }
        
        
        if (med_label==1){
          sub_ind<-c(2:size)
        } else if (med_label==size) {
          sub_ind<-c(1:(med_label-1))
        } else {
          sub_ind<-c(1:(med_label-1),(med_label+1):size)
        }
        
      } else {
        IsContinueMed<-F
      }
    }
    #Calculate the MAD sequences lie in interval: [theta1_max, pi) with theta1_max < pi
    beta0<-theta0
    MADXMat<-MADX[,sub_ind2[1:sizeS2]]
    MADXMatSubMADXj<-cbind(MADXMat-MADX[,mid_ind],MADXMat-MADX[,mid_indP1])
    slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
    MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
    MADtheta[which(MADtheta<=beta0)]=2*pi
    beta1<-min(MADtheta)
    SubindMAD<-which.min(MADtheta)
    while (beta1<theta1) {
      subIndx1MAD<-sub_ind2[SubindMAD]
      subIndx2MAD<-mid_ind+(SubindMAD>sizeS2)%>%as.numeric()
      
      tmpVecMAD<-MADX[,subIndx1MAD]
      MADX[,subIndx1MAD]<-MADX[,subIndx2MAD]
      MADX[,subIndx2MAD]<-tmpVecMAD
      
      MADXj = (MADX[,mid_ind] + MADX[,mid_indP1]) / 2
      #Update arguments u, MedVal and MADVal
      u=cbind(u,c(cos(beta1),sin(beta1)))
      Med=cbind(Med,c(t(u[,dim(u)[2]])%*%sortX[,med_label]))
      MAD=cbind(MAD,c(t(u[,dim(u)[2]])%*%MADXj))
      beta0=beta1
      
      MADXMat<-MADX[,sub_ind2[1:sizeS2]]
      MADXMatSubMADXj<-cbind(MADXMat-MADX[,mid_ind],MADXMat-MADX[,mid_indP1])
      slopemad<--MADXMatSubMADXj[1,]/MADXMatSubMADXj[2,]
      MADtheta<-atan(slopemad)+pi*((slopemad<0)%>%as.numeric())
      MADtheta[which(MADtheta<=beta0)]=2*pi
      beta1<-min(MADtheta)
      SubindMAD<-which.min(MADtheta)
    }
  }
  output<-list(u=cbind(u,-u),Med=c(Med,-Med),MAD=c(MAD,MAD),NumDirVecU=2*length(Med))
  return(output)
}

#combine SDO and PDV, for convenience
ex2dprojdepth<-function(X,x,vec){
  SDOvalue<-SDO(x,vec)
  PD_value<-PDV(SDOvalue$SDO)
  return(PD_value)
}


#Use top two dimensions of sim_data

test.x<-sim_data[,1:2]%>%as.matrix()

#one example
sort_var<-sort.list(test.x[,1],decreasing = T)
med<-test.x[sort_var[2],]
fixedvec2<-PD_vec_fixmed(test.x[1:1000,],med)
SDO_fixed2<-SDO(test.x,fixedvec2)
PDV_fixed2<-PDV(SDO_fixed2$SDO)

#depth.contours
library(metR)
library(ddalpha)
tsdata<-expand.grid(seq(-1.05,0,length=200),seq(-0.54,0.55,length=200))%>%as.matrix()
z<-SDO(tsdata,fixedvec2)
zd<-PDV(z$SDO)

med2<-test.x[which.min(test.x[,2]),]
fixedvec3<-PD_vec_fixmed(test.x,med2)
z2<-SDO(tsdata,fixedvec3)
zd2<-PDV(z2$SDO)

z3<-depth.projection(test.x,test.x)
zd3<-depth.projection(tsdata,test.x)

med3<-test.x[which.max(z3),]
fixedvec4<-PD_vec_fixmed(test.x,med3)
z4<-SDO(tsdata,fixedvec4)
zd4<-PDV(z4$SDO)
tsdata<-cbind(tsdata,zd,zd2,zd3,zd4)%>%as.data.frame()
ggplot()+geom_point(data = test.x%>%as.data.frame(),aes(x=V1,y=V2))+
  geom_contour(data=tsdata,colour='purple',aes(x = Var1, y = Var2, z = zd))+
  geom_contour(data=tsdata,colour='blue',aes(x = Var1, y = Var2, z = zd2))+geom_text_contour(data=tsdata,aes(x = Var1, y = Var2, z = zd))+geom_text_contour(data=tsdata,aes(x = Var1, y = Var2, z = zd2))

ggplot()+geom_point(data = test.x%>%as.data.frame(),aes(x=V1,y=V2))+
  geom_contour(data=tsdata,colour='purple',aes(x = Var1, y = Var2, z = zd3))+
  geom_text_contour(data=tsdata,aes(x = Var1, y = Var2, z = zd3))+geom_contour(data=tsdata,colour='green',aes(x = Var1, y = Var2, z = zd4))+
  geom_text_contour(data=tsdata,aes(x = Var1, y = Var2, z = zd4))

