#### Packages
rm(list=ls(all=TRUE))
source('MultiRankCpp.R')
library("ecp")
library("mvtnorm")
library("hdbinseg")
library("InspectChangepoint")
library("factorcpt")
###functions
simulation_setting = function(mean1,m,p,R,omega,n,signal,df=5){
  Theta<-matrix(rnorm(m*p, mean1, 1),m,p)
  Theta_SVD<-svd(Theta)
  d1<-length(Theta_SVD$d)
  Theta_SVD$d[(R+1):d1]<-0
  Theta<-Theta_SVD$u%*%diag(Theta_SVD$d)%*%t(Theta_SVD$v)
  
  #X=rmvnorm(n=n,  mean=rep(0,(m*p)), sigma=omega)
  X=rmvnorm(n=(n*m),  mean=rep(0,p), sigma=omega)
  
  E=rt((n*m), df)
  E=(E-mean(E))/sd(E)
  E=matrix(E,n,m)
  
  Y=matrix(0,n,m)
  
  for (j in 1:n){
    Xj=X[((j-1)*m+1):(j*m),]
    #print(Xj)
    #Y[j,]=rowSums(Theta*t(matrix(X[j,],p,m)))+signal*E[j,]
    Y[j,]=rowSums(Theta*Xj)+signal*E[j,]
  }
  result=list("","")
  result[[1]]=X
  result[[2]]=Y
  return(result)}



evaluate=function(Selected_Index,true_cpt_pos){
  cpt_dist=matrix(0,length(true_cpt_pos),1)
  sel_dist=matrix(0,length(Selected_Index),1)
  for (i in 1:length(true_cpt_pos)){
    dist=abs(Selected_Index-true_cpt_pos[i])
    cpt_dist[i,]=min(dist)
  }
  cpt_dist_max=max(cpt_dist)
  for (j in 1:length(Selected_Index)){
    dist=abs(Selected_Index[j]-true_cpt_pos)
    sel_dist[j,]=min(dist)
  }
  sel_dist_max=max(sel_dist)
  result=list("","")
  result[[1]]=cpt_dist_max
  result[[2]]=sel_dist_max
  return(result) }




####Simulation start
p=8
m=8
#p_total=p*m
rho=0.5
#omega=rho^(abs(matrix(rep(c(1:p_total),each=p_total),ncol=p_total,byrow=T)
#               -matrix(rep(c(1:p_total),each=p_total),ncol=p_total)))

omega=rho^(abs(matrix(rep(c(1:p),each=p),ncol=p,byrow=T)
               -matrix(rep(c(1:p),each=p),ncol=p)))

R1=3

T1=100
T2=75
T3=150
T4=75
T5=100
n=T1+T2+T3+T4+T5
true_cpt_pos=c(T1,(T1+T2),(T1+T2+T3),(T1+T2+T3+T4));
signal=0.5
#mu=0
#sigma=1
Itermax=500
###################
ecp_num=matrix(NA,Itermax,1)
ecp_under=matrix(NA,Itermax,1)
ecp_over=matrix(NA,Itermax,1)

mul_num=matrix(NA,Itermax,1)
mul_under=matrix(NA,Itermax,1)
mul_over=matrix(NA,Itermax,1)

ins_change=matrix(NA,Itermax,500)
ins_num=matrix(NA,Itermax,1)
#ins_num1=matrix(NA,Itermax,1)
ins_under=matrix(NA,Itermax,1)
ins_over=matrix(NA,Itermax,1)

sbs1_num=matrix(NA,Itermax,1)
sbs1_under=matrix(NA,Itermax,1)
sbs1_over=matrix(NA,Itermax,1)

sbs2_change=matrix(NA,Itermax,500)
sbs2_num=matrix(NA,Itermax,1)
sbs2_under=matrix(NA,Itermax,1)
sbs2_over=matrix(NA,Itermax,1)

factor_change=matrix(NA,Itermax,500)
factor_num=matrix(NA,Itermax,1)
factor_under=matrix(NA,Itermax,1)
factor_over=matrix(NA,Itermax,1)

###############

set.seed(2019)
for (Iter in 1:Itermax){
  print(Iter)
  
  X1=simulation_setting(-1,m,p,R1,omega,T1,signal)[[1]]
  Y1=simulation_setting(-1,m,p,R1,omega,T1,signal)[[2]]
  
  X2=simulation_setting(0,m,p,R1,omega,T2,signal)[[1]]
  Y2=simulation_setting(0,m,p,R1,omega,T2,signal)[[2]]
  
  
  X3=simulation_setting(1,m,p,R1,omega,T3,signal)[[1]]
  Y3=simulation_setting(1,m,p,R1,omega,T3,signal)[[2]]
  
  X4=simulation_setting(-1,m,p,R1,omega,T4,signal)[[1]]
  Y4=simulation_setting(-1,m,p,R1,omega,T4,signal)[[2]]
  
  X5=simulation_setting(0,m,p,R1,omega,T5,signal)[[1]]
  Y5=simulation_setting(0,m,p,R1,omega,T5,signal)[[2]]

Y=rbind(Y1,Y2,Y3,Y4,Y5)

####ecp
ecp = e.divisive(Y,R=499,alpha=1)
ecp_num[Iter]=ecp$k.hat-1
 if (ecp_num[Iter]!=0){
   last_index=length(ecp$estimates)
   Selected_Index_ecp=ecp$estimates[-c(1,last_index)]-1
   res=evaluate(Selected_Index_ecp,true_cpt_pos)
   ecp_over[Iter]=res[[1]]
   ecp_under[Iter]=res[[2]]
 }else{
   Selected_Index_ecp=ecp$estimates
   res=evaluate(Selected_Index_ecp,true_cpt_pos)
   ecp_over[Iter]=res[[1]]
   ecp_under[Iter]=res[[2]]
 }
#print(Selected_Index_ecp)


###MultiRank
mul = MultiRank(t(Y),n/30)
mul_num[Iter]=mul$cluster_number-1
if (mul_num[Iter]!=0){
  last_index=length(mul$list)
  Selected_Index_mul=mul$list[-c(1,last_index)]-1
  res=evaluate(Selected_Index_mul,true_cpt_pos)
  mul_over[Iter]=res[[1]]
  mul_under[Iter]=res[[2]]
}else{
  Selected_Index_mul=mul$list
  res=evaluate(Selected_Index_mul,true_cpt_pos)
  mul_over[Iter]=res[[1]]
  mul_under[Iter]=res[[2]]
}
#print(Selected_Index_mul)



###Tengyao Wang and Richard Samworth
n=dim(t(Y))[2]
p=dim(t(Y))[1]
threshold <- compute.threshold(n,p)
ins_res <- inspect(t(Y), threshold = threshold)
#ins_num[Iter]=length(ins_res$changepoints)/3
ins_num[Iter]=length(ins_res$changepoints[,1])
if (ins_num[Iter]!=0){
  ins_change[Iter,1:ins_num[Iter]]=ins_res$changepoints[,1]
  res=evaluate(ins_res$changepoints[,1],true_cpt_pos)
  ins_over[Iter]=res[[1]]
  ins_under[Iter]=res[[2]]
}else{
  ins_se=c(1,n)
  res=evaluate(ins_se,true_cpt_pos)
  ins_over[Iter]=res[[1]]
  ins_under[Iter]=res[[2]]
}


###Haeran Cho and Piotr Fryzlewicz
sbs1=sbs.alg(t(Y), cp.type=1, do.parallel=4)$ecp ##change point in mean
sbs2=sbs.alg(t(Y), cp.type=2, do.parallel=4)$ecp ##change point in second-order
sbs1_num[Iter]=length(sbs1)
if (sbs1_num[Iter]!=0){
  res=evaluate(sbs1,true_cpt_pos)
  sbs1_over[Iter]=res[[1]]
  sbs1_under[Iter]=res[[2]]
}else{
  sbs1_se=c(1,n)
  res=evaluate(sbs1_se,true_cpt_pos)
  sbs1_over[Iter]=res[[1]]
  sbs1_under[Iter]=res[[2]]
}

sbs2_num[Iter]=length(sbs2)
if (sbs2_num[Iter]!=0){
  sbs2_change[Iter,1:sbs2_num[Iter]]=sbs2
  res=evaluate(sbs2,true_cpt_pos)
  sbs2_over[Iter,]=res[[1]]
  sbs2_under[Iter]=res[[2]]
}else{
  sbs2_se=c(1,n)
  res=evaluate(sbs2_se,true_cpt_pos)
  sbs2_over[Iter]=res[[1]]
  sbs2_under[Iter]=res[[2]]
}


###factor
n=dim(t(Y))[2]
p=dim(t(Y))[1]
factor_res=factor.seg.alg(t(Y),max.q=sqrt(min(n,p)))
factor_num[Iter]=length(factor_res$common.est.cps)
if (factor_num[Iter]!=0){
  factor_change[Iter,1:factor_num[Iter]]=factor_res$common.est.cps
  res=evaluate(factor_res$common.est.cps,true_cpt_pos)
  factor_over[Iter]=res[[1]]
  factor_under[Iter]=res[[2]]
}else{
  factor_se=c(1,n)
  res=evaluate(factor_se,true_cpt_pos)
  factor_over[Iter]=res[[1]]
  factor_under[Iter]=res[[2]]
}

#print(factor_res$common.est.cps)
}
save(list=ls(all=TRUE),file="cptmeanT5.Rdata")