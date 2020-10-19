ecp1=length(which(ecp_num<2))
ecp2=length(which(ecp_num==2))
ecp3=length(which(ecp_num==3))
ecp4=length(which(ecp_num==4))
ecp5=length(which(ecp_num==5))
ecp6=length(which(ecp_num==6))
ecp7=length(which(ecp_num>6))

Index1=which(ecp_num==0)

if (length(Index1)!=0){
  ecp8=mean(ecp_over[-Index1,])
  ecp9=sd(ecp_over[-Index1,])
  ecp10=mean(ecp_under[-Index1,])
  ecp11=sd(ecp_under[-Index1,])
}else{
  ecp8=mean(ecp_over)
  ecp9=sd(ecp_over)
  ecp10=mean(ecp_under)
  ecp11=sd(ecp_under)
}



mul1=length(which(mul_num<2))
mul2=length(which(mul_num==2))
mul3=length(which(mul_num==3))
mul4=length(which(mul_num==4))
mul5=length(which(mul_num==5))
mul6=length(which(mul_num==6))
mul7=length(which(mul_num>6))

Index1=which(mul_num==0)

if (length(Index1)!=0){
  mul8=mean(mul_over[-Index1,])
  mul9=sd(mul_over[-Index1,])
  mul10=mean(mul_under[-Index1,])
  mul11=sd(mul_under[-Index1,])
}else{
  mul8=mean(mul_over)
  mul9=sd(mul_over)
  mul10=mean(mul_under)
  mul11=sd(mul_under)
}






ins1=length(which(ins_num<2))
ins2=length(which(ins_num==2))
ins3=length(which(ins_num==3))
ins4=length(which(ins_num==4))
ins5=length(which(ins_num==5))
ins6=length(which(ins_num==6))
ins7=length(which(ins_num>6))

Index1=which(ins_num==0)

if (length(Index1)!=0){
  ins8=mean(ins_over[-Index1,])
  ins9=sd(ins_over[-Index1,])
  ins10=mean(ins_under[-Index1,])
  ins11=sd(ins_under[-Index1,])
}else{
  ins8=mean(ins_over)
  ins9=sd(ins_over)
  ins10=mean(ins_under)
  ins11=sd(ins_under)
}


sbs11=length(which(sbs1_num<2))
sbs12=length(which(sbs1_num==2))
sbs13=length(which(sbs1_num==3))
sbs14=length(which(sbs1_num==4))
sbs15=length(which(sbs1_num==5))
sbs16=length(which(sbs1_num==6))
sbs17=length(which(sbs1_num>6))

Index1=which(sbs1_num==0)

if (length(Index1)!=0){
  sbs18=mean(sbs1_over[-Index1,])
  sbs19=sd(sbs1_over[-Index1,])
  sbs110=mean(sbs1_under[-Index1,])
  sbs111=sd(sbs1_under[-Index1,])
}else{
  sbs18=mean(sbs1_over)
  sbs19=sd(sbs1_over)
  sbs110=mean(sbs1_under)
  sbs111=sd(sbs1_under)
}



sbs21=length(which(sbs2_num<2))
sbs22=length(which(sbs2_num==2))
sbs23=length(which(sbs2_num==3))
sbs24=length(which(sbs2_num==4))
sbs25=length(which(sbs2_num==5))
sbs26=length(which(sbs2_num==6))
sbs27=length(which(sbs2_num>6))


Index1=which(sbs2_num==0)

if (length(Index1)!=0){
  sbs28=mean(sbs2_over[-Index1,])
  sbs29=sd(sbs2_over[-Index1,])
  sbs210=mean(sbs2_under[-Index1,])
  sbs211=sd(sbs2_under[-Index1,])
}else{
  sbs28=mean(sbs2_over)
  sbs29=sd(sbs2_over)
  sbs210=mean(sbs2_under)
  sbs211=sd(sbs2_under)
}




factor1=length(which(factor_num<2))
factor2=length(which(factor_num==2))
factor3=length(which(factor_num==3))
factor4=length(which(factor_num==4))
factor5=length(which(factor_num==5))
factor6=length(which(factor_num==6))
factor7=length(which(factor_num>6))


Index1=which(factor_num==0)
if (length(Index1)!=0){
  factor8=mean(factor_over[-Index1,])
  factor9=sd(factor_over[-Index1,])
  factor10=mean(factor_under[-Index1,])
  factor11=sd(factor_under[-Index1,])
}else{
  factor8=mean(factor_over)
  factor9=sd(factor_over)
  factor10=mean(factor_under)
  factor11=sd(factor_under)
}



###summary

library(pracma)
Table=matrix(0,5,11)
Table[1,]=c(ecp1,ecp2,ecp3,ecp4,ecp5,ecp6,ecp7,ecp8,ecp9,ecp10,ecp11)
Table[2,]=c(mul1,mul2,mul3,mul4,mul5,mul6,mul7,mul8,mul9,mul10,mul11)
Table[3,]=c(ins1,ins2,ins3,ins4,ins5,ins6,ins7,ins8,ins9,ins10,ins11)
#Table[4,]=c(sbs11,sbs12,sbs13,sbs14,sbs15,sbs16,sbs17,sbs18,sbs19,sbs110,sbs111)
Table[4,]=c(sbs21,sbs22,sbs23,sbs24,sbs25,sbs26,sbs27,sbs28,sbs29,sbs210,sbs211)
Table[5,]=c(factor1,factor2,factor3,factor4,factor5,factor6,factor7,factor8,factor9,factor10,factor11)

for (i in 1:5) {
  for (j in 1:7){
    fprintf("%4d & ",Table[i,j] )
  }
  fprintf("%4.2f (%4.2f)&%4.2f (%4.2f)\\\\  ",Table[i,8],Table[i,9],Table[i,10],Table[i,11] )
  fprintf("\n")
}