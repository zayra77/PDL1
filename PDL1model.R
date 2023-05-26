setwd("X:/wanghua/PDL1")
rm(list=ls())
library(mlbench)
library(elasticnet)
library(caret)
library(dplyr)
library(glmnet)
library(pROC)
library(survival)
library(dbplyr)
library(devtools)

data.clinic<-read.csv("./csv/data.clinical.csv")
data.spectral<-read.csv("./csv/data.spectral.csv")
data.CT<-read.csv("./csv/data.radiomicCT.csv")
data.IM<-read.csv("./csv/data.radiomicIM.csv")
data.VNC<-read.csv("./csv/data.radiomicVNC.csv")
na.row<-rownames(data.clinic[!complete.cases(data.clinic),])  ##display na rows
which(rownames(data.clinic)=="610289")#74
which(rownames(data.clinic)=="625271")#137
PDL1<-data.CT$PD_L1
PDL1.191<-PDL1[-c(74,137)]
data.clinic.pdl1<-read.csv("./csv/data.clinical.csv")

####reader agreement
library(irr)
reader2.CT<-read.csv("./csv/Radiomics_observation1.csv")
reader2.IM<-read.csv("./csv/Radiomics_observation2.csv")
reader2.VNC<-read.csv("./csv/Radiomics_observation3.csv")
rownames(reader2.CT)<-reader2.CT$ID
rownames(reader2.IM)<-reader2.IM$ID
rownames(reader2.VNC)<-reader2.VNC$ID
reader2.CT<-reader2.CT[,-1]
reader2.IM<-reader2.IM[,-1]
reader2.VNC<-reader2.VNC[,-1]

rownames(data.CT)<-data.CT$ID
rownames(data.IM)<-data.IM$ID
rownames(data.VNC)<-data.VNC$ID

reader1.CT<-data.CT[rownames(reader2.CT),]
reader1.IM<-data.IM[rownames(reader2.IM),]
reader1.VNC<-data.VNC[rownames(reader2.VNC),]
reader1.CT<-reader1.CT[,-(1:2)]
reader1.IM<-reader1.IM[,-(1:2)]
reader1.VNC<-reader1.VNC[,-(1:2)]
icc_val<-c()
for (i in 1:length(reader1.IM)){
  ratings<-cbind(reader1.IM[,i],reader2.IM[,i])
  icc<-icc(ratings,model = "twoway",
           type="agreement",unit="single",r0 = 0,conf.level = 0.95)
  icc_val[i]<-icc$value
  }

index<-which(icc_val>0.75) #ICC >0.75 

#####radiomics features with agreement>0.75 will be used for further analysis
data.CT<-data.CT[,index]
data.IM<-data.IM[,index]
data.VNC<-data.VNC[,index]
sapply(data.CT, class) #check the data format: numeric
data.CT<- apply(data.CT,2,            # Specify own function within apply
                function(x) as.numeric(as.character(x)))

data.IM<- apply(data.IM,2,            # Specify own function within apply
                function(x) as.numeric(as.character(x)))

data.VNC<- apply(data.VNC,2,            # Specify own function within apply
                 function(x) as.numeric(as.character(x)))
data.spectral<- apply(data.spectral,2,            # Specify own function within apply
                      function(x) as.numeric(as.character(x)))
data.clinic<- apply(data.clinic,2,            # Specify own function within apply
                function(x) as.numeric(as.character(x)))

data.CT<-data.frame(data.CT)
data.IM<-data.frame(data.IM)
data.VNC<-data.frame(data.VNC)
data.spectral<-data.frame(data.spectral)
data.clinic<-data.frame(data.clinic)

rownames(data.CT)<-data.CT$ID
rownames(data.IM)<-data.IM$ID
rownames(data.VNC)<-data.VNC$ID
rownames(data.spectral)<-data.spectral$ID
rownames(data.clinic)<-data.clinic$ID
data.CT<-data.CT[,-(1:2)]
data.IM<-data.IM[,-(1:2)]
data.VNC<-data.VNC[,-(1:2)]
data.spectral<-data.spectral[,-(1:2)]
data.clinic<-data.clinic[,-(1:2)]
#########combine CT IM and VNC into one dataframe
colnames(data.CT)<-paste("CT",colnames(data.CT),sep = "_")
colnames(data.IM)<-paste("IM",colnames(data.IM),sep = "_")
colnames(data.VNC)<-paste("VNC",colnames(data.VNC),sep = "_")
data.radiomic<-cbind(data.CT,data.IM,data.VNC)
#######scale data
data.CT.scale<-as.data.frame(scale(data.CT),center=TRUE,scale=TRUE)
data.IM.scale<-as.data.frame(scale(data.IM),center=TRUE,scale=TRUE)
data.VNC.scale<-as.data.frame(scale(data.VNC),center=TRUE,scale=TRUE)
data.radiomic.scale<-as.data.frame(scale(data.radiomic),center=TRUE,scale=TRUE) #default are TRUE for both
data.spectral.scale<-as.data.frame(scale(data.spectral))
######remove high correlated features###

highcorremove<-function(df){
  library(Hmisc)
  df.nona<-df[,colSums(is.na(df))==0]##remove NA columns
  cor.radiomic<-cor(df.nona)
  cor_matrix_zero <- cor.radiomic                 # Modify correlation matrix
  cor_matrix_zero[upper.tri(cor_matrix_zero)] <- 0
  diag(cor_matrix_zero) <- 0
  df.ind<- df.nona[,!apply(cor_matrix_zero,2,function(x) any(abs(x)> 0.9))]
  return(df.ind)
}
data.CT.ind<-highcorremove(data.CT.scale)
data.radiomic.ind<-highcorremove(data.radiomic.scale)
data.IM.ind<-highcorremove(data.IM.scale)
data.VNC.ind<-highcorremove(data.VNC.scale)
data.CT.191<-data.CT.ind[-c(74,137),]
data.IM.191<-data.IM.ind[-c(74,137),]
data.VNC.191<-data.VNC.ind[-c(74,137),]
data.radiomic.191<-data.radiomic.ind[-c(74,137),]
data.clinic.191<-data.clinic[-c(74,137),]
data.spectral.191<-data.spectral.scale[-c(74,137),]

########sample data to create train and test

set.seed(6397)

Samples<-sample(c(rep(0, 0.8*nrow(data.radiomic.191)), rep(1, 0.2*nrow(data.radiomic.191))))

Train.CT<-data.CT.191[Samples==0,]
Test.CT<-data.CT.191[Samples==1,]

Train.IM<-data.IM.191[Samples==0,]
Test.IM<-data.IM.191[Samples==1,]

Train.VNC<-data.VNC.191[Samples==0,]
Test.VNC<-data.VNC.191[Samples==1,]

Train.radiomic<-data.radiomic.191[Samples==0,]
Test.radiomic<-data.radiomic.191[Samples==1,]

Train.spectral<-data.spectral.191[Samples==0,]
Test.spectral<-data.spectral.191[Samples==1,]
Train.clinic<-data.clinic.191[Samples==0,]
Test.clinic<-data.clinic.191[Samples==1,]
y_train = PDL1.191[Samples==0]
y_test=PDL1.191[Samples==1]

######t test to select the features between pdl expression
source("X:/yzc_r_library.R")

ttest.CT<-try(yzc.ttest(Train.CT,0.05,y_train),TRUE)
if("try-error" %in% class(ttest.CT))
{
  next
}
ttest.VNC<-try(yzc.ttest(Train.VNC,0.05,y_train),TRUE)
if("try-error" %in% class(ttest.VNC))
{
  next
}
ttest.IM<-try(yzc.ttest(Train.IM,0.05,y_train),TRUE)
if("try-error" %in% class(ttest.IM))
{
  next
}
ttest.radiomic<-try(yzc.ttest(Train.radiomic,0.05,y_train),TRUE)
if("try-error" %in% class(ttest.radiomic))
{
  next
}
ttest.spectral<-yzc.ttest(Train.spectral,0.05,y_train)#NA

#####radiomic risk score base on lasso
set.seed(10)
mydf<-Train.VNC
mydf.test<-Test.VNC
y<-y_train
sfeature<-ttest.VNC$ID
#yzc.lasso.rad<-function(mydf,mydf.test,y,sfeature){
y_train<-y
fitCV<- cv.glmnet(data.matrix(mydf[,sfeature]),y_train,alpha=1,family = "binomial",type.measure = "mse",nfolds = 10)
pdf("./nomogramMethod/Spectralmse.pdf")
plot(fitCV)
dev.off()
pdf("./nomogramMethod/spectrallasso.pdf")
plot(fitCV$glmnet.fit,"lambda",label = FALSE)
dev.off()
best_lambda <- fitCV$lambda.min
best_model <- glmnet(data.matrix(mydf[,sfeature]),y_train, alpha = 1, lambda = best_lambda)
coef<-coef(best_model)[which(coef(best_model)!=0),]
print(names(coef))
lasso.df.train<-mydf[,names(coef)[-1]]
lasso.df.test<-mydf.test[,names(coef)[-1]]
temp.test<-mapply(`*`, lasso.df.test, coef[-1])
radscore.test<-rowSums(temp.test)+coef[1]
temp.train<-mapply(`*`, lasso.df.train, coef[-1])
radscore.train<-rowSums(temp.train)+coef[1]
rad.df<-list(radtest=radscore.test,radtrain=radscore.train)
radscore.radiomic<-try(yzc.lasso.rad(Train.radiomic,Test.radiomic,y_train,ttest.radiomic$ID),TRUE)
if("try-error" %in% class(radscore.radiomic))
{next}
radscore.IM<-try(yzc.lasso.rad(Train.IM,Test.IM,y_train,ttest.IM$ID),TRUE)
if("try-error" %in% class(radscore.IM))
{next}
radscore.VNC<-try(yzc.lasso.rad(Train.VNC,Test.VNC, y_train,ttest.VNC$ID),TRUE)
if("try-error" %in% class(radscore.VNC))
{next}

radscore.CT<-try(yzc.lasso.rad(Train.CT,Test.CT,y_train,ttest.CT$ID),TRUE)
if("try-error" %in% class(radscore.CT))
{next}
radscore.Spectral<-try(yzc.lasso.rad(Train.spectral,Test.spectral,y_train,ttest.spectral$ID),TRUE)#or colnames(Train.spectral)
if("try-error" %in% class(radscore.Spectral))
{next}

#############################
Test.clinic$radscore_IM<-radscore.IM[[1]]
Test.clinic$radscore_CT<-radscore.CT[[1]]
Test.clinic$radscore_VNC<-radscore.VNC[[1]]
Train.clinic$radscore_IM<-radscore.IM[[2]]
Train.clinic$radscore_PCI<-radscore.CT[[2]]
Train.clinic$radscore_VNC<-radscore.VNC[[2]]
Train.clinic<- lapply(Train.clinic,factor)
Train.clinic<-as.data.frame(Train.clinic)
sapply(Train.clinic, class)
Train.clinic$radscore_radiomic<-as.numeric(as.character(Train.clinic$radscore_radiomic))
Train.clinic$Age<-as.numeric(as.character(Train.clinic$Age))
Train.clinic$radscore_IM<-as.numeric(as.character(Train.clinic$radscore_IM))
Train.clinic$radscore_PCI<-as.numeric(as.character(Train.clinic$radscore_PCI))
Train.clinic$radscore_VNC<-as.numeric(as.character(Train.clinic$radscore_VNC))
Train.clinic$Papillary<-as.numeric(as.character(Train.clinic$Papillary))
Train.clinic$Acinar<-as.numeric(as.character(Train.clinic$Acinar))
Train.clinic<-within(Train.clinic,rm(Mucinous))
Train.clinic<-within(Train.clinic,rm(radscore_radiomic))
#univariate
Train.clinic$radscore_PCI<-Train.clinic$radscore_PCI*5
Train.clinic$radscore_IM<-Train.clinic$radscore_IM*5
Train.clinic$radscore_VNC<-Train.clinic$radscore_VNC*5
Train.clinic<-within(Train.clinic,rm(Solid))
univariate.result.sig.all<-yzc.univariate(Train.clinic,y_train,0.05)
multivariate.analysis.IM<-glm(y_train~pTNM+Heterogeneity+radscore_IM, data=Train.clinic,family = binomial)
multivariate.analysis.CT<-glm(y_train~pTNM+Heterogeneity+radscore_CT, data=Train.clinic,family = binomial)
multivariate.analysis.VNC<-glm(y_train~pTNM+Heterogeneity+radscore_VNC, data=Train.clinic,family = binomial)
multivariate.analysis.Spectral<-glm(y_train~pTNM+Heterogeneity+radscore_Spectral, data=Train.clinic,family = binomial)
multivariate.analysis.radiomic<-glm(y_train~pTNM+radscore_radiomic, data=Train.clinic,family = binomial)
multivariate.analysis.CT.VNC.IM.Spe<-glm(y_train~pTNM+Heterogeneity+radscore_VNC+radscore_CT+radscore_IM+radscore_Spectral, data=Train.clinic,family = binomial)
multivariate.analysis.CT.IM<-glm(y_train~pTNM+radscore_CT+radscore_IM, data=Train.clinic,family = binomial)
multivariate.analysis.CT.VNC<-glm(y_train~pTNM+radscore_CT+radscore_VNC, data=Train.clinic,family = binomial)
multivariate.analysis.clinical<-glm(y_train~pTNM+Heterogeneity, data=Train.clinic,family = binomial)

multivariate.analysis.IM.VNC<-glm(y_train~pTNM+radscore_VNC+radscore_IM,data=Train.clinic)
multivariate.analysis.IM.VNC.CT<-glm(y_train~pTNM+radscore_VNC+radscore_IM+radscore_CT,data=Train.clinic)
multivariate.analysis.IM.spectral.CT<-glm(y_train~pTNM+radscore_Spectral+radscore_IM+radscore_CT,data=Train.clinic)
multivariate.analysis.spectral.radiomic<-glm(y_train~pTNM+radscore_Spectral+radscore_radiomic,data=Train.clinic)
multivariate.analysis.IM.spectral<-glm(y_train~pTNM+radscore_Spectral+radscore_IM,data=Train.clinic)
multivariate.analysis.CT.spectral<-glm(y_train~pTNM+radscore_Spectral+radscore_CT,data=Train.clinic)
multivariate.analysis.VNC.spectral<-glm(y_train~pTNM+radscore_Spectral+radscore_VNC,data=Train.clinic)

summary(multivariate.analysis.IM)
summary(multivariate.analysis.CT)
summary(multivariate.analysis.VNC)
summary(multivariate.analysis.clinical)
summary(multivariate.analysis.Spectral)
summary(multivariate.analysis.radiomic)
summary(multivariate.analysis.CT.VNC.IM.Spe)
summary(multivariate.analysis.CT.IM)
summary(multivariate.analysis.IM.VNC)
summary(multivariate.analysis.IM.VNC.CT)
summary(multivariate.analysis.IM.spectral.CT)
summary(multivariate.analysis.spectral.radiomic)
summary(multivariate.analysis.IM.spectral)
summary(multivariate.analysis.CT.spectral)
summary(multivariate.analysis.VNC.spectral)
summary(multivariate.analysis.CT.VNC)
#####forest plot
library(questionr)
forestplot.clinical<-odds.ratio(multivariate.analysis.clinical)
forestplot.clinical$Variable<-rownames(forestplot.clinical)
forestplot.clinical<-forestplot.clinical[-1,]
clinical.pvalue<-round(forestplot.clinical$p,3)
colnames(forestplot.clinical)<-c("OR_mean","Lower","Upper","P value","Variable")
forestplot.clinical<-data.frame(lapply(forestplot.clinical, function(x)if(is.numeric(x))round(x,3) else x))
forestplot.clinical$P.value<-clinical.pvalue
forestplot.clinical$OR<-paste(forestplot.clinical$OR_mean,"[",forestplot.clinical$Lower,"-",forestplot.clinical$Upper,"]")
write.csv(forestplot.clinical,file="./nomogramMethod/clinical.multi.forestplot.csv")
library(forestplot)
colnames(univariate.result.sig.all)[1]<-"OR_mean"
colnames(univariate.result.sig.all)[2]<-"Lower"
colnames(univariate.result.sig.all)[3]<-"Upper"
colnames(univariate.result.sig.all)[4]<-"P value"
univariate.result.sig.all.round <- data.frame(lapply(univariate.result.sig.all,    # Using Base R functions
                                 function(x) if(is.numeric(x)) round(x, 3) else x))
univariate.result.sig.all.round$OR<-paste(univariate.result.sig.all.round$OR_mean,"[",univariate.result.sig.all.round$Lower,"-",univariate.result.sig.all.round$Upper,"]")
univariate.result.sig.all.round$Variable<-rownames(univariate.result.sig.all)
write.csv(univariate.result.sig.all.round,file="./nomogramMethod/uniresult.csv")
trellis.device(device="windows", height = 6, width = 25, color=TRUE)

df.forestplot<-read.csv("./nomogramMethod/PCI.multi.forestplot.csv")
df.forestplot<-df.forestplot[,-1]

forestplot(labeltext=as.matrix(df.forestplot[,c(4,5,3)]),
           mean=df.forestplot$OR_mean,
           lower=df.forestplot$Lower,
           upper=df.forestplot$Upper,
           #is.summary=c(T,T,T,F,F,T,F,F,T,T,T,T,T,T,F,F,T,F,F,T,F,F,T,T,T,T,T),
           graph.pos=2, 
           zero=1, 
           graphwidth=unit(40,"mm"),
           lwd.ci = 2,# 误差条的线的宽度
           ci.vertices.height = 0.05, # # 误差条末端的长度
           lineheight="auto", #线的高度eg:unit(20,"cm")
           line.margin=0.1,
           boxsize=0.3,  ##误差条中的圆心点大小
          # xlog=TRUE,#x轴的坐标取对数
           #   xlim=c(0,50),
           lty.ci = 7,# 误差条的线的线型
           xticks.digits = 1,
           xticks=c(1,2,4,6,8),
           fn.ci_norm = fpDrawCircleCI,#误差条显示方式
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.6), xlab = gpar(cex = 0.8), cex = 0.9), #文本大小
           col= fpColors(line = "#CC79A7", #误差条的线的线型
                         box="#D55E00")#误差条的线的宽度
           )

#dev.off()      
####静态nomogram
library(rms)
library(stringr)
class.name<-c("IM","VNC","CT","Spectral","radiomic")
classRad.test<-c(radscore.IM[1],radscore.VNC[1],radscore.CT[1],radscore.Spectral[1],radscore.radiomic[1])
classRad.train<-c(radscore.IM[2],radscore.VNC[2],radscore.CT[2],radscore.Spectral[2],radscore.radiomic[2])
result.allclass<-vector("list",5)

radscore.test<-radscore.VNC[[1]]
radscore.train<-radscore.VNC[[2]]
#single
data.nom.train.VNC<-data.frame(PDL=y_train,Radscore_VNC=unlist(radscore.VNC[[2]]),pTNM=Train.clinic$pTNM)#TNM=Train.clinic$pTNM
data.nom.test.VNC<-data.frame(PDL=y_test,Radscore_VNC=unlist(radscore.VNC[[1]]),pTNM=Test.clinic$pTNM)#TNM=Test.clinic$pTNM
data.nom.train.PCI<-data.frame(PDL=y_train,Radscore_PCI=unlist(radscore.CT[[2]]),pTNM=Train.clinic$pTNM)#TNM=Train.clinic$pTNM
data.nom.test.PCI<-data.frame(PDL=y_test,Radscore_PCI=unlist(radscore.CT[[1]]),pTNM=Test.clinic$pTNM)#TNM=Test.clinic$pTNM
data.nom.train.IM<-data.frame(PDL=y_train,Radscore_IM=unlist(radscore.IM[[2]]),pTNM=Train.clinic$pTNM)#TNM=Train.clinic$pTNM
data.nom.test.IM<-data.frame(PDL=y_test,Radscore_IM=unlist(radscore.IM[[1]]),pTNM=Test.clinic$pTNM)#TNM=Test.clinic$pTNM
data.nom.train.clinical<-data.frame(PDL=y_train,pTNM=Train.clinic$pTNM,Heterogeneity=Train.clinic$Heterogeneity)#TNM=Train.clinic$pTNM
data.nom.test.clinical<-data.frame(PDL=y_test,pTNM=Test.clinic$pTNM,Heterogeneity=Test.clinic$Heterogeneity)#TNM=Train.clinic$pTNM

#multiple

data.nom.test1<-data.nom.test.IM

##使用Platt法和Isotonic Regression一起对预测概率进行校正
library(platt)

isoreg.data<-isoreg(data.nom.test$Radscore_IM) #isotonic choosed in this study
platt.data<-plattScaling(isoreg.data$yf,data.nom.test$PDL)

data.nom.test.IM$Radscore_IM<-platt.data$pred
ddist <- datadist(data.nom.train.VNC)
options(datadist = 'ddist')
logit_lrm<-lrm(PDL~Radscore_VNC+pTNM,data = data.nom.train.VNC)
vif(logit_lrm)  # calculate the VIF

pdf("./nomogramMethod/PCInomo.pdf")
nom <- nomogram(logit_lrm,
                #fun= funCTion(x)1/(1+exp(-x)), # or 
                fun=plogis,
                #cex.var=3,
                fun.at = c(0.1,seq(0.2,1.0,by=0.2),1.0),
                funlabel="Risk"
                )
plot(nom, cex.var=2.5, cex.axis=2)
dev.off()
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
points.clinical.train<-points_cal(formula = results$formula,rd=data.nom.train.clinical)
auc.train<-auc(y_train,points.clinical.train)
roc.clinical.train<-roc(y_train,points.clinical.train)
points.clinical.test<-points_cal(formula = results$formula,rd=data.nom.test.clinical)

auc.test<-auc(y_test,points.clinical.test)
roc.clinical.test<-roc(y_test,points.clinical.test)

###calculate c index
cindex.train<-rcorrcens(y_train~ predict(logit_lrm), data =  Train.clinic)

###calibration curve
logit_lrm<-lrm(PDL~Radscore_IM+pTNM,x=T,y=T,data = data.nom.train.IM)
logit_lrm<-lrm(PDL~Radscore_IM+pTNM,x=T,y=T,data = data.nom.test.IM)

cal.IM.test <- calibrate(logit_lrm, 
                  cmethod='KM', 
                  method="boot", 
                  u=365, # u需要与之前模型中定义好的time.inc一致，即365或730；
                  m=76, #每次抽样的样本量，
                  B=1000) #抽样次数
library(lattice)
trellis.device(device="windows", height = 20, width = 30, color=TRUE)

plot(1,type ="n",
     xlim =c(0,1),
     ylim =c(0,1),
     xlab = "Nomogram-Predicted Probability of PD-L1 positive expression",
     ylab = "Actual Probability",
     legend = FALSE,
     cex.lab=1.2,
     cex.axis=1.8,#坐标轴字体大小
     subtitles = FALSE)
abline(0,1,col = "black",lty = 2,lwd = 2)
lines(cal.IM.test[,c("predy","calibrated.orig")], type = "l",lwd = 2,col="red",pch =16)
legend(0.7,0.15,   #x,y 定位 或者“bottomright”
       c("Ideal","Bias-corrected"),
       lty = c(2,1),
       lwd = c(2,1),
       cex=1.4,
       col = c("black","red"),
       bty = "n") # "o"为加边框
######Hosmer-Lemeshow test
library(ResourceSelection)
hoslem.test(as.numeric(y_train),attr(cal.IM.train,"predicted"))
hoslem.test(as.numeric(y_test),attr(cal.IM.test,"predicted"))

##############roc curve test
roc.list=list("PCI AUC=0.647"=roc.CT.test,"IM AUC=0.737"=roc.IM.test,"VNC AUC=0.670"=roc.VNC.test,"Clinical AUC=0.564"=roc.clinical.test)
roc.gg<-ggroc(roc.list,size=1,alpha=.9,aes =c("linetype","color")) #size 粗细，alpha透明度 aes按照属性区分
roc.gg+ggsci::scale_color_lancet()+theme(legend.position = c(0.87, 0.15),legend.background = element_rect(fill = "white", color = "black"),legend.text = element_text (size = 15))+
  theme (axis.text=element_text (size=20))+
  theme (axis.title.x = element_text (color = "black", size = 20, face = "bold"), axis.title.y = element_text (color = "black", size = 20, face = "bold"))
##train
roc.list=list("PCI AUC=0.7227"=roc.CT.train,"IM AUC=0.791"=roc.IM.train,"VNC AUC=0.724"=roc.VNC.train,"Clinical AUC=0.687"=roc.clinical.train)
roc.gg<-ggroc(roc.list,size=1.2,alpha=.6,aes =c("linetype","color")) #size 粗细，alpha透明度 aes按照属性区分
roc.gg+ggsci::scale_color_lancet()+theme(legend.position = c(0.87, 0.15),legend.background = element_rect(fill = "white", color = "black"),legend.text = element_text (size = 15))+
  theme (axis.text=element_text (size=20))+
  theme (axis.title.x = element_text (color = "black", size = 20, face = "bold"), axis.title.y = element_text (color = "black", size = 20, face = "bold"))

########compare roc curve
roc.combine.train<-data.frame(PCI=roc.CT.train$original.predictor,IM=roc.IM.train$original.predictor,VNC=roc.VNC.train$original.predictor,Clinical=roc.clinical.train$original.predictor)
roc.combine.test<-data.frame(PCI=roc.CT.test$original.predictor,IM=roc.IM.test$original.predictor,VNC=roc.VNC.test$original.predictor,Clinical=roc.clinical.test$original.predictor)
rocTestResult<- matrix(nrow = ncol(roc.combine.test),ncol = ncol(roc.combine.test))
rocTrainResult<- matrix(nrow = ncol(roc.combine.train),ncol = ncol(roc.combine.train))

for (a in 1:ncol(roc.combine.train)){
  for(b in 1:ncol(roc.combine.train)){
    temp<-roc.test(roc(y_train,roc.combine.train[,a]),roc(y_train,roc.combine.train[,b]),method="delong")
    rocTrainResult[a,b]<-temp$p.value
  }
}
for (a in 1:ncol(roc.combine.test)){
  for(b in 1:ncol(roc.combine.test)){
    temp<-roc.test(roc(y_test,roc.combine.test[,a]),roc(y_test,roc.combine.test[,b]),method="venkatraman") #venkatraman
    rocTestResult[a,b]<-temp$p.value
  }
}
##reorder 函数
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}
rownames(rocTestResult)<-colnames(roc.combine.test)
colnames(rocTestResult)<-colnames(roc.combine.test)
cormat <- reorder_cormat(rocTestResult)
rownames(rocTrainResult)<-colnames(roc.combine.train)
colnames(rocTrainResult)<-colnames(roc.combine.train)
cormat <- reorder_cormat(rocTrainResult)

####
library(RColorBrewer)
library(corrplot)
trellis.device(device="windows", height = 20, width = 30, color=TRUE)
#win.graph(width=4.875, height=2.5,pointsize=8)# to solve the error "figure margin too large"

corrplot::corrplot.mixed(rocTestResult, lower = "number", upper = "circle", number.digits = 3,tl.col = "black",
               lower.col = brewer.pal(5, "Set2"), 
               number.cex=2.5,tl.cex = 2.5,cl.cex = 1.5,
               upper.col = brewer.pal(5, "Set2"))
corrplot::corrplot.mixed(rocTrainResult,number.digits = 3, lower = "number", upper = "circle", tl.col = "black",
                         lower.col = brewer.pal(5, "Set2"), 
                         number.cex=2.5,tl.cex = 2.5,cl.cex = 1.5, 
                        upper.col = brewer.pal(5, "Set2"))

####decision curve
library(rmda)
decisionCurve.test.Clinical<-decision_curve(PDL~Heterogeneity+pTNM,family = binomial(link ='logit'),data = data.nom.test.clinical1,
                                       thresholds = seq(0,1, by = 0.01),
                                      # confidence.intervals= 0.95,
                                       study.design = 'case-control',
                                       population.prevalence= 0.55)
List<- list(decisionCurve.test.Clinical,decisionCurve.test.IM,decisionCurve.test.PCI,decisionCurve.test.VNC)
plot_decision_curve(List,
                    curve.names=c('Clinical','IM','PCI','VNC'),
                    cost.benefit.axis =FALSE,col= c("green",'red','blue',"pink"),
                    confidence.intervals=FALSE,
                    xlim=c(0.24,0.6),
                    standardize = TRUE,
                    cex.axis=1.5,
                    legend.position = "topright"
                    )