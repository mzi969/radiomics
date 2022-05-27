rm(list=ls())
library(glmnet)
library(corrplot)  
library(ggplot2)
library(mRMRe)
library(ResourceSelection)
library(RColorBrewer) 
library(rmda)
library(MASS)
getwd()  
setwd("D:/muzi/Desktop/study")  
train <- read.csv("primary_radiomic.csv") 
test <- read.csv("validation_radiomic.csv") 
train <- read.csv("primary_clinical.csv") 
test <- read.csv("validation_clinical.csv") 
train <- train[,-1]
test <- test[,-1]
train_label <- train$group
test_label  <- test$group

#Z-score normalization
standardization <- function(data0, data1 = NULL, data2 = NULL,
                            type = c("minmax", "zscore")){

  if ( class(data1) == "NULL" & class(data2) == "NULL") {
    if( type == "minmax" ){
      
      max_col = apply(data0, 2, function(x) max(x))
      min_col = apply(data0, 2, function(x) min(x)) 
      max_min = max_col - min_col
      
      sub_col0 = sweep(data0, 2, min_col, "-")  
      mer_col0 = sweep(sub_col0, 2, max_min, "/") 
      
      data_result = list(data0 = mer_col0)
    }
    if( type == "zscore" ){
      
      mean_col = apply(data0, 2, function(x) mean(x))
      sd_col   = apply(data0, 2, function(x) sd(x)) 
      
      sub_col0 = sweep(data0, 2, mean_col, "-")
      mer_col0 = sweep(sub_col0, 2, sd_col, "/")
      
      data_result = list(data0 = mer_col0)
      
    }
  }
  
  if ( class(data1) != "NULL" & class(data2) == "NULL" ) {
    if( type == "minmax" ){
      
      max_col = apply(data0, 2, function(x) max(x))
      min_col = apply(data0, 2, function(x) min(x))
      max_min = max_col - min_col
      
      sub_col0 = sweep(data0, 2, min_col, "-")
      sub_col1 = sweep(data1, 2, min_col, "-")
      
      mer_col0 = sweep(sub_col0, 2, max_min, "/")
      mer_col1 = sweep(sub_col1, 2, max_min, "/")
      
      data_result = list(data0 = mer_col0,
                         data1 = mer_col1)
    }
    if( type == "zscore" ){
      
      mean_col = apply(data0, 2, function(x) mean(x))
      sd_col   = apply(data0, 2, function(x) sd(x))
      
      sub_col0 = sweep(data0, 2, mean_col, "-")
      sub_col1 = sweep(data1, 2, mean_col, "-")
      
      mer_col0 = sweep(sub_col0, 2, sd_col, "/")
      mer_col1 = sweep(sub_col1, 2, sd_col, "/")
      
      data_result = list(data0 = mer_col0,
                         data1 = mer_col1)
      
    }
  }
  
  if ( class(data2) != "NULL" ){
    if( type == "minmax" ){
      
      max_col = apply(data0, 2, function(x) max(x))
      min_col = apply(data0, 2, function(x) min(x))
      max_min = max_col - min_col
      
      sub_col0 = sweep(data0, 2, min_col, "-")
      sub_col1 = sweep(data1, 2, min_col, "-")
      sub_col2 = sweep(data2, 2, min_col, "-")
      
      mer_col0 = sweep(sub_col0, 2, max_min, "/")
      mer_col1 = sweep(sub_col1, 2, max_min, "/")
      mer_col2 = sweep(sub_col2, 2, max_min, "/")
      
      data_result = list(data0 = mer_col0,
                         data1 = mer_col1,
                         data2 = mer_col2)
    }
    if( type == "zscore" ){
      
      mean_col = apply(data0, 2, function(x) mean(x))
      sd_col   = apply(data0, 2, function(x) sd(x))
      
      sub_col0 = sweep(data0, 2, mean_col, "-")
      sub_col1 = sweep(data1, 2, mean_col, "-")
      sub_col2 = sweep(data2, 2, mean_col, "-")
      
      mer_col0 = sweep(sub_col0, 2, sd_col, "/")
      mer_col1 = sweep(sub_col1, 2, sd_col, "/")
      mer_col2 = sweep(sub_col2, 2, sd_col, "/")
      
      data_result = list(data0 = mer_col0,
                         data1 = mer_col1,
                         data2 = mer_col2)
      
    }
  }
  
  
  return(data_result)
  }
df_standed <- standardization(data0 = train[,-1], data1 = test[,-1] , type = "zscore")
train2 <- df_standed$data0
test2 <- df_standed$data1

t2<-as.data.frame(train2)
g1<-as.factor(train$group)
t3<-data.frame(g1)
head(t3)
Pvaluekw<-c(rep(0,ncol(t2)))
for(i in 1:ncol(t2))
{
  ab<-as.numeric(t2[1:nrow(t2),i])
  b<-t3$g1
  aa<-data.frame(ab,b)
  y1=kruskal.test(ab~b,data=aa)
  Pvaluekw[i]<-y1$p.value
}
Pvaluekw_p_Index <- which(do.call("rbind", lapply(Pvaluekw, as.data.frame)) < 0.05)
train3 <- t2[,Pvaluekw_p_Index]
dim(train3) 
#Spearman’s correlation
reduce_redundency <- function (dat, threshold = 0.9, method = "spearman") 
{
  if (!("data.frame" %in% class(dat))) {
    stop("Input data must be class data.frame")
  }
  if (sum(is.na(dat)) != 0) {
    stop("Input data contain missing values")
  }
  feature_names = colnames(dat)
  types = sapply(dat, class)
  dataIndex = which(types == "numeric" | types == "integer")
  despIndex = which(types != "numeric" & types != "integer")
  data_numeric = dat[, dataIndex]
  if (length(despIndex) > 0) {
    desp = dat[, despIndex]
  } else {
    desp = NULL
  }
  cor = cor(data_numeric, method = "spearman")
  cor[upper.tri(cor)] = 0
  diag(cor) = 0
  dat.redd = data_numeric[, !apply(cor, 2, function(x) any(abs(x) >
                                                             threshold))]
  features_selected = colnames(dat.redd)
  if (length(despIndex) > 0) {
    dat.redd = data.frame(desp, dat.redd)
  }
  return(list(names = features_selected, dat.redd = dat.redd))
}
train4 <- data.frame(reduce_redundency(train3, threshold = 0.9)$dat.redd)
dim(train4)   
                                  
#mRMR
mrmr_feature<-train4
mrmr_feature$y <-train_label
target_indices = which(names(mrmr_feature)=='y')
for (m in which(sapply(mrmr_feature, class)!="numeric")){
  mrmr_feature[,m]=as.numeric(mrmr_feature[,m])
}

Data <- mRMR.data(data = data.frame(mrmr_feature))
mrmr=mRMR.ensemble(data = Data, target_indices = target_indices, 
                   feature_count = 10, solution_count = 1)

mRME_index=mrmr@filters[[as.character(mrmr@target_indices)]]
train5<- mrmr_feature[,mRME_index]
dim(train5) 

model <- polr(as.factor(train$group)~., data = train5, Hess=TRUE)
tstep <- step(model)
summary(tstep)
drop1(tstep)
View(train5)
train6 <- train2[,c(1,3,7,10,11)]
train6 <- train5[,c(2,3,7,8,9,10)]
library(MASS)
model <- polr(as.factor(train$group) ~., data = train6, Hess=TRUE)
summary(model)
 
#P, OR and 95%CI                                  
(ctable <- coef(summary(model)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))
(ci <- confint(model))
exp(cbind(OR = coef(model), ci))

#brant test                                   
library("brant")
brant(model)

#radscore                                
coefPara <- coef(model)
beta_matrix <- as.matrix(coefPara)
Rad_matrix_train <- as.matrix((train6))
Radscore_train <- Rad_matrix_train %*% beta_matrix
test6 <- test2[,names(train6)]
Rad_matrix_test <- as.matrix((test6))
Radscore_test <- Rad_matrix_test %*% beta_matrix
                                   
#prediction
pred_train_c <- predict(model,newdata=train2,type="class")
(tab <- table(pred_train_c,train$group))
sum(diag(tab))/sum(tab)
pred_train_p <- predict(model,newdata=train2,type="probs")
pred_test_c <- predict(model,newdata=test2,type="class")
(tab <- table(pred_test_c,test$group))
sum(diag(tab))/sum(tab)
pred_test_p <- predict(model,newdata=test2,type="probs")
write.csv(cbind(train6,Radscore_train,pred_train_p,pred_train_c),"result_train.csv")
write.csv(cbind(test6,Radscore_test,pred_test_p,pred_test_c),"result_test.csv")
                                   
#calibration plots
dt_train <- data.frame(Rad_score = Radscore_train) 
ddist_train <- datadist(dt_train)
options(datadist="ddist_train")
lrm_train <- lrm(train_label ~ .,x=T,y=T,data=dt_train)
lrm_train    
summary(lrm_train)  
va.train <- validate(lrm_train,method="boot",B=1000,dxy=T)
cal.train <- calibrate(lrm_train,method="boot",B=1000)
write.csv(cal.train,"cal_train.csv")
Dxy_train = va.train[rownames(va.train)=="Dxy", colnames(va.train)=="index.corrected"]
orig_Dxy_train = va.train[rownames(va.train)=="Dxy", colnames(va.train)=="index.orig"]
bias_corrected_c_index_train  <- abs(Dxy_train)/2+0.5
orig_c_index_train <- abs(orig_Dxy_train)/2+0.5
orig_c_index_train            
bias_corrected_c_index_train  
c_train <- rcorrcens(train_label~predict(lrm_train,newdata=dt_train),data=dt_train) 
dt_test <- data.frame(Rad_score = Radscore_test) 
ddist_test <- datadist(dt_test)
options(datadist="ddist_test")
lrm_test <- lrm(test_label ~ .,data=dt_test,x=T,y=T)
lrm_test    
summary(lrm_test) 
va.test <- validate(lrm_test,method="boot",B=1000,dxy=T) 
cal.test <- calibrate(lrm_test,method="boot",B=1000)
write.csv(cal.test,"cal.test.csv")
Dxy_test = va.test[rownames(va.test)=="Dxy", colnames(va.test)=="index.corrected"]
orig_Dxy_test = va.test[rownames(va.test)=="Dxy", colnames(va.test)=="index.orig"]
bias_corrected_c_index_test <- abs(Dxy_test)/2+0.5
orig_c_index_test <- abs(orig_Dxy_test)/2+0.5
orig_c_index_test
bias_corrected_c_index_test
c_test <- rcorrcens(test_label~predict(lrm_test,newdata=dt_test),data=dt_test) 
tiff(file = "CT_primary_calibration.tiff", res =600, width =4800, height =3600, compression = "lzw")
plot(cal.train,xlim = c(0,1),ylim= c(0,1),scat1d.opts=list(nhistSpike=240,side=1,frac=0),xlab="Nomogram predicted probability of adenocarcinoma")
dev.off()
tiff(file = "CT_validation_calibration.tiff", res =600, width =4800, height =3600, compression = "lzw")
plot(cal.test,xlim = c(0,1),ylim= c(0,1),scat1d.opts=list(nhistSpike=240,side=1,frac=0),xlab="Validation nomogram predicted probability of adenocarcinoma")
dev.off()
tiff(file = "CT_calibration.tiff", res =600, width =4800, height =3600, compression = "lzw")
plot(cal.train,xlab="Predicted probability",
     ylab="Actual probability",
     xlim=c(0,1),ylim=c(0,1),legend,
     lines(cal.train,lty=1,lwd=2,type="o",pch=20,cex=0.7,col=my.col7[2])
     +lines(cal.test,lty=1,lwd=2,type="o",pch=20,cex=0.7,col=my.col7[5])
     +abline(0,1,lty=1,lwd=2,col="grey")
     +legend(x=0.75,y=0.2,col=c(my.col7[2],my.col7[5]),
             legend=c("Primary","Validation"),
             lty=c(1,1),bty="n"))
dev.off()


