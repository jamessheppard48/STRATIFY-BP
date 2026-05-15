###########################################################
##  STRATIFY - CVD risk                                  ##
##  This programme is used for the final CVDrisk model   ##
###########################################################

## Packages install
install.packages("mice")
install.packages("dplyr")
install.packages("VIM")
install.packages("naniar")
install.packages("questionr")
install.packages("data.table")
install.packages("plyr")
install.packages("flextable")
install.packages("sjPlot")
install.packages("foreach")
install.packages("dynpred")
install.packages("fastcmprsk")
install.packages("foreign")
install.packages("tidyr")
install.packages("matrixStats")
install.packages("mfp")
install.packages("visreg")
install.packages("varhandle")

library(mice)
library(plyr)
library(dplyr)
library(VIM)
library(naniar)
library(questionr)
library(data.table)
library(flextable)
library(sjPlot)
library(foreach)
library(dynpred)
library(fastcmprsk)
library(foreign)
library(tidyr)
library(matrixStats)
library(mfp)
library(visreg)
library(varhandle)

## Competing risk model (Fine & Gray)

##############################
## Multivariate
##############################

## Linearity in multivariate models 
## FP1 transformation for Age selected (based on cox model)
## Reselect transformation after imputation

#####################################################################
#####################################################################

#### Imputation procedure for CVD model

# Variables requiring imputation
# 1) BMI (ordered logistic regression)  
# 2) Ethnicity (multinomial logistic regression)
# 3) Smoking (multinomial logistic regression)
# 4) TC_HDL_value (Predictive mean matching)

CVD_variables<-select(mydata, BMIgroups, Ethnicity_recoded, Smoking, TC_HDL_value,
                      CVDCox, CVD_comp, TimetoEvent, age, Gender,
                      patid, imd2015_5, sysbp_inclusion, TreatedHypertension,
                      RheumatoidArthritis, Diabetes, Qrisk2FH, CKD, Atrialfib, 
                      prevStroke, prevMI, HF, statins )


# make sure categorical variables are registered as factor

CVD_variables$Smoking<-as.factor(CVD_variables$Smoking)
CVD_variables$imd2015_5<-as.factor(CVD_variables$imd2015_5)
CVD_variables$Ethnicity_recoded<-as.factor(CVD_variables$Ethnicity_recoded)
CVD_variables$BMIgroups<-as.factor(CVD_variables$BMIgroups)
CVD_variables$Gender<-as.factor(CVD_variables$Gender)

## Estimate the NA estimator for main event and competing event

CVD_variables$Competing<-CVD_variables$CVD_comp
CVD_variables$Competing[CVD_variables$CVD_comp==2] <- 1
CVD_variables$Competing[CVD_variables$CVD_comp==1] <- 0

## estimate the NA estimators for the main event and competing risks
CVD_variables$nelsonalen_main<-nelsonaalen(CVD_variables, TimetoEvent, CVDCox )
CVD_variables$nelsonalen_competing<-nelsonaalen(CVD_variables, TimetoEvent, Competing)

# Start imputation by generating the predictor matrix
ini <- mice(CVD_variables, maxit = 0)
mat<-ini$predictorMatrix

mat[,"nelsonalen_competing"]<-1 # set NA comp to 1 in order to be used as predictor for imp model
mat["nelsonalen_competing",]<-0 # set NA comp to 0 in order not be considererd for imputation
mat["TimetoEvent",]<-0  # set TimetoEvent to 0
mat[,"TimetoEvent"]<-0  # set TimetoEvent to 0

# set patid variable to 0 for both rows and columns as we don't want to use it all in imputation
mat["patid",]<-0
mat[,"patid"]<-0

# set the variables with no missing data to have 0s as rows as we don't want to impute them

mat["CVDCox" ,]<-0
mat["CVD_comp" ,]<-0
mat[,"CVD_comp"]<-0 # not to be included as predictor due to collinearity 
mat["age",]<-0
mat["Gender" ,]<-0
mat["imd2015_5" ,]<-0
mat["sysbp_inclusion" ,]<-0
mat["TreatedHypertension",]<-0
mat["RheumatoidArthritis",]<-0
mat["Diabetes",]<-0
mat["Qrisk2FH",]<-0
mat["CKD",]<-0
mat["Atrialfib",]<-0
mat["prevStroke",]<-0
mat["prevMI",]<-0
mat["HF",]<-0
mat["statins",]<-0
mat["nelsonalen_main",]<-0 
mat["Competing",]<-0 

imp <- mice(CVD_variables, m=10, pred=mat, method=c("polr","polyreg", "polyreg", "pmm", 
                                                    "", "", "", "", "", "", "", "", "", "",
                                                    "", "", "", "", "", "", "", "", "", "", ""), 
            seed=123, maxit=10)

com3<-mice::complete(imp, "broad", include=TRUE)  


imp1$age_1<-(imp1$age/100)^3
summary(imp1$age_1)
sum(is.na(imp1$age_1))
imp1$age_1<-(imp1$age_1-0.24188)
summary(imp1$age_1)

### First imputation ####
fit1<-fastCrr(Crisk(imp1$TimetoEvent, imp1$CVD_comp)~
                imp1$age_1 + 
                as.factor(imp1$Gender) +
                as.factor(imp1$imd2015_5) +
                relevel(factor(imp1$BMIgroups), ref = "1") +
                as.factor(imp1$Ethnicity_recoded)+
                as.factor(imp1$Smoking)+
                imp1$TC_HDL_value+
                imp1$sysbp_inclusion+
                imp1$TreatedHypertension+
                imp1$RheumatoidArthritis+
                imp1$Diabetes+
                imp1$Qrisk2FH+                   
                imp1$CKD+
                imp1$Atrialfib+
                imp1$prevStroke+
                imp1$prevMI+                   
                imp1$HF+
                imp1$statins,  variance =TRUE, max.iter = 1000, 
              var.control = varianceControl(B = 100))
summary(fit1)

### Repeat 9 imputations ####
#### Combine coefficients
models <- mget(paste0("fit", 1:10))
m <- length(models)
coef_list <- lapply(models, coef)
var_list  <- lapply(models, function(x) diag(vcov(x)))
all(sapply(coef_list, function(x) length(x)) == length(coef_list[[1]]))
coef_mat <- do.call(cbind, coef_list)  
var_mat  <- do.call(cbind, var_list)

Beta <- rowMeans(coef_mat)
Covar <- rowMeans(var_mat)
B <- apply(coef_mat, 1, var)
Varian <- Covar + (1 + 1/m) * B
SE_pooled <- sqrt(Varian)

lower <- Beta - 1.96 * SE_pooled
upper <- Beta + 1.96 * SE_pooled
SHR  <- exp(Beta)
LCI <- exp(lower)
UCI <- exp(upper)

Final_estimates <- data.frame(Beta = Beta, SHR = SHR, LowerCI = LCI, UpperCI = UCI)

rownames(Final_estimates) <- names(coef_list[[1]])
options(scipen = 999)
Final_estimates

## Baseline CIF 

cov2 <- rep(0, 27)

CIF.hat1 <- cumsum(exp(sum(cov2 * fit1$coef)) * fit1$breslowJump[, 2]) 
test<-cbind(fit1$breslowJump, CIF.hat1)
imp1$time<-imp1$TimetoEvent
test <- merge(imp1 , test , by="time", all= TRUE)
test<-test %>% fill(CIF.hat1)
test$CIF.hat1[1]<-0
test$CIF.hat2 <- 1 - exp(-test$CIF.hat1)
imp1_bs<-cbind(test$time, test$CIF.hat2)
colnames(imp1_bs)<-c("time1", "CIFimp1")
sum(is.na(imp1_bs))

### Repeat 9 imputations ####
cif_matrix<-cbind(imp1_bs, imp2_bs, imp3_bs, imp4_bs, imp5_bs, imp6_bs, imp7_bs, imp8_bs,
                  imp9_bs, imp10_bs)
all(cif_matrix[,1] == cif_matrix[,3])
all(cif_matrix[,3] == cif_matrix[,5])
all(cif_matrix[,5] == cif_matrix[,7])
all(cif_matrix[,7] == cif_matrix[,9])
all(cif_matrix[,9] == cif_matrix[,11])
all(cif_matrix[,11] == cif_matrix[,13])
all(cif_matrix[,13] == cif_matrix[,15])
all(cif_matrix[,15] == cif_matrix[,17])
all(cif_matrix[,17] == cif_matrix[,19])
cif_matrix2<-cbind(imp1_bs, imp2_bs[,2], imp3_bs[,2], imp4_bs[,2], 
                   imp5_bs[,2], imp6_bs[,2], imp7_bs[,2], imp8_bs[,2],
                   imp9_bs[,2], imp10_bs[,2])
colnames(cif_matrix2)<-c("time", "CIF1", "CIF2", "CIF3", "CIF4", "CIF5",
                         "CIF6", "CIF7", "CIF8", "CIF9", "CIF10")
cif_matrix2<-as.data.frame(cif_matrix2)
CIF3650<- cif_matrix2[which(cif_matrix2$time==3650), ]
CIF3650<-c(mean(CIF3650$CIF1), mean(CIF3650$CIF2), mean(CIF3650$CIF3), mean(CIF3650$CIF4),
           mean(CIF3650$CIF5), mean(CIF3650$CIF6), mean(CIF3650$CIF7), mean(CIF3650$CIF8),
           mean(CIF3650$CIF9), mean(CIF3650$CIF10))
loglogt<-log(-log(1-CIF3650))
mean(loglogt)
backtransf<-1-exp(-exp(mean(loglogt)))
CIF3650<-cbind(CIF3650,median(CIF3650), mean(CIF3650), backtransf)
colnames(CIF3650)<-c("CIF10y", "median CIF", "mean CIF", "log log")


CIF1825<- cif_matrix2[which(cif_matrix2$time==1825), ]
CIF1825<-c(mean(CIF1825$CIF1), mean(CIF1825$CIF2), mean(CIF1825$CIF3), mean(CIF1825$CIF4),
           mean(CIF1825$CIF5), mean(CIF1825$CIF6), mean(CIF1825$CIF7), mean(CIF1825$CIF8),
           mean(CIF1825$CIF9), mean(CIF1825$CIF10))
loglogt<-log(-log(1-CIF1825))
mean(loglogt)
backtransf<-1-exp(-exp(mean(loglogt)))
CIF1825<-cbind(CIF1825,median(CIF1825), mean(CIF1825), backtransf)
colnames(CIF1825)<-c("CIF5y", "median CIF", "mean CIF", "log log")


CIF365<- cif_matrix2[which(cif_matrix2$time==365), ]
CIF365<-c(mean(CIF365$CIF1), mean(CIF365$CIF2), mean(CIF365$CIF3), mean(CIF365$CIF4),
          mean(CIF365$CIF5), mean(CIF365$CIF6), mean(CIF365$CIF7), mean(CIF365$CIF8),
          mean(CIF365$CIF9), mean(CIF365$CIF10))
loglogt<-log(-log(1-CIF365))
mean(loglogt)
backtransf<-1-exp(-exp(mean(loglogt)))
CIF365<-cbind(CIF365,median(CIF365), mean(CIF365), backtransf)
colnames(CIF365)<-c("CIF1y", "median CIF", "mean CIF", "log log")

######## predictions ########
##### Imputation 1 ##### 

imp1$Underweight[imp1$BMIgroups==0]<-1
imp1$Underweight[imp1$BMIgroups!=0]<-0
imp1$Normal[imp1$BMIgroups==1]<-1
imp1$Normal[imp1$BMIgroups!=1]<-0
imp1$Overweight[imp1$BMIgroups==2]<-1
imp1$Overweight[imp1$BMIgroups!=2]<-0
imp1$Obese[imp1$BMIgroups==3]<-1
imp1$Obese[imp1$BMIgroups!=3]<-0
imp1$Mobese[imp1$BMIgroups==4]<-1
imp1$Mobese[imp1$BMIgroups!=4]<-0

imp1$White[imp1$Ethnicity_recoded==0]<-1
imp1$White[imp1$Ethnicity_recoded!=0]<-0
imp1$Black[imp1$Ethnicity_recoded==1]<-1
imp1$Black[imp1$Ethnicity_recoded!=1]<-0
imp1$Asian[imp1$Ethnicity_recoded==2]<-1
imp1$Asian[imp1$Ethnicity_recoded!=2]<-0
imp1$Other[imp1$Ethnicity_recoded==3]<-1
imp1$Other[imp1$Ethnicity_recoded!=3]<-0

imp1$Non_smoker[imp1$Smoking==0]<-1
imp1$Non_smoker[imp1$Smoking!=0]<-0
imp1$Ex_smoker[imp1$Smoking==1]<-1
imp1$Ex_smoker[imp1$Smoking!=1]<-0
imp1$Smoker[imp1$Smoking==2]<-1
imp1$Smoker[imp1$Smoking!=2]<-0

imp1$IMD1[imp1$imd2015_5==0]<-1
imp1$IMD1[imp1$imd2015_5!=0]<-0
imp1$IMD2[imp1$imd2015_5==1]<-1
imp1$IMD2[imp1$imd2015_5!=1]<-0
imp1$IMD3[imp1$imd2015_5==2]<-1
imp1$IMD3[imp1$imd2015_5!=2]<-0
imp1$IMD4[imp1$imd2015_5==3]<-1
imp1$IMD4[imp1$imd2015_5!=3]<-0
imp1$IMD5[imp1$imd2015_5==4]<-1
imp1$IMD5[imp1$imd2015_5!=4]<-0

imp1$Female[imp1$Gender==0] <- 0
imp1$Female[imp1$Gender==1] <- 1

PI1<-((imp1$age_1) * (Final_estimates[1,1])) + 
  (imp1$Female * (Final_estimates[2,1])) +
  (imp1$IMD2 *(Final_estimates[3,1])) +
  (imp1$IMD3 *(Final_estimates[4,1])) +
  (imp1$IMD4 * (Final_estimates[5,1])) + 
  (imp1$IMD5 * (Final_estimates[6,1])) + 
  (imp1$Underweight * (Final_estimates[7,1])) +
  (imp1$Overweight * (Final_estimates[8,1])) + 
  (imp1$Obese * (Final_estimates[9,1])) +
  (imp1$Mobese * (Final_estimates[10,1])) +
  (imp1$Black * (Final_estimates[11,1])) +  
  (imp1$Asian * (Final_estimates[12,1])) +
  (imp1$Other * (Final_estimates[13,1])) + 
  (imp1$Ex_smoker * (Final_estimates[14,1])) +
  (imp1$Smoker * (Final_estimates[15,1])) + 
  (imp1$TC_HDL_value * (Final_estimates[16,1])) +
  (imp1$sysbp_inclusion * (Final_estimates[17,1])) +
  (imp1$TreatedHypertension * (Final_estimates[18,1])) +
  (imp1$RheumatoidArthritis * (Final_estimates[19,1])) +
  (imp1$Diabetes * (Final_estimates[20,1])) +
  (imp1$Qrisk2FH * (Final_estimates[21,1])) +
  (imp1$CKD *  (Final_estimates[22,1])) +
  (imp1$Atrialfib *  (Final_estimates[23,1])) +
  (imp1$prevStroke * (Final_estimates[24,1])) +
  (imp1$prevMI * (Final_estimates[25,1])) +
  (imp1$HF *  (Final_estimates[26,1])) +
  (imp1$statins *  (Final_estimates[27,1])) 

### Repeat 9 imputations ####
#### combine and save

PIcomb<-as.data.frame(cbind(imp1$TimetoEvent, imp1$CVD_comp, PI1, PI2, PI3, PI4, PI5, PI6,
                            PI7, PI8, PI9, PI10))

PIcomb$combPI<-(PI1+PI2+PI3+PI4+PI5+PI6+PI7+PI8+PI9+PI10)/10

risk10y<-1-(((1- 0.0146165))^exp(PIcomb$combPI))
risk5y<-1-(((1- 0.00699506))^exp(PIcomb$combPI))
risk1y<-1-(((1- 0.001500368))^exp(PIcomb$combPI))


PIcomb<-as.data.frame(cbind(PIcomb, risk1y, risk5y, risk10y, mydata$patid))


