options(max.print=10000000)

#Packages
library(ggpubr)
library(lmSupport)
library(lm.beta)
library(ggplot2)
library(ggfortify)
library(dplyr)

##### Load in data ########
dataSet <- read.csv("/Users/JasminePlows/Documents/GitHub/breast_milk_HMOs_and_infant_appetite/Data/hmos_appetite_data.csv",header=TRUE,check.names = FALSE)

##View
# head(dataSet); names(dataSet); dim(dataSet); #View(dataSet)

####Split by timepoint and secretor status####
dataSet1mo <- dataSet[dataSet$timepoint %in% "1",]
dataSet6mo <- dataSet[dataSet$timepoint %in% "6",]


######Add SES score to 6 month dataset
merge_data <- dataSet1mo[,names(dataSet1mo) %in% c("study_id","SES_index_final")]
names(merge_data) <- c("study_id","SES_index_final")
dataSet6mo <- merge(dataSet6mo,merge_data,by="study_id")
names(dataSet6mo)[names(dataSet6mo) == "SES_index_final.y"] <- "SES_index_final"
head(dataSet6mo)

dataSet1moSecretors <- dataSet1mo[dataSet1mo$Secretor %in% 1,]
dataSet1moNonSecretors <- dataSet1mo[dataSet1mo$Secretor %in% 0,]
dataSet6moSecretors <- dataSet6mo[dataSet6mo$Secretor %in% 1,]
dataSet6moNonSecretors <- dataSet6mo[dataSet6mo$Secretor %in% 0,]


# #Make data matrix for PCA analysis
which( colnames(dataSet1moSecretors)=="mom_age_at_birth" )
myvars <- c(1, 3:30) 
dataSet_reduced <- dataSet1moSecretors[,myvars]

clean <- function(x){
  as.numeric(gsub(",","", x))
}

dataSet_reduced[]<- sapply(dataSet_reduced,clean)
dataSet_reduced <- na.omit(dataSet_reduced)
rownames(dataSet_reduced) <- dataSet_reduced[,1]
dataSet_reduced <- dataSet_reduced[-1]
dataSet_reduced
dataSet_reduced$Secretor <- as.factor(dataSet_reduced$Secretor)

pca <- stats::prcomp(dataSet_reduced[-c(1:9)], scale = TRUE, center = TRUE)

# PCA plot
plot <- autoplot(pca, data = dataSet_reduced, colour = "food_responsiveness",
                 main = "1 month secretors only")
         # loadings = TRUE, loadings.colour = 'blue',
         # loadings.label = TRUE, loadings.label.size = 3)

print(plot + labs(colour = "Food Responsiveness"))

# Plotting and extracting PCs
pca <- prcomp(dataSet_reduced[-c(1:9)], scale = TRUE, center = TRUE)
pca

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("1 month PCA Plot")

loading_scores <- pca$rotation[,1]
hmo_scores <- abs(loading_scores)
hmo_scores_ranked <- sort(hmo_scores, decreasing=TRUE)
top_hmos <- names(hmo_scores_ranked[1:19])

top_hmos

pca$rotation[top_hmos,1]

df.projected <- as.data.frame(predict(pca, dataSet_reduced), stringsAsFactors = FALSE)
df.projected$food_responsiveness <- dataSet_reduced$food_responsiveness
df.projected$mom_age_at_birth <- dataSet_reduced$mom_age_at_birth
df.projected$prepreg_bmi_kgm2 <- dataSet_reduced$prepreg_bmi_kgm2
df.projected$baby_gender <- dataSet_reduced$baby_gender
df.projected$birthweight_kg <- dataSet_reduced$birthweight_kg
df.projected$breastmilk_per_day <- dataSet_reduced$breastmilk_per_day
df.projected$mode_of_delivery <- dataSet_reduced$mode_of_delivery
df.projected$SES_index_final <- dataSet_reduced$SES_index_final
ncomp = 2
regvar = paste(paste("PC", 1:ncomp,sep=""), collapse="+")
fmla = paste("food_responsiveness ~", regvar)
mod.com <- glm(fmla, data=df.projected) 
exp(coef(mod.com))
exp(confint(mod.com))

Model_simple <- lm( food_responsiveness ~ mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
                    + birthweight_kg  + mode_of_delivery + breastmilk_per_day + SES_index_final,
                    data=df.projected )
Model_1 <- lm( food_responsiveness~ PC1 + mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
               + birthweight_kg  + mode_of_delivery + breastmilk_per_day + SES_index_final,
               data=df.projected )
Model_2 <- lm( food_responsiveness~ PC2 + mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
               + birthweight_kg  + mode_of_delivery + breastmilk_per_day + SES_index_final,
               data=df.projected )

summary(Model_simple)
summary(Model_1)
lm.beta(Model_1)
modelCompare(Model_simple, Model_1)
anova(Model_simple,Model_1)
summary(Model_2)
lm.beta(Model_2)
modelCompare(Model_simple, Model_2)
anova(Model_simple,Model_2)

##### Start of coding -- select under here and ctrl+f IN SELECTION to replace 15 occurences with the dataSet of interest.
hmos <- c(12:31)
hmoList <- names(dataSet1moSecretors)[hmos]

hmoNamesList <- c("HMO");
allhmoList <- hmoList;

corrList <- numeric(0)
corrpList <- numeric(0)
betaList <- numeric(0)
betaStandardizedList <- numeric(0)
pValList <- numeric(0)
RsqList <- numeric(0)
deltaRsqList <- numeric(0)
FStatList <- numeric(0)
numdfList <- numeric(0)
dendfList <- numeric(0)
pdf(file = paste("/Users/JasminePlows/Documents/GitHub/breast_milk_HMOs_and_infant_appetite/Plots/hmo_foodresponsiveness.pdf",sep=""))
par(mfrow=c(2,2))

for(i in seq(1:length(hmoList)))
{
  print(hmoList[i])
  exposureValues <- as.numeric(gsub(",","",dataSet1moSecretors[,names(dataSet1moSecretors) %in% hmoList[i]]))
  foodResponsiveness <- dataSet1moSecretors$food_responsiveness
  thisDataInstance <- na.omit(data.frame(food_responsiveness=dataSet1moSecretors$food_responsiveness,
                                         hmo=as.numeric(gsub(",","",dataSet1moSecretors[,names(dataSet1moSecretors) %in% hmoList[i]])), 
                                         mom_age_at_birth=dataSet1moSecretors$mom_age_at_birth,
                                         prepreg_bmi_kgm2=dataSet1moSecretors$prepreg_bmi_kgm2,
                                         baby_gender=dataSet1moSecretors$baby_gender,
                                         birthweight_kg=dataSet1moSecretors$birthweight_kg,
                                         breastmilk_per_day=dataSet1moSecretors$breastmilk_per_day,
                                         mode_of_delivery=dataSet1moSecretors$mode_of_delivery,
                                         ses=dataSet1moSecretors$SES_index_final))
  
  #spearman correlations
  corr <- cor.test(formula = ~ food_responsiveness + hmo,
                   data = thisDataInstance)
  #linear models
  Model_simple <- lm( food_responsiveness ~ mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
                      + birthweight_kg  + mode_of_delivery + breastmilk_per_day + ses, 
                      data=thisDataInstance )
  Model_1 <- lm( food_responsiveness~ hmo + mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
                 + birthweight_kg  + mode_of_delivery + breastmilk_per_day + ses, 
                 data=thisDataInstance )
  modelFormForPlot=as.formula(paste("food_responsiveness","~","hmo"));
  statModForPlot <- lm(modelFormForPlot, data=thisDataInstance)
  
  pearcorr <- corr$estimate
  pearp <- corr$p.value 
  beta <- summary(Model_1)$coefficients[2,1]
  betaStandard <- lm.beta(Model_1)$standardized.coefficients[2] #Extract standardized beta
  deltaR2 <- modelCompare(Model_simple, Model_1)$DeltaR2
  pVal <- anova(Model_1,Model_simple)$"Pr(>F)"[2]
  Rsq <- summary(Model_1)$adj.r.squared[1]
  FStat <- anova(Model_1,Model_simple)$"F"[2]
  numdf <- summary(Model_1)$fstatistic[2]
  dendf <- summary(Model_1)$fstatistic[3]
  
  if(pVal <= 0.05)
  {
    title<-paste("Beta = ", format(beta,digits=3),", p-value=", format.pval(anova(Model_1,Model_simple)$"Pr(>F)"[2],digits=3), sep="")
    
    plot(exposureValues, foodResponsiveness,main=title,
         ylab = "Food Responsiveness", 
         xlab = hmoList[i], 
         pch = 19, cex = 1.25, col = "grey39",
         sub = paste("n= ",table(is.na(exposureValues))[1],sep = "")
    )
    abline(statModForPlot, lwd = 2, col="grey39");
  }
  
  corrList[[length(corrList)+1]] <- pearcorr
  corrpList[[length(corrpList)+1]] <- pearp
  pValList[[length(pValList)+1]] <- pVal
  betaList[[length(betaList)+1]] <- beta
  betaStandardizedList[[length(betaStandardizedList)+1]] <- betaStandard
  RsqList[[length(RsqList)+1]] <- Rsq
  deltaRsqList[[length(deltaRsqList)+1]] <- deltaR2
  FStatList[[length(FStatList)+1]] <- FStat
  numdfList[[length(numdfList)+1]] <- numdf
  dendfList[[length(dendfList)+1]] <- dendf
}

adjPvals <- p.adjust(pValList,method = "BH") #Adjusting p-value for multiple comparisons

allhmoList <- cbind(allhmoList, corrList, corrpList, betaList,betaStandardizedList,FStatList, numdfList, dendfList, RsqList, deltaRsqList, pValList, adjPvals);
hmoNamesList <- cbind(hmoNamesList,paste("food_responsiveness","pearson_coefficient",sep="_"),paste("food_responsiveness","pearson_p_value",sep="_"),paste("food_responsiveness","beta",sep="_"),paste("food_responsiveness","standardized_beta",sep="_"),paste("food_responsiveness","f_stat",sep="_"), paste("food_responsiveness","numdf",sep="_"), paste("food_responsiveness","dendf",sep="_"), paste("food_responsiveness","rsq",sep="_"),paste("food_responsiveness","deltarsq",sep="_"),paste("food_responsiveness","pval",sep="_"),paste("food_responsiveness","adj_pval",sep="_"))
graphics.off()

allhmoList<- data.frame(allhmoList)
names(allhmoList) <- hmoNamesList
write.table(allhmoList,paste("/Users/JasminePlows/Documents/GitHub/breast_milk_HMOs_and_infant_appetite/Output/effectSizes_pValues_hmo_foodResponsiveness_dataSet1moSecretors.txt",sep=""),quote=FALSE, sep="\t",append=FALSE, row.names=FALSE, col.names=TRUE)

# FULL MODEL (PCA cluster) at 1month
thisDataInstance <- na.omit(data.frame(food_responsiveness=dataSet1mo$food_responsiveness,
                                       sum=as.numeric(gsub(",","",dataSet1mo$`SUM [ug/mL]`)),
                                       dslnt=as.numeric(gsub(",","",dataSet1mo$`DSLNT [ug/mL]`)),
                                       dflnt=as.numeric(gsub(",","",dataSet1mo$`DFLNT [ug/mL]`)),
                                       mom_age_at_birth=dataSet1mo$mom_age_at_birth,
                                       prepreg_bmi_kgm2=dataSet1mo$prepreg_bmi_kgm2,
                                       baby_gender=dataSet1mo$baby_gender,
                                       birthweight_kg=dataSet1mo$birthweight_kg,
                                       breastmilk_per_day=dataSet1mo$breastmilk_per_day,
                                       mode_of_delivery=dataSet1mo$mode_of_delivery,
                                       ses=dataSet1mo$SES_index_final))

Model_base <- lm( food_responsiveness ~ mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
                    + birthweight_kg  + mode_of_delivery + breastmilk_per_day + ses,
                    data=thisDataInstance )
Model_sum <- lm( food_responsiveness~ sum + mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
               + birthweight_kg  + mode_of_delivery + breastmilk_per_day + ses,
               data=thisDataInstance )
Model_DSLNT <- lm( food_responsiveness~ sum + dslnt + mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
                   + birthweight_kg  + mode_of_delivery + breastmilk_per_day + ses,
                   data=thisDataInstance )
Model_DFLNT <- lm( food_responsiveness~ sum + dslnt + dflnt + mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
                   + birthweight_kg  + mode_of_delivery + breastmilk_per_day + ses,
                   data=thisDataInstance )

summary(Model_base)
summary(Model_sum)
lm.beta(Model_sum)
modelCompare(Model_base, Model_sum)
anova(Model_sum,Model_base)
summary(Model_DSLNT)
lm.beta(Model_DSLNT)
modelCompare(Model_sum, Model_DSLNT)
anova(Model_DSLNT,Model_sum)
summary(Model_DFLNT)
lm.beta(Model_DFLNT)
modelCompare(Model_DSLNT, Model_DFLNT)
anova(Model_DFLNT,Model_DSLNT)

# FULL MODEL (PCA cluster) 6 months
thisDataInstance <- na.omit(data.frame(food_responsiveness=dataSet6mo$food_responsiveness,
                                       dslnh=as.numeric(gsub(",","",dataSet6mo$`DSLNH [ug/mL]`)),
                                       lnh=as.numeric(gsub(",","",dataSet6mo$`LNH [ug/mL]`)),
                                       lstc=as.numeric(gsub(",","",dataSet6mo$`LSTc [ug/mL]`)),
                                       flnh=as.numeric(gsub(",","",dataSet6mo$`FLNH [ug/mL]`)),
                                       mom_age_at_birth=dataSet6mo$mom_age_at_birth,
                                       prepreg_bmi_kgm2=dataSet6mo$prepreg_bmi_kgm2,
                                       baby_gender=dataSet6mo$baby_gender,
                                       birthweight_kg=dataSet6mo$birthweight_kg,
                                       breastmilk_per_day=dataSet6mo$breastmilk_per_day,
                                       mode_of_delivery=dataSet6mo$mode_of_delivery,
                                       ses=dataSet6mo$SES_index_final))

Model_base <- lm( food_responsiveness ~ mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
                  + birthweight_kg  + mode_of_delivery + breastmilk_per_day + ses,
                  data=thisDataInstance )
Model_DSLNH <- lm( food_responsiveness~ dslnh + mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
                 + birthweight_kg  + mode_of_delivery + breastmilk_per_day + ses,
                 data=thisDataInstance )
Model_LNH <- lm( food_responsiveness~ dslnh + lnh + mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
                   + birthweight_kg  + mode_of_delivery + breastmilk_per_day + ses,
                   data=thisDataInstance )
Model_LSTc <- lm( food_responsiveness~ dslnh + lnh + lstc + mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
                   + birthweight_kg  + mode_of_delivery + breastmilk_per_day + ses,
                   data=thisDataInstance )
Model_FLNH <- lm( food_responsiveness~ dslnh + lnh + lstc + flnh + mom_age_at_birth + prepreg_bmi_kgm2 + baby_gender
                   + birthweight_kg  + mode_of_delivery + breastmilk_per_day + ses,
                   data=thisDataInstance )

summary(Model_base)
summary(Model_DSLNH)
lm.beta(Model_DSLNH)
modelCompare(Model_base, Model_DSLNH)
anova(Model_base,Model_DSLNH)
summary(Model_LNH)
lm.beta(Model_LNH)
modelCompare(Model_DSLNH, Model_LNH)
anova(Model_LNH,Model_DSLNH)
summary(Model_LSTc)
lm.beta(Model_LSTc)
modelCompare(Model_LNH, Model_LSTc)
anova(Model_LSTc,Model_LNH)
summary(Model_FLNH)
lm.beta(Model_FLNH)
modelCompare(Model_LSTc, Model_FLNH)
anova(Model_FLNH,Model_LSTc)

#####Correlation matrix
cor.dataSet <- dataSet1mo[c(3:31)]
cor.dataSet[] <- lapply(cor.dataSet, gsub, pattern = ",", replacement ="")
cor.dataSet[] <- sapply(cor.dataSet, as.numeric)
cor.dataSet <- na.omit(cor.dataSet)

cor.matrix_whole <- cor(cor.dataSet)

cor.matrix.data.frame <- as.data.frame(cor.matrix_whole)

cor.data.frame.food.responsiveness.only <- cor.matrix.data.frame[2]

cor.matrix.food.responsiveness.only <- as.matrix(cor.data.frame.food.responsiveness.only)

cor.matrix.food.responsiveness.only

library(corrplot)
corrplot(cor.matrix.food.responsiveness.only, method = "color")



