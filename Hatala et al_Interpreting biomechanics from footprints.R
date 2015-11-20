## Article citation information:
## Interpreting locomotor biomechanics from the morphology of human footprints
## Kevin G. Hatala, Roshna E. Wunderlich, Heather L. Dingwall, Brian G. Richmond
## Journal of Human Evolution 90, 38-48 (2016).



## Load data and required packages

## Load table containing measurements of Daasanach footprints with associated kinematic and kinetic data for each trial
daasData<-read.csv(file.choose())
str(daasData)
## For information regarding this dataset, contact me at kevin.g.hatala@gmail.com

## Load packages 'dplyr','MASS','nlme','ggplot2','grid','gridExtra'
library(dplyr)
library(MASS)
library(nlme)
library(ggplot2)
library(grid)
library(gridExtra)
## Load functions from Zuur et al. 2009 AED package (http://www.highstat.com/book2.htm) - used to calculate variance inflation factors for mixed effects models (corvif function)
source("~/HighstatLibV6.R")



##Analysis - This follows (more or less) the sequence of Section 2.5 in JHE publication

## Use PCA to identify principal axes of variation in human footprint topography
## Columns 5-18 in Daasanach data set are normalized regional depth measurements from footprints produced in each trial
daasPrintPCA<-prcomp(daasData[,5:18],retx=TRUE)
summary(daasPrintPCA)

## Create scree plot to show PCA results (Figure 4)
## Calculate proportion of variance explained by first 9 PC axes (together explain over 95% of total variance)
propVar<-as.data.frame((daasPrintPCA$sdev^2/sum(daasPrintPCA$sdev^2))[1:9])
## Format data for desired use in ggplot2
screeData<-data.frame(row.names(propVar),propVar)
colnames(screeData)<-c("Component","propVar")
screeData[,1]<-c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9")
screeData$Component<-factor(screeData$Component)
## Draw scree plot
screePlot<-ggplot(screeData,aes(x=Component,y=propVar))+
  geom_bar(stat="identity",colour="black")+
  scale_y_continuous(breaks=round(seq(0,0.40,by = 0.05),2))+
  geom_text(aes(label=round(propVar,digits=3)),vjust=1.15,colour="white",size=3)+
  labs(y="Proportion of total variance", x="Component")+
  theme(axis.text=element_text(size=10,colour="black"),axis.title=element_text(size=12,face="bold",vjust=0.2))
print(screePlot)
pdf(file="screePlot.pdf")
print(screePlot)
dev.off()

## Plot and examine rotations for the first five PC axes - these all individually describe at least 5% of total variance (>85% together)
footRegions<-c("Lateral heel","Medial heel","Lateral midfoot","Medial midfoot","Metatarsal 1","Metatarsal 2","Metatarsal 3","Metatarsal 4","Metatarsal 5","Hallux","Toe 2","Toe 3","Toe 4","Toe 5")
pcaRotations<-as.data.frame(daasPrintPCA$rotation)
rownames(pcaRotations)<-c(1:14)
pcaRotations<-cbind(footRegions,pcaRotations)
colnames(pcaRotations)[1]<-"Region"
str(pcaRotations)
## Create logical columns assessing whether rotations on each of PCs 1-5 are positive or negative - these are used to color positive and negative loadings differently in the plot 
pcaRotations$pc1Pos<-pcaRotations$PC1>=0
pcaRotations$pc2Pos<-pcaRotations$PC2>=0
pcaRotations$pc3Pos<-pcaRotations$PC3>=0
pcaRotations$pc4Pos<-pcaRotations$PC4>=0
pcaRotations$pc5Pos<-pcaRotations$PC5>=0
## Draw a separate graph for loadings on each of PCs 1-5 (Figure 5)
pc1Plot<-ggplot(pcaRotations,aes(x=Region,y=PC1,fill=pc1Pos))+
  geom_bar(stat="identity",position="identity",colour="black",size=0.25)+
  scale_y_continuous(breaks=round(seq(min(pcaRotations$PC1),max(pcaRotations$PC1),by = 0.2),1))+
  scale_fill_manual(values=c("#CCCCCC","#333333"),guide=FALSE)+
  labs(y="PC1 rotation")+
  theme(axis.text.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),plot.margin=unit(c(0.05,0.5,0.05,0.1),"cm"))
print(pc1Plot)
pc2Plot<-ggplot(pcaRotations,aes(x=Region,y=PC2,fill=pc2Pos))+
  geom_bar(stat="identity",position="identity",colour="black",size=0.25)+
  scale_y_continuous(breaks=round(seq(min(pcaRotations$PC2),max(pcaRotations$PC2),by = 0.2),1))+
  scale_fill_manual(values=c("#CCCCCC","#333333"),guide=FALSE)+
  labs(y="PC2 rotation")+
  theme(axis.text.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),plot.margin=unit(c(0.05,0.5,0.05,0.1),"cm"))
print(pc2Plot)
pc3Plot<-ggplot(pcaRotations,aes(x=Region,y=PC3,fill=pc3Pos))+
  geom_bar(stat="identity",position="identity",colour="black",size=0.25)+
  scale_y_continuous(breaks=round(seq(min(pcaRotations$PC3),max(pcaRotations$PC3),by = 0.2),1))+
  scale_fill_manual(values=c("#CCCCCC","#333333"),guide=FALSE)+
  labs(y="PC3 rotation")+
  theme(axis.text.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),plot.margin=unit(c(0.05,0.5,0.05,0.1),"cm"))
print(pc3Plot)
pc4Plot<-ggplot(pcaRotations,aes(x=Region,y=PC4,fill=pc4Pos))+
  geom_bar(stat="identity",position="identity",colour="black",size=0.25)+
  scale_y_continuous(breaks=round(seq(min(pcaRotations$PC4),max(pcaRotations$PC4),by = 0.2),1))+
  scale_fill_manual(values=c("#CCCCCC","#333333"),guide=FALSE)+
  labs(y="PC4 rotation")+
  theme(axis.text.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),plot.margin=unit(c(0.05,0.5,0.05,0.1),"cm"))
print(pc4Plot)
pc5Plot<-ggplot(pcaRotations,aes(x=Region,y=PC5,fill=pc5Pos))+
  geom_bar(stat="identity",position="identity",colour="black",size=0.25)+
  scale_y_continuous(breaks=round(seq(min(pcaRotations$PC5),max(pcaRotations$PC5),by = 0.2),1))+
  scale_fill_manual(values=c("#CCCCCC","#333333"),guide=FALSE)+
  labs(y="PC5 rotation")+
  theme(axis.text.y=element_text(size=10),axis.text.x=element_text(angle=90,vjust=0.1,size=10),plot.margin=unit(c(0.05,0.5,0.05,0.1),"cm"))
print(pc5Plot)
## Create a stacked grid combining all 5 PC plots - only PC 5 plot, at the bottom, has a labeled x-axis
pdf(file="pcPlot.pdf")
pcPlot<-grid.arrange(pc1Plot,pc2Plot,pc3Plot,pc4Plot,pc5Plot,nrow=5,ncol=1,heights=c(1,1,1,1,1.75))
dev.off()

## Extract PC scores for first five PC axes (those which individually describe at least 5% of total topographic variance)
daasPCScores<-daasPrintPCA$x
daasPCs<-daasPCScores[,1:5]
## Bind PC scores to original data set
daasData<-cbind(daasData,daasPCs)

## Linear mixed effects modeling - analyzing effects of biomechanical variables on PC axes of footprint shape variation
## What is the shape change described by PC1?
pcaRotations[,c(1,2)]
## PC1 is described by deep heel vs deep forefoot
## Biomechanical hypothesis: Is this influenced by relative forefoot versus heel pressure? Total movement of COM (hip excursion)? Gait type?
## Does substrate compliance (overall depth) also have a significant effect?
## Collect subset of data to test these hypotheses (observations where all variables were measured)
pc1Data<-daasData[complete.cases(daasData[,c("PC1","Gait","Max.forefoot.heel.P","Mean.depth","First.half.hip.excursion","Second.half.hip.excursion")]),]
## Check variance inflation factors of continuous predictor variables (test for collinearity)
corvif(daasData[,c("Max.forefoot.heel.P","Mean.depth","First.half.hip.excursion","Second.half.hip.excursion")])
## Gait type and second half hip excursion show evidence of some collinearity (VIFs of about 4) when both are included in this test - potential problem if both are included in model
## When either is excluded, all VIFs are less than 3 (suggests no collinearity among other predictors)
## Create linear mixed effects model of PC1 and various hypothesized predictors - use random intercept model with Subject as random effect
model1PC1<-lme(PC1~Gait+Max.forefoot.heel.P+Mean.depth+First.half.hip.excursion+Second.half.hip.excursion,random=~1|Subject,data=pc1Data,na.action=na.omit,method="ML")
summary(model1PC1)
## First half hip extension has least significant effect - remove from model
model2PC1<-update(model1PC1,.~.-First.half.hip.excursion)
anova(model1PC1,model2PC1)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model.
summary(model2PC1)
## Second half hip extension has least significant effect - remove from model
model3PC1<-update(model2PC1,.~.-Second.half.hip.excursion)
anova(model2PC1,model3PC1)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model.
summary(model3PC1)
## Difference between max forefoot and heel pressures also has a nonsignificant effect - remove from model.
model4PC1<-update(model3PC1,.~.-Max.forefoot.heel.P)
anova(model3PC1,model4PC1)
## Removal does not significantly affect model fit - keep out of model. 
## Note: In this case, AIC was very slightly increased with removal (from 252.7016 to 252.7427) but the scale of this increase is incredibly small. Variable was still removed for model simplification because its effect was non-significant.
summary(model4PC1)
## All remaining effects are significant
## Check structure of residuals
plot(pc1Data$Mean.depth,model4PC1$residuals[,1])
plot(pc1Data$Gait,model4PC1$residuals[,1])
plot(model4PC1)
plot(density(residuals(model4PC1)))
qqnorm(model4PC1)
## Extract Pearson (standardized, divided by standard errors) residuals
PC1resids<-residuals(model4PC1,type="pearson")
mean(PC1resids^2)
## Plot example of gait effects for one subject (Figure 6)
dcData<-subset(daasData,Subject=="DC")
dcPC1Plot<-ggplot(dcData,aes(x=Gait,y=PC1))+
  geom_boxplot()+
  labs(y="PC1")+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(size=8,color="black"),panel.background=element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),
        axis.text.y=element_text(size=8,color="black"),axis.title.y=element_text(size=10,color="black"),axis.ticks=element_line(colour="black"))
print(dcPC1Plot)
jpeg(filename="dcPC1Plot.jpg",width=3,height=3,units="in",quality=100,res=300)
dcPC1Plot
dev.off()

## What is the shape change described by PC2?
pcaRotations[,c(1,3)]
## PC2 described by relatively deep hallux and 2nd toe impressions
## Biomechanical hypothesis: Is this influenced by toe-off mechanism (peak pressures on toes)? Or hip retraction? Gait type?
## Does substrate compliance (overall depth) also have a significant effect?
## Collect subset of data to test these hypotheses (observations where all variables were measured)
pc2Data<-daasData[complete.cases(daasData[,c("PC2","Gait","Hallux.peak.P","Lat.toes.peak.P","Second.half.hip.excursion","Mean.depth")]),]
## Check variance inflation factors of continuous predictor variables (test for collinearity)
corvif(daasData[,c("Gait","Hallux.peak.P","Lat.toes.peak.P","Second.half.hip.excursion","Mean.depth")])
## Once again, if gait type is included, it appears somewhat collinear with second half hip excursion
## If either is removed all VIFs are less than 3, suggesting no collinearity among other predictors
## Create linear mixed effects model of PC2 and various hypothesized predictors - use random intercept model with Subject as random effect
model1PC2<-lme(PC2~Gait+Hallux.peak.P+Lat.toes.peak.P+Second.half.hip.excursion+Mean.depth,random=~1|Subject,data=pc2Data,na.action=na.omit,method="ML")
summary(model1PC2)
## Gait type has least significant effect - remove from model
model2PC2<-update(model1PC2,.~.-Gait)
anova(model1PC2,model2PC2)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model.
summary(model2PC2)
## Mean depth has least significant effect - remove from model
model3PC2<-update(model2PC2,.~.-Mean.depth)
anova(model2PC2,model3PC2)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model.
summary(model3PC2)
## Hallux peak P has least significant effect - remove from model
model4PC2<-update(model3PC2,.~.-Hallux.peak.P)
anova(model3PC2,model4PC2)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model.
summary(model4PC2)
## Lateral toes peak P has non-significant effect - remove from model
model5PC2<-update(model4PC2,.~.-Lat.toes.peak.P)
anova(model4PC2,model5PC2)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model.
summary(model5PC2)
## Check structure of residuals
plot(pc2Data$Second.half.hip.excursion,model5PC2$residuals[,1])
plot(model5PC2)
plot(density(residuals(model5PC2)))
qqnorm(model5PC2)
## Extract Pearson (standardized, divided by standard errors) residuals
pc2ModelResiduals<-residuals(model5PC2,type="pearson")
mean(pc2ModelResiduals^2)
## Plot example of biomechanical relationship (PC2 to hip excursion) for one subject (see Figure 7)
bkData<-subset(daasData,Subject=="BK")
bkPC2Data<-as.matrix(cbind(bkData$Second.half.hip.excursion,bkData$PC2))
bkPC2Data<-as.data.frame(bkPC2Data)
bkPC2Plot<-ggplot(bkPC2Data,aes(x=bkPC2Data[,1],y=bkPC2Data[,2]))+
  geom_point(shape=19,size=2,alpha=.9)+
  labs(x="Hip retraction (radians)",y="PC2")+
  theme(axis.text.x=element_text(size=10,colour="black"),axis.text.y=element_text(size=10,colour="black"),axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),legend.position=c(1,0),legend.justification=c(1,0),
        legend.background=element_rect(colour="black",fill="white"),panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),axis.ticks=element_line(colour="black"))
print(bkPC2Plot)
jpeg(file="bkPC2Plot.jpg",width=3.25,height=4.25,units="in",quality=100,res=300)
print(bkPC2Plot)
dev.off()

## What is the shape change described by PC3?
pcaRotations[,c(1,4)]
## PC3 described by deep heel and toes relative to metatarsal heads and midfoot 
## Biomechanical hypotheses: Related to amount of hip or ankle flexion/extension? Gait type?
## Substrate effect?
## Collect subset of data to test these hypotheses (observations where all variables were measured)
pc3Data<-daasData[complete.cases(daasData[,c("PC3","Gait","Footprint.strike.hip","Footprint.toeoff.hip","Footprint.strike.ankle","Footprint.toeoff.ankle","Mean.depth")]),]
## Check variance inflation factors of continuous predictor variables (test for collinearity)
corvif(daasData[,c("Gait","Footprint.strike.hip","Footprint.toeoff.hip","Footprint.strike.ankle","Footprint.toeoff.ankle","Mean.depth")])
## All VIFs are less than 3 (suggests no collinearity).
## Create linear mixed effects model of PC3 and various hypothesized predictors - use random intercept model with Subject as random effect
model1PC3<-lme(PC3~Gait+Footprint.strike.hip+Footprint.toeoff.hip+Footprint.strike.ankle+Footprint.toeoff.ankle+Mean.depth,random=~1|Subject,data=pc3Data,na.action=na.omit,method="ML")
summary(model1PC3)
## Hip angle at toe-off has least significant effect - remove from model
model2PC3<-update(model1PC3,.~.-Footprint.toeoff.hip)
anova(model1PC3,model2PC3)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model
summary(model2PC3)
## Gait type has least significant effect - remove from model
model3PC3<-update(model2PC3,.~.-Gait)
anova(model2PC3,model3PC3)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model
summary(model3PC3)
## Ankle angle at toe-off has least significant effect - remove from model
model4PC3<-update(model3PC3,.~.-Footprint.toeoff.ankle)
anova(model3PC3,model4PC3)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model
summary(model4PC3)
## Ankle angle at strike has least significant effect - remove from model
model5PC3<-update(model4PC3,.~.-Footprint.strike.ankle)
anova(model4PC3,model5PC3)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model
summary(model5PC3)
## Hip angle at strike has least significant effect - remove from model
model6PC3<-update(model5PC3,.~.-Footprint.strike.hip)
anova(model5PC3,model6PC3)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model
summary(model6PC3)
## Check structure of residuals
plot(pc3Data$Mean.depth,model6PC3$residuals[,1])
plot(model6PC3)
plot(density(residuals(model6PC3)))
qqnorm(model6PC3)
## Plot example of relationship (PC3 to mean depth) for one subject (see Figure 7)
dhData<-subset(daasData,Subject=="DH")
dhPC3Plot<-ggplot(dhData,aes(x=Mean.depth,y=PC3))+
  geom_point(shape=19,size=2,alpha=.9)+
  labs(x="Substrate deformation (cm)",y="PC3")+
  theme(axis.text.x=element_text(size=10,colour="black"),axis.text.y=element_text(size=10,colour="black"),axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),legend.position=c(1,0),legend.justification=c(1,0),legend.background=element_rect(colour="black",
        fill="white"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),axis.ticks=element_line(colour="black"))
print(dhPC3Plot)
jpeg(file="dhPC3Plot.jpg",width=3.25,height=4.25,units="in",quality=100,res=300)
print(dhPC3Plot)
dev.off()

## What is the shape change described by PC4?
pcaRotations[,c(1,5)]
## PC4 described by deep toes 4+5 and metatarsal 5 relative to metatarsal 1 and hallux
## Biomechanical hypotheses: related to medial transfer of pressure? Degree of plantarflexion at toeoff? Gait type?
## Does amount of substrate deformation have an effect?
## Collect subset of data to test these hypotheses (observations where all variables were measured)
pc4Data<-daasData[complete.cases(daasData[,c("PC4","Gait","Transverse.oblique.P","Footprint.toeoff.ankle","Mean.depth")]),]
## Check variance inflation factors of predictor variables (test for collinearity)
corvif(daasData[,c("Gait","Transverse.oblique.P","Footprint.toeoff.ankle","Mean.depth")])
## All VIFs are less than 3 (suggests no collinearity)
## Create linear mixed effects model of PC4 and various hypothesized predictors - use random intercept model with Subject as random effect
model1PC4<-lme(PC4~Gait+Transverse.oblique.P+Footprint.toeoff.ankle+Mean.depth,random=~1|Subject,data=pc4Data,na.action=na.omit,method="ML")
summary(model1PC4)
## Gait type has least significant effect - remove from model
model2PC4<-update(model1PC4,.~.-Gait)
anova(model1PC4,model2PC4)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model
summary(model2PC4)
## Ankle angle at toe-off has least significant effect - remove from model
model3PC4<-update(model2PC4,.~.-Footprint.toeoff.ankle)
anova(model2PC4,model3PC4)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model
summary(model3PC4)
## Remainder of variables have a significant effect at the p<=0.05 level
## Check structure of residuals
plot(pc4Data$Transverse.oblique.P,model3PC4$residuals[,1])
plot(pc4Data$Mean.depth,model3PC4$residuals[,1])
plot(model3PC4)
plot(density(residuals(model3PC4)))
qqnorm(model3PC4)
## Plot example of biomechanical relationship (PC4 to medial pressure) for one subject (see Figure 7)
bsData<-subset(daasData,Subject=="BS")
bsPC4Plot<-ggplot(bsData,aes(x=Transverse.oblique.P,y=PC4))+
  geom_point(shape=19,size=2,alpha=.9)+
  labs(x="Relative medial pressure",y="PC4")+
  theme(axis.text.x=element_text(size=10,colour="black"),axis.text.y=element_text(size=10,colour="black"),axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),legend.position=c(1,0),legend.justification=c(1,0),legend.background=element_rect(colour="black",
        fill="white"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),axis.ticks=element_line(colour="black"))
print(bsPC4Plot)
jpeg(file="bsPC4Plot.jpg",width=3.25,height=4.25,units="in",quality=100,res=300)
print(bsPC4Plot)
dev.off()

## What is the shape change described by PC5?
pcaRotations[,c(1,6)]
## PC5 described by deep second toe and lateral midfoot relative to fifth toe
## Biomechanical hypotheses: related to medial transfer of pressure for toeoff? pressure on hallux? ankle extension at toeoff? Gait type?
## Does amount of substrate deformation have an effect?
## Collect subset of data to test these hypotheses (observations where all variables were measured)
pc5Data<-daasData[complete.cases(daasData[,c("PC5","Gait","Transverse.oblique.P","Hallux.peak.P","Footprint.toeoff.knee","Footprint.toeoff.ankle","Mean.depth")]),]
## Check variance inflation factors of continuous predictor variables (test for collinearity)
corvif(daasData[,c("Gait","Transverse.oblique.P","Hallux.peak.P","Footprint.toeoff.knee","Footprint.toeoff.ankle","Mean.depth")])
## Gait type appears collinear with knee angle at toeoff (VIF of about 4). If either of these variables is excluded, then all VIFs are less than 3 (suggests no collinearity among other predictors).
## Create linear mixed effects model of PC5 and various hypothesized predictors - use random intercept model with Subject as random effect
model1PC5<-lme(PC5~Gait+Transverse.oblique.P+Hallux.peak.P+Footprint.toeoff.knee+Footprint.toeoff.ankle+Mean.depth,random=~1|Subject,data=pc5Data,na.action=na.omit,method="ML")
summary(model1PC5)
## Knee angle at toeoff has least significant effect - remove from model
model2PC5<-update(model1PC5,.~.-Footprint.toeoff.knee)
anova(model1PC5,model2PC5)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model
summary(model2PC5)
## Relative medial pressure has least significant effect - remove from model
model3PC5<-update(model2PC5,.~.-Transverse.oblique.P)
anova(model2PC5,model3PC5)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model
summary(model3PC5)
## Hallux peak pressure has least significant effect - remove from model
model4PC5<-update(model3PC5,.~.-Hallux.peak.P)
anova(model3PC5,model4PC5)
## Removal does not significantly affect model fit and AIC is reduced - keep out of model
summary(model4PC5)
## All remaining variables have effects that are significant at the p<0.05 level
## Check structure of residuals
plot(pc5Data$Gait,model4PC5$residuals[,1])
plot(pc5Data$Footprint.toeoff.ankle,model4PC5$residuals[,1])
plot(pc5Data$Mean.depth,model4PC5$residuals[,1])
plot(model4PC5)
plot(density(residuals(model4PC5)))
qqnorm(model4PC5)
## Plot example of significant biomechanical relationship (PC5 to ankle angle at toeoff) for one subject (see Figure 7)
bjData<-subset(daasData,Subject=="BJ")
bjPC5Plot<-ggplot(bjData,aes(x=Footprint.toeoff.ankle,y=PC5))+
  geom_point(shape=19,size=2,alpha=.9)+
  labs(x="Ankle angle at toeoff (radians)",y="PC5")+
  theme(axis.text.x=element_text(size=10,colour="black"),axis.text.y=element_text(size=10,colour="black"),axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),legend.position=c(1,0),legend.justification=c(1,0),legend.background=element_rect(colour="black",
        fill="white"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),axis.ticks=element_line(colour="black"))
print(bjPC5Plot)
jpeg(file="bjPC5Plot.jpg",width=3.25,height=4.25,units="in",quality=100,res=300)
print(bjPC5Plot)
dev.off()

## Create grid of plots showing biomechanical effects on PC2-5 axes (Figure 7)
jpeg(file="pcModelPlot.jpg",width=6.5,height=8.5,units="in",quality=100,res=300)
pcModelPlot<-grid.arrange(bkPC2Plot,dhPC3Plot,bsPC4Plot,bjPC5Plot,nrow=2,ncol=2)
dev.off()
