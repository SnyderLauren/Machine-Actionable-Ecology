####################STATISTICAL ANALYSIS################################

# Contrasting effects of landscape composition on crop yield mediated by specialist herbivores
# https://doi.org/10.1002/eap.1695
# Authors: Ricardo Perez-Alvarez, Brian Nault, and Katja Poveda
# Analysis conducted with R version 4.0.4 (2021-02-15)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)

########### LOAD DATA  ##############
#####This script is for Fig. 4c - Effect of meadows (1000m) on lepidotera incidence####
##There are three key sets of data to extract:
#1 - Output from ANOVA (Type III sum of squares) 
#2 - Output of summary function on lme (fixed effects)
#3 - Predicted model values to recreate figure 4c (based on model predictions, not the raw data)

# https://datadryad.org/stash/dataset/doi:10.5061/dryad.484tt
Landscape.affects.pest.and.crop.yield_2023<-read.csv("https://datadryad.org/stash/downloads/file_stream/28280", na.strings=c("NA", "", "N/A"))

#shorten dataset name
LandscapeData <- Landscape.affects.pest.and.crop.yield_2023
#remove the long name dataset, so we accidentally don't end up 
#working with the wrong data
rm(Landscape.affects.pest.and.crop.yield_2023)
#check dataset: should yield 43 obs. of 20 variables
head(LandscapeData)
#Dataset looks fine! Let's check the summary
summary(LandscapeData)
#year is read as a numeric variable. same with Plot_ID
#converting Year and Plot_ID from numeric to character
LandscapeData$Year <- as.character(LandscapeData$Year)
LandscapeData$Plot_ID <- as.character(LandscapeData$Plot_ID)
#check summary again
summary(LandscapeData)
#looks good now!

############ LOAD PACKAGES ##################### 
library(nlme)#nlme_3.1-152  
library(MuMIn)#MuMIn_1.43.17  
library (ggplot2)#ggplot2_3.3.6  
library (visreg)#visreg_2.7.0 
library(effects)#effects_4.2-0
library(orkg)
library (ggpubr)#vesion 0.4.0
library(dplyr)#version 1.08

######## REMOVING MISSING VALUES ########
#Removing the rows with missing data from each variable individually (this way we don't delete more data than we need to)

names(LandscapeData)
FleaBeetlesIncidence <- LandscapeData[!is.na(LandscapeData$FleaBeetles_incidence), ]
FleaBeetlesAbundance <- LandscapeData
LepidopteranIncidence <- LandscapeData[!is.na(LandscapeData$Lepidoptera_incidence), ]
LepidopteranAbundance <-LandscapeData[!is.na(LandscapeData$Lepidoptera_abundance), ]
AphidsIncidence <- LandscapeData[!is.na(LandscapeData$Aphid_incidence), ]
ParasitoidHostRatio <-LandscapeData[!is.na(LandscapeData$ParasitoidHostRatio), ]
PlantDamageIndex <-LandscapeData[!is.na(LandscapeData$Plant_damage), ]


#######CREATING PLOTS BASED ON LME MODEL PREDICTIONS ########
#making model predictions and calculating confidence intervals

####Lepidoptera incidence data - Figure 4c ##############

fitlme.Li <- lme(sqrt(Lepidoptera_incidence)~  mead_1000+Year, data = LepidopteranIncidence, random=~1|Farm_ID/Plot_ID)
LMMOutput <- data.frame(fitlme.Li$coefficients$fixed)
colnames(LMMOutput)[1] <- 'Value'

#First set of data to EXTRACT from anova. Here, we run an ANOVA (Type III sum of squares) and would like to capture: numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), the F-value (test statistic from the F test), and the associated p-value for all three rows - Intercept, Mead_1000 (proportion of meadows within a 1000 meter radius of the field plot) and Year (sampling year)
anovaOutput <- data.frame(anova(fitlme.Li,type='marginal'))
anovaOutput

#Second set of data to EXTRACT from summary output. Here we would like to extract the information associated with the fixed effects: Value (slope estimates), Std.Error (approximate standard error of the slope estimates), DF (denominator degrees of freedom), t- value (ratios between slope estimates and their standard errors), p-value (associated p-value from a t-distribution)
sum1 <- data.frame(summary(fitlme.Li)$tTable, check.names=FALSE)
sum1 

newdat.lme.Li  = data.frame(Year = LepidopteranIncidence$Year,
                            mead_1000 = LepidopteranIncidence$mead_1000,
                            Lepidoptera_incidence=LepidopteranIncidence$Lepidoptera_incidence)
head(newdat.lme.Li)
newdat.lme.Li$predlme.Li = predict(fitlme.Li, newdata = newdat.lme.Li, level = 0)
#ggplot(LepidopteranIncidence, aes(x = mead_1000, y = Lepidoptera_incidence, color = Year) ) +
# geom_rug(sides = "b", size = 1) +
#geom_line(data = newdat.lme.Li, aes(y = predlme.Li), size = 1)
des.Li = model.matrix(formula(fitlme.Li)[-2], newdat.lme.Li)
predvar.Li = diag( des.Li %*% vcov(fitlme.Li) %*% t(des.Li) )

newdat.lme.Li$lower = with(newdat.lme.Li, predlme.Li - 2*sqrt(predvar.Li) )
newdat.lme.Li$upper = with(newdat.lme.Li, predlme.Li + 2*sqrt(predvar.Li) )

p1 <- ggplot(LepidopteranIncidence, aes(x = mead_1000, y = Lepidoptera_incidence, color = Year) ) +
  geom_point()+
  geom_rug(sides = "b", size = 1) +
  geom_ribbon(data = newdat.lme.Li, aes(y = NULL, ymin = lower, ymax = upper,
                                        color = NULL, fill = Year),
              alpha = .15) +
  geom_line(data = newdat.lme.Li, aes(y = predlme.Li), size = .75)+
  theme_classic()+xlab('Proportion of meadows at 1000-m')+ylab('Proportion of plants infested by lepidoptera')

p1
#ggsave("Fig.4c.png", plot = p1, scale=0.5)


############Lepidoptera incidence- FIGURE 4c #############################
#Figure 4c--Keeping the default settings in ggplot
Lepidoptera_incidence.default<- ggplot(LandscapeData, aes(x=mead_1000, y=Lepidoptera_incidence, colour = Year)) +   geom_point(shape=1) + geom_smooth (aes (x=mead_1000, y=Lepidoptera_incidence, colour=factor(Year)), method=lm, se=TRUE, fullrange=TRUE) + theme(panel.background = element_rect(fill='white', colour='black'))+scale_y_continuous(breaks=c(0,0.3, 0.9, 1.2))+xlab('Proportion of meadows at 1000-m')+ylab('Proportion of plants infested by Lepidoptera')
Lepidoptera_incidence.default#we can ignore the warning. It is just telling us that there is missing data
#Figure4c--Reproducing the plot theme selected in the paper (the way the plot looks like in the paper)
Lepidoptera_incidence.papertheme<- ggplot(LandscapeData, aes(x=mead_1000, y=Lepidoptera_incidence, colour = Year, linetype = Year)) +  geom_point(aes(shape = factor(Year)), size = 2) + geom_smooth (aes (x=mead_1000, y=Lepidoptera_incidence, colour=factor(Year)), method=lm, se=TRUE, fullrange=TRUE,show.legend = NA) +
  scale_color_grey (start=0.7, end=0.2)+ scale_linetype_manual(values = c("dashed", "solid"))+ scale_shape_manual(name = "Year",values = c(16, 17))+theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(legend.position = "top")+xlab('Proportion of meadows at 1000-m')+ ylab('Proportion of plants infested by Lepidoptera')
Lepidoptera_incidence.papertheme
#add equation to the plot---
Lepidoptera_incidence.papertheme+ stat_regline_equation()
# Save ggplot figure as png
#ggsave("Fig.4c.2.png", plot = Lepidoptera_incidence.papertheme, scale=0.5)


#not sure how to extract equation from the ggplot figure. Alternatively, we could run the individual models
#first, we need to split the dataset per year
Data2014<-filter (LandscapeData, Year=="2014")
Data2015 <- filter(LandscapeData, Year=="2015")
#now we can create individual models (lm)

Lepidoptera_incidence2014 <-lm (Lepidoptera_incidence~mead_1000, data=Data2014 )
sum2014 <- data.frame(summary(Lepidoptera_incidence2014)$coefficients, check.names=FALSE) #the estimate of mead_250 represents the slope 
rownames(sum2014) <- c("(Intercept)(2014)","mead_1000 (2014)")
Lepidoptera_incidence2015 <-lm (Lepidoptera_incidence~mead_1000, data=Data2015 )
sum2015 <- data.frame(summary(Lepidoptera_incidence2015)$coefficients, check.names=FALSE) #the estimate of mead_250 represents the slope
rownames(sum2015) <- c("(Intercept)(2015)","mead_1000 (2015)")
sumLR <- rbind(sum2014,sum2015)
#keep in mind that the model significance could be different from the significance
#obtained from the linear mixed-effect models (lme)


#Third set of data to EXTRACT. This would allow someone to recreate FIG. 4c based on the model predictions.
PredictedValuesLepidopteraIncidence <- newdat.lme.Li
PredictedValuesLepidopteraIncidence <- subset(PredictedValuesLepidopteraIncidence, select = -c(Lepidoptera_incidence))
PredictedValuesLepidopteraIncidence
# We would like to extract the data from all rows and columns, EXCEPT for the column labeled Lepidoptera_Incidence as this is just the raw data. Description of data in each column: 
#Year: sampling year (either 2014 or 2015)
#mead_1000: same as described above
#predlme.Li: predicted (based on model) percentage of plants in a field that would have one or more lepidoteran larva on it 
#lower: lower limit of the 95% CI
#upper: Upper limit of the 95% CI

###End of script

#Input data set
inputDF <- LepidopteranIncidence[, c("Year", "Farm_ID", "Plot_ID", "Lepidoptera_incidence", "mead_1000")]
