####################STATISTICAL ANALYSIS################################

# Contrasting effects of landscape composition on crop yield mediated by specialist herbivores
# https://doi.org/10.1002/eap.1695
# Authors: Ricardo Perez-Alvarez, Brian Nault, and Katja Poveda
# Analysis conducted with R version 4.0.4 (2021-02-15)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)

########### LOAD DATA  ##############
#####This script is for Fig. 4d - Effect of meadows (500 m radius) on parasitoid-host ratio####
##There are three key sets of data to extract:
#1 - Output from ANOVA (Type III sum of squares) 
#2 - Output of summary function on lme (fixed effects)
#3 - Predicted model values to recreate figure 4a (based on model predictions, not the raw data)

#read in following CSV file:'Landscape.affects.pest.and.crop.yield_2023.cvs'
Landscape.affects.pest.and.crop.yield_2023<-read.csv("Landscape affects pest and crop yield_2023.csv", na.strings=c("NA", ""))

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


####Parasitoid Host-ratio - figure 4d#########
fitlme.Ph <- lme(ParasitoidHostRatio~  mead_500, data = ParasitoidHostRatio, random=~1|Farm_ID/Plot_ID)
LMMOutput <- data.frame(fitlme.Ph$coefficients$fixed)
colnames(LMMOutput)[1] <- 'Value'


#First set of data to EXTRACT from anova. Here, we run an ANOVA (Type III sum of squares) and would like to capture: numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), the F-value (test statistic from the F test), and the associated p-value for all three rows - Intercept, Mead_250 (proportion of meadows within a 500 meter radius of the field plot) and Year (sampling year)
anovaOutput <- data.frame(anova(fitlme.Ph,type='marginal'))
anovaOutput


#Second set of data to EXTRACT from summary output. Here we would like to extract the information associated with the fixed effects: Value (slope estimates), Std.Error (approximate standard error of the slope estimates), DF (denominator degrees of freedom), t- value (ratios between slope estimates and their standard errors), p-value (associated p-value from a t-distribution)
sum1 <- data.frame(summary(fitlme.Ph)$tTable, check.names=FALSE)
sum1 

newdat.lme.Ph  = data.frame(mead_500 = ParasitoidHostRatio$mead_500,
                            ParasitoidHostRatio=ParasitoidHostRatio$ParasitoidHostRatio
)
head(newdat.lme.Ph )
newdat.lme.Ph$predlme.Ph = predict(fitlme.Ph, newdata = newdat.lme.Ph, level = 0)
#ggplot(ParasitoidHostRatio, aes(x = mead_500, y = ParasitoidHostRatio) ) +
 # geom_rug(sides = "b", size = 1) +
  #geom_line(data = newdat.lme.Ph, aes(y = predlme.Ph), size = 1)
des.Ph = model.matrix(formula(fitlme.Ph)[-2], newdat.lme.Ph)
predvar.Ph = diag( des.Ph %*% vcov(fitlme.Ph) %*% t(des.Ph) )
newdat.lme.Ph$lower = with(newdat.lme.Ph, predlme.Ph - 2*sqrt(predvar.Ph) )
newdat.lme.Ph$upper = with(newdat.lme.Ph, predlme.Ph + 2*sqrt(predvar.Ph) )

p1 <- ggplot(ParasitoidHostRatio, aes(x = mead_500, y = ParasitoidHostRatio) ) +
geom_point()+
geom_rug(sides = "b", size = 1) +
geom_ribbon(data = newdat.lme.Ph, aes(y = NULL, ymin = lower, ymax = upper,
                                    color = NULL),
         alpha = .15) +
geom_line(data = newdat.lme.Ph, aes(y = predlme.Ph), size = .75)+
theme_classic()+xlab('Proportion of meadows at 500-m')+ylab('Parasitoid Host-ratio')
p1
ggsave("Fig.4d.png", plot = p1, scale=0.5)


############Parasitoid-host ratios- FIGURE 4d #############################

#Figure 4d--Keeping the default settings in ggplot
ParasitoidHostRatio.default<- ggplot(LandscapeData, aes(x=mead_500, y=ParasitoidHostRatio)) +   geom_point(shape=1) + geom_smooth (aes (x=mead_500, y=ParasitoidHostRatio), method=lm, se=TRUE, fullrange=TRUE) + theme(panel.background = element_rect(fill='white', colour='black'))+xlab('Proportion of meadows at 500-m')+ylab('ParasitoidHostRatio')
ParasitoidHostRatio.default#we can ignore the warning. It is just telling us that there is missing data
#Figure4c--Reproducing the plot theme selected in the paper (the way the plot looks like in the paper)
ParasitoidHostRatio.papertheme<- ggplot(LandscapeData, aes(x=mead_500, y=ParasitoidHostRatio)) +  geom_point(size = 2) + geom_smooth (aes (x=mead_500, y=ParasitoidHostRatio), method=lm, se=TRUE, fullrange=TRUE,color = "black") +
  scale_color_grey (start=0.7, end=0.2)+theme(panel.background = element_rect(fill='white', colour='black'))+
  xlab('Proportion of meadows at 500-m')+ ylab('ParasitoidHostRatio')
ParasitoidHostRatio.papertheme
#add equation to the plot---

ParasitoidHostRatio.papertheme+ stat_regline_equation()
ggsave("Fig.4d.2.png", plot = ParasitoidHostRatio.papertheme, scale=0.5)

#not sure how to extract equation from the ggplot figure. Alternatively, we could run the individual models
ParasitoidHostRatio.model <-lm (ParasitoidHostRatio~mead_500, data=LandscapeData,na.action=na.omit )
sumLR <- data.frame(summary(ParasitoidHostRatio.model)$coefficients, check.names=FALSE) #the estimate of mead_250 represents the slope

#Third set of data to EXTRACT. This would allow someone to recreate FIG. 4d based on the model predictions.
PredictedValuesParasitoidHostRatio <- newdat.lme.Ph
PredictedValuesParasitoidHostRatio <- subset(PredictedValuesParasitoidHostRatio, select = -c(ParasitoidHostRatio))
# We would like to extract the data from all rows and columns, EXCEPT for the column labeled ParasitoidHostRatio as this is just the raw data. Description of data in each column: 
#Year: sampling year (either 2014 or 2015)
#mead_500: same as described above
#predlme.Ph: predicted (based on model) parasitoid-host ratio of a field
#lower: lower limit of the 95% CI
#upper: Upper limit of the 95% CI



###End of script

#Input data set
inputDF <- ParasitoidHostRatio[, c("Year", "Farm_ID", "Plot_ID", "ParasitoidHostRatio", "mead_500")]

####################################### 
############### ORKG ##################
####################################### 

orkg <- ORKG(host="https://incubating.orkg.org")
#Template 'LMM Planned Process'
orkg$templates$materialize_template(template_id = "R492225")
tp = orkg$templates$list_templates()
keys(tp)
tp$linear_mixed_model_fitting(text= 'doc')
tp$linear_mixed_model(text= 'doc')

##################################
######### Definitions ############
##################################
Meter <- tp$qudt_unit(label="m", qudtucumcode="m", same_as="http://qudt.org/vocab/unit/M")
Radius <- tp$quantity_value(label="500 m radius", qudtnumericvalue= 250, qudtunit=Meter)
Parasitoid <- tp$entity(label="Parasitoid", same_as="http://purl.obolibrary.org/obo/ECOCORE_00000087")
ParasitoidHost <- tp$entity(label="Parasitoid and their host")
Host <- tp$entity(label="Parasitoid", same_as="http://purl.obolibrary.org/obo/ECOCORE_00000087")
Plot <- tp$entity(label="Agricultural experimental plot", same_as="http://purl.obolibrary.org/obo/AGRO_00000301")
Ratio <- tp$property(label="Incidence", same_as="http://purl.obolibrary.org/obo/NCIT_C16726")
Meadow <- tp$entity(label="Meadow", same_as="http://purl.obolibrary.org/obo/ENVO_00000108")
Region <- tp$entity(label="Region", same_as="http://purl.obolibrary.org/obo/NCIT_C41129", is_constrained_by=Radius)
Proportion <- tp$property(label="Proportion", same_as="http://purl.obolibrary.org/obo/NCIT_C49159")
Year <- tp$entity(label="Year", same_as="http://purl.obolibrary.org/obo/NCIT_C29848")
Farm <- tp$entity(label="Farm", same_as="http://purl.obolibrary.org/obo/NCIT_C48953")


################################
######## LMM Variables #########
################################
var_parasitoid_host_ratio <- tp$variable(
  label="Parasitoid-host ratio",
  has_object_of_interest_= ParasitoidHost,
  has_matrix = Plot,
  has_property = Ratio
)

var_mead_500 <- tp$variable(
  label="Meadow proportion in a region with 500 m radius",
  has_object_of_interest_= Meadow,
  has_matrix = Region,
  has_property = Proportion
)

var_Plot_Farm <- tp$variable(
  label="Plot and farm",
  has_object_of_interest_ = Plot,
  has_matrix = Farm,
  
)


################################
############# LR  #############
################################
lr <- tp$linear_regression(
  label="A linear regression (LR) with parasitoid-host ratio (ParasitoidHostRatio) as the dependent variable and the proportion of meadows within a 500 meter radius of the experimental 
  farm plot (mead_500) as the independent variable",
  has_input_dataset=tuple(inputDF, "Raw field data with parasitoid-host ratio"),
  has_dependent_variable = var_parasitoid_host_ratio,
  has_independent_variable = var_mead_500,
  has_output_figure="https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig.4d.2.png",
  has_output_statement= "Relationship between the proportion of meadows around the experimental fields (500 m radius) and parasitoid-host ratio. Lines are the fixed-effect predictions 
  from the best models without covariables and shading represents the associated 95% confidence intervals",
  has_output_dataset = tuple(sumLR, "Results of LR with ParasitoidHostRatio as the dependent variable and mead_500 as the independent variable")
)


################################
############# LMM  #############
################################
lmm <- tp$linear_mixed_model(
  label="A linear mixed model (LMM) with parasitoid-host ratio (ParasitoidHostRatio) as the response variable,  the proportion of meadows within a 500 meter radius of the experimental 
  farm plot (mead_500) as fixed effects, and farm (Farm_ID) and plot identity (Plot_ID) as random effects",
  has_response_variable = var_parasitoid_host_ratio,
  has_fixed_effect_term_i = var_mead_500,
  has_random_effect_term_ = var_Plot_Farm,
)


################################
######### LMM Fitting  #########
################################
lmmFitting <- tp$linear_mixed_model_fitting(
  label="A linear mixed model (LMM) fitting with parasitoid-host ratio (ParasitoidHostRatio) as the response variable,the proportion of meadows within a 500 meter radius of the experimental 
  farm plot (mead_500) as fixed effects, and farm (Farm_ID) and plot identity (Plot_ID) as random effects",
  has_input_dataset= tuple(inputDF, "Raw field data on parasitoid-host ratio"),
  has_input_model=lmm,
  has_output_dataset= tuple(LMMOutput, 'Results of LMM fitting with ParasitoidHostRatio as the response variable and mead_500 as a fixed effect'),
)

################################
## LMM Significance Testing  ###
################################
LMMSignificanceTesting <- tp$lmm_significance_testing(
  label="Significance testing to attain p-values for Intercept and mead_500 (proportion of meadows within a 500 meter radius of the field plot)",
  has_input_dataset= tuple(LMMOutput, "Fitted LMM with ParasitoidHostRatio as the response variable and mead_500 as a fixed effect"),
  has_output_dataset= tuple(sum1, 'Value (slope estimates), Std.Error, DF, t- value and p-value for fixed effect'),
)


################################
############ Anova  ############
################################
ANOVA <- tp$anova(
  label="ANOVA to attain F-values for Intercept and mead_500 (proportion of meadows within a 500 meter radius of the field plot)",
  has_input_dataset= tuple(LMMOutput, "Fitted LMM with ParasitoidHostRatio as the response variable and mead_500 as a fixed effect"),
  has_output_dataset= tuple(anovaOutput, 'numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), F-value, and the associated p-value for fixed effects'),
)


################################
######## LMM Prediction  #######
################################
LMMPrediction <- tp$lmm_prediction(
  label="Prediction using the fitted LMM with ParasitoidHostRatio as the response variable and mead_500 as a fixed effect.",
  has_input_dataset= tuple(LMMOutput, "Fitted LMM with ParasitoidHostRatio as the response variable and mead_500 as a fixed effect."),
  has_output_dataset= tuple(PredictedValuesParasitoidHostRatio, 'Predicted results of the LMM with ParasitoidHostRatio as the response variable and mead_500 as a fixed effect'),
  has_output_figure = "https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig.4d.png",
)


################################
#### LMM Planned Process #######
################################
instance <- tp$lmm_planned_process(
  has_implementation= "https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig4d.snippet.R",
  label="Parasitoid-host ratio in experimental farm plots evaluated by a LMM with the proportion of meadows within a 500 meter radius as a fixed effect", 
  has_lmm_fitting= lmmFitting,
  has_anova = ANOVA,
  has_lmm_significance_testing = LMMSignificanceTesting,
  has_lmm_prediction = LMMPrediction,
  has_linear_regression= lr
)
instance$serialize_to_file("article.contribution.4.json", format="json-ld")
#instance$pretty_print()