####################STATISTICAL ANALYSIS################################

# Contrasting effects of landscape composition on crop yield mediated by specialist herbivores
# https://doi.org/10.1002/eap.1695
# Authors: Ricardo Perez-Alvarez, Brian Nault, and Katja Poveda
# Analysis conducted with R version 4.0.4 (2021-02-15)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)

########### LOAD DATA  ##############
#####This script is for Fig. 4b - Effect of meadows (250 m radius) on flea beetle abundance####
##There are three key sets of data to extract:
#1 - Output from ANOVA (Type III sum of squares) 
#2 - Output of summary function on lme (fixed effects)
#3 - Predicted model values to recreate figure 4b (based on model predictions, not the raw data)

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

#######CREATING PLOTS BASED ON MODEL PREDICTIONS ########
#making model predictions and calculating confidence intervals

###### Flea beetle abundance Fig 4b###############
fitlme.fb <- lme(FleaBeetles_abundance~  mead_250+Year, data = LandscapeData, random=~1|Farm_ID/Plot_ID)
LMMOutput <- data.frame(fitlme.fb$coefficients$fixed)
colnames(LMMOutput)[1] <- 'Value'

#First set of data to EXTRACT from anova. Here, we run an ANOVA (Type III sum of squares) and would like to capture: numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), the F-value (test statistic from the F test), and the associated p-value for all three rows - Intercept, Mead_250 (proportion of meadows within a 250 meter radius of the field plot) and Year (sampling year)
anovaOutput <- data.frame(anova(fitlme.fb,type='marginal'))
anovaOutput

#Second set of data to EXTRACT from summary output. Here we would like to extract the information associated with the fixed effects: Value (slope estimates), Std.Error (approximate standard error of the slope estimates), DF (denominator degrees of freedom), t- value (ratios between slope estimates and their standard errors), p-value (associated p-value from a t-distribution)
sum1 <- data.frame(summary(fitlme.fb)$tTable, check.names=FALSE)
sum1

newdat.lme.fb  = data.frame(Year = LandscapeData$Year,
                            mead_250 = LandscapeData$mead_250,
                            FleaBeetles_abundance=LandscapeData$FleaBeetles_abundance
)
head(newdat.lme.fb )
newdat.lme.fb$predlme.fb = predict(fitlme.fb, newdata = newdat.lme.fb, level = 0)
#ggplot(LandscapeData, aes(x = mead_250, y = FleaBeetles_abundance, color = Year) ) +
 # geom_rug(sides = "b", size = 1) +
 # geom_line(data = newdat.lme.fb, aes(y = predlme.fb), size = 1)

des.fb = model.matrix(formula(fitlme.fb)[-2], newdat.lme.fb)
predvar.fb = diag( des.fb %*% vcov(fitlme.fb) %*% t(des.fb) )
newdat.lme.fb$lower = with(newdat.lme.fb, predlme.fb - 2*sqrt(predvar.fb) )
newdat.lme.fb$upper = with(newdat.lme.fb, predlme.fb + 2*sqrt(predvar.fb) )

p1 <- ggplot(LandscapeData, aes(x = mead_250, y = FleaBeetles_abundance, color = Year) ) +
 geom_point()+
 geom_rug(sides = "b", size = 1) +
 geom_ribbon(data = newdat.lme.fb, aes(y = NULL, ymin = lower, ymax = upper,
                                       color = NULL, fill = Year),
             alpha = .15) +
 geom_line(data = newdat.lme.fb, aes(y = predlme.fb), size = .75)+
 theme_classic()+xlab('Proportion of meadows at 250-m')+ylab('Mean flea beetles abundance/trap')

p1
# Save ggplot figure as png
ggsave("Fig.4b.png", plot = p1, scale=0.5)


###########FIGURE 4b---Flea beetles abundance#############################
#Figure 4b--Keeping the default settings in ggplot
FleaBeetles_abundance.default<- ggplot(LandscapeData, aes(x=mead_250, y=FleaBeetles_abundance, colour = Year)) +   geom_point(shape=1) + geom_smooth (aes (x=mead_250, y=FleaBeetles_abundance, colour=factor(Year)), method=lm, se=TRUE, fullrange=TRUE) + theme(panel.background = element_rect(fill='white', colour='black'))+scale_y_continuous(breaks=c(0,5, 10, 15))+xlab('Proportion of meadows at 250-m')+ylab('Mean flea beetle abundance/trap')
FleaBeetles_abundance.default#we can ignore the warning. It is just telling us that there is missing data
#Figure4b--Reproducing the plot theme selected in the paper (the way the plot looks like in the paper)
FleaBeetles_abundance.papertheme<- ggplot(LandscapeData, aes(x=mead_250, y=FleaBeetles_abundance, colour = Year, linetype = Year)) +  geom_point(aes(shape = factor(Year)), size = 2) + geom_smooth (aes (x=mead_250, y=FleaBeetles_abundance, colour=factor(Year)), method=lm, se=TRUE, fullrange=TRUE,show.legend = NA) +
  scale_color_grey (start=0.7, end=0.2)+ scale_linetype_manual(values = c("dashed", "solid"))+ scale_shape_manual(name = "Year",values = c(16, 17))+theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(legend.position = "top")+xlab('Proportion of meadows at 250-m')+ ylab('Mean flea beetle abundance/trap')
FleaBeetles_abundance.papertheme
#add equation to the plot---
FleaBeetles_abundance.papertheme+ stat_regline_equation()
# Save ggplot figure as png
ggsave("Fig.4b.2.png", plot = FleaBeetles_abundance.papertheme, scale=0.5)

#not sure how to extract equation from the ggplot figure. Alternatively, we could run the individual models
#first, we need to split the dataset per year
Data2014<-filter (LandscapeData, Year=="2014")
Data2015 <- filter(LandscapeData, Year=="2015")
#now we can create individual models (lm)

FleaBeetles_abundance2014 <-lm (FleaBeetles_abundance~mead_1000, data=Data2014 )
sum2014 <- data.frame(summary(FleaBeetles_abundance2014)$coefficients, check.names=FALSE) #the estimate of mead_250 represents the slope 
rownames(sum2014) <- c("(Intercept)(2014)","mead_1000 (2014)")
FleaBeetles_abundance2015 <-lm (FleaBeetles_abundance~mead_1000, data=Data2015 )
sum2015 <- data.frame(summary(FleaBeetles_abundance2015)$coefficients, check.names=FALSE) #the estimate of mead_250 represents the slope
rownames(sum2015) <- c("(Intercept)(2015)","mead_1000 (2015)")
sumLR <- rbind(sum2014,sum2015)
#keep in mind that the model significance could be different from the significance
#obtained from the linear mixed-effect models (lme)

#Third set of data to EXTRACT. This would allow someone to recreate FIG. 4b based on the model predictions.
PredictedValuesFleaBeetlesAbundance <- newdat.lme.fb
PredictedValuesFleaBeetlesAbundance <- subset(PredictedValuesFleaBeetlesAbundance, select = -c(FleaBeetles_abundance))
PredictedValuesFleaBeetlesAbundance

# We would like to extract the data from all rows and columns, EXCEPT for the column labeled FleaBeetles_Abundance as this is just the raw data. Description of data in each column: 
#Year: sampling year (either 2014 or 2015)
#mead_250: same as described above
#predlme.fb: predicted (based on model) number of flea beetle adults/pitfall trip  
#lower: lower limit of the 95% CI
#upper: Upper limit of the 95% CI

###End of script

#Input data set
inputDF <- FleaBeetlesAbundance[, c("Year", "Farm_ID", "Plot_ID", "FleaBeetles_abundance", "mead_250")]

####################################### 
############### ORKG ##################
####################################### 

orkg <- ORKG(host="https://incubating.orkg.org")
# Template 'Model Fitting 3'
orkg$templates$materialize_template(template_id = "R492225")
tp = orkg$templates$list_templates()
keys(tp)
tp$linear_mixed_model_fitting(text= 'doc')
tp$linear_mixed_model(text= 'doc')

##################################
######### Definitions ############
##################################
Meter <- tp$qudt_unit(label="m", qudtucumcode="m", same_as="http://qudt.org/vocab/unit/M")
Radius <- tp$quantity_value(label="250 m radius", qudtnumericvalue= 250, qudtunit=Meter)
Fleabeetle <- tp$entity(label="Flea beetle", same_as="http://purl.obolibrary.org/obo/OMIT_0009010")
Plot <- tp$entity(label="Agricultural experimental plot", same_as="http://purl.obolibrary.org/obo/AGRO_00000301")
Abundance <- tp$property(label="Abundance", same_as="http://purl.obolibrary.org/obo/NCIT_C70589")
Meadow <- tp$entity(label="Meadow", same_as="http://purl.obolibrary.org/obo/ENVO_00000108", is_constrained_by=Radius)
Region <- tp$entity(label="Region", same_as="http://purl.obolibrary.org/obo/NCIT_C41129")
Proportion <- tp$property(label="Proportion", same_as="http://purl.obolibrary.org/obo/NCIT_C49159")
Year <- tp$entity(label="Year", same_as="http://purl.obolibrary.org/obo/NCIT_C29848")
Farm <- tp$entity(label="Farm", same_as="http://purl.obolibrary.org/obo/NCIT_C48953")



################################
######## LMM Variables #########
################################
var_FleaBeetles_abundance <- tp$variable(
  label="Flea beetle abundance in an agricultural experimental plot",
  has_object_of_interest_= Fleabeetle,
  has_matrix = Plot,
  has_property = Abundance
)

var_mead_250 <- tp$variable(
  label="Meadow proportion in a region with 250 m radius",
  has_object_of_interest_= Meadow,
  has_matrix = Region,
  has_property = Proportion
)

var_Year <- tp$variable(
  label="Year",
  has_object_of_interest= Year,
  
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
  label="A linear regression (LR) with  flea beetle abundance (FleaBeetles_abundance) as the dependent variable and the proportion of meadows within a 250 meter radius of the experimental 
  farm plot (mead_250) as the independent variable, grouped by year (Year)",
  has_input_dataset=tuple(inputDF, "Raw field data on flea beetle abundance"),
  has_dependent_variable = var_FleaBeetles_abundance,
  has_independent_variable = var_mead_250,
  has_output_figure="https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig.4b.2.png",
  has_output_statement= "Relationship between the proportion of meadows around the experimental fields (250 m radius) and  flea beetle abundance. Lines are the fixed-effect predictions 
  from the best models without covariables and shading represents the associated 95% confidence intervals.",
  has_output_dataset = tuple(sumLR, "Results of LR with FleaBeetles_abundance as the dependent variable and mead_250 as the independent variable")
)


################################
############# LMM  #############
################################
lmm <- tp$linear_mixed_model(
  label="A linear mixed model (LMM) with flea beetle abundance (FleaBeetles_abundance) as the response variable, study year (Year) and
  the proportion of meadows within a 250 meter radius of the experimental farm plot (mead_250) as fixed effects, and farm (Farm_ID) and plot identity (Plot_ID) as random effects.",
  has_response_variable = var_FleaBeetles_abundance,
  has_fixed_effect_term_i = var_mead_250,
  has_fixed_effect_term_ii = var_Year,
  has_random_effect_term_ = var_Plot_Farm,
)

################################
######### LMM Fitting  #########
################################
lmmFitting <- tp$linear_mixed_model_fitting(
  label="A linear mixed model (LMM) with flea beetle abundance (FleaBeetles_abundance) as the response variable, study year (Year) and the proportion of meadows within a 250 meter radius of the experimental 
  farm plot (mead_250) as fixed effects, and farm (Farm_ID) and plot identity (Plot_ID) as random effects.",
  has_input_dataset= tuple(inputDF, "Raw field data on flea beetle abundance"),
  has_input_model=lmm,
  has_output_dataset= tuple(LMMOutput, 'Results of LMM with FleaBeetles_abundance as the response variable, and Year and mead_250 as fixed effects.'),
)

################################
## LMM Significance Testing  ###
################################
LMMSignificanceTesting <- tp$lmm_significance_testing(
  label="Significance testing to attain p-values for Intercept, mead_250 (proportion of meadows within a 250 meter radius of the field plot) and Year (sampling year)",
  has_input_dataset= tuple(LMMOutput, "Fitted LMM with FleaBeetles_abundance as the response variable, and Year and mead_250 as fixed effects."),
  has_output_dataset= tuple(sum1, 'Value (slope estimates), Std.Error, DF, t- value and p-value for fixed effects'),
)

################################
############ Anova  ############
################################
ANOVA <- tp$anova(
  label="ANOVA to attain F-values for Intercept, mead_250 (proportion of meadows within a 250 meter radius of the field plot) and Year (sampling year)",
  has_input_dataset= tuple(LMMOutput, "Fitted LMM with FleaBeetles_abundance as the response variable, and Year and mead_250 as fixed effects."),
  has_output_dataset= tuple(anovaOutput, 'numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), F-value, and the associated p-value for fixed effects.'),
)

################################
######## LMM Prediction  #######
################################
LMMPrediction <- tp$lmm_prediction(
  label="Prediction using the fitted LMM with FleaBeetles_abundance as the response variable, and Year and mead_250 as fixed effects.",
  has_input_dataset= tuple(LMMOutput, "Fitted LMM with FleaBeetles_abundance as the response variable, and Year and mead_250 as fixed effects."),
  has_output_dataset= tuple(PredictedValuesFleaBeetlesAbundance, 'Predicted results of the LMM with FleaBeetles_abundance as the response variable, and Year and mead_250 as fixed effects.'),
  has_output_figure = "https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig.4b.png",
)

################################
#### LMM Planned Process #######
################################
instance <- tp$lmm_planned_process(
  has_implementation= "https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig4b.snippet.R",
  label="Flea beetle abundance in experimental farm plots evaluated by a LMM with study year and the proportion of meadows within a 250 meter radius as fixed effects.", 
  has_lmm_fitting= lmmFitting,
  has_anova = ANOVA,
  has_lmm_significance_testing = LMMSignificanceTesting,
  has_lmm_prediction = LMMPrediction,
  has_linear_regression= lr
)
instance$serialize_to_file("article.contribution.2.json", format="json-ld")
#instance$pretty_print()