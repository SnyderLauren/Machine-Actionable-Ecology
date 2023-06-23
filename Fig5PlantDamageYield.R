####################STATISTICAL ANALYSIS################################

# Contrasting effects of landscape composition on crop yield mediated by specialist herbivores
# https://doi.org/10.1002/eap.1695
# Authors: Ricardo Perez-Alvarez, Brian Nault, and Katja Poveda
# Analysis conducted with R version 4.0.4 (2021-02-15)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)

########### LOAD DATA  ##############
#####This script is for Fig. 5 - Relationship between plant damage and crop yield####
##There are three key sets of data to extract:
#1 - Output from ANOVA (Type III sum of squares) 
#2 - Output of summary function on lme (fixed effects)
#3 - Predicted model values to recreate figure 5 (based on model predictions, not the raw data)

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

#####Relationship between plant damage and crop yield - Figure 5##################
fitlme.Cr <- lme(sqrt(FreshHead_biomass)~  Plant_damage+Year, data = PlantDamageIndex , random=~1|Farm_ID/Plot_ID)
LMMOutput <- data.frame(fitlme.Cr$coefficients$fixed)
colnames(LMMOutput)[1] <- 'Value'

#First set of data to EXTRACT from anova. Here, we run an ANOVA (Type III sum of squares) and would like to capture: numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), the F-value (test statistic from the F test), and the associated p-value for all three rows - Intercept, Plant_damage (estimation of the percent of plant damage caused by insect pests) and Year (sampling year)
anovaOutput <- data.frame(anova(fitlme.Cr,type='marginal'))
anovaOutput


#Second set of data to EXTRACT from summary output. Here we would like to extract the information associated with the fixed effects: Value (slope estimates), Std.Error (approximate standard error of the slope estimates), DF (denominator degrees of freedom), t- value (ratios between slope estimates and their standard errors), p-value (associated p-value from a t-distribution)
sum1 <- data.frame(summary(fitlme.Cr)$tTable, check.names=FALSE)
sum1 

newdat.lme.Cr  = data.frame(Year = PlantDamageIndex$Year,
                            Plant_damage = PlantDamageIndex$Plant_damage,
                            FreshHead_biomass=PlantDamageIndex$FreshHead_biomass,
                            Plot_ID= PlantDamageIndex$Plot_ID
)
head(newdat.lme.Cr )
newdat.lme.Cr$predlme.Cr = predict(fitlme.Cr, newdata = newdat.lme.Cr, level = 0)
#ggplot(PlantDamageIndex , aes(x = Plant_damage, y = sqrt(FreshHead_biomass), color = Year) ) +
 # geom_rug(sides = "b", size = 1) +
 # geom_line(data = newdat.lme.Cr, aes(y = predlme.Cr), size = 1)
des.Cr = model.matrix(formula(fitlme.Cr)[-2], newdat.lme.Cr)
predvar.Cr = diag( des.Cr %*% vcov(fitlme.Cr) %*% t(des.Cr) )
newdat.lme.Cr$lower = with(newdat.lme.Cr, predlme.Cr - 2*sqrt(predvar.Cr) )
newdat.lme.Cr$upper = with(newdat.lme.Cr, predlme.Cr + 2*sqrt(predvar.Cr) )

p1 <- ggplot(PlantDamageIndex , aes(x = Plant_damage, y = sqrt(FreshHead_biomass), color = Year) ) +
geom_point()+
geom_rug(sides = "b", size = 1) +
 geom_ribbon(data = newdat.lme.Cr, aes(y = NULL, ymin = lower, ymax = upper,
         color = NULL, fill = Year),
  alpha = .15) +
geom_line(data = newdat.lme.Cr, aes(y = predlme.Cr), size = .75)+
theme_classic()+xlab('Plant Damage Index (%) ')+ylab('Cabbage head mass (g)')
p1
ggsave("Fig.5.png", plot = p1, scale=0.5)

#### Fig. 5 ####

PlantDamageIndex.papertheme<- ggplot(LandscapeData, aes(x=Plant_damage, y=FreshHead_biomass, colour = Year, linetype = Year)) +  geom_point(aes(shape = factor(Year)), size = 2) + geom_smooth (aes (x=Plant_damage, y=FreshHead_biomass, colour=factor(Year)), method=lm, se=TRUE, fullrange=TRUE,show.legend = NA) + 
scale_color_grey (start=0.7, end=0.2)+ scale_linetype_manual(values = c("dashed", "solid"))+ scale_shape_manual(name = "Year",values = c(16, 17))+theme(panel.background = element_rect(fill='white', colour='black'))+
theme(legend.position = "top")+xlab('Plant Damage Index (%)')+ ylab('Cabbage head mass (g)')+coord_cartesian(ylim = c(0,2500))
PlantDamageIndex.papertheme

ggsave("Fig.5.2.png", plot = PlantDamageIndex.papertheme, scale=0.5)

#not sure how to extract equation from the ggplot figure. Alternatively, we could run the individual models
#First, we need to split the dataset per year
Data2014<-filter (LandscapeData, Year=="2014")
Data2015 <- filter(LandscapeData, Year=="2015")
#now we can create individual models (lm)

PlantDamageIndex2014 <-lm (FreshHead_biomass~Plant_damage, data=Data2014 )
sum2014 <- data.frame(summary(PlantDamageIndex2014)$coefficients, check.names=FALSE) #the estimate of mead_250 represents the slope 
rownames(sum2014) <- c("(Intercept)(2014)","Plant_damage (2014)")
PlantDamageIndex2015 <-lm (FreshHead_biomass~Plant_damage, data=Data2015 )
sum2015 <- data.frame(summary(PlantDamageIndex2015)$coefficients, check.names=FALSE) #the estimate of mead_250 represents the slope
rownames(sum2015) <- c("(Intercept)(2015)","Plant_damage (2015)")
sumLR <- rbind(sum2014,sum2015)

#Third set of data to EXTRACT. This would allow someone to recreate FIG. 5 based on the model predictions.
PredictedValuesCropYield <- newdat.lme.Cr
PredictedValuesCropYield
PredictedValuesCropYield <- subset(PredictedValuesCropYield, select = -c(FreshHead_biomass))
# We would like to extract the data from all rows and columns, EXCEPT for the column labeled "FreshHead_biomass" as this is just the raw crop yield data. Description of data in each column: 
#Year: sampling year (either 2014 or 2015)
#Plant_damage: same as described above
#Plot_ID: identification of the plot where the sample was collected
#predlme.Cr: predicted (based on model) biomass of the crop in grams
#lower: lower limit of the 95% CI
#upper: Upper limit of the 95% CI

###End of script


#Input data set
inputDF <- PlantDamageIndex[, c("Year", "Farm_ID", "Plot_ID", "FreshHead_biomass", "Plant_damage")]

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
CabbageHead <- tp$entity(label="Cabbage Head", same_as="http://purl.obolibrary.org/obo/FOODON_00003406")
CabbageLeaf <- tp$entity(label="Cabbage Leaf", same_as="http://purl.obolibrary.org/obo/FOODON_00003514")
Damage <- tp$property(label="Damage (percentage of leaf area removed by herbivores)")
MeanMass <- tp$property(label="Mean Mass (g)")
Plot <- tp$entity(label="Agricultural experimental plot", same_as="http://purl.obolibrary.org/obo/AGRO_00000301")
Mean <- tp$property(label="Incidence", same_as="http://purl.obolibrary.org/obo/NCIT_C16726")
Meadow <- tp$entity(label="Meadow", same_as="http://purl.obolibrary.org/obo/ENVO_00000108")
Region <- tp$entity(label="Region", same_as="http://purl.obolibrary.org/obo/NCIT_C41129", is_constrained_by=Radius)
Proportion <- tp$property(label="Proportion", same_as="http://purl.obolibrary.org/obo/NCIT_C49159")
Year <- tp$entity(label="Year", same_as="http://purl.obolibrary.org/obo/NCIT_C29848")
Farm <- tp$entity(label="Farm", same_as="http://purl.obolibrary.org/obo/NCIT_C48953")


################################
######## LMM Variables #########
################################
var_fresh_head_biomass <- tp$variable(
  label="Mean mass (g) of marketable cabbage heads",
  has_object_of_interest_= CabbageHead,
  has_matrix = Plot,
  has_property = MeanMass
)

var_plant_damage <- tp$variable(
  label="Plant damage (percentage of leaf area removed by herbivores)",
  has_object_of_interest_= CabbageLeaf,
  has_property = Damage
)
var_Plot_Farm <- tp$variable(
  label="Plot and farm",
  has_object_of_interest_ = Plot,
  has_matrix = Farm,
  
)

var_Year <- tp$variable(
  label="Year",
  has_object_of_interest= Year,
  
)


################################
############# LR  #############
################################
lr <- tp$linear_regression(
  label="A linear regression (LR) with mean mass (g) of marketable cabbage heads (FreshHead_biomass) as the dependent variable and the percentage of leaf area removed by herbivores (Plant_damage) as the independent variable",
  has_input_dataset=tuple(inputDF, "Raw field data on cabbage head mean mass and plant damage"),
  has_dependent_variable = var_fresh_head_biomass,
  has_independent_variable = var_plant_damage,
  has_output_figure="https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig.5.2.png",
  has_output_statement= "Relationship between plant damage index (percentage of leaf area removed by herbivores) and cabbage yield (mean mass of marketable cabbage heads). Lines are the fixed-effect predictions from the best model without covariables and associated 95% confidence intervals (gray shaded)",
  has_output_dataset = tuple(sumLR, "Results of LR with FreshHead_biomass as the dependent variable and Plant_damage as the independent variable")
)
################################
############# LMM  #############
################################
lmm <- tp$linear_mixed_model(
  label="A linear mixed model (LMM) with mean mass (g) of marketable cabbage heads (FreshHead_biomass) as the response variable, study year (Year) and the percentage of leaf area removed by herbivores (Plant_damage) as fixed effects, and farm (Farm_ID) and plot identity (Plot_ID) as random effects",
  has_response_variable = var_fresh_head_biomass,
  has_fixed_effect_term_i = var_plant_damage,
  has_fixed_effect_term_ii = var_Year,
  has_random_effect_term_ = var_Plot_Farm,
)

################################
######### LMM Fitting  #########
################################
lmmFitting <- tp$linear_mixed_model_fitting(
  label="A linear mixed model (LMM) fitting with mean mass (g) of marketable cabbage heads (FreshHead_biomass) as the response variable, study year (Year) and the percentage of leaf area removed by herbivores (Plant_damage) as fixed effects, and farm (Farm_ID) and plot identity (Plot_ID) as random effects",
  has_input_dataset= tuple(inputDF, "Raw field data on cabbage heads"),
  has_input_model=lmm,
  has_output_dataset= tuple(LMMOutput, 'Results of LMM fitting with FreshHead_biomass as the response variable, and Year and Plant_damage as fixed effects'),
)

################################
## LMM Significance Testing  ###
################################
LMMSignificanceTesting <- tp$lmm_significance_testing(
  label="Significance testing to attain p-values for Intercept, Plant_damage (the percentage of leaf area removed by herbivores) and Year (sampling year)",
  has_input_dataset= tuple(LMMOutput, "Fitted LMM with FreshHead_biomass as the response variable, and Year and Plant_damage as fixed effects"),
  has_output_dataset= tuple(sum1, 'Value (slope estimates), Std.Error, DF, t- value and p-value for fixed effects'),
)

################################
############ Anova  ############
################################
ANOVA <- tp$anova(
  label="ANOVA to attain F-values for Intercept, Plant_damage (the percentage of leaf area removed by herbivores) and Year (sampling year)",
  has_input_dataset= tuple(LMMOutput, "Fitted LMM with FreshHead_biomass as the response variable, and Plant_damage and Year as fixed effects"),
  has_output_dataset= tuple(anovaOutput, 'numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), F-value, and the associated p-value for fixed effects'),
)


################################
######## LMM Prediction  #######
################################
LMMPrediction <- tp$lmm_prediction(
  label="Prediction using the fitted LMM with FreshHead_biomass as the response variable, and Year and Plant_damage as fixed effects",
  has_input_dataset= tuple(LMMOutput, "Fitted LMM with FreshHead_biomass as the response variable, and Year and Plant_damage as fixed effects"),
  has_output_dataset= tuple(PredictedValuesCropYield, 'Predicted results of the LMM with FreshHead_biomass as the response variable, and Year and Plant_damage as fixed effects'),
  has_output_figure = "https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig.5.png",
)

################################
#### LMM Planned Process #######
################################
instance <- tp$lmm_planned_process(
  has_implementation= "https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig5.snippet.R",
  label="Mean mass of marketable cabbage heads in experimental farm plots evaluated by a LMM with study year and the percentage of leaf area removed by herbivores as fixed effects", 
  has_lmm_fitting= lmmFitting,
  has_anova = ANOVA,
  has_lmm_significance_testing = LMMSignificanceTesting,
  has_lmm_prediction = LMMPrediction,
  has_linear_regression= lr
)
instance$serialize_to_file("article.contribution.5.json", format="json-ld")
#instance$pretty_print()