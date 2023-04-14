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


#First set of data to EXTRACT from anova. Here, we run an ANOVA (Type III sum of squares) and would like to capture: numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), the F-value (test statistic from the F test), and the associated p-value for all three rows - Intercept, Plant_damage (estimation of the percent of plant damage caused by insect pests) and Year (sampling year)
anova(fitlme.Cr,type='marginal')


#Second set of data to EXTRACT from summary output. Here we would like to extract the information associated with the fixed effects: Value (slope estimates), Std.Error (approximate standard error of the slope estimates), DF (denominator degrees of freedom), t- value (ratios between slope estimates and their standard errors), p-value (associated p-value from a t-distribution)
sum1 <- summary(fitlme.Cr)
sum1$tTable

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


# Note the y axis transformation. 
p1 <- ggplot(PlantDamageIndex , aes(x = Plant_damage, y = sqrt(FreshHead_biomass), color = Year) ) +
geom_point()+
geom_rug(sides = "b", size = 1) +
 geom_ribbon(data = newdat.lme.Cr, aes(y = NULL, ymin = lower, ymax = upper,
         color = NULL, fill = Year),
  alpha = .15) +
geom_line(data = newdat.lme.Cr, aes(y = predlme.Cr), size = .75)+
theme_classic()+xlab('Plant Damage Index (%) ')+ylab('Cabbage head mass (g)')

p1

# Save ggplot figure as png
ggsave("Fig.5.png", plot = p1, scale=0.5)

#Third set of data to EXTRACT. This would allow someone to recreate FIG. 5 based on the model predictions.
PredictedValuesCropYield <- newdat.lme.Cr
PredictedValuesCropYield <- subset(PredictedValuesCropYield, select = -c(FreshHead_biomass))
PredictedValuesCropYield
# We would like to extract the data from all rows and columns, EXCEPT for the column labeled "FreshHead_biomass" as this is just the raw crop yield data. Description of data in each column: 
#Year: sampling year (either 2014 or 2015)
#Plant_damage: same as described above
#Plot_ID: identification of the plot where the sample was collected
#predlme.Cr: predicted (based on model) biomass of the crop in grams
#lower: lower limit of the 95% CI
#upper: Upper limit of the 95% CI

###End of script


####################################### 
############### ORKG ##################
####################################### 

orkg <- ORKG(host="https://incubating.orkg.org/")

# Template 'Model Fitting 4'
orkg$templates$materialize_template(template_id = "R479769")
tp = orkg$templates$list_templates()

instance <- tp$model_fitting_4(
  label="Plant damage index and cabbage yield", 
  
  
  has_input_dataset="https://doi.org/10.5061/dryad.484tt",
  
  # LandscapeData can be used instead of URI
  #has_input_dataset=tuple(LandscapeData, "Landscape affects pest and crop yield" ),
  
  # Description of the statistical model used
  has_input_model=tp$statistical_model(
    label="Mixed-effect model ",
    is_denoted_by=tp$formula(
      label="The formula for the mixed-effect model",
      
      has_value_specification=tp$value_specification(
        label="formula",
        has_specified_value="formula"
      )
    )
  ),
  
  # Predicted model values to recreate figure 5 (based on model predictions, not the raw data)
  has_output_dataset= tuple(PredictedValuesAphid_incidence, 'Relationship between plant damage index and cabbage yield'),
  
  # Output from ANOVA (Type III sum of squares) 
  has_output_dataset_2= tuple(anova(fitlme.Cr,type='marginal'), 'ANOVA (Type III sum of squares)'),
  
  # Output of summary function on lme (fixed effects)
  has_output_dataset_3= tuple(sum1$tTable, 'Summary of LME (fixed effects)'),
  
  # PNG output from ggplot - Git Repo is currently set to private.
  has_output_figure="https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig.5.png",
  
  # Output statement if applicable.
  has_output_statement= "Relationship between plant damage index (percentage of leaf area removed by herbivores) and cabbage yield (mean mass of
  marketable cabbage heads). Lines are the fixed-effect predictions from the best model without covariables and associated 95% confidence intervals (shaded).",
  
  # Snippet is essentially a concise version of this script with redundant code removed.
  # Git Repo is currently set to private.
  has_implementation="https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig5.snippet.R"
  
)
instance$serialize_to_file("article.contribution.6.json", format="json-ld")
