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

#First set of data to EXTRACT from anova. Here, we run an ANOVA (Type III sum of squares) and would like to capture: numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), the F-value (test statistic from the F test), and the associated p-value for all three rows - Intercept, Mead_250 (proportion of meadows within a 250 meter radius of the field plot) and Year (sampling year)
anova(fitlme.fb,type='marginal')


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

####################################### 
############### ORKG ##################
####################################### 

orkg <- ORKG(host="https://incubating.orkg.org/")

# Template 'Model Fitting 3'
orkg$templates$materialize_template(template_id = "R474043")
tp = orkg$templates$list_templates()

instance <- tp$model_fitting_3(
  label="LMM evaluating the effect of study year and the proportion of meadows on flea beetle abundance in experimental farm plots.", 
  
  
  has_input_dataset="https://doi.org/10.5061/dryad.484tt",
  
  
  # Description of the statistical model used
  has_input_model=tp$statistical_model(
    label="A linear mixed model (LMM) with flea beetle abundance (FleaBeetles_abundance) as the response variable, study year (Year) and the proportion of meadows within a 250 meter
    radius of the experimental farm plot (mead_250) as fixed effects, and farm (Farm_ID) and plot identity (Plot_ID) as random effects.",
    is_denoted_by=tp$formula(
      label="The formula of the LMM with FleaBeetles_abundance as the response variable, Year and mead_250 as fixed effects, and Farm_ID and Plot_ID as random effects.",
      
      has_value_specification=tp$value_specification(
        label="FleaBeetles_abundance ~  mead_250 + Year + (1|Farm_ID/Plot_ID)",
        has_specified_value="FleaBeetles_abundance ~  mead_250 + Year + (1|Farm_ID/Plot_ID)"
      )
    )
  ),
  
  # Output of summary function on lme (fixed effects)
  has_output_dataset= tuple(sum1, 'Results of LMM with FleaBeetles_abundance as the response variable, and Year and mead_250 as fixed effects.'),
  
  
  # PNG output from ggplot - Git Repo is currently set to private.
  has_output_figure="https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig.4b.png",
  
  # Output statement if applicable.
  has_output_statement= "Relationship between the proportion of meadows around the experimental fields (250 m radius) and flea beetle abundance. Lines are the fixed-effect predictions 
  from the best models without covariables and shading represents the associated 95% confidence intervals.",
  
  # Snippet is essentially a concise version of this script with redundant code removed.
  # Git Repo is currently set to private.
  has_implementation="https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig4b.snippet.R"
  
)
instance$serialize_to_file("article.contribution.2.json", format="json-ld")