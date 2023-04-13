####################STATISTICAL ANALYSIS################################

# Contrasting effects of landscape composition on crop yield mediated by specialist herbivores
# https://doi.org/10.1002/eap.1695
# Authors: Ricardo Perez-Alvarez, Brian Nault, and Katja Poveda
# Analysis conducted with R version 4.0.4 (2021-02-15)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)

########### LOAD DATA  ##############
#####This script is for Fig. 4a - Effect of meadows (250 m radius) on aphid incidence####
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

##########Aphid Incidence -Figure 4a##############

#The first step is to fit a simplified mixed-effect model (nlme)
fitlme <- lme(sqrt(Aphid_incidence)~  mead_250+Year, data = AphidsIncidence, random=~1|Farm_ID/Plot_ID)


#First set of data to EXTRACT from anova. Here, we run an ANOVA (Type III sum of squares) and would like to capture: numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), the F-value (test statistic from the F test), and the associated p-value for all three rows - Intercept, Mead_250 (proportion of meadows within a 250 meter radius of the field plot) and Year (sampling year)
anova (fitlme,type='marginal')


#Second set of data to EXTRACT from summary output. Here we would like to extract the information associated with the fixed effects: Value (slope estimates), Std.Error (approximate standard error of the slope estimates), DF (denominator degrees of freedom), t- value (ratios between slope estimates and their standard errors), p-value (associated p-value from a t-distribution)
sum1 <- summary(fitlme)
sum1$tTable

#Next, create a new dataset from which we will make predictions
newdat.lme = data.frame(Year = AphidsIncidence$Year,
                        mead_250 = AphidsIncidence$mead_250,
                        Aphid_incidence=AphidsIncidence$Aphid_incidence
)
head(newdat.lme)

#I can add the predicted values to the dataset.
newdat.lme$predlme = predict(fitlme, newdata = newdat.lme, level = 0)
#And then use these in geom_line() to add fitted lines based on the new predlm variable. (figure code just provided for reference)
#ggplot(AphidsIncidence, aes(x = mead_250, y = Aphid_incidence, color = Year) ) +
 # geom_rug(sides = "b", size = 1) +
  #geom_line(data = newdat.lme, aes(y = predlme), size = 1)

#create confidence intervals for lme objects
des = model.matrix(formula(fitlme)[-2], newdat.lme)
#Then we use matrix multiplication on the model matrix and variance-covariance matrix
#extracted from the model with vcov(). We pull out the values on the diagonal, which are the
#variances of the predicted values.
predvar = diag( des %*% vcov(fitlme) %*% t(des) )
#I add the confidence interval limits to the dataset for plotting.
newdat.lme$lower = with(newdat.lme, predlme - 2*sqrt(predvar) )
newdat.lme$upper = with(newdat.lme, predlme + 2*sqrt(predvar) )

#Here's the code to generate the figure, with a confidence envelope for each line. (Again, figure just provided for reference, but not integral to code to generate data for extraction)
p1 <- ggplot(AphidsIncidence, aes(x = mead_250, y = Aphid_incidence, color = Year) ) +
 geom_point()+
geom_rug(sides = "b", size = 1) +
geom_ribbon(data = newdat.lme, aes(y = NULL, ymin = lower, ymax = upper,
           color = NULL, fill = Year),
alpha = .15) +
geom_line(data = newdat.lme, aes(y = predlme), size = .75)+
theme_classic()+xlab('Proportion of meadows at 250-m')+ylab('Proportion of plants infested by aphids')

# Save ggplot figure as png
ggsave("Fig.4a.png", plot = p1, scale=0.5)

#Third set of data to EXTRACT. This would allow someone to recreate FIG. 4a based on the model predictions.
PredictedValuesAphid_incidence <- newdat.lme
PredictedValuesAphid_incidence 
# We would like to extract the data from all rows and columns, EXCEPT for the column labeled Aphid_Incidence as this is just the raw data. Description of data in each column: 
#Year: sampling year (either 2014 or 2015)
#mead_250: same as described above
#predlme-: predicted (based on model) percentage of plants in a field that would have 10 or more aphids on it 
#lower: lower limit of the 95% CI
#upper: Upper limit of the 95% CI

###End of script


####################################### 
############### ORKG ##################
####################################### 

orkg <- ORKG(host="https://incubating.orkg.org/")

# Template 'Model Fitting 4' allows for has_input_dataset as Table as opposed to URI
orkg$templates$materialize_template(template_id = "R479769")
tp = orkg$templates$list_templates()

instance <- tp$model_fitting_4(
  label="Relationship between the proportion of meadows around the experimental fields and aphid incidence (250 m radius)", 
  
  # Git Repo is currently set to private
  #has_input_dataset="https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Landscape%20affects%20pest%20and%20crop%20yield_2023.csv?",
  
  # LandscapeData can be used instead of URI
  has_input_dataset=tuple(LandscapeData, "Suitable title for dataset" ),
  
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
  

  #Predicted model values to recreate figure 4a (based on model predictions, not the raw data)
  has_output_dataset= tuple(PredictedValuesAphid_incidence, 'Dataset title 1'),

  #Output from ANOVA (Type III sum of squares) 
  has_output_dataset_2= tuple(anova (fitlme,type='marginal'), 'Dataset title 2'),
  
  # Output of summary function on lme (fixed effects)
  has_output_dataset_3= tuple(sum1$tTable, 'Dataset title 3'),
  
  # PNG output from ggplot - Git Repo is currently set to private.
  has_output_figure="https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main/Fig.4a.png",
  
  # Output statement if applicable.
  has_output_statement= "Relationship between the proportion of meadows around the experimental fields and aphid incidence (250 m radius)",
  
  # Snippet is essentially a concise version of this script with redundant code removed.
  # Git Repo is currently set to private.
  has_implementation="https://raw.githubusercontent.com/SnyderLauren/Machine-Actionable-Ecology/main//Fig4a.snippet.R"
  
)
instance$serialize_to_file("article.contribution.1.json", format="json-ld")