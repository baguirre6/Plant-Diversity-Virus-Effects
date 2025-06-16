#Total hosts infected with Virus Analyses
#June 09, 2025
#Author: Beatriz A. Aguirre

######################################
#Clear environment
rm(list= ls())

library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(glmmTMB)
library(ggpubr)
library(forcats)
library(glmmTMB)

##########################################################################################
#Import BYDV infection data 
##########################################################################################

virus.data <- read.csv("1_Data/virus.data.csv", stringsAsFactors =
                         TRUE, strip.white = TRUE, na.strings =
                         c("NA", ""))

#fix data structure
virus.data$Year<-as.factor(virus.data$Year)

#add fg_richness column to infection data
virus.data %>% 
  mutate(fg_richness = case_when(
    FunctionalDiversity == 'G' ~ '1', 
    FunctionalDiversity == 'GF' ~ '2',
    FunctionalDiversity == 'GL' ~ '2',
    FunctionalDiversity == 'GLF' ~ '3',
  )) ->
  virus.data

virus.data$fg_richness <- as.factor(virus.data$fg_richness)

#Calculate virus prevalence (proportion infected) for each species by plot:
virus.data %>% 
  group_by(Year, Block, Plot, FunctionalDiversity, GrassType, fg_richness, Species)%>% 
  summarise(tot.infected = sum(PAV_Infection, na.rm=TRUE)) -> tot.infected.df

###################################################################################
#Total infected with virus by Species Analyses 
###################################################################################
#### Virus Model By Species & Functional Diversity
sp.virus.prev.model = glmmTMB(tot.infected ~ FunctionalDiversity * Species + Year #removed grass type variable because it's redundant with species
                              + (1|Year:Block) + (1|Year:Block:Plot),
                              data=tot.infected.df)
summary(sp.virus.prev.model)
car::Anova(sp.virus.prev.model, type=3) #Same results as infection probability and virus prevalence (proportion infected)

emmeans(sp.virus.prev.model, pairwise~Species, type="response")
emm <-emmeans(sp.virus.prev.model, pairwise~Species, type="response")

multcomp::cld(emm, Letters = LETTERS) #compact letter display for species

#############################################################################
# Total infected with virus Analyses 
#############################################################################
virus.prev.model = glmmTMB(tot.infected ~GrassType * FunctionalDiversity + Year
                           + (1|Year:Block) + (1|Year:Block:Plot) + (1|Species),
                           data=tot.infected.df)
summary(virus.prev.model)
car::Anova(virus.prev.model,type=3) #Same results as infection probability and virus prevalence (proportion infected)


# #Post-hoct tests for effects of functional diversity on BYDV Prevalence (proportion infected)
emm_FD <-emmeans(virus.prev.model, ~FunctionalDiversity, type="response")
emm_FD
multcomp::cld(emm_FD, Letters = LETTERS) #compact letter display for functional diversity treatments
contrast(emm_FD, "trt.vs.ctrl", ref = "G")
#Post-hoc tests indicate no difference between functional diversity treatments

#############################################################################
# Total infected with virus -- FG Richness Analyses 
#############################################################################

virus.prev.richness.model = glmmTMB(tot.infected ~GrassType * fg_richness + Year
                                    + (1|Year:Block) + (1|Year:Block:Plot) + (1|Species),
                                    data=tot.infected.df)

summary(virus.prev.richness.model)
car::Anova(virus.prev.richness.model, type=3) #Same results as infection probability and virus prevalence (proportion infected)


