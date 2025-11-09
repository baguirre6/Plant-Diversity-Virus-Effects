#Virus Prevalence (proportion infected) Analyses
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
library(googledrive)

##########################################################################################
#Import BYDV infection data 
##########################################################################################

#drive authentication
drive_auth()

#create data folder
dir.create(file.path("data"), showWarnings = F)

#identify desired file
focal_file <- "virus.data.csv"

#download data file
googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/1dpjkQh9MKMO8LInw4A1IaRXhddenbFzq")) %>% 
  dplyr::filter(name == focal_file) %>% 
  googledrive::drive_download(file = .$id, overwrite = T,
                              path = file.path("data", .$name))

# Read in infection data
virus.data <- read.csv(file = file.path("data", "virus.data.csv"))

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
  summarise(mean.plot.infection = mean(PAV_Infection, na.rm=TRUE)) -> mean.spp.BYDV

#Add weights column
mean.spp.BYDV$weights <- 10

#Adjust weights for species with less than 10 samples :
mean.spp.BYDV$weights[90] <- 7
mean.spp.BYDV$weights[127] <- 7
mean.spp.BYDV$weights[133] <- 3
mean.spp.BYDV$weights[141] <- 8
mean.spp.BYDV$weights[160] <- 3

###################################################################################
#Virus prevalence (proportion infected) by Species Analyses 
###################################################################################
#### Virus Model By Species & Functional Diversity
sp.virus.prev.model = glmmTMB(mean.plot.infection ~ FunctionalDiversity * Species + Year #removed grass type variable because it's redundant with species
                         + (1|Year:Block) + (1|Year:Block:Plot),
                          data=mean.spp.BYDV)
summary(sp.virus.prev.model)
car::Anova(sp.virus.prev.model, type=3) #TABLE S2 (BYDV-PAV INFECTION PREVALENCE)

# emmeans(sp.virus.prev.model, pairwise~Species, type="response")
# emm <-emmeans(sp.virus.prev.model, pairwise~Species, type="response")
# 
# multcomp::cld(emm, Letters = LETTERS) #compact letter display for species


#############################################################################
# Virus Prevalence (proportion infected) Analyses 
#############################################################################
virus.prev.model = glmmTMB(mean.plot.infection ~GrassType * FunctionalDiversity + Year
                      + (1|Year:Block) + (1|Year:Block:Plot) + (1|Year:Block:Plot:Species) + (1|Species), 
                      weights=weights,
                      family= 'binomial', data=mean.spp.BYDV)
summary(virus.prev.model)
car::Anova(virus.prev.model,type=3) #TABLE S3 (BYDV-PAV INFECTION PREVALENCE)


# #Post-hoct tests for effects of functional diversity on BYDV Prevalence (proportion infected)
emm_FD <-emmeans(virus.prev.model, ~FunctionalDiversity, type="response")
emm_FD
multcomp::cld(emm_FD, Letters = LETTERS) #compact letter display for functional diversity treatments
contrast(emm_FD, "trt.vs.ctrl", ref = "G")
#Both post-hoc tests indicate no difference between functional diversity treatments 
#(this is consistent with infection probability analysis)

#############################################################################
# Virus Prevalence (proportion infected) X FG Richness Analyses 
#############################################################################

virus.prev.richness.model = glmmTMB(mean.plot.infection ~GrassType * fg_richness + Year
                               + (1|Year:Block) + (1|Year:Block:Plot) + (1|Year:Block:Plot:Species) + (1|Species),
                               weights=weights,
                               family= 'binomial', data=mean.spp.BYDV)

summary(virus.prev.richness.model)
car::Anova(virus.prev.richness.model, type=3) #TABLE S4 (BYDV-PAV INFECTION PREVALENCE)

