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
  summarise(tot.infected = sum(PAV_Infection, na.rm=TRUE)) -> tot.infected.df

###################################################################################
#Total infected with virus by Species Analyses 
###################################################################################
#### Virus Model By Species & Functional Diversity
sp.virus.ct.model = glmmTMB(tot.infected ~ FunctionalDiversity * Species + Year #removed grass type variable because it's redundant with species
                              + (1|Year:Block) + (1|Year:Block:Plot),
                              family="poisson",
                              data=tot.infected.df) #TABLE S2 (TOTAL HOSTS INFECTED)
summary(sp.virus.ct.model)
car::Anova(sp.virus.ct.model, type=3) #Similar results as infection probability and virus prevalence (proportion infected)

# emmeans(sp.virus.ct.model, pairwise~Species, type="response")
# emm <-emmeans(sp.virus.ct.model, pairwise~Species, type="response")
# 
# multcomp::cld(emm, Letters = LETTERS) #compact letter display for species

#############################################################################
# Total infected with virus Analyses 
#############################################################################
virus.ct.model = glmmTMB(tot.infected ~GrassType * FunctionalDiversity + Year
                           + (1|Year:Block) + (1|Year:Block:Plot) + (1|Year:Block:Plot:Species) + (1|Species),
                           family="poisson",
                           data=tot.infected.df)
summary(virus.ct.model)
car::Anova(virus.ct.model, type=3) #TABLE S3 (TOTAL HOSTS INFECTED)
#Similar results as infection probability and virus prevalence (proportion infected)


#############################################################################
# Total infected with virus -- FG Richness Analyses 
#############################################################################

virus.ct.richness.model = glmmTMB(tot.infected ~GrassType * fg_richness + Year
                                    + (1|Year:Block) + (1|Year:Block:Plot) +  (1|Year:Block:Plot:Species) + (1|Species),
                                    family="poisson",
                                    data=tot.infected.df)

summary(virus.ct.richness.model)
car::Anova(virus.ct.richness.model, type=3) #TABLE S4 (TOTAL HOSTS INFECTED)


