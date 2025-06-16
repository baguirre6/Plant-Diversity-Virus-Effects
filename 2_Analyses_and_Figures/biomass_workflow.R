#Biomass Analyses
#Beatriz A. Aguirre

#Clear environment
rm(list= ls())

#Load required libraries
library(tidyverse)
library(lsmeans)
library(forcats)
library(lme4)
library(MuMIn)
library(dplyr)
library(lmerTest)

##########################################################################################
# Import Biomass Data
##########################################################################################

# #Set working directory and import final biomass data with both years
# setwd("~/Desktop/Dissertation/1_SBF/Manuscipt/Final_SBF_Data")

biomass.data <- read.csv("1_Data/biomass.data.csv", stringsAsFactors =
                           TRUE, strip.white = TRUE, na.strings =
                           c("NA", ""))
#Check data structure:
biomass.data$Year=as.factor(biomass.data$Year)

###########################
#Subset Grass biomass data:

biomass.data %>% 
  filter(PlantType =='grass') ->grass_biomass_data

#Adjust biomass for seeding density based on Loreau and Hector (2001): 
grass_biomass_data %>%
  mutate(Adj_biomass = case_when(
    FuncDiversity == 'G' ~ grass_biomass_data$Biomass_grams.m2,
    FuncDiversity == 'GF' ~ grass_biomass_data$Biomass_grams.m2*2, # multiply x2 because grasses in this plot are planted at half the density
    FuncDiversity == 'GL' ~ grass_biomass_data$Biomass_grams.m2*2, # multiply x2 because grasses in this plot are planted at half the density
    FuncDiversity == 'GLF' ~ grass_biomass_data$Biomass_grams.m2*3 # multiply x3 because grasses in this plot are planted at a third of the density
  )) -> grass_biomass_data

##########################################################################################
# Import virus data:
##########################################################################################
# setwd("~/Desktop/Dissertation/1_SBF/Manuscipt/Final_SBF_Data")

virus.data <- read.csv("1_Data/virus.data.csv", stringsAsFactors =
                         TRUE, strip.white = TRUE, na.strings =
                         c("NA", ""))

virus.data$Year <- as.factor(virus.data$Year)

#Calculate mean BYDV prevalence for each plot:
virus.data %>% 
  group_by(Year, Block, Plot, FunctionalDiversity, GrassType)%>% 
  summarise(mean.plot.infection = mean(PAV_Infection, na.rm=TRUE)) -> mean.plot.BYDV

#Calculate mean grass biomass for each plot:
grass_biomass_data %>% 
  group_by(Year, Block, Plot, FuncDiversity, planted_fg_richness, GrassType)%>% 
  summarise(tot.plot.adj = sum(Adj_biomass)) ->plot.total.grass.biomass

#join mean plot infection data to mean plot grass biomass data by the "Year" and "Plot" columns:
plot.total.grass.biomass %>%
  left_join(mean.plot.BYDV %>% select(Year, Plot, mean.plot.infection), by = c("Year", "Plot")) %>% 
  select(-c(Block.x, FunctionalDiversity)) ->grass.regression.df

# tidy up df: 
grass.regression.df %>% 
  rename(Block = Block.y) %>% 
  select(Year, Block, Plot, FuncDiversity, GrassType, planted_fg_richness, tot.plot.adj, mean.plot.infection) -> grass.regression.df

rm(plot.total.grass.biomass)

# Look at the distribution of the grass biomass data:
hist(grass.regression.df$tot.plot.adj)
# The data does not look normally distributed --> Log transform

hist(log(grass.regression.df$tot.plot.adj))
# Log transformed data is normally distributed

#######################################################################################################
#### MODEL GRASS BIOMASS DATA -- TABLE 3 
#######################################################################################################

#Grass Productivity
grass.model <- lmer(log(tot.plot.adj) ~  mean.plot.infection * GrassType * FuncDiversity + Year + 
                           (1|Year:Block), data= grass.regression.df)

summary(grass.model)
anova(grass.model, type=3) 

r.squaredGLMM(grass.model) # Grass Productivity Marginal R2 provided in Table 3


#######################################################################################################
#### MODEL COMMUNITY BIOMASS DATA -- TABLE 3 
#######################################################################################################

#Calculate total plot biomass: 
biomass.data %>% 
  filter(Species!= 'other') %>%  #exclude weeds
  group_by(Year, Block, Plot, FuncDiversity, GrassType, planted_fg_richness) %>% 
  summarize(
    total_plot_biomass = sum(Biomass_grams.m2)) %>%
  ungroup() -> plot_biomass

#join mean plot infection data to mean plot grass biomass data by the "Year" and "Plot" columns:
plot_biomass %>%
  left_join(mean.plot.BYDV %>% select(Year, Plot, mean.plot.infection), by = c("Year", "Plot")) %>% 
  select(-c(Block.x, FunctionalDiversity)) %>% 
  rename(Block = Block.y)->comm.regression.df

#Check data frame structure: make planted fg richness factor
comm.regression.df$planted_fg_richness <- as.factor(comm.regression.df$planted_fg_richness)

################################

# Look at the distribution of the community biomass data:
hist(comm.regression.df$total_plot_biomass)
# The data does not look normally distributed --> Log transform

hist(log(comm.regression.df$total_plot_biomass))
# Log transformed data is normally distributed

#Table 3 -- Community Productivity
community.model <- lmer(log(total_plot_biomass) ~  mean.plot.infection * GrassType * FuncDiversity + Year + 
                               (1|Year:Block), data= comm.regression.df)
summary(community.model)
anova(community.model, type=3)

r.squaredGLMM(community.model) # Community Productivity Marginal R2 provided in Table 3

################################

#Post-hoc test for 3-way interaction (mean.plot infection x Functional group diversity x Grass Type)
emtrends(community.model, ~ FuncDiversity * GrassType, var = "mean.plot.infection", infer=TRUE)


# FIGURE 4

ggplot(comm.regression.df, aes(x=mean.plot.infection, y=log(total_plot_biomass))) +
  geom_point() +
  xlab("BYDV-PAV Prevalence in Grasses") +
  ylab(bquote('Log Community Productivity' ~ (g/m^-2))) +
  theme_bw() + 
  theme(axis.title = element_text(size=16),
        axis.text  = element_text(size=13)) +
  facet_grid(GrassType~FuncDiversity) -> comm.regression.fig

comm.regression.fig


# Add a regression line to only one specific facet panel, using subset() to filter data
comm.regression.fig +
  geom_smooth(data = subset(comm.regression.df, GrassType == "C3" & FuncDiversity == "G"),
              color = "#E69F00", method = "lm", se = TRUE) +
  geom_smooth(data = subset(comm.regression.df, GrassType == "C3" & FuncDiversity == "GL"),
              linetype = "dashed", color = "#E69F00", method = "lm", se = TRUE) +
  geom_text(data = data.frame(
    mean.plot.infection = c(0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.5, 0.5,
                            0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.5, 0.5),
    total_plot_biomass = c(7.50, 7.30, 7.50, 7.30, 7.50, 7.30, 7.5, 7.3,
                           7.50, 7.30, 7.50, 7.30, 7.50, 7.30, 7.5, 7.3),
    GrassType = c("C3", "C3", "C3", "C3", "C3", "C3", "C3", "C3",
                  "C4", "C4", "C4", "C4", "C4", "C4", "C4", "C4"),
    FuncDiversity =c("G", "G", "GL", "GL", "GF", "GF", "GLF", 'GLF',
                     "G", "G", "GL", "GL", "GF", "GF", "GLF", 'GLF'),
    label = c("~'t=4.340'", "~ 'P<0.0001***'", "'t=1.941'", "~'P=0.06'", "t== -0.251", "~'P=0.80'","t== -0.729", "~'P=0.46'",
              "~'t= -0.319'", "~ 'P=0.75'", "'t= -1.081'", "~'P=0.28'", "t== 0.246", "~'P=0.81'","t== -0.356", "~'P=0.72'")
  ),
  aes(x = mean.plot.infection, y = total_plot_biomass, label = label),
  parse = TRUE) -> figure.4

figure.4
#ggsave("~/Desktop/Dissertation/1_SBF/Manuscipt/2025_SubmissionMaterials/Analyses_March2025/Figures/Figure3_CommProductivity.pdf", width = 7, height = 5)

#######################################################################################################
#### VISUALIZE GRASS BIOMASS DATA -- FIGURE 3
#######################################################################################################

#Post-hoc test for mean.plot infection x Grass Type interaction
emtrends(grass.model, pairwise~GrassType, var = "mean.plot.infection", infer=TRUE)

# FIGURE 4

ggplot(grass.regression.df, aes(x=mean.plot.infection, y=log(tot.plot.adj))) +
  geom_point() +
  xlab("BYDV-PAV Prevalence in Grasses") +
  ylab(bquote('Log Grass Productivity' ~ (g/m^-2))) +
  theme_bw() + 
  theme(axis.title = element_text(size=16),
           axis.text  = element_text(size=13)) + 
  facet_grid(~GrassType) -> grass.regression.fig


# Add a regression line to only one specific facet panel, using subset() to filter data
grass.regression.fig +
  geom_smooth(data = subset(grass.regression.df, GrassType == "C3"),
              color = "#E69F00", method = "lm", se = TRUE) +
  geom_text(data = data.frame(
    mean.plot.infection = c(0.59, 0.59, 0.59, 0.59, 0.59),
    tot.plot.adj = c(7.50, 7.30, 7.10, 7.30, 7.10),
    GrassType = c("C3", "C3", "C3", "C4", "C4"),
    label = c("R^2==0.58", "~'t=3.578'", "~ 'P=0.0005***'","~'t= -1.174'", "~ 'P=0.245'")
  ),
  aes(x = mean.plot.infection, y = tot.plot.adj, label = label),
  parse = TRUE) -> figure.5

figure.5
#ggsave("~/Desktop/Dissertation/1_SBF/Manuscipt/2025_SubmissionMaterials/Analyses_March2025/Figures/Figure3_GrassProductivity.pdf", width = 6, height = 5)


#######################################################################################################
#### TOTAL FORB AND LEGUME BIOMASS ANALYSES 
#######################################################################################################

#Calculate total plant type biomass at the plot level
biomass.data %>% 
  filter(Species!= 'other') %>%  #remove other species
  group_by(Year, Block, Plot, GrassType, FuncDiversity, PlantType) %>% 
  summarize(plot.biomass = sum(Biomass_grams.m2)) -> plant.type.plot.biomass

#Add mean plot infection data to total.plot.community.biomass df
plant.type.biomass <- left_join(plant.type.plot.biomass, mean.plot.BYDV, by = c("Year", "Block", "Plot"))

plant.type.biomass %>% 
  select(-GrassType.x) %>% 
  rename(GrassType = GrassType.y) -> plant.type.biomass

#Clean up dataframe
plant.type.biomass %>% 
  select(Year, Block, Plot, GrassType, FuncDiversity, PlantType, plot.biomass, mean.plot.infection) -> plant.type.biomass

#########################################################
#### Subset net forb data
plant.type.biomass %>% 
  filter(PlantType == 'forb') -> net.forb.biomass

net.forb.biomass %>% 
  mutate(Adj_biomass = case_when(
    FuncDiversity == 'GF' ~ plot.biomass,
    FuncDiversity == 'GLF' ~ plot.biomass*1.5 # multiply x1.5 to standardize forb seeding density in GLF with GF
  )) -> adjusted.forb.data

hist(log(adjusted.forb.data$plot.biomass))

forb.model.4.add <- lmer(log(plot.biomass) ~  mean.plot.infection * GrassType * FuncDiversity + Year + 
                           (1|Year:Block), data= adjusted.forb.data)

anova(forb.model.4.add)
emmeans(forb.model.4.add, ~ Year, type="response")

r.squaredGLMM(forb.model.4.add)

#Run model diagnostics:
#Plot residuals vs. fitted values to check for heteroskedasticity
plot(fitted(forb.model.4.add), resid(forb.model.4.add), 
     xlab = "Fitted Values", ylab = "Residuals", 
     main = "Residuals vs. Fitted Values")

abline(h = 0, col = "red")


#########################################################
#### Subset net legume data

plant.type.biomass %>% 
  filter(PlantType == 'legume') -> net.legume.biomass

net.legume.biomass %>% 
  mutate(Adj_biomass = case_when(
    FuncDiversity == 'GL' ~ plot.biomass,
    FuncDiversity == 'GLF' ~ plot.biomass*1.5 # multiply x1.5 to standardize seeding density in GLF with GL
  )) -> adjusted.legume.data

adjusted.legume.data %>% 
  mutate(log.plot.biomass = log1p(plot.biomass), # log1p(x) computes the log (x +1) this is necessary since there are zeros in the data
         sq.plot.biomass = sqrt(plot.biomass)
  ) -> adjusted.legume.data

hist(adjusted.legume.data$plot.biomass)
hist(adjusted.legume.data$log.plot.biomass)
hist(adjusted.legume.data$sq.plot.biomass) #Data follows normal distribution; use this for the model

legume.model.sq.trans <- lmer(sq.plot.biomass ~  mean.plot.infection * GrassType * FuncDiversity + Year + 
                       (1|Year:Block), data= adjusted.legume.data)

summary(legume.model.sq.trans)
anova(legume.model.sq.trans, type=3)

#Run model diagnostics:
#Plot residuals vs. fitted values to check for heteroskedasticity
plot(fitted(legume.model.sq.trans), resid(legume.model.sq.trans), 
     xlab = "Fitted Values", ylab = "Residuals", 
     main = "Residuals vs. Fitted Values")

abline(h = 0, col = "red")

# legume.model.sq.trans is a good model. 
r.squaredGLMM(legume.model.sq.trans)

emmeans(legume.model.sq.trans, pairwise~ FuncDiversity, type="response") #significant


#####################################################################
# VISUALIZE LEGUME PRODUCTIVITY 
#####################################################################
# legume.emm =as.data.frame(summary(emmeans(legume.model.sq.trans, ~FuncDiversity, type="response")))
# legume.emm
# 
# legume.figure<-ggplot(legume.emm,
#                 aes(x=FuncDiversity, y=emmean), size=3) +
#   geom_point(size=5) +
#   geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.05, position=position_dodge(0.01)) +
#   ylab("Net Legume Productivity") +
#   xlab("Functional Diversity") +
#   theme_bw() +
#   theme(axis.title   = element_text(face = "bold", size=15),
#         axis.text    = element_text(face = "bold", size=13))
# 
# legume.figure


