#Virus Analyses and Figures

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

###################################################
#Virus by Species Analyses - Included in Results
###################################################
#### Virus Model By Species & Functional Diversity
sp.virus.model = glmmTMB(PAV_Infection ~ FunctionalDiversity * Species + Year #removed grass type variable because it's redundant with species
                         + (1|Year:Block) + (1|Year:Block:Plot),
                         family='binomial', data=virus.data)
summary(sp.virus.model)
car::Anova(sp.virus.model, type=3)

emmeans(sp.virus.model, pairwise~Species, type="response")
emm <-emmeans(sp.virus.model, pairwise~Species, type="response")

multcomp::cld(emm, Letters = LETTERS) #compact letter display for figure 1

#######################
# BYDV X Species -- FIGURE 1
species.info0=as.data.frame(summary(emmeans(sp.virus.model, ~Species, type="response")))
species.info0

#add grass type column
species.info0 %>% 
  mutate(GrassType =
           case_when(
             Species == 'Avena fatua' ~ 'C3',
             Species == 'Lolium multiflorum' ~ 'C3',
             Species == 'Secale cereale' ~ 'C3',
             Species == 'Echinochloa crusgalli' ~ 'C4',
             Species == 'Setaria italica' ~ 'C4',
             Species == 'Panicum miliaceum' ~ 'C4'
           )
  ) -> species.info0

species.info0$GrassType <- as.factor(species.info0$GrassType)

species.info0 %>% 
  mutate(Species = fct_relevel(Species,
                               "Avena fatua", "Lolium multiflorum", "Secale cereale", 
                               "Echinochloa crusgalli", "Panicum miliaceum", "Setaria italica")) -> species.info0

levels(species.info0$Species)

library(wesanderson)

Figure1_bydv.species <-ggplot(species.info0,
                       aes(x=Species, y=prob)) +
  geom_point(aes(color=GrassType), size=5) +
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE), width=.05, position=position_dodge(0.01)) +
  ylab("BYDV-PAV Prevalence") +
  xlab("Grass Species") +
  ylim(c(0,.5)) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold", size=15),
        axis.text.y    = element_text(face = "bold", size=13),
        axis.text.x = element_text(face = "bold.italic", size=13, angle = 45, hjust = 1),
        legend.position = "top") +
  annotate("text", x="Avena fatua", y=.36, label= 'c',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Lolium multiflorum", y=.15, label= 'a',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Secale cereale", y=.38, label= 'c',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Echinochloa crusgalli", y=.20, label= 'ab',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Panicum miliaceum", y=.15, label= 'a',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Setaria italica", y=.23, label= 'b',
           col="black", size=6, parse=TRUE)

Figure1_bydv.species +scale_color_manual(guide = guide_legend(title = "Grass Type"), values=wes_palette(n=2, name="GrandBudapest1"))


###################################################
#Virus Analyses - Table 2
###################################################
library(glmmTMB)

virus.model = glmmTMB(PAV_Infection ~GrassType * FunctionalDiversity + Year
                                + (1|Year:Block) + (1|Year:Block:Plot) + (1|Species),
                                family='binomial', data=virus.data)
summary(virus.model)
car::Anova(virus.model, type=3)

###################################################
#FG Richness X Virus Analyses - Table S2
###################################################

virus.richness.model = glmmTMB(PAV_Infection ~GrassType * fg_richness + Year 
                                  + (1|Year:Block) + (1|Year:Block:Plot) + (1|Species),
                                  family='binomial', data=virus.data)

summary(virus.richness.model)
car::Anova(virus.richness.model, type=3)

######################################
#Virus Figure 2
######################################

#Visualize BYDV x year:
yearly.bydv =emmeans(virus.model, ~Year, type="response")
yearly.bydv

result.annual.bydv = as.data.frame(summary(yearly.bydv)) #Export emmeans data into a data frame for plotting w/ggplot

panel.1<-ggplot(result.annual.bydv,
                aes(x=Year, y=prob, group=Year)) +
  geom_point(aes(shape=Year, color=Year), size=5) +
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE, color=Year), width=.05, position=position_dodge(0.8)) +
  ylab("BYDV-PAV\n Prevalence") +
  xlab("Year") +
  ylim(c(0,.5)) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold", size=15),
        axis.text    = element_text(face = "bold", size=13),
        legend.position = "none",
        plot.caption = element_text(hjust = 0, size=10)) +
  geom_signif(comparisons = list(c("2021", "2022")), annotations = "*", map_signif_level=TRUE,
              textsize=8, y_position = 0.35, tip_length = 0.05, vjust=0.4) +
  annotate("text", x=1.5, y=0.42, label= '"P=0.006"',
           col="black", size=4, parse=TRUE) +
  annotate("text", 0.55, y=.50, label= '"(a)"',
           col="black", size=6, parse=TRUE)

panel.1

#######################

#Create panel figure for Grass Type

panel_gt=as.data.frame(summary(emmeans(virus.model, ~GrassType+Year, type="response")))
panel_gt

panel.2<-ggplot(panel_gt,
                aes(x=GrassType, y=prob, group=Year), size=3) +
  geom_point(aes(shape=Year, color=Year), size=5) +
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE, color=Year), width=.05, position=position_dodge(0.01)) +
  ylab("Probability of Infection") +
  xlab("Grass Type") +
  ylim(c(0,.5)) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold", size=15),
        axis.text    = element_text(face = "bold", size=13),
        plot.caption = element_text(hjust = 0, size=10),
        axis.title.y = element_blank(),
        legend.position = c(.85,.82)) +
  annotate("text", x=0.55, y=.50, label= '"(b)"',
           col="black", size=6, parse=TRUE) +
  geom_signif(comparisons = list(c("C3", "C4")), annotations = "", map_signif_level=TRUE,
              textsize=10, y_position = 0.40, tip_length = 0.05, vjust=0.4) +
  annotate("text", x=1.5, y=0.45, label= 'NS',
           col="black", size=4, parse=TRUE)

panel.2

#######################
#Visualize BYDV x functional group richness

marginal.2 =as.data.frame(summary(emmeans(virus.richness.model, ~fg_richness+Year, type="response")))
marginal.2

panel.3<-ggplot(marginal.2,
                aes(x=as.factor(fg_richness), y=prob, group=Year)) +
  geom_point(aes(shape=Year, color=Year), size=5) +
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE, color=Year), width=.05, position=position_dodge(0)) +
  ylab("BYDV-PAV\n Prevalence") +
  xlab("Functional Group Richness") +
  ylim(c(0,.5)) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold", size=15),
        axis.text    = element_text(face = "bold", size=13),
        legend.position = "none") +
  annotate("text", x=.60, y=.50, label= '"(c)"',
           col="black", size=6, parse=TRUE)  +
  geom_signif(comparisons = list(c("1", "3")), annotations = "", map_signif_level=TRUE,
              textsize=10, y_position = 0.40, tip_length = 0.05, vjust=0.4) +
  annotate("text", x=2, y=0.45, label= 'NS',
           col="black", size=4, parse=TRUE)

panel.3

#######################

#Visualize BYDV x functional diversity
panel_fd=as.data.frame(summary(emmeans(virus.model, ~FunctionalDiversity+Year, type="response")))
panel_fd

panel.4 <-ggplot(panel_fd,
                 aes(x=FunctionalDiversity, y=prob, group=Year)) +
  geom_point(aes(shape=Year, color=Year), size=5) +
  geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE, color=Year), width=.05, position=position_dodge(0.01)) +
  ylab("Probability of Infection") +
  xlab("Functional Group Composition") +
  ylim(c(0,.5)) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold", size=15),
        axis.text    = element_text(face = "bold", size=13),
        plot.caption = element_text(hjust = 0, size=10),
        axis.title.y = element_blank(),
        legend.position = "none") +
  annotate("text", x=0.65, y=.50, label= '"(d)"',
           col="black", size=6, parse=TRUE) +
  geom_signif(comparisons = list(c("G", "GLF")), annotations = "", map_signif_level=TRUE,
              textsize=10, y_position = 0.40, tip_length = 0.05, vjust=0.4) +
  annotate("text", x=2.5, y=0.45, label= 'NS',
           col="black", size=4, parse=TRUE)

panel.4

#####################################################################
# Join panels to make final BYDV Figure -- Figure 2

BYDV_figure2 <- ggarrange(panel.1, panel.2, panel.3, panel.4,
                                     ncol = 2, nrow = 2)
BYDV_figure2

#####################################################################
