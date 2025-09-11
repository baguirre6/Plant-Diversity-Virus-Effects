# Supplemental Analyses and Figures
# Beatriz A. Aguirre

rm(list = ls())

library(dplyr)
library(ggplot2)

#####################################################################
# Figure S2 -Seed density by FG Richness
#####################################################################

seed.density.data2 <- data.frame(
  fg_richness = c(1, 2, 3), 
  seed_density = c(2499, 1251, 834),
  stringsAsFactors = FALSE
)


seed.density.data2 %>%
  ggplot(aes(x = fg_richness, y = seed_density)) +
  geom_col() +
  xlab("Functional Group Richness") +
  ylab(bquote('Grass seed density per' ~ m^-2)) +
  theme_bw() +
  theme(axis.title   = element_text(size=17),
        axis.text    = element_text(size=13)) -> figure.S2

figure.S2

# ggsave("4_Figures/FigureS2.pdf", figure.S2, 
#        width = 6, height = 5)

#####################################################################
# Table S2 and Figure S3 (MODEL GRASS SPECIES BIOMASS DATA)
#####################################################################
source("2_Analyses_and_Figures/biomass_workflow.R")

#calculate mean infection for each species in each plot:
virus.data %>% 
  group_by(Year, Block, Plot, Species) %>% 
  summarise(mean.sp.infection = mean(PAV_Infection, na.rm=TRUE)) -> mean.sp.plot.BYDV

#join mean plot infection data to mean plot grass biomass data by the "Year" and "Plot" columns:
grass_biomass_data %>%
  left_join(mean.sp.plot.BYDV %>%
              select(Year, Plot, Species, mean.sp.infection), 
            by = c("Year",  "Block", "Plot", "Species")) -> sp.grass.data

hist(sp.grass.data$Adj_biomass) #right skewed
hist(log(sp.grass.data$Adj_biomass)) #log transformed data distribution is normal

grass.spp.model <- lmer(log(Adj_biomass) ~  mean.sp.infection * Species * FuncDiversity + Year + 
                          (1|Year:Block), data= sp.grass.data)

summary(grass.spp.model)
anova(grass.spp.model, type=3) 

r.squaredGLMM(grass.spp.model)

# #Post-hoc trend for Functional Diversity
# emmeans(grass.spp.model, pairwise~FuncDiversity, type="response")
emtrends(grass.spp.model, ~ Species, var="mean.sp.infection", infer=TRUE)

#Post-hoc trend for mean.sp.infection X Species
emtrends(grass.spp.model, ~ Species, var="mean.sp.infection", infer=TRUE)

species_labels <- c(
  "Echinochloa crusgalli" = "italic('Echinochloa crusgalli')",
  "Panicum miliaceum" = "italic('Panicum miliaceum')",
  "Secale cereale" = "italic('Secale cereale')",
  "Avena fatua" = "italic('Avena fatua')",
  "Lolium multiflorum" = "italic('Lolium multiflorum')",
  "Setaria italica" = "italic('Setaria italica')",
  "C3" = "C3",
  "C4" = "C4"
)

library(wesanderson)

# Plot regressions of species biomass by infection
sp.grass.data %>% 
  mutate(Species = fct_relevel(Species,
                               "Avena fatua", "Lolium multiflorum", "Secale cereale", 
                               "Echinochloa crusgalli", "Panicum miliaceum", "Setaria italica")) %>% 
  ggplot(aes(x=mean.sp.infection, y=log(Adj_biomass), group=Species)) +
  geom_point(aes(color=GrassType), alpha=0.3) +
  xlab("BYDV-PAV Prevalence") +
  ylab(bquote('Log Productivity' ~ (g/m^-2))) +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=13),
        legend.position = "none") +
  facet_wrap(~Species, labeller = as_labeller(species_labels, label_parsed))-> figure5c.base

figure5c.base

figure5c.base +
  geom_text(data = data.frame(
    mean.sp.infection = c(0.75, 0.75, 0.75, 0.75, 0.75, 0.75),
    Adj_biomass = c(1500, 1500, 1500, 1500, 1500, 1500),
    GrassType = factor(c("C4", "C4", "C3", "C4", "C3", "C3"), levels = c("C3", "C4")),
    Species = factor(c("Echinochloa crusgalli", "Panicum miliaceum", "Secale cereale", "Setaria italica", "Avena fatua", "Lolium multiflorum"), 
                     levels = c("Avena fatua", "Lolium multiflorum", "Secale cereale",
                                "Echinochloa crusgalli", "Panicum miliaceum", "Setaria italica")),
    label = c("~'t= -1.784 P=0.07'", "~'t= 3.451 P=0.0007*'", "~'t= 2.736 P=0.006*'","~'t= 1.443 P=0.15'",
              "~'t= 0.049 P=0.96'", "~'t=0.516 P=0.60'")
  ),
  aes(x = mean.sp.infection, y = log(Adj_biomass), label = label), 
  parse = TRUE) +
  geom_smooth(data = subset(sp.grass.data, GrassType == "C4"  & Species == "Echinochloa crusgalli"),
              linetype = "dashed", color = "#FD6467", method = "lm", se = TRUE) +
  geom_smooth(data = subset(sp.grass.data, GrassType == "C4" & Species == "Panicum miliaceum"), 
              color = "#FD6467", method = "lm", se = TRUE) +
  geom_smooth(data = subset(sp.grass.data, GrassType == "C3"  & Species == "Secale cereale"),
              color = "#FFB90F", method = "lm", se = TRUE) +
              scale_color_manual(guide = guide_legend(title = "Grass Type"), values=wes_palette(n=2, name="GrandBudapest1")) ->figure5c 
figure5c

###################################################################################################
########## Plot biomass by species (Figure 5b): 
emmeans(grass.spp.model, pairwise~Species)
emm_sp_biomass <-emmeans(grass.spp.model, pairwise~Species)

multcomp::cld(emm_sp_biomass, Letters = LETTERS) #compact letter display for figure 5b

#######################

calc_SE <- function(x) {
  x2<-na.omit(x)
  sd(x2)/sqrt(length(x2))
}

sp.grass.data %>% 
  group_by(Species, GrassType) %>% 
  mutate(Species = fct_relevel(Species,
                               "Avena fatua", "Lolium multiflorum", "Secale cereale", 
                               "Echinochloa crusgalli", "Panicum miliaceum", "Setaria italica")) %>% 
  summarise(mean.prod = mean(log(Adj_biomass)), se.prod= calc_SE(log(Adj_biomass))) %>% 
  ggplot(aes(x=Species, y=mean.prod)) +
  geom_point(aes(color=GrassType), size=4) +
  geom_errorbar(aes(ymin=mean.prod-se.prod, ymax=mean.prod+se.prod), width=.05, position=position_dodge(0.01)) +
  xlab("Grass Species") +
  ylab(bquote('Log Productivity' ~ (g/m^-2))) + 
  ylim(3.5, 6.0) +
  theme_bw() +
  theme(axis.title   = element_text(size=16),
          axis.text.y    = element_text(size=13),
          axis.text.x = element_text(face = "italic", size=13, angle = 45, hjust = 1),
          legend.position = "top") +
  annotate("text", x="Avena fatua", y=5.5, label= 'd',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Lolium multiflorum", y=4.2, label= 'a',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Secale cereale", y=4.8, label= 'ab',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Echinochloa crusgalli", y=5.25, label= 'bc',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Panicum miliaceum", y=4.8, label= 'b',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Setaria italica", y=5.4, label= 'cd',
           col="black", size=6, parse=TRUE) -> Figure5b_base

Figure5b_base +
  scale_color_manual(guide = guide_legend(title = "Grass Type"), values=wes_palette(n=2, name="GrandBudapest1")) -> figure5b

figure5b

#import species infection figure:
source("2_Analyses_and_Figures/virus_workflow.R")

#####################################################################
# Join panels to make final Figure 5 
Figure5 <- ggarrange(
  ggarrange(Fig2.panela, figure5b, 
            ncol = 2, labels = c("(a)", "(b)"), 
            font.label = list(size = 16), 
            common.legend = TRUE, legend = "top"),
  figure5c,
  ncol = 1, heights = c(1, 1),
  labels = c("", "(c)"),
  font.label = list(size = 16))

Figure5

# ggsave("4_Figures/Figure5.pdf", Figure5,
#        width = 12, height = 12)

#####################################################################
# Plot species biomass by infection
sp.grass.data %>% 
  group_by(Species, GrassType, planted_fg_richness) %>% 
  mutate(Species = fct_relevel(Species,
                               "Avena fatua", "Lolium multiflorum", "Secale cereale", 
                               "Echinochloa crusgalli", "Panicum miliaceum", "Setaria italica")) %>% 
  summarise(mean.prod = mean(log(Adj_biomass)), se.prod= calc_SE(log(Adj_biomass))) %>% 
  ggplot(aes(x=as.factor(planted_fg_richness), y=mean.prod, group=Species)) +
  geom_point(aes(color=GrassType)) +
  geom_errorbar(aes(ymin=mean.prod-se.prod, ymax=mean.prod+se.prod, color=GrassType), width=.05, position=position_dodge(0.01)) +
  xlab("Functional Group Richness") +
  ylab(bquote('Log Productivity' ~ (g/m^-2))) +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=13),
        legend.position = "top") +
  facet_wrap(~Species, labeller = as_labeller(species_labels, label_parsed))-> figure.S4

figure.S4 +
  scale_color_manual(guide = guide_legend(title = "Grass Type"), values=wes_palette(n=2, name="GrandBudapest1")) -> figure.S4

figure.S4

# ggsave("4_Figures/FigureS4.pdf", figure.S4, 
#        width = 6, height = 5)

#####################################################################
# Supplemental analyses to check if biomass is positively 
# correlated with infection: 
#####################################################################
grass.model.sup <- lmer(mean.plot.infection ~ log(tot.plot.adj) * GrassType * FuncDiversity + Year + 
                          (1|Year:Block), data= grass.regression.df)

anova(grass.model.sup, type=3) 

#####################################################################
# Biomass by Richness Analyses (TABLE S5)
#####################################################################
#Table S5 -- Net Community Productivity
community.model.richness <- lmer(log(total_plot_biomass) ~  mean.plot.infection * GrassType * planted_fg_richness  + Year + 
                          (1|Year:Block), data= comm.regression.df)
summary(community.model.richness)
anova(community.model.richness, type=3) #Table S5 

r.squaredGLMM(community.model.richness) # Community Productivity Marginal R2 provided in Table S5

#Post-hoc tests for significant 3-way interaction (mean.plot infection x Functional group diversity x Grass Type):
emtrends(community.model.richness, ~ planted_fg_richness | GrassType, var = "mean.plot.infection", infer=TRUE)

##############
#Table S5 -- Net Grass Productivity
grass.model.richness <- lmer(log(tot.plot.adj) ~  mean.plot.infection * GrassType * planted_fg_richness + Year + 
                      (1|Year:Block), data= grass.regression.df)

summary(grass.model.richness)
anova(grass.model.richness, type=3) #Table S5 

r.squaredGLMM(grass.model.richness) # Grass Productivity Marginal R2 provided in Table 3

#Post-hoc tests for significant effect of Grass Type):
emmeans(grass.model.richness, ~GrassType)
