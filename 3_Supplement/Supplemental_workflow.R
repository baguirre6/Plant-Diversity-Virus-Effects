# Supplemental Analyses and Figures
# Beatriz A. Aguirre

rm(list = ls())

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


# Plot regressions of species biomass by infection
ggplot(sp.grass.data, aes(x=mean.sp.infection, y=log(Adj_biomass), group=Species)) +
  geom_point() +
  xlab("BYDV-PAV Prevalence in Grasses") +
  ylab(bquote('Log Grass Species Productivity' ~ (g/m^-2))) +
  theme_bw() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=13)) +
  facet_wrap(GrassType~Species, labeller = as_labeller(species_labels, label_parsed))-> sp.figure.1


sp.figure.1 +
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
              linetype = "dashed", color = "#00BFC4", method = "lm", se = TRUE) +
  geom_smooth(data = subset(sp.grass.data, GrassType == "C4" & Species == "Panicum miliaceum"), 
              color = "#00BFC4", method = "lm", se = TRUE) +
  geom_smooth(data = subset(sp.grass.data, GrassType == "C3"  & Species == "Secale cereale"),
              color = "#00BFC4", method = "lm", se = TRUE)

###################################################################################################
########## Plot biomass by species: 
emmeans(grass.spp.model, pairwise~Species, type="response")
emm_sp_biomass <-emmeans(grass.spp.model, pairwise~Species, type="response")

multcomp::cld(emm_sp_biomass, Letters = LETTERS) #compact letter display for figure 1

#######################
# BYDV infection probability X Species -- FIGURE 1
species.biomass.info =as.data.frame(summary(emmeans(grass.spp.model, ~Species, type="response")))
species.biomass.info

#add grass type column
species.biomass.info %>% 
  mutate(GrassType =
           case_when(
             Species == 'Avena fatua' ~ 'C3',
             Species == 'Lolium multiflorum' ~ 'C3',
             Species == 'Secale cereale' ~ 'C3',
             Species == 'Echinochloa crusgalli' ~ 'C4',
             Species == 'Setaria italica' ~ 'C4',
             Species == 'Panicum miliaceum' ~ 'C4'
           )
  ) -> species.biomass.info

species.biomass.info$GrassType <- as.factor(species.biomass.info$GrassType)

species.biomass.info %>% 
  mutate(Species = fct_relevel(Species,
                               "Avena fatua", "Lolium multiflorum", "Secale cereale", 
                               "Echinochloa crusgalli", "Panicum miliaceum", "Setaria italica")) -> species.biomass.info

levels(species.biomass.info$Species)

library(wesanderson)

Figure1b_species_biomass <-ggplot(species.biomass.info,
                       aes(x=Species, y=response)) +
  geom_point(aes(color=GrassType), size=4) +
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=.05, position=position_dodge(0.01)) +
  xlab("Grass Species") +
  ylab(bquote('Grass Species Productivity' ~ (g/m^-2))) +
  ylim(40, 350) +
  theme_bw() +
  theme(axis.title   = element_text(size=16),
        axis.text.y    = element_text(size=13),
        axis.text.x = element_text(face = "italic", size=13, angle = 45, hjust = 1),
        legend.position = "top") +
  annotate("text", x="Avena fatua", y=330, label= 'd',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Lolium multiflorum", y=115, label= 'a',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Secale cereale", y=130, label= 'ab',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Echinochloa crusgalli", y=175, label= 'bc',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Panicum miliaceum", y=150, label= 'b',
           col="black", size=6, parse=TRUE) +
  annotate("text", x="Setaria italica", y=250, label= 'cd',
           col="black", size=6, parse=TRUE)  +
  annotate("text", x=0.65, y=345, label= '"(b)"',
           col="black", size=6, parse=TRUE)

Figure1b_species_biomass +scale_color_manual(guide = guide_legend(title = "Grass Type"), values=wes_palette(n=2, name="GrandBudapest1")) ->Fig1.panelB
Fig1.panelB

#import species infection figure:
source("virus_workflow.R")

