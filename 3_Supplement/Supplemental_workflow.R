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

