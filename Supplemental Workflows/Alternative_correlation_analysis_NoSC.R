# Is the positive correlation between C3 grass productivity and community-level BYDV prevalence driven by S. cereale?

source("3_Supplement/Supplemental_workflow.R")

view(sp.grass.data)
names(sp.grass.data)

#Calculate mean grass biomass for each plot:
sp.grass.data %>% 
  filter(!Species == "Secale cereale") %>%  #exclude S. cereale
  group_by(Year, Block, Plot, FuncDiversity, planted_fg_richness, GrassType)%>% 
  summarise(tot.plot.adj = sum(Adj_biomass),
            mean.plot.infection = mean(mean.sp.infection, na.rm=TRUE)) -> plot.total.grass.biomass


grass.model.no.sc <- lmer(log(tot.plot.adj) ~  mean.plot.infection * GrassType * FuncDiversity + Year + 
                      (1|Year:Block), data= plot.total.grass.biomass)

summary(grass.model.no.sc)
anova(grass.model.no.sc, type=3) 

#Post-hoc test for mean.plot infection x Grass Type interaction WITHOUT SECALE CEREALE
emtrends(grass.model.no.sc, pairwise~GrassType, var = "mean.plot.infection", infer=TRUE)

# VISUALIZE C3 PRODUCTIVITY - BYDV RELATIONSHIP WITHOUT SECALE CEREALE

ggplot(plot.total.grass.biomass, aes(x=mean.plot.infection, y=log(tot.plot.adj))) +
  geom_point() +
  xlab("BYDV-PAV Prevalence in Grasses") +
  ylab(bquote('Log Grass Productivity' ~ (g/m^-2))) +
  theme_bw() + 
  theme(axis.title = element_text(size=16),
        axis.text  = element_text(size=13)) + 
  facet_grid(~GrassType) -> grass.regression.fig2

grass.regression.fig2


# Add a regression line to only one specific facet panel, using subset() to filter data
grass.regression.fig2 +
  geom_smooth(data = subset(grass.regression.df, GrassType == "C3"),
              color = "#E69F00", method = "lm", se = TRUE) +
  geom_text(data = data.frame(
    mean.plot.infection = c(0.59, 0.59, 0.59, 0.59),
    tot.plot.adj = c(7.30, 7.10, 7.30, 7.10),
    GrassType = c("C3", "C3", "C4", "C4"),
    label = c("~'t=1.958'", "~ 'P=0.05'","~'t= -1.331'", "~ 'P=0.187'")
  ),
  aes(x = mean.plot.infection, y = tot.plot.adj, label = label),
  parse = TRUE) -> supp.fig_no_SC

supp.fig_no_SC

##############################################################
# Is the positive correlation between C3 grass productivity and community-level BYDV prevalence 
#driven by A. fatua?

#Calculate mean grass biomass for each plot:
sp.grass.data %>% 
  filter(!Species == "Avena fatua") %>%  #exclude A. fatua
  group_by(Year, Block, Plot, FuncDiversity, planted_fg_richness, GrassType)%>% 
  summarise(tot.plot.adj = sum(Adj_biomass),
            mean.plot.infection = mean(mean.sp.infection, na.rm=TRUE)) -> plot.total.grass.biomass.no.af


grass.model.no.af <- lmer(log(tot.plot.adj) ~  mean.plot.infection * GrassType * FuncDiversity + Year + 
                            (1|Year:Block), data= plot.total.grass.biomass.no.af)

summary(grass.model.no.af)
anova(grass.model.no.af, type=3) 

#Post-hoc test for mean.plot infection x Grass Type interaction WITHOUT AVENA FATUA
emtrends(grass.model.no.af, pairwise~GrassType, var = "mean.plot.infection", infer=TRUE)

##############################################################
# Is the positive correlation between net community productivity and community-level BYDV prevalence 
#driven by A. fatua?

source("2_Analyses_and_Figures/biomass_workflow.R")

#Calculate total plot biomass: 
biomass.data %>% 
  filter(Species!= c('other', 'Avena fatua')) %>%  #exclude weeds
  group_by(Year, Block, Plot, FuncDiversity, GrassType, planted_fg_richness) %>% 
  summarize(
    total_plot_biomass = sum(Biomass_grams.m2)) %>%
  ungroup() -> plot_biomass.no.af

#Calculate mean plot infection without A. fatua
sp.grass.data %>%
  filter(!Species == "Avena fatua") %>%  #exclude S. cereale
  group_by(Year, Block, Plot, FuncDiversity, planted_fg_richness, GrassType)%>% 
  summarise(mean.plot.infection = mean(mean.sp.infection, na.rm=TRUE)) -> plot.inf.no.af

#join mean plot infection data to mean plot grass biomass data by the "Year" and "Plot" columns:
plot_biomass.no.af %>% 
  select(-c(FuncDiversity, planted_fg_richness, Block)) %>% 
  left_join(plot.inf.no.af %>% select(Year, Plot, mean.plot.infection), by = c("Year", "Plot")) ->comm.regression.df.no.af


#Check data frame structure: make planted fg richness factor
comm.regression.df.no.af$planted_fg_richness <- as.factor(comm.regression.df.no.af$planted_fg_richness)

################################

# Look at the distribution of the community biomass data:
hist(comm.regression.df.no.af$total_plot_biomass)
# The data does not look normally distributed --> Log transform

hist(log(comm.regression.df.no.af$total_plot_biomass))
# Log transformed data is normally distributed

#Community Productivity (WITHOUT AVENA FATUA)
community.model.no.af <- lmer(log(total_plot_biomass) ~  mean.plot.infection * GrassType * FuncDiversity + Year + 
                          (1|Year:Block), data= comm.regression.df.no.af)
summary(community.model.no.af)
anova(community.model.no.af, type=3)


