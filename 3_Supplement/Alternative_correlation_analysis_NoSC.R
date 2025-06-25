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