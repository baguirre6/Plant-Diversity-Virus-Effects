#Supplemental Analyses of Climate Factors
#Beatriz A. Aguirre
#February 22, 2024

#Clear environment
rm(list= ls())

#Load required libraries
library(dplyr)
library(ggplot2)
library(emmeans)
library(ggsignif)
library(ggpubr)
library(googledrive)

#################################################################
# Import data
#################################################################

#drive authentication
drive_auth()

#create data folder
dir.create(file.path("data"), showWarnings = F)

#identify desired file
focal_file <- "climate_data.csv"

#download data file
googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/1dpjkQh9MKMO8LInw4A1IaRXhddenbFzq")) %>% 
  dplyr::filter(name == focal_file) %>% 
  googledrive::drive_download(file = .$id, overwrite = T,
                              path = file.path("data", .$name))

# Read in infection data
clim.data <- read.csv(file = file.path("data", "climate_data.csv"))

# Check data structure
clim.data$year <- as.factor(clim.data$year)

#################################################################
#Use day-time (6AM -6PM) data for growing/experiment season (May-August) by year
#Mean, max, min day time temperature
#Mean, max, min day time RH
#Mean, max, min day time VPD
#################################################################

#Calculate daytime daily mean, max and min temp, RH and VPD for each date
clim.data %>% 
  filter(hour %in% c('6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18'),
         month %in% c('may', 'jun', 'jul', 'aug')) %>% 
  group_by(date, year, month, date_num) %>% 
  summarize(mean_day_daily_temp = mean(Temp_c), 
            mean_day_daily_RH = mean(RH),
            mean_day_daily_vpd = mean(VPD_kpa),
            max_day_daily_temp = max(Temp_c), 
            max_day_daily_RH = max(RH),
            max_day_daily_vpd = max(VPD_kpa),
            min_day_daily_temp = min(Temp_c), 
            min_day_daily_RH = min(RH),
            min_day_daily_vpd = min(VPD_kpa)
  ) -> daytime.daily.means

#################################################################
# Run Statistics on mean daily environmental factors - Table S1
#################################################################

#Mean Temp
m1 <- lm(mean_day_daily_temp ~ year, data=daytime.daily.means)
summary(m1)
anova(m1)
#No difference in daytime mean daily temperatures between years

#Max temp
m1.max <- lm(max_day_daily_temp ~ year, data=daytime.daily.means)
summary(m1.max)
emmeans(m1.max, ~ year, type="response")
#Marginal difference in daytime max daily temperatures between years

#Min temp
m1.min <- lm(min_day_daily_temp ~ year, data=daytime.daily.means)
summary(m1.min)
#No difference in daytime min daily temperatures between years
emmeans(m1.min, ~year, type="response")


######################################################
# Relative Humidity - Table S1
######################################################
#Mean RH
m2 <- lm(mean_day_daily_RH ~ year, data=daytime.daily.means)
summary(m2)
#Highly significant difference in daytime mean relative humidity between years****
emmeans(m2, ~ year, type="response")

#Max RH
m2.max <- lm(max_day_daily_RH ~ year, data=daytime.daily.means)
summary(m2.max)
#No significant difference in daytime max relative humidity between years
emmeans(m2.max, ~ year, type="response")

#Min RH
m2.min <- lm(min_day_daily_RH ~ year, data=daytime.daily.means)
summary(m2.min)
#Significant difference in daytime min relative humidity between years ****
emmeans(m2.min, ~ year, type="response")

######################################################
# VPD - Table S1
######################################################
#Mean VPD
m3 <- lm(mean_day_daily_vpd ~ year, data=daytime.daily.means)
summary(m3)
#Highly significant difference in daytime mean VPD between years ****

#Max VPD
m3.max <- lm(max_day_daily_vpd ~ year, data=daytime.daily.means)
summary(m3.max)
#Highly significant difference in daytime max VPD between years ****
emmeans(m3.max, ~ year, type="response")

#Min VPD
m3.min <- lm(min_day_daily_vpd ~ year, data=daytime.daily.means)
summary(m3.min)
#No significant difference in daytime min VPD between years
emmeans(m3.min, ~ year, type="response")

#################################################################
#Visualize data
#################################################################

##visualize daytime mean daily RH data

#Export emmeans data into a data frame for plotting w/ggplot
m1.result =as.data.frame(summary(emmeans(m1, ~year, type="response")))
m1.result

panel.1<-ggplot(m1.result,
                aes(x=year, y=emmean)) +
  geom_point(aes(color=year, shape=year), size=5) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=year), width=.05, position=position_dodge(0.8)) +
  ylab("Mean daytime\n daily temperature (C)") +
  xlab("Year") +
  ylim(c(20,23.5)) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold", size=15),
        axis.text    = element_text(face = "bold", size=13),
        plot.caption = element_text(hjust = 0, size=10),
        legend.position = "none") +
  geom_signif(comparisons = list(c("2021", "2022")), annotations = "", map_signif_level=TRUE,
              textsize=10, y_position = 22.5, tip_length = 0.05, vjust=0.4) +
  annotate("text", x=1.5, y=22.75, label= 'italic("NS")',
           col="black", size=4, parse=TRUE) +
  annotate("text", x=0.60, y=23.5, label= '"(b)"',
           col="black", size=6, parse=TRUE) 

panel.1

#Export emmeans data into a data frame for plotting w/ggplot
m2.result =as.data.frame(summary(emmeans(m2, ~year, type="response")))
m2.result

panel.2<-ggplot(m2.result,
                aes(x=year, y=emmean)) +
  geom_point(aes(color=year, shape=year), size=5) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=year), width=.05, position=position_dodge(0.8)) +
  ylab("Mean daytime\n daily RH (%)") +
  xlab("Year") +
  ylim(c(61,73)) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold", size=15),
        axis.text    = element_text(face = "bold", size=13),
        plot.caption = element_text(hjust = 0, size=10),
        legend.position = "none") +
  geom_signif(comparisons = list(c("2021", "2022")), annotations = "*", map_signif_level=TRUE,
              textsize=10, y_position = 71, tip_length = 0.05, vjust=0.4) +
  annotate("text", x=1.5, y=72.75, label= 'italic("P<0.01")',
           col="black", size=4, parse=TRUE) +
  annotate("text", x=0.60, y=73, label= '"(c)"',
           col="black", size=6, parse=TRUE) 

panel.2

##visualize daytime mean daily VPD data

#Export emmeans data into a data frame for plotting w/ggplot
m3.result =as.data.frame(summary(emmeans(m3, ~year, type="response")))
m3.result

panel.3<-ggplot(m3.result,
                aes(x=year, y=emmean)) +
  geom_point(aes(color=year, shape=year), size=5) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=year), width=.05, position=position_dodge(0.8)) +
  ylab("Mean daytime\n daily VPD (kPa)") +
  xlab("Year") +
  ylim(c(.7, 1.3)) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold", size=15),
        axis.text    = element_text(face = "bold", size=13),
        plot.caption = element_text(hjust = 0, size=10), 
        legend.position="none") +
  geom_signif(comparisons = list(c("2021", "2022")), annotations = "*", map_signif_level=TRUE,
              textsize=10, y_position = 1.2, tip_length = 0.05, vjust=0.4) +
  annotate("text", x=1.5, y=1.275, label= 'italic("P<0.01")',
           col="black", size=4, parse=TRUE) +
  annotate("text", x=0.60, y=1.30, label= '"(d)"',
           col="black", size=6, parse=TRUE) 

panel.3
#ggsave("~/Desktop/Dissertation/1_SBF/Ithaca Airport Climate Data/Figures/Mean_VPD_figure.pdf", width =3, height = 5)

##################################################################################
#Calculate daily total precipitation by year or growing season precip

##################################################################################

#identify desired file
focal_file <- "Precip_data.csv"

#download data file
googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/u/0/folders/1dpjkQh9MKMO8LInw4A1IaRXhddenbFzq")) %>% 
  dplyr::filter(name == focal_file) %>% 
  googledrive::drive_download(file = .$id, overwrite = T,
                              path = file.path("data", .$name))

# Read in infection data
precip.data <- read.csv(file = file.path("data", "Precip_data.csv"))

# check data structure
precip.data$Year <- as.factor(precip.data$Year)

#################################################################
# Visualize precip data 
#################################################################

# reorder month levels
month_order <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')

# Convert Month column to a factor with ordered levels
precip.data$Month <- factor(precip.data$Month, levels = month_order)

precip.data %>% 
  mutate(monthly_precip_mm = MonthlySum_inches * 25.4,
         Norm_monthly_precip_mm = Normal_1991to2020 * 25.4, 
         DepartureFromNormal_mm = DepartureFromNormal * 25.4) %>% 
  select(Year, Month, monthly_precip_mm, Norm_monthly_precip_mm, DepartureFromNormal_mm)->precip.data

#################################################################
# Precipitation Statistics
#################################################################

precip.data %>% 
  filter(Month %in% c('May', 'Jun', 'Jul', 'Aug')) ->gs.data

model.a <- lm(monthly_precip_mm ~ Year, data = gs.data) # Table S1 - monthly growing season precip
summary(model.a)
anova(model.a)
#Marginally significant difference in **growing season** monthly precipitation between 2021 and 2022
emmeans(model.a, ~Year, type="response")

model.a0 <- lm(DepartureFromNormal_mm ~ Year, data = gs.data) # Table S1 - growing season precip deviation from 30 yr avg
summary(model.a0)
anova(model.a0)
#Significant difference (P=0.03) in departure from normal for growing season

# get averages for precip departure: 
methods.text =as.data.frame(summary(emmeans(model.a0, ~Year, type="response")))
methods.text

result = as.data.frame(summary(methods.text)) #Export emmeans data into a data frame for plotting w/ggplot

calcSE<-function(x){
  x2<-na.omit(x)
  sd(x2)/sqrt(length(x2))
}

gs.data %>% 
  group_by(Year) %>% 
  summarise(gs.precip.yr = mean(monthly_precip_mm),
            gs.precip.yr.se = calcSE(monthly_precip_mm)) ->methods.text.final
# 133.5 mm in 2021 growing season - 94 mm (30 year mean) = + 39.5 mm difference
# 39.5 mm is 42.02%  of 94 mm  (42% increase from the 30 year average in 2021)


# 75.2 mm in 2021 growing season - 94 mm (30 year mean) = - 18.8 mm difference
# 18.8 mm is % of 94 mm ( 20% reduction from the 30 year average in 2022)

#################################################################
# Figure S1 - Panel A 
#################################################################

delta.normal.figure <-
  gs.data %>%
  group_by(Year) %>%
  summarise(
    mean.pcp = mean(DepartureFromNormal_mm),
    se.pcp = calcSE(DepartureFromNormal_mm)) %>%
  ggplot(aes(x=as.factor(Year), y=mean.pcp)) +
  geom_point(aes(color= as.factor(Year), shape= as.factor(Year)), size=5) +
  geom_errorbar(aes(ymin= mean.pcp- se.pcp, ymax= mean.pcp + se.pcp, color=as.factor(Year)), width=.05, position=position_dodge(0.8)) +
  ylab("Precipitation departure\n from 30 year average (mm)") +
  xlab("Year") +
  theme_bw() +
  ylim(-60, 85) +
  theme(axis.title   = element_text(face = "bold", size=15),
        axis.text    = element_text(face = "bold", size=13),
        plot.caption = element_text(hjust = 0, size=10),
        legend.position="none") +
  geom_signif(comparisons = list(c("2021", "2022")), annotations = "", map_signif_level=TRUE,
              textsize=10, y_position = 65, tip_length = 0.05, vjust=0.4) +
  annotate("text", x=1.5, y=75, label= '"*"',
           col="black", size=7, parse=TRUE) +
  annotate("text", x=1.5, y=84, label= 'italic("P=0.03")',
           col="black", size=4, parse=TRUE) +
  annotate("text", x=0.60, y=85, label= '"(a)"',
           col="black", size=6, parse=TRUE)

delta.normal.figure

delta.normal.figure +
  geom_hline(yintercept = as.numeric("0"),
             linetype = "dashed", color = "black") -> final.delta.normal.figure

final.delta.normal.figure

####################### #

manuscript.figure.S1 <- ggarrange(final.delta.normal.figure, panel.1, panel.2, panel.3,
                                  ncol = 2, nrow = 2)
manuscript.figure.S1
#ggsave("Figure PDFs/FigureS1.pdf", manuscript.figure.S1, width = 12, height = 12)
