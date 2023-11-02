#### Jijón S, Czuppon P, Blanquart F & Débarre F (2023). 
#### Using early detection data to estimate the date of emergence of an epidemic outbreak.
#### https://github.com/sjijon/estimate-emergence-from-data
####
#### Plot epicurve from data
####
#### 0. SETUP ################################################################
####

## CLEAR ALL
rm(list = ls())

## PACKAGES 
library(ggplot2)
library(tidyverse)
library(lubridate)
library(ggh4x)


## USER-DETERMINED OPTIONS

## Colors
MyBlue = "#00648C"      # Blue
MyLightBlue = "#BDCED6" # LightBlue

## Save results?
SAVE_RES = "YES";
# SAVE_RES = "NO";

####
####
#### 1. UK Alpha cases #########################
####

##
## Read data collected by November 30, 2020 (main results)
##
Data_UK_SubmNov30 = read.csv(file = "Data/CaseData_Alpha_UK.csv") %>%
    as_tibble()
# Convert to date
Data_UK_SubmNov30$Date = Data_UK_SubmNov30$Date %>% as.Date()

##
## Read data collected by December 15, 2020
##
Data_UK_SubmDec15 = read.csv(file = "Figs_Appendix/FigS1-S2_CaseData/SupplementaryData/CaseData_Alpha_UK_SubmittedByDec15.csv") %>%
    as_tibble()
## Rename columns
colnames(Data_UK_SubmDec15) = c("Date","Cases")
## Convert to date
Data_UK_SubmDec15$Date = Data_UK_SubmDec15$Date %>% as.Date()
## Add cumulative cases
Data_UK_SubmDec15 = Data_UK_SubmDec15 %>%
    mutate(Cumul = cumsum(Cases))

# Data_UK_SubmNov30
# Data_UK_SubmDec15

## Dates
Dates_UK=c(as.Date("2020-09-01"),
           seq(as.Date("2020-09-07"),as.Date("2020-09-28"),by="7 days"),
           as.Date("2020-10-01"),
           seq(as.Date("2020-10-07"),as.Date("2020-10-28"),by="7 days"),
           as.Date("2020-11-01"),
           seq(as.Date("2020-11-07"),as.Date("2020-11-28"),by="7 days"),
           as.Date("2020-12-01"),
           seq(as.Date("2020-12-07"),as.Date("2020-12-28"),by="7 days"))

Dates_UK_minor=seq(as.Date("2020-09-01"),as.Date("2020-12-30"),by="1 days")

coeff_uk = 2.5

p_uk = ggplot(data=Data_UK_SubmNov30, 
              aes(x=Date, y=Cases)) + 
    geom_rect(aes(xmin=as.Date("2020-09-20"), xmax=as.Date("2020-11-11"), 
                  ymin=0, ymax=165),
              fill=MyLightBlue) +
    ## Submitted up to Dec 15
    geom_bar(data=Data_UK_SubmDec15, 
             aes(x=Date, y=Cases),
             stat="identity", 
             fill='white',
             colour='black', 
             alpha=1) +
    ## Submitted up to Nov 30
    geom_bar(stat="identity", 
             fill="black", 
             alpha=1) +
    ## Cumulative cases up to Nov 11
    geom_line(data=Data_UK_SubmNov30[Data_UK_SubmNov30$Date<="2020-11-11",],
              aes(y=Cumul/coeff_uk),
              color=MyBlue,
              size = 1,
              alpha=0.7) +
    geom_segment(aes(x=as.Date("2020-11-11"), xend=as.Date("2020-12-05"), 
                     y = 406/coeff_uk, yend = 406/coeff_uk),
                 color=MyBlue,
                 size = 0.3,
                 # size = 1,
                 # alpha=0.02,
                 linetype="dashed") +
    ## Labels
    scale_x_date(name="Date of sample collection (2020)",
                 limits=as.Date(c("2020-09-19","2020-12-05")),
                 breaks=Dates_UK,
                 minor_breaks=Dates_UK_minor,
                 guide="axis_minor",
                 date_labels="%b %d",
                 expand=c(0,0)) +
    labs(# title = "Alpha variant in UK",
        y= "Daily number of cases") +
    scale_y_continuous(expand = c(0, 3),
                       breaks=seq(0,165,20),
                       sec.axis=sec_axis(~.*coeff_uk,name="Cumulative cases")) +
    ## Theme
    theme_classic() +
    theme(axis.title.y.right=element_text(color=MyBlue,angle = 90),
          axis.text.y.right=element_text(color=MyBlue),
          axis.line.y.right = element_line(color=MyBlue),
          axis.ticks.y.right = element_line(color=MyBlue),
          axis.text.x=element_text(angle=90, hjust=1, vjust=1),
          axis.ticks.length=unit(6, "pt"),
          ggh4x.axis.ticks.length.minor=rel(0.5)) +
    NULL


p_uk

####
#### 2. Wuhan COVID-19 cases #########################
####

## Read csv file
Data_Wuhan = read.csv(file="Data/CaseData_COVID-19_Wuhan.csv") %>%
    as_tibble() 
## Convert to date the Dates column
Data_Wuhan$Date = Data_Wuhan$Date %>% as.Date()
## Up to Jan 31
Data_Wuhan_Jan31 = Data_Wuhan[Data_Wuhan$Date<=as.Date("2020-01-31"),]
## Up to Dec 31 or Jan 
Data_Wuhan = Data_Wuhan[Data_Wuhan$Date<=as.Date("2020-01-19"),]
## Ad cumulative cases
Data_Wuhan = Data_Wuhan %>%
    mutate(Cumul = cumsum(Cases))
Data_Wuhan

N_cases_Wuhan = Data_Wuhan$Cases %>% sum()
N_cases_Wuhan

##
## Dates
##
Dates_WU=c(as.Date("2019-12-01"), as.Date("2019-12-10"),
           seq(as.Date("2019-12-07"),as.Date("2019-12-28"),by="7 days"),
           as.Date("2020-01-01"),
           seq(as.Date("2020-01-07"),as.Date("2020-01-28"),by="7 days"))

Dates_WU_minor=seq(as.Date("2019-12-01"),as.Date("2020-01-28"),by="1 days")

coeff_wu = 2.0

p_wu = ggplot(data=Data_Wuhan, 
              aes(x=Date, y=Cases)) +
    geom_rect(aes(xmin=as.Date("2019-12-10"), xmax=as.Date("2020-01-19"),
                  ymin=0, ymax=1600),
              fill=MyLightBlue) +
    ## Number of cases up to Jan 31
    geom_bar(data=Data_Wuhan_Jan31, 
             aes(x=Date, y=Cases),
             stat="identity", 
             color= "black",
             fill="white", 
             alpha=1) +
    ## Number of cases up to Jan 19
    geom_bar(stat="identity", 
             fill="black", 
             alpha=1) +
    # Cumulative number of cases
    geom_line(data=Data_Wuhan,
              aes(y=Cumul/coeff_wu),
              color=MyBlue,
              size = 1,
              alpha=0.7) +
    geom_segment(aes(x=as.Date("2020-01-19"), xend=as.Date("2020-01-31"),
                     y = 3072/coeff_wu, yend = 3072/coeff_wu),
                 color=MyBlue,
                 size = 0.3,
                 linetype="dashed") +
    scale_x_date(name="\nDate of symptoms onset (2019)",
                 limits=as.Date(c("2019-12-09","2020-01-31")),
                 breaks=Dates_WU,
                 minor_breaks=Dates_WU_minor,
                 guide="axis_minor",
                 date_labels="%b %d",
                 expand = c(0, 0)) +
    labs(
        # title = "COVID-19 in Wuhan",
        y= "Daily number of cases") +
    scale_y_continuous(expand = c(0, 20),
                       limits=c(0,1600),
                       sec.axis=sec_axis(~.*coeff_wu,name="Cumulative cases")) +
    theme_classic() +
    theme(axis.title.y.right=element_text(color=MyBlue,angle = 90),
          axis.text.y.right=element_text(color=MyBlue),
          axis.line.y.right = element_line(color=MyBlue),
          axis.ticks.y.right = element_line(color=MyBlue),
          axis.text.x=element_text(angle=90, hjust=1, vjust=1),
          axis.ticks.length=unit(6, "pt"),
          ggh4x.axis.ticks.length.minor=rel(0.5)) +
    NULL

p_wu


####
#### 3. SAVE FIGURES ############################################################
####
dir.create("Figs_Appendix/FigS1-S2_CaseData/Output",showWarnings = FALSE)

if (SAVE_RES == "YES"){
    ggsave("Figs_Appendix/FigS1-S2_CaseData/Output/FigS1_EpiCurve_Alpha_UK.pdf",
           plot=p_uk, height=12, width=16, units=c("cm"), dpi=600)
    ggsave("Figs_Appendix/FigS1-S2_CaseData/Output/FigS2_EpiCurve_COVID_Wuhan.pdf",
           plot=p_wu, height=12, width=16, units=c("cm"), dpi=600)
}