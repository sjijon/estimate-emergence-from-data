## Estimating the time distribution between the first infection and N cases
##
## Jijon S, Czuppon P, Blanqart F., Débarre F
## iEES, 2022
##
## Plot epicurves and results
##
#################################################
####
#### 0. SETUP ################################################################
####

## CLEAR ALL
rm(list=ls())

## PACKAGES 
library(ggplot2)
library(tidyverse)
library(lubridate)
library(cowplot)
library(ggh4x)
library(viridis) # Color the violin plot

## USER-DETERMINED OPTIONS

## Colors
MyBlue =	 "#00648c"
MyGray =	 "#646464"
MyOrange =	 "#cc5e00"
MyViolet =	 "#843fe5"

## Save results?
SAVE_RES="YES";
# SAVE_RES="NO";

AllEstim = NULL

####
####
#### 1. ABC #########################
####


##
## 2.1 Observed epi curve #######
##

##
## Read data
##
Data_Wuhan = read.csv(file="../Data/CaseData_COVID-19_Wuhan.csv") %>%
    as_tibble() 

## Convert to date the Dates column
Data_Wuhan$Date=Data_Wuhan$Date %>% as.Date()
## Up to Jan 19
Data_Wuhan = Data_Wuhan[Data_Wuhan$Date<=as.Date("2020-01-19"),]
Data_Wuhan

## Number of cases
N_cases_Wuhan=Data_Wuhan$Cases %>% sum()
N_cases_Wuhan

##
## Dates
##
Dates_WU=c(seq(as.Date("2019-09-14"),as.Date("2019-09-28"),by="7 days"),
           as.Date("2019-10-01"),
           seq(as.Date("2019-10-07"),as.Date("2019-10-28"),by="7 days"),
           as.Date("2019-11-01"),
           seq(as.Date("2019-11-07"),as.Date("2019-11-28"),by="7 days"),
           as.Date("2019-12-01"),
           seq(as.Date("2019-12-07"),as.Date("2019-12-28"),by="7 days"),
           as.Date("2020-01-01"),
           seq(as.Date("2020-01-07"),as.Date("2020-01-28"),by="7 days"))

Dates_WU_minor=seq(as.Date("2019-09-10"),as.Date("2020-01-28"),by="1 days")

##
## Plot ###############
##

##
## Our estimates
##
Date_N = as.Date("2020-01-19")


d=read.csv(file="ABCresults.csv") 
dates = seq(as.Date("2019-12-10"),as.Date("2019-09-02"),by ="-1 day")
df = data.frame(date = dates, freq=d[,1]*10000)

aux = matrix(nrow=cumsum(df$freq)[length(cumsum(df$freq))],ncol=2)

j = 1
for (i in 1:length(df$date)) {
  while (j <= cumsum(df$freq[1:i])[i]){
    aux[j,1] = as.character(df$date[i])
    aux[j,2] = "ABC estimate"
    j = j + 1
  }
}

Results_ABC <- as_tibble(data.frame(Date=as.Date(aux[,1]),Study = as.factor(aux[,2])))

## IqR
Results_ABC_IqR = Results_ABC[Results_ABC$Date>=quantile(Results_ABC$Date,0.025,type=1) & Results_ABC$Date<=quantile(Results_ABC$Date,0.975,type=1),] 

## Save all estimates
AllEstim = tibble(EpiContext = "COVID-19_Wuhan",
                  Study = "Jijon et al. (2022)", 
                  Median = median(Results_ABC$Date),
                  P025 = quantile(Results_ABC$Date,0.025,type=1),
                  P975 = quantile(Results_ABC$Date,0.975,type=1),
                  Earliest = sort(Results_ABC$Date)[1])

## Pekar et al 2022
Results_Pekar2022 = read.csv(file="../Data/Emergence_Pekar2022.csv") %>%
    mutate(Study=as.factor("Pekar et al. 2022")) %>%
    as_tibble() 
    
Results_Pekar2022$Date = Results_Pekar2022$Date %>% as.Date()
# Results_Pekar2022

AllEstim = rbind(AllEstim, 
                 tibble(EpiContext = "COVID-19_Wuhan",
                        Study = "Pekar et al. (2022)", 
                        Median = median(Results_Pekar2022$Date),
                        P025 = quantile(Results_Pekar2022$Date,0.025,type=1),
                        P975 = quantile(Results_Pekar2022$Date,0.975,type=1),
                        Earliest = sort(Results_Pekar2022$Date)[1]))

##
## Plot Estimates
##
Estimates_WU = rbind(Results_ABC,Results_Pekar2022)
#levels(Estimates_WU$Study) <- as.factor(c("Pekar et al. (2022)","ABC estimate"))

ggplot(data=Estimates_WU,
                  aes(x=Date,y=Study,fill=Study,color=Study)) +
    geom_violin(width=0.5, size=0.2, alpha=0.3, trim=TRUE) +
    stat_summary(fun = "mean",
                 geom = "point",
                 shape=5,
                 size=3) +
    stat_summary(fun = "median",
                 geom = "point",
                 shape=3,
                 size=3) +
    ## Add IqR
    geom_segment(x=min(Results_ABC_IqR$Date), xend=min(Results_ABC_IqR$Date),
                 y=1.7,yend=2.3,
                 size=0.2, color=MyBlue) + 
    geom_segment(x=max(Results_ABC_IqR$Date), xend=max(Results_ABC_IqR$Date),
                 y=1.7,yend=2.3,
                 size=0.2, color=MyBlue)  +
    ## Annotate estimates
    annotate(geom="text", 
             x=as.Date("2019-12-12"), y=2, 
             hjust=0,
             label="ABC Estimate",
             size=3,
             color=MyBlue) +
    annotate(geom="text", 
             x=as.Date("2019-12-12"), y=1, 
             hjust=0,
             label="Pekar et al. (2022)",
             size=3,
             color=MyViolet) +
    ## Add stats to legend
    annotate("point", shape=5, x=as.Date("2020-01-11"), y=2.2,size=3) +
    annotate("text",x=as.Date("2020-01-13"), y=2.2,label="Mean",hjust = 0,size=2.75) +
    annotate("point", shape=3, x=as.Date("2020-01-11"), y=1.8,size=3) +
    annotate("text",x=as.Date("2020-01-13"), y=1.8,label="Median",hjust = 0,size=2.75) +
    ## Colors
    scale_fill_manual(name="", values= c(MyBlue,MyViolet)) +
    scale_color_manual(name="", values= c(MyBlue,MyViolet))  +
    ## Axis and Labels
    scale_y_discrete(name="Distributions",
                     label=c("   ", ""), 
                     limits=rev(levels(Estimates_WU$Study)),
                     expand=c(0,0.5)) +
    scale_x_date(name="Date of emergence",
               limits=as.Date(c("2019-09-11","2020-01-22")),
               breaks=Dates_WU,
               minor_breaks=Dates_WU_minor,
               guide="axis_minor",
               date_labels="%b %d",
               expand=c(0,0)) +
    theme_classic() +
    theme(legend.position="none",
          axis.text.x=element_text(angle=90, hjust=1, vjust=1),
          axis.ticks.length=unit(6, "pt"),
          ggh4x.axis.ticks.length.minor=rel(0.5),
          plot.margin = unit(c(t=1,r=0.2,b=0,l=0.2), "cm")) +
    guides(fill = guide_legend(byrow = TRUE))  # Important to increase space between legend elements
    

 

####
#### 3. PRINT ALL ESTIMATES ############
####

AllEstim

####
#### 4. SAVE RESULTS ###############
####
#if (SAVE_RES == "YES"){
#    ggsave("FigABC/Output/ABC_Pekar_Comp.pdf",
#           plot=p_wu, height=12.5, width=16, units=c("cm"))
#    
#}