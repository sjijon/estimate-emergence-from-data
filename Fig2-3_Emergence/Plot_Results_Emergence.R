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
MyGreen=     "#1d771d"
MyOrange =	 "#cc5e00"
MyViolet =	 "#843fe5"

## Save results?
SAVE_RES="YES";
# SAVE_RES="NO";

AllEstim = NULL

####
####
#### 1. UK Alpha cases #########################
####

##
## 1.1 Observed epi curve #######
##
Date_N = as.Date("2020-11-11")

## Read data
Data_UK=read.csv(file="Data/CaseData_Alpha_UK.csv",head=TRUE) %>%
    as_tibble()

## Convert to date the Dates column
Data_UK$Date=Data_UK$Date %>% as.Date()

## Manipulate dataset
Data_UK=Data_UK %>%
    filter(Date <= Date_N) # Use data of samples collected up to Nov 11
# Data_UK

## Dates for ticks
Dates_UK=c(seq(as.Date("2020-05-07"),as.Date("2020-05-28"),by="7 days"),
           as.Date("2020-06-01"),
           seq(as.Date("2020-06-07"),as.Date("2020-06-28"),by="7 days"),
           as.Date("2020-07-01"),
           seq(as.Date("2020-07-07"),as.Date("2020-07-28"),by="7 days"),
           as.Date("2020-08-01"),
           seq(as.Date("2020-08-07"),as.Date("2020-08-28"),by="7 days"),
           as.Date("2020-09-01"),
           seq(as.Date("2020-09-07"),as.Date("2020-09-28"),by="7 days"),
           as.Date("2020-10-01"),
           seq(as.Date("2020-10-07"),as.Date("2020-10-28"),by="7 days"),
           as.Date("2020-11-01"),
           as.Date("2020-11-07"),
           as.Date("2020-11-11"))
Dates_UK_minor=seq(as.Date("2020-05-01"),as.Date("2020-11-20"),by="1 days")

## Plot epicurve
epi_uk = ggplot(data=Data_UK, 
                aes(x=Date, y=Cases)) +
    ## 1. Epicurve
    geom_bar(stat="identity", 
             width=1, 
             fill="black",
             color="white",
             size=0.3,
             alpha=1,
             position=position_nudge(x=0)) +
    ## 2.First case
    annotate("segment", x=as.Date("2020-09-20"), y=25, xend=as.Date("2020-09-20"), yend=2,
             arrow = arrow(length=unit(6,"pt"))) +
    annotate("text", x=as.Date("2020-09-20"), y=35, label="1st case", size=3) +
    scale_x_date(name="Date of sample collection (2020)",
                 limits=as.Date(c("2020-05-26","2020-11-14")),
                 breaks=Dates_UK,
                 minor_breaks=Dates_UK_minor,
                 guide="axis_minor",
                 date_labels="%b %d",
                 expand=c(0,0)) +
    labs(y="Daily number of Alpha cases") +
    scale_y_continuous(expand=c(0, 2),
    ) +
    theme_classic() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1),
        axis.ticks.length=unit(6, "pt"),
        ggh4x.axis.ticks.length.minor=rel(0.5),
        plot.margin = unit(c(t=1,r=0.2,b=0,l=0.2), "cm")) +
    NULL
# epi_uk

##
## 1.2 Estimates ###############
##

##
## Our estimates
##
Results_Alpha=read.csv(file="Fig1-2_Emergence/Output/Alpha_UK/Cases_EpiSize_Time_406cases.csv") %>% 
    as_tibble()

Results_Alpha = Results_Alpha %>%
    mutate(Date=as.Date(Date_N-MinTime),
           Study=as.factor("Estimates")) %>% 
    select(Date,Study)

## IqR
Results_UK_IqR = Results_Alpha[Results_Alpha$Date>=quantile(Results_Alpha$Date,0.025,type=1) & Results_Alpha$Date<=quantile(Results_Alpha$Date,0.975,type=1),] 


AllEstim = rbind(AllEstim, 
                 tibble(EpiContext = "Alpha_UK",
                        Study = "Jijon et al. (2022)", 
                        Median = median(Results_Alpha$Date),
                        P025 = quantile(Results_Alpha$Date,0.025,type=1),
                        P975 = quantile(Results_Alpha$Date,0.975,type=1),
                        Earliest = sort(Results_Alpha$Date)[1]))

##
## ABC results
##
dates_ABC = seq(as.Date("2019-12-10"),as.Date("2019-09-01"),by ="-1 day") # Set dates

Results_ABC = read.csv(file="ABC_COVID-19_Wuhan/ABCresults.csv",,head=FALSE) %>% # Read results
    as_tibble() %>%
    mutate(Date = dates_ABC,
           Freq = V1 * 10000,
           Study=as.factor("ABC estimates"))  %>%
    select(-V1)
# Results_ABC

# Expand Results_ABC frequencies
aux = matrix(nrow=cumsum(Results_ABC$Freq)[length(cumsum(Results_ABC$Freq))],ncol=2)

j = 1
for (i in 1:length(Results_ABC$Date)) {
    while (j <= cumsum(Results_ABC$Freq[1:i])[i]){
        aux[j,1] = as.character(Results_ABC$Date[i])
        aux[j,2] = "ABC estimates"
        j = j + 1
    }
}
Res_ABC <- as_tibble(data.frame(Date=as.Date(aux[,1]),Study = as.factor(aux[,2])))

## IqR
IqR_ABC = Res_ABC[Res_ABC$Date>=quantile(Res_ABC$Date,0.025,type=1) & Res_ABC$Date<=quantile(Res_ABC$Date,0.975,type=1),] 


AllEstim = rbind(AllEstim, 
                 tibble(EpiContext = "COVID-19_Wuhan",
                        Study = "ABC estimates", 
                        Median = median(Res_ABC$Date),
                        P025 = quantile(Res_ABC$Date,0.025,type=1),
                        P975 = quantile(Res_ABC$Date,0.975,type=1),
                        Earliest = sort(Res_ABC$Date)[1]))
# AllEstim

##
## Czuppon et al. 2021 with negBinom and R=1.9
##
Results_Czuppon2021new = read.csv(file="Data/Emergence_Czuppon2021updated.csv") %>% 
    as_tibble()

Results_Czuppon2021new$Date = Results_Czuppon2021new$Date %>% as.Date()

Results_Czuppon2021new = Results_Czuppon2021new %>%
    mutate(Study="Czuppon et al. 2021 (update)") 

AllEstim = rbind(AllEstim, 
                 tibble(EpiContext = "Alpha_UK",
                        Study = "Czuppon et al. (2021, updated)", 
                        Median = median(Results_Czuppon2021new$Date),
                        P025 = quantile(Results_Czuppon2021new$Date,0.025,type=1),
                        P975 = quantile(Results_Czuppon2021new$Date,0.975,type=1),
                        Earliest = sort(Results_Czuppon2021new$Date)[1]))

##
## Hill et al 2022
##
Results_Hill2022 = read.csv(file="Data/Emergence_Hill2022.csv",sep="\t") %>% 
    as_tibble()

Results_Hill2022$Date = Results_Hill2022$Date %>% as.Date()

Results_Hill2022 = Results_Hill2022 %>%
    mutate(Study="Hill et al. 2022")
Results_Hill2022

AllEstim = rbind(AllEstim, 
                 tibble(EpiContext = "Alpha_UK",
                        Study = "Hill et al. (2022)", 
                        Median = median(Results_Hill2022$Date),
                        P025 = quantile(Results_Hill2022$Date,0.025,type=1),
                        P975 = quantile(Results_Hill2022$Date,0.975,type=1),
                        Earliest = sort(Results_Hill2022$Date)[1]))

## Build data frame
Estimates_UK = rbind(Results_Alpha,Results_Czuppon2021new,Results_Hill2022)

## Plot ##########
estim_uk = ggplot(data=Estimates_UK, 
                  aes(x=Date, y=Study, fill=Study,color=Study)) +
    geom_violin(width=1, size=0.2, alpha=0.3,trim=TRUE) + 
    # geom_boxplot(width=0.25,size=0.2,alpha=0,outlier.shape=NA) +
    stat_summary(fun = "mean",
                 geom = "point",
                 shape=5, # Diamond
                 size=3) +
    stat_summary(fun = "median",
                 geom = "point",
                 shape=3, # Cross
                 size=3) +
    ## Add IqR
    geom_segment(x=min(Results_UK_IqR$Date), xend=min(Results_UK_IqR$Date),
                 y=2.7,yend=3.3,
                 size=0.2, color=MyBlue) + 
    geom_segment(x=max(Results_UK_IqR$Date), xend=max(Results_UK_IqR$Date),
                 y=2.7,yend=3.3,
                 size=0.2, color=MyBlue) +
    ## Annotate estimates
    annotate(geom="text", 
             x=as.Date("2020-09-23"), y=3, 
             hjust=0,
             label="Estimates",
             size=3,
             color=MyBlue) +
    annotate(geom="text", 
             x=as.Date("2020-09-23"), y=2, 
             hjust=0,
             label="Czuppon et al. (2021, updated)",
             size=3,
             color=MyOrange) +
    annotate(geom="text", 
             x=as.Date("2020-09-23"), y=1, 
             hjust=0,
             label="Hill et al. (2022)",
             size=3,
             color=MyGray) + 
    ## Add stats to legend
    annotate("point", shape=5, x=as.Date("2020-10-25"), y=3.2,size=3) +
    annotate("text",x=as.Date("2020-10-28"), y=3.2,label="Mean",hjust=0,size=2.75) +
    annotate("point", shape=3, x=as.Date("2020-10-25"), y=2.8,size=3) +
    annotate("text",x=as.Date("2020-10-28"), y=2.8,label="Median",hjust=0,size=2.75) +
    ## Colors
    scale_fill_manual(name="", values= c(MyBlue,MyOrange,MyGray)) +
    scale_color_manual(name="", values= c(MyBlue,MyOrange,MyGray)) +
    ## Axis and Labels
    scale_x_date(name="\nDate of emergence/tMRCA",
                 limits=as.Date(c("2020-05-26","2020-11-14")),
                 breaks=Dates_UK,
                 minor_breaks=Dates_UK_minor,
                 guide="axis_minor",
                 date_labels="%b %d",
                 expand=c(0,0)) +
    scale_y_discrete(name="Distributions",
                     label=c("    ", "", ""), 
                     expand=c(0,0.7),
                     limits=c("Hill et al. 2022","Czuppon et al. 2021 (update)","Estimates")) +
    theme_classic() +
    theme(legend.position="none",
          axis.text.x=element_text(angle=90, hjust=1, vjust=1),
          axis.ticks.length=unit(6, "pt"),
          ggh4x.axis.ticks.length.minor=rel(0.5)) +
    guides(fill = guide_legend(byrow = TRUE)) + # Important to increase space between legend elements
    NULL

p_uk=plot_grid(epi_uk, estim_uk,
               ncol=1, nrow=2,
               rel_heights=c(1,1))
p_uk

####
#### 2. Wuhan COVID-19 cases #############################################
####

##
## 2.1 Observed epi curve #######
##

##
## Read data
##
Data_Wuhan = read.csv(file="Data/CaseData_COVID-19_Wuhan.csv") %>%
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
## Plot Epicurve
##
epi_wu = ggplot(data=Data_Wuhan, 
                aes(x=Date, y=Cases)) +
    ## 1.Epicurve
    geom_bar(stat="identity",
             width=1,
             fill="black",
             color="white",
             size=0.3,
             alpha=1,
             position=position_nudge(x=0)) +
    ## 3.First case
    annotate("segment", x=as.Date("2019-12-10"), xend=as.Date("2019-12-10"), 
             y=75, yend=15,
             arrow = arrow(length=unit(6,"pt"))) +
    annotate("text", x=as.Date("2019-12-10"), y=100, label="1st case", size=3) +
    ## 4. Axis, legends et al.
    scale_x_date(name="Date of symptoms onset (2019-2020)",
                 limits=as.Date(c("2019-09-11","2020-01-22")),
                 breaks=Dates_WU,
                 minor_breaks=Dates_WU_minor,
                 guide="axis_minor",
                 date_labels="%b %d",
                 expand=c(0, 0)) +
    labs(y="Daily number of COVID-19 cases") +
    scale_y_continuous(expand=c(0, 5),
                       # sec.axis=sec_axis(~.*coeff,name="Cumulative cases")
    ) +
    theme_classic() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1),
        axis.ticks.length=unit(6, "pt"),
        ggh4x.axis.ticks.length.minor=rel(0.5),
        plot.margin = unit(c(t=1,r=0.2,b=0,l=0.2), "cm")) +
    NULL
# epi_wu

##
## 2.2 Estimates ###############
##

##
## Our estimates
##
Date_N = as.Date("2020-01-19")

Results_COVID=read.csv(file="Fig1-2_Emergence/Output/COVID-19_Wuhan/Cases_EpiSize_Time_3072cases.csv") %>% 
    as_tibble() %>%
    mutate(Date=as.Date(Date_N-MinTime),
           Study=as.factor("Estimates")) %>%
    select(Date,Study)

## IqR
Results_COVID_IqR = Results_COVID[Results_COVID$Date>=quantile(Results_COVID$Date,0.025,type=1) & Results_COVID$Date<=quantile(Results_COVID$Date,0.975,type=1),] 

## Save all estimates
AllEstim = tibble(EpiContext = "COVID-19_Wuhan",
                  Study = "Jijon et al. (2022)", 
                  Median = median(Results_COVID$Date),
                  P025 = quantile(Results_COVID$Date,0.025,type=1),
                  P975 = quantile(Results_COVID$Date,0.975,type=1),
                  Earliest = sort(Results_COVID$Date)[1])

## Pekar et al 2022
Results_Pekar2022 = read.csv(file="Data/Emergence_Pekar2022.csv") %>%
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
Estimates_WU = rbind(Results_COVID,Results_Pekar2022)

estim_wu = ggplot(data=Estimates_WU,
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
    geom_segment(x=min(Results_COVID_IqR$Date), xend=min(Results_COVID_IqR$Date),
                 y=1.7,yend=2.3,
                 size=0.2, color=MyBlue) + 
    geom_segment(x=max(Results_COVID_IqR$Date), xend=max(Results_COVID_IqR$Date),
                 y=1.7,yend=2.3,
                 size=0.2, color=MyBlue) +
    geom_segment(x=min(IqR_ABC$Date), xend=min(IqR_ABC$Date),
                 y=1.8,yend=2.2,
                 size=0.2, color=MyGreen) + 
    geom_segment(x=max(IqR_ABC$Date), xend=max(IqR_ABC$Date),
                 y=1.8,yend=2.2,
                 size=0.2, color=MyGreen) + 
    geom_segment(x=min(IqR_Pekar2022$Date), xend=min(IqR_Pekar2022$Date),
                 y=0.8,yend=1.2,
                 size=0.2, color=MyViolet) + 
    geom_segment(x=max(IqR_Pekar2022$Date), xend=max(IqR_Pekar2022$Date),
                 y=0.8,yend=1.2,
                 size=0.2, color=MyViolet) +
    ## Annotate estimates
    annotate(geom="text", 
             x=as.Date("2019-12-12"), y=3, 
             hjust=0,
             label="Estimates",
             size=3,
             color=MyBlue) +
    annotate(geom="text", 
             x=as.Date("2019-12-12"), y=2, 
             hjust=0,
             label="ABC estimates",
             size=3,
             color=MyGreen) +
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
    scale_fill_manual(name="", values= c(MyBlue,MyGreen,MyViolet)) +
    scale_color_manual(name="", values= c(MyBlue,MyGreen,MyViolet))  +
    ## Axis and Labels
    scale_x_date(name="Date of emergence",
                 limits=as.Date(c("2019-09-11","2020-01-22")),
                 breaks=Dates_WU,
                 minor_breaks=Dates_WU_minor,
                 guide="axis_minor",
                 date_labels="%b %d",
                 expand=c(0,0)) +
    scale_y_discrete(name="Distributions",
                     label=c("   ", "", ""),
                     limits=c("Pekar et al. 2022","ABC estimates","Estimates"),
                     expand=c(0,0.5)) +
    theme_classic() +
    theme(legend.position="none",
          axis.text.x=element_text(angle=90, hjust=1, vjust=1),
          axis.ticks.length=unit(6, "pt"),
          ggh4x.axis.ticks.length.minor=rel(0.5),
          plot.margin = unit(c(t=1,r=0.2,b=0,l=0.2), "cm")) +
    guides(fill = guide_legend(byrow = TRUE)) + # Important to increase space between legend elements
    NULL

# Create one figure
p_wu=plot_grid(epi_wu, estim_wu,
               ncol=1, nrow=2,
               rel_heights=c(1,1))
p_wu 
 

####
#### 3. PRINT ALL ESTIMATES ############
####

AllEstim

####
#### 4. SAVE RESULTS ###############
####
if (SAVE_RES == "YES"){
    ggsave("Fig1-2_Emergence/Output/Alpha_UK/Fig1_Emergence_Alpha_UK.pdf",
           plot=p_uk, height=15, width=16, units=c("cm"))
    
    ggsave("Fig1-2_Emergence/Output/COVID-19_Wuhan/Fig2_Emergence_COVID-19_Wuhan.pdf",
           plot=p_wu, height=15, width=16, units=c("cm"))
    
}