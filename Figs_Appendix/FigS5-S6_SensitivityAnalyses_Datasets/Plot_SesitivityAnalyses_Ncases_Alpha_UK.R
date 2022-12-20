#### Estimating the time distribution between the first infection and N cases
####
#### Jijon S, Czuppon P, Blanqart F. and Débarre F
#### iEES, 2022
####
#### 0. SETUP ################################################################
####

## PACKAGES 
library(ggplot2)
library(tidyverse)
library(lubridate)
library(plyr)
library(latex2exp)
library(cowplot)


## Colors
ColorOne = "#00648C" # Blue
ColorTwo = "#A01E18" # Red
ColorFour = "#6C6CEB" # Violet

EpiContext="Alpha_UK"

####
#### 1. Aplha #####################
####

Ytitle="Alpha"
xAxisTitle = "\nDate of emergence (2020)"

Ncases_var=c(1,406)
LastDate = as.Date(c("2020-09-20","2020-11-11"))
MyLabels = c("N=1","N=406")

x_lims=as.Date(c("2020-05-28", "2020-09-30"))

MyBreaks=c(seq(as.Date("2020-04-07"),as.Date("2020-04-28"),by="7 days"),
           as.Date("2020-05-01"),
           seq(as.Date("2020-05-07"),as.Date("2020-05-28"),by="7 days"),
            as.Date("2020-06-01"),
            seq(as.Date("2020-06-07"),as.Date("2020-06-28"),by="7 days"),
            as.Date("2020-07-01"),
            seq(as.Date("2020-07-07"),as.Date("2020-07-28"),by="7 days"),
            # as.Date("2020-07-26"),
            as.Date("2020-08-01"),
            seq(as.Date("2020-08-07"),as.Date("2020-08-28"),by="7 days"),
            # as.Date("2020-08-26"),
            as.Date("2020-09-01"),
            seq(as.Date("2020-09-07"),as.Date("2020-09-28"),by="7 days"))

SA_Ncases_long=NULL
SA_tol_delay_long=NULL
Results_SA = NULL


####
#### 1. Read results #####################
####
k=0
for (Ncases in Ncases_var){
    k=k+1
    # Strings
    Ncases_str=as.character(Ncases)
    
    # Read data
    if (Ncases==1){
        SA_Ncases=read.csv("Figs_Appendix/FigS5-S6_SensitivityAnalyses_Datasets/Output/Alpha_UK/Cases_EpiSize_Time_1cases.csv") %>%
            as_tibble()
        
        Date_N = as.Date("2020-09-20")
    }else if (Ncases==406){
        SA_Ncases=read.csv("Fig1-2_Emergence/Output/Alpha_UK/Cases_EpiSize_Time_406cases.csv") %>% 
            as_tibble()
        
        Date_N = as.Date("2020-11-11")
    }
    
    colnames(SA_Ncases)=c("Time","Cases","EpiSize","Case1")
    
    
    SA_Ncases=SA_Ncases%>%
        mutate(EpiContext=EpiContext,
               Study = paste0("Estimates (N=",Ncases,")"),
               Ncases=Ncases,
               Date=Date_N-Time,
               DateN = Date_N)
    
    SA_Ncases = SA_Ncases %>%
        relocate(EpiContext, .before=Time) %>%
        relocate(Ncases, .before=Time )
    
    # Build dataframe of stats
    Results_SA = rbind(Results_SA,
                       tibble(EpiContext=EpiContext,
                              Ncases = Ncases,
                              Study = paste0("Estimates (N=",Ncases,")"),
                              Mean = mean(SA_Ncases$Date),
                              Median = median(SA_Ncases$Date),
                              P025 = quantile(SA_Ncases$Date,0.025,type=1)[[1]],
                              P975 = quantile(SA_Ncases$Date,0.975,type=1)[[1]],
                              Earliest = min(SA_Ncases$Date)))
    
    # Build dataframe with all values
    SA_Ncases_long=rbind(SA_Ncases_long,SA_Ncases)
}
SA_Ncases_long
SA_Ncases_long$Ncases=as_factor(SA_Ncases_long$Ncases)

SA_Ncases_long = SA_Ncases_long %>%
    select(Date,Time,Ncases,Study)

##
## Czuppon et al. 2021
##
Res_Czuppon2021 = read.csv(file="Data/Emergence_Czuppon2021updated.csv") %>% 
    as_tibble()
Res_Czuppon2021$Date = as.Date(Res_Czuppon2021$Date)
Date_N = as.Date("2020-09-20")
Ncases = 1

Res_Czuppon2021 = Res_Czuppon2021 %>%
    mutate(Time=as.numeric(Date_N-Date),
           Ncases=1,
           Study="Czuppon et al. 2021 (updated, N=1)") %>% 
    select(Date,Time,Ncases,Study)

Results_SA = rbind(Results_SA,
                   tibble(EpiContext=EpiContext,
                          Ncases = Ncases,
                          Study = "Czuppon et al. 2021 (updated, N=1)",
                          Mean = mean(Res_Czuppon2021$Date),
                          Median = median(Res_Czuppon2021$Date),
                          P025 = quantile(Res_Czuppon2021$Date,0.025,type=1)[[1]],
                          P975 = quantile(Res_Czuppon2021$Date,0.975,type=1)[[1]],
                          Earliest = min(Res_Czuppon2021$Date)))

SA_Ncases_long = rbind(SA_Ncases_long,Res_Czuppon2021)
SA_Ncases_long$Study=as_factor(SA_Ncases_long$Study)

SA_Ncases_long


####
#### 2. Plot  Violins  #####################
####
p_Ncases = ggplot(data=SA_Ncases_long, 
                          aes(x=Date, y=Study, fill=Study,color=Study)) +
    geom_violin(width=0.8, size=0.2, alpha=0.3) + 
    stat_summary(fun = "mean",
                 geom = "point",
                 shape=5, # Diamond
                 size=3) +
    stat_summary(fun = "median",
                 geom = "point",
                 shape=3, # Cross
                 size=3) +
    ## Annotate estimates
    annotate(geom="text", 
             x=as.Date("2020-06-04"), y=1.3, 
             hjust=0,
             label="Estimates (N=406)",
             size=3,
             color=ColorOne) +
    annotate(geom="text", 
             x=as.Date("2020-06-04"), y=3.3, 
             hjust=0,
             label="Estimates (N=1)",
             size=3,
             color=ColorFour) +
    annotate(geom="text", 
             x=as.Date("2020-06-04"), y=2.3, 
             hjust=0,
             label="Czuppon et al., 2021 (updated; N=1)",
             size=3,
             color=ColorTwo) + 
    ## Add stats to legend
    annotate("point", shape=5, x=as.Date("2020-06-26"), y=4.2,size=3) +
    annotate("text",x=as.Date("2020-06-30"), y=4.2,label="Mean",hjust=0,size=2.75) +
    annotate("point", shape=3, x=as.Date("2020-08-03"), y=4.2,size=3) +
    annotate("text",x=as.Date("2020-08-07"), y=4.2,label="Median",hjust=0,size=2.75) +
    ## Colors
    scale_fill_manual(name="", values= c(ColorFour,ColorOne,ColorTwo)) +
    scale_color_manual(name="", values= c(ColorFour,ColorOne,ColorTwo)) +
    ## Axis and Labels
    scale_x_date(name="\nDate of emergence (2019)",
                 limits=x_lims,
                 breaks=MyBreaks,
                 date_labels="%b %d",
                 expand=c(0,0)) +
    scale_y_discrete(name="Distributions",
                     label=c("","","",""),
                     limits = c("Estimates (N=406)","Czuppon et al. 2021 (updated, N=1)","Estimates (N=1)"," "),
                     expand=c(0,0.7)) +
    theme_classic() +
    theme(legend.position="none",
          axis.text.x=element_text(angle=90, hjust=1, vjust=1),
          axis.ticks.length=unit(6, "pt"),
          axis.ticks.y = element_blank()) +
    NULL

p_Ncases

ggsave("Figs_Appendix/FigS5-S6_SensitivityAnalyses_Datasets/Output/Alpha_UK/FigS5_SensitivityAnalyses_Datasets_Alpha_UK.pdf",
       plot=p_Ncases, height=10, width=12, units=c("cm"))

####
#### 4.Build Results table #####################
####
Table_Results_SA = Results_SA %>%
    select(-Mean) %>%
    mutate(EarliestDate = format(Earliest, "%b %d, %Y"),
           Results = paste0(format(Median, "%b %d")," (",format(P025, "%b %d"),'--',format(P975, "%b %d"),")", format(P975, " %Y"))) %>%
    select(-c(Median,Ncases,Earliest,P025,P975))
Table_Results_SA

