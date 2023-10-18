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

EpiContext="COVID-19_Wuhan"

####
#### 1. COVID-19 #####################
####

Ytitle="COVID-19"
xAxisTitle = "\nDate of emergence (2019)"

Ncases_var=c(169,202,3072) # Baseline: 3072
LastDate = as.Date(c("2019-12-31","2019-12-31","2020-01-19"))
DataSets = c("Pekar 2022","WHO 2020","Pekar 2022")
MyLabels = c("Jan 19, 2020 (Pekar et al. 2022)","Dec 31, 2019 (Pekar et al. 2022)","Dec 31, 2019 (WHO 2020)")

x_lims=as.Date(c("2019-10-21", "2019-12-14"))
MyLimits = as.Date(c("2019-11-14","2019-12-28"))

MyBreaks=c(as.Date("2019-07-01"),
           seq(as.Date("2019-07-07"),as.Date("2019-07-28"),by="7 days"),
           as.Date("2019-08-01"),
           seq(as.Date("2019-08-07"),as.Date("2019-08-28"),by="7 days"),
           as.Date("2019-09-01"),
           seq(as.Date("2019-09-07"),as.Date("2019-09-28"),by="7 days"),
           as.Date("2019-10-01"),
           seq(as.Date("2019-10-07"),as.Date("2019-10-28"),by="7 days"),
           as.Date("2019-11-01"),
           seq(as.Date("2019-11-07"),as.Date("2019-11-28"),by="7 days"),
           as.Date("2019-12-01"),
           as.Date("2019-12-07"),
           as.Date("2019-12-10"))

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
    if (Ncases==169){ # (Pekar et al. 2020)
        SA_Ncases=read.csv(paste0("Figs_Appendix/FigS5-S6_SensitivityAnalyses_Datasets/Output/COVID-19_Wuhan/Cases_EpiSize_Time_",Ncases,"cases.csv")) %>%
            as_tibble()
        
        Date_N = as.Date("2019-12-31") 
    } else if (Ncases==202){# (WHO 2020)
        SA_Ncases=read.csv(paste0("Figs_Appendix/FigS5-S6_SensitivityAnalyses_Datasets/Output/COVID-19_Wuhan/Cases_EpiSize_Time_",Ncases,"cases.csv")) %>%
            as_tibble()
        
        Date_N = as.Date("2019-12-31")
    }else if (Ncases==3072){ # (Pekar et al. 2020)
        SA_Ncases=read.csv("Fig1-2_Emergence/Output/COVID-19_Wuhan/Cases_EpiSize_Time_3072cases.csv") %>% 
            as_tibble()
        
        Date_N = as.Date("2020-01-19")
    }

    colnames(SA_Ncases)=c("Time","Cases","EpiSize","Case1")
    
    
    SA_Ncases=SA_Ncases%>%
        mutate(EpiContext=EpiContext,
               DataSet = DataSets[k],
               Ncases=Ncases,
               Date=Date_N-Time)
    
    SA_Ncases = SA_Ncases %>%
        relocate(EpiContext, .before=Time) %>%
        relocate(Ncases, .before=Time )
    
    # Build dataframe of stats
    Results_SA = rbind(Results_SA,
                       tibble(EpiContext=EpiContext,
                              DataSet = DataSets[k],
                              Ncases = Ncases,
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

## Mean time by Ncases
mu_Ncases = ddply(SA_Ncases_long, "Ncases", summarise, grp.mean=mean(Time))
me_Ncases = ddply(SA_Ncases_long, "Ncases", summarise, grp.median=median(Time))

me_Ncases = me_Ncases %>% 
    mutate(Date_N = LastDate,
           Date_1sCase = Date_N-grp.median)

####
#### 2. Plot  Violins  #####################
####
p_Ncases = ggplot(data=SA_Ncases_long, 
                          aes(x=Date, y=Ncases, fill=Ncases,color=Ncases)) +
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
             x=as.Date("2019-09-28"), y=3.3, 
             hjust=0,
             label="(Pekar et al. 2022) data up to Jan 19, 2020",
             size=3,
             color=ColorOne) +
    annotate(geom="text", 
             x=as.Date("2019-09-28"), y=2.3, 
             hjust=0,
             label="(Pekar et al. 2022) data up to Dec 31, 2019",
             size=3,
             color=ColorFour) +
    annotate(geom="text", 
             x=as.Date("2019-09-28"), y=1.3, 
             hjust=0,
             label="(WHO 2020) data up to Dec 31, 2019",
             size=3,
             color=ColorTwo) + 
    ## First cases
    geom_segment(x=as.Date("2019-12-10"), xend=as.Date("2019-12-10"), y=0,yend=3,
                 color=ColorOne,linetype="dashed",size=0.3) +
    geom_segment(x=as.Date("2019-12-10"), xend=as.Date("2019-12-10"), y=0,yend=2,
                 color=ColorFour,linetype="dashed",size=0.3) +
    geom_segment(x=as.Date("2019-12-02"), xend=as.Date("2019-12-02"), y=0,yend=1,
                 color=ColorTwo,linetype="dashed",size=0.3) +
    ## Add stats to legend
    annotate("point", shape=5, x=as.Date("2019-10-10"), y=4.2,size=3) +
    annotate("text",x=as.Date("2019-10-14"), y=4.2,label="Mean",hjust=0,size=2.75) +
    annotate("point", shape=3, x=as.Date("2019-11-05"), y=4.2,size=3) +
    annotate("text",x=as.Date("2019-11-09"), y=4.2,label="Median",hjust=0,size=2.75) +
    ## Colors
    scale_fill_manual(name="", values= c(ColorFour,ColorTwo,ColorOne)) +
    scale_color_manual(name="", values= c(ColorFour,ColorTwo,ColorOne)) +
    ## Axis and Labels
    scale_x_date(name="\nDate of emergence (2019)",
                 limits=as.Date(c("2019-09-24","2019-12-15")),
                 breaks=MyBreaks,
                 date_labels="%b %d",
                 expand=c(0,0)) +
    scale_y_discrete(name="Distributions",
                     label=c("","","",""),
                     expand=c(0,0.7),
                     limits=c("202","169","3072","")) +
    theme_classic() +
    theme(legend.position="none",
          axis.text.x=element_text(angle=90, hjust=1, vjust=1),
          axis.ticks.length=unit(6, "pt"),
          axis.ticks.y = element_blank()) +
    NULL

p_Ncases

ggsave("Figs_Appendix/FigS5-S6_SensitivityAnalyses_Datasets/Output/COVID-19_Wuhan/FigS6_SensitivityAnalyses_Datasets_COVID-19_Wuhan.pdf",
       plot=p_Ncases, height=10, width=12, units=c("cm"))

####
#### 4.Build TeX table #####################
####
Table_Results_SA = Results_SA %>%
    select(-Mean) %>%
    mutate(EarliestDate = format(Earliest, "%b %d, %Y"),
           Results = paste0(format(Median, "%b %d")," (",format(P025, "%b %d"),'--',format(P975, "%b %d"),")", format(P975, " %Y"))) %>%
    select(-c(Median,Earliest,P025,P975))
Table_Results_SA


