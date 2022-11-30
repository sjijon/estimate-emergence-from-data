####
#### Early epidemics
#### iEES
####
#### November, 2022
####
#### Plot results from the sensitivity analyses
#### Varying tol_epi and tol_delay
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
ColorOne = "#FFAA00" # Yellow
ColorTwo = "#00648C" # Blue
ColorThree = "#A01E18" # Red


EpiContext="COVID-19_Wuhan"
# EpiContext="Alpha_UK"

Cond_Name="Cond_Cumul_Delay"

if (EpiContext == "Alpha_UK"){
    Date_N = as.Date("2020-11-11")
    Ytitle="Alpha"
    xAxisTitle = "\nDate of emergence (2020)"
    MyPanelLabels=c("A", "C")
    tol_epi_var=c(0.3,0.5) # Baseline: 0.3
    tol_delay_var=c(0.8,0.9,1.0) # Baseline: 0.9
    
    x_lims=as.Date(c("2020-06-21", "2020-09-07"))
    MyLimits = as.Date(c("2020-04-14","2020-07-28"))
    
    MyBreaks= c(seq(as.Date("2020-04-07"),as.Date("2020-04-28"),by="7 days"),
                as.Date("2020-05-01"),
                seq(as.Date("2020-05-07"),as.Date("2020-05-28"),by="7 days"),
                as.Date("2020-06-01"),
                seq(as.Date("2020-06-07"),as.Date("2020-06-28"),by="7 days"),
                as.Date("2020-07-01"),
                seq(as.Date("2020-07-07"),as.Date("2020-07-28"),by="7 days"),
                as.Date("2020-08-01"),
                seq(as.Date("2020-08-07"),as.Date("2020-08-28"),by="7 days"),
                as.Date("2020-09-01"),
                seq(as.Date("2020-09-07"),as.Date("2020-09-28"),by="7 days"))
} else if (EpiContext == "COVID-19_Wuhan"){
    Date_N = as.Date("2020-01-19")
    
    Ytitle="COVID-19"
    xAxisTitle = "\nDate of emergence (2019)"
    MyPanelLabels=c("B", "D")
    tol_epi_var=c(0.1,0.3,0.5) # Baseline: 0.3
    tol_delay_var=c(0.8,0.9,1.0) # Baseline: 0.9
    
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
}
    
SA_tol_epi_long=NULL
SA_tol_delay_long=NULL
Results_SA = NULL

####
#### 1. Varying tol_epi ################################################################
####
tol_delay=0.9
for (tol_epi in tol_epi_var){
    # Strings
    tol_delay_str=as.character(10*tol_delay)
    tol_epi_str=as.character(10*tol_epi)

    # Read data
    SA_tol_epi=read_csv(paste0("Fig_3_SensitivityAnalyses_Calibration/Output/",EpiContext,"/Varying_tol_epi/Cases_EpiSize_Time_tol_epi_0",tol_epi_str,"_tol_delay_0",tol_delay_str,".csv"),
                        col_names=FALSE,
                        col_types = "ddddd") 
    colnames(SA_tol_epi)=c("Time","Cases","EpiSize","Case1")
    
    SA_tol_epi=SA_tol_epi%>%
        mutate(EpiContext=EpiContext,
               Date=Date_N-Time,
               tol_epi=tol_epi,
               tol_delay=tol_delay) %>%
        as_tibble()
    
    SA_tol_epi = SA_tol_epi %>%
        relocate(EpiContext, .before=Time) %>%
        relocate(tol_epi, .before=Time )%>%
        relocate(tol_delay, .before=Time)

    # Build dataframe of stats
    Results_SA = rbind(Results_SA,
                       tibble(EpiContext=EpiContext,
                              tol_epi = tol_epi,
                              tol_delay=tol_delay,
                              Mean = mean(SA_tol_epi$Date),
                              Median = median(SA_tol_epi$Date),
                              P025 = quantile(SA_tol_epi$Date,0.025,type=1)[[1]],
                              P975 = quantile(SA_tol_epi$Date,0.975,type=1)[[1]]))
    
    # Build dataframe with all values
    SA_tol_epi_long=rbind(SA_tol_epi_long,SA_tol_epi)
}
SA_tol_epi_long

SA_tol_epi_long$tol_epi=as_factor(SA_tol_epi_long$tol_epi)

SA_tol_epi_long=SA_tol_epi_long %>%
    mutate(Date=Date_N-Time)

## Mean time by tol_epi
mu_tol_epi = ddply(SA_tol_epi_long, "tol_epi", summarise, grp.mean=mean(Time))
me_tol_epi = ddply(SA_tol_epi_long, "tol_epi", summarise, grp.median=median(Time))

#### Plot
p_tol_epi = ggplot(SA_tol_epi_long, aes(x=Date, y=..density..,fill=tol_epi)) +
    geom_histogram(alpha=0.4, 
                   binwidth=1, # each bar represents 1 day
                   position='identity') +
    ## When two values match
    geom_vline(data=me_tol_epi[2,], aes(xintercept=Date_N-grp.median),
               color=ColorTwo) +
    geom_vline(data=me_tol_epi, aes(xintercept=Date_N-grp.median, color=tol_epi),
               linetype="dashed") +
    scale_fill_manual(name=TeX("$delta_{tol}^Y$"), breaks=as.character(rev(tol_epi_var)),values=c(ColorOne,ColorTwo,ColorThree)) +
    scale_color_manual(name=TeX("$delta_{tol}^Y$"), breaks=as.character(rev(tol_epi_var)),values=c(ColorOne,ColorTwo,ColorThree)) +
    scale_y_continuous(expand=c(0.01,0)) +
    scale_x_date(name="\nDate of symptoms onset (2019)",
                 breaks=MyBreaks,
                 date_labels="%b %d") +
    labs(title="",
         x ="Days",
         y=" ") +
    theme_classic() +
    theme(axis.title.x=element_text(vjust=+2),
          axis.text.x=element_text(angle=90, hjust=-0.5, vjust=0.5)) +
    NULL

p_tol_epi

####
#### 2. Varying tol_delay ################################################################
####
tol_epi=0.3

for (tol_delay in tol_delay_var){
    # Strings
    tol_delay_str=as.character(10*tol_delay)
    tol_epi_str=as.character(10*tol_epi)

    # Read data
    SA_tol_delay = read_csv(paste0("Fig_3_SensitivityAnalyses_Calibration/Output/",EpiContext,"/Varying_tol_delay/Cases_EpiSize_Time_tol_epi_0",tol_epi_str,"_tol_delay_0",tol_delay_str,".csv"),
                          col_names=FALSE,
                          col_types = "ddddd") 
    colnames(SA_tol_delay)=c("Time","Cases","EpiSize","Case1")
    
    SA_tol_delay=SA_tol_delay %>%
        mutate(Date=Date_N-Time,
               tol_epi=tol_epi,
               tol_delay=tol_delay) %>%
        as_tibble()
    
    # Build dataframe of stats
    Results_SA = rbind(Results_SA,
                       tibble(EpiContext=EpiContext,
                              tol_epi = tol_epi,
                              tol_delay=tol_delay,
                              Mean = mean(SA_tol_delay$Date),
                              Median = median(SA_tol_delay$Date),
                              P025 = quantile(SA_tol_delay$Date,0.025,type=1)[[1]],
                              P975 = quantile(SA_tol_delay$Date,0.975,type=1)[[1]]))

    # Build dataframe
    SA_tol_delay_long=rbind(SA_tol_delay_long,SA_tol_delay)
}
SA_tol_delay_long

SA_tol_delay_long$tol_delay=as_factor(SA_tol_delay_long$tol_delay)

## Mean time by tol_epi
mu_tol_delay <- ddply(SA_tol_delay_long, "tol_delay", summarise, grp.mean=mean(Time))
me_tol_delay <- ddply(SA_tol_delay_long, "tol_delay", summarise, grp.median=median(Time))

SA_tol_delay_long=SA_tol_delay_long %>%
    mutate(Date=Date_N-Time)

##
## Plot
##
p_tol_delay = ggplot(SA_tol_delay_long, aes(x=Date, y=..density..,fill=tol_delay)) +
    geom_histogram(alpha=0.4,
                   binwidth=1, # each bar represents 1 day
                   position='identity') +
    geom_vline(data=me_tol_delay, aes(xintercept=Date_N-grp.median, color=tol_delay),
               linetype="dashed") +
    scale_fill_manual(name=TeX("$delta_{tol}^{tau}$"), breaks=as.character(rev(tol_delay_var)),values=c(ColorOne,ColorTwo,ColorThree)) +
    scale_color_manual(name=TeX("$delta_{tol}^{tau}$"), breaks=as.character(rev(tol_delay_var)),values=c(ColorOne,ColorTwo,ColorThree)) +
    scale_y_continuous(expand=c(0.01,0)) +
    scale_x_date(name=xAxisTitle,
                 limits=MyLimits,
                 breaks=MyBreaks,
                 date_labels="%b %d") +
    labs(title="",
         x ="Days",
         y=" ") +
    theme_classic() +
    theme(axis.title.x=element_text(vjust=+2),
          axis.text.x=element_text(angle=90, hjust=-0.5, vjust=0.5)) +
    NULL
p_tol_delay


##
## Put vertically
##

p_tol_delay_V = p_tol_delay +
    scale_y_continuous(expand=c(0.01,0)) +
    scale_x_date(name="\n",
                 limits=x_lims,
                 breaks=MyBreaks,
                 date_labels="%b %d") +
    ggtitle(Ytitle) +
    theme_classic() +
    theme(plot.title=element_text(size=14), # face="bold"
          axis.text.x=element_text(angle=90, hjust=-0.5, vjust=0.5)) +
    NULL

tol_epi_V = p_tol_epi + 
    scale_y_continuous(expand=c(0.01,0)) +
    scale_x_date(name=xAxisTitle,
                 limits=x_lims,
                 breaks=MyBreaks,
                 date_labels="%b %d") +
    ggtitle(Ytitle) +
    theme_classic() +
    theme(plot.title=element_text(size=14), #face="bold"
          axis.text.x=element_text(angle=90, hjust=-0.5, vjust=0.5)) +
    NULL

p_SA=plot_grid(p_tol_delay_V, tol_epi_V,
                  labels=MyPanelLabels,
                  ncol=1, nrow=2)
p_SA


####
#### Display results ####################
####

ggsave(paste0("Fig_3_SensitivityAnalyses_Calibration/Output/",EpiContext,"/SensitivityAnalyses_Conditions_",EpiContext,"_V.pdf"),
       plot=p_SA, height=10, width=12, units=c("cm"))

Results_SA 

