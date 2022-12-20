####
#### Early epidemics
#### iEES
####
#### June, 2022
####
#### Read and plot results
####
#### 0. SETUP ################################################################
####

## CLEAR ALL
rm(list = ls())

## PACKAGES 
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(ggh4x)

## Colors
MyGray = 	"#A0A0A0" # Gray

## Simulations
# num_sims = 5000
num_sims = 10

# sim_num = sample(1:num_sims, 1) # One random sim
sim_num = 1

dir.create("Figs_Appendix/FigS3-S4_CumulativeCases/Output",showWarnings = FALSE)
###
#### 1. READ DATA ################################################################
####
for (EpiContext in c("Alpha_UK","COVID-19_Wuhan")){
    ####
    #### Observations
    ####
    if (EpiContext == "Alpha_UK"){
        DateN = as.Date("2020-11-11")
        date_breaks = seq(as.Date("2020-07-01"), as.Date("2020-11-30"), by = "10 days")
        
        Obs_Cases = read.csv("Data/CaseData_Alpha_UK.csv") %>%
            as_tibble()
        
        Obs_Cases = Obs_Cases[Obs_Cases$Date<="2020-11-11",]
        
        Dates_Ticks=c(seq(as.Date("2020-05-07"),as.Date("2020-05-28"),by="7 days"),
                      as.Date("2020-06-01"),
                      seq(as.Date("2020-06-07"),as.Date("2020-06-28"),by="7 days"),
                      as.Date("2020-07-01"),
                      seq(as.Date("2020-07-07"),as.Date("2020-07-28"),by="7 days"),
                      # as.Date("2020-07-26"),
                      as.Date("2020-08-01"),
                      seq(as.Date("2020-08-07"),as.Date("2020-08-28"),by="7 days"),
                      # as.Date("2020-08-26"),
                      as.Date("2020-09-01"),
                      seq(as.Date("2020-09-07"),as.Date("2020-09-28"),by="7 days"),
                      # as.Date("2020-09-26"),
                      as.Date("2020-10-01"),
                      seq(as.Date("2020-10-07"),as.Date("2020-10-28"),by="7 days"),
                      # as.Date("2020-10-26"),
                      as.Date("2020-11-01"),
                      as.Date("2020-11-07"),
                      as.Date("2020-11-11"))
        
        Dates_Ticks_minor=seq(as.Date("2020-05-01"),as.Date("2020-11-20"),by="1 days")
        
        DatesLimits = as.Date(c("2020-08-15","2020-11-30"))
        Xtitle = "\nDate of sample collection (2020)"
        
    }else if (EpiContext == "COVID-19_Wuhan"){
        DateN = as.Date("2020-01-19")
        date_breaks = seq(as.Date("2019-11-15"), as.Date("2020-01-10"), by = "5 days")
        
        Obs_Cases = read.csv("Data/CaseData_COVID-19_Wuhan.csv") %>%
            as_tibble()
        Obs_Cases = Obs_Cases[Obs_Cases$Date<="2020-01-19",]
        
        Dates_Ticks=c(seq(as.Date("2019-09-14"),as.Date("2019-09-28"),by="7 days"),
                      as.Date("2019-10-01"),
                      seq(as.Date("2019-10-07"),as.Date("2019-10-28"),by="7 days"),
                      as.Date("2019-11-01"),
                      seq(as.Date("2019-11-07"),as.Date("2019-11-28"),by="7 days"),
                      as.Date("2019-12-01"),
                      seq(as.Date("2019-12-07"),as.Date("2019-12-28"),by="7 days"),
                      as.Date("2020-01-01"),
                      seq(as.Date("2020-01-07"),as.Date("2020-01-28"),by="7 days"))
        
        Dates_Ticks_minor=seq(as.Date("2019-09-10"),as.Date("2020-01-28"),by="1 days")
        
        DatesLimits = as.Date(c("2019-11-21","2020-01-22"))
        Xtitle = "\nDate of symptoms onset (2019--2020)"
        
    }
    
    Obs_Cases$Date = as.Date(Obs_Cases$Date)
    
    ## Last observed case
    Date_Obs_NthCase = tail(Obs_Cases$Date,1)
    ## Number of cases
    N_Cases = sum(Obs_Cases$Cases)
    
    ####
    #### 2. READ RESULTS ################################################################
    ####
    data_allsims_cases = read.csv(paste0("Fig1-2_Emergence/Output/",EpiContext,"/cumul_",N_Cases,"cases.csv"),
                                  header = TRUE)
    data_allsims_cases
    
    ####
    #### 4. PLOT ALL ################################################################
    ####
    # Indx_Sims = seq(1,num_sims,20)        # Figure in the Appendix
    Indx_Sims = seq(1,num_sims,500)       # 100 sims
    # Indx_Sims = sample(1:num_sims, 1) # One random sim
    
    p_Nsims <- ggplot()
    
    for (sim_num in Indx_Sims){
        ## read data
        data_sim_cases = t(data_allsims_cases[sim_num,])
        data_sim_cases = data_sim_cases[!is.na(data_sim_cases)]
        data_sim_cases
        
        sim_cases_daynum = seq(1,length(data_sim_cases),1)
        sim_cases_daynum = sim_cases_daynum[data_sim_cases>0]
        data_sim_cases = data_sim_cases[data_sim_cases>0]
        
        ## Build table
        Sim_Cases = cbind(data_sim_cases,sim_cases_daynum) %>%
            as_tibble()
        colnames(Sim_Cases) = c("CumulCases","DayNumber")
        
        Sim_Dates = NULL
        n = length(sim_cases_daynum)
        for (d in seq(1,n,1)){
            Sim_Dates = rbind(Sim_Dates,
                              Date_Obs_NthCase - Sim_Cases$DayNumber[n] + Sim_Cases$DayNumber[d])
        }
        
        Sim_Cases =  Sim_Cases %>%
            mutate(Date = as.Date(Sim_Dates,origin="1970-01-01")) %>%
            relocate(Date,.before="CumulCases")
        
        ## Plot line
        p_Nsims = p_Nsims + 
            geom_line(data=Sim_Cases, aes(x=Date, y=CumulCases, color="Simulated"), lwd=0.2) +
            # geom_point(data=Sim_Cases, aes(x=Date, y=CumulCases, color="Simulated"), shape=1) +
            NULL
    }
    
    p_Nsims = p_Nsims + 
        geom_line(data=Obs_Cases, aes(x=Date, y=Cumul, color="Observed"), lwd=0.5) +
        geom_point(data=Obs_Cases, aes(x=Date, y=Cumul, color="Observed"),size=1) +
        scale_color_manual(name="", values = c("black",MyGray)) + 
        scale_x_date(name=Xtitle,
                     limits=DatesLimits,
                     breaks=Dates_Ticks,
                     minor_breaks=Dates_Ticks_minor,
                     guide="axis_minor",
                     date_labels="%b %d",
                     expand=c(0, 0)) +
        labs(title="",
             y = "Cumulative number of cases") +
        theme_classic() +
        theme(legend.position=c(0.15,0.8), 
              axis.text.x=element_text(angle=90, hjust=1, vjust=1),
              axis.ticks.length=unit(6, "pt")) +
        NULL
    

    p_Nsims
    
    #### save figure
    ggsave(paste0("Figs_Appendix/FigS3-S4_CumulativeCases/Output/CumulativeCases_",EpiContext,".pdf"),
           plot=p_Nsims, width=14, height=10, units=c("cm"), dpi=600)
}