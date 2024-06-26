#### Jijón S, Czuppon P, Blanquart F & Débarre F (2023). 
#### Using early detection data to estimate the date of emergence of an epidemic outbreak.
#### https://github.com/sjijon/estimate-emergence-from-data
####
####
#### Plot all results from the sensitivity analyses
####
#### 0. SETUP ################################################################
####

## CLEAR ALL
rm(list=ls())

## PACKAGES 
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(plyr)
library(latex2exp)
library(gridExtra)
library(cowplot)

## Colors
ColorOne = "#00648C" # Blue
ColorTwo = "#FFAA00" # Yellow
ColorThree = "#A01E18" # Red

Cond_Name = "Cond_Cumul"

Results_SA = NULL
str_aux=10000

for(EpiContext in c("Alpha_UK","COVID-19_Wuhan")){
    # Context-specific parameters
    if(EpiContext == "Alpha_UK"){
        # Date of N-th case
        Date_N = as.Date("2020-11-11")
        
        ## Baseline
        p_detect=0.0105
        R0=1.9
        kappa=0.57
        theta_tau=12
        tol_epi=0.3
        
        ## Varying parameters
        Pdetect_var=c(0.005,0.0105,0.025) # Baseline: 0.0105
        R0_var=c(1.7,1.9,2.1)             # Baseline: 1.9
        kappa_var=c(0.35,0.57,0.75)       # Baseline: 0.57
        theta_tau_var=c(7,12,15)          # Baseline: 12
        tol_epi_var=c(0.3,0.5)            # Baseline: 0.3
        
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
        
        x_lims=as.Date(c("2020-05-21", "2020-09-14"))
        y_lims=c(0,0.07)
        Ytitle="Alpha"
        MyXLabel="\nDate of emergence (2020)"
        MyPanelLabels=c(" ","A","C","E","G","I","K")
        StatPos = as.Date(c("2020-05-10","2020-05-14","2020-06-05","2020-06-09"))
    }else if (EpiContext == "COVID-19_Wuhan"){
        # Date of N-th case
        Date_N = as.Date("2020-01-19")
        
        # Baseline
        p_detect=0.15
        R0=2.5
        kappa = 0.1
        theta_tau=6.25
        tol_epi=0.3
        
        # Varying parameters
        Pdetect_var=c(0.10,0.15,0.25)
        R0_var=c(2.0,2.5,3.5)
        kappa_var=c(0.05,0.1,0.25)
        theta_tau_var=c(2.8,6.25,11.6)    # Baseline: 6.25
        tol_epi_var=c(0.1,0.3,0.5)        # Baseline: 0.3
        
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
        
        x_lims=as.Date(c("2019-09-21", "2019-12-10"))
        y_lims=c(0,0.08)
        Ytitle="COVID-19"
        MyXLabel="\nDate of emergence (2019)"
        MyPanelLabels=c(" ","B","D","F","H","J","L")
        StatPos = as.Date(c("2019-10-10","2019-10-14","2019-11-05","2019-11-09"))
    }
    
    ####
    #### 1. Varying R0 ################################################################
    ####
    if(EpiContext == "Alpha_UK"){
        p_detect=0.0105
        kappa=0.57
        theta_tau=12
    }else if (EpiContext == "COVID-19_Wuhan"){
        kappa=0.10
        p_detect=0.15
        theta_tau=6.25
    }
    
    SA_R0_long=NULL
    
    for (R0 in R0_var){
        # Strings
        R0_str=as.character(10*R0)
        kappa_str=as.character(100*kappa)
        p_detect_str=as.character(str_aux*p_detect)
        theta_tau_str=as.character(100*theta_tau)
        
        # Read data
        SA_R0=read.csv(paste0("Fig4_SensitivityAnalyses/Output/",EpiContext,"/Varying_R0/Cases_EpiSize_Time_R0_",R0_str,"_kappa_0",kappa_str,"_p_detect_0",p_detect_str,"_theta_tau",theta_tau_str,".csv"), header=FALSE) 
        colnames(SA_R0)=c("Time","Cases","EpiSize","Case1")
        
        SA_R0=SA_R0 %>%
            mutate(EpiContext=EpiContext,
                   Date=Date_N-Time,
                   R0=R0) %>%
            as_tibble()
        
        # Build dataframe
        SA_R0_long=rbind(SA_R0_long,SA_R0)
        
        # Build dataframe of stats
        Results_SA = rbind(Results_SA,
                           tibble(EpiContext=EpiContext,
                                  R0 = R0,
                                  kappa = kappa,
                                  p_detect=p_detect,
                                  theta_tau=theta_tau,
                                  tol_epi = tol_epi,
                                  Mean = mean(SA_R0$Date),
                                  Median = median(SA_R0$Date),
                                  P025 = quantile(SA_R0$Date,0.025,type=1)[[1]],
                                  P975 = quantile(SA_R0$Date,0.975,type=1)[[1]]))
    }
    
    SA_R0_long$R0=as_factor(SA_R0_long$R0)
    
    ## Mean time by R0
    mu_R0 <- ddply(SA_R0_long, "R0", summarise, grp.mean=mean(Time))
    me_R0 <- ddply(SA_R0_long, "R0", summarise, grp.median=median(Time))
    
    #### Plot
    p_R0 = ggplot(data=SA_R0_long, 
                aes(x=Date, y=R0, fill=R0,color=R0)) +
        geom_violin(width=0.6, size=0.2, alpha=0.3)  + 
        # geom_boxplot(width=0.25,size=0.2,alpha=0,outlier.shape=NA) +
        stat_summary(fun = "mean",
                     geom = "point",
                     shape=5, # Diamond
                     size=3) +
        stat_summary(fun = "median",
                     geom = "point",
                     shape=3, # Cross
                     size=3) +
        ## Add stats to legend
        # annotate("point", shape=5, x=as.Date("2020-05-15"), y=2.4,size=3) +
        # annotate("text",x=as.Date("2020-05-19"), y=2.4,label="Mean",hjust=0,size=2.75) +
        # annotate("point", shape=3, x=as.Date("2020-05-15"), y=2.2,size=3) +
        # annotate("text",x=as.Date("2020-05-19"), y=2.2,label="Median",hjust=0,size=2.75) +
        ## Axis and Labels
        scale_fill_manual(name=TeX("$R$"), breaks=as.character(rev(R0_var)),values=c(ColorTwo,ColorOne,ColorThree)) +
        scale_color_manual(name=TeX("$R$"), breaks=as.character(rev(R0_var)),values=c(ColorTwo,ColorOne,ColorThree)) +
        scale_y_discrete(name="Distributions",
                         label=c("","",""),
                         expand=c(0,0.7)) +
        scale_x_date(name="",
                     breaks=MyBreaks,
                     date_labels="%b %d") +
        ggtitle("") +
        theme_classic() +
        theme(legend.position=c(0.12,0.8),
              plot.title=element_text(size=14), # face="bold"
              axis.text.x=element_text(angle=90, hjust=-0.5, vjust=0.5)) +
        NULL
    # print(p_R0)
    
    ####
    #### 2. Varying kappa ################################################################
    ####
    if(EpiContext == "COVID-19_Wuhan"){
        p_detect=0.15
        R0=2.5
        theta_tau=6.25
    }else if (EpiContext == "Alpha_UK"){
        p_detect=0.0105
        R0=1.9
        theta_tau=12
    }

    SA_kappa_long=NULL
    for (kappa in kappa_var){
        # Strings
        R0_str=as.character(10*R0)
        kappa_str=as.character(100*kappa)
        p_detect_str=as.character(str_aux*p_detect)
        theta_tau_str=as.character(100*theta_tau)

        # Read data
        SA_kappa=read.csv(paste0("Fig4_SensitivityAnalyses/Output/",EpiContext,"/Varying_kappa/Cases_EpiSize_Time_R0_",R0_str,"_kappa_0",kappa_str,"_p_detect_0",p_detect_str,"_theta_tau",theta_tau_str,".csv"),
                          header=FALSE)
        colnames(SA_kappa)=c("Time","Cases","EpiSize","Case1")

        SA_kappa=SA_kappa %>%
            mutate(EpiContext=EpiContext,
                   Date=Date_N-Time,
                   kappa = kappa) %>%
            as_tibble()

        # Build dataframe of stats
        Results_SA = rbind(Results_SA,
                           tibble(EpiContext=EpiContext,
                                  R0 = R0,
                                  kappa = kappa,
                                  p_detect=p_detect,
                                  theta_tau=theta_tau,
                                  tol_epi = tol_epi,
                                  Mean = mean(SA_kappa$Date),
                                  Median = median(SA_kappa$Date),
                                  P025 = quantile(SA_kappa$Date,0.025,type=1)[[1]],
                                  P975 = quantile(SA_kappa$Date,0.975,type=1)[[1]]))

        # Build dataframe
        SA_kappa_long=rbind(SA_kappa_long,SA_kappa)
    }

    SA_kappa_long$kappa=as_factor(SA_kappa_long$kappa)

    ## Mean time by kappa
    mu_kappa <- ddply(SA_kappa_long, "kappa", summarise, grp.mean=mean(Time))
    me_kappa <- ddply(SA_kappa_long, "kappa", summarise, grp.median=median(Time))

    #### Plot
    p_kappa = ggplot(data=SA_kappa_long, 
                   aes(x=Date, y=kappa, fill=kappa,color=kappa)) +
        geom_violin(width=0.6, size=0.2, alpha=0.3)  + 
        # geom_boxplot(width=0.25,size=0.2,alpha=0,outlier.shape=NA) +
        stat_summary(fun = "mean",
                     geom = "point",
                     shape=5, # Diamond
                     size=3) +
        stat_summary(fun = "median",
                     geom = "point",
                     shape=3, # Cross
                     size=3) +
        ## Add stats to legend
        # annotate("point", shape=5, x=as.Date("2020-05-17"), y=2.4,size=3) +
        # annotate("text",x=as.Date("2020-05-21"), y=2.4,label="Mean",hjust=0,size=2.75) +
        # annotate("point", shape=3, x=as.Date("2020-05-17"), y=2.2,size=3) +
        # annotate("text",x=as.Date("2020-05-21"), y=2.2,label="Median",hjust=0,size=2.75) +
        ## Axis and Labels
        scale_fill_manual(name=TeX("$kappa$"), breaks=as.character(rev(kappa_var)),values=c(ColorTwo,ColorOne,ColorThree)) +
        scale_color_manual(name=TeX("$kappa$"), breaks=as.character(rev(kappa_var)),values=c(ColorTwo,ColorOne,ColorThree)) +
        scale_y_discrete(name="Distributions",
                         label=c("","",""),
                         expand=c(0,0.7)) +
        scale_x_date(name="",
                     breaks=MyBreaks,
                     date_labels="%b %d") +
        ggtitle("") +
        theme_classic() +
        theme(legend.position=c(0.12,0.8),
              plot.title=element_text(size=14), # face="bold"
              axis.text.x=element_text(angle=90, hjust=-0.5, vjust=0.5)) +
        NULL
    # print(p_kappa)

    ####
    #### 3. Varying p_detect ################################################################
    ####

    if(EpiContext == "COVID-19_Wuhan"){
        R0=2.5
        kappa = 0.1
        theta_tau=6.25
    }else if (EpiContext == "Alpha_UK"){
        R0=1.9
        kappa=0.57
        theta_tau=12
    }

    SA_p_detect_long=NULL

    for (p_detect in Pdetect_var){
        # Strings
        R0_str=as.character(10*R0)
        kappa_str=as.character(100*kappa)
        p_detect_str=as.character(str_aux*p_detect)
        theta_tau_str=as.character(100*theta_tau)

        # Read data
        SA_p_detect=read.csv(paste0("Fig4_SensitivityAnalyses/Output/",EpiContext,"/Varying_pdetect/Cases_EpiSize_Time_R0_",R0_str,"_kappa_0",kappa_str,"_p_detect_0",p_detect_str,"_theta_tau",theta_tau_str,".csv"),
                             header=FALSE)
        colnames(SA_p_detect)=c("Time","Cases","EpiSize","Case1")

        SA_p_detect=SA_p_detect%>%
            mutate(ProbDetect=p_detect,
                   EpiContext=EpiContext,
                   Date=Date_N-Time) %>%
            as_tibble()

        # Build dataframe with all values
        SA_p_detect_long=rbind(SA_p_detect_long,SA_p_detect)

        # Build dataframe of stats
        Results_SA = rbind(Results_SA,
                           tibble(EpiContext=EpiContext,
                                  R0 = R0,
                                  kappa = kappa,
                                  p_detect=p_detect,
                                  theta_tau=theta_tau,
                                  tol_epi = tol_epi,
                                  Mean = mean(SA_p_detect$Date),
                                  Median = median(SA_p_detect$Date),
                                  P025 = quantile(SA_p_detect$Date,0.025,type=1)[[1]],
                                  P975 = quantile(SA_p_detect$Date,0.975,type=1)[[1]]))

    }

    SA_p_detect_long$ProbDetect=as_factor(SA_p_detect_long$ProbDetect)

    ## Mean time by p_detect
    mu_pdet <- ddply(SA_p_detect_long, "ProbDetect", summarise, grp.mean=mean(Time))
    me_pdet <- ddply(SA_p_detect_long, "ProbDetect", summarise, grp.median=median(Time))


    #### Plot
    p_pdetect = ggplot(data=SA_p_detect_long,
                     aes(x=Date, y=ProbDetect, fill=ProbDetect,color=ProbDetect)) +
        geom_violin(width=0.6, size=0.2, alpha=0.3)  +
        # geom_boxplot(width=0.25,size=0.2,alpha=0,outlier.shape=NA) +
        stat_summary(fun = "mean",
                     geom = "point",
                     shape=5, # Diamond
                     size=3) +
        stat_summary(fun = "median",
                     geom = "point",
                     shape=3, # Cross
                     size=3) +
        ## Add stats to legend
        # annotate("point", shape=5, x=as.Date("2020-06-02"), y=2.4,size=3) +
        # annotate("text",x=as.Date("2020-06-04"), y=2.4,label="Mean",hjust=0,size=2.75) +
        # annotate("point", shape=3, x=as.Date("2020-06-02"), y=2.2,size=3) +
        # annotate("text",x=as.Date("2020-06-04"), y=2.2,label="Median",hjust=0,size=2.75) +
        ## Axis and Labels
            scale_fill_manual(name=TeX("$p_{detect}$"), breaks=as.character(rev(Pdetect_var)),values=c(ColorTwo,ColorOne,ColorThree)) +
            scale_color_manual(name=TeX("$p_{detect}$"), breaks=as.character(rev(Pdetect_var)),values=c(ColorTwo,ColorOne,ColorThree)) +
        scale_y_discrete(name="Distributions",
                         label=c("","",""),
                         expand=c(0,0.7)) +
        scale_x_date(name="",
                     breaks=MyBreaks,
                     date_labels="%b %d") +
        ggtitle("") +
        theme_classic() +
        theme(legend.position=c(0.12,0.8),
              plot.title=element_text(size=14), # face="bold"
              axis.text.x=element_text(angle=90, hjust=-0.5, vjust=0.5)) +
        NULL

    # print(p_pdetect)

    ####
    #### 4. Varying theta_tau ################################################################
    ####

    if(EpiContext == "COVID-19_Wuhan"){
      R0=2.5
      kappa = 0.1
      p_detect=0.15
    }else if (EpiContext == "Alpha_UK"){
      R0=1.9
      kappa=0.57
      p_detect=0.0105
    }

    SA_theta_tau_long=NULL

    for (theta_tau in theta_tau_var){
      # Strings
      R0_str=as.character(10*R0)
      kappa_str=as.character(100*kappa)
      p_detect_str=as.character(str_aux*p_detect)
      theta_tau_str=as.character(100*theta_tau)

      # Read data
      SA_theta_tau=read.csv(paste0("Fig4_SensitivityAnalyses/Output/",EpiContext,"/Varying_theta_tau/Cases_EpiSize_Time_R0_",R0_str,"_kappa_0",kappa_str,"_p_detect_0",p_detect_str,"_theta_tau",theta_tau_str,".csv"),
                             header=FALSE)
      colnames(SA_theta_tau)=c("Time","Cases","EpiSize","Case1")

      SA_theta_tau=SA_theta_tau%>%
        mutate(TimeToDetect=theta_tau,
               EpiContext=EpiContext,
               Date=Date_N-Time) %>%
        as_tibble()

      # Build dataframe with all values
      SA_theta_tau_long=rbind(SA_theta_tau_long,SA_theta_tau)

      # Build dataframe of stats
      Results_SA = rbind(Results_SA,
                         tibble(EpiContext=EpiContext,
                                R0 = R0,
                                kappa = kappa,
                                p_detect=p_detect,
                                theta_tau=theta_tau,
                                tol_epi = tol_epi,
                                Mean = mean(SA_theta_tau$Date),
                                Median = median(SA_theta_tau$Date),
                                P025 = quantile(SA_theta_tau$Date,0.025,type=1)[[1]],
                                P975 = quantile(SA_theta_tau$Date,0.975,type=1)[[1]]))

    }

    SA_theta_tau_long$TimeToDetect=as_factor(SA_theta_tau_long$TimeToDetect)

    ## Mean time by theta_tau
    mu_theta_tau <- ddply(SA_theta_tau_long, "TimeToDetect", summarise, grp.mean=mean(Time))
    me_theta_tau <- ddply(SA_theta_tau_long, "TimeToDetect", summarise, grp.median=median(Time))


    #### Plot
    p_theta_tau = ggplot(data=SA_theta_tau_long,
                       aes(x=Date, y=TimeToDetect, fill=TimeToDetect,color=TimeToDetect)) +
      geom_violin(width=0.6, size=0.2, alpha=0.3)  +
      # geom_boxplot(width=0.25,size=0.2,alpha=0,outlier.shape=NA) +
      stat_summary(fun = "mean",
                   geom = "point",
                   shape=5, # Diamond
                   size=3) +
      stat_summary(fun = "median",
                   geom = "point",
                   shape=3, # Cross
                   size=3) +
      ## Add stats to legend
      # annotate("point", shape=5, x=as.Date("2020-06-02"), y=2.4,size=3) +
      # annotate("text",x=as.Date("2020-06-04"), y=2.4,label="Mean",hjust=0,size=2.75) +
      # annotate("point", shape=3, x=as.Date("2020-06-02"), y=2.2,size=3) +
      # annotate("text",x=as.Date("2020-06-04"), y=2.2,label="Median",hjust=0,size=2.75) +
      ## Axis and Labels
      scale_fill_manual(name=TeX("$\\theta_\\tau$"), breaks=as.character(rev(theta_tau_var)),values=c(ColorTwo,ColorOne,ColorThree)) +
      scale_color_manual(name=TeX("$\\theta_\\tau$"), breaks=as.character(rev(theta_tau_var)),values=c(ColorTwo,ColorOne,ColorThree)) +
      scale_y_discrete(name="Distributions",
                       label=c("","",""),
                       expand=c(0,0.7)) +
      scale_x_date(name="",
                   breaks=MyBreaks,
                   date_labels="%b %d") +
      ggtitle("") +
      theme_classic() +
      theme(legend.position=c(0.12,0.8),
            plot.title=element_text(size=14), # face="bold"
            axis.text.x=element_text(angle=90, hjust=-0.5, vjust=0.5)) +
      NULL

    # print(p_theta_tau)

    ####
    #### 5. Varying tol_epi ################################################################
    ####
    if(EpiContext == "Alpha_UK"){
        R0=1.9
        kappa=0.57
        p_detect=0.0105
        theta_tau=12
    }else if (EpiContext == "COVID-19_Wuhan"){
        R0=2.5
        kappa = 0.1
        p_detect=0.15
        theta_tau=6.25
    }

    SA_tol_epi_long=NULL

    for (tol_epi in tol_epi_var){
        # Strings
        tol_epi_str=as.character(10*tol_epi)

        # Read data
        SA_tol_epi=read_csv(paste0("Fig4_SensitivityAnalyses/Output/",EpiContext,"/Varying_tol_epi/Cases_EpiSize_Time_tol_epi_0",tol_epi_str,".csv"),
                            col_names=FALSE,
                            col_types = "ddddd")
        colnames(SA_tol_epi)=c("Time","Cases","EpiSize","Case1")

         if (tol_epi==0.3){
            SA_tol_epi = SA_tol_epi[2:5000,] # ignore first row
        }

        SA_tol_epi=SA_tol_epi%>%
            mutate(EpiContext=EpiContext,
                   Date=Date_N-Time,
                   tol_epi=tol_epi) %>%
            as_tibble()

        SA_tol_epi = SA_tol_epi %>%
            relocate(EpiContext, .before=Time) %>%
            relocate(tol_epi, .before=Time )

        # Build dataframe of stats
        Results_SA = rbind(Results_SA,
                           tibble(EpiContext=EpiContext,
                                  R0 = R0,
                                  kappa = kappa,
                                  p_detect=p_detect,
                                  theta_tau=theta_tau,
                                  tol_epi = tol_epi,
                                  Mean = mean(SA_tol_epi$Date),
                                  Median = median(SA_tol_epi$Date),
                                  P025 = quantile(SA_tol_epi$Date,0.025,type=1,na.rm=TRUE)[[1]],
                                  P975 = quantile(SA_tol_epi$Date,0.975,type=1,na.rm=TRUE)[[1]]))

        # Build dataframe with all values
        SA_tol_epi_long=rbind(SA_tol_epi_long,SA_tol_epi)
    }
    SA_tol_epi_long

    SA_tol_epi_long$tol_epi=as_factor(SA_tol_epi_long$tol_epi)

    ## Mean time by tol_epi
    mu_tol_epi = ddply(SA_tol_epi_long, "tol_epi", summarise, grp.mean=mean(Time))
    me_tol_epi = ddply(SA_tol_epi_long, "tol_epi", summarise, grp.median=median(Time))

    #### Plot
    p_tol_epi = ggplot(data=SA_tol_epi_long,
                       aes(x=Date, y=tol_epi, fill=tol_epi,color=tol_epi)) +
        geom_violin(width=0.6, size=0.2, alpha=0.3)  +
        # geom_boxplot(width=0.25,size=0.2,alpha=0,outlier.shape=NA) +
        stat_summary(fun = "mean",
                     geom = "point",
                     shape=5, # Diamond
                     size=3) +
        stat_summary(fun = "median",
                     geom = "point",
                     shape=3, # Cross
                     size=3) +
        ## Add stats to legend
        # annotate("point", shape=5, x=as.Date("2020-06-02"), y=2.4,size=3) +
        # annotate("text",x=as.Date("2020-06-04"), y=2.4,label="Mean",hjust=0,size=2.75) +
        # annotate("point", shape=3, x=as.Date("2020-06-02"), y=2.2,size=3) +
        # annotate("text",x=as.Date("2020-06-04"), y=2.2,label="Median",hjust=0,size=2.75) +
        ## Axis and Labels
        scale_fill_manual(name=TeX("$delta_{tol}^Y$"), breaks=as.character(rev(tol_epi_var)),values=c(ColorTwo,ColorOne,ColorThree)) +
        scale_color_manual(name=TeX("$delta_{tol}^Y$"), breaks=as.character(rev(tol_epi_var)),values=c(ColorTwo,ColorOne,ColorThree)) +
        scale_y_discrete(name="Distributions",
                         label=c("","",""),
                         expand=c(0,0.7)) +
        scale_x_date(name="",
                     breaks=MyBreaks,
                     date_labels="%b %d") +
        ggtitle("") +
        theme_classic() +
        theme(legend.position=c(0.12,0.8),
              plot.title=element_text(size=14), # face="bold"
              axis.text.x=element_text(angle=90, hjust=-0.5, vjust=0.5)) +
        NULL

    # print(p_tol_epi)

    
    ####
    #### Put in one column and save ################################################################
    ####
    title = ggdraw() +
        draw_label(Ytitle,
            fontface='bold',x=0,hjust=0) +
        theme(plot.margin = margin(0,0,0,24)) # align title

    p_draft=plot_grid(title,p_R0,p_kappa,p_pdetect,p_theta_tau,p_tol_epi,
                      labels=MyPanelLabels,
                      ncol=1, nrow=7,
                      rel_heights = c(.1,1,1,1,1,1,1))
    print(p_draft)

    ggsave(paste0("Fig4_SensitivityAnalyses/Output/Fig4_SensitivityAnalyses_",EpiContext,".pdf"),
           plot=p_draft, height=30, width=10, units=c("cm"))

    ####
    #### Save individual figs ################### ########################################
    ####
    filedir=paste0("Fig4_SensitivityAnalyses/Output/",EpiContext,"/")
    ggsave(paste0(filedir,"Varying_R0/Fig4_SensitivityAnalyses_R0.pdf"),
           plot=p_R0, height=10, width=12, units=c("cm"))
    ggsave(paste0(filedir,"Varying_kappa/Fig4_SensitivityAnalyses_kappa.pdf"),
           plot=p_kappa, height=10, width=12, units=c("cm"))
    ggsave(paste0(filedir,"Varying_pdetect/Fig4_SensitivityAnalyses_p_detect.pdf"),
           plot=p_pdetect, height=10, width=12, units=c("cm"))
    ggsave(paste0(filedir,"Varying_tol_epi/Fig4_SensitivityAnalyses_tol_epi.pdf"),
           plot=p_tol_epi, height=10, width=12, units=c("cm"))
}