#### Jijón S, Czuppon P, Blanquart F & Débarre F (2023). 
#### Using early detection data to estimate the date of emergence of an epidemic outbreak.
#### https://github.com/sjijon/estimate-emergence-from-data
####
#### Plot Estims_Allults from running the model on synthetic data
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

####
#### 1. Read results ################################################################
####
EpiContext="SimulatedData"
Cond_Name = "Cond_Cumul_Delay"

##
## 1.1 Read simulated datasets
##
Datasets = tibble(Dataset=numeric(),
                  DateInf1=date(),
                  DateCase1=date(),
                  DateCase10=date(),
                  DateCase100=date(),
                  DateCase1000=date())

for (i in seq(1,100,1)){
    ## Infections
    data_inf = read.csv(paste0("Data/SimulatedData/COVID-19_Wuhan/sims/simdata_inf_sim", i,".csv"), header=TRUE, sep="\t") %>%
        as_tibble()
    ## Detections
    data_det = read.csv(paste0("Data/SimulatedData/COVID-19_Wuhan/sims/simdata_det_sim", i,".csv"), header=FALSE, sep="\t") %>%
        as_tibble()
    colnames(data_det) = c("Date","Cases")
    data_det = data_det %>%
        mutate(Cumul = cumsum(Cases))
    
    Datasets = Datasets %>%
        add_row(Dataset=i,
                DateInf1=data_inf$Date[1],
                DateCase1=data_det$Date[1],
                DateCase10=data_det[data_det$Cumul>=10,][[1,1]],
                DateCase100=data_det[data_det$Cumul>=100,][[1,1]],
                DateCase1000=data_det[data_det$Cumul>=1000,][[1,1]])
}

# Convert dates to dates
Datasets$DateInf1 =  Datasets$DateInf1 %>% 
    as.Date()
Datasets$DateCase1 =  Datasets$DateCase1 %>% 
    as.Date()
Datasets$DateCase10 =  Datasets$DateCase10 %>% 
    as.Date()
Datasets$DateCase100 =  Datasets$DateCase100 %>% 
    as.Date()
Datasets$DateCase1000 =  Datasets$DateCase1000 %>% 
    as.Date()

# Compute number of days between 1st infection and Nth case
Datasets = Datasets %>%
    mutate(DaysInf1Case10 = as.numeric(Datasets$DateCase10-Datasets$DateInf1),
           DaysInf1Case100 = as.numeric(Datasets$DateCase100-Datasets$DateInf1),
           DaysInf1Case1000 = as.numeric(Datasets$DateCase1000-Datasets$DateInf1))

Datasets

###
### 1.Z Compute error by dataset (median)
###

##
## Using first 1000 cases
##
Estims_All_Case1000_ = read.csv(paste0("Time distribution for N cases/RunSims_to_Ncases/Output/",EpiContext,"/",Cond_Name,"/MinTime_N_EpiSize_Case1000_300sims_1000datasets.csv"), header=TRUE)
Estims_All_Case1000_ = tibble(Estims_All_Case1000_)
Estims_All_Case1000_$Dataset=as_factor(Estims_All_Case1000_$Dataset)

## Keep only first 100 estimates of the first 100 simulated datasets
Estims_All_Case1000 = NULL
for (i in seq(1,100,1)){
    AUX = Estims_All_Case1000_ %>%
        filter(Dataset == paste0("Sim ",i))
    
    Estims_All_Case1000 = rbind(Estims_All_Case1000,AUX[1:100,])
}

Estims_All_Case1000 = Estims_All_Case1000 %>%
    mutate(ObsDaysInf1CaseN=NA,
           ErrorMedian=NA,
           FirstCases = 1000)

for (i in seq(1,100,1)){
    AUX = Estims_All_Case1000 %>%
        filter(Dataset == paste0("Sim ",i))
    AUX$ObsDaysInf1CaseN =Datasets$DaysInf1Case1000[i]
    AUX$ErrorMedian = abs(AUX$MinTime - Datasets$DaysInf1Case1000[i])
    
    Estims_All_Case1000 = Estims_All_Case1000 %>% 
        mutate_at(vars(c("ObsDaysInf1CaseN")), ~ifelse(Dataset==paste0("Sim ",i), AUX$ObsDaysInf1CaseN, .)) %>% 
        mutate_at(vars(c("ErrorMedian")), ~ifelse(Dataset==paste0("Sim ",i), AUX$ErrorMedian, .))
}
Estims_All_Case1000

##
## Using data on first 100 cases
##
Estims_All_Case100 = read.csv(paste0("Time distribution for N cases/RunSims_to_Ncases/Output/",EpiContext,"/",Cond_Name,"/MinTime_N_EpiSize_Case100_100sims_100datasets.csv"), header=TRUE)
Estims_All_Case100 = tibble(Estims_All_Case100)%>%
    mutate(ObsDaysInf1CaseN=NA,
           ErrorMedian=NA,
           FirstCases = 100)
Estims_All_Case100$Dataset=as_factor(Estims_All_Case100$Dataset)

for (i in seq(1,100,1)){
    AUX = Estims_All_Case100 %>%
        filter(Dataset == paste0("Sim ",i))
    AUX$ObsDaysInf1CaseN =Datasets$DaysInf1Case100[i]
    AUX$ErrorMedian = abs(AUX$MinTime - Datasets$DaysInf1Case100[i])
    
    Estims_All_Case100 = Estims_All_Case100 %>% 
        mutate_at(vars(c("ObsDaysInf1CaseN")), ~ifelse(Dataset==paste0("Sim ",i), AUX$ObsDaysInf1CaseN, .)) %>% 
        mutate_at(vars(c("ErrorMedian")), ~ifelse(Dataset==paste0("Sim ",i), AUX$ErrorMedian, .))
}
Estims_All_Case100

##
## Using data on first 10 cases
##
Estims_All_Case10 = read.csv(paste0("Time distribution for N cases/RunSims_to_Ncases/Output/",EpiContext,"/",Cond_Name,"/MinTime_N_EpiSize_Case10_100sims_100datasets.csv"), header=TRUE)
Estims_All_Case10 = tibble(Estims_All_Case10)%>%
    mutate(ObsDaysInf1CaseN=NA,
           ErrorMedian=NA,
           FirstCases = 10)
Estims_All_Case10$Dataset=as_factor(Estims_All_Case10$Dataset)

for (i in seq(1,100,1)){
    AUX = Estims_All_Case10 %>%
        filter(Dataset == paste0("Sim ",i))
    AUX$ObsDaysInf1CaseN =Datasets$DaysInf1Case10[i]
    AUX$ErrorMedian = abs(AUX$MinTime - Datasets$DaysInf1Case10[i])
    
    Estims_All_Case10 = Estims_All_Case10 %>% 
        mutate_at(vars(c("ObsDaysInf1CaseN")), ~ifelse(Dataset==paste0("Sim ",i), AUX$ObsDaysInf1CaseN, .)) %>% 
        mutate_at(vars(c("ErrorMedian")), ~ifelse(Dataset==paste0("Sim ",i), AUX$ErrorMedian, .))
}
Estims_All_Case10

# ##
# ## Using data on first case
# ##
# Estims_All_Case1 = read.csv(paste0("Time distribution for N cases/RunSims_to_Ncases/Output/",EpiContext,"/",Cond_Name,"/MinTime_N_EpiSize_Case1_100sims_100datasets.csv"), header=TRUE)
# Estims_All_Case1 = tibble(Estims_All_Case1)%>%
#     mutate(ObsDaysInf1CaseN=NA,
#            ErrorMedian=NA,
#            FirstCases = 10)
# Estims_All_Case10$Dataset=as_factor(Estims_All_Case10$Dataset)
# 
# for (i in seq(1,100,1)){
#     AUX = Estims_All_Case1 %>%
#         filter(Dataset == paste0("Sim ",i))
#     AUX$ObsDaysInf1CaseN =Datasets$DaysInf1Case1[i]
#     AUX$ErrorMedian = abs(AUX$MinTime - Datasets$DaysInf1Case10[i])
#     
#     Estims_All_Case1 = Estims_All_Case1 %>% 
#         mutate_at(vars(c("ObsDaysInf1CaseN")), ~ifelse(Dataset==paste0("Sim ",i), AUX$ObsDaysInf1CaseN, .)) %>% 
#         mutate_at(vars(c("ErrorMedian")), ~ifelse(Dataset==paste0("Sim ",i), AUX$ErrorMedian, .))
# }
# Estims_All_Case1
##
## Build array
##
All_Estims = rbind(#Estims_All_Case1,
                   Estims_All_Case10,
                   # Estims_All_Case100,
                   Estims_All_Case1000)
All_Estims$FirstCases = as.factor(All_Estims$FirstCases)
All_Estims = All_Estims %>%
    mutate(RelError=ErrorMedian/ObsDaysInf1CaseN)
All_Estims

####
#### 3. PLOT ERROR ##########################
####


p_error = ## Plot asistogram
    # ggplot(data=All_Estims, aes(x=RelError, y=..density.., color = FirstCases)) +
    # geom_histogram(alpha=0.5, position = 'identity',binwidth = 0.1, fill=NA) +
    ## Step function
    ggplot(data=All_Estims, aes(x=RelError, color = FirstCases)) + 
    stat_bin(aes(y=..density..),geom="step", position = 'identity',binwidth = 0.1, origin=0) +
    scale_color_manual(values=c(ColorOne,ColorThree,ColorTwo),
                      name="Number of cases (N)",
                      breaks=c("1000","100","10"),
                      labels=c("~1000","~100","~10")) +
    scale_x_continuous(name="Relative error (median of estimates vs. observation)",
                       limits = c(0,3),
                       expand=c(0,0)) +
    scale_y_continuous(name="",
                       expand=c(0,0)) +
    theme_classic() +
    NULL

p_error = ggplot(data=All_Estims, aes(x=RelError, color = FirstCases)) + 
    stat_bin(aes(y=..density..),geom="step", position = 'identity',binwidth = 0.1, origin=0) +
    scale_color_manual(values=c(ColorOne,ColorThree),
                       name="Number of cases (N)",
                       breaks=c("1000","10"),
                       labels=c("~1000","~10")) +
    scale_x_continuous(name="Relative error (median of estimates vs. observation)",
                       limits = c(0,2),
                       expand=c(0,0)) +
    scale_y_continuous(name="",
                       limits = c(0,3.1),
                       expand=c(0,0)) +
    theme_classic() +
    NULL

p_error

####
#### 4. SAVE ################################################################
####

print(p_error)

ggsave("Communicating results/Articles/Figures/SupplementaryFigures/RunOnSyntheticData_Error.pdf",
       plot=p_error, height=9, width=14, units=c("cm"))
