## Jijón S, Czuppon P, Blanquart F & Débarre F (2023). 
## Using early detection data to estimate the date of emergence of an epidemic outbreak.
## https://github.com/sjijon/estimate-emergence-from-data
##
##
## Parametrizing the model
## 
## Infectious disease: COVID-19
## Location: Wuhan
## Period: December 8-31, 2019
## Data: Simulated data
#################################################

##
## 1. Read simulated data #################################
##
data = readdlm(string("Data/SimulatedData/COVID-19_Wuhan/sims/simdata_det_sim", dataset_num_,".csv"),'\t',skipstart=0)

## Add cumsum column
data = [data cumsum(data[:,2])]

# data = data[data[:,3].<=1100 ,:]## First 1000 obs_cases_cumul only
# data = data[data[:,3].<=120 ,:] ## First 100 obs_cases_cumul only
# data = data[data[:,3].<=12 ,:] ## First ~10 obs_cases_cumul only
# data = data[data[:,3].<2 ,:]  ## First case only (selecting the first row only)


# Number of cases
N_cases = sum(data[:,2])

Date_1 = Dates.Date(data[1,1])
Date_N = Dates.Date(data[end,1])

##
## 2. Context-specific parameters #################################
##
## Same as in the COVID-19 @ Wuhan application
## All times are given in days

##
## 2.1. Case detection (symptoms onset)
##

## Max. time for detecting the first case
t_detect_max = 90

##  Time from infection to case declaration (clinically)
##  ~ Gamma distribution (mean=6.5,std=2.6) (Backer et al. Eurosurveillance, 2020)
shape_detect = 6.5  # kappa = mean^2/std^2
scale_detect = 1.04 # theta = std^2/mean

##  Probability of detection
p_detect = 0.15             # Acertainement rate (Hao et al., Nature, 2020)

##
## 2.2. Offspring generation
##

## Time to secondary infection
## ~ Gamma distribution (mean = 5.5 days)
shape_inf = 6.6
scale_inf = 0.833

### Number of secondary infections
## ~ Negative binomial distribution
## Mean number of secondary infections
R0 = 2.5
kappa = 0.1                 # overdispersion parameter (Endo et al., Wellcome Open Research, 2020)
p_neg = kappa/(kappa+R0)    # success probability for the negative-binomial

## Stopping criteria
min_n_infect = 5*N_cases/p_detect
# min_n_infect = 2*N_cases/p_detect
# min_n_infect = 100*N_cases/p_detect
max_t_infect = 365                  # Max. time (days) for epidemic developing