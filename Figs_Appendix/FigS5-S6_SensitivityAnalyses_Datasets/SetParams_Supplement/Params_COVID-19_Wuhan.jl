## Parametrizing the model
## 
## Infectious disease: COVID-19
## Location: Wuhan
## Period: Dec 8-Jan 19, 2019
## Data: date of symptoms onset (N=3072)
#################################################

##
## 1. Read data #################################
##

## Oserved data
## Date and number of cases (symptoms onset)
data = readdlm("Data/CaseData_COVID-19_Wuhan.csv",',',skipstart=1)
data = data = data[data[:,1].<="2020-01-19",:] # Up to Jan 19 2020

## Select a dataset:
#### i) Run on data up to Dec 31, 2019
data = data[data[:,1].<="2019-12-31",:] 

#### ii) Run on an outdated dataset (WHO, 2020)
# data = readdlm("Figs_Appendix/FigS5-S6_SensitivityAnalyses_Datasets/SetParams_Supplement/CaseData_Wuhan_WHO2020.csv",',',skipstart=1)


# Number of first N cases (symptoms onset)
N_cases = sum(data[:,2])

# Date of N-th case
Date_1 = Dates.Date(data[1,1])
Date_N = Dates.Date(data[end,1])

##
## 2. Context-specific parameters #################################
##

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
# R0 = 3.54                 # (Hao et al., Nature, 2020)
R0 = 2.5
kappa = 0.1                 # overdispersion parameter (Endo et al., Wellcome Open Research, 2020)
p_neg = kappa/(kappa+R0)    # success probability for the negative-binomial

## Stopping criteria
min_n_infect = 5*N_cases/p_detect   # = 102400 
max_t_infect = 365                  # Max. time (days) for epidemic developing
