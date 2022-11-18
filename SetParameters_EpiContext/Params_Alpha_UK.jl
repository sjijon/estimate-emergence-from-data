## Parametrizing the model
## 
## Infectious disease: Alpha variant
## Location: UK
## Period: Sept. 20-Nov. 30, 2021
## Data: date of sequenced samples of the B.1.1.7 variant
##       cases submitted up to Nov 30th (Gisaid)
##       (N=457)
#################################################

##
## 1. Read data #################################
##

## Read data on Alpha cases in the UK, where a case is defined as a sample:
##   i) collected in the UK between Sept 20th and Nov. 30th, 2020, 
##  ii) sequenced at unknown date, determining infection with the Alpha variant,
## iii) submitted to Gisaid by November 30th, 2020.

## Read data file
data = readdlm("Data/CaseData_Alpha_UK.csv",',',skipstart=1)
data = data[data[:,1].<="2020-11-11",:] ## Use data of sample collected up to Nov 11

## Date and number of Confirmed B.1.1.7 cases in the UK
Date_1 = Dates.Date(data[1,1])
Date_N = Dates.Date(data[end,1])
N_cases = sum(data[:,2])

##
## 2. Context-specific parameters #################################
##
## All times are given in days

##
## 2.1. Detection (sampling and sequencing)
##
## Time from infection to sampling (sample collection)
## ~ Gamma distribution, mean of 7 days
shape_detect = 12
scale_detect = 7/12

## Probability of detection 
p_detect = 0.0105           # = p_sampling * p_sequencing =  0.25 * 0.042

##
## 2.2. Offspring generation
##
## Time to secondary infection
## ~ Gamma distribution (mean = 5.5 days)
shape_inf = 6.6
scale_inf = 0.833

## Number of secondary infections 
## ~ NegativeBinomial(kappa=successes,p_neg=success rate)
kappa = 0.57                # dispersion parameter (Salje et al., Science, 2020)
R0 = 1.9      # Hill et al.
p_neg = kappa/(kappa+R0)    # success probability for the negative-binomial      

## Stopping criteria
min_n_infect = 10*N_cases/p_detect      ## Min. number of infections
max_t_infect = 365      ## Max. time (days) for epidemic developing
