## Estimating the time distribution between the first infection and N cases
##
## Jijon S, Czuppon P, Blanqart F., Débarre F
## iEES, 2022/2023
##
#################################################

##
## 0. Setup #####################################
##
## Load packages
using DelimitedFiles
using StatsBase, Distributions, Random
using Dates

## Load user functions and other setup scripts
include("../Routines/DefineStructs.jl")

## Display new run in the REPL
println("\n\n.................ABC Figure.................\n")

##
## 1a. User-determined parameters ########################d
##
## Epidemiological context
##
EpiContext = "COVID-19_Wuhan"

## Select wether or not to save results
SaveResults = "Yes";
# SaveResults = "No";

##
## 1b. Read data and set parameters ########################
##

## Read parameters for chosen epi context
include(string("../SetParameters_EpiContext/Params_",EpiContext,".jl"))
println("\nEpidemiological context: $EpiContext")
println("(N=$N_cases cases reported by $(Date(Date_N)))\n")

## Set tolerances
# abs(Difference in daily cases) <= tol_epi
tol_epi = 0.3   

## Create object of observations
ObsCases = ObsDetect(Date.(data[:,1]),data[:,2],N_cases,Date_N)

## Data
obs_cases_cumul = Array{Int64}(cumsum(ObsCases.cases)) # Cumulative number of cases
obs_cases_daynum = Dates.value.(ObsCases.dates- ObsCases.dates[1]) .+ 1 # Day number

##
## 2. Define function to evaluate condition (C4) ##########
##
function maxDiff(obs,sim)
    local Dist_pw = abs.(sim .- obs)
    if (maximum(Dist_pw)<(tol_epi*N_cases)) 
        return(1)
    else 
        return(0)
    end
end

##
## 3a. Read in and evaluate ABC simulations  ##########
##
## Number of days evaluated
repeatDays = 101
startDay = 40
simreps = 9999
keepTraj = zeros(repeatDays)

println("\n\n.................Reading in data.................\n")
for i in 1:repeatDays
    println(string("Day_", startDay+i,"\n"))
    for j in 0:simreps
        # read in the simulation
        a = readdlm(string("../ABC_Simulations/Days_",startDay+i,"/ABCSim_no_",j,".txt"),',')
        local cumvec_sim = cumsum(a[:,3])       # cumulative sum of detections
        # add zeros to observations and align at the last date (=Nth detection)
        local obs_cases_cumul_adj = [zeros(Int,length(cumvec_sim)-length(obs_cases_cumul));obs_cases_cumul]
        # keep simulated trajectory if condition (C4) is satisfied 
        if (maxDiff(obs_cases_cumul_adj,cumvec_sim) == 1)
            keepTraj[i] += 1/(simreps+1)
        end
    end
end

# store result in txt file
writedlm("../FigABC/ABCresults.csv", keepTraj)