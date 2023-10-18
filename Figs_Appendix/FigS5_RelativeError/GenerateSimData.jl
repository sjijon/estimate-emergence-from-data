## Estimating the time distribution between the first infection and N cases
##
## Script to generate synthetic data
##
## Jijon S, Czuppon P, Blanqart F., Débarre F
## iEES, 2023
##
#################################################

##
## 0. Setup #####################################
##
## Load packages
using DelimitedFiles
using StatsBase, Distributions, Random
using Dates
using CSV, DataFrames

## Create and start the timer
to = TimerOutput()

## Load user functions
include("../../Routines/InfectionProcess.jl")
include("../../Routines/DetectionProcess.jl")

## Defining the structures
include("../../Routines/DefineStructs.jl")

## Display new run in the REPL
println("\n\n.................NEW RUN.................\n")

##
## 1a. User-determined parameters ########################d
## (uncomment the line corresponding to the selection)
##
## Epidemiological context
##
# EpiContext = "Alpha_UK"
EpiContext = "COVID-19_Wuhan"

## Number of generated datasets
num_sims = 1000

## Select wether or not to save results
SaveResults = "Yes";
# SaveResults = "No";

##
## 1b. Read data and set parameters ########################
##
## Read parameters for chosen epi context
include(string("../../SetParameters_EpiContext/Params_",EpiContext,".jl"))

## Set number of cases
N_cases = 1000

println("\nGenerating synthetic data\n\nEpidemiological context: $EpiContext")
println("(N=$N_cases cases reported by $(Date(Date_N)))\n")

## Time parameters (in days)
dt = 1/10 # time step in the simulation (only<1 day alllowed)
length_inf_vec = trunc(Int,max_t_infect/dt) # final vector length

##
## 1c. Define filepaths ########################
##
## Create path
dir_output = string("Data/SimulatedData/",EpiContext,"/sims")
mkpath(dir_output)

##
## 2. Run simulations ################################
##

## Initialization of counters
run_num = 0                                 # repetition index
successes = 0                               # successes counter

## Initialization of results arrays
# SIMS = [1:max_t_infect;]                    # All simulations
SecInf_stats = Array{Int64}(undef,0,6)      # All number of secondary infections (stats for each SimEpi)

###
### Distributions
###
## Time of secondary infection ~ Gamma(α=shape,θ=scale)
global InfTimeDist = Gamma(shape_inf,scale_inf)
## Number of secondary infections ~ NegativeBinomial(r=successes,p=success rate)
global TransDist = NegativeBinomial(kappa,p_neg)
## Time of detection ~ Gamma(α=shape,θ=scale)
global SampTimeDist = Gamma(shape_detect, scale_detect)

## Run simulations until reaching "num_sims" epidemics
while (successes < num_sims)
    ## Run stochastic simulations
    global run_num += 1 # Increase run index

    ##
    ## 2.1. Run transmission model 
    ##
    Random.seed!(run_num)             # setting the seed
    
    global SimEpi = InfectionProcess()

    ##
    ## 2.2 Run detection model detection on successful epidemics
    ##
    if SimEpi.epidemic == true # Epidemic? Yes (true) or No (false)
        # Update number of successes
        global successes += 1

        ##
        ## Attribute dates
        ##
        ## Random date of first infection
        local Date_1stInf = rand(Date("2019-09-01"):Day(1):Date("2019-12-10"))

        ## Get cases from simulated cumulative curve
        local num_days = length(SimEpi.daily_infec)+1
        global SimData_Inf = Array{Any}(undef,num_days,2)
        ## header
        SimData_Inf[1,1] = "Date"
        SimData_Inf[1,2] = "Infections"
        ## 1st row
        SimData_Inf[2,1] = Date_1stInf
        for i in 3:num_days
            ## Date
            SimData_Inf[i,1] = SimData_Inf[i-1,1] + Day(1)
        end
        SimData_Inf[2:end,2] = SimEpi.daily_infec

        ##
        ## Detected cases
        ##
        global SimCases = DetectionProcess(SimEpi)

        ##
        ## Attribute dates to cases
        ##
        local Dates_SimCases = Array{Any}(undef,length(SimCases.d_detect),1)
        for i in 1:length(SimCases.d_detect)
            Dates_SimCases[i] = Date_1stInf + Day(SimCases.d_detect[i]- 1)
        end

        
        local Dates_SimCases_list = unique(Dates_SimCases)
        local NumDet = Array{Int64}(undef,length(Dates_SimCases_list),1)
        for i in 1:length(NumDet)
            NumDet[i] = count(x->x==Dates_SimCases_list[i],Dates_SimCases)
        end
        global SimData_Det = [Dates_SimCases_list NumDet]

        ##
        ## Print progress
        ##
        if successes == 1
            println("\nPROGRESS\n")
            println("set_num: $successes, run_num: $run_num, infections: $(sum(SimEpi.infections)), cases: $(SimCases.cumul[end])")
        elseif successes % (num_sims/10) == 0 || successes==num_sims
            println("set_num: $successes, run_num: $run_num, infections: $(sum(SimEpi.infections)), cases: $(SimCases.cumul[end])")
        end

        ##
        ## 2.3. Save results
        ##
        local file_sim_data_inf = string(dir_output,"/simdata_inf_sim",successes,".csv")
        open(file_sim_data_inf, "w") do io
            writedlm(io, SimData_Inf)
        end

        local file_sim_data_det = string(dir_output,"/simdata_det_sim",successes,".csv")
        open(file_sim_data_det, "w") do io
            writedlm(io, SimData_Det)
        end                
    end # if (epidemic)
end # while (iterations)