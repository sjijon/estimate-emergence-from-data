## Estimating the time distribution between the first infection and N cases
##
## Jijon S, Czuppon P, Blanqart F., Débarre F
## iEES, 2022
##
#################################################

##
## 0. Setup #####################################
##
## Load packages
using DelimitedFiles
using StatsBase, Distributions, Random
using Dates

## Load user functions
include("../Routines/InfectionProcess.jl")
include("../Routines/DetectionProcess.jl")
include("../Routines/DefineStructs.jl")

## Display new run in the REPL
println("\n\n.................NEW RUN.................\n")

##
## 1a. User-determined parameters ########################d
## (uncomment the line corresponding to the selection)
##
## Epidemiological context
##
EpiContext = "Alpha_UK"
# EpiContext = "COVID-19_Wuhan"

## Varying parameter 
# Var_str = "Varying_tol_epi"
Var_str = "Varying_tol_delay"

## Number of repetitions
repeats = 5000


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
# length(Simulated cases) >= tol_delay.8*length(Observed cases)
global tol_delay = 0.9
# abs(Difference in daily cases) <= tol_epi
global tol_epi = 0.3   

## Create object of observations
ObsCases = ObsDetect(Date.(data[:,1]),data[:,2],N_cases,Date_N)

## Data
obs_cases_cumul = Array{Int64}(cumsum(ObsCases.cases))               # Cumulative number of cases
obs_cases_daynum = Dates.value.(ObsCases.dates- ObsCases.dates[1]) .+ 1   # Day number

## Time parameters (in days)
dt = 1/10 # time step in the simulation (only<1 day alllowed)
length_inf_vec = trunc(Int,max_t_infect/dt) # final vector length

##
## 2. Run simulations ################################
##
println("\n", string(split(Var_str,"_")[1], " ", split(Var_str,"_")[2], "_", split(Var_str,"_")[3]))

# Varying values of tolerances
if EpiContext == "Alpha_UK"
    Var_tol_epi = [0.3 0.5]     # Baseline: 0.3
elseif EpiContext == "COVID-19_Wuhan"
    Var_tol_epi = [0.1 0.3 0.5]     # Baseline: 0.3
end
Var_tol_delay = [0.8 0.9 1.0]   # Baseline: 0.9

## Run for different parameters
if Var_str == "Varying_tol_epi"
    VarParam = Var_tol_epi
elseif Var_str == "Varying_tol_delay"
    VarParam = Var_tol_delay
end

for param in VarParam
    ## Set values
    if Var_str == "Varying_tol_epi"
        global tol_epi = param
        tol_delay = 0.9
    elseif Var_str == "Varying_tol_delay"
        tol_epi = 0.3
        global tol_delay = param
    end

    ## Print params values
    println("\n\ntol_epi: $tol_epi")
    println("tol_delay: $tol_delay")

    if SaveResults == "Yes"
        ## Filepaths
        global tol_epi_str = string(Int(10*tol_epi))
        global tol_delay_str = string(Int(10*tol_delay))
        global dir_output = string("Fig3_SensitivityAnalyses/Output/",EpiContext,"/",Var_str)
        mkpath(dir_output)

        ## Time to N cases
        global file_Cases_EpiSize_Time = string(dir_output,"/Cases_EpiSize_Time_tol_epi_0",tol_epi_str,"_tol_delay_0",tol_delay_str,".csv")
    end
    
    ## Initialization
    global run_num = 0                                 # repetition index
    local successes = 0                                # successes counter
    global Cases_EpiSize_Time = Array{Int64}(undef,0,4) # Time of detection, number of cases and size of the epidemic
    global SIMS = [1:max_t_infect;]                    # All simulations
    
    ###
    ### Distributions
    ###
    ## Time of secondary infection ~ Gamma(α=shape,θ=scale)
    global InfTimeDist = Gamma(shape_inf,scale_inf)
    ## Number of secondary infections ~ NegativeBinomial(r=successes,p=success rate)
    global TransDist = NegativeBinomial(kappa,p_neg)
    ## Time of detection ~ Gamma(α=shape,θ=scale)
    global SampTimeDist = Gamma(shape_detect, scale_detect)

    while (successes < repeats)
        ## Run stochastic simulations
        global run_num += 1 # Increase run index
        
        Random.seed!(run_num)             # setting the seed
        
        global SimEpi = InfectionProcess()
        
        if SimEpi.epidemic == true # Epidemic? Yes (true) or No (false)
            ##
            ## 2.1.2 Detection model
            ##
            global SimCases = DetectionProcess(SimEpi)

            # println("infections: $(sum(SimEpi.infections)), cases: $(SimCases.cumul[end])")

            if SimCases.cumul[end]>=N_cases
                ##
                ## 2.2. Selecting the simulations
                ##
                ## Selection conditioned to the cumulative number of cases & the length of time period (Con_Cumul_Delay)
                ## b) Time period
                ## b.1.) The time period between the 1st infection and the last observed case
                local obs_num_days = Dates.value(Date_N - Date_1)
                local sim_num_days = length(SimCases.cumul[SimCases.cumul.>0])
                local Diff_SimInf1_ObsCas1 = Int(SimCases.d_detect[end] - obs_num_days)
                ## b.2.) The length of the time period where cases occur
                local Diff_NumDays = abs(obs_num_days - sim_num_days)

                ## c) The similarity with the observed cumulative number of cases
                ## Add zeros if&where needed and align the cumulative curves at right
                local obs_cases_cumul_ = [zeros(Int,maximum([length(SimCases.cumul),length(obs_cases_cumul)])-length(obs_cases_cumul));obs_cases_cumul]
                local sim_cases_cumul_ = [zeros(Int,maximum([length(SimCases.cumul),length(obs_cases_cumul)])-length(SimCases.cumul));SimCases.cumul]
                ## Compute the maximum of the absolute difference between the simulated and the observed cumulative cases
                local Dist_pw = abs.(sim_cases_cumul_.-obs_cases_cumul_)
                local Dist = maximum(Dist_pw)

                ## Condition
                VerifiesCondition = (Diff_SimInf1_ObsCas1>=0 &&  sim_num_days>= tol_delay*obs_num_days && Dist<(tol_epi*N_cases))

                ## Select the simulation if the selected condition is verified
                if VerifiesCondition
                    ## Total epidemic size at the time of infection of the N-th case
                    global EpiSize_NthCaseInfect = SimEpi.daily_infec[1:Int(SimCases.d_detect[N_cases])] |> sum
                
                    ## Time from the first infection and epidemic size at time where the N-th case is determined
                    global Cases_EpiSize_Time = [Cases_EpiSize_Time;[SimCases.d_NthCase SimCases.cumul[end] EpiSize_NthCaseInfect] SimCases.d_detect[1]]
                    # println("Cases_EpiSize_Time: $Cases_EpiSize_Time")

                    ## Save
                    if SaveResults == "Yes"
                        writedlm(file_Cases_EpiSize_Time, Cases_EpiSize_Time, ',')
                    end
                    
                    ##
                    ## 2.3. Success ! (No condition to select simulations)
                    ##
                    successes += 1

                    # Print progress
                    if successes == 1
                        println("\nPROGRESS\n")
                        println("successes No. $successes of $repeats \t(epi_size: $EpiSize_NthCaseInfect, num_cases: $(SimCases.cumul[end]), run_num: $run_num)")
                    elseif successes % (repeats/10) == 0 || successes==repeats
                        println("successes No. $successes of $repeats \t(epi_size: $EpiSize_NthCaseInfect, num_cases: $(SimCases.cumul[end]), run_num: $run_num)")
                    end
                end # if (comparison)
            end # if (cases >= N_cases)
        end # if (epidemic)
    end # while (iterations)
end # for (tol_epi or tol_delay) 