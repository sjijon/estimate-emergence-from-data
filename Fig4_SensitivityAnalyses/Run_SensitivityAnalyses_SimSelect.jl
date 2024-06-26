## Jijón S, Czuppon P, Blanquart F & Débarre F (2023). 
## Using early detection data to estimate the date of emergence of an epidemic outbreak.
## https://github.com/sjijon/estimate-emergence-from-data
##
##
## Varying the parameter values on the mean simulated epidemic
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
Var_str = "Varying_tol_epi"

## Number of repetitions
# repeats = 5000
repeats = 100


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
VarParam = Var_tol_epi


for param in VarParam
    ## Set values
    global tol_epi = param

    ## Print params values
    println("\n\ntol_epi: $tol_epi")

    if SaveResults == "Yes"
        ## Filepaths
        global tol_epi_str = string(Int(10*tol_epi))
        global dir_output = string("Fig4_SensitivityAnalyses/Output/",EpiContext,"/",Var_str)
        mkpath(dir_output)

        ## Time to N cases
        global file_Cases_EpiSize_Time = string(dir_output,"/Cases_EpiSize_Time_tol_epi_0",tol_epi_str,".csv")
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
                ## i) The time period between the 1st infection and the last observed case
                local obs_num_days = Dates.value(Date_N - Date_1)
                local sim_num_days = length(SimCases.cumul[SimCases.cumul.>0])
                local Diff_SimInf1_ObsCas1 = Int(SimCases.d_detect[end] - obs_num_days)
                
                ## ii) The similarity with the observed cumulative number of cases
                ## Add zeros if&where needed and align the cumulative curves at right
                local obs_cases_cumul_ = [zeros(Int,maximum([length(SimCases.cumul),length(obs_cases_cumul)])-length(obs_cases_cumul));obs_cases_cumul]
                local sim_cases_cumul_ = [zeros(Int,maximum([length(SimCases.cumul),length(obs_cases_cumul)])-length(SimCases.cumul));SimCases.cumul]
                ## Compute the maximum of the absolute difference between the simulated and the observed cumulative cases
                local Dist_pw = abs.(sim_cases_cumul_.-obs_cases_cumul_)
                local Dist = maximum(Dist_pw)

                ## Condition
                VerifiesCondition = (Diff_SimInf1_ObsCas1>=0 && Dist<(tol_epi*N_cases))

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