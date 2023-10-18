## Estimating the time distribution between the first infection and N samples
## Jijon S, Czuppon P, Blanqart F., Débarre F
## iEES, 2022
##
## Varying th parameter values on the mean simulated epidemic
## Under no condition
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

## Defining the structures
include("../Routines/DefineStructs.jl")

## Display new run in the REPL
println("\n\n.................NEW RUN.................\n")

##
## 1a. User-determined parameters ########################
## (uncomment the line corresponding to the selection)
##
## Epidemiological context
##
EpiContext = "Alpha_UK"
# EpiContext = "COVID-19_Wuhan"

## Varying parameter 
Var_str = "Varying_R0"
# Var_str = "Varying_kappa"
# Var_str = "Varying_pdetect"
# Var_str = "Varying_theta_tau"

## Number of repetitions
repeats = 5000
# repeats = 100

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
tol_epi = 0.30

## Create object of observations
ObsCases = ObsDetect(Date.(data[:,1]),data[:,2],N_cases,Date_N)

## Data
obs_cases_cumul = Array{Int64}(cumsum(ObsCases.cases))               # Cumulative number of cases
obs_cases_daynum = Dates.value.(ObsCases.dates- ObsCases.dates[1]) .+ 1   # Day number

## Time parameters (in days)
dt = 1/10 # time step in the simulation (only<1 day alllowed)

##
## 2. Run simulations ##############################
##
println("\n", string(split(Var_str,"_")[1], " ", split(Var_str,"_")[2]))

# Set varying values depending on EpiContext
if EpiContext == "COVID-19_Wuhan"
    Var_p_detect=[0.1,0.15,0.25]        # Baseline = 0.15
    Var_R0 = [1.5,2.0,2.5,3.5]          # Baseline = 2.5
    Var_kappa = [0.05,0.1,0.25]         # Baseline = 0.1
    Var_theta_tau = [2.8,6.25,11.6]     # Baseline = 6.25
elseif EpiContext == "Alpha_UK"
    Var_p_detect=[0.005,0.0105,0.025]   # Baseline = 0.0105
    Var_R0 = [1.7,1.9,2.1]              # Baseline = 1.9
    Var_kappa = [0.35,0.57,0.75]        # Baseline = 0.57
    Var_theta_tau = [7,12,15]           # Baseline = 12
end

## Run for different parameters
if Var_str == "Varying_R0"
    VarParam = Var_R0
elseif Var_str == "Varying_kappa"
    VarParam = Var_kappa
elseif Var_str == "Varying_pdetect"
    VarParam = Var_p_detect
elseif Var_str == "Varying_theta_tau"
    VarParam = Var_p_detect
end

for param in VarParam
    if Var_str == "Varying_R0"
        global R0=param
        global p_neg = kappa/(kappa+R0)
    elseif Var_str == "Varying_kappa"
        global kappa=param
        global p_neg = kappa/(kappa+R0)
    elseif Var_str == "Varying_pdetect"
        global p_detect = param
    elseif Var_str == "Varying_theta_tau"
        global shape_detect = param
        global SampTimeDist = Gamma(shape_detect, scale_detect)
    end
    
    ## Print params values
    println("\n\np_detect: $p_detect")
    println("R0: $R0")
    println("kappa: $kappa")
    println("theta_tau: $shape_detect")

    ## Increase the minimum number of infections to reach N cases
    if EpiContext == "Alpha_UK" && p_detect<0.01
        global min_n_infect = 20*N_cases/p_detect
    elseif EpiContext == "Alpha_UK" && R0<2
        global min_n_infect = 50*N_cases/p_detect
    elseif EpiContext == "COVID-19_Wuhan" && p_detect<=0.1
        global max_t_infect = 365        
        global min_n_infect = 20*N_cases/p_detect
    end

    if SaveResults == "Yes"
        ## Filepaths
        global R0_str = string(Int(10*R0))
        global p_detect_str = string(Int(10000*p_detect))
        global kappa_str = string(Int(round(100*kappa)))
        global theta_str = string(Int(round(100*shape_detect)))

        global dir_output = string("Fig4_SensitivityAnalyses/Output/",EpiContext,"/",Var_str)
        mkpath(dir_output)

        ## Cumulative cases
        global file_cumul_cases = string(dir_output,"/CumulCases_R0_",R0_str,"_kappa",kappa_str,"_p_detect_0",p_detect_str,"_theta_tau",theta_tau_str,".csv")            
        open(file_cumul_cases, "w")
        ## Time to N cases
        global file_Cases_EpiSize_Time = string(dir_output,"/Cases_EpiSize_Time_R0_",R0_str,"_kappa_0",kappa_str,"_p_detect_0",p_detect_str,"_theta_tau",theta_tau_str,".csv")
    end

    ## Initialization
    global run_num = 0                                      # repetition index
    local successes = 0                                     # successes counter
    global Cases_EpiSize_Time = Array{Float64}(undef,0,4)    # time and size of the epidemic, and param

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
        run_num += 1 # Increase run index

        Random.seed!(run_num)             # setting the seed
        
        global SimEpi = InfectionProcess()
        
        if SimEpi.epidemic==true # Epidemic? Yes (true) or No (false)
            ##
            ## 2.1.2 Detection model
            ##
            global SimCases = DetectionProcess(SimEpi)

            # println("num_cases: $(sum(SimEpi.infections)), num_cases: $(SimCases.cumul[end])")

            if SimCases.num>N_cases
                ## Total epidemic size at the time of infection of the N-th case
                global EpiSize_NthCaseInfect = SimEpi.daily_infec[1:Int(SimCases.d_detect[N_cases])] |> sum

                ##
                ## 2.2. Selecting the simulations
                ##
                ## i) The time period between the 1st infection and the first observed case
                global obs_num_days = Dates.value(Date_N - Date_1)
                global sim_num_days = length(SimCases.cumul[SimCases.cumul.>0])
                global Diff_SimInf1_ObsCas1 = Int(SimCases.d_detect[end] - obs_num_days)
                
                ## ii) The similarity with the observed cumulative number of cases
                ## Add zeros if&where needed and align the cumulative curves at right
                global obs_cases_cumul_ = [zeros(Int,maximum([length(SimCases.cumul),length(obs_cases_cumul)])-length(obs_cases_cumul));obs_cases_cumul]
                global sim_cases_cumul_ = [zeros(Int,maximum([length(SimCases.cumul),length(obs_cases_cumul)])-length(SimCases.cumul));SimCases.cumul]
                ## Compute the maximum of the absolute difference between the simulated and the observed cumulative cases
                global Dist_pw = abs.(sim_cases_cumul_.-obs_cases_cumul_)
                global Dist = maximum(Dist_pw)

                ## Select the simulation if the selected condition is verified
                if (Diff_SimInf1_ObsCas1>=0 && Dist<(tol_epi*N_cases))
                
                    # Update number of successes
                    successes += 1
            
                    ## Time from the first infection and epidemic size at time where the N-th case is determined
                    global Cases_EpiSize_Time = [Cases_EpiSize_Time;[SimCases.d_NthCase SimCases.cumul[end] EpiSize_NthCaseInfect] SimCases.d_detect[1]]

                    ## Save
                    if SaveResults == "Yes"
                        ## All cumulative cases                            
                        open(file_cumul_cases, "a") do io
                        writedlm(io,SimCases.cumul',",")
                        end

                        ## Time from 1st infection to N-th case
                        writedlm(file_Cases_EpiSize_Time, Cases_EpiSize_Time, ',')
                    end
                
                    # Print progress
                    if successes == 1
                        println("\nProgress:\nSuccess No. $successes of $repeats")
                    elseif successes % (repeats/10) == 0
                        # println("Success No. $successes of $repeats")
                        # println("Success No. $successes of $repeats (epi_size: $EpiSize_NthCaseInfect, num_cases: $(SimCases.cumul[end]), day: $(SimCases.d_NthCase), run_num: $run_num)")
                        println("Success No. $successes of $repeats (epi_size: $EpiSize_NthCaseInfect, num_cases: $(SimCases.cumul[end]), day: $(SimCases.d_NthCase))")
                    end # if (progress)
                end # if (sim selected)
            end # if (epidemic)
        end # if (Cases >= N_cases)
    end # while (iterations)
end # for (p_detect or R0) 