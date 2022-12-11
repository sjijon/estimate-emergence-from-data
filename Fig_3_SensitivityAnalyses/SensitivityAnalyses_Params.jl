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
using LinearAlgebra, SparseArrays
using StatsBase, Distributions, Random
using Printf, DelimitedFiles, CSV
using Dates, TimerOutputs

## Create and start the timer
to = TimerOutput()

## Load user functions
include("../Processes/InfectionProcess.jl")
include("../Processes/DetectionProcess.jl")

## Defining the structures
include("../Processes/DefineStructs.jl")

## Display new run in the REPL
@printf "\n\n.................NEW RUN.................\n"

##
## 1a. User-determined parameters ########################
## (uncomment the line corresponding to the selection)
##
## Epidemiological context
##
# EpiContext = "Alpha_UK"
EpiContext = "COVID-19_Wuhan"

## Varying parameter 
# Var_str = "Varying_R0"
Var_str = "Varying_pdetect"

## Number of repetitions
repeats = 5000
# repeats = 1

## Select wether or not to save results
SaveResults = "Yes";
# SaveResults = "No";

##
## 1b. Read data and set parameters ########################
##
@timeit to "Parameters set" begin # timer
    ## Read parameters for chosen epi context
    include(string("../SetParameters/Params_",EpiContext,".jl"))
    println("\nEpidemiological context: $EpiContext")
    println("(N=$N_cases cases reported by $(Date(Date_N)))\n")

    ## Create object of observations
    ObsCases = ObsDetect(Date.(data[:,1]),data[:,2],N_cases,Date_N)

    ## Data
    obs_cases_cumul = Array{Int64}(cumsum(ObsCases.cases))               # Cumulative number of cases
    obs_cases_daynum = Dates.value.(ObsCases.dates- ObsCases.dates[1]) .+ 1   # Day number

    ## Time parameters (in days)
    dt = 1/10 # time step in the simulation (only<1 day alllowed)
end # timer (parameters)

##
## 2. Run simulations ##############################
##

@timeit to "Iterations" begin # TimerOutput
    println("\n", string(split(Var_str,"_")[1], " ", split(Var_str,"_")[2]))
    
    # Set varying values depending on EpiContext
    if EpiContext == "COVID-19_Wuhan"
        # Var_p_detect=[0.1,0.15,0.25]
        Var_p_detect=[0.05]
        Var_R0 = [1.5,2.0,2.5,3.5]
    elseif EpiContext == "Alpha_UK"
        Var_p_detect=[0.005,0.0105,0.025]
        Var_R0 = [1.5,1.9,2.5]
    end

    ## Run for different parameters
    if Var_str == "Varying_R0"
        VarParam = Var_R0
    elseif Var_str == "Varying_pdetect"
        VarParam = Var_p_detect
    end

    for param in VarParam
        if Var_str == "Varying_R0"
            global R0=param
            global p_neg = kappa/(kappa+R0)
        elseif Var_str == "Varying_pdetect"
            global p_detect = param
        end
        
        ## Print params values
        println("\n\np_detect: $p_detect")
        println("R0: $R0")

        ## Increase the minimum number of infections to reach N cases
        if EpiContext == "Alpha_UK" && p_detect==p_detect<0.01
            global min_n_infect = 20*N_cases/p_detect
        elseif EpiContext == "COVID-19_Wuhan" && p_detect<=0.1
            global max_t_infect = 365        
            global min_n_infect = 20*N_cases/p_detect
        end

        if SaveResults == "Yes"
            ## Filepaths
            global R0_str = string(Int(10*R0))
            global p_detect_str = string(Int(10000*p_detect))
            global dir_output = string("Time distribution for N cases/RunSims_to_Ncases/Output/",EpiContext,"/SensitivityAnalyses/",Var_str)
            mkpath(dir_output)

            ## Cumulative cases
            file_cumul_cases = string(dir_output,"/CumulCases_R0_",R0_str,"_p_detect_0",p_detect_str,".csv")            
            open(file_cumul_cases, "w")
            ## Time to N cases
            file_MinTime_N_EpiSize = string(dir_output,"/MinTime_N_EpiSize_R0_",R0_str,"_p_detect_0",p_detect_str,".csv")
        end

        ## Initialization
        global run_num = 0                               # repetition index
        local successes = 0                             # successes counter
        global MinTime_N_EpiSize = Array{Float64}(undef,0,4)   # time and size of the epidemic, and param

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
                
                    ## Time from the first infection and epidemic size at time where the N-th case is determined
                    global MinTime_N_EpiSize = [MinTime_N_EpiSize;[SimCases.d_NthCase SimCases.cumul[end] EpiSize_NthCaseInfect] SimCases.d_detect[1]]

                    ## Save
                    if SaveResults == "Yes"
                        ## All cumulative cases                            
                        open(file_cumul_cases, "a") do io
                        writedlm(io,SimCases.cumul',",")
                        end

                        ## Time from 1st infection to N-th case
                        writedlm(file_MinTime_N_EpiSize, MinTime_N_EpiSize, ',')
                    end
                    
                    ##
                    ## 2.3. Success ! (No condition to select simulations)
                    ##
                    successes += 1
                    
                    # Print progress
                    if successes == 1
                        println("\nProgress:\nSuccess No. $successes of $repeats")
                    elseif successes % (repeats/10) == 0
                        # println("Success No. $successes of $repeats")
                        # println("Success No. $successes of $repeats (epi_size: $EpiSize_NthCaseInfect, num_cases: $(SimCases.cumul[end]), day: $(SimCases.d_NthCase), run_num: $run_num)")
                        println("Success No. $successes of $repeats (epi_size: $EpiSize_NthCaseInfect, num_cases: $(SimCases.cumul[end]), day: $(SimCases.d_NthCase))")
                    end
                end # if (epidemic)
            end # if (Cases >= N_cases)
        end # while (iterations)
    end # for (p_detect or R0) 
end # timer (iterations)

#
# 5. Display results and save ####################
#


@printf "\n\nComputation time:\n"
show(to; allocations = false)