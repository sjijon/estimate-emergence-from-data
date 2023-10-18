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
using TimerOutputs
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
## 1a. User-determined parameters ########################
EpiContext = "SimulatedData"

# global num_datasets = 500 
# global num_datasets = 100 # (Results in the manuscript)
# global num_datasets = 10
global num_datasets = 2
# global num_datasets = 1

## Number of repetitions per dataset
# repeats = 1000
# repeats = 100 # (Results in the manuscript)
repeats = 10
# repeats = 2
# repeats = 1

N_cases_aux = 1
# N_cases_aux = 10
# N_cases_aux = 100
# N_cases_aux = 1000 # (Results in the manuscript)

## Select wether or not to save results
SaveResults = "Yes";
# SaveResults = "No";

##
## 1b. Define filepaths ########################
##
if SaveResults == "Yes"
    ## Create path
    dir_output = string("Time distribution for N cases/RunSims_to_Ncases/Output/",EpiContext)
    mkpath(dir_output)

    ## Time and epi size
    file_MinTime_N_EpiSize_All = string(dir_output,"/MinTime_N_EpiSize_Case",N_cases_aux,"_",repeats,"sims_",num_datasets,"datasets.csv")
end

##
## 2. Run simulations for each dataset ########################
##
println("\nEpidemiological context: $EpiContext (first ~$N_cases_aux cases)")

println("\n\n\tSIMULATIONS")
println("______________________________________________\n")
# dataset_num = 1;

MinTime_N_EpiSize_All = Array{Any}(undef,0,5) # Stats per dataset

@timeit to "Iterations" begin # timer
    global SelectedSims = Array{Int64}(undef,0,3)

    for dataset_num in 1:num_datasets

        ##
        ## 2a. Read data and set parameters ########################
        ##
        ## Read parameters for chosen epi context
        global dataset_num_ = dataset_num
        include(string("../../SetParameters_EpiContext/Params_",EpiContext,".jl"))

        ## Print progress
        if dataset_num == 1
            println("dataset_num: $dataset_num (N=$N_cases in $(length(data)) days)")
        elseif dataset_num % 50 == 0 || dataset_num==num_datasets
            println("dataset_num: $dataset_num (N=$N_cases in $(length(data)) days)")
            # println("dataset_num: $dataset_num (N=$N_cases reported by $(Date(Date_N)))")
        end

        ## Set tolerances
        # length(Simulated cases) >= tol_delay*length(Observed cases)
        global tol_delay = 0.90
        # abs(Difference in daily cases) <= tol_epi
        global tol_epi = 0.30
        
        ## Create object of observations
        global ObsCases = ObsDetect(Date.(data[:,1]),data[:,2],N_cases,Date_N)

        ## Data
        global obs_cases_cumul = Array{Int64}(cumsum(ObsCases.cases))               # Cumulative number of cases
        global obs_cases_daynum = Dates.value.(ObsCases.dates- ObsCases.dates[1]) .+ 1   # Day number

        ## Time parameters (in days)
        global dt = 1/10 # time step in the simulation (only<1 day alllowed)
        global length_inf_vec = trunc(Int,max_t_infect/dt) # final vector length

        ##
        ## 2b. Run simulations ################################
        ##

        ## Initialization of counters
        local run_num = 0                                 # repetition index
        local successes = 0                               # successes counter

        ## Initialization of results arrays
        global MinTime_N_EpiSize = Array{Int64}(undef,0,4) # Time of detection, number of cases and size of the epidemic
        ###
        ### Distributions
        ###
        ## Time of secondary infection ~ Gamma(α=shape,θ=scale)
        global InfTimeDist = Gamma(shape_inf,scale_inf)
        ## Number of secondary infections ~ NegativeBinomial(r=successes,p=success rate)
        global TransDist = NegativeBinomial(kappa,p_neg)
        ## Time of detection ~ Gamma(α=shape,θ=scale)
        global SampTimeDist = Gamma(shape_detect, scale_detect)

        ## Run simulations until reaching "repeats" epidemics
        while (successes < repeats)
            ## Run stochastic simulations
            run_num += 1 # Increase run index
            # println("run_num: ",run_num)

            try
                ##
                ## 2.1. Run the model
                ##                
                ##
                ## 2.1.1 Transmission model 
                ##
                global SimEpi = InfectionProcess()

                ##
                ## Keep simulation
                ##
                if SimEpi.epidemic == true # Epidemic? Yes (true) or No (false)

                    ##
                    ## 2.1.2 Detection model
                    ##
                    global SimCases = DetectionProcess(SimEpi)

                        ## Continue if SimCases >= N_cases
                        if SimCases.cumul[end]>=N_cases
                            ## Total epidemic size at the time of infection of the N-th case
                            global EpiSize_NthCaseInfect = SimEpi.daily_infec[1:Int(SimCases.d_detect[N_cases])] |> sum

                            ##
                            ## 2.2. Selecting the simulations
                            ##
                            ## i) The time period between the 1st infection and the first observed case
                            global obs_num_days = Dates.value(Date_N - Date_1)
                            global dataset_num_days = length(SimCases.cumul[SimCases.cumul.>0])
                            
                            global Diff_SimInf1_ObsCas1 = Int(SimCases.d_detect[end] - obs_num_days)
                            ## ii) The length of the time period where cases occur
                            global Diff_NumDays = abs(obs_num_days - dataset_num_days)
                                
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
                                
                                ##
                                ## 2.4. Build results arrays
                                ##
                                ## Time from the first infection and epidemic size at time where the N-th case is determined
                                global MinTime_N_EpiSize = [MinTime_N_EpiSize; SimCases.d_NthCase SimCases.cumul[end] EpiSize_NthCaseInfect SimCases.d_detect[1]]
                            end # if (success)
                        end # if (cases >= N_cases)
                end # if (epidemic)
            catch
                println("\t(some error ocurred, excluding run_num=$run_num and running next)")
            end
        end # while (iterations)
        
        ##
        ## 4. Buid all results array
        ##
        global MinTime_N_EpiSize_All = [MinTime_N_EpiSize_All;repeat([string("Sim ",dataset_num)],size(MinTime_N_EpiSize,1)) MinTime_N_EpiSize]

        global SelectedSims = [SelectedSims; [repeats run_num round(100*repeats/run_num,digits=2)]]
    end # for (datasets)
end # timer (iterations)

println("______________________________________________\n")

##
## 5. Display results and save ####################
##
@timeit to "Display and save" begin # timer
    ##
    ## Display
    ##
    print_results("Retained simulations", SelectedSims[:,3]) 
    println("______________________________________________\n")

    ##
    ## Save
    ##
    if SaveResults == "Yes"
        ## Save
        writedlm(file_MinTime_N_EpiSize_All, [["Dataset" "MinTime" "Cases" "EpiSize" "Case1"]; MinTime_N_EpiSize_All], ',')
    end
end     