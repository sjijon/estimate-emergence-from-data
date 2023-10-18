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

## Create and start the timer
to = TimerOutput()

## Load user functions
include("../Time distribution for N cases/RunSims_to_Ncases/InfectionProcess.jl")
include("../Time distribution for N cases/RunSims_to_Ncases/DetectionProcess.jl")

## Defining the structures
include("../Time distribution for N cases/DefineStructs.jl")

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

## Number of repetitions
repeats = 5000 # number of repetitions
# repeats = 1000
# repeats = 500
# repeats = 10
# repeats = 2
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

    ## Time parameters (in days)
    dt = 1/10 # time step in the simulation (only<1 day alllowed)
    length_inf_vec = trunc(Int,max_t_infect/dt) # final vector length
end # timer (parameters)

##
## 1c. Define filepaths ########################
##
if SaveResults == "Yes"
    ## Create path
    dir_output = string("Time distribution for N cases/RunSims_to_Ncases/Output/",EpiContext,"/",Cond_Name)
    mkpath(dir_output)

    ## Per-day cumulative number of cases (per SimEpi)
    SimData = string(dir_output,"/SimData",N_cases,"cases_",repeats,"sims.csv")

    ## Open file                
    open(file_cumul_cases, "w")
    open(file_cumul_cases, "a") do io
    writedlm(io,[1:1:4*30],",")
    end
    open(file_dist_pw, "w")
end

##
## 2. Run simulations ################################
##

## Initialization of counters
run_num = 0                                 # repetition index
successes = 0                               # successes counter

## Initialization of results arrays
SIMS = [1:max_t_infect;]                    # All simulations
SecondaryInfec_all = Array{Int64}(undef,0,2)# All number of secondary infections
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


## Run simulations until reaching "repeats" epidemics
@timeit to "Iterations" begin # timer
    while (successes < repeats)
        ## Run stochastic simulations
        global run_num += 1 # Increase run index
        ## Print run number
        # if run_num == 1
        #     println("\nRun number (run_num)\nsim_num: $run_num")
        # elseif run_num % 10000 == 0
            # println("run_num: $run_num")
        # end

        ##
        ## 2.1. Run the model
        ##
        Random.seed!(run_num)             # setting the seed
        
        ##
        ## 2.1.1 Transmission model 
        ##
        global SimEpi = InfectionProcess()
        
        ## Print progress
        # println("run_num: $run_num, infections: $(sum(SimEpi.infections))")
    
        ##
        ## Keep simulation
        ##
        if SimEpi.epidemic == true # Epidemic? Yes (true) or No (false)
            # println("\nEpidemic reached !") 

            ##
            ## 2.1.2 Detection model
            ##
            global SimCases = DetectionProcess(SimEpi)

            # println("infections: $(sum(SimEpi.infections)), cases: $(SimCases.cumul[end])")

            if epidemic
                # Update number of successes
                global successes += 1
                
                ##
                ## 2.4. Build results arrays
                ##
                ## Save all simulated epidemics
                global SIMS = hcat(SIMS,SimEpi.daily_infec)
                
                ## Save all number of secondary infections 
                global SecondaryInfec_all = vcat(SecondaryInfec_all,SecInf_dNthCase)
                
                ## Save
                if SaveResults == "Yes"
                    ## All cumulative cases
                    open(file_cumul_cases, "a") do io
                        writedlm(io,SimCases.cumul',",")
                    end

                    if Cond_Name == "Cond_Cumul_Delay"
                        # All distances (pointwise) to observations
                        open(file_dist_pw, "a") do io
                            writedlm(io,Dist_pw',",")
                        end
                    end
                end

                # Print progress
                if successes == 1
                    println("\nPROGRESS\n")
                    println("successes No. $successes of $repeats \t(epi_size: $EpiSize_NthCaseInfect, num_cases: $(SimCases.cumul[end]), run_num: $run_num)")
                elseif successes % (repeats/10) == 0 || successes==repeats
                    println("successes No. $successes of $repeats \t(epi_size: $EpiSize_NthCaseInfect, num_cases: $(SimCases.cumul[end]), run_num: $run_num)")
                end # if (progress)
            end # if (cases >= N_cases)
        end # if (epidemic)
    end # while (iterations)
end # timer (iterations)

Prop_SelectedSims = round(100*repeats/run_num,digits=2)
println("\n\nNumber of runs (run_num): $run_num\nSuccesses = $Prop_SelectedSims % of runs")

##
## 5. Display results and save ####################
##

@timeit to "Display and save" begin # timer
    if SaveResults == "Yes"
        ## Save
        writedlm(file_SecondaryInfec_all, SecondaryInfec_all, ',')
        writedlm(file_SecInf_stats, [["Mean" "Q025" "Q25" "Q50" "Q75" "Q975"];SecInf_stats], ',')
        writedlm(file_all_sim_inf, SIMS, ',')
    end

    ##
    ## Display results
    ##
    # include("../Display_and_plot_results/TimeDist_to_Ncases/DisplayAndPlot_results.jl")
end # timer (save)

##
## Computation time
##
# Display timer outputs
println("\n\n\tCOMPUTATION TIME\n")
show(to; allocations=false)
println(" ")