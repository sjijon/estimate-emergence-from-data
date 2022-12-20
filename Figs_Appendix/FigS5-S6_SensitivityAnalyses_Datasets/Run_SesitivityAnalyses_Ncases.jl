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

## Load user functions and other setup scripts
include("../Processes/InfectionProcess.jl")
include("../Processes/DetectionProcess.jl")
include("../Processes/DefineStructs.jl")

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
repeats = 5000
# repeats = 1000
# repeats = 10
# repeats = 1

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
tol_delay = 0.9
# abs(Difference in daily cases) <= tol_epi
tol_epi = 0.3   

## Create object of observations
ObsCases = ObsDetect(Date.(data[:,1]),data[:,2],N_cases,Date_N)

## Data
obs_cases_cumul = Array{Int64}(cumsum(ObsCases.cases)) # Cumulative number of cases
obs_cases_daynum = Dates.value.(ObsCases.dates- ObsCases.dates[1]) .+ 1 # Day number

## Time parameters (in days)
dt = 1/10  # time step in the simulation (only<1 day alllowed)
length_inf_vec = trunc(Int,max_t_infect/dt) # final vector length

##
## 1c. Define filepaths ########################
##
if SaveResults == "Yes"
    ## Create path
    dir_output = string("Fig1-2_Emergence/Output/",EpiContext)
    mkpath(dir_output)

    ## Per-day cumulative number of cases (per SimEpi)
    file_cumul_cases = string(dir_output,"/cumul_",N_cases,"cases_",repeats,"sims.csv")
    ## Per-day distance to observations
    file_dist_pw = string(dir_output,"/dist_pw_",N_cases,"cases_",repeats,"sims.csv")
    ## Open file                
    open(file_cumul_cases, "w")
    open(file_cumul_cases, "a") do io
    writedlm(io,[1:1:90],",")
    end
    open(file_dist_pw, "w")

    ## Time and epi size
    file_Cases_EpiSize_Time = string(dir_output,"/Cases_EpiSize_Time_",N_cases,"cases_",repeats,"sims.csv")
    ## Simulated epidemics
    file_all_sim_inf = string(dir_output,"/SimInf_",N_cases,"cases_",repeats,"sims.csv")
    ## Secondary infections (stats of means)
    file_SecInf_stats = string(dir_output,"/SecInf_stats_",N_cases,"cases_",repeats,"sims.csv")
end

##
## 2. Run simulations ################################
##

## Initialization of counters
run_num = 0                                 # repetition index
successes = 0                               # successes counter

## Initialization of results arrays
Cases_EpiSize_Time = Array{Int64}(undef,0,4) # Time of detection, number of cases and size of the epidemic
SIMS = [1:max_t_infect;]                    # All simulations
SecondaryInfec_all = Array{Int64}(undef,0,2)# All number of secondary infections
SecInf_stats = Array{Int64}(undef,0,4)      # All number of secondary infections (stats for each SimEpi)

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

        ## Continue if SimCases >= N_cases
        if SimCases.cumul[end]>=N_cases
            ## Total epidemic size at the time of infection of the N-th case
            global EpiSize_NthCaseInfect = SimEpi.daily_infec[1:Int(SimCases.d_detect[N_cases])] |> sum
        
            ##
            ## 2.2. Model calibration
            ##
            ## Constraints for selecting a simulation:
            ## a) The time period between the 1st infection and the first observed case
            global obs_num_days = Dates.value(Date_N - Date_1)
            global sim_num_days = length(SimCases.cumul[SimCases.cumul.>0])
            global Diff_SimInf1_ObsCas1 = Int(SimCases.d_detect[end] - obs_num_days)
            ## b) The length of the time period where cases occur
            global Diff_NumDays = abs(obs_num_days - sim_num_days)
            ## c) The similarity with the observed cumulative number of cases
            ## Add zeros if&where needed and align the cumulative curves at right (N-th case)
            global obs_cases_cumul_ = [zeros(Int,maximum([length(SimCases.cumul),length(obs_cases_cumul)])-length(obs_cases_cumul));obs_cases_cumul]
            global sim_cases_cumul_ = [zeros(Int,maximum([length(SimCases.cumul),length(obs_cases_cumul)])-length(SimCases.cumul));SimCases.cumul]
            ## Compute the maximum of the absolute difference between the simulated and the observed cumulative cases
            global Dist_pw = abs.(sim_cases_cumul_.-obs_cases_cumul_)
            global Dist = maximum(Dist_pw)

            ## Select the simulation if all conditions are verified
            if (Diff_SimInf1_ObsCas1>=0 &&  sim_num_days>= tol_delay*obs_num_days && Dist<(tol_epi*N_cases))
                # Update number of successes
                global successes += 1
                
                ##
                ## 2.4. Build results arrays
                ##
                ## Time from the first infection and epidemic size at time where the N-th case is determined
                global Cases_EpiSize_Time = [Cases_EpiSize_Time;[SimCases.d_NthCase SimCases.cumul[end] EpiSize_NthCaseInfect] SimCases.d_detect[1]]
                # println("Cases_EpiSize_Time: $Cases_EpiSize_Time")

                ## Save all simulated epidemics
                global SIMS = hcat(SIMS,SimEpi.daily_infec)
                
                ## Save all number of secondary infections 
                ## up to the date of the N-th case
                global SecInf_dNthCase= SimEpi.second_infec[SimEpi.second_infec[:,1].<SimCases.d_NthCase,:]
                local SecInf_mean = mean(SecInf_dNthCase[:,2])
                local SecInf_IQR = quantile(SecInf_dNthCase[:,2], [0.025 0.5 0.975])
                global SecInf_stats = vcat(SecInf_stats,[SecInf_mean SecInf_IQR])
                
                ## Save
                if SaveResults == "Yes"
                    ## All cumulative cases
                    open(file_cumul_cases, "a") do io
                        writedlm(io,SimCases.cumul',",")
                    end

                    # All distances (pointwise) to observations
                    open(file_dist_pw, "a") do io
                        writedlm(io,Dist_pw',",")
                    end
                end

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

Prop_SelectedSims = round(100*repeats/run_num,digits=2)
println("\n\nNumber of runs (run_num): $run_num\nSuccesses = $Prop_SelectedSims % of runs")

##
## 5. Display results and save ####################
##

if SaveResults == "Yes"
    ## Save
    writedlm(file_Cases_EpiSize_Time, [["MinTime" "Cases" "EpiSize" "Case1"]; Cases_EpiSize_Time], ',')
    writedlm(file_SecInf_stats, [["Mean" "Q025" "Q50" "Q975"];SecInf_stats], ',')
    writedlm(file_all_sim_inf, SIMS, ',')
end

##
## Display results
##
include("../Processes/DisplayResults.jl")
