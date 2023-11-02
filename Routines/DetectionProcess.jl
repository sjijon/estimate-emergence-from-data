## Jijón S, Czuppon P, Blanquart F & Débarre F (2023). 
## Using early detection data to estimate the date of emergence of an epidemic outbreak.
## https://github.com/sjijon/estimate-emergence-from-data
##
## DETECTION PROCESS
##
## Function building the cases time series
#################################################
function DetectionProcess(SimEpi)
    ## Initialize the number of cases
    global case_num = 0                     # Number of detected cases
    global d_detect_all = Vector{Int64}() # All detection times

    Infections = SimEpi.daily_infec

    ## Draw the per-day number of detections among infections
    for d in range(start=1,stop=length(Infections))
        ## Number of detections ~ Binom(n=Infected(d),p=prob. of detection)
        SamplingDist = Binomial(Infections[d], p_detect)

        ## Number of detected among infected I(t)
        n_detect = rand(SamplingDist,1)[1]
        # println("n_detect: $n_detect")

        ## Detection time
        if n_detect>0
            ## Compute detection time for each individual during [t,t+dt]
            d_detect = trunc.(Int, d .+ rand(SampTimeDist,n_detect))
            # println(d_detect,"\n")
            
            ## Add to vector of detection times (before max_t_detect)
            d_detect_all = append!(d_detect_all,d_detect)
            d_detect_all = sort(d_detect_all)
            # println("Detection times: $d_detect_all\n")

            ## Increase the number of cases
            global case_num+=n_detect
            # println("case_num: $case_num")

            ## Run until the same day of the N-th observed case
            if case_num >= N_cases
                ## Day of ocurrence of the N-th case
                global d_NthCase = d_detect_all[N_cases]
                # println("day min to detect N cases = $(d_NthCase)")

                ## Ensure to keep the earliest cases
                if d >= d_NthCase
                    ## Keep only the events occuring on the same day as N-th case
                    d_detect_all = d_detect_all[d_detect_all.<=d_NthCase]
                    ## Cumulative cases
                    global cum_sum = d_detect_all |> counts |> cumsum
                    ## Exit the for loop
                    break
                end # if (date)
                
                ## Cumulative cases
                global cum_sum = d_detect_all |> counts |> cumsum
            end # if (number of cases)
        end # if (detection time)
    end # for
        
    return SimDetect(cum_sum[end],d_NthCase,d_detect_all,cum_sum)
end # Function
    