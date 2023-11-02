## Jijón S, Czuppon P, Blanquart F & Débarre F (2023). 
## Using early detection data to estimate the date of emergence of an epidemic outbreak.
## https://github.com/sjijon/estimate-emergence-from-data
##
## TRANSMISSION
##
## Function that builds the epidemic:
##      - Epidemic              (boolean)
##      - NewInfections            (vector)
##      - Secondary infections  (vector)
#################################################
function InfectionProcess()
    ## Initialization
    global t = 0                                # Time (increases by dt)
    global t_indx = 1                           # Working index (increases by 1)
    length_inf_vec = trunc(Int,max_t_infect/dt) # final vector length
    global NewInfections = Vector{Int64}(zeros(length_inf_vec))   # Number of infected
    NewInfections[1] = 1                        # 1 infected individual at time 0
    global SecondaryInfec = Array{Any}(undef,0,2)     # All secondary infections
    global DailyInfections = Vector{Int64}(zeros(max_t_infect))   # Daily number of infected

    ## Run until max time
    while (sum(NewInfections)<=min_n_infect && t<max_t_infect)
        ## Display t
        # println("\nt: $t")
        # if round(t,digits=2) % 10 == 0
        #     println("\nt: $(round(t,digits=2))")
        # end

        ## Number of secondary infections
        offsp = rand(TransDist,NewInfections[t_indx])
        n_offsp = sum(offsp)
        # println("n_offspring: $n_offsp")

        ## Loop over n_offsp from individuals that were added at time t
        if n_offsp>0
            ## Time of secondary infection occurrence
            t_infec = t .+ rand(InfTimeDist,n_offsp)
        
            ## Get time index
            # IMPORTANT THAT dt = 10^(-k)! otherwise adapt this rounding formula!
            t_infec_indx = trunc.(Int,t_infec/dt)
            # println("t_infec_indx: $t_infec_indx")
        
            # Increase number of infected at time t_infec_indx
            for i in 1:n_offsp    
                if t_infec_indx[i] < length_inf_vec
                    global NewInfections[t_infec_indx[i]] += 1
                end
            end

            ## Save the number of secondary infections
            SecondaryInfec = vcat(SecondaryInfec,hcat(repeat([t_indx*dt],length(offsp)),offsp))
        end # if (offspring)
        
        # Increase time
        t += dt        
        t_indx += 1
        # println("t_indx: $t_indx")
    end # while (time)
    
    ## Daily infections
    for d in 1:Int(length_inf_vec*dt)
        di=Int((d-1)/dt)+1
        df=Int((d/dt))
        DailyInfections[d] = sum(NewInfections[di:df])
    end

    ## Epi size
    if sum(NewInfections) >= min_n_infect
        ## Epidemic!
        epidemic = true

        return SimInfect(epidemic,NewInfections,DailyInfections,SecondaryInfec)
    else
        return SimInfect(false,[],[],[])
    end # Success
end # Function