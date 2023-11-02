## Jijón S, Czuppon P, Blanquart F & Débarre F (2023). 
## Using early detection data to estimate the date of emergence of an epidemic outbreak.
## https://github.com/sjijon/estimate-emergence-from-data
##
## DEFINE STRUCTS
##
## Define the structures for the simulated and observed epidemics
#################################################
## Observed cases
struct ObsDetect
    dates::Vector{Date}   # Dates of observations
    cases::Array{Float64} # Observed cases (>=1)
    N::Int64              # Number of first cases
    Date_N::Dates.Date    # Date of N-th case
end

## Simulated epidemic
struct SimInfect
    epidemic::Bool              # Epidemic? Yes (true) or No (false)
    infections::Vector{Int64}   # Vector of infections (time)
    daily_infec:: Vector{Int64} # Infections (days)
    second_infec::Array{Any}    # Secondary infections
end

## Simulated case-declarations
struct SimDetect
    num::Int64                  # Number of cases
    d_NthCase::Int64            # Day of detection of N-th case
    d_detect::Vector{Int64}     # Day of detection
    cumul::Vector{Int64}        # Cumulative number of cases
end