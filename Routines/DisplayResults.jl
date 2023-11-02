## Jijón S, Czuppon P, Blanquart F & Débarre F (2023). 
## Using early detection data to estimate the date of emergence of an epidemic outbreak.
## https://github.com/sjijon/estimate-emergence-from-data
##
## DISPLAY RESULTS
##
## Median and 95% Inter quantile ranges (IqR; values between 2.5th and 97.5th percentiles)
#################################################
println("\n\n.............Display and plot results...............")

##
## 0. EPI CONTEXT #########################################################
##
println("\n\nEpidemiological context: $EpiContext") 

println("R0=$R0; kappa=$kappa; p_detect=$p_detect; delta_tol_epi=$tol_epi\t\n")

##
## 0. DATA #########################################################
##
println("\n\n1. OBSERVED DATA")
println("\n1st case detected on $(Date(Date_1)), \nN-th case detected on $(Date(Date_N))\nthat is, $(Date_N-Date_1+Day(1)) between 1st and N-th cases\n")

println("Total number of cases: $N_cases")

println("\n\n2. MODEL CALIBRATION\n")
##
## 1. DISPLAY RESULTS #########################################################
## The results array us of the following form:
##
##                      col1                col2                col3                    col4
##                      Day of N-th case    Number of cases     Epidemic size           Day of 1st Case
## Cases_EpiSize_Time = [SimCases.d_NthCase  SimCases.cumul[end] EpiSize_NthCaseInfect   SimCases.d_detect[1]]

Delay_Inf1CaseN = Cases_EpiSize_Time[:,1]    # Number of days between 1syt infection and N-th case
NumCases_tN = Cases_EpiSize_Time[:,2]        # Number of cases the day of the N-th detection (d_K)
EpiSize_dK = Cases_EpiSize_Time[:,3]         # Epi size the day of infection of the N-th case (t_N)
Day_Case1 = Cases_EpiSize_Time[:,4]          # Time period between case 1 and case N

##
## Display results
##
println("\n\n3. RESULTS")

##
## Time period between 1st and N-th case
##
Delay_Case1CaseN = Delay_Inf1CaseN .- Day_Case1 .+ 1

println("\nDays between 1st and the N-th cases\n")
println("\tMedian: $(trunc.(Int,median(Delay_Case1CaseN))) (95%IqR: $(trunc.(Int,quantile(Delay_Case1CaseN,0.025)))-$(trunc.(Int,quantile(Delay_Case1CaseN,0.975))))")

##
## Date of 1st case
##
Date_Case1 = Dates.Date(Date_N) .- Dates.Day.(trunc.(Int,quantile(Delay_Case1CaseN, [0.50 0.025 0.975]))) 
println("\twhich dates the 1st case at : $(Date_Case1[1]) (95%IqR: $(Date_Case1[3]), $(Date_Case1[2]))") 

##
## Time period before 1st detection
##
Delay_Inf1Case1 = trunc.(Int,quantile(Delay_Inf1CaseN.-Delay_Case1CaseN, [0.50 0.025 0.975]))
## = trunc.(Int,quantile(Day_Case1, [0.50 0.025 0.975])) .- 1
println("\nTime elapsed until 1st case: \n\n\t$(Delay_Inf1Case1[1]) (95%IqR: $(Delay_Inf1Case1[2])-$(Delay_Inf1Case1[3]))") 

##
## Time between 1st infection and  N-th case
##

# Table of estimates
Estims = [  "Median"    round(quantile(Cases_EpiSize_Time[:,1], 0.50),digits=2);
            "Q25"       round(quantile(Cases_EpiSize_Time[:,1], 0.25),digits=2);
            "Q75"       round(quantile(Cases_EpiSize_Time[:,1], 0.75),digits=2);
            "Q025"      round(quantile(Cases_EpiSize_Time[:,1], 0.025),digits=2);
            "Q975"      round(quantile(Cases_EpiSize_Time[:,1], 0.975),digits=2);
            "Mean"      round(mean(Cases_EpiSize_Time[:,1]),digits=2);
            "StD"       round(std(Cases_EpiSize_Time[:,1]),digits=2)]

println("\nNumber of days between first infection and the time of occurrence of $N_cases cases\n")
            
## Median, IqR
println("\tMedian: $(Estims[1,2]) (IqR: $(Estims[2,2])-$(Estims[3,2])),")
## date (median):
## assuming Date_N is the date where the N-th case ocurred, 
## we go back in time using the distribution of the number of days
## between first infection and the N-th case
OriginDate = Dates.Date(Date_N) .- Dates.Day.(trunc.(Int,quantile(Cases_EpiSize_Time[:,1], [0.50 0.25 0.75])))
println("\twhich dates origin at : $(OriginDate[1]) ($(OriginDate[3]), $(OriginDate[2]))") 

## Median, 95%IqR
println("\n\tMedian: $(Estims[1,2]) (95%IqR: $(Estims[4,2])-$(Estims[5,2])),")
## date (median):

## Date_N is the date where the N-th case ocurred, 
## we go back in time using the distribution of the number of days
## between first infection and the N-th case
OriginDate = Dates.Date(Date_N) .- Dates.Day.(trunc.(Int,quantile(Cases_EpiSize_Time[:,1], [0.50 0.025 0.975])))
println("\twhich dates origin at : $(OriginDate[1]) (95%IqR: $(OriginDate[3]) to $(OriginDate[2]))") 
## Min
println("\tand not earlier than $(Dates.Date(Date_N) .- Dates.Day.(maximum(Cases_EpiSize_Time[:,1])))") 

##
## Number of cases at day where N-th case occurs
##
println("\nNumber of cases at day where N-th case occurs\n")
## Rounded to multiples of 100
println("\tMedian: $(Int(round(median(NumCases_tN),sigdigits=2))) (95%IqR: $(Int(round(quantile(NumCases_tN,0.025),sigdigits=3)))-$(Int(round(quantile(NumCases_tN,0.975),sigdigits=2))))")



##
## Epidemic
##
println("\nEpidemic size at day of infection of N-th case\n")
## Rounded to multiples of 100
println("\tMedian: $(Int(round(median(Cases_EpiSize_Time[:,3]),sigdigits=3))) (95%IqR: $(Int(round(quantile(Cases_EpiSize_Time[:,3],0.025),sigdigits=3)))-$(Int(round(quantile(Cases_EpiSize_Time[:,3],0.975),sigdigits=3))))")

##
## Hiddden epidemic 
##
println("\nHidden epidemic \n")
println("\tMedian (95%IqR): $(round(100*(1-median(NumCases_tN./EpiSize_dK)),digits=2)) ($(round(100*(1-quantile(NumCases_tN./EpiSize_dK,0.025)),digits=2))-$(round(100*(1-quantile(NumCases_tN./EpiSize_dK,0.975)),digits=2)))")

##
## Proportion of detected infections
##
println("\n\tProportion of detected:")
println("\tMedian (95%IqR): $(round(100*(median(NumCases_tN./EpiSize_dK)),digits=2)) ($(round(100*(quantile(NumCases_tN./EpiSize_dK,0.025)),digits=2))-$(round(100*(quantile(NumCases_tN./EpiSize_dK,0.975)),digits=2)))")