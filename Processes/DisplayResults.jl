## Estimating the time distribution between the first infection and N cases
##
## Jijon S, Czuppon P, Blanqart F., Débarre F
## iEES, 2022
##
##
## DISPLAY AND PLOT RESULTS
##
## Display results
#################################################
println("\n\n.............Display and plot results...............")

##
## 0. EPI CONTEXT #########################################################
##
println("\n\nEpidemiological context: $EpiContext") 

println("\n\n1. OBSERVED DATA")
println("\n1st case detected on $(Date(Date_1)), \nN-th case detected on $(Date(Date_N))\nthat is, $(Date_N-Date_1+Day(1)) between 1st and N-th cases\n")

println("Total number of cases: $N_cases")

##
## 1. RESULTS #########################################################
println("\n\n3. RESULTS")

##
## Time period between 1st and N-th case
##
println("\nDays between 1st and the N-th cases\n")
println("\tMedian: $(median(Cases_EpiSize_Time[:,4])) (95%CrI: $(quantile(Cases_EpiSize_Time[:,4],0.025))-$(quantile(Cases_EpiSize_Time[:,4],0.975)))")
Case1Date = Dates.Date(Date_N) .- Dates.Day.(trunc.(Int,quantile(Cases_EpiSize_Time[:,1].-Cases_EpiSize_Time[:,4], [0.50 0.025 0.925]))) 
println("\twhich dates the 1st case at : $(Case1Date[1]) (95%CrI: $(Case1Date[3]), $(Case1Date[2]))") 

##
## Minimum time until N-th case
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
            
## Median, IQR
println("\tMedian: $(Estims[1,2]) (IQR: $(Estims[2,2])-$(Estims[3,2])),")
## date (median):
## assuming Date_N is the date where the N-th case ocurred, 
## we go back in time using the distribution of the number of days
## between first infection and the N-th case
OriginDate = Dates.Date(Date_N) .- Dates.Day.(trunc.(Int,quantile(Cases_EpiSize_Time[:,1], [0.50 0.25 0.75])))
println("\twhich dates origin at : $(OriginDate[1]) ($(OriginDate[3]), $(OriginDate[2]))") 
## Mean
println("\n\tMean: $(Estims[6,2]) (sd: $(Estims[7,2]))")  
OriginDate_mean = Dates.Date(Date_N) - Dates.Day(trunc.(Int,mean(Cases_EpiSize_Time[:,1])))
println("\twhich dates origin at : $OriginDate_mean") 
## Median, 95%CrI
println("\n\tMedian: $(Estims[1,2]) (95%CrI: $(Estims[4,2])-$(Estims[5,2])),")
## date (median):

## Date_N is the date where the N-th case ocurred, 
## we go back in time using the distribution of the number of days
## between first infection and the N-th case
OriginDate = Dates.Date(Date_N) .- Dates.Day.(trunc.(Int,quantile(Cases_EpiSize_Time[:,1], [0.50 0.025 0.975])))
println("\twhich dates origin at : $(OriginDate[1]) (95%CrI: $(OriginDate[3]) to $(OriginDate[2]))") 
## Min
println("\tand not earlier than $(Dates.Date(Date_N) .- Dates.Day.(maximum(Cases_EpiSize_Time[:,1])))") 

##
## Number of cases at day where N-th case occurs
##
println("\nNumber of cases at day where N-th case occurs\n")
## Rounded to multiples of 100
println("\tMedian: $(Int(round(median(Cases_EpiSize_Time[:,2]),sigdigits=2))) (95%CrI: $(Int(round(quantile(Cases_EpiSize_Time[:,2],0.025),sigdigits=3)))-$(Int(round(quantile(Cases_EpiSize_Time[:,2],0.975),sigdigits=2))))")


##
## Epidemic
##
println("\nEpidemic size at day of infection of N-th case\n")
## Rounded to multiples of 100
println("\tMedian: $(Int(round(median(Cases_EpiSize_Time[:,3]),sigdigits=3))) (95%CrI: $(Int(round(quantile(Cases_EpiSize_Time[:,3],0.025),sigdigits=3)))-$(Int(round(quantile(Cases_EpiSize_Time[:,3],0.975),sigdigits=3))))")

##
## Hiddden epidemic 
##
println("\nHidden epidemic \n")
println("\tMedian (95%CrI): $(round(100*(1-median(Cases_EpiSize_Time[:,2]./Cases_EpiSize_Time[:,3])),digits=2)) ($(round(100*(1-quantile(Cases_EpiSize_Time[:,2]./Cases_EpiSize_Time[:,3],0.025)),digits=2))-$(round(100*(1-quantile(Cases_EpiSize_Time[:,2]./Cases_EpiSize_Time[:,3],0.975)),digits=2)))")

##
## Proportion of detected infections
##
println("\n\tProportion of detected:")
println("\tMedian (95%CrI): $(round(100*(median(Cases_EpiSize_Time[:,2]./Cases_EpiSize_Time[:,3])),digits=2)) ($(round(100*(quantile(Cases_EpiSize_Time[:,2]./Cases_EpiSize_Time[:,3],0.025)),digits=2))-$(round(100*(quantile(Cases_EpiSize_Time[:,2]./Cases_EpiSize_Time[:,3],0.975)),digits=2)))")

##
## Number of secondary infections
##
println("\nSecondary infections\n")
SecondaryInfec_all_g0 = SecondaryInfec_all[SecondaryInfec_all[:,2].>0,:]
println("\tMean: $(round(mean(SecondaryInfec_all_g0[:,2]),digits=3)) (95%CrI: $(round(quantile(SecondaryInfec_all_g0[:,2],0.025),digits=3))-$(round(quantile(SecondaryInfec_all_g0[:,2],0.975),digits=3)))")