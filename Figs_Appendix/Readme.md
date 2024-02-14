## Supplementary figures (Appendix)

- `FigS1-S2_CaseData` contains the script (`.R`) to plot the epicurves.
- `FigS3-Reanalysis_Roberts_etal_2021` contains the script (`.R`) to reproduce and extend the analysis presented in Roberts et al. (2021)
- `FigS4-S5_CumulativeCases` contains the script (`.R`) to plot a subset of curves of cumulative cases of the retained epidemics.
- `FigS6-S7_Error` it contains the scripts to run the model on simulated data (generated using `GenerateSinData.jl`). It also contains the scripts (`.R`) to analyze the error in the estimates (difference in dates) by number of cases considered.
- `FigS8-S89_SensitivityAnalyses_Datasets` contains
    - `SetParams_Supplement` where there are other datasets on which we ran our model 
    - `Run_SesitivityAnalyses_Ncases.jl` runs our model on different datasets other than those presented in the main text. Please select the epi context (Alpha variant in the UK or COVID-19 in Wuhan in lines 30-31).
    - `Plot_SesitivityAnalyses_Ncases_Alpha_UK.R` plots and compares the results (Fig S5).
    - `Plot_SesitivityAnalyses_Ncases_COVID-19_Wuhan.R` plots and compares the results (Fig S6).
