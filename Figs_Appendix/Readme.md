## Supplementary figures (Appendix)

- `FigS1-S2_CaseData` contains the script (.R) to plot the epicurves.
- `FigS3-S4_CumulativeCases` contains the script (.R) to plot a set of curves of cumulative cases os simulated epidemics.
- `FigS5-S6_SensitivityAnalyses_Datasets` contains
    - `SetParams_Supplement` where there are other datasets on which we ran out model 
    - `Run_SesitivityAnalyses_Ncases.jl` runs our model on different datasets other than those presented in the main text. Please select the epi context (Alpha variant in the UK or COVID-19 in Wuhan in lines 30-31).
    - `Plot_SesitivityAnalyses_Ncases_Alpha_UK.R` plots and compares the results (Fig S5).
    - `Plot_SesitivityAnalyses_Ncases_COVID-19_Wuhan.R` plots and compares the results (Fig S6).