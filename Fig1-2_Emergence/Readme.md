## Figures 1 and 2: Emergence of the epidemic outbreak

Code to reproduce Figs 1 and 2 of the main text.

- `Run_EstimateEmergence` 
    <br> Main script running the infection and detection model on each epidemiological context (see the pseudo-algorithm in the Appendix).
    <br> User-determined options:</br>
    - Select the epidemiological context (`Alpha_UK` or `COVID-19_Wuhan`; lines 30-31)
    - Select how many accepted simulations to retrieve (line 34)
    - Select wether or not to save the results (lines 40-41)


- `Plot_Results_Emergence` 
    <br> Script plotting Figures 1 (Alpha) and 2 (COVID-19):
    - Epidemic curve from the case data and the 
    - Violin plots for the estimates of the emergence date
    - Violin plots for the results from previous studies: 
        - <strong>Alpha variant in the UK:</strong> We compare our results to i) an updated version of the results on emergence presented in (Czuppon et al., 2021), by setting R=1.9 (as in Hill et al., instead of 2.5) and using a negative-binomial distribution for the secondary infections (instead of Poisson) and ii) the tMRCA estimated by (Hill et al., 2021; personal communication).
        - <strong>COVID-19 cases in Wuhan:</strong> We compare our results to the estimates of emergence presented in (Pekar et al., 2022; personal communication)


### References
Czuppon P, Schertzer E, Blanquart F and Débarre F (2021). “The Stochastic Dynamics of Early Epidemics: Probability of Establishment, Initial Growth Rate, and Infection Cluster Size at First Detection.” Journal of The Royal Society Interface 18 (184): 20210575. <a href="https://doi.org/10.1098/rsif.2021.0575" rel="_blank">10.1098/rsif.2021.0575</a>.

Hill V, Du Plessis L, Peacock TP, Aggarwal D, Colquhoun R, Carabelli AM, Ellaby N et al. (2022). “The Origins and Molecular Evolution of SARS-CoV-2 Lineage B.1.1.7 in the UK.” Virus Evolution, August, veac080. DOI: <a href="https://doi.org/10.1093/ve/veac080" rel="_blank">10.1093/ve/veac080</a>.

Pekar JE, Magee A, Parker E, Moshiri N, Izhikevich K, Havens JL, Gangavarapu K, et al. (2022). “The Molecular Epidemiology of Multiple Zoonotic Origins of SARS-CoV-2.” Science, July, eabp8337. DOI: <a href="https://doi.org/10.1126/science.abp8337" rel="_blank">10.1126/science.abp8337</a>.


