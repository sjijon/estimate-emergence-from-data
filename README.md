# Estimate emergence from data
Code for reproducing the results from the preprint: 

<strong>Using early detection data to estimate the date of emergence of an epidemic outbreak</strong>

Sofía Jijón<sup>1</sup>, Peter Czuppon<sup>2</sup>, François Blanquart<sup>3</sup> and Florence Débarre<sup>4</sup>

<sup>1</sup>Institute of ecology and environmental sciences of Paris (iEES Paris, UMR 7618), Sorbonne Université, CNRS, UPEC, IRD, INRAE, Paris 75005, France

<sup>2</sup>Institute for Evolution and Biodiversity, University of Münster, Germany

<sup>3</sup>Center for Interdisciplinary Research in Biology, CNRS, Collège de France, PSL Research University

<strong>URL:</strong> 
<a href="https://www.medrxiv.org/content/10.1101/2023.01.09.23284284v1" >Preprint</a> | 
<a href="https://doi.org/10.5281/zenodo.10623614" >Static code</a> | 
<!-- <a href="" >Article</a> -->

## Abstract

While the first infection by an emerging disease is often unknown, information on early cases can be used to date it, which is of great interest to trace the disease's origin and understand early infection dynamics. Here, we extend the method presented in (Czuppon et al., 2021) to estimate the time series from emergence (i.e., infection of the first human case by the focal disease or variant) to the detection of N cases. We run numerical simulations of infectious disease spreading from a single infectious individual and calibrate the model to reproduce the observed cases in two main epidemiological contexts (the emergence of the Alpha SARS-CoV-2 variant in the UK and the early cases of COVID-19 in Wuhan), but the code was built in a generic form to facilitate applications to other datasets. The main outcome of our model is the expected time to N detected cases.

## Contents

- `Data/` contains the files containing the time series of the early cases we use to calibrate the model for different epidemiological contexts, as well as the emergence or tMRCA estimates from previously published studies.

- `Routines/` contains the functions used repeatedly in each application.

- `SetParameters_EpiContext/` contains the parameterizations for each application of the model

- `Fig1_FirstCaseDiagrama` generates a diagram illustrating the notions of extinct trees, first infection of the ongoing tee and the MRCA in an epidemic outbreak.

- `Fig2-3_Emergence/` contains the codes to reproduce the results (`.jl`) and figures (`.R`) presented in the main text. The subfolder `ABC_Simulations/` contains the files containing the time series of the early cases to estimate the first infection in the COVID-19 example. 

- `Fig4_Sensitivity Analyses/` contains a variation of the previous scripts to run the simulation over varying values of the main parameters.

- `Fig5_ModelDiagram/` plots a simplified diagram of the model (Ti<emph>k</emph>Z, `.tex`).

- `Figs_Appendix/` contain scripts to reproduce the figures presented in the article's appendix

See the `Readme.md` files contained in each folder for more details.

## Implementation details

Simulations are built in `Julia` version 1.8.0. Simulations for the Approximate Bayesian Computation are built in `C++` version 17. 

Figures are built in `R` version 4.1.2, using the `ggplot2` package version 3.3.

