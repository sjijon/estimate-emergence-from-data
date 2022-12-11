# estimate-emergence-from-data
Code for reproducing the results from the preprint : 

<strong>Estimating the date of emergence of an epidemic outbreak from detection data: Applications to COVID-19</strong>

Sofía Jijón<sup>1</sup>, Peter Czuppon<sup>2</sup>, François Blanquart<sup>3</sup> and Florence Débarre<sup>4</sup>

<sup>1</sup>Institute of ecology and environmental sciences of Paris (iEES Paris, UMR 7618), Sorbonne Université, CNRS, UPEC, IRD, INRAE, Paris 75005, France

<sup>2</sup>Institute for Evolution and Biodiversity, University of Münster

<sup>3</sup>Center for Interdisciplinary Research in Biology, CNRS, Collège de France, PSL Research University

<strong>URL:</strong> 
<a href="" >Preprint</a> | <a href="" >Static code</a>

## Abstract

While the first infection by an emerging disease is often unknown, information on early cases can be used to date it, which is of great interest to trace the disease's origin and understand early infection dynamics. Here, we extend the method presented in Czuppon et al. 2021 to estimate the time series from emergence (i.e., infection of the first human case by the focal disease or variant) to the detection of N cases. We run numerical simulations of infectious disease spreading from a single infectious individual and calibrate the model to reproduce the observed cases in two main epidemiological contexts (the emergence of the Alpha SARS-CoV-2 variant in the UK and the early cases of COVID-19 in Wuhan), but the code was built in a generic form to facilitate applications to other datasets. The main outcome of our model is the expected time to the detection of N cases.

## Contents

- `Data` contains the files containing the time series of the early cases we use to calibrate the model for different epidemiological contexts, as well as the emergence or tMRCA estimates from previously published studies.

- `Fig 1-2 and Fig3` contain the codes to reproduce the results and figures presented in the article. See more details in the readme.md files contained in each folder. 

- `Processes` contains the functions used repeatedly in each application.

- `SetParameters_EpiContext` contains the parametrizations for each application of the model