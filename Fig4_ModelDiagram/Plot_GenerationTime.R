library(ggplot2)
library(tidyverse)

##
## Data ############
##
GammaDist = tibble(x=seq(0, 12, by = 0.05),
              y=dgamma(x, shape=6.6, scale=0.8))


write_csv(GammaDist, "Intro figures/ModelDiagram_TikZ/GammaDist.csv")

##
## Plot ############
##
p = ggplot(data=GammaDist, 
           aes(x=x, y=y)) +
    geom_line() + 
    labs(x="", y="Generation time distribution") + 
    theme_classic() +
    NULL
p