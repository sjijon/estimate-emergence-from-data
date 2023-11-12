#### Jijón S, Czuppon P, Blanquart F & Débarre F (2023). 
#### Using early detection data to estimate the date of emergence of an epidemic outbreak.
#### https://github.com/sjijon/estimate-emergence-from-data
####
#### Plot results from running the model on synthetic data
####
#### 0. SETUP ################################################################
####

## CLEAR ALL
rm(list=ls())

# Load data
allestims <- read.csv("results/All_Estims.csv")

# Error without taking the absolute value
allestims$DiffDelay <- allestims$ObsDaysInf1CaseN - allestims$MinTime

# Values of N
Ns <- sort(unique(allestims$FirstCases))

# Colors for plotting
colsN <- c("#D81B60", "#1E88E5", "#004D40")
stopifnot(length(colsN) == length(Ns))
names(colsN) <- as.character(Ns)

# Function to make colos transparent
t_col <- function(colors, opacity) {
  vapply(colors, function(color){
    rgb.val <- col2rgb(color)
    rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (opacity) * 255)
  }, FUN.VALUE = "x")
}

# Define breaks for histograms (so that they look the same)
brks <- seq(min(allestims$DiffDelay)-0.5, max(allestims$DiffDelay) + 0.5, by = 1)

# Plot the result
par(las = 1, mgp = c(2, 0.5, 0), tck = -0.025)
layout(matrix(1:3, ncol = 1))
for(i in 1:3){
  # Subset of the data for this N
  tmp <- allestims[allestims$FirstCases == Ns[i], ]
  # Plot histogram
  hist(tmp$DiffDelay, breaks = brks, 
       col = t_col(colsN[i], 0.5), border = "white", lwd = 0.1,
       xlab = "Difference from real date of first infection", 
       main = paste0("N = ", Ns[i]), axes = FALSE)
  # All line for median
  abline(v = median(tmp$DiffDelay), col = colsN[i], lwd = 2)
  # Add axes
  axis(1, pos = 0)
  axis(2, pos = brks[1])
}
