####
#### Early epidemics
#### iEES
####
#### January, 2024
####
#### Plot Results from running the model on synthetic data
####
################################################################

## CLEAR ALL
rm(list=ls())

# Load data
folder_name <- "Figs_Appendix/FigS6-S7_Error"
allestims <- read.csv(paste0(folder_name,"/All_Estims.csv"))

# Error without taking the absolute value
allestims$DiffDelay <- allestims$ObsDaysInf1CaseN - allestims$MinTime

# Sample size
table(allestims$Dataset, allestims$FirstCases)

# Values of N
Ns <- sort(unique(allestims$FirstCases))

# Colors for plotting
colsN <- c("#D81B60", "#1E88E5", "#004D40")
stopifnot(length(colsN) == length(Ns))
names(colsN) <- as.character(Ns)


# Simulation number
allestims$Sim <- vapply(allestims$Dataset, function(x) as.numeric(substr(x, 5, nchar(x))), 1)

allq <- aggregate(allestims$DiffDelay, by = list(N = allestims$FirstCases, Sim = allestims$Sim), FUN = function(x) quantile(x, c(0.025, 0.5, 0.975)))

xr <- range(allq$x)

# One plot for each N ####
# Plot the result
for(j in seq_along(Ns)){
  pdf(file=paste0(folder_name,"/distribution_N_bySim_", Ns[j], ".pdf"), width = 7, height = 7)
  
  theN <- Ns[j]
  
  par(las = 1)
  sub <- allq[allq$N == theN, ]
  includes0 <- vapply(seq_len(nrow(sub)), function(i){
    ifelse((0 >= sub$x[i, "2.5%"]) & (0 <= sub$x[i, "97.5%"]), 1, 0)
  }, FUN.VALUE = 1)
  
  plot(sub$x[, "50%"], sub$Sim, xlim = xr, pch = 1 + includes0 * 15, col = colsN[j], 
       xlab = "Difference between estimated date and real date (median and 95% interval)", ylab = "", axes = FALSE)
  axis(1)
  segments(x0 = sub$x[, "2.5%"], x1 = sub$x[, "97.5%"], 
           y0 = sub$Sim, y1 = sub$Sim, col = colsN[j])
  abline(v = 0, lty = 2)
  title(main = paste0("N = ", theN))
  dev.off()
}


# All in the same plot ####

# Other plot: all in the same plot
pdf(paste0(folder_name,"/distribution_N_bySim_", "allN", ".pdf"), width = 7, height = 14)

# Initialize plot
dx <- 0.2
plot(0, 0, xlim = xr, ylim = c(1, (nrow(allq) / length(Ns))) + c(0, 2*dx), pch = 1 + includes0 * 15, col = colsN[j], 
     xlab = "Difference between estimated date and real date (median and 95% interval)", ylab = "", axes = FALSE, 
     type = "n")

for(j in seq_along(Ns)){
  theN <- Ns[j]
  par(las = 1)
  sub <- allq[allq$N == theN, ]
  includes0 <- vapply(seq_len(nrow(sub)), function(i){
    ifelse((0 >= sub$x[i, "2.5%"]) & (0 <= sub$x[i, "97.5%"]), 1, 0)
  }, FUN.VALUE = 1)
  
  points(sub$x[, "50%"], sub$Sim + dx * (j-1), xlim = xr, pch = 1 + includes0 * 15, col = colsN[j])
  segments(x0 = sub$x[, "2.5%"], x1 = sub$x[, "97.5%"], 
           y0 = sub$Sim + dx * (j-1), y1 = sub$Sim + dx * (j-1), col = colsN[j])
  #title(main = paste0("N = ", theN))
}
lines(rep(0, 2), c(1, (nrow(allq) / length(Ns)) + 2*dx), lty = 2)
#abline(v = 0, lty = 2)
axis(1, pos = 0.5)
par(xpd = TRUE)
legend("top", inset = c(0, -0.02), col = colsN, legend = paste("N =", names(colsN)), horiz = TRUE, pch = 15)
par(xpd = FALSE)
dev.off()


