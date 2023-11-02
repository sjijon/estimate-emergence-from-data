## Florence Débare
## iEES 
##
## ----------------------------------------------------------------------------------
## Reproducing and extending the analysis presented in Roberts et al. (2021):
##  Roberts DL, Rossman JS, Jarić I (2021) 
##  Dating first cases of COVID-19. PLoS Pathog 17(6): e1009620. 
##  https://doi.org/10.1371/journal.ppat.1009620
##
## -> Take-Home message: 
## The estimated dates of origin change a lot with changes in the data. 
## With corrected datasets, the estimated origin is much later. 
## ----------------------------------------------------------------------------------

## Data ####
# dates at which there is at least one onset date

# Robert et al. used 10 case days (cf. their Table S2)

# Date used by Roberts et al. 
datesHuang <- c("2019-12-01", "2019-12-10", "2019-12-15", "2019-12-17", "2019-12-18", 
                "2019-12-19", "2019-12-20", "2019-12-21", "2019-12-22", "2019-12-23")

# WHO 2021 data, after Worobey correction
datesWHOWorobey <- c("2019-12-10", "2019-12-11", "2019-12-12", "2019-12-13", "2019-12-15", 
                     "2019-12-16", "2019-12-17", "2019-12-18", "2019-12-19", "2019-12-20")

# WHO 2021 data
datesWHOraw <- c("2019-12-08", "2019-12-11", "2019-12-12", "2019-12-13", "2019-12-15", 
                 "2019-12-16", "2019-12-17", "2019-12-18", "2019-12-19", "2019-12-20")

# Other way to write the Huang et al. dates, for checks
tt <- rev(c(1, 10, 15, 17, 18, 19, 20, 21, 22, 23))
# Needs to be ordered from most recent to latest date
k <- length(tt) # Number of dates

## -------------------------------------------------------------------

# Functions ####

OLE.CC <- function(sights, alpha){
  # sights: vector needs to be ordered first
  # alpha: significance level, 0.1 for one-sided
  
  # This is part of the OLE.func function 
  # taken from https://github.com/cran/sExtinct
  k <- length(sights) # Number of observations
  # Compute nu
  v <- (1/(k-1)) * sum(log((sights[1] - sights[k])/(sights[1] - sights[2:(k-1)])))
  e <- matrix(rep(1, k), ncol=1)
  # 1/c(alpha)
  SL <- (-log(1-alpha/2)/length(sights))^-v
  SU <- (-log(alpha/2)/length(sights))^-v
  # Compute lambda matrix
  myfun <- function(i,j,v){(gamma(2*v+i)*gamma(v+j))/(gamma(v+i)*gamma(j))}
  lambda <- outer(1:k, 1:k, myfun, v=v)
  lambda <- ifelse(lower.tri(lambda), lambda, t(lambda)) # lambda is symmetrical
  # Compute a
  a <- as.vector(solve(t(e)%*%solve(lambda)%*%e)) * solve(lambda)%*%e
  # Compute CI
  lowerCI <- max(sights) + ((max(sights)-min(sights))/(SL-1))
  upperCI <- max(sights) + ((max(sights)-min(sights))/(SU-1))
  # Estimated date
  extest <- sum(t(a)%*%sights)
  res <- data.frame(Estimate=extest, lowerCI=lowerCI, upperCI=upperCI)
  res
}

# Function for OLE on emergence dates
OLE.emergence <- function(sights, endDate = as.Date("2019-12-31"), alpha){
  # sights: dates of the sights, YYYY-MM-DD
  # endDate: final date from which dates are calculated
  # alpha: significance level (multiply by 2 for one-sided)
  
  estim <- OLE.CC(rev(sort(as.numeric(endDate - as.Date(sights)))), alpha)
  c(Estimate = endDate - estim[1, 1], 
             lowerCI = endDate - estim[1, 2], 
             upperCI = endDate - estim[1, 3])
}

# Checks
if(FALSE){
  # No need to evaluate this, these checks are just kept in case one wants to see them
  OLE.CC(sort(tt), 0.1)
  OLE.CC(rev(sort(tt)), 0.1)
  OLE.CC(sort(tt-max(tt)), 0.1)
  OLE.CC(rev(sort(max(tt)-tt)), 0.1) # -> correct one
  endDate <- as.Date("2019-12-31") 
  estim <- OLE.CC(rev(sort(as.numeric(endDate - as.Date(datesHuang)))), 0.1)
  endDate - estim[1, 1]
  endDate - estim[1, 3]
}

## -------------------------------------------------------------------

# Analysis ####

# Compute emergence dates
# With Huang et al. data
# alpha=0.1 because one-sided (they only reported earliest date)
resHuang <- OLE.emergence(datesHuang, alpha = 0.1)
resHuang
# >     Estimate    lowerCI    upperCI
# > 1 2019-11-17 2019-11-23 2019-10-04
# -> same dates as given in the paper

# Now with the other datasets
resWHOWorobey <- OLE.emergence(datesWHOWorobey, alpha = 0.1)
resWHOWorobey
# >     Estimate    lowerCI    upperCI
# > 1 2019-12-08 2019-12-09 2019-12-03

resWHOraw <- OLE.emergence(datesWHOraw, alpha = 0.1)
resWHOraw
# >     Estimate    lowerCI    upperCI
# > 1 2019-12-04 2019-12-07 2019-11-24

# Plots ####

alldates2 <- sort(as.Date(unique(c(as.Date(c(datesHuang, datesWHOWorobey, datesWHOraw)), resHuang, resWHOWorobey, resWHOraw))))
fname <- "plotRoberts.pdf"

# Initialize plot
pdf(file = fname, width = 8, height = 3)
par(las = 1, mai = c(0.85, 1.65, 0.65, 0.1))

ypos <- rev(seq(0.5, 2, length.out = 3)) # Positionts of the different datasets/results
xvals <- c(alldates2, as.Date("2019-10-01"), max(alldates2)+1) # xvalues (to have rounder numbers)
plot(x = xvals, y = rep(0, length(xvals)), type = "n", axes = FALSE, 
     ylim = c(0, 2.5), 
     xlab = "date (2019)", ylab = "")
# x axis
tcks <- -0.05 # x axis tick length
xx <- seq(min(c(alldates2), as.Date("2019-10-01")), max(alldates2), by = "day") # x axis tick positions
axis(1, at = xx, labels = rep("", length(xx)), tck = tcks, pos = 0)
xl <- c("2019-10-01", "2019-10-15", "2019-11-01", "2019-11-15", "2019-12-01", "2019-12-15") # x axis labels positions
axis(1, at = as.Date(xl), labels = rep("", length(xl)), tck = tcks - 0.02, pos = 0, lwd = 0, lwd.ticks = 2, lend = "butt") # ticks for the labels
par(xpd = TRUE)
text(x = as.Date(xl), y = 0, labels = paste0(format(as.Date(xl), "%b %d"), "   "), srt = 90, cex = 0.9, adj = c(1, 0.5)) # labels

# Add results
suffix <- c("Huang", "WHOraw", "WHOWorobey") # suffixes of the vectors we will use
nms <- c("Huang et al. (2020)*", "WHO (2021)", "WHO (2021), corrected") # Corresponding names to print

# Plotting parameters
colDates <- gray(0.4) # color case dates
colEstimates <- "#E69F00" # color estimates
pchDates <- 16
pchEstimates <- 15
lwdCI <- 2

for(i in seq_along(suffix)){
  # Case dates
  x <- get(paste0("dates", suffix[i]))
  points(as.Date(x), rep(ypos[i], length(x)), pch = pchDates, cex = 1, col = colDates)
  # Estimates
  estim <- get(paste0("res", suffix[i]))
  points(as.Date(estim[1]), ypos[i], col = colEstimates, pch = pchEstimates, cex = 1.2) # estimate
  lines(as.Date(estim[2:3]), rep(ypos[i], 2), col = colEstimates, lwd = lwdCI, lend = "butt") # CI
  # Text label for the dataset
  text(x = min(as.Date(xvals)), y = ypos[i], adj = c(1, 0.5), labels = nms[i])
}

# Add legend
yleg <- max(ypos) + 0.5 # y position of the legend
xleg <- min(xvals) - 25
legend(x = xleg, y = yleg, 
       yjust = 0, xjust = 0,
       legend = c("case dates", "estimated 1st infection", "95% confidence interval"), 
       col = c(colDates, colEstimates, colEstimates), 
       pch = c(pchDates, pchEstimates, NA), 
       lwd = c(0, 0, lwdCI), bty = "n", lty = c(0, 0, 1), 
       horiz = FALSE)
# Add title for data names
text(x = xleg, y = yleg, adj = c(0, 0.51), labels = "Dataset         ")
dev.off()
system(paste("open", fname))
