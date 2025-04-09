###########################################################################################
#
#         Design 1: Case-Control
#         See Lammertse, Brunner, van Berkel et al. (2022)
#         For information on SIMR package, see:
#         https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12504
#
###########################################################################################

################# User input #################

### Set working directory to folder where results should be saved (change to your own directory)
setwd("")

# Specify filename for output
filename <- "ipsc_results"

# specify effect size (cohen's d) to be modeled
x <- 0.2

# specify ICC based on own pilot data
ICC <- 0.2

# Determine total number of independent iPSC-lines (i.e. 2 cases vs. 2 controls = 4 lines)
lines <- 30

# Set significance level for the statistical test
alpha <- 0.05

################## Running the power simulation ##################

### Load libraries
library(tidyverse)
library(lme4)
library(simr)

# Create simulated dataset
# Start with 2 lines per condition and 100 observations per line
# Randomly sample standardized data from a normal distribution with a mean of 0 and standard deviation of 1
set.seed(1234)
X <- data.frame(
  Condition = rep(0:1, 200),
  Line = rep(1:4, 100),
  y = rnorm(400, mean = 0, sd = 1)
)

# Create a multilevel model
mod1 <- lmer(y ~ Condition + (1 | Line), data = X)

# Set residual variance
s <- 1 # Standardized
# Formula to calculate ICC: ICC = b/(b+s), therefore:
# Formula to get b for a given ICC and a given s
bcalc <- (ICC * s) / (1 - ICC)
# specify bcalc as random variance of line
VarCorr(mod1)["Line"] <- bcalc
# specify st.Dev of residual variance
sigma(mod1) <- s # Standardized
# Set effect size
fixef(mod1)["Condition"] <- x

# Increase number of lines to user settings
mod2 <- extend(mod1, along = "Line", n = lines)

# Statistical test is: F-test and degrees of freedom based on Kenward-Roger approximation "kr" (default)
# Significance threshold: alpha = 0.5 (default)
# Create powerCurve with 1000 simulations (default)
PC1 <- powerCurve(mod2, within = "Line+Condition", seed = 1234, alpha = alpha, nsim = 1000, fixed("Condition", "kr"), progress = F)

# Create dataframe of powercurve and save
PC1df <- as.data.frame(summary(PC1))
csv <- paste0(filename, ".csv")
write.csv(PC1df, file = csv, row.names = F)

# Save original plot
pdf <- paste0(filename, ".pdf")
pdf(file = pdf)
plot(PC1)
dev.off()

# Save PowerCurve
rdata <- paste0(filename, ".Rdata")
save(PC1, file = rdata)
