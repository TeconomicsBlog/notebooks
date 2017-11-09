# Evan Magnusson

# Clear memory
rm(list=ls())
gc()

# Load Packages
library(simstudy)
library(glmnet)
library(rdd)

# # # # # # # # #
# Simulate Data #
# # # # # # # # #

set.seed(1)
# Number of Observations
N <- 1e3

total.covar <- 50 + 1e3
# Number of covariates (excluding treatment variable)
p <- total.covar - 2

# Simulate Data
mu.vector <- rep(0, total.covar)
variance.vector <- abs(rnorm(total.covar, mean = 1, sd = .5))
#variance.vector <- rep(1, total.covar)
simulated.data <- as.data.frame.matrix(genCorGen(n = N, nvars = total.covar, params1 = mu.vector, params2 = variance.vector, dist = 'normal',  rho = .5,
                            corstr = 'ar1', wide='True'))[2:(total.covar+1)]
colnames(simulated.data)[1] <- 'W'
colnames(simulated.data)[total.covar] <- 'C'


X <- simulated.data[, 2:(total.covar-1)]

covariate.names <- colnames(X)

error <- rnorm(n = N)

# Make W a function of the X's and unobservable
simulated.data$W <- simulated.data$W + .5 * simulated.data$C + 3 * simulated.data$V80 - 6 * simulated.data$V81
# Assign treatment
treated <- (simulated.data$W > 0) * 1.0


# Generate Y as a function of treatment, W, X's, and unobservable C
beta.true.linear <- rnorm(p, mean = 5, sd = 5)
beta.true.linear[30:p] <- 0
Y <- 1.2 * treated - 4 * simulated.data$W  + data.matrix(X) %*% beta.true.linear + .6 * simulated.data$C + error

df <- cbind(Y, simulated.data)
colnames(df)[1] <- 'Y'

# Equations to be used
sumx <- paste(covariate.names, collapse=" + ")
eq.propensity <- paste("W", sumx, sep=" ~ ")
eq.propensity <- as.formula(eq.propensity)

eq.outcome <- paste("Y", sumx, sep=" ~ ")
eq.outcome <- as.formula(eq.outcome)


# LASSO for outcome variables
lasso.fit.outcome <- cv.glmnet(data.matrix(X), df$Y, alpha=1)

coef <- predict(lasso.fit.outcome, type = "nonzero")
colnames <- colnames(X)
H <- colnames[unlist(coef)]
# Vars selected by LASSO:
H

# LASSO for propensity variables
lasso.fit.propensity <- cv.glmnet(data.matrix(X), df$W, alpha=1)

coef <- predict(lasso.fit.propensity, type = "nonzero")
K <- colnames[unlist(coef)]
# Vars selected by LASSO:
K

# Union of selected variables:
doubleselection.names <- unique(c(H, K))
doubleselection.names
sum.doubleselection <- paste(doubleselection.names, collapse = " + ")
eq.doubleselection <- paste("Y ~ W | ", sum.doubleselection)

# RDD, using all covariates selected by double selection
fit <- RDestimate(eq.doubleselection, data = df)
summary(fit)
