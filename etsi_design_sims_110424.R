library(etsi)
library(quantreg)

# only change setting and run everything below
# also need to load est.delta.B, est.delta.AB, get.truth, get.parameters from the master sims file
setting <- 2

# read in study A from simulations
setwd("C:/Users/rkkno/Documents/University of Texas at Austin/etsi")
study.A <- read.table(paste("etsi.studyA",setting,"_110424.txt", sep=""), header = T)

# empirical power
output <- read.table(paste("etsi.output", setting, "_110424",".txt", sep=""), header = T)
emp.power <- apply(output[,c(3,6,9,12,15)], 2, function(x) mean(x < 0.05))

kappa <- c(0.5, 0.6, 0.7)
thresholds <- c(0,1, kappa)

power.sim <- rep(NA, length(thresholds))

# only want to compare versions for n with the pooled versions
n.sim <- rep(NA, length(thresholds) - 2)

for (k in 1:length(thresholds)) {
  set.seed(1) # should start with the same seed for each threshold
  if (k == 1) {
    iterations <- 100
    holdout <- 0.5
    power.vec <- rep(NA, iterations)
    for (i in 1:iterations) {
      train <- sample(1:nrow(study.A), round(nrow(study.A)*holdout))
      study.A.train <- study.A[train,]
      study.A.test <- study.A[-train,]
      
      # treat everyone as weak surrogate
      delta.B <- est.delta.B(y1 = study.A.test$Y[study.A.test$A==1],
                             y0 = study.A.test$Y[study.A.test$A==0])
      power.vec[i] <- 1 - pnorm(1.96 - delta.B$delta.B / delta.B$se.delta.B)
    }
    power.sim[k] <- mean(power.vec)
  } else if (k == 2) {
    iterations <- 100
    holdout <- 0.5
    power.vec <- rep(NA, iterations)
    for (i in 1:iterations) {
      train <- sample(1:nrow(study.A), round(nrow(study.A)*holdout))
      study.A.train <- study.A[train,]
      study.A.test <- study.A[-train,]
      
      # treat everyone as strong surrogate
      control.A = study.A.train
      delta.AB <- est.delta.AB(s1 = study.A.test$S[study.A.test$A==1],
                               s0 = study.A.test$S[study.A.test$A==0])
      power.vec[i] <- 1 - pnorm(1.96 - delta.AB$delta.AB / delta.AB$se.delta.AB)
    }
    power.sim[k] <- mean(power.vec)
  } else {
    power.sim[k] <- etsi.design(Study.A = study.A, n.b0 = 400, n.b1 = 500, kappa = thresholds[k])
    n.sim[k-2] <- etsi.design(Study.A = study.A, kappa = thresholds[k], desired.power = emp.power[k])
  }
}

# results 
results <- matrix(ncol = 5, nrow = 3)
rownames(results) <- c("Estimated Power", "Empirical Power","pi_B")  
colnames(results) <- c("Delta.B", "Delta.AB", "Delta.P: k > 0.5", "Delta.P: k > 0.6", "Delta.P: k > 0.7")

results[1,] <- unlist(power.sim)
results[2,] <- emp.power
set.seed(1)
results[3,] <- c(0, 1, get.truth(setting)$pi)
results.latex <- format(round(results ,3),nsmall=3)
latex.table(results.latex, paste0("designsimres", setting), caption = "", dcolumn = T)
print(round(results,3))

# sample size for pooled estimated only
n.results <- matrix(ncol = 3, nrow = 1)
rownames(n.results) <- c("Estimated n")  
colnames(n.results) <- c("Delta.P: k > 0.5", "Delta.P: k > 0.6", "Delta.P: k > 0.7")
n.results[1,] <- unlist(n.sim)
n.results.latex <- format(round(n.results ,3),nsmall=3)
latex.table(n.results.latex, paste0("n_designsimres", setting), caption = "", dcolumn = T)
print(round(n.results,3))
