library(etsi)
library(quantreg)

# only change setting and run everything below
# also need to load get.truth, get.parameters from the master sims file
setting <- 2

# read in study A from simulations
setwd("C:/Users/rkkno/Documents/University of Texas at Austin/etsi")
study.A <- read.table(paste("etsi.studyA",setting,"_110424.txt", sep=""), header = T)

# empirical power
output <- read.table(paste("etsi.output", setting, "_110424",".txt", sep=""), header = T)
emp.power <- apply(output[,c(9,12,15)], 2, function(x) mean(x < 0.05))

thresholds <- c(0.5, 0.6, 0.7)

power.sim <- rep(NA, length(thresholds))
n.emp.sim <- rep(NA, length(thresholds))
n.est.sim <- rep(NA, length(thresholds))

for (k in 1:length(thresholds)) {
  set.seed(1) # Should start with the same seed each time we call etsi.design, so train/test splits are the same
  power.sim[k] <- etsi.design(Study.A = study.A, n.b0 = 400, n.b1 = 500, kappa = thresholds[k])
  set.seed(1)
  n.emp.sim[k] <- etsi.design(Study.A = study.A, kappa = thresholds[k], desired.power = emp.power[k])
  set.seed(1)
  n.est.sim[k] <- etsi.design(Study.A = study.A, kappa = thresholds[k], desired.power = unlist(power.sim[k]))
}

# results 
results <- matrix(ncol = 3, nrow = 3)
rownames(results) <- c("Estimated Power", "Empirical Power","pi_B")  
colnames(results) <- c("Delta.P: k = 0.5", "Delta.P: k = 0.6", "Delta.P: k = 0.7")

results[1,] <- unlist(power.sim)
results[2,] <- emp.power
set.seed(1)
results[3,] <- get.truth(setting)$pi
results.latex <- format(round(results ,3),nsmall=3)
latex.table(results.latex, paste0("designsimres", setting), caption = "", dcolumn = T)
print(round(results,3))

# sample size results
n.results <- matrix(ncol = 3, nrow = 2)
rownames(n.results) <- c("Estimated n (from empirical power)", "Estimated n (from estimated power)")  
colnames(n.results) <- c("Delta.P: k = 0.5", "Delta.P: k = 0.6", "Delta.P: k = 0.7")
n.results[1,] <- unlist(n.emp.sim)
n.results[2,] <- unlist(n.est.sim)
n.results.latex <- format(round(n.results ,3),nsmall=3)
latex.table(n.results.latex, paste0("n_designsimres", setting), caption = "", dcolumn = T)
print(round(n.results,3))
