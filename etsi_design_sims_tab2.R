library(etsi)
library(quantreg)

design_function <- function(setting) {
  # read in study A from simulations
  study.A <- read.table(paste("sim_output/etsi.studyA_setting",setting,".txt", sep=""), header = T)
  
  # empirical power
  output <- read.table(paste("sim_output/etsi.output_setting", setting,".txt", sep=""), header = T)
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
  
  # sample size results
  n.results <- matrix(ncol = 3, nrow = 2)
  rownames(n.results) <- c("Estimated n (from empirical power)", "Estimated n (from estimated power)")  
  colnames(n.results) <- c("Delta.P: k = 0.5", "Delta.P: k = 0.6", "Delta.P: k = 0.7")
  n.results[1,] <- unlist(n.emp.sim)
  n.results[2,] <- unlist(n.est.sim)

  return(list("results" = results,
           "n.results" = n.results))
}


tab2_setting1 <- design_function(1)$results
tab2_setting2 <- design_function(2)$results

tab2 <- rbind(tab2_setting1, tab2_setting2)
saveRDS(tab2, file = "results/tab2.rds")
