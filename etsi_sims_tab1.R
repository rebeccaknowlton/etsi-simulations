library(etsi)
library(hetsurr)
library(matrixStats)
library(quantreg)
source("etsi_helper_functions.R")

# set simulation parameters
n1.a <- 1000 
n0.a <- 1100
n1.b <- 500
n0.b <- 400
kappa <- c(0.5, 0.6, 0.7)
num.sim <- 1000

sim_function <- function(setting) {
  
  gen.data <- function(n1, n0, setting, study="A"){
    params <- get.parameters(setting=setting)
    if(setting == 1) {
      data.temp <- data.frame(A = c(rep(1, n1), rep(0, n0)),
                              W = runif(n1 + n0, params$w.lower, params$w.upper))
      
      data.temp$S[data.temp$A==1] <- rgamma(n1, shape = params$s1.shape, scale = params$s1.scale)
      data.temp$S[data.temp$A==0] <- rgamma(n0, shape = params$s0.shape, scale = params$s0.scale)
      data.temp$Y[data.temp$W < params$w.cutoff] <- params$beta0 + params$beta1 * data.temp$A[data.temp$W < params$w.cutoff] + 
        params$beta2 * data.temp$S[data.temp$W < params$w.cutoff] + params$beta3 * data.temp$A[data.temp$W < params$w.cutoff] * data.temp$S[data.temp$W < params$w.cutoff] + rnorm(sum(data.temp$W < params$w.cutoff), 0, params$sd.y)
      data.temp$Y[data.temp$W >= params$w.cutoff] <- params$beta4 * data.temp$S[data.temp$W >= params$w.cutoff] +
        params$beta5 * data.temp$A[data.temp$W >= params$w.cutoff] * data.temp$S[data.temp$W >= params$w.cutoff] + rnorm(sum(data.temp$W >= params$w.cutoff), 0, params$sd.y)
      return(data.temp)
    }
    
    if(setting == 2) {
      data.temp <- data.frame(A = c(rep(1, n1), rep(0, n0)),
                              W = runif(n1 + n0, params$w.lower, params$w.upper))
      data.temp$S[data.temp$A==1] <- rgamma(n1, shape = params$s1.shape, scale = params$s1.scale)
      data.temp$S[data.temp$A==0] <- rgamma(n0, shape = params$s0.shape, scale = params$s0.scale)
      
      #W1
      data.temp$Y[data.temp$W < params$w.cutoff.1] <- params$beta0[1] + params$beta1[1] * data.temp$A[data.temp$W < params$w.cutoff.1] +
        params$beta4[1] * data.temp$S[data.temp$W < params$w.cutoff.1] +
        params$beta5[1] * data.temp$A[data.temp$W < params$w.cutoff.1] * data.temp$S[data.temp$W < params$w.cutoff.1] +
        rnorm(sum(data.temp$W < params$w.cutoff.1), 0, params$sd.y)
      #W2
      data.temp$Y[(data.temp$W >= params$w.cutoff.1) & (data.temp$W < params$w.cutoff.2)] <- params$beta0[2] + params$beta1[2] * data.temp$A[(data.temp$W >= params$w.cutoff.1) & (data.temp$W < params$w.cutoff.2)] +
        params$beta4[2] * data.temp$S[(data.temp$W >= params$w.cutoff.1) & (data.temp$W < params$w.cutoff.2)] +
        params$beta5[2] * data.temp$A[(data.temp$W >= params$w.cutoff.1) & (data.temp$W < params$w.cutoff.2)] * data.temp$S[(data.temp$W >= params$w.cutoff.1) & (data.temp$W < params$w.cutoff.2)] +
        rnorm(sum((data.temp$W >= params$w.cutoff.1) & (data.temp$W < params$w.cutoff.2)), 0, params$sd.y)
      #W3
      data.temp$Y[(data.temp$W >= params$w.cutoff.2) & (data.temp$W < params$w.cutoff.3)] <- params$beta0[3] + params$beta1[3] * data.temp$A[(data.temp$W >= params$w.cutoff.2) & (data.temp$W < params$w.cutoff.3)] +
        params$beta4[3] * data.temp$S[(data.temp$W >= params$w.cutoff.2) & (data.temp$W < params$w.cutoff.3)] +
        params$beta5[3] * data.temp$A[(data.temp$W >= params$w.cutoff.2) & (data.temp$W < params$w.cutoff.3)] * data.temp$S[(data.temp$W >= params$w.cutoff.2) & (data.temp$W < params$w.cutoff.3)] +
        rnorm(sum((data.temp$W >= params$w.cutoff.2) & (data.temp$W < params$w.cutoff.3)), 0, params$sd.y)
      #W4
      data.temp$Y[data.temp$W >= params$w.cutoff.3] <- params$beta0[4] + params$beta1[4] * data.temp$A[data.temp$W >= params$w.cutoff.3] +
        params$beta4[4] * data.temp$S[data.temp$W >= params$w.cutoff.3] +
        params$beta5[4] * data.temp$A[data.temp$W >= params$w.cutoff.3] * data.temp$S[data.temp$W >= params$w.cutoff.3] +
        rnorm(sum(data.temp$W >= params$w.cutoff.3), 0, params$sd.y)
      
      return(data.temp)
    }
    
    if(setting == 3) {
      if (study == "A") {
        data.temp <- data.frame(A = c(rep(1,n1),rep(0,n0)),
                                W = runif(n1 + n0, params$w.lower, params$w.upper))  
      } else if (study == "B") {
        data.temp <- data.frame(A = c(rep(1,n1),rep(0,n0)),
                                W = runif(n1 + n0, params$w.lower, params$w.upper))
      }
      data.temp$S[data.temp$A==1] <- rnorm(n1, mean = params$s.mean.treat, sd = params$s.sd.treat)
      data.temp$S[data.temp$A==0] <- rnorm(n0, mean = params$s.mean.control, sd = params$s.sd.control)
      data.temp$Y <- params$beta2 * data.temp$S + params$beta4 * data.temp$W + params$beta5 * data.temp$A * data.temp$W + rnorm(n1+n0, mean = 0, sd = params$sd.y)
      return(data.temp)
    }
    if(setting == 4) {
    	if (study == "A") {
    		params <- get.parameters(setting=1)
      data.temp <- data.frame(A = c(rep(1, n1), rep(0, n0)),
                              W = runif(n1 + n0, params$w.lower, params$w.upper))
      
      data.temp$S[data.temp$A==1] <- rgamma(n1, shape = params$s1.shape, scale = params$s1.scale)
      data.temp$S[data.temp$A==0] <- rgamma(n0, shape = params$s0.shape, scale = params$s0.scale)
      data.temp$Y[data.temp$W < params$w.cutoff] <- params$beta0 + params$beta1 * data.temp$A[data.temp$W < params$w.cutoff] + 
        params$beta2 * data.temp$S[data.temp$W < params$w.cutoff] + params$beta3 * data.temp$A[data.temp$W < params$w.cutoff] * data.temp$S[data.temp$W < params$w.cutoff] + rnorm(sum(data.temp$W < params$w.cutoff), 0, params$sd.y)
      data.temp$Y[data.temp$W >= params$w.cutoff] <- params$beta4 * data.temp$S[data.temp$W >= params$w.cutoff] +
        params$beta5 * data.temp$A[data.temp$W >= params$w.cutoff] * data.temp$S[data.temp$W >= params$w.cutoff] + rnorm(sum(data.temp$W >= params$w.cutoff), 0, params$sd.y)}
      if (study == "B") {
    		params <- get.parameters(setting=4)
      data.temp <- data.frame(A = c(rep(1, n1), rep(0, n0)),
                              W = runif(n1 + n0, params$w.lower, params$w.upper))
      
      data.temp$S[data.temp$A==1] <- rgamma(n1, shape = params$s1.shape, scale = params$s1.scale)
      data.temp$S[data.temp$A==0] <- rgamma(n0, shape = params$s0.shape, scale = params$s0.scale)
      data.temp$Y[data.temp$W < params$w.cutoff] <- params$beta0 + params$beta1 * data.temp$A[data.temp$W < params$w.cutoff] + 
        params$beta2 * data.temp$S[data.temp$W < params$w.cutoff] + params$beta3 * data.temp$A[data.temp$W < params$w.cutoff] * data.temp$S[data.temp$W < params$w.cutoff] + rnorm(sum(data.temp$W < params$w.cutoff), 0, params$sd.y)
      data.temp$Y[data.temp$W >= params$w.cutoff] <- params$beta4 * data.temp$S[data.temp$W >= params$w.cutoff] +
        params$beta5 * data.temp$A[data.temp$W >= params$w.cutoff] * data.temp$S[data.temp$W >= params$w.cutoff] + rnorm(sum(data.temp$W >= params$w.cutoff), 0, params$sd.y)}
      return(data.temp)
    }
    if(setting == 5) {
    	if (study == "A") {
    		params <- get.parameters(setting=1)
      data.temp <- data.frame(A = c(rep(1, n1), rep(0, n0)),
                              W = runif(n1 + n0, params$w.lower, params$w.upper))
      
      data.temp$S[data.temp$A==1] <- rgamma(n1, shape = params$s1.shape, scale = params$s1.scale)
      data.temp$S[data.temp$A==0] <- rgamma(n0, shape = params$s0.shape, scale = params$s0.scale)
      data.temp$Y[data.temp$W < params$w.cutoff] <- params$beta0 + params$beta1 * data.temp$A[data.temp$W < params$w.cutoff] + 
        params$beta2 * data.temp$S[data.temp$W < params$w.cutoff] + params$beta3 * data.temp$A[data.temp$W < params$w.cutoff] * data.temp$S[data.temp$W < params$w.cutoff] + rnorm(sum(data.temp$W < params$w.cutoff), 0, params$sd.y)
      data.temp$Y[data.temp$W >= params$w.cutoff] <- params$beta4 * data.temp$S[data.temp$W >= params$w.cutoff] +
        params$beta5 * data.temp$A[data.temp$W >= params$w.cutoff] * data.temp$S[data.temp$W >= params$w.cutoff] + rnorm(sum(data.temp$W >= params$w.cutoff), 0, params$sd.y)}
      if (study == "B") {
    		params <- get.parameters(setting=5)
      data.temp <- data.frame(A = c(rep(1, n1), rep(0, n0)),
                              W = runif(n1 + n0, params$w.lower, params$w.upper))
      
      data.temp$S[data.temp$A==1] <- rgamma(n1, shape = params$s1.shape, scale = params$s1.scale)
      data.temp$S[data.temp$A==0] <- rgamma(n0, shape = params$s0.shape, scale = params$s0.scale)
      
      data.temp$Y<- params$beta0 + params$beta1 * data.temp$A + 
        params$beta2 * data.temp$S + params$beta3 * data.temp$A * data.temp$S+ rnorm(length(data.temp$A), 0, params$sd.y)}
      return(data.temp)
    }
  }
  
  # Generate study A
  set.seed(1)
  study.A <- gen.data(n1.a, n0.a, setting, "A")
  # Export study A
  write.table(study.A, paste("sim_output/etsi.studyA_setting",setting,".txt",sep=""), quote = FALSE, row.names = FALSE)
  
  control.A <- study.A[study.A$A == 0,]
  
  # calculate R on a grid of W for study A
  PTE.results <- hetsurr.fun(y1 = study.A$Y[study.A$A==1],
                             y0 = study.A$Y[study.A$A==0], 
                             s1 = study.A$S[study.A$A==1], 
                             s0 = study.A$S[study.A$A==0], 
                             w1 = study.A$W[study.A$A==1], 
                             w0 = study.A$W[study.A$A==0])
  
  
  # Generate study B with same data generation setup, and test H_0: \delta_B = 0
  
  delta.B.estimates <- data.frame("delta.B" = rep(NA, num.sim),
                                  "se.delta.B" = rep(NA, num.sim),
                                  "p.value" = rep(NA, num.sim))
  delta.AB.estimates <- data.frame("delta.AB" = rep(NA, num.sim),
                                   "se.delta.AB" = rep(NA, num.sim),
                                   "p.value" = rep(NA, num.sim))
  delta.P.estimates.k1 <- data.frame("delta.P" = rep(NA, num.sim),
                                     "se.delta.P" = rep(NA, num.sim),
                                     "p.value" = rep(NA, num.sim))
  delta.P.estimates.k2 <- data.frame("delta.P" = rep(NA, num.sim),
                                     "se.delta.P" = rep(NA, num.sim),
                                     "p.value" = rep(NA, num.sim))
  delta.P.estimates.k3 <- data.frame("delta.P" = rep(NA, num.sim),
                                     "se.delta.P" = rep(NA, num.sim),
                                     "p.value" = rep(NA, num.sim))
  
  for (i in 1:num.sim) {
    study.B <- gen.data(n1.b, n0.b, setting, "B")
    
    # Export one iteration for example data, setting 1 only
    if (i == 1 & setting == 1) {
      write.table(study.B, paste("sim_output/etsi.studyB.example_setting",setting,".txt",sep=""), quote = FALSE, row.names = FALSE)
    }
    
    # delta.B
    results.temp <- est.delta.B(y1 = study.B$Y[study.B$A==1], y0 = study.B$Y[study.B$A==0])
    delta.B.estimates$delta.B[i] <- results.temp$delta.B
    delta.B.estimates$se.delta.B[i] <- results.temp$se.delta.B
    delta.B.estimates$p.value[i] <- 2 * (1 - pnorm(abs(results.temp$delta.B / results.temp$se.delta.B)))
    
    # delta.AB
    results.temp <- est.delta.AB(s1 = study.B$S[study.B$A == 1], s0 = study.B$S[study.B$A == 0], var.want = TRUE, control.A)
    delta.AB.estimates$delta.AB[i] <- results.temp$delta.AB
    delta.AB.estimates$se.delta.AB[i] <- results.temp$se.delta.AB
    delta.AB.estimates$p.value[i] <- 2 * (1 - pnorm(abs(results.temp$delta.AB / results.temp$se.delta.AB)))
    
    # delta.P - kappa 1
    study.B$delta <- get.delta(study.B, kappa[1], PTE.results)
    study.A$delta <- get.delta(study.A, kappa[1], PTE.results)
    results.temp <- etsi.main(study.A, study.B)
    delta.P.estimates.k1$delta.P[i] <- results.temp$delta.P
    delta.P.estimates.k1$se.delta.P[i] <- results.temp$se.delta.P
    delta.P.estimates.k1$p.value[i] <- results.temp$p.value
    
    # delta.P - kappa 2
    study.B$delta <- get.delta(study.B, kappa[2], PTE.results)
    study.A$delta <- get.delta(study.A, kappa[2], PTE.results)
    results.temp <- etsi.main(study.A, study.B)
    delta.P.estimates.k2$delta.P[i] <- results.temp$delta.P
    delta.P.estimates.k2$se.delta.P[i] <- results.temp$se.delta.P
    delta.P.estimates.k2$p.value[i] <- results.temp$p.value
    
    # delta.P - kappa 3
    study.B$delta <- get.delta(study.B, kappa[3], PTE.results)
    study.A$delta <- get.delta(study.A, kappa[3], PTE.results)
    results.temp <- etsi.main(study.A, study.B)
    delta.P.estimates.k3$delta.P[i] <- results.temp$delta.P
    delta.P.estimates.k3$se.delta.P[i] <- results.temp$se.delta.P
    delta.P.estimates.k3$p.value[i] <- results.temp$p.value
    
    print(i)
  }
  
  
  sim.results <- cbind(delta.B.estimates$delta.B, delta.B.estimates$se.delta.B, delta.B.estimates$p.value,
                       delta.AB.estimates$delta.AB, delta.AB.estimates$se.delta.AB, delta.AB.estimates$p.value,
                       delta.P.estimates.k1$delta.P, delta.P.estimates.k1$se.delta.P, delta.P.estimates.k1$p.value,
                       delta.P.estimates.k2$delta.P, delta.P.estimates.k2$se.delta.P, delta.P.estimates.k2$p.value,
                       delta.P.estimates.k3$delta.P, delta.P.estimates.k3$se.delta.P, delta.P.estimates.k3$p.value)
  colnames(sim.results) <- c("delta.B", "se.delta.B", "p.value.B",
                             "delta.AB", "se.delta.AB", "p.value.AB",
                             "delta.P.k1", "se.delta.P.k1", "p.value.P.k1",
                             "delta.P.k2", "se.delta.P.k2", "p.value.P.k2",
                             "delta.P.k3", "se.delta.P.k3", "p.value.P.k3")
  write.table(sim.results, paste("sim_output/etsi.output_setting",setting,".txt",sep=""), quote = FALSE, row.names = FALSE)
  
  results <- matrix(ncol = 5, nrow = 6)
  if (setting == 1 | setting == 2 | setting == 4 |setting == 5) {
    rownames(results) <- c("Point Estimate", "ESE", "ASE", "Effect Size", "Power", "pi_B")  
  } else if (setting == 3) {
    rownames(results) <- c("Point Estimate", "ESE", "ASE", "Effect Size", "Type 1 Error", "pi_B")
  }
  colnames(results) <- c("Delta.B", "Delta.AB", "Delta.P: k = 0.5", "Delta.P: k = 0.6", "Delta.P: k = 0.7")
  
  
  results[1,] <- colMeans(sim.results[,c(1,4,7,10,13)])
  results[2,] <- sqrt(colVars(as.matrix(sim.results[,c(1,4,7,10,13)])))
  results[3,] <- colMeans(sim.results[,c(2,5,8,11,14)])
  results[4,] <- colMeans(sim.results[,c(1,4,7,10,13)] / sim.results[,c(2,5,8,11,14)])
  results[5,] <- apply(sim.results[,c(3,6,9,12,15)], 2, function(x) mean(x < 0.05))
  set.seed(1)
  results[6,] <- c(0, 1, get.truth(setting)$pi)
  
  print(round(results,3))
  return(results)
}

tab1_setting1 <- sim_function(1)
tab1_setting2 <- sim_function(2)
tab1_setting3 <- sim_function(3)

tab1 <- rbind(tab1_setting1, tab1_setting2, tab1_setting3)
saveRDS(tab1, file = "results/tab1.rds")

