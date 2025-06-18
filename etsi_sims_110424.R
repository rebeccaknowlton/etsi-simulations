library(etsi)
library(hetsurr)

# only change setting and run everything below
setting <- 5

# set simulation parameters
n1.a <- 1000 
n0.a <- 1100
n1.b <- 500
n0.b <- 400
kappa <- c(0.5, 0.6, 0.7)
num.sim <- 1000

get.parameters <- function(setting){
  if(setting ==1){
    w.lower <- 0
    w.upper <- 10
    w.cutoff <- 5
    s0.scale <- 2.4  
    s0.shape <- 2.4
    s1.scale <- 2.55
    s1.shape <- 2.55
    beta0 <- 1
    beta1 <- 1.8
    beta2 <- 0
    beta3 <- 0
    beta4 <- 2.8
    beta5 <- 0.1
    sd.y <- 1
    return(list("w.lower" = w.lower, "w.upper" = w.upper, "w.cutoff" = w.cutoff, "s0.scale" = s0.scale, "s0.shape" = s0.shape, "s1.scale" = s1.scale, "s1.shape" = s1.shape, "beta0" = beta0, "beta1" = beta1, "beta2" = beta2, "beta3" = beta3, "beta4" = beta4, "beta5" = beta5, "sd.y" = sd.y))
  }
  if(setting ==2){
    w.lower <- 0
    w.upper <- 10
    w.cutoff.1 <- 2.5
    w.cutoff.2 <- 5
    w.cutoff.3 <- 7.5
    s0.scale <- 2.4  
    s0.shape <- 2.4
    s1.scale <- 2.55
    s1.shape <- 2.55
    beta0 <- c(1, 0.8, 1, 0)
    beta1 <- c(1.8, 0.3, 0.5, 0)
    beta4 <- c(0, .3, 1.5, 1.8)
    beta5 <- c(0, .1, 0.1, 0.05)
    sd.y <- 3
    return(list("w.lower" = w.lower, "w.upper" = w.upper, "w.cutoff.1" = w.cutoff.1, "w.cutoff.2" = w.cutoff.2, "w.cutoff.3" = w.cutoff.3, "s0.scale" = s0.scale, "s0.shape" = s0.shape, "s1.scale" = s1.scale, "s1.shape" = s1.shape, "beta0" = beta0, "beta1" = beta1, "beta4" = beta4, "beta5" = beta5, "sd.y" = sd.y))
  }
  if(setting == 3) {
    w.lower = 0
    w.upper = 12
    s.mean.control = 2
    s.mean.treat = 2
    s.sd.control = 3
    s.sd.treat = 3
    beta0 = 0
    beta1 = 0
    beta2 = 2
    beta3 = 0
    beta4 = 1
    beta5 = 0
    sd.y = 6
    return(list("w.lower" = w.lower, "w.upper" = w.upper, "s.mean.control" = s.mean.control, "s.mean.treat" = s.mean.treat, "s.sd.control" = s.sd.control, "s.sd.treat" = s.sd.treat,"beta0"=beta0,"beta1"=beta1, "beta2"=beta2,"beta3"=beta3, "beta4"=beta4, "beta5"=beta5, "sd.y"=sd.y))
  }
   if(setting ==4){
    w.lower <- 0
    w.upper <- 10
    w.cutoff <- 5
    s0.scale <- 2.4  
    s0.shape <- 2.4
    s1.scale <- 2.55
    s1.shape <- 2.55
    beta0 <- 1
    beta1 <- 1.8
    beta2 <- 0
    beta3 <- 0
    beta4 <- 2.1
    beta5 <- 0.2
    sd.y <- 1
    return(list("w.lower" = w.lower, "w.upper" = w.upper, "w.cutoff" = w.cutoff, "s0.scale" = s0.scale, "s0.shape" = s0.shape, "s1.scale" = s1.scale, "s1.shape" = s1.shape, "beta0" = beta0, "beta1" = beta1, "beta2" = beta2, "beta3" = beta3, "beta4" = beta4, "beta5" = beta5, "sd.y" = sd.y))
  }
  if(setting ==5){
    w.lower <- 0
    w.upper <- 10
    w.cutoff <- 5
    s0.scale <- 2.4  
    s0.shape <- 2.4
    s1.scale <- 2.55
    s1.shape <- 2.55
    beta0 <- 1
    beta1 <- 0.9
    beta2 <- 0
    beta3 <- 0
    beta4 <- 2.8
    beta5 <- 0.1
    sd.y <- 1
    return(list("w.lower" = w.lower, "w.upper" = w.upper, "w.cutoff" = w.cutoff, "s0.scale" = s0.scale, "s0.shape" = s0.shape, "s1.scale" = s1.scale, "s1.shape" = s1.shape, "beta0" = beta0, "beta1" = beta1, "beta2" = beta2, "beta3" = beta3, "beta4" = beta4, "beta5" = beta5, "sd.y" = sd.y))
  }
}

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

# calculate true R and pi
get.truth <- function(setting) {
  params <- get.parameters(setting)
  if (setting == 1 | setting == 4 | setting == 5) {
    W.grid <- runif(100000, params$w.lower, params$w.upper)
    delta.s <- rep(NA, length(grid))
    delta <- rep(NA, length(grid))
    delta.s[W.grid < params$w.cutoff] <- params$beta1 + params$beta3 * params$s0.shape * params$s0.scale 
    delta.s[W.grid >= params$w.cutoff] <- params$beta5 * params$s0.shape * params$s0.scale
    delta[W.grid < params$w.cutoff] <- params$beta1 + (params$beta2 + params$beta3) * (params$s1.shape * params$s1.scale) - 
      (params$beta2) * (params$s0.shape * params$s0.scale)
    delta[W.grid >= params$w.cutoff] <- (params$beta4 + params$beta5) * (params$s1.shape * params$s1.scale) -
      (params$beta4) * (params$s0.shape * params$s0.scale)
    R <- 1-delta.s/delta
    R.mean = mean(R)
    pi <- rep(NA, length(kappa))
    for (k in 1:length(kappa)) {
      pi[k] <- sum(R > kappa[k]) / length(R)
    }
  }
  
  if (setting == 2) {
    W.grid <- runif(100000, params$w.lower, params$w.upper)
    beta0 <- rep(NA, length(W.grid))
    for (i in 1:length(W.grid)) {
      if (W.grid[i] < params$w.cutoff.1) {beta0[i] <- params$beta0[1]}
      else if (W.grid[i] < params$w.cutoff.2) {beta0[i] <- params$beta0[2]}
      else if (W.grid[i] < params$w.cutoff.3) {beta0[i] <- params$beta0[3]}
      else {beta0[i] <- params$beta0[4]}
    }
    beta1 <- rep(NA, length(W.grid))
    for (i in 1:length(W.grid)) {
      if (W.grid[i] < params$w.cutoff.1) {beta1[i] <- params$beta1[1]}
      else if (W.grid[i] < params$w.cutoff.2) {beta1[i] <- params$beta1[2]}
      else if (W.grid[i] < params$w.cutoff.3) {beta1[i] <- params$beta1[3]}
      else {beta1[i] <- params$beta1[4]}
    }
    beta4 <- rep(NA, length(W.grid))
    for (i in 1:length(W.grid)) {
      if (W.grid[i] < params$w.cutoff.1) {beta4[i] <- params$beta4[1]}
      else if (W.grid[i] < params$w.cutoff.2) {beta4[i] <- params$beta4[2]}
      else if (W.grid[i] < params$w.cutoff.3) {beta4[i] <- params$beta4[3]}
      else {beta4[i] <- params$beta4[4]}
    }
    beta5 <- rep(NA, length(W.grid))
    for (i in 1:length(W.grid)) {
      if (W.grid[i] < params$w.cutoff.1) {beta5[i] <- params$beta5[1]}
      else if (W.grid[i] < params$w.cutoff.2) {beta5[i] <- params$beta5[2]}
      else if (W.grid[i] < params$w.cutoff.3) {beta5[i] <- params$beta5[3]}
      else {beta5[i] <- params$beta5[4]}
    }
    
    delta.s <- beta1 + beta5 * params$s0.shape * params$s0.scale
    delta <- beta1 + (beta4 + beta5) * (params$s1.shape * params$s1.scale) - beta4 * (params$s0.shape * params$s0.scale) 
    R <- 1 - delta.s / delta 
    pi <- rep(NA, length(kappa)) 
    for (k in 1:length(kappa)) {
      pi[k] <- sum(R > kappa[k]) / length(R)
    }
  }
  
  if (setting == 3) {
    W.grid <- runif(100000, params$w.lower, params$w.upper)
    delta.s <- params$beta1+params$beta3*params$s.mean.control + W.grid * params$beta5
    delta <- params$beta1+(params$beta2+params$beta3)*(params$s.mean.treat) - (params$beta2)*(params$s.mean.control) + W.grid * params$beta5
    R <- 1-delta.s/delta
    pi <- rep(NA, length(kappa))
    for (k in 1:length(kappa)) {
      pi[k] <- sum(R > kappa[k]) / length(R)
    }
  }
  return(list("W.grid" = W.grid, "R" = R, "pi" = pi))
}

# Generate study A
set.seed(1)
study.A <- gen.data(n1.a, n0.a, setting, "A")
# Export study A
write.table(study.A, paste("etsi.studyA",setting,"_110424.txt",sep=""), quote = FALSE, row.names = FALSE)

control.A <- study.A[study.A$A == 0,]

# calculate R on a grid of W for study A
PTE.results <- hetsurr.fun(y1 = study.A$Y[study.A$A==1],
                           y0 = study.A$Y[study.A$A==0], 
                           s1 = study.A$S[study.A$A==1], 
                           s0 = study.A$S[study.A$A==0], 
                           w1 = study.A$W[study.A$A==1], 
                           w0 = study.A$W[study.A$A==0])

est.delta.B <- function(y1, y0) {
  delta.B <- mean(y1) - mean(y0)
  se.delta.B <- sqrt((1 / length(y1)) * var(y1) + (1 / length(y0)) * var(y0))    
  return(list("delta.B" = delta.B, "se.delta.B" = se.delta.B))
}

est.delta.AB <- function(s1, s0, var.want = TRUE) {
  # Kernel smoothing
  kernel <- function(x, h) {return(dnorm(x / h) / h)}
  get.mu.hat.0 <- function(s, A.s0, A.y0) {
    h.0 <- bw.nrd(A.s0)*length(A.s0)^(-0.2)
    return(sum(kernel(A.s0 - s, h.0) * A.y0 / sum(kernel(A.s0 - s, h.0))))
  }
  
  # Predict y values based on control group from study A
  y1 <- unlist(lapply(s1, get.mu.hat.0, A.s0 = control.A$S, A.y0 = control.A$Y))
  y0 <- unlist(lapply(s0, get.mu.hat.0, A.s0 = control.A$S, A.y0 = control.A$Y))
  if (sum(is.na(y1)) != 0) {
    ind = which(is.na(y1))
    sub.nona = cbind(s1, y1)[-ind,]
    for (i in ind) {
      mat <- cbind(abs(sub.nona[,1] - s1[i]), sub.nona[,2])
      mmm <- which(mat[,1] == min(mat[,1]))[1]
      y1[i] <- sub.nona[mmm,2]
    }
  }
  if (sum(is.na(y0)) != 0) {
    ind = which(is.na(y0))
    sub.nona = cbind(s0, y0)[-ind,]
    for (i in ind) {
      mat <- cbind(abs(sub.nona[,1] - s0[i]), sub.nona[,2])
      mmm <- which(mat[,1] == min(mat[,1]))[1]
      y0[i] <- sub.nona[mmm,2]
    }
  }
  
  delta.AB <- mean(y1) - mean(y0)
  return.list <- list("delta.AB" = delta.AB)
  
  if (var.want) {
    se.delta.AB <- sqrt((1 / length(y1)) * var(y1) + (1 / length(y0)) * var(y0)) 
    return.list <- c(return.list, list("se.delta.AB" = se.delta.AB))
  }
  
  return(return.list)
}


# Function that takes in study B dataframe and a kappa threshold, and returns indicator whether that individual has strong surrogate based on study A
get.delta <- function(df, k) {
  closest.index <- sapply(df$W, function(w) {which.min(abs(PTE.results$w.values - w))})
  return(PTE.results$R.w.s[closest.index] > k)
}

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
    write.table(study.B, paste("etsi.studyB.example",setting,"_110424.txt",sep=""), quote = FALSE, row.names = FALSE)
  }
  
  # delta.B
  results.temp <- est.delta.B(y1 = study.B$Y[study.B$A==1], y0 = study.B$Y[study.B$A==0])
  delta.B.estimates$delta.B[i] <- results.temp$delta.B
  delta.B.estimates$se.delta.B[i] <- results.temp$se.delta.B
  delta.B.estimates$p.value[i] <- 2 * (1 - pnorm(abs(results.temp$delta.B / results.temp$se.delta.B)))
  
  # delta.AB
  results.temp <- est.delta.AB(s1 = study.B$S[study.B$A == 1], s0 = study.B$S[study.B$A == 0], var.want = TRUE)
  delta.AB.estimates$delta.AB[i] <- results.temp$delta.AB
  delta.AB.estimates$se.delta.AB[i] <- results.temp$se.delta.AB
  delta.AB.estimates$p.value[i] <- 2 * (1 - pnorm(abs(results.temp$delta.AB / results.temp$se.delta.AB)))
  
  # delta.P - kappa 1
  study.B$delta <- get.delta(study.B, kappa[1])
  study.A$delta <- get.delta(study.A, kappa[1])
  results.temp <- etsi.main(study.A, study.B)
  delta.P.estimates.k1$delta.P[i] <- results.temp$delta.P
  delta.P.estimates.k1$se.delta.P[i] <- results.temp$se.delta.P
  delta.P.estimates.k1$p.value[i] <- results.temp$p.value
  
  # delta.P - kappa 2
  study.B$delta <- get.delta(study.B, kappa[2])
  study.A$delta <- get.delta(study.A, kappa[2])
  results.temp <- etsi.main(study.A, study.B)
  delta.P.estimates.k2$delta.P[i] <- results.temp$delta.P
  delta.P.estimates.k2$se.delta.P[i] <- results.temp$se.delta.P
  delta.P.estimates.k2$p.value[i] <- results.temp$p.value
  
  # delta.P - kappa 3
  study.B$delta <- get.delta(study.B, kappa[3])
  study.A$delta <- get.delta(study.A, kappa[3])
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

write.table(sim.results, paste("etsi.output",setting, "_110424",".txt",sep=""), quote = FALSE, row.names = FALSE)

