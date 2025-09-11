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

est.delta.B <- function(y1, y0) {
  delta.B <- mean(y1) - mean(y0)
  se.delta.B <- sqrt((1 / length(y1)) * var(y1) + (1 / length(y0)) * var(y0))    
  return(list("delta.B" = delta.B, "se.delta.B" = se.delta.B))
}

est.delta.AB <- function(s1, s0, var.want = TRUE, control.A) {
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
get.delta <- function(df, k, PTE.results) {
  closest.index <- sapply(df$W, function(w) {which.min(abs(PTE.results$w.values - w))})
  return(PTE.results$R.w.s[closest.index] > k)
}