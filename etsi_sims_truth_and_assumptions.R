set.seed(0)
# Calculate true R, true pi, and check assumptions

############################ SETTING 1 ############################
params <- get.parameters(1)
truth.1 <- get.truth(1)
range(truth.1$R)
truth.1$pi

### C1: conditional mean function is monotone increasing in S ###
# this is true because E(Y1) increases as S1 increases, and same for control group 
# (specifically this is true within the strong surrogacy group, not for the weak group)

### C2: P(S1 > s) > P(S0 > s) for all s ###

s.grid <- seq(0, 30, by = 0.01)
surv1 <- 1 - pgamma(s.grid, shape = params$s1.shape, scale = params$s1.scale)
surv0 <- 1 - pgamma(s.grid, shape = params$s0.shape, scale = params$s0.scale)

sum(surv1 > surv0) / length(s.grid)
# holds for almost all values

# survival plot to visually check
plot(s.grid, surv1, type = "l", col = "blue", lwd = 2, 
     ylab = "P(S > s)", xlab = "s", main = "Checking Assumption (C2) (blue should be above red)")
lines(s.grid, surv0, col = "red", lwd = 2)
legend("topright", legend = c("S(1)", "S(0)"), 
       col = c("blue", "red"), lwd = 2)

### C3: E(Y1 |s) > E(Y0 | s) for all s ###
# yes this is true because
params$beta5 > 0

### C4: mu_{A0}(s) = mu_{B0}(s) ###
# this is true because data generating process is the same for A and B

### C5: support of S1 contained within support of S0 ###

pdf1 <- dgamma(s.grid, shape = params$s1.shape, scale = params$s1.scale)
pdf0 <- dgamma(s.grid, shape = params$s0.shape, scale = params$s0.scale)

# plot the pdfs 
plot(s.grid, pdf1, type = "l", col = "blue", lwd = 2, 
     ylab = "Density", xlab = "s", main = "Checking Assumption (C5) (blue should be contained within red)")
lines(s.grid, pdf0, col = "red", lwd = 2)
legend("topright", legend = c("S(1)", "S(0)"), 
       col = c("blue", "red"), lwd = 2)

# yes, blue is contained within red

############################ SETTING 2 ############################
params <- get.parameters(2)
truth.2 <- get.truth(2)
range(truth.2$R)
truth.2$pi

### C1: conditional mean function is monotone increasing in S ###
# this is true because E(Y1) increases as S1 increases, and same for control group 

### C2: P(S1 > s) > P(S0 > s) for all s ###

s.grid <- seq(0, 30, by = 0.01)
surv1 <- 1 - pgamma(s.grid, shape = params$s1.shape, scale = params$s1.scale)
surv0 <- 1 - pgamma(s.grid, shape = params$s0.shape, scale = params$s0.scale)

sum(surv1 > surv0) / length(s.grid)
# holds for nearly all values

# survival plot to visually check
plot(s.grid, surv1, type = "l", col = "blue", lwd = 2, 
     ylab = "P(S > s)", xlab = "s", main = "Checking Assumption (C2) (blue should be above red)")
lines(s.grid, surv0, col = "red", lwd = 2)
legend("topright", legend = c("S(1)", "S(0)"), 
       col = c("blue", "red"), lwd = 2)

### C3: E(Y1 |s) > E(Y0 | s) for all s ###
# true by how we defined the data generating process


### C4: mu_{A0}(s) = mu_{B0}(s) ###
# this is true because data generating process is the same for A and B

### C5: support of S1 contained within support of S0 ###

pdf1 <- dgamma(s.grid, shape = params$s1.shape, scale = params$s1.scale)
pdf0 <- dgamma(s.grid, shape = params$s0.shape, scale = params$s0.scale)

# plot the pdfs 
plot(s.grid, pdf1, type = "l", col = "blue", lwd = 2, 
     ylab = "Density", xlab = "s", main = "Checking Assumption (C5) (blue should be contained within red)")
lines(s.grid, pdf0, col = "red", lwd = 2)
legend("topright", legend = c("S(1)", "S(0)"), 
       col = c("blue", "red"), lwd = 2)

# yes blue is contained within red

############################ SETTING 3 ############################
params <- get.parameters(3)

# remember - truth won't work for the null setting because R undefined
truth.3 <- get.truth(3)
range(truth.3$R)
truth.3$pi

# all assumptions are met by definition since there's no treatment effect and the data generating process is the same in both groups