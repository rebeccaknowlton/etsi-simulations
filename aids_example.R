library(etsi)
library(quantreg)

############################################
#### USING ACTG 320 AS STUDY A
################################################

# Set working directory to data location
setwd("C:/Users/rkkno/Documents/University of Texas at Austin/etsi/AIDS Data/ACTG320-Data")
tmpg <- function(xx){exp(xx)}; tmpginv <- function(xx){log(xx)}
tmpg <- function(xx){xx^2}; tmpdg <- function(xx){2*xx}; tmpginv <- function(xx){sqrt(xx)}
log.link <- list(tmpg, tmpg, tmpginv); rm(tmpg, tmpginv)
tmpdir <- "./"
aids.base <- read.table(paste(tmpdir, "actg320.base.dat", sep=""))
aids.rna <- read.table(paste(tmpdir, "actg320.rna", sep=""))
aids.cd4 <- read.table(paste(tmpdir, "actg320.cd4", sep=""))
aids.trt <- read.table(paste(tmpdir, "actg320.trt", sep=""))

n.total = dim(aids.base)[1]
aids.base$TREAT = vector(length = n.total)
aids.base$CD4BASE = vector(length = n.total)
aids.base$CD424 = vector(length = n.total)
aids.base$RNABASE = vector(length = n.total)
aids.base$RNA24 = vector(length = n.total)

for(i in 1:n.total) {
	aids.base$TREAT[i] = aids.trt[aids.trt$PATID == aids.base$PATID[i],2]
	if(dim(aids.cd4[aids.cd4$PATID == aids.base$PATID[i] & aids.cd4$WEEK == 0,])[1] > 0) {aids.base$CD4BASE[i] = aids.cd4[aids.cd4$PATID == aids.base$PATID[i] & aids.cd4$WEEK == 0,]$CD4} else {aids.base$CD4BASE[i] = NA}
	if(dim(aids.cd4[aids.cd4$PATID == aids.base$PATID[i] & aids.cd4$WEEK == 24,])[1] > 0) {aids.base$CD424[i] = aids.cd4[aids.cd4$PATID == aids.base$PATID[i] & aids.cd4$WEEK == 24,]$CD4} else {aids.base$CD424[i] = NA}
	if(dim(aids.rna[aids.rna$PATID == aids.base$PATID[i] & aids.rna$WEEK == 0,])[1] >0) {aids.base$RNABASE[i] = aids.rna[aids.rna$PATID == aids.base$PATID[i] & aids.rna$WEEK == 0,]$ULOGRNA} else {aids.base$RNABASE[i] = NA}
	if(dim(aids.rna[aids.rna$PATID == aids.base$PATID[i] & aids.rna$WEEK == 24,])[1] >0) {aids.base$RNA24[i] = aids.rna[aids.rna$PATID == aids.base$PATID[i] & aids.rna$WEEK == 24,]$ULOGRNA} else {aids.base$RNA24[i] = NA}
}

aids.base$CD4CHANGE = aids.base$CD424 - aids.base$CD4BASE
aids.base$RNACHANGE = aids.base$RNA24 - aids.base$RNABASE

library(hetsurr)
library(Rsurrogate)

aids.base = aids.base[!(is.na(aids.base$RNACHANGE) | is.na(aids.base$CD4CHANGE)),]

# Multiply outcome by -1 so that higher values are better, matching our assumptions
study.A <- data.frame(A = aids.base$TREAT - 1,
                      Y = -1*aids.base$RNACHANGE,
                      S = aids.base$CD4CHANGE,
                      W = aids.base$AVECD4)

sum(study.A$A == 1) # n1 = 418
sum(study.A$A == 0) # n0 = 412

# test for evidence of heterogeneity
set.seed(1)
het.ob = hetsurr.fun(y1 = study.A$Y[study.A$A == 1],
                     y0 = study.A$Y[study.A$A == 0],
                     s1 = study.A$S[study.A$A == 1],
                     s0 = study.A$S[study.A$A == 0],
                     w1 = study.A$W[study.A$A == 1],
                     w0 = study.A$W[study.A$A == 0],
                     c.adj=2, var.want =TRUE,  test.want = TRUE, type = "cont")
het.ob
# Evidence of significant heterogeneity

# Plot estimated PTE for Study A
library(ggplot2)
het_df <- data.frame(het.ob)
ggplot(het_df, aes(x = w.values)) +
  geom_line(aes(y = R.w.s), linewidth = 1.2, color = "black") +
  geom_ribbon(aes(ymin = band.R.w.s.lower, ymax = band.R.w.s.upper), 
              fill = "grey40", alpha = 0.4) +
  geom_hline(yintercept = c(0.5, 0.6), linetype = "dashed", size = 1) +
  labs(x = "Baseline CD4", y = expression(R[S])) +
  xlim(9, 185) +
  theme_minimal(base_size = 15) 

############################################
#### USING ACTG 193A AS STUDY B
################################################

# Set working directory to data location
setwd("C:/Users/rkkno/Documents/University of Texas at Austin/etsi/AIDS Data/ACTG 193A Data")
# Control group: Zidovudine and didanosine (2 NRTIs)
# Treatment group: Zidovudine and didanosine and nevirapine (2 NRTIs plus NNRTI)
# So let GROUP 0 = ZDV+ddI; GROUP 1 = ZDV+ddI+NVP

aids.193 = read.csv("193A_Base.csv", header = T)
aids.193 = aids.193[aids.193$TRT == "ZDV+ddI" | aids.193$TRT == "ZDV+ddI+NVP",]
aids.193$TRT = 1*(aids.193$TRT == "ZDV+ddI+NVP")
aids.193.cd4 = read.csv("193A_CD4.csv", header = T)
aids.193.rna = read.csv("193A_RNA.csv", header = T)

n.total = dim(aids.193)[1]
aids.193$CD4BASE = vector(length = n.total)
aids.193$CD424 = vector(length = n.total)
aids.193$RNABASE = vector(length = n.total)
aids.193$RNA24 = vector(length = n.total)

# Get baseline and 24 week CD4 and logRNA
for(i in 1:n.total) {
	if(dim(aids.193.cd4[aids.193.cd4$PATID == aids.193$PATID[i] & aids.193.cd4$VISIT == "Baseline",])[1] > 0) {aids.193$CD4BASE[i] = aids.193.cd4[aids.193.cd4$PATID == aids.193$PATID[i] & aids.193.cd4$VISIT == "Baseline",]$CD4} else {aids.193$CD4BASE[i] = NA}
	if(dim(aids.193.cd4[aids.193.cd4$PATID == aids.193$PATID[i] & aids.193.cd4$VISIT == "Week 24",])[1] > 0) {aids.193$CD424[i] = aids.193.cd4[aids.193.cd4$PATID == aids.193$PATID[i] & aids.193.cd4$VISIT == "Week 24",]$CD4} else {aids.193$CD424[i] = NA}
	if(dim(aids.193.rna[aids.193.rna$PATID == aids.193$PATID[i] & aids.193.rna$VISIT == "Baseline",])[1] > 0) {aids.193$RNABASE[i] = aids.193.rna[aids.193.rna$PATID == aids.193$PATID[i] & aids.193.rna$VISIT == "Baseline",]$LOGRNA} else {aids.193$RNABASE[i] = NA}
	if(dim(aids.193.rna[aids.193.rna$PATID == aids.193$PATID[i] & aids.193.rna$VISIT == "Week 24",])[1] > 0) {aids.193$RNA24[i] = aids.193.rna[aids.193.rna$PATID == aids.193$PATID[i] & aids.193.rna$VISIT == "Week 24",]$LOGRNA} else {aids.193$RNA24[i] = NA}
}

# Calculate change
aids.193$CD4CHANGE = aids.193$CD424 - aids.193$CD4BASE
aids.193$RNACHANGE = aids.193$RNA24 - aids.193$RNABASE

# Only use the observations that include RNA measurements
# Note this is a significant cut in sample size
studyb = aids.193[!is.na(aids.193$CD4CHANGE) & !is.na(aids.193$RNACHANGE),]

# Multiply outcome by -1 so that higher values are better, matching our assumptions
study.B <- data.frame(A = studyb$TRT,
                      Y = -1*studyb$RNACHANGE,
                      S = studyb$CD4CHANGE,
                      W = studyb$CD4BASE)

sum(study.B$A == 1) # n1 = 28
sum(study.B$A == 0) # n0 = 37

# Helper function to calculate delta col for strong surrogacy
get.delta.aids <- function(df, k) {
  closest.index <- sapply(df$W, function(w) {which.min(abs(het.ob$w.values - w))})
  return(het.ob$R.w.s[closest.index] > k)
}

# Export study A dataframe for study design purposes with delta based on k = 0.5
study.A$delta <- 1*get.delta.aids(study.A, 0.5)
setwd("C:/Users/rkkno/Documents/University of Texas at Austin/etsi")
write.table(study.A, paste("aids.studyA.txt",sep=""), quote = FALSE, row.names = FALSE)

# First run est.delta.B and est.delta.AB functions from main sim file
results.B <- est.delta.B(y1 = study.B$Y[study.B$A == 1], 
                         y0 = study.B$Y[study.B$A == 0])

control.A <- study.A[study.A$A == 0, ]

results.AB <- est.delta.AB(s1 = study.B$S[study.B$A == 1],
                           s0 = study.B$S[study.B$A == 0])

thresholds <- c(0.5, 0.6)
results.P <- matrix(nrow = 2, ncol = 3)
colnames(results.P) <- c("delta.P", "se.delta.P", "p.value")
rownames(results.P) <- c("k = 0.5", "k = 0.6")

for (k in 1:length(thresholds)) {
  study.A$delta <- 1*get.delta.aids(study.A, thresholds[k])
  study.B$delta <- 1*get.delta.aids(study.B, thresholds[k])
  results.P[k,] <- unlist(etsi.main(study.A, study.B))
}

results <- matrix(ncol = 4, nrow = 4)
rownames(results) <- c("Estimate", "SE", "Effect Size", "p")
colnames(results) <- c("Delta.B", "Delta.AB", "Delta.P: k = 0.5", "Delta.P: k = 0.6")

# fill in B results
results[1,1] <- results.B$delta.B
results[2,1] <- results.B$se.delta.B
# AB results
results[1,2] <- results.AB$delta.AB
results[2,2] <- results.AB$se.delta.AB
# P results
results[1:2,3] <- results.P[1,1:2]
results[1:2,4] <- results.P[2,1:2]

results[3,] <- results[1,] / results[2,]
results[4,] <- 2* (1 - pnorm(abs(results[3,])))

setwd("C:/Users/rkkno/Documents/University of Texas at Austin/etsi")
results.latex <- format(round(results ,3),nsmall=3)
latex.table(results.latex, "aidsres", caption = "", dcolumn = T)
print(round(results,3))



#################### CHECK ASSUMPTIONS #################### 

# Check assumptions in region of strong surrogacy, where estimated PTE > 0.5
threshold <- 0.5

# Strong surrogacy in Study A:
y1.strong <- y1[A.w1.closest.R >= threshold]
y0.strong <- y0[A.w0.closest.R >= threshold]
s1.strong <- s1[A.w1.closest.R >= threshold]
s0.strong <- s0[A.w0.closest.R >= threshold]
w1.strong <- w1[A.w1.closest.R >= threshold]
w0.strong <- w0[A.w0.closest.R >= threshold]

# Strong surrogacy in Study B:
studyb.y1.strong <- studyb.y1[R.studyb.w1 >= threshold]
studyb.y0.strong <- studyb.y0[R.studyb.w0 >= threshold]
studyb.s1.strong <- studyb.s1[R.studyb.w1 >= threshold]
studyb.s0.strong <- studyb.s0[R.studyb.w0 >= threshold]
studyb.w1.strong <- studyb.w1[R.studyb.w1 >= threshold]
studyb.w0.strong <- studyb.w0[R.studyb.w0 >= threshold]

### C1: conditional mean function is monotone DECREASING in S for A1, A0, B1, B0 ###

check.c1 <- function(y, s) {
  fit_smooth <- loess(y ~ s)
  preds <- predict(fit_smooth, newdata = data.frame(s = s), se = TRUE)
  plot.data <- data.frame(s = s,
                          y = preds$fit,
                          lower.ci = preds$fit - 1.96*preds$se.fit,
                          upper.ci = preds$fit + 1.96*preds$se.fit)
  ggplot(plot.data, aes(x = s, y = y)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.2, fill = "blue") +
    title("Estimated E(y | s)") +
    ylab("E(y | s)") +
    xlab("s")
}

check.c1(y1.strong, s1.strong) # OK
check.c1(y0.strong, s0.strong)# OK
check.c1(studyb.y1.strong, studyb.s1.strong) # lightly violated for negative s, but remember small sample size
check.c1(studyb.y0.strong, studyb.s0.strong) # OK



### C2: P(S1 > s) > P(S0 > s) for all s, in both studies A and B ###

check.c2 <- function(s1, s0) {
  # Combine both vectors to get a range of s values to test
  all_s <- sort(unique(c(s1, s0)))
  
  # Define empirical CDF functions for s1 and s0
  ecdf_s1 <- ecdf(s1)
  ecdf_s0 <- ecdf(s0)
  
  # Compute P(s1 > s) and P(s0 > s) for all s
  p_s1_greater <- 1 - ecdf_s1(all_s)
  p_s0_greater <- 1 - ecdf_s0(all_s)
  
  print(paste0("Proportion of s where P(S1 > s) > P(S0 > s): ", sum(p_s1_greater > p_s0_greater) / length(all_s), sep = ""))
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    s = all_s,
    p_s1_greater = p_s1_greater,
    p_s0_greater = p_s0_greater
  )
  ggplot(plot_data) +
    geom_line(aes(x = s, y = p_s1_greater, color = "s1"), size = 1.2) +
    geom_line(aes(x = s, y = p_s0_greater, color = "s0"), size = 1.2) +
    ggtitle("Comparison of P(s1 > s) and P(s0 > s)") +
    xlab("s") +
    ylab("P(S > s)") +
    scale_color_manual(values = c("s1" = "blue", "s0" = "red")) +
    theme_minimal() +
    theme(legend.title = element_blank())
}

check.c2(s1.strong, s0.strong) # OK
check.c2(studyb.s1.strong, studyb.s0.strong) # Mostly OK

### C3: E(Y1 |s) > E(Y0 | s) for all s, in both studies A and B ###

check.c3 <- function(s1, y1, s0, y0) {
  # Fit LOESS models for E(Y1 | s1) and E(Y0 | s0)
  fit_y1 <- loess(y1 ~ s1)
  fit_y0 <- loess(y0 ~ s0)
  
  all_s <- sort(unique(c(s1, s0)))
  
  # Predict E(Y1 | s1) and E(Y0 | s0)
  preds1 <- predict(fit_y1, newdata = data.frame(s1 = all_s), se = TRUE)
  preds0 <- predict(fit_y0, newdata = data.frame(s0 = all_s), se = TRUE)
  
  idx <- (!is.na(preds1$fit) & !is.na(preds0$fit))
  all_s <- all_s[idx]
  predicted_y1 <- preds1$fit[idx]
  predicted_y1_lower <- preds1$fit[idx] - 1.96 * preds1$se.fit[idx]
  predicted_y1_upper <- preds1$fit[idx] + 1.96 * preds1$se.fit[idx]
  predicted_y0 <- preds0$fit[idx]
  predicted_y0_lower <- preds0$fit[idx] - 1.96 * preds0$se.fit[idx]
  predicted_y0_upper <- preds0$fit[idx] + 1.96 * preds0$se.fit[idx]
  
  print(paste0("Proportion of s where E(y1 > s) < E(y0 > s): ", sum(predicted_y1 < predicted_y0) / length(all_s), sep = ""))
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    s = all_s,
    preds = c(predicted_y1, predicted_y0),
    lower = c(predicted_y1_lower, predicted_y0_lower),
    upper = c(predicted_y1_upper, predicted_y0_upper),
    group = rep(c("Y1", "Y0"), c(length(predicted_y1), length(predicted_y0)))
  )
  
  # Create the plot
  ggplot(plot_data, aes(x = s, y = preds, color = group, fill = group)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    ggtitle("Comparison of E(Y1 | s1) and E(Y0 | s0)") +
    xlab("s") +
    ylab("E(Y | s)") +
    scale_color_manual(values = c("Y1" = "blue", "Y0" = "red")) +
    theme_minimal() +
    theme(legend.title = element_blank())
}

check.c3(s1 = s1.strong, y1 = y1.strong, s0 = s0.strong, y0 = y0.strong) # OK
check.c3(s1 = studyb.s1.strong, y1 = studyb.y1.strong, s0 = studyb.s0.strong, y0 = studyb.y0.strong) # lightly violated


### C4: mu_{A0}(s) = mu_{B0}(s) ###

fit_y0.A <- loess(y0.strong ~ s0.strong)
fit_y0.B <- loess(studyb.y0.strong ~ studyb.s0.strong)

s <- sort(unique(c(s0.strong, studyb.s0.strong)))
preds.A <- predict(fit_y0.A, newdata = data.frame(s0.strong = s), se= TRUE)
preds.B <- predict(fit_y0.B, newdata = data.frame(studyb.s0.strong = s), se = TRUE)
idx <- (!is.na(preds.A$fit) & !is.na(preds.B$fit))

s <- s[idx]
pred.y.A <- preds.A$fit[idx]
pred.y.A.lower <- preds.A$fit[idx] - 1.96 * preds.A$se.fit[idx]
pred.y.A.upper <- preds.A$fit[idx] + 1.96 * preds.A$se.fit[idx]
pred.y.B <- preds.B$fit[idx]
pred.y.B.lower <- preds.B$fit[idx] - 1.96 * preds.B$se.fit[idx]
pred.y.B.upper <- preds.B$fit[idx] + 1.96 * preds.B$se.fit[idx]

plot_data <- data.frame(
  s = s,
  preds = c(pred.y.A, pred.y.B),
  lower = c(pred.y.A.lower, pred.y.B.lower),
  upper = c(pred.y.A.upper, pred.y.B.upper),
  group = rep(c("A", "B"), c(length(pred.y.A), length(pred.y.B)))
)

# Create the plot
ggplot(plot_data, aes(x = s, y = preds, color = group, fill = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  ggtitle("Comparison of E(Y0.A | s) and E(Y0.B | s)") +
  xlab("s") +
  ylab("E(Y | s)") +
  scale_color_manual(values = c("A" = "blue", "B" = "red")) +
  scale_fill_manual(values = c("A" = "blue", "B" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank())

# this assumption doesn't really seem to be met

### C5: support of S0_B, S1_B, and S1_A contained within support of S0_A ###

c5.data <- data.frame(
  value = c(s0.strong, s1.strong, studyb.s0.strong, studyb.s1.strong),
  group = rep(c("s0.strong", "s1.strong", "studyb.s0.strong", "studyb.s1.strong"),
              times = c(length(s0.strong), length(s1.strong), length(studyb.s0.strong), length(studyb.s1.strong)))
)

# Create the plot
ggplot(c5.data, aes(x = value, color = group, fill = group)) +
  geom_density(alpha = 0.3, adjust = 1) +  # Use alpha to make the areas transparent, adjust smooths the density estimate
  ggtitle("Empirical PDFs of s0.strong, s1.strong, studyb.s0.strong, and studyb.s1.strong") +
  xlab("Value") +
  ylab("Density") +
  scale_color_manual(values = c("s0.strong" = "blue", "s1.strong" = "red", "studyb.s0.strong" = "green", "studyb.s1.strong" = "purple")) +
  scale_fill_manual(values = c("s0.strong" = "blue", "s1.strong" = "red", "studyb.s0.strong" = "green", "studyb.s1.strong" = "purple")) +
  theme_minimal() +
  theme(legend.title = element_blank())

# everything should be contained within blue density
# seems mostly true, but on the right tail red isn't really contained, so slightly violated for large values of s1.strong