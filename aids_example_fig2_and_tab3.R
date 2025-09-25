library(etsi)
library(quantreg)
source("etsi_helper_functions.R")

############################################
#### USING ACTG 320 AS STUDY A
################################################


tmpg <- function(xx){exp(xx)}; tmpginv <- function(xx){log(xx)}
tmpg <- function(xx){xx^2}; tmpdg <- function(xx){2*xx}; tmpginv <- function(xx){sqrt(xx)}
log.link <- list(tmpg, tmpg, tmpginv); rm(tmpg, tmpginv)
tmpdir <- "./"
aids.base <- read.table(paste(tmpdir, "aids_data/ACTG320-Data/actg320.base.dat", sep=""))
aids.rna <- read.table(paste(tmpdir, "aids_data/ACTG320-Data/actg320.rna", sep=""))
aids.cd4 <- read.table(paste(tmpdir, "aids_data/ACTG320-Data/actg320.cd4", sep=""))
aids.trt <- read.table(paste(tmpdir, "aids_data/ACTG320-Data/actg320.trt", sep=""))

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
fig2 <- ggplot(het_df, aes(x = w.values)) +
  geom_line(aes(y = R.w.s), linewidth = 1.2, color = "black") +
  geom_ribbon(aes(ymin = band.R.w.s.lower, ymax = band.R.w.s.upper), 
              fill = "grey40", alpha = 0.4) +
  geom_hline(yintercept = c(0.5, 0.6), linetype = "dashed", size = 1) +
  labs(x = "Baseline CD4", y = expression(R[S])) +
  xlim(9, 185) +
  theme_minimal(base_size = 15) 
ggsave("results/fig2.png", fig2)

############################################
#### USING ACTG 193A AS STUDY B
################################################

# Control group: Zidovudine and didanosine (2 NRTIs)
# Treatment group: Zidovudine and didanosine and nevirapine (2 NRTIs plus NNRTI)
# So let GROUP 0 = ZDV+ddI; GROUP 1 = ZDV+ddI+NVP

aids.193 = read.csv("aids_data/ACTG 193A Data/193A_Base.csv", header = T)
aids.193 = aids.193[aids.193$TRT == "ZDV+ddI" | aids.193$TRT == "ZDV+ddI+NVP",]
aids.193$TRT = 1*(aids.193$TRT == "ZDV+ddI+NVP")
aids.193.cd4 = read.csv("aids_data/ACTG 193A Data/193A_CD4.csv", header = T)
aids.193.rna = read.csv("aids_data/ACTG 193A Data/193A_RNA.csv", header = T)

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
write.table(study.A, paste("aids_output/aids.studyA.txt",sep=""), quote = FALSE, row.names = FALSE)

results.B <- est.delta.B(y1 = study.B$Y[study.B$A == 1], 
                         y0 = study.B$Y[study.B$A == 0])

control.A <- study.A[study.A$A == 0, ]

results.AB <- est.delta.AB(s1 = study.B$S[study.B$A == 1],
                           s0 = study.B$S[study.B$A == 0],
                           control.A = control.A)

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

tab3 <- results
saveRDS(tab3, file = "results/tab3.rds")
print(round(tab3,3))


