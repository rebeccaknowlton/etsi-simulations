library(etsi)

# options for n, kappa, and psi
n.vec <- c(25, 50, 75, 100, 125, 150, 175)
kappa.vec <- c(0.5, 0.6)
psi.vec <- c(0.5, 0.75, 1.0)

options <- expand.grid(n.vec, kappa.vec, psi.vec)
colnames(options) <- c("n","kappa", "psi")

# read in Study A
study.A <- read.table(paste("aids_output/aids.studyA.txt", sep=""), header = T)
study.A <- study.A[,-5]

power.res <- rep(NA, nrow(options))
for (i in 1:nrow(options)) {
  n <- options[i,1]
  kappa <- options[i,2]
  psi <- options[i,3]
  set.seed(1) # set seed before each run of etsi.design so they have same train/test splits
  power.res[i] <- etsi.design(study.A, n.b0 = n, n.b1 = n, psi = psi, kappa = kappa)$power 
  print(i)
}

results.df <- data.frame("sample.size" = options[,1]*2, # multiply by 2 for n1 + n0 = n
                         "kappa" = options[,2],
                         "psi" = options[,3],
                         "power" = power.res)

# plot 
library(ggplot2)
# Plot sample size vs power with lines for different combinations of kappa and psi
fig3 <- ggplot(results.df, aes(x = power, y = sample.size, color = factor(psi), linetype = factor(kappa))) +
  geom_line(size = 0.9) +
  geom_point() +
  labs(
    x = "Power",
    y = "Sample Size",
    color = "Psi",
    linetype = "Kappa") +
  theme_minimal() +
  theme(legend.position = "right") 
ggsave("results/fig3.png", fig3)
