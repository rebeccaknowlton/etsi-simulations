library(matrixStats)
library(quantreg)

#set the working directory to where the txt files are
setwd("C:/Users/rkkno/Documents/University of Texas at Austin/etsi")
setting <- 3


output <- read.table(paste("etsi.output", setting, "_110424",".txt", sep=""), header = T)


colMeans(output)

# proportion rejected 
apply(output[,c(3,6,9,12,15)], 2, function(x) mean(x < 0.05))

results <- matrix(ncol = 5, nrow = 6)
if (setting == 1 | setting == 2) {
  rownames(results) <- c("Point Estimate", "ESE", "ASE", "Effect Size", "Power", "pi_B")  
} else if (setting == 3) {
  rownames(results) <- c("Point Estimate", "ESE", "ASE", "Effect Size", "Type 1 Error", "pi_B")
}
colnames(results) <- c("Delta.B", "Delta.AB", "Delta.P: k > 0.5", "Delta.P: k > 0.6", "Delta.P: k > 0.7")


results[1,] <- colMeans(output[,c(1,4,7,10,13)])
results[2,] <- sqrt(colVars(as.matrix(output[,c(1,4,7,10,13)])))
results[3,] <- colMeans(output[,c(2,5,8,11,14)])
results[4,] <- colMeans(output[,c(1,4,7,10,13)] / output[,c(2,5,8,11,14)])
results[5,] <- apply(output[,c(3,6,9,12,15)], 2, function(x) mean(x < 0.05))
set.seed(1)
results[6,] <- c(0, 1, get.truth(setting)$pi)

results.latex <- format(round(results ,3),nsmall=3)
latex.table(results.latex, paste0("simres", setting), caption = "", dcolumn = T)
print(round(results,3))

