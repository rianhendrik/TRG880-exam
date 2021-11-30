#Getting table values from DD Drone's simulations
setwd("/home/rian/Dropbox/2. TRG880/EXAM/Q2 - MC and Coefficients")
options(scipen = 999)

true_betas = c(0.12, 0.25, 0.76)

t = matrix(true_betas, nrow = 1, ncol = 18 )

#VIF<5
betas_0.75 = readRDS("VIF less than 5/betas_0.75.rds")
se_0.75 = readRDS("VIF less than 5/vars_0.75.rds")
ts_0.75 = readRDS("VIF less than 5/ts_0.75")

average_betas = data.frame(colMeans(betas_0.75))
average_se = data.frame(colMeans(se_0.75))
average_t = data.frame(colMeans(ts_0.75))
bias_0.75 = data.frame(average_betas - true_betas)

#5<=VIF<10
betas_0.85 = readRDS("VIF between 5 and 10/betas_0.85.rds")
se_0.85 = readRDS("VIF between 5 and 10/vars_0.85.rds")
ts_0.85 = readRDS("VIF between 5 and 10/ts_0.85")

average_betas = data.frame(colMeans(betas_0.85))
average_se = data.frame(colMeans(se_0.85))
average_t = data.frame(colMeans(ts_0.85))
bias_0.85 = data.frame(average_betas - true_betas)

#10<=VIF<15
betas_0.915 = readRDS("VIF between 10 and 15/betas_0.915.rds")
se_0.915 = readRDS("VIF between 10 and 15/vars_0.915.rds")
ts_0.915 = readRDS("VIF between 10 and 15/ts_0.915")

average_betas = data.frame(colMeans(betas_0.915))
average_se = data.frame(colMeans(se_0.915))
average_t = data.frame(colMeans(ts_0.915))
bias_0.915 = data.frame(average_betas - true_betas)

#VIF>100
betas_0.995 = readRDS("VIF greater then 100 (much higher)/betas_0.995_higher.rds")
se_0.995 = readRDS("VIF greater then 100 (much higher)/vars_0.995_higher.rds")
ts_0.995 = readRDS("VIF greater then 100 (much higher)/ts_0_higher.995")


average_betas = data.frame(colMeans(betas_0.995))
average_se = data.frame(colMeans(se_0.995))
average_t = data.frame(colMeans(ts_0.995))
bias_0.995 = data.frame(average_betas - true_betas)

iter1 = seq(1, 18, 3)
iter2 = seq(2, 18, 3)
iter3 = seq(3, 18, 3)
b1s = abs(bias_0.995[iter1,])
b2s = abs(bias_0.995[iter2,])
b3s = abs(bias_0.995[iter3,])
n = c(20, 100, 1000, 10000, 100000, 1000000)

seb1 = average_se[iter1,]
seb2 = average_se[iter2,]
seb3 = average_se[iter3,]


par(mfrow = c(1,1))
op = par(cex = 2.5)
png("se_plot.png", width = 700)
plot(n, seb1, col = "red", type = "b", lwd = 2, pch = 16,
     main = "Regression parameter standard errors against sample size",
     xlab = "Sample size",
     ylab = "Standard errors of betas",
     ylim = c(0, max(seb1, seb2, seb3)))
lines(n, seb2, col = "green", type = "b", lwd = 5, pch = 15)
lines(n, seb3, col = "blue", type = "b", lwd = 2, pch = 19)
legend(x ="topright", legend=c("b1", "b2", "b3"),
       col=c("red", "green", "blue"), 
       lty=1, cex=1.5,
       text.font = 2)
dev.off()

par(mfrow = c(1,3))
plot(n, seb1, col = "red", type = "b", lwd = 2, pch = 19,
     main = "Beta 1 bias for different sample sizes",
     xlab = "Sample size",
     ylab = "Beta 1 bias")
plot(n, seb2, col = "green", type = "b", lwd = 2, pch = 19,
     main = "Beta 2 bias for different sample sizes",
     xlab = "Sample size",
     ylab = "Beta 2 bias")
plot(n, seb3, col = "blue", type = "b", lwd = 2, pch = 19,
     main = "Beta 3 bias for different sample sizes",
     xlab = "Sample size",
     ylab = "Beta 3 bias")


