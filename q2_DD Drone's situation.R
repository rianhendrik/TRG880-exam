#DD Drone's situation
library(JWileymisc)
library(MASS)


n = c(20, 100, 1000, 10000, 100000, 1000000)
rs = c(0.75, 0.85, 0.915, 0.995)
true_betas = c(100, 0.12, 0.25, 0.76)

names = c("b1_20", "b2_20", "b3_20",
  "b1_100", "b2_100", "b3_100",
  "b1_1000", "b2_1000", "b3_1000",
  "b1_10k", "b2_10k", "b3_10k",
  "b1_100k", "b2_100k", "b3_100k",
  "b1_1m", "b2_1m", "b3_1m")


nsim = 1000

vifs_0.75 = matrix(NA, nsim, length(n)*3) #to ensure that we are (mostly) sticking within our desired VIF
vifs_0.85 = matrix(NA, nsim, length(n)*3)
vifs_0.915 = matrix(NA, nsim, length(n)*3)
vifs_0.995 = matrix(NA, nsim, length(n)*3)

betas_0.75 = matrix(NA, nsim, length(n)*3) #the estimates themselves. We can use these to evaluate the bias
betas_0.85 = matrix(NA, nsim, length(n)*3)
betas_0.915 = matrix(NA, nsim, length(n)*3)
betas_0.995 = matrix(NA, nsim, length(n)*3)

vars_0.75 = matrix(NA, nsim, length(n)*3)  #to see the variances of the estimates
vars_0.85 = matrix(NA, nsim, length(n)*3)
vars_0.915 = matrix(NA, nsim, length(n)*3)
vars_0.995 = matrix(NA, nsim, length(n)*3)

ts_0.75 = matrix(NA, nsim, length(n)*3)    #to verify if the estimates are significant
ts_0.85 = matrix(NA, nsim, length(n)*3)
ts_0.915 = matrix(NA, nsim, length(n)*3)
ts_0.995 = matrix(NA, nsim, length(n)*3)



colnames(betas_0.75) = names
colnames(betas_0.85) = names
colnames(betas_0.915) = names
colnames(betas_0.995) = names

colnames(vifs_0.75) = names
colnames(vifs_0.85) = names
colnames(vifs_0.915) = names
colnames(vifs_0.995) = names

colnames(vars_0.75) = names
colnames(vars_0.85) = names
colnames(vars_0.915) = names
colnames(vars_0.995) = names

colnames(ts_0.75) = names
colnames(ts_0.85) = names
colnames(ts_0.915) = names
colnames(ts_0.995) = names




#for VIF<5, x3 error var = 900, and x2 error = 900
#for 5<= VIF <10, x3 error = 450, and x2 error = 500
#for 10 <= VIF <10, x3 error =   , and x2 error = 



for( j in rs){
  
    V = matrix(c(1, j,
                j, 1),
             2, 2)
  
  #vector of standard deviations for variables
  sigma = c(10000, 1000)  #R10 000 sd for Revenue, and R1000 sd for Mine eq. spending
  mu = c(100000, 10000)   # R100 000 mean for Revenue, and R10 000 mean for Mine eq. spending
  
  Sigma = cor2cov(V, sigma)
  
  counter = 0
  for  (i in n){
    counter = counter + 1
    nums = c(0, 2, 4, 6, 8, 10)
    inum = nums[counter]
    
    for(iter in 1:nsim){
      
      d = data.frame(mvrnorm(i, mu, Sigma, 2, 2))
      x1 = d[,1]             #revenue
#      x2 = d[,2]             #mine eq. spend
      
      #create spending on aerial data service, 
      #as linear combination of other two, with correct adj. R2
      
      x2 = 0.1*x1 +rnorm(i, 0, 30)
      
      x3 = 0.074*x1 + 0.19*x2 + rnorm(i, 0, 30) 
          cor(cbind(x1, x2, x3))
      
      y = 100 + 0.12*x1 + 0.25*x2 + 0.76*x3 + rnorm(i, 0, 120) #these are the true betas
      
      d = data.frame(y, x1, x2, x3)
      #    cor(d)
      
      reg = lm(y ~ ., d)
      car::vif(reg)
      sum = summary(reg)
      
      
           
      
      if (j == rs[1]){vifs_0.75[iter, (counter+inum):(3*counter)] = car::vif(reg)
                  betas_0.75[iter, (counter+inum):(3*counter)] = reg$coefficients[2:4]
                  vars_0.75[iter, (counter+inum):(3*counter)] = sum[["coefficients"]][2:4,2]
                  ts_0.75[iter, (counter+inum):(3*counter)] = sum[["coefficients"]][2:4,3]}
      if (j == rs[2]){vifs_0.85[iter, (counter+inum):(3*counter)] = car::vif(reg)
                  betas_0.85[iter, (counter+inum):(3*counter)] = reg$coefficients[2:4]
                  vars_0.85[iter, (counter+inum):(3*counter)] = sum[["coefficients"]][2:4,2]
                  ts_0.85[iter, (counter+inum):(3*counter)] = sum[["coefficients"]][2:4,3]}
      if (j == rs[3]){vifs_0.915[iter, (counter+inum):(3*counter)] = car::vif(reg)
                  betas_0.915[iter, (counter+inum):(3*counter)] = reg$coefficients[2:4]
                  vars_0.915[iter, (counter+inum):(3*counter)] = sum[["coefficients"]][2:4,2]
                  ts_0.915[iter, (counter+inum):(3*counter)] = sum[["coefficients"]][2:4,3]}
      if (j == rs[4]){vifs_0.995[iter, (counter+inum):(3*counter)] = car::vif(reg)
                  betas_0.995[iter, (counter+inum):(3*counter)] = reg$coefficients[2:4]
                  vars_0.995[iter, (counter+inum):(3*counter)] = sum[["coefficients"]][2:4,2]
                  ts_0.995[iter, (counter+inum):(3*counter)] = sum[["coefficients"]][2:4,3]}
     
       #progress = paste0((iter/nsim)*100,"%"))
      #print(progress)

    }
    print(cat(j,i))
  }

  if (j == rs[1]){saveRDS(vifs_0.75, "vifs_0.75.rds")
                  saveRDS(betas_0.75, "betas_0.75.rds")
                  saveRDS(vars_0.75, "vars_0.75.rds")
                  saveRDS(ts_0.75, "ts_0.75")}
  
  if (j == rs[2]){saveRDS(vifs_0.85, "vifs_0.85.rds")
                  saveRDS(betas_0.85, "betas_0.85.rds")
                  saveRDS(vars_0.85, "vars_0.85.rds")
                  saveRDS(ts_0.85, "ts_0.85")}
  
  if (j == rs[3]){saveRDS(vifs_0.915, "vifs_0.915.rds")
                  saveRDS(betas_0.915, "betas_0.915.rds")
                  saveRDS(vars_0.915, "vars_0.915.rds")
                  saveRDS(ts_0.915, "ts_0.915")}
  
  if (j == rs[4]){saveRDS(vifs_0.995, "vifs_0.995.rds")
                  saveRDS(betas_0.995, "betas_0.995.rds")
                  saveRDS(vars_0.995, "vars_0.995.rds")
                  saveRDS(ts_0.995, "ts_0.995")}
}






