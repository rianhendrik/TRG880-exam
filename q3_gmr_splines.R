#Question 3
library(pracma)
d = read.csv("mix.csv")
n = nrow(d)
data_plot = plot(d,
                 main = "A scatterplot of mix.csv",
                 ylab = "y",
                 xlab = "x",
                 pch = 19)

#there are 2 components (clusters) from the looks of it.
K = 2

write(data_plot, "mix_plot.png")
("mix_plot.png")

#Using kmeans to obtain initial values
km = kmeans(d, centers = 2)
kmeans_plot = plot(d, col = km$cluster,
                   main = "A scatterplot of mix.csv \n
             Kmeans clusters denoted in black and pink.",
                   xlab = "x",
                   ylab = "y",
                   pch = 19)

clusters = cbind(km$cluster==1, km$cluster==2)

#Testing for equality of variance - H0 variances are not sig. different.
F = km$withinss[1]/km$withinss[2]
crit = qf(0.9, km$size[1], km$size[2])
paste0("The F test statistic ",F, "is less than the relevant 90th percentile of
       the F distribution, which is ", crit,".")
#Do not reject H0 at 10% level of significance.
#Variances can be assumed equal at this stage (with current cluster specification)

#initial parameters
pis.o = c(km$size)
#pooled variance (assuming equality of variance)
sigma.o = sqrt(sum(km$withinss)/n)


# x1_opt = cbind(int = matrix(1, km$size[1] , 1), d[,1],
#                ifelse(d[clusters[,1],1]<optimal_knots.n[1,1], 0, d[clusters[,1],1]-optimal_knots.n[1,2]),
#                ifelse(d[clusters[,1],1]<optimal_knots.n[1,2], 0, d[clusters[,1],1]-optimal_knots.n[1,2]))
# x2_opt = cbind(int = matrix(1, km$size[2], 1), d[,1],
#                ifelse(d[clusters[,2],1]<optimal_knots.n[2,1], 0, d[,1]-optimal_knots.n[2,1]),
#                ifelse(d[clusters[,2],1]<optimal_knots.n[2,2], 0, d[clusters[,2],1]-optimal_knots.n[2,2]))
# sigma = sqrt(sum(gammas[,1]*(d[,2] - x1_opt%*%betas.n[,1])^2 + gammas[,2]*(d[,2] - x2_opt%*%betas.n[,2])^2)/n)


#We will have to do a grid search to determine the optimal knot values
knots = seq(min(d[,1]), max(d[,1]), 20)
counter = 0
p = 4 #intercept, slope, and two hockey stick functions.
betas.o = matrix(NA, p, K)
optimal_knots.o = matrix(NA, K, 2)
startt = Sys.time()
for (k in 1:K){
  min_sse = sum(d[,2]^2)
  for (knot1 in knots){
    for (knot2 in knots[which(knots == knot1):length(knots)]){
      counter = counter + 1
      x = cbind(int = matrix(1, km$size[k], 1), d[clusters[,k],1],
                ifelse(d[clusters[,k],1]<knot1, 0, d[clusters[,k],1]-knot1),
                ifelse(d[clusters[,k],1]<knot2, 0, d[clusters[,k],1]-knot2))
      if(det(t(x)%*%x)>0){
        yhat = x%*%pinv(t(x)%*%x)%*%t(x)%*%d[clusters[,k],2]
        sse = sum(d[clusters[,k],2]-yhat)^2
        if (sse<min_sse){
          betas.o[,k] = pinv(t(x)%*%x)%*%t(x)%*%d[clusters[,k],2]
          min_sse = sse
          optimal_knots.o[k,] = c(knot1, knot2)
          print(counter)
          print(betas.o)
          print(optimal_knots.o)
          plot(d, col = km$cluster, pch = 19)
          print(sse)
          if (k == 1){lines(d[clusters[,1],1], yhat, lwd = 2, col = "blue", main = paste0(sse))}
          if (k == 2){lines(d[clusters[,2],1], yhat, lwd =2, col = "red", main = paste0(sse))}
        }
      }
    }
  }
}
duration = Sys.time() - startt
duration

#Plotting the initial regression
x1 = cbind(int = matrix(1, km$size[1], 1), d[clusters[,1],1],
           ifelse(d[clusters[,1],1]<optimal_knots.o[1,1], 0, d[clusters[,1],1]-optimal_knots.o[1,1]),
           ifelse(d[clusters[,1],1]<optimal_knots.o[1,2], 0, d[clusters[,1],1]-optimal_knots.o[1,2]))
x2 = cbind(int = matrix(1, km$size[2], 1), d[clusters[,2],1],
           ifelse(d[clusters[,2],1]<optimal_knots.o[2,1], 0, d[clusters[,2],1]-optimal_knots.o[2,1]),
           ifelse(d[clusters[,2],1]<optimal_knots.o[2,2], 0, d[clusters[,2],1]-optimal_knots.o[2,2]))

yhat1 = x1%*%pinv(t(x1)%*%x1)%*%t(x1)%*%d[clusters[,1],2]
yhat2 = x2%*%pinv(t(x2)%*%x2)%*%t(x2)%*%d[clusters[,2],2]

png("initialplot.png")
plot(d, pch = 19, col = km$cluster)
lines(d[clusters[,1],1], yhat1, lwd = 2, col = "blue")
lines(d[clusters[,2],1], yhat2, lwd =2, col = "red")
dev.off()

plot(d, pch = 19, col = km$cluster)
lines(d[clusters[,1],1], yhat1, lwd = 2, col = "blue")
lines(d[clusters[,2],1], yhat2, lwd =2, col = "red")

#The EM Algorithm

diffs = c()
tol = 0.001
iteration = 0
hard_clus = matrix(NA, n, 1)
diff = 1
while(diff > tol){
  iteration = iteration + 1
  
  #calculating belongings (expectation)
  nums = matrix(NA, n, 2)
  for (k in 1:K){
    x = cbind(int = matrix(1, n, 1), d[,1],
              ifelse(d[,1]<optimal_knots.o[k,1], 0, d[,1]-optimal_knots.o[k,1]),
              ifelse(d[,1]<optimal_knots.o[k,2], 0, d[,1]-optimal_knots.o[k,2]))
    nums[,k] = pis.o[k]*dnorm(d[,2], x%*%betas.o[,k], sigma.o)
  }
  gammas = nums/rowSums(nums)
  
  #determining the hard cluster solution & optimising pi
  
  for (i in 1:n){hard_clus[i] = which(gammas[i,] == max(gammas[i,]))}
  
  plot(d, col = hard_clus, pch = 19)
  
  hc = matrix(NA, n, 1)
  c1 = which(hard_clus == 1) #do not delete these c's!
  c2 = which(hard_clus == 2)
  hc[c1,] = 1
  hc[c2,] = 2
  
  
  #updating parameters (maximisation)
  
  pis.n = colSums(gammas)/n
  
  counter = 0
  p = 4 #intercept, slope, and two hockey stick functions.
  betas.n = matrix(NA, p, K)
  optimal_knots.n = matrix(NA, K, 2)
  for (k in 1:K){
    min_sse = 10000000000
    for (knot1 in knots){
      for (knot2 in knots[which(knots == knot1):length(knots)]){
        counter = counter + 1
        progress = round(counter/(length(knots)*(length(knots)+1))*100,4)
        print(paste0(progress,"%"))
        x = cbind(int = matrix(1, n, 1), d[,1],
                  ifelse(d[,1]<knot1, 0, d[,1]-knot1),
                  ifelse(d[,1]<knot2, 0, d[,1]-knot2))
        if(det(t(x)%*%x)>0){
          yhat = x%*%pinv(t(x)%*%diag(gammas[,k])%*%x)%*%t(x)%*%diag(gammas[,k])%*%d[,2]
          sse = sum(d[,2]-yhat)^2
          if (sse<min_sse){
            betas.n[,k] = pinv(t(x)%*%diag(gammas[,k])%*%x)%*%t(x)%*%diag(gammas[,k])%*%d[,2]
            min_sse = sse
            optimal_knots.n[k,] = c(knot1, knot2)
            plot(d, col = hc, pch = 19)
            if (k == 1){lines(d[,1], yhat, lwd = 2, col = "red")}
            if (k == 2){lines(d[,1], yhat, lwd = 2, col = "blue")}
          }
        }
      }
    }
  }

  
  
  #updating sigma
  x1_opt = cbind(int = matrix(1, n, 1), d[,1],
             ifelse(d[,1]<optimal_knots.n[1,1], 0, d[,1]-optimal_knots.n[1,1]),
             ifelse(d[,1]<optimal_knots.n[1,2], 0, d[,1]-optimal_knots.n[1,2]))
  x2_opt = cbind(int = matrix(1, n, 1), d[,1],
             ifelse(d[,1]<optimal_knots.n[2,1], 0, d[,1]-optimal_knots.n[2,1]),
             ifelse(d[,1]<optimal_knots.n[2,2], 0, d[,1]-optimal_knots.n[2,2]))
  sigma2 = sum(gammas[,1]*(d[,2] - x1_opt%*%betas.n[,1])^2 + gammas[,2]*(d[,2] - x2_opt%*%betas.n[,2])^2)/n
  sigma.n = sqrt(sigma2)
  
  #Plotting current results at ith iteration
  
  x1 = cbind(int = matrix(1, length(c1), 1), d[c1,1],
             ifelse(d[c1,1]<optimal_knots.n[1,1], 0, d[c1,1]-optimal_knots.n[1,1]),
             ifelse(d[c1,1]<optimal_knots.n[1,2], 0, d[c1,1]-optimal_knots.n[1,2]))
  x2 = cbind(int = matrix(1, length(c2), 1), d[c2,1],
             ifelse(d[c2,1]<optimal_knots.n[2,1], 0, d[c2,1]-optimal_knots.n[2,1]),
             ifelse(d[c2,1]<optimal_knots.n[2,2], 0, d[c2,1]-optimal_knots.n[2,2]))
  
  yhat1 = x1%*%pinv(t(x1)%*%x1)%*%t(x1)%*%d[c1,2]
  yhat2 = x2%*%pinv(t(x2)%*%x2)%*%t(x2)%*%d[c2,2]
  
  png(paste0("results", iteration,".png"))
  plot(d, pch = 19, col = hc)
  lines(d[c1,1], yhat1, lwd = 2, col = "blue")
  lines(d[c2,1], yhat2, lwd = 2, col = "red")
  dev.off()
  
  plot(d, pch = 19, col = hc)
  lines(d[c1,1], yhat1, lwd = 2, col = "blue")
  lines(d[c2,1], yhat2, lwd = 2, col = "red")
  
  diff = sum(abs(abs(betas.n)-abs(betas.o)))
  
  pis.o = pis.n
  betas.o = betas.n
  sigma.o = sigma.n
  
  #updating knots
  optimal_knots.o = optimal_knots.n
  
  diffs = rbind(diffs, diff)
  

}
