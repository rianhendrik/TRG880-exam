#credit is given to Mike Crowson and his YouTube video titled
#"Steps for simulating multivariate normal data in R" in assisting
#with the below code. His video may be viewed here:
#https://www.youtube.com/watch?v=Wk5wlkWY3rw&ab_channel=MikeCrowson

library(JWileymisc)
library(MASS)

n = 1000

#V=symmetric matrix of correlations among variables
V = matrix(c(1, 0.9567,
             0.9567, 1),
           2, 2)

#vector of standard deviations for variables
sigma = c(5.578206, 3.177563)
mu = c(21.28791, 14.87100)

Sigma = cor2cov(V, sigma)

d = data.frame(mvrnorm(n, mu, Sigma, 5, 5))

#Creating multicollinear variables
x1 = d[,1]
x2 = d[,2]
x3 = 0.5*x1 + 1.1*x2 + rnorm(n, 0, 0.002)
x4 = x1 + 0.9*x2 - x3 + rnorm(n, 0, 0.002)
t = cbind(x1, x2, x3, x4)
#confirming high correlations
cor(t)

#creating the dependent variables
y = 104 - 20*x1 + 13*x2 +         #creating the dependent variable
  3.1*x3 + 40*x4 +rnorm(n, 0, ) 
d = data.frame(y, x1, x2, x3, x4)   #creating a dataset
cor(d)                              #checking correlation of dataset

mod = lm(y ~ x1 + x2 + x3 + x4, d)  #fitting model
summary(mod)                        #note the insignificant variable estimates
car::vif(mod)                       #vif shows beyond a doubt that we are dealing wiht MC.

write.csv(d, "mc_data2.csv")




