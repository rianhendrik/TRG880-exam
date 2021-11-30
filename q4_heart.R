library(readxl)
library(dplyr)
library(tidyverse)
library(splines)
library(ggplot2)
library(pracma)
library(splines)

setwd("/home/rian/Dropbox/2. TRG880/EXAM/Q4 - B-splines/")


d = read_xls("Heart.xls")


n = nrow(d)

d = d[,2:ncol(d)] %>%
  dplyr::select(-adiposity,-typea, -famhist, -alcohol)

# d01 = data.frame(d)
# for (i in 1:ncol(d)){
#   d01[,i] = d[,i]/max(d[,i])
# }

sbp = scale(d[,1])
tobacco = scale(d[,2])
ldl = scale(d[,3])
obesity = scale(d[,4])
age = scale(d[,5])
y = as.numeric(unlist(d[,6]))




## This function is a (simplified) R implementation of the bs()
## function in the splines library and illustrates how the Cox-de Boor
## recursion formula is used to construct B-splines.
basis <- function(x, degree, i, knots) {
  if(degree == 0){
    B <- ifelse((x >= knots[i]) & (x < knots[i+1]), 1, 0)
  } else {
    if((knots[degree+i] - knots[i]) == 0) {
      alpha1 <- 0
    } else {
      alpha1 <- (x - knots[i])/(knots[degree+i] - knots[i])
    }
    if((knots[i+degree+1] - knots[i+1]) == 0) {
      alpha2 <- 0
    } else {
      alpha2 <- (knots[i+degree+1] - x)/(knots[i+degree+1] - knots[i+1])
    }
    B <- alpha1*basis(x, (degree-1), i, knots) + alpha2*basis(x, (degree-1), (i+1), knots)
  }
  return(B)
}


bs <- function(x, degree=3, interior.knots=NULL, intercept=FALSE, Boundary.knots = c(0,1)) {
  if(missing(x)) stop("You must provide x")
  if(degree < 1) stop("The spline degree must be at least 1")
  Boundary.knots <- sort(Boundary.knots)
  interior.knots.sorted <- NULL
  if(!is.null(interior.knots)) interior.knots.sorted <- sort(interior.knots)
  knots <- c(rep(Boundary.knots[1], (degree+1)), interior.knots.sorted, rep(Boundary.knots[2], (degree+1)))
  K <- length(interior.knots) + degree + 1
  B.mat <- matrix(0,length(x),K)
  for(j in 1:K) B.mat[,j] <- basis(x, degree, j, knots)
  if(any(x == Boundary.knots[2])) B.mat[x == Boundary.knots[2], K] <- 1
  if(intercept == FALSE) {
    return(B.mat[,-1])
  } else {
    return(B.mat)
  }
}

#Building sbp basis
ps = c(0.25, 0.5, 0.75)
ps_string = c("0.25, 0.5, 0.75")
degree = 1

k_sbp = unname(quantile(sbp, probs = ps))
B_sbp = bs(sbp, degree=degree, interior.knots = k_sbp, 
           intercept = T, Boundary.knots=c(min(sbp), max(sbp)))


t = splines::bs(sbp, knots = k_sbp, degree = 1)
matplot(sbp, t, type="l", lwd=2,
        main = "1st Degree B-spline fitted to sbp",
        xlab = "sbp",
        ylab = "B(sbp)")
#Building tobacco basis
k_tob = unname(quantile(tobacco, probs = ps))
B_tob = bs(tobacco, degree=degree, interior.knots = k_tob,
           intercept = F, Boundary.knots=c(min(tobacco), max(tobacco)))

#Building ldl basis
k_ldl = unname(quantile(ldl, probs = ps))
B_ldl = bs(ldl, degree=degree, interior.knots = k_ldl,
           intercept = F, Boundary.knots=c(min(ldl), max(ldl)))

#Building obesity basis
k_ob = unname(quantile(obesity, probs = ps))
B_ob = bs(obesity, degree=degree, interior.knots = k_ob,
          intercept = F, Boundary.knots=c(min(obesity), max(obesity)))

#Building age basis
k_age = unname(quantile(age, probs = ps))
B_age = bs(age, degree=degree, interior.knots = k_age,
           intercept = F, Boundary.knots=c(min(age), max(age)))

#Now that we have built our 5 bases (each with m basis functions),
#we can use them to fit our logistic regression model.
x = data.matrix(data.frame(B_sbp, B_tob, B_ldl, B_ob, B_age))
x = scale(x)
x = cbind(intercept = matrix(1, n, 1), x) #standardised X

#Question 4.2 - 4.5
beta = matrix(0, nrow = ncol(x), ncol = 1)
prob = function(x, beta){
  return(1/(1+exp(-x%*%beta)))
}

tol = 0.0001
crit = 1
loglik = c()
counter = 0
while (crit > tol){
  counter = counter+1
  probs = prob(x, beta)
  W = diag(diag(probs%*%t((1-probs))), n, n)
  #score = t(x)%*%(y - probs)
  #hessian = -t(x)%*%W%*%x
  #beta_new = beta + solve(hessian)%*%score
  z = x%*%beta + solve(W)%*%(y-probs)
  #beta_new = solve(t(x)%*%W%*%x)%*%t(x)%*%W%*%z
  beta_new = pinv(t(x)%*%W%*%x)%*%t(x)%*%W%*%z
  crit = max(beta-beta_new)
  beta = beta_new
  progress = paste0("Iteration #",counter,"\n",
                    "Covergence = ",crit)
  loglik = cbind(loglik, sum(matrix(y)*x%*%beta - log(1+exp(x%*%beta))))
  print(progress)
}

BIC = -2*loglik[length(loglik)] + ncol(x)*log(nrow(x))
AIC = 2*ncol(x)-2*loglik[length(loglik)]
BIC
AIC
####Covariance matrix of betas (SIGMA) ####

#inv(H`*W*H)

#sigma = solve(t(x)%*%W%*%x)
sigma = pinv(t(x)%*%W%*%x)

#For indexing
basis_funcs = ncol(B_sbp)

#####Plots####

#sbp plot:


x_sbp = x[,2:basis_funcs+1]
x_sbp = x_sbp[order(x_sbp[,1]),]

#sbp variance plot:
betas_sbp = as.matrix(unname(beta[2:basis_funcs+1,]))
yhat_sbp = x[,2:basis_funcs+1]%*%betas_sbp

sigma_sbp = sigma[2:basis_funcs+1, 2:basis_funcs+1]

var_sbp = x[,2:basis_funcs+1]%*%sigma_sbp%*%t(x[,2:basis_funcs+1])
pointwise_var_sbp = diag(var_sbp)
ci_upper_sbp = yhat_sbp + 2*sqrt(pointwise_var_sbp)
ci_lower_sbp = yhat_sbp - 2*sqrt(pointwise_var_sbp)

sbp_plot = data.frame(d$sbp, yhat_sbp, ci_lower_sbp, ci_upper_sbp)
sbp_plot = sbp_plot[order(sbp_plot$d.sbp),]


plot(sbp_plot[,1], sbp_plot[,2], type = "l", ylim = c(-3, 3))
lines(sbp_plot[,1], sbp_plot[,4], col = "red")
lines(sbp_plot[,1], sbp_plot[,3], col = "blue")

x11()

the_sbp_plot = ggplot(sbp_plot, aes(sbp_plot[,1])) + ggtitle(paste0("Plot of sbp (degree", degree, " B-Spline) \n Interior knots at ", ps_string)) + theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="f_hat(sbp)", x = "sbp") +
  geom_line(aes(y=sbp_plot[,2]), colour="green") + 
  geom_ribbon(aes(ymin=sbp_plot[,3], ymax=sbp_plot[,4]), alpha=0.2)

the_sbp_plot

ggsave(
  paste0("sbp_d",degree,".png"),
  plot = last_plot())

#tobacco plot

betas_tob = as.matrix(unname(beta[(basis_funcs+2):(basis_funcs*2+1)]))
yhat_tob = x[,(basis_funcs+2):(basis_funcs*2+1)]%*%betas_tob

sigma_tob = sigma[(basis_funcs+2):(basis_funcs*2+1), 
                  (basis_funcs+2):(basis_funcs*2+1)]

var_tob = x[,(basis_funcs+2):(basis_funcs*2+1)]%*%sigma_tob%*%t(x[,(basis_funcs+2):(basis_funcs*2+1)])
pointwise_var_tob = diag(var_tob)
ci_upper_tob = yhat_tob + 2*sqrt(pointwise_var_tob)
ci_lower_tob = yhat_tob - 2*sqrt(pointwise_var_tob)

tob_plot = data.frame(d$tobacco, yhat_tob, ci_lower_tob, ci_upper_tob)
tob_plot = tob_plot[order(tob_plot$d.tobacco),]

plot(tob_plot[,1], tob_plot[,2], type = "l")
lines(tob_plot[,1], tob_plot[,4], col = "red")
lines(tob_plot[,1], tob_plot[,3], col = "blue")

x11()

the_tob_plot = ggplot(tob_plot, aes(tob_plot[,1])) + 
  ggtitle(paste0("Plot of tobacco (degree", degree, " B-Spline) \n Interior knots at ", ps_string)) + theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="f_hat(tob)", x = "tob") +
  geom_line(aes(y=tob_plot[,2]), colour="green") + 
  geom_ribbon(aes(ymin=tob_plot[,3], ymax=tob_plot[,4]), alpha=0.2)

the_tob_plot

ggsave(
  paste0("tob_d",degree,".png"),
  plot = last_plot())

#ldl plot

betas_ldl = as.matrix(unname(beta[(2*basis_funcs+2):(basis_funcs*3+1)]))
yhat_ldl = x[,(2*basis_funcs+2):(basis_funcs*3+1)]%*%betas_ldl

sigma_ldl = sigma[(2*basis_funcs+2):(basis_funcs*3+1), (2*basis_funcs+2):(basis_funcs*3+1)]

var_ldl = x[,(2*basis_funcs+2):(basis_funcs*3+1)]%*%sigma_ldl%*%t(x[,(2*basis_funcs+2):(basis_funcs*3+1)])
pointwise_var_ldl = diag(var_ldl)
ci_upper_ldl = yhat_ldl + 2*sqrt(pointwise_var_ldl)
ci_lower_ldl = yhat_ldl - 2*sqrt(pointwise_var_ldl)

ldl_plot = data.frame(d$ldl, yhat_ldl, ci_lower_ldl, ci_upper_ldl)
ldl_plot = ldl_plot[order(ldl_plot$d.ldl),]


plot(ldl_plot[,1], ldl_plot[,2], type = "l")
lines(ldl_plot[,1], ldl_plot[,4], col = "red")
lines(ldl_plot[,1], ldl_plot[,3], col = "blue")

x11()

the_ldl_plot = ggplot(ldl_plot, aes(ldl_plot[,1])) + 
  ggtitle(paste0("Plot of ldl (degree", degree, " B-Spline) \n Interior knots at ", ps_string)) + theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="f_hat(ldl)", x = "ldl") +
  geom_line(aes(y=ldl_plot[,2]), colour="green") + 
  geom_ribbon(aes(ymin=ldl_plot[,3], ymax=ldl_plot[,4]), alpha=0.2)

the_ldl_plot

ggsave(
  paste0("ldl_d",degree,".png"),
  plot = last_plot())

# plot obesity
betas_ob = as.matrix(unname(beta[(3*basis_funcs+2):(basis_funcs*4+1),]))
yhat_ob = x[,(3*basis_funcs+2):(basis_funcs*4+1)]%*%betas_ob

sigma_ob = sigma[(3*basis_funcs+2):(basis_funcs*4+1), (3*basis_funcs+2):(basis_funcs*4+1)]

var_ob = x[,(3*basis_funcs+2):(basis_funcs*4+1)]%*%sigma_ob%*%t(x[,(3*basis_funcs+2):(basis_funcs*4+1)])
pointwise_var_ob = diag(var_ob)
ci_upper_ob = yhat_ob + 2*sqrt(pointwise_var_ob)
ci_lower_ob = yhat_ob - 2*sqrt(pointwise_var_ob)

ob_plot = data.frame(d$obesity, yhat_ob, ci_lower_ob, ci_upper_ob)
ob_plot = ob_plot[order(ob_plot$d.obesity),]


plot(ob_plot[,1], ob_plot[,2], type = "l")
lines(ob_plot[,1], ob_plot[,4], col = "red")
lines(ob_plot[,1], ob_plot[,3], col = "blue")

x11()

the_ob_plot = ggplot(ob_plot, aes(ob_plot[,1])) + 
  ggtitle(paste0("Plot of obesity (degree", degree, " B-Spline) \n Interior knots at ", ps_string)) + theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="f_hat(Obesity)", x = "Obesity") +
  geom_line(aes(y=ob_plot[,2]), colour="green") + 
  geom_ribbon(aes(ymin=ob_plot[,3], ymax=ob_plot[,4]), alpha=0.3)

the_ob_plot

ggsave(
  paste0("ob_d",degree,".png"),
  plot = last_plot())

# plot age
betas_age = as.matrix(unname(beta[(4*basis_funcs+2):(basis_funcs*5+1),]))
yhat_age = x[,(4*basis_funcs+2):(basis_funcs*5+1)]%*%betas_age

sigma_age = sigma[(4*basis_funcs+2):(basis_funcs*5+1), 
                  (4*basis_funcs+2):(basis_funcs*5+1)]

var_age = x[,(4*basis_funcs+2):(basis_funcs*5+1)]%*%sigma_age%*%t(x[,(4*basis_funcs+2):(basis_funcs*5+1)])
pointwise_var_age = diag(var_age)
ci_upper_age = yhat_age + 2*sqrt(pointwise_var_age)
ci_lower_age = yhat_age - 2*sqrt(pointwise_var_age)

age_plot = data.frame(d$age, yhat_age, ci_lower_age, ci_upper_age)
age_plot = age_plot[order(age_plot$d.age),]


plot(age_plot[,1], age_plot[,2], type = "l")
lines(age_plot[,1], age_plot[,4], col = "red")
lines(age_plot[,1], age_plot[,3], col = "blue")

x11()

the_age_plot = ggplot(age_plot, aes(age_plot[,1])) + 
  ggtitle(paste0("Plot of age (degree", degree, " B-Spline) \n Interior knots at ", ps_string)) + theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="f_hat(Age)", x = "Age") +
  geom_line(aes(y=age_plot[,2]), colour="green") + 
  geom_ribbon(aes(ymin=age_plot[,3], ymax=age_plot[,4]), alpha=0.3)

the_age_plot

ggsave(
  paste0("age_d",degree,".png"),
  plot = last_plot())



