#### Overview ####

# This file reproduces the simulation results in "Censored panel 
# quantile regression with fixed effects via an asymmetric link 
# function" by Fulden Komuryakan & Selahattin Guris.

# CORRESPONDING AUTHOR:   Fulden Komuryakan
# AFFILIATION:            Bandırma Onyedi Eylül University
# EMAIL:                  fkomuryakan@bandirma.edu.tr
# R version 4.2.2

rm(list=ls())

#### Libraries ####

library(quantreg)
library(MASS)
library(splines)
library(gam)
library(mfx)

#### Data Generating Process ####

beta1 <- 10
beta2 <- -2
beta <- c(beta1,beta2)
gamma <- 0.5
Ci <- -0.95
# Ci <- -10
# tau <- 0.10
# tau <- 0.25
tau <- 0.5
# tau <- 0.75
# tau <- 0.9
const <- 0.05
deltan <- const
Rep <- 500 
iter <- 0

datagen <- function(n,t,beta1,beta2,gamma,Ci) {
  s <- rep(1:n,rep(t,n))
  X <-array(0,c(n*t,2))
  for (i in 1:(n*t))
  {
    signal=0
    while (signal==0){
      x<-rnorm(2)
      if (max(abs(x))<=2) {X[i,]<-x;signal=1}
    }
  }
  X1 <- X[,1]
  X2 <- X[,2]
  X1s <- X1^2
  X2s <- X2^2
  X<-cbind(X1,X2,X1s,X2s)
  eta.aux <- rnorm(n)
  medX <- (X[,1]+X[,2])/2
  eta<-array(0,c(n,1))
  eta[1]<-eta.aux[1]+(1/sqrt(1))*sum(medX[1:t])
  for (k in 1:(n-1)){
    eta[k+1]<-eta.aux[k+1]+(1/sqrt(1))*sum(medX[(k*t+1):((k+1)*t)])
  }
  eta <- rep(eta,rep(t,n))
  u <- rnorm(n*t)
  ystar <- eta + beta1*X1+ beta2*X2 + (1+gamma*(X1+X2+X1s+X2s))*u
  y <-replace(ystar, ystar < Ci, Ci)
  delta <- 1-((y==Ci)*1)
  return(cbind(y,X,s,delta,ystar))
}

#### Simulation ####

# t <- 5
t <- 15
# t <- 50
n <- 300
set.seed(14)

for (j in 1:Rep)
{
  
  iter <- iter+1
  print(iter)
  data <- datagen(n,t,beta1,beta2,gamma,Ci)
  y <- data[,1]
  X1 <- data[,2]
  X2 <- data[,3]
  X <- data[,2:5]
  s <- data[,6]
  delta <- data[,7]
  ystar <- data[,8]
  yc <- rep(Ci,n*t)
  Z <- as.matrix(model.matrix(~as.factor(s)-1))
  tuning <- (n*t)^(-1/5) * tau
  
  #### Two-step estimator ####
  
  #### Generalized linear models #### 
  
  #### Generalized linear models - Logit #### 

  lm.fit.p2.step.1.logit <- glm(delta ~ X1 + X2 + Z,family=binomial(link = "logit"))
  results.p2step.1.logit <- summary(lm.fit.p2.step.1.logit,se="ker")$coef[1:2,]
  logit.fitted.p2.step.1 <- fitted(lm.fit.p2.step.1.logit)
  logit.aic.p2.step.1 <- lm.fit.p2.step.1.logit$aic
  c.logit.aic <- c(logit.aic.p2.step.1)
  logit.dev.p2.step.1 <- lm.fit.p2.step.1.logit$deviance
  c.logit.deviance <- c(logit.dev.p2.step.1)
  logit.null.dev.p2.step.1 <- lm.fit.p2.step.1.logit$null.deviance
  c.logit.null.deviance <- c(logit.null.dev.p2.step.1)
  stat.p1step.logit <- cbind(c.logit.aic, c.logit.deviance, c.logit.null.deviance)
  
  JTplogit<-logit.fitted.p2.step.1>(1-tau+tuning)
  
  dataJTplogit<-as.data.frame(data[JTplogit,])
  y.s2s<-dataJTplogit[,1]
  x1.s2s<-dataJTplogit[,2]
  x2.s2s<-dataJTplogit[,3]
  x1s.s2s<-dataJTplogit[,4]
  x2s.s2s<-dataJTplogit[,5]
  s.s2s<-dataJTplogit[,6]
  X.s2s<-cbind(x1.s2s,x2.s2s)
  Z.s2s <- as.matrix(model.matrix(~as.factor(s.s2s)-1))
  n.s2s<-ncol(Z.s2s)
  
  as.matrix.csr(cbind(X.s2s,Z.s2s)) -> sM.s2s
  coefs.step2s <- rq.fit.sfn(sM.s2s,y.s2s,tau=tau)$coef
  
  fit.p2.step.2.logit <- rq(y.s2s ~ X.s2s + Z.s2s - 1,tau=tau)
  results.p2step.logit <- summary(fit.p2.step.2.logit,se="ker")$coef[1:2,]
  residuals.p2step.logit <- summary(fit.p2.step.2.logit,se="ker")$residuals
  res.logit <- c(residuals.p2step.logit)
  RMSE.p.2.step.logit <- sqrt(mean((res.logit)^2))
  bias.logit <- c(beta-results.p2step.logit[1:2,1])
  stat.p2step.logit <- cbind(RMSE.p.2.step.logit, bias.logit)
  
  #### Generalized linear models - Cloglog ####
  
  lm.fit.p2.step.1.cloglog <- glm(delta ~ X1 + X2 + Z,family=binomial(link = "cloglog"))
  results.p2step.1.cloglog <- summary(lm.fit.p2.step.1.cloglog,se="ker")$coef[1:2,]
  cloglog.fitted.p2.step.1 <- fitted(lm.fit.p2.step.1.cloglog)
  cloglog.aic.p2.step.1 <- lm.fit.p2.step.1.cloglog$aic
  c.cloglog.aic <- c(cloglog.aic.p2.step.1)
  cloglog.dev.p2.step.1 <- lm.fit.p2.step.1.cloglog$deviance
  c.cloglog.dev <- c(cloglog.dev.p2.step.1)
  cloglog.null.dev.p2.step.1 <- lm.fit.p2.step.1.cloglog$null.deviance
  ccloglog.null.deviance <- c(cloglog.null.dev.p2.step.1)
  stat.p1step.cloglog <- cbind(c.cloglog.aic, c.cloglog.dev, ccloglog.null.deviance)
  
  JTpcloglog<-cloglog.fitted.p2.step.1>(1-tau+tuning)
  
  dataJTpcloglog<-as.data.frame(data[JTpcloglog,])
  y.s5s<-dataJTpcloglog[,1]
  x1.s5s<-dataJTpcloglog[,2]
  x2.s5s<-dataJTpcloglog[,3]
  x1s.s5s<-dataJTpcloglog[,4]
  x2s.s5s<-dataJTpcloglog[,5]
  s.s5s<-dataJTpcloglog[,6]
  X.s5s<-cbind(x1.s5s,x2.s5s)
  Z.s5s <- as.matrix(model.matrix(~as.factor(s.s5s)-1))
  n.s5s<-ncol(Z.s5s)
  
  as.matrix.csr(cbind(X.s5s,Z.s5s)) -> sM.s5s
  coefs.step5s <- rq.fit.sfn(sM.s5s,y.s5s,tau=tau)$coef
  
  fit.p2.step.2.cloglog <- rq(y.s5s ~ X.s5s + Z.s5s - 1,tau=tau)
  results.p2step.cloglog <- summary(fit.p2.step.2.cloglog,se="ker")$coef[1:2,]
  residuals.p2step.cloglog <- summary(fit.p2.step.2.cloglog,se="ker")$residuals
  res.cloglog <- c(residuals.p2step.cloglog)
  RMSE.p.2.step.cloglog <- sqrt(mean((res.cloglog)^2))
  bias.cloglog <- c(beta-results.p2step.cloglog[1:2,1])
  stat.p2step.cloglog <- cbind(RMSE.p.2.step.cloglog, bias.cloglog)
  
  #### Generalized additive models #### 
  
  #### Generalized additive models - Logit #### 

  tuning <- (n*t)^(-1/5) * tau
  gam.fit.logit <- gam(delta ~ s(X1) + s(X2) + Z, family=binomial(link = "logit"))
  gam.fitted.logit <- fitted(gam.fit.logit)
  logit.dev.n2.step.1 <- gam.fit.logit$deviance
  clogit.dev.n2.step.1 <- c(logit.dev.n2.step.1)
  logit.null.dev.n2.step.1 <- gam.fit.logit$null.deviance
  clogit.null.dev.n2.step.1 <- c(logit.null.dev.n2.step.1)
  stat.n1step.logit <- cbind(clogit.dev.n2.step.1, clogit.null.dev.n2.step.1)
  
  JTnlogit<-gam.fitted.logit>(1-tau+tuning)
  
  dataJTnlogit<-as.data.frame(data[JTnlogit,])
  y.s1s<-dataJTnlogit[,1]
  x1.s1s<-dataJTnlogit[,2]
  x2.s1s<-dataJTnlogit[,3]
  x1s.s1s<-dataJTnlogit[,4]
  x2s.s1s<-dataJTnlogit[,5]
  s.s1s<-dataJTnlogit[,6]
  X.s1s<-cbind(x1.s1s,x2.s1s)
  Z.s1s <- as.matrix(model.matrix(~as.factor(s.s1s)-1))
  n.s1s<-ncol(Z.s1s)
  
  as.matrix.csr(cbind(X.s1s,Z.s1s)) -> sM.s1s
  coefs.step1s <- rq.fit.sfn(sM.s1s,y.s1s,tau=tau)$coef
  
  fitn2steplogit <- rq(y.s1s ~ X.s1s + Z.s1s - 1,tau=tau)
  results.n2step.logit <- summary(fitn2steplogit,se="ker")$coef[1:2,]
  residuals.n2step.logit <- summary(fitn2steplogit,se="ker")$residuals
  res.n2step.logit <- c(residuals.n2step.logit)
  RMSE.n2step.logit <- sqrt(mean((res.n2step.logit)^2))
  bias.n2step.logit <- c(beta-results.n2step.logit[1:2,1])
  stat.n2step.logit <- cbind(RMSE.n2step.logit, bias.n2step.logit)
  
  #### Generalized additive models - Cloglog #####
  
  gam.fit.cloglog <- gam(delta ~ s(X1) + s(X2) + Z, family=binomial(link = "cloglog"))
  gam.fitted.cloglog <- fitted(gam.fit.cloglog)
  cloglog.dev.n2.step.1 <- gam.fit.cloglog$deviance
  ccloglog.dev.n2.step.1 <- c(cloglog.dev.n2.step.1)
  cloglog.null.dev.n2.step.1 <- gam.fit.cloglog$null.deviance
  ccloglog.null.dev.n2.step.1 <- c(cloglog.null.dev.n2.step.1)
  stat.n1step.cloglog <- cbind(ccloglog.dev.n2.step.1, ccloglog.null.dev.n2.step.1)
  
  JTncloglog<-gam.fitted.cloglog>(1-tau+tuning)
  
  dataJTncloglog<-as.data.frame(data[JTncloglog,])
  y.s13s<-dataJTncloglog[,1]
  x1.s13s<-dataJTncloglog[,2]
  x2.s13s<-dataJTncloglog[,3]
  x1s.s13s<-dataJTncloglog[,4]
  x2s.s13s<-dataJTncloglog[,5]
  s.s13s<-dataJTncloglog[,6]
  X.s13s<-cbind(x1.s13s,x2.s13s)
  Z.s13s <- as.matrix(model.matrix(~as.factor(s.s13s)-1))
  n.s13s<-ncol(Z.s13s)
  
  as.matrix.csr(cbind(X.s13s,Z.s13s)) -> sM.s13s
  coefs.step13s <- rq.fit.sfn(sM.s13s,y.s13s,tau=tau)$coef
  
  fitn2stepcloglog <- rq(y.s13s ~ X.s13s + Z.s13s - 1,tau=tau)
  results.n2step.cloglog <- summary(fitn2stepcloglog,se="ker")$coef[1:2,]
  residuals.n2step.cloglog <- summary(fitn2stepcloglog,se="ker")$residuals
  res.n2step.cloglog <- c(residuals.n2step.cloglog)
  RMSE.n2step.cloglog <- sqrt(mean((res.n2step.cloglog)^2))
  bias.n2step.cloglog <- c(beta-results.n2step.cloglog[1:2,1])
  stat.n2step.cloglog <- cbind(RMSE.n2step.cloglog, bias.n2step.cloglog)
  
#### Three-step estimator #### 
  
  #### Three-step - Logit #### 
  
  lm.fit <- glm(delta ~ X + Z,family=binomial(link = "logit"))
  logit.fitted <- fitted(lm.fit)
  
  logit.aic.p3.step.1 <- lm.fit$aic
  c.logit.p3.aic <- c(logit.aic.p3.step.1)
  logit.dev.p3.step.1 <- lm.fit$deviance
  c.logit.p3.deviance <- c(logit.dev.p3.step.1)
  logit.null.dev.p3.step.1 <- lm.fit$null.deviance
  c.logit.p3.null.deviance <- c(logit.null.dev.p3.step.1)
  stat.p3step.logit <- cbind(c.logit.p3.aic, c.logit.p3.deviance, c.logit.p3.null.deviance)
  
  d <- min(quantile(logit.fitted,.1),const)
  J0 <- logit.fitted > 1-tau + d
  
  dataJ0 <- as.data.frame(data[J0,])
  y.s1   <- dataJ0[,1]
  x1.s1  <- dataJ0[,2]
  x2.s1  <- dataJ0[,3]
  x1s.s1 <- dataJ0[,4]
  x2s.s1 <- dataJ0[,5]
  s.s1   <- dataJ0[,6]
  X.s1   <- cbind(x1.s1,x2.s1)
  Z.s1 <- as.matrix(model.matrix(~as.factor(s.s1)-1))
  n.s1<-ncol(Z.s1)
  
  as.matrix.csr(cbind(X.s1,Z.s1)) -> sM.s1
  coefs.step1 <- rq.fit.sfn(sM.s1,y.s1,tau=tau)$coef
  
  alphas <- coefs.step1[-c(1:2)]
  Zas <- matrix(0,n*t,1)
  s.s1[match(s,s.s1,nomatch = 999999)] -> new.S
  for (k in 1:length(unique(s.s1)))
  {
    J <- unique(s.s1)[k]
    Zas[which(new.S==J),1] <- rep(alphas[k],length(which(new.S==J)))
  }
  fit.rq.fitted <- X[,1:2]%*%coefs.step1[1:2] + Zas
  
  deltan <- quantile(fit.rq.fitted[fit.rq.fitted>0],((1/3)*(n*t)^(-1/3)))
  J1 <- fit.rq.fitted>(Ci+deltan)
  
  dataJ1<-as.data.frame(data[J1,])
  y.s2<-dataJ1[,1]
  x1.s2<-dataJ1[,2]
  x2.s2<-dataJ1[,3]
  x1s.s2<-dataJ1[,4]
  x2s.s2<-dataJ1[,5]
  s.s2<-dataJ1[,6]
  X.s2<-cbind(x1.s2,x2.s2)
  Z.s2 <- as.matrix(model.matrix(~as.factor(s.s2)-1))
  n.s2<-ncol(Z.s2)
  
  as.matrix.csr(cbind(X.s2,Z.s2)) -> sM.s2
  
  fit <- rq(y.s2 ~ X.s2 + Z.s2 - 1,tau=tau)
  results.3step <- summary(fit,se="ker")$coef[1:2,]
  residuals.3step <- summary(fit,se="ker")$residuals
  res.3step <- c(residuals.3step)
  RMSE.3step <- sqrt(mean((res.3step)^2))
  bias.3step <- c(beta-results.3step[1:2,1])
  stat.3step <- cbind(RMSE.3step, bias.3step)
  
  #### Three-step - Cloglog #### 
  
  lm.fitc <- glm(delta ~ X + Z, family=binomial(link = "cloglog"))
  logit.fittedc <- fitted(lm.fitc)
  
  clog.aic.p3.step.1 <- lm.fitc$aic
  c.clog.p3.aic <- c(clog.aic.p3.step.1)
  clog.dev.p3.step.1 <- lm.fitc$deviance
  c.clog.p3.deviance <- c(clog.dev.p3.step.1)
  clog.null.dev.p3.step.1 <- lm.fitc$null.deviance
  c.clog.p3.null.deviance <- c(clog.null.dev.p3.step.1)
  stat.p3step.clog <- cbind(c.clog.p3.aic, c.clog.p3.deviance, c.clog.p3.null.deviance)
  
  d <- min(quantile(logit.fittedc,.1),const)
  J0c <- logit.fittedc > 1-tau + d
  
  dataJ0c <- as.data.frame(data[J0c,])
  y.s1c   <- dataJ0c[,1]
  x1.s1c  <- dataJ0c[,2]
  x2.s1c  <- dataJ0c[,3]
  x1s.s1c <- dataJ0c[,4]
  x2s.s1c <- dataJ0c[,5]
  s.s1c   <- dataJ0c[,6]
  X.s1c   <- cbind(x1.s1c,x2.s1c)
  Z.s1c <- as.matrix(model.matrix(~as.factor(s.s1c)-1))
  n.s1c <-ncol(Z.s1c)
  
  as.matrix.csr(cbind(X.s1c,Z.s1c)) -> sM.s1c
  coefs.step1c <- rq(y.s1c ~ X.s1c + Z.s1c - 1,tau=tau)$coef
  coefs.s1c <- coefs.step1c[1:2]
  
  alphas <- coefs.step1c[-c(1:2)]
  Zasc <- matrix(0,n*t,1)
  s.s1c[match(s,s.s1c,nomatch = 999999)] -> new.Sc
  for (k in 1:length(unique(s.s1c)))
  {
    Jc <- unique(s.s1c)[k]
    Zasc[which(new.Sc==Jc),1] <- rep(alphas[k],length(which(new.Sc==Jc)))
  }
  
  fit.rq.fittedc <- X[,1:2]%*%coefs.step1c[1:2] + Zasc
  
  deltan <- quantile(fit.rq.fittedc[fit.rq.fittedc>0],((1/3)*(n*t)^(-1/3)))
  J1c <- fit.rq.fittedc>(Ci+deltan)
  
  dataJ1c<-as.data.frame(data[J1c,])
  y.s2c<-dataJ1c[,1]
  x1.s2c<-dataJ1c[,2]
  x2.s2c<-dataJ1c[,3]
  x1s.s2c<-dataJ1c[,4]
  x2s.s2c<-dataJ1c[,5]
  s.s2c<-dataJ1c[,6]
  X.s2c<-cbind(x1.s2c,x2.s2c)
  Z.s2c <- as.matrix(model.matrix(~as.factor(s.s2c)-1))
  n.s2c<-ncol(Z.s2c)
  
  as.matrix.csr(cbind(X.s2c,Z.s2c)) -> sM.s2c
  
  fit3stepc <- rq(y.s2c ~ X.s2c + Z.s2c - 1,tau=tau)
  results.3stepc <- summary(fit3stepc,se="ker")$coef[1:2,]
  residuals.3stepc <- summary(fit3stepc,se="ker")$residuals
  res.3stepc <- c(residuals.3stepc)
  RMSE.3stepc <- sqrt(mean((res.3stepc)^2))
  bias.3stepc <- c(beta-results.3stepc[1:2,1])
  stat.3stepc <- cbind(RMSE.3stepc, bias.3stepc)
  
  #### Omniscient and Naive #### 
  
  as.matrix.csr(cbind(X,Z)) -> sM
  coefs.omni <- rq.fit.sfn(sM,ystar,tau=tau)$coef[1:2]
  ccoefs.omni <- c(coefs.omni)
  fitcoef.omni <- rq.fit.sfn(sM,ystar,tau=tau)
  res.omni <- fitcoef.omni[["residuals"]]
  RMSE.omni <- sqrt(mean((res.omni)^2))
  bias.omni <- c(beta-ccoefs.omni)
  stat.omni <- cbind(RMSE.omni, bias.omni)
  
  as.matrix.csr(cbind(X,Z)) -> sM
  coefs.naive <- rq.fit.sfn(sM,y,tau=tau)$coef[1:2]
  ccoefs.naive <- c(coefs.naive)
  fitcoef.naive <- rq.fit.sfn(sM,y,tau=tau)
  res.naive <- fitcoef.naive[["residuals"]]
  RMSE.naive <- sqrt(mean((res.naive)^2))
  bias.naive <- c(beta-ccoefs.naive)
  stat.naive <- cbind(RMSE.naive, bias.naive)
  
}
