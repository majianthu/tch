library(copent)
library(copula)
library(mvtnorm)
library(gofCopula)
library(latex2exp)

tGaussianCopula = function(x){
  corx = cor(x)
  mu1 = mean(x[,1]); mu2 = mean(x[,2])
  std1 = sqrt(var(x[,1])); std2 = sqrt(var(x[,2]))
  u = construct_empirical_copula(x) - 0.5 / dim(x)[1]
  uq1 = qnorm(u[,1],mu1,std1); uq2 = qnorm(u[,2],mu2,std2)
  uq = cbind(uq1,uq2)
  copent(x) - mean( -0.5 * log(det(corx)) - 0.5 * diag(uq %*% (solve(corx) - diag(dim(corx)[1])) %*% t(uq)) )
}

tGumbelCopula = function(x){
  u = construct_empirical_copula(x) - 0.5 / dim(x)[1]
  fit.tau = fitCopula(gumbelCopula(), u, method = "mpl")
  tau = fit.tau@estimate
  gCopula = gumbelCopula(param=tau)
  copent(x) - mean(dCopula(u,gCopula,log = TRUE))
}

tFrankCopula = function(x){
  u = construct_empirical_copula(x) - 0.5 / dim(x)[1]
  fit.tau = fitCopula(frankCopula(), u, method = "mpl")
  tau = fit.tau@estimate
  fCopula = frankCopula(param=tau)
  copent(x) - mean(dCopula(u,fCopula,log = TRUE))
}

tClaytonCopula = function(x){
  u = construct_empirical_copula(x) - 0.5 / dim(x)[1]
  fit.tau = fitCopula(claytonCopula(), u, method = "mpl")
  tau = fit.tau@estimate
  cCopula = claytonCopula(param=tau)
  copent(x) - mean(dCopula(u,cCopula,log = TRUE))
}

M1 = 1
# simulation for testing Gaussian copula
tgc1a = tgc1b = tgc1c = tgc1d = 0
cvm1a = kcvm1a = kks1a = ks1a = piosrn1a = piostn1a = rcq1a = rg1a = rsnb1a = rsnc1a = white1a = 0
cvm1b = kcvm1b = kks1b = ks1b = piosrn1b = piostn1b = rcq1b = rg1b = rsnb1b = rsnc1b = white1b = 0
cvm1c = kcvm1c = kks1c = ks1c = piosrn1c = piostn1c = rcq1c = rg1c = rsnb1c = rsnc1c = white1c = 0
cvm1d = kcvm1d = kks1d = ks1d = piosrn1d = piostn1d = rcq1d = rg1d = rsnb1d = rsnc1d = white1d = 0
for(i in 1:9){
  mu = c(0,0)
  rho = i * 0.1
  sigma = matrix(c(1,rho,rho,1),2,2)
  n = 300
  x1 = rmvnorm(n,mu,sigma)
  tgc1a[i] = tGaussianCopula(x1)
  tgc1b[i] = tGumbelCopula(x1)
  tgc1c[i] = tFrankCopula(x1)
  tgc1d[i] = tClaytonCopula(x1)
  
  cvm1a[i] = gofCvM("normal",x1,M=M1)$normal$res.tests[2]
  kcvm1a[i] = gofKendallCvM("normal",x1,M=M1)$normal$res.tests[2]
  kks1a[i] = gofKendallKS("normal",x1,M=M1)$normal$res.tests[2]
  ks1a[i] = gofKS("normal",x1,M=M1)$normal$res.tests[2]
  piosrn1a[i] = gofPIOSRn("normal",x1,M=M1)$normal$res.tests[2]
  piostn1a[i] = gofPIOSTn("normal",x1,M=M1)$normal$res.tests[2]
  rcq1a[i] = gofRosenblattChisq("normal",x1,M=M1)$normal$res.tests[2]
  rg1a[i] = gofRosenblattGamma("normal",x1,M=M1)$normal$res.tests[2]
  rsnb1a[i] = gofRosenblattSnB("normal",x1,M=M1)$normal$res.tests[2]
  rsnc1a[i] = gofRosenblattSnC("normal",x1,M=M1)$normal$res.tests[2]
  white1a[i] = gofWhite("normal",x1,M=M1)$normal$res.tests[2]

  cvm1b[i] = gofCvM("gumbel",x1,M=M1)$gumbel$res.tests[2]
  kcvm1b[i] = gofKendallCvM("gumbel",x1,M=M1)$gumbel$res.tests[2]
  kks1b[i] = gofKendallKS("gumbel",x1,M=M1)$gumbel$res.tests[2]
  ks1b[i] = gofKS("gumbel",x1,M=M1)$gumbel$res.tests[2]
  piosrn1b[i] = gofPIOSRn("gumbel",x1,M=M1)$gumbel$res.tests[2]
  piostn1b[i] = gofPIOSTn("gumbel",x1,M=M1)$gumbel$res.tests[2]
  rcq1b[i] = gofRosenblattChisq("gumbel",x1,M=M1)$gumbel$res.tests[2]
  rg1b[i] = gofRosenblattGamma("gumbel",x1,M=M1)$gumbel$res.tests[2]
  rsnb1b[i] = gofRosenblattSnB("gumbel",x1,M=M1)$gumbel$res.tests[2]
  rsnc1b[i] = gofRosenblattSnC("gumbel",x1,M=M1)$gumbel$res.tests[2]
  white1b[i] = gofWhite("gumbel",x1,M=M1)$gumbel$res.tests[2]

  cvm1c[i] = gofCvM("frank",x1,M=M1)$frank$res.tests[2]
  kcvm1c[i] = gofKendallCvM("frank",x1,M=M1)$frank$res.tests[2]
  kks1c[i] = gofKendallKS("frank",x1,M=M1)$frank$res.tests[2]
  ks1c[i] = gofKS("frank",x1,M=M1)$frank$res.tests[2]
  piosrn1c[i] = gofPIOSRn("frank",x1,M=M1)$frank$res.tests[2]
  piostn1c[i] = gofPIOSTn("frank",x1,M=M1)$frank$res.tests[2]
  rcq1c[i] = gofRosenblattChisq("frank",x1,M=M1)$frank$res.tests[2]
  rg1c[i] = gofRosenblattGamma("frank",x1,M=M1)$frank$res.tests[2]
  rsnb1c[i] = gofRosenblattSnB("frank",x1,M=M1)$frank$res.tests[2]
  rsnc1c[i] = gofRosenblattSnC("frank",x1,M=M1)$frank$res.tests[2]
  white1c[i] = gofWhite("frank",x1,M=M1)$frank$res.tests[2]

  cvm1d[i] = gofCvM("clayton",x1,M=M1)$clayton$res.tests[2]
  kcvm1d[i] = gofKendallCvM("clayton",x1,M=M1)$clayton$res.tests[2]
  kks1d[i] = gofKendallKS("clayton",x1,M=M1)$clayton$res.tests[2]
  ks1d[i] = gofKS("clayton",x1,M=M1)$clayton$res.tests[2]
  piosrn1d[i] = gofPIOSRn("clayton",x1,M=M1)$clayton$res.tests[2]
  piostn1d[i] = gofPIOSTn("clayton",x1,M=M1)$clayton$res.tests[2]
  rcq1d[i] = gofRosenblattChisq("clayton",x1,M=M1)$clayton$res.tests[2]
  rg1d[i] = gofRosenblattGamma("clayton",x1,M=M1)$clayton$res.tests[2]
  rsnb1d[i] = gofRosenblattSnB("clayton",x1,M=M1)$clayton$res.tests[2]
  rsnc1d[i] = gofRosenblattSnC("clayton",x1,M=M1)$clayton$res.tests[2]
  white1d[i] = gofWhite("clayton",x1,M=M1)$clayton$res.tests[2]
}

x11(width = 12, height = 7)
par(mfrow = c(3,4))
xlab1 = TeX(r'($\rho$)') 
rho1 = seq(0.1,0.9,by = 0.1)
cex1 = 0.7

plot(rho1,tgc1a, xlab = xlab1, ylab = 'statistic', main = "CE", ylim = c(min(tgc1a,tgc1b,tgc1c,tgc1d),max(tgc1a,tgc1b,tgc1c,tgc1d)), col = 'red', pch = 1);
lines(rho1,tgc1a, col = 'red')
points(rho1,tgc1b, col = 'blue', pch = 2); 
lines(rho1,tgc1b, col = 'blue')
points(rho1,tgc1c, col = 'green', pch = 3); 
lines(rho1,tgc1c, col = 'green')
points(rho1,tgc1d, col = 'black', pch = 4); 
lines(rho1,tgc1d, col = 'black')
legend(0.1,max(tgc1a,tgc1b,tgc1d,tgc1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(rho1,cvm1a, xlab = xlab1, ylab = 'statistic', main = "CvM", ylim = c(min(cvm1a,cvm1b,cvm1c,cvm1d),max(cvm1a,cvm1b,cvm1c,cvm1d)), col = 'red', pch = 1);
lines(rho1,cvm1a, col = 'red')
points(rho1,cvm1b, col = 'blue', pch = 2); 
lines(rho1,cvm1b, col = 'blue')
points(rho1,cvm1c, col = 'green', pch = 3); 
lines(rho1,cvm1c, col = 'green')
points(rho1,cvm1d, col = 'black', pch = 4); 
lines(rho1,cvm1d, col = 'black')
legend(0.1,max(cvm1a,cvm1b,cvm1d,cvm1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(rho1,kcvm1a, xlab = xlab1, ylab = 'statistic', main = "KendallCvM", ylim = c(min(kcvm1a,kcvm1b,kcvm1c,kcvm1d),max(kcvm1a,kcvm1b,kcvm1c,kcvm1d)), col = 'red', pch = 1);
lines(rho1,kcvm1a, col = 'red')
points(rho1,kcvm1b, col = 'blue', pch = 2); 
lines(rho1,kcvm1b, col = 'blue')
points(rho1,kcvm1c, col = 'green', pch = 3); 
lines(rho1,kcvm1c, col = 'green')
points(rho1,kcvm1d, col = 'black', pch = 4); 
lines(rho1,kcvm1d, col = 'black')
legend(0.1,max(kcvm1a,kcvm1b,kcvm1d,kcvm1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(rho1,kks1a, xlab = xlab1, ylab = 'statistic', main = "KendallKS", ylim = c(min(kks1a,kks1b,kks1c,kks1d),max(kks1a,kks1b,kks1c,kks1d)), col = 'red', pch = 1);
lines(rho1,kks1a, col = 'red')
points(rho1,kks1b, col = 'blue', pch = 2); 
lines(rho1,kks1b, col = 'blue')
points(rho1,kks1c, col = 'green', pch = 3); 
lines(rho1,kks1c, col = 'green')
points(rho1,kks1d, col = 'black', pch = 4); 
lines(rho1,kks1d, col = 'black')
legend(0.1,max(kks1a,kks1b,kks1d,kks1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(rho1,ks1a, xlab = xlab1, ylab = 'statistic', main = "KS", ylim = c(min(ks1a,ks1b,ks1c,ks1d),max(ks1a,ks1b,ks1c,ks1d)), col = 'red', pch = 1);
lines(rho1,ks1a, col = 'red')
points(rho1,ks1b, col = 'blue', pch = 2); 
lines(rho1,ks1b, col = 'blue')
points(rho1,ks1c, col = 'green', pch = 3); 
lines(rho1,ks1c, col = 'green')
points(rho1,ks1d, col = 'black', pch = 4); 
lines(rho1,ks1d, col = 'black')
legend(0.1,max(ks1a,ks1b,ks1d,ks1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(rho1,piosrn1a, xlab = xlab1, ylab = 'statistic', main = "PIOSRn", ylim = c(min(piosrn1a,piosrn1b,piosrn1c,piosrn1d),max(piosrn1a,piosrn1b,piosrn1c,piosrn1d)), col = 'red', pch = 1);
lines(rho1,piosrn1a, col = 'red')
points(rho1,piosrn1b, col = 'blue', pch = 2); 
lines(rho1,piosrn1b, col = 'blue')
points(rho1,piosrn1c, col = 'green', pch = 3); 
lines(rho1,piosrn1c, col = 'green')
points(rho1,piosrn1d, col = 'black', pch = 4); 
lines(rho1,piosrn1d, col = 'black')
legend(0.1,max(piosrn1a,piosrn1b,piosrn1d,piosrn1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(rho1,piostn1a, xlab = xlab1, ylab = 'statistic', main = "PIOSTn", ylim = c(min(piostn1a,piostn1b,piostn1c,piostn1d),max(piostn1a,piostn1b,piostn1c,piostn1d)), col = 'red', pch = 1);
lines(rho1,piostn1a, col = 'red')
points(rho1,piostn1b, col = 'blue', pch = 2); 
lines(rho1,piostn1b, col = 'blue')
points(rho1,piostn1c, col = 'green', pch = 3); 
lines(rho1,piostn1c, col = 'green')
points(rho1,piostn1d, col = 'black', pch = 4); 
lines(rho1,piostn1d, col = 'black')
legend(0.1,max(piostn1a,piostn1b,piostn1d,piostn1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(rho1,rcq1a, xlab = xlab1, ylab = 'statistic', main = "RosenblattChisq", ylim = c(min(rcq1a,rcq1b,rcq1c,rcq1d),max(rcq1a,rcq1b,rcq1c,rcq1d)), col = 'red', pch = 1);
lines(rho1,rcq1a, col = 'red')
points(rho1,rcq1b, col = 'blue', pch = 2); 
lines(rho1,rcq1b, col = 'blue')
points(rho1,rcq1c, col = 'green', pch = 3); 
lines(rho1,rcq1c, col = 'green')
points(rho1,rcq1d, col = 'black', pch = 4); 
lines(rho1,rcq1d, col = 'black')
legend(0.1,max(rcq1a,rcq1b,rcq1d,rcq1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(rho1,rg1a, xlab = xlab1, ylab = 'statistic', main = "RosenblattGamma", ylim = c(min(rg1a,rg1b,rg1c,rg1d),max(rg1a,rg1b,rg1c,rg1d)), col = 'red', pch = 1);
lines(rho1,rg1a, col = 'red')
points(rho1,rg1b, col = 'blue', pch = 2); 
lines(rho1,rg1b, col = 'blue')
points(rho1,rg1c, col = 'green', pch = 3); 
lines(rho1,rg1c, col = 'green')
points(rho1,rg1d, col = 'black', pch = 4); 
lines(rho1,rg1d, col = 'black')
legend(0.1,max(rg1a,rg1b,rg1d,rg1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(rho1,rsnb1a, xlab = xlab1, ylab = 'statistic', main = "RosenblattSnB", ylim = c(min(rsnb1a,rsnb1b,rsnb1c,rsnb1d),max(rsnb1a,rsnb1b,rsnb1c,rsnb1d)), col = 'red', pch = 1);
lines(rho1,rsnb1a, col = 'red')
points(rho1,rsnb1b, col = 'blue', pch = 2); 
lines(rho1,rsnb1b, col = 'blue')
points(rho1,rsnb1c, col = 'green', pch = 3); 
lines(rho1,rsnb1c, col = 'green')
points(rho1,rsnb1d, col = 'black', pch = 4); 
lines(rho1,rsnb1d, col = 'black')
legend(0.1,max(rsnb1a,rsnb1b,rsnb1d,rsnb1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(rho1,rsnc1a, xlab = xlab1, ylab = 'statistic', main = "RosenblattSnC", ylim = c(min(rsnc1a,rsnc1b,rsnc1c,rsnc1d),max(rsnc1a,rsnc1b,rsnc1c,rsnc1d)), col = 'red', pch = 1);
lines(rho1,rsnc1a, col = 'red')
points(rho1,rsnc1b, col = 'blue', pch = 2); 
lines(rho1,rsnc1b, col = 'blue')
points(rho1,rsnc1c, col = 'green', pch = 3); 
lines(rho1,rsnc1c, col = 'green')
points(rho1,rsnc1d, col = 'black', pch = 4); 
lines(rho1,rsnc1d, col = 'black')
legend(0.1,max(rsnc1a,rsnc1b,rsnc1d,rsnc1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(rho1,white1a, xlab = xlab1, ylab = 'statistic', main = "White", ylim = c(min(white1a,white1b,white1c,white1d),max(white1a,white1b,white1c,white1d)), col = 'red', pch = 1);
lines(rho1,white1a, col = 'red')
points(rho1,white1b, col = 'blue', pch = 2); 
lines(rho1,white1b, col = 'blue')
points(rho1,white1c, col = 'green', pch = 3); 
lines(rho1,white1c, col = 'green')
points(rho1,white1d, col = 'black', pch = 4); 
lines(rho1,white1d, col = 'black')
legend(0.1,max(white1a,white1b,white1c,white1d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)


# simulation for testing Archimedean copula
tgc2a = tgc2b = tgc2c = tgc2d = 0
cvm2a = kcvm2a = kks2a = ks2a = piosrn2a = piostn2a = rcq2a = rg2a = rsnb2a = rsnc2a = white2a = 0
cvm2b = kcvm2b = kks2b = ks2b = piosrn2b = piostn2b = rcq2b = rg2b = rsnb2b = rsnc2b = white2b = 0
cvm2c = kcvm2c = kks2c = ks2c = piosrn2c = piostn2c = rcq2c = rg2c = rsnb2c = rsnc2c = white2c = 0
cvm2d = kcvm2d = kks2d = ks2d = piosrn2d = piostn2d = rcq2d = rg2d = rsnb2d = rsnc2d = white2d = 0

M1 = 1
for(t in 2:10){
  # cop2 = gumbelCopula(t)
  # cop2 = frankCopula(t)
  cop2 = claytonCopula(t)
  mvdc2 <- mvdc(cop2, c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = 2)))
  n = 300
  x2 = rMvdc(n,mvdc2)
  tgc2a[t-1] = tGaussianCopula(x2)
  tgc2b[t-1] = tGumbelCopula(x2)
  tgc2c[t-1] = tFrankCopula(x2)
  tgc2d[t-1] = tClaytonCopula(x2)
  
  cvm2a[t-1] = gofCvM("normal",x2,M=M1)$normal$res.tests[2]
  kcvm2a[t-1] = gofKendallCvM("normal",x2,M=M1)$normal$res.tests[2]
  kks2a[t-1] = gofKendallKS("normal",x2,M=M1)$normal$res.tests[2]
  ks2a[t-1] = gofKS("normal",x2,M=M1)$normal$res.tests[2]
  piosrn2a[t-1] = gofPIOSRn("normal",x2,M=M1)$normal$res.tests[2]
  piostn2a[t-1] = gofPIOSTn("normal",x2,M=M1)$normal$res.tests[2]
  rcq2a[t-1] = gofRosenblattChisq("normal",x2,M=M1)$normal$res.tests[2]
  rg2a[t-1] = gofRosenblattGamma("normal",x2,M=M1)$normal$res.tests[2]
  rsnb2a[t-1] = gofRosenblattSnB("normal",x2,M=M1)$normal$res.tests[2]
  rsnc2a[t-1] = gofRosenblattSnC("normal",x2,M=M1)$normal$res.tests[2]
  white2a[t-1] = gofWhite("normal",x2,M=M1)$normal$res.tests[2]
  
  cvm2b[t-1] = gofCvM("gumbel",x2,M=M1)$gumbel$res.tests[2]
  kcvm2b[t-1] = gofKendallCvM("gumbel",x2,M=M1)$gumbel$res.tests[2]
  kks2b[t-1] = gofKendallKS("gumbel",x2,M=M1)$gumbel$res.tests[2]
  ks2b[t-1] = gofKS("gumbel",x2,M=M1)$gumbel$res.tests[2]
  piosrn2b[t-1] = gofPIOSRn("gumbel",x2,M=M1)$gumbel$res.tests[2]
  piostn2b[t-1] = gofPIOSTn("gumbel",x2,M=M1)$gumbel$res.tests[2]
  rcq2b[t-1] = gofRosenblattChisq("gumbel",x2,M=M1)$gumbel$res.tests[2]
  rg2b[t-1] = gofRosenblattGamma("gumbel",x2,M=M1)$gumbel$res.tests[2]
  rsnb2b[t-1] = gofRosenblattSnB("gumbel",x2,M=M1)$gumbel$res.tests[2]
  rsnc2b[t-1] = gofRosenblattSnC("gumbel",x2,M=M1)$gumbel$res.tests[2]
  white2b[t-1] = gofWhite("gumbel",x2,M=M1)$gumbel$res.tests[2]
  
  cvm2c[t-1] = gofCvM("frank",x2,M=M1)$frank$res.tests[2]
  kcvm2c[t-1] = gofKendallCvM("frank",x2,M=M1)$frank$res.tests[2]
  kks2c[t-1] = gofKendallKS("frank",x2,M=M1)$frank$res.tests[2]
  ks2c[t-1] = gofKS("frank",x2,M=M1)$frank$res.tests[2]
  piosrn2c[t-1] = gofPIOSRn("frank",x2,M=M1)$frank$res.tests[2]
  piostn2c[t-1] = gofPIOSTn("frank",x2,M=M1)$frank$res.tests[2]
  rcq2c[t-1] = gofRosenblattChisq("frank",x2,M=M1)$frank$res.tests[2]
  rg2c[t-1] = gofRosenblattGamma("frank",x2,M=M1)$frank$res.tests[2]
  rsnb2c[t-1] = gofRosenblattSnB("frank",x2,M=M1)$frank$res.tests[2]
  rsnc2c[t-1] = gofRosenblattSnC("frank",x2,M=M1)$frank$res.tests[2]
  white2c[t-1] = gofWhite("frank",x2,M=M1)$frank$res.tests[2]
  
  cvm2d[t-1] = gofCvM("clayton",x2,M=M1)$clayton$res.tests[2]
  kcvm2d[t-1] = gofKendallCvM("clayton",x2,M=M1)$clayton$res.tests[2]
  kks2d[t-1] = gofKendallKS("clayton",x2,M=M1)$clayton$res.tests[2]
  ks2d[t-1] = gofKS("clayton",x2,M=M1)$clayton$res.tests[2]
  piosrn2d[t-1] = gofPIOSRn("clayton",x2,M=M1)$clayton$res.tests[2]
  piostn2d[t-1] = gofPIOSTn("clayton",x2,M=M1)$clayton$res.tests[2]
  rcq2d[t-1] = gofRosenblattChisq("clayton",x2,M=M1)$clayton$res.tests[2]
  rg2d[t-1] = gofRosenblattGamma("clayton",x2,M=M1)$clayton$res.tests[2]
  rsnb2d[t-1] = gofRosenblattSnB("clayton",x2,M=M1)$clayton$res.tests[2]
  rsnc2d[t-1] = gofRosenblattSnC("clayton",x2,M=M1)$clayton$res.tests[2]
  white2d[t-1] = gofWhite("clayton",x2,M=M1)$clayton$res.tests[2]
}

x11(width = 12, height = 7)
par(mfrow = c(3,4))
xlab2 = TeX(r'($\alpha$)') 
alpha2 = 2:10
cex1 = 0.7

plot(alpha2,tgc2a, xlab = xlab2, ylab = 'statistic', main = "CE", ylim = c(min(tgc2a,tgc2b,tgc2c,tgc2d),max(tgc2a,tgc2b,tgc2c,tgc2d)), col = 'red', pch = 1);
lines(alpha2,tgc2a, col = 'red')
points(alpha2,tgc2b, col = 'blue', pch = 2); 
lines(alpha2,tgc2b, col = 'blue')
points(alpha2,tgc2c, col = 'green', pch = 3); 
lines(alpha2,tgc2c, col = 'green')
points(alpha2,tgc2d, col = 'black', pch = 4); 
lines(alpha2,tgc2d, col = 'black')
legend(2,max(tgc2a,tgc2b,tgc2d,tgc2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(alpha2,cvm2a, xlab = xlab2, ylab = 'statistic', main = "CvM", ylim = c(min(cvm2a,cvm2b,cvm2c,cvm2d),max(cvm2a,cvm2b,cvm2c,cvm2d)), col = 'red', pch = 1);
lines(alpha2,cvm2a, col = 'red')
points(alpha2,cvm2b, col = 'blue', pch = 2); 
lines(alpha2,cvm2b, col = 'blue')
points(alpha2,cvm2c, col = 'green', pch = 3); 
lines(alpha2,cvm2c, col = 'green')
points(alpha2,cvm2d, col = 'black', pch = 4); 
lines(alpha2,cvm2d, col = 'black')
legend(2,max(cvm2a,cvm2b,cvm2d,cvm2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(alpha2,kcvm2a, xlab = xlab2, ylab = 'statistic', main = "KendallCvM", ylim = c(min(kcvm2a,kcvm2b,kcvm2c,kcvm2d),max(kcvm2a,kcvm2b,kcvm2c,kcvm2d)), col = 'red', pch = 1);
lines(alpha2,kcvm2a, col = 'red')
points(alpha2,kcvm2b, col = 'blue', pch = 2); 
lines(alpha2,kcvm2b, col = 'blue')
points(alpha2,kcvm2c, col = 'green', pch = 3); 
lines(alpha2,kcvm2c, col = 'green')
points(alpha2,kcvm2d, col = 'black', pch = 4); 
lines(alpha2,kcvm2d, col = 'black')
legend(2,max(kcvm2a,kcvm2b,kcvm2d,kcvm2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(alpha2,kks2a, xlab = xlab2, ylab = 'statistic', main = "KendallKS", ylim = c(min(kks2a,kks2b,kks2c,kks2d),max(kks2a,kks2b,kks2c,kks2d)), col = 'red', pch = 1);
lines(alpha2,kks2a, col = 'red')
points(alpha2,kks2b, col = 'blue', pch = 2); 
lines(alpha2,kks2b, col = 'blue')
points(alpha2,kks2c, col = 'green', pch = 3); 
lines(alpha2,kks2c, col = 'green')
points(alpha2,kks2d, col = 'black', pch = 4); 
lines(alpha2,kks2d, col = 'black')
legend(2,max(kks2a,kks2b,kks2d,kks2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(alpha2,ks2a, xlab = xlab2, ylab = 'statistic', main = "KS", ylim = c(min(ks2a,ks2b,ks2c,ks2d),max(ks2a,ks2b,ks2c,ks2d)), col = 'red', pch = 1);
lines(alpha2,ks2a, col = 'red')
points(alpha2,ks2b, col = 'blue', pch = 2); 
lines(alpha2,ks2b, col = 'blue')
points(alpha2,ks2c, col = 'green', pch = 3); 
lines(alpha2,ks2c, col = 'green')
points(alpha2,ks2d, col = 'black', pch = 4); 
lines(alpha2,ks2d, col = 'black')
legend(2,max(ks2a,ks2b,ks2d,ks2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(alpha2,piosrn2a, xlab = xlab2, ylab = 'statistic', main = "PIOSRn", ylim = c(min(piosrn2a,piosrn2b,piosrn2c,piosrn2d),max(piosrn2a,piosrn2b,piosrn2c,piosrn2d)), col = 'red', pch = 1);
lines(alpha2,piosrn2a, col = 'red')
points(alpha2,piosrn2b, col = 'blue', pch = 2); 
lines(alpha2,piosrn2b, col = 'blue')
points(alpha2,piosrn2c, col = 'green', pch = 3); 
lines(alpha2,piosrn2c, col = 'green')
points(alpha2,piosrn2d, col = 'black', pch = 4); 
lines(alpha2,piosrn2d, col = 'black')
legend(2,max(piosrn2a,piosrn2b,piosrn2d,piosrn2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(alpha2,piostn2a, xlab = xlab2, ylab = 'statistic', main = "PIOSTn", ylim = c(min(piostn2a,piostn2b,piostn2c,piostn2d),max(piostn2a,piostn2b,piostn2c,piostn2d)), col = 'red', pch = 1);
lines(alpha2,piostn2a, col = 'red')
points(alpha2,piostn2b, col = 'blue', pch = 2); 
lines(alpha2,piostn2b, col = 'blue')
points(alpha2,piostn2c, col = 'green', pch = 3); 
lines(alpha2,piostn2c, col = 'green')
points(alpha2,piostn2d, col = 'black', pch = 4); 
lines(alpha2,piostn2d, col = 'black')
legend(2,max(piostn2a,piostn2b,piostn2d,piostn2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(alpha2,rcq2a, xlab = xlab2, ylab = 'statistic', main = "RosenblattChisq", ylim = c(min(rcq2a,rcq2b,rcq2c,rcq2d),max(rcq2a,rcq2b,rcq2c,rcq2d)), col = 'red', pch = 1);
lines(alpha2,rcq2a, col = 'red')
points(alpha2,rcq2b, col = 'blue', pch = 2); 
lines(alpha2,rcq2b, col = 'blue')
points(alpha2,rcq2c, col = 'green', pch = 3); 
lines(alpha2,rcq2c, col = 'green')
points(alpha2,rcq2d, col = 'black', pch = 4); 
lines(alpha2,rcq2d, col = 'black')
legend(2,max(rcq2a,rcq2b,rcq2d,rcq2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(alpha2,rg2a, xlab = xlab2, ylab = 'statistic', main = "RosenblattGamma", ylim = c(min(rg2a,rg2b,rg2c,rg2d),max(rg2a,rg2b,rg2c,rg2d)), col = 'red', pch = 1);
lines(alpha2,rg2a, col = 'red')
points(alpha2,rg2b, col = 'blue', pch = 2); 
lines(alpha2,rg2b, col = 'blue')
points(alpha2,rg2c, col = 'green', pch = 3); 
lines(alpha2,rg2c, col = 'green')
points(alpha2,rg2d, col = 'black', pch = 4); 
lines(alpha2,rg2d, col = 'black')
legend(2,max(rg2a,rg2b,rg2d,rg2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(alpha2,rsnb2a, xlab = xlab2, ylab = 'statistic', main = "RosenblattSnB", ylim = c(min(rsnb2a,rsnb2b,rsnb2c,rsnb2d),max(rsnb2a,rsnb2b,rsnb2c,rsnb2d)), col = 'red', pch = 1);
lines(alpha2,rsnb2a, col = 'red')
points(alpha2,rsnb2b, col = 'blue', pch = 2); 
lines(alpha2,rsnb2b, col = 'blue')
points(alpha2,rsnb2c, col = 'green', pch = 3); 
lines(alpha2,rsnb2c, col = 'green')
points(alpha2,rsnb2d, col = 'black', pch = 4); 
lines(alpha2,rsnb2d, col = 'black')
legend(2,max(rsnb2a,rsnb2b,rsnb2d,rsnb2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(alpha2,rsnc2a, xlab = xlab2, ylab = 'statistic', main = "RosenblattSnC", ylim = c(min(rsnc2a,rsnc2b,rsnc2c,rsnc2d),max(rsnc2a,rsnc2b,rsnc2c,rsnc2d)), col = 'red', pch = 1);
lines(alpha2,rsnc2a, col = 'red')
points(alpha2,rsnc2b, col = 'blue', pch = 2); 
lines(alpha2,rsnc2b, col = 'blue')
points(alpha2,rsnc2c, col = 'green', pch = 3); 
lines(alpha2,rsnc2c, col = 'green')
points(alpha2,rsnc2d, col = 'black', pch = 4); 
lines(alpha2,rsnc2d, col = 'black')
legend(2,max(rsnc2a,rsnc2b,rsnc2d,rsnc2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)

plot(alpha2,white2a, xlab = xlab2, ylab = 'statistic', main = "White", ylim = c(min(white2a,white2b,white2c,white2d),max(white2a,white2b,white2c,white2d)), col = 'red', pch = 1);
lines(alpha2,white2a, col = 'red')
points(alpha2,white2b, col = 'blue', pch = 2); 
lines(alpha2,white2b, col = 'blue')
points(alpha2,white2c, col = 'green', pch = 3); 
lines(alpha2,white2c, col = 'green')
points(alpha2,white2d, col = 'black', pch = 4); 
lines(alpha2,white2d, col = 'black')
legend(2,max(white2a,white2b,white2c,white2d),legend = c("Gaussian","Gumbel","Frank","Clayton"), 
       col = c("red","blue","green","black"), pch = 1:4, cex = cex1)
