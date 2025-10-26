library(copent)
library(copula)
library(mvtnorm)
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

tgc1a = tgc1b = 0
for(i in 1:9){
  mu = c(0,0)
  rho = i * 0.1
  sigma = matrix(c(1,rho,rho,1),2,2)
  n = 300
  x1 = rmvnorm(n,mu,sigma)
  tgc1a[i] = tGaussianCopula(x1)
  tgc1b[i] = tGumbelCopula(x1)
}

x11()
xlab1 = TeX(r'($\rho$)') 
rho1 = seq(0.1,0.9,by = 0.1)
plot(rho1,tgc1a, xlab = xlab1, ylab = 'statistic', ylim = c(min(tgc1a,tgc1b),max(tgc1a,tgc1b)), col = 'red', pch = 1);
lines(rho1,tgc1a, col = 'red')
points(rho1,tgc1b, col = 'blue', pch = 2); 
lines(rho1,tgc1b, col = 'blue')
legend(0.1,max(tgc1a,tgc1b),legend = c("Gaussian","Gumbel"), col = c("red","blue"), pch = 1:2)

tgc2a = tgc2b = 0
for(t in 2:10){
  cop2 = gumbelCopula(t)
  mvdc2 <- mvdc(cop2, c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = 2)))
  n = 300
  x2 = rMvdc(n,mvdc2)
  tgc2a[t-1] = tGaussianCopula(x2)
  tgc2b[t-1] = tGumbelCopula(x2)
}

x11()
xlab2 = TeX(r'($\alpha$)') 
alpha1 = 2:10
plot(alpha1,tgc2a, xlab = xlab2, ylab = 'statistic', ylim = c(min(tgc2a,tgc2b),max(tgc2a,tgc2b)), col = 'red', pch =1);
lines(alpha1,tgc2a, col = 'red')
points(alpha1,tgc2b, col = 'blue', pch = 2); 
lines(alpha1,tgc2b, col = 'blue')
legend(2,max(tgc2a,tgc2b),legend = c("Gaussian","Gumbel"), col = c("red","blue"), pch = 1:2)
