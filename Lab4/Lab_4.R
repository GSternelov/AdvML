library(ggplot2)
library(dlm)
library(gridExtra)
library(mvtnorm)
# Assignment 1

# a) - Simulate the model
ssm1 <- function(C=matrix(0), A=matrix(1), B=matrix(0), u=matrix(0), x0=0,
                 omegaE=0, omegaV=0, TimeS=0) {
  xt <- 0
  yt <- 0
  vt <- rnorm(1, 0, omegaV)
  xt[1] <- A%*%x0 + vt
  et <- rnorm(1, 0, omegaE)
  yt[1] <- C%*%xt + B%*%u + et
  for(i in 2:TimeS){
    vt <- rnorm(1, 0, omegaV)
    xt[i] <- A%*%xt[i-1] + vt
    et <- rnorm(1, 0, omegaE)
    yt[i] <- C%*%xt[i]+ B%*%u + et
  }
  res <- data.frame(Val=c(xt,yt), t= rep(1:TimeS,2), 
                    Value=rep(c("State","Observation"), each=TimeS))
  return(res)
}
set.seed(1234)
resA <- ssm1(C = matrix(0.93), x0 = 0, omegaE = sqrt(0.5), omegaV = sqrt(0.1), TimeS = 100)
plotA <- ggplot(resA, aes(x=t, y=Val, col=Value)) + geom_line(size=1) + theme_classic() + 
  ylim(-4,4) + scale_color_manual(values=c("skyblue", "royalblue")) 
plotA + labs(caption = "Figure 1")
KalmanF <- function(yt, C=matrix(0), A=matrix(1),
                    B=matrix(0), u=matrix(0), omegaE=0, omegaV=0, TimeS=0){
  KalmanMat <- matrix(c(0,0), ncol=2, nrow=101)
  for(i in 1:TimeS){
    # Prediction update
    muBar_t <- A %*%  KalmanMat[i,1] + B%*%u
    sigmaBar_t <- A%*%KalmanMat[i,2]%*%t(A) + omegaV
    # Measurement update
    K_t <- sigmaBar_t%*%t(C) %*% solve(C%*%sigmaBar_t%*%t(C) + omegaE)
    KalmanMat[i+1,1] <- muBar_t + K_t %*% (yt[i] - C%*%muBar_t)
    KalmanMat[i+1,2] <- (diag(ncol(K_t))-K_t%*%C) %*% sigmaBar_t
  }
  KalmanFrame <- data.frame(KalmanMat[-1, ], Time=1:TimeS)
  return(KalmanFrame)
}

resB <- KalmanF(yt = resA[101:200,1], C = matrix(0.93), omegaE = 0.5, omegaV = 0.1, TimeS=100)
resB_1 <- data.frame(Val=c(resB[,1],resB[,1]+1.96*sqrt(resB[,2]),resB[,1]-1.96*sqrt(resB[,2])), Time=rep(1:100, 3), Value=rep(c("Kalman mean", "Upper", "Lower"), each=100))
                     
plotA + geom_line(data=resB_1, aes(x=Time, y=Val, col=Value), size=1) +ylim(-5,4)+
  scale_colour_manual(values=c("red3", "darkorange","skyblue","royalblue","darkorange"))+
  labs(x="Time", title="Kalman filter estimates with 0.95 probability intervals", caption="Figure 2")
filteringB <- dlmFilter(y=resA[101:200,1], dlm(m0=0, V=0.5, W=0.1, C0=10, FF=0.93, GG=1) )

filterB <- data.frame(y=(filteringB$m)[-1])
filterVar <- unlist(dlmSvd2var(u = filteringB$U.C, d = filteringB$D.C))
filterB[101:300,1] <- c(filterB$y + 1.96*sqrt(filterVar[-1]), 
                        filterB$y - 1.96*sqrt(filterVar[-1]))
filterB$Time <- rep(1:100, 3)
filterB$Value <- rep(c("dlm Kalman mean", "dlm Upper", "dlm Lower"), each=100)

ggplot(resB_1, aes(x=Time, y=Val, col=Value,linetype=Value)) + geom_line(size=1) + theme_classic()+ ylim(-5,4)+ geom_line(data=filterB, aes(x=Time, y=y, col=Value, linetype=Value),size=1) +scale_colour_manual(values=c("purple","seagreen","seagreen","red3","darkorange","darkorange"))+scale_linetype_manual(values=c("dashed","dashed","dashed","solid","solid","solid"))+labs(title="Comparsion between dlmFilter and KalmanF function", caption="Figure 3")
PredFunc <- function(C=matrix(0), A=matrix(1), B=matrix(0), u=matrix(0), x0=0,
                     omegaE=0, omegaV=0, TimeS=0){
  xt <- 0
  yt <- 0
  vt <- rnorm(1, 0, omegaV)
  xt[1] <- x0 + vt
  et <- rnorm(1, 0, omegaE)
  yt[1] <- C%*%xt + B%*%u + et
  for(i in 2:TimeS){
    vt <- rnorm(1, 0, omegaV)
    xt[i] <- A%*%xt[i-1] + vt
    et <- rnorm(1, 0, omegaE)
    yt[i] <- C%*%xt[i]+ B%*%u + et
  }
  res <- data.frame(x=c(xt, y=yt))
  return(res)
}

PredC <- as.data.frame(matrix(nrow=10,ncol=10000))
for(i in 1:10000){
  PredC[,i] <- PredFunc(C = matrix(0.93), x0 = rnorm(1, mean=resB$X1[100], 
                        sd=sqrt(resB$X2[100])),omegaE = sqrt(0.5), omegaV = sqrt(0.1),
                        TimeS = 5)
}
quanTs <- t(apply(PredC, 1, quantile, probs=c(0.025, 0.975)))
quanTsFrame <- data.frame(x=c(quanTs[,1], quanTs[,2]), Time= 100+rep(1:5, 4),
                 Var=rep(c("x", "y"), 2,each=5), InterV=rep(c("Lower", "Upper"), each=10))
quanTsFrame$Interval <- interaction(quanTsFrame$Var, quanTsFrame$InterV)
quanTsFrame$Interval<- factor(quanTsFrame$Interval, levels=rev(levels(quanTsFrame$Interval)))

ggplot(resA, aes(x=t, y=Val, col=Value)) + geom_line(size=1) + theme_classic() + ylim(-5,4) +
  geom_line(data=resB_1[1:100,], aes(x=Time, y=Val), col="seagreen", size=1.1) +
  geom_line(data=resB_1[101:200,], aes(x=Time, y=Val), col="darkorange", size=1) +
  geom_line(data=resB_1[201:300,], aes(x=Time, y=Val), col="darkorange", size=1) +
  geom_line(data=quanTsFrame, aes(x=Time, y=x, col=Interval), size=1, linetype="dashed") +
  scale_color_manual(values=c("skyblue", "royalblue","darkmagenta","darkmagenta","red3", "red3"))+ labs(x="Time", title="0.95 probability intervals for k=5", caption="Figure 4")
KalmanF_2 <- function(yt, C=matrix(0), A=matrix(1),
                    B=matrix(0), u=matrix(0), omegaE=0, omegaV=0, TimeS=0){
  KalmanMat <- matrix(c(0,0,0,0), ncol=4, nrow=101)
  for(i in 1:TimeS){
    # Prediction update
    KalmanMat[i+1,3] <- A %*%  KalmanMat[i,1] + B%*%u #muBar_t 
    KalmanMat[i+1,4] <- A%*%KalmanMat[i,2]%*%t(A) + omegaV #sigmaBar_t
    # Measurement update
    K_t <- KalmanMat[i+1,4] %*%t(C) %*% solve(C%*%KalmanMat[i+1,4] %*%t(C) + omegaE)
    KalmanMat[i+1,1] <- KalmanMat[i+1,3] + K_t %*% (yt[i] - C%*%KalmanMat[i+1,3]) #mu_t
    KalmanMat[i+1,2] <- (diag(ncol(K_t))-K_t%*%C) %*% KalmanMat[i+1,4]  # sigma_t
  }
  KalmanFrame <- data.frame(KalmanMat[-1,], Time=1:TimeS)
  return(KalmanFrame)
}

BS <- function(A, Kalman, x){
  x_t <- data.frame(Val=0)
  for(k in (nrow(Kalman)-1):1){
    ht <- Kalman[k, 1] + Kalman[k, 2]%*%t(A)%*%solve(Kalman[k+1, 4]) %*% (x[k+1] - Kalman[k+1, 3])
    Ht <- Kalman[k,2]-Kalman[k,2]%*%t(A)%*%solve(Kalman[k+1,4])%*%A%*%Kalman[k,2]
    x_t[k,1] <- rnorm(1, mean=ht, sd=sqrt(Ht))
  }
  return(x_t)
}

bayesReg <- function(beta0, Xt, y, omega0, cov0){
  betaHat <- solve(t(Xt)%*%Xt)%*% t(Xt) %*% as.matrix(y)
  sigma0 <- diag(omega0/cov0,2)
  betaNew <- solve(t(Xt)%*%Xt + sigma0) %*% (t(Xt)%*%Xt%*%betaHat + sigma0%*%beta0)
  sigmaNew <- t(Xt)%*%Xt + sigma0
  Vari <- omega0 * solve(sigmaNew)
  coef <- rmvnorm(n=length(y), mean=betaNew, sigma=sqrt(diag(sqrt(diag(Vari)),length(betaNew))))
  return(coef)
}

# Gibbs sampling
set.seed(311015)
theta <- data.frame(matrix(ncol=2, nrow=5000)) #nrow sets number of sims
theta[1,2] <- 1  #start value for theta
Xt= matrix(data = c(rep(1,100),resA[1:100,1]), ncol=2)
for(j in 1:nrow(theta)){
  LGSS <- ssm1(C = matrix(theta[j,2]), x0 = 0, omegaE = sqrt(0.5), omegaV = sqrt(0.1),
              TimeS = 100)
  y = LGSS[101:200, 1]
  x = LGSS[1:100, 1]
  KalmanRes <- KalmanF_2(yt = y, C = matrix(theta[j,2]), omegaE = 0.5, omegaV = 0.1,
                         TimeS = 100)
  NewX <- (BS(A = matrix(1), Kalman = KalmanRes, x=x))
  Xt[1:99,2] <- as.numeric(NewX[,1])
  Xt[100,2] <- KalmanRes[100,1]
  coefs <- bayesReg(beta0 = matrix(data=c(0,0), ncol=1),Xt=Xt, y=y, omega0=0.5,cov0=100)
  theta[j+1,] <- colMeans(coefs) 
}
ggplot(theta, aes(y=X2, x=1:nrow(theta))) + geom_line() + theme_classic() + labs(caption="Figure 5", title="Trace plot", x="1:1000")
ggplot(theta, aes(x=X2)) + geom_histogram(binwidth = 0.1) + labs(caption="Figure 6", title=expression(Posterior~distribution~of~theta)) +  theme_classic()
mean(theta$X2)
# Assignment 2
load("C:\\Users\\Gustav\\Documents\\AdvML\\Lab4\\CAPM.Rda")

# a) - Estimating the variance components by MLE

CAPMv <- CAPM_data[,c(10,17,18)]
Zt <- matrix(CAPM_data$MARKET - CAPM_data$RKFREE)
Yt <- matrix(CAPMv$IBM - CAPMv$RKFREE)

buildLocalTrend <- function(x,data){
  V = exp(x[1])
  W = diag(exp(x[2:3]),2,2)
  return(dlm(
    m0 = c(0,1),
    C0 = diag(c(100,0.5),2),
    FF = matrix(c(1,1),1,2),     
    GG = diag(2),
    JFF = matrix(c(0,1),1,2),
    V = V,
    W = W,
    X=data))
}

initVal <- c(1,1,1) # Initial values for optim on the estimated parameters 
MLEs <- dlmMLE(Yt, parm = initVal, build = buildLocalTrend, data = Zt)
dlmWithMLEs <- buildLocalTrend(MLEs$par, data = Zt)

## Filtering ##
Filter_2b <- dlmFilter(y = Yt, dlmWithMLEs)

# Variances
u1 <- matrix((t(data.frame(Filter_2b$U.C)))[c(seq(1,242,2)),1])
u2 <- matrix((t(data.frame(Filter_2b$U.C)))[c(seq(2,242,2)),2])
u_1 <-list()
u_2 <- list()
for(i in 1:120){
  u_1[i]<-c(u1[i+1])
  u_2[i]<-c(u2[i+1])
}
Filter_y <- data.frame(alpha=c(Filter_2b$m)[2:121],beta=c(Filter_2b$m)[123:242])
Filter_y$alphaVar <- unlist(dlmSvd2var(u = u_1, d = matrix(Filter_2b$D.C[-1,1])))
Filter_y$betaVar <- unlist(dlmSvd2var(u = u_2, d = matrix(Filter_2b$D.C[-1,2])))
# Intervals
Filter_y$alphaLower <- Filter_y$alpha - 1.96*sqrt(Filter_y$alphaVar)
Filter_y$alphaUpper <- Filter_y$alpha + 1.96*sqrt(Filter_y$alphaVar)
Filter_y$betaLower <- Filter_y$beta - 1.96*sqrt(Filter_y$betaVar)
Filter_y$betaUpper <- Filter_y$beta + 1.96*sqrt(Filter_y$betaVar)

alphaFilter <- ggplot(Filter_y, aes(x=1:120, y=alpha)) + geom_line(col="seagreen", size=1) + theme_classic() + labs(title="Alpha filtered", caption="Figure 7") +
  geom_line(data=Filter_y, aes(x=1:120, y=alphaLower), col="darkorange", size=1) +
  geom_line(data=Filter_y, aes(x=1:120, y=alphaUpper), col="darkorange", size=1) 

betaFilter <- ggplot(Filter_y, aes(x=1:120, y=beta)) + geom_line( col="seagreen", size=1) + theme_classic()+ labs(title="Beta filtered", caption="Figure 9") +
  geom_line(data=Filter_y, aes(x=1:120, y=betaLower), col="darkorange", size=1) +
  geom_line(data=Filter_y, aes(x=1:120, y=betaUpper), col="darkorange", size=1)

#### Smoothed ####
Smooth_2b <- dlmSmooth(y = Yt, dlmWithMLEs)
Smooth_frame <- data.frame(alpha=Smooth_2b$s[-1,1], beta=Smooth_2b$s[-1,2])

covList = dlmSvd2var(Smooth_2b$U.S, Smooth_2b$D.S)
varMat = t(sapply(covList, FUN=function(x) sqrt(diag(x))))

Smooth_frame$alphaSd <- varMat[-1,1]
Smooth_frame$betaSd <- varMat[-1,2]
# Intervals
Smooth_frame$alphaLower <- Smooth_frame$alpha - 1.96*Smooth_frame$alphaSd
Smooth_frame$alphaUpper <- Smooth_frame$alpha + 1.96*Smooth_frame$alphaSd
Smooth_frame$betaLower <- Smooth_frame$beta - 1.96*Smooth_frame$betaSd
Smooth_frame$betaUpper <- Smooth_frame$beta + 1.96*Smooth_frame$betaSd

alphaSmooth <- ggplot(Smooth_frame, aes(x=1:120, y=alpha)) + geom_line(size=1, col="seagreen") + theme_classic() +labs(title="Alpha smoothed", caption="Figure 8") +
  geom_line(data=Smooth_frame, aes(x=1:120, y=alphaLower), size=1, col="darkorange") +
  geom_line(data=Smooth_frame, aes(x=1:120, y=alphaUpper), size=1, col="darkorange") 

betaSmooth <- ggplot(Smooth_frame, aes(x=1:120, y=beta)) + geom_line(size=1, col="seagreen")+theme_classic()+labs(title="Beta smoothed", caption="Figure 10") +
  geom_line(data=Smooth_frame, aes(x=1:120, y=betaLower),size=1 ,col="darkorange") +
  geom_line(data=Smooth_frame, aes(x=1:120, y=betaUpper),size=1, col="darkorange") 

grid.arrange(alphaFilter, alphaSmooth, ncol=2)
grid.arrange(betaFilter, betaSmooth, ncol=2)
# c)
TimeB <- data.frame(Beta=(dlmSmooth(y = Yt[1:48], dlmWithMLEs)$s)[,2])
TimeA <- data.frame(Beta=(dlmSmooth(y = Yt[1:120], dlmWithMLEs)$s)[,2])

Sensitivity<- data.frame(Beta=c(rnorm(10000, mean(TimeB$Beta), sd(TimeB$Beta)),
rnorm(10000, mean(TimeA$Beta), sd(TimeA$Beta))), Period=rep(c("Before 82", "After 82"),
                                                            each=10000))
Sensitivity$Period <- factor(Sensitivity$Period, levels=rev(levels(Sensitivity$Period)))

ggplot(Sensitivity, aes(x=Beta))+geom_density(aes(y=..density..,fill=Period),alpha=.7)+
  scale_fill_manual(values = c("seagreen", "darkorange"))+theme_classic() + labs(caption="Figure 11")
## NA
