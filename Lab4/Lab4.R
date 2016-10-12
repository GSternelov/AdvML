library(ggplot2)
library(scales)
library(dlm)

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
load(file ="Lab4/simulated_data.RDA")
ggplot(simulated_data, aes(x=time, y=y)) + geom_line(col="royalblue") + theme_classic() + 
  geom_line(data=simulated_data, aes(x=time, y=state), col="red3") + ylim(-4,4)


# b) - Filtering: Sequential state inference

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

resB <- KalmanF(yt = resA[101:200,1], C = matrix(0.93), omegaE = 0.5, omegaV = 0.1, TimeS = 100)
resB_1 <- data.frame(Val=c(resB[,1],resB[,1]+1.96*sqrt(resB[,2]),resB[,1]-1.96*sqrt(resB[,2])),
                     Time=rep(1:100, 3), Value=rep(c("Kalman mean", "Upper", "Lower"), each=100))
resB_1$Value <- factor(resB_1$Value, levels=rev(levels(resB_1$Value)))
ggplot(resB_1, aes(x=Time, y=Val, col=Value))+geom_line(size=1) + theme_classic() +ylim(-5,4)+
  scale_colour_manual(values=c("darkorange","darkorange","seagreen"))

plotA + geom_line(data=resB_1, aes(x=Time, y=Val, col=Value), size=1) +ylim(-5,4)+
  scale_colour_manual(values=c("seagreen", "darkorange","skyblue","royalblue","darkorange"))+
  labs(x="Time", title="Kalman filter estimates with \n0.95 probability intervals")

model <- dlm(m0=0, V=0.5, W=0.1, C0=10, FF=0.93, GG=1)
filteringB <- dlmFilter(y=resA[101:200,1], model )

filterB <- data.frame(y=(filteringB$m)[-1])
filterVar <- unlist(dlmSvd2var(u = filteringB$U.C, d = filteringB$D.C))
filterB[101:300,1] <- c(filterB$y + 1.96*sqrt(filterVar[-1]), 
                        filterB$y - 1.96*sqrt(filterVar[-1]))
filterB$Time <- rep(1:100, 3)
filterB$Value <- rep(c("dlm Kalman mean", "dlm Upper", "dlm Lower"), each=100)

ggplot(resB_1, aes(x=Time, y=Val, col=Value, linetype=Value))+geom_line(size=1)+theme_classic()+
  ylim(-5,4)+ geom_line(data=filterB, aes(x=Time, y=y, col=Value, linetype=Value),size=1) +
  scale_colour_manual(values=c("yellow","red3","red3","seagreen","darkorange","darkorange"))+
  scale_linetype_manual(values=c("dashed","dashed","dashed","solid","solid","solid"))

  
plotA +
  ggplot() + geom_line(data=resB, aes(x=1:100, y=V1), col="seagreen", size=1.1) +
  geom_line(data=resB, aes(x=1:100, y=Upper), col="darkorange", size=1) +
  geom_line(data=resB, aes(x=1:100, y=Lower), col="darkorange", size=1) + ylim(-4.2,4)  +
  geom_line(data=filterB, aes(x=1:100, y=y), col="pink", size=1) +

ggplot(filterB, aes(x=1:100, y=y)) + geom_line(col="seagreen") + theme_classic() +
  geom_line(data=filterB, aes(x=1:100, y=Upper), col="darkorange") +
  geom_line(data=filterB, aes(x=1:100, y=Lower), col="darkorange") + ylim(-5,4)
  

# c) - Prediction of state and data by simulation
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

ggplot(resA[1:100,], aes(y=Val, x=t)) + geom_line() + theme_classic() + ylim(-5,4) +
  geom_line(data=resB_1[1:100,], aes(x=Time, y=Val), col="seagreen", size=1.1) +
  geom_line(data=resB_1[101:200,], aes(x=Time, y=Val), col="darkorange", size=1) +
  geom_line(data=resB_1[201:300,], aes(x=Time, y=Val), col="darkorange", size=1) +
  geom_line(data=quanTsFrame, aes(x=Time, y=x, col=Interval), size=1, linetype="dashed") +
  scale_color_manual(values=c("royalblue","red3","royalblue", "red3"))+
  labs(x="Time", title="0.95 probability intervals for k=5")

dlmFc <- dlmForecast(mod = model, nAhead = 5, sampleNew = 10000)
newStates <- t(apply(data.frame(dlmFc$newStates), 1, quantile, probs=c(0.025, 0.975)))
newObs <- t(apply(data.frame(dlmFc$newObs), 1, quantile, probs=c(0.025, 0.975)))

quanTsFrame2 <- data.frame(x=c(newStates[,1],newStates[,2],c(newObs[,1],newObs[,2])),
                           Time= 100+rep(1:5, 4),Var=rep(c("x", "y"),each=10),
                           InterV=rep(c("Lower", "Upper"),2, each=5))
quanTsFrame2$Interval <- interaction(quanTsFrame2$Var, quanTsFrame2$InterV)
quanTsFrame2$Interval<- factor(quanTsFrame2$Interval, levels=rev(levels(quanTsFrame2$Interval)))

ggplot(resA[1:100,], aes(y=Val, x=t)) + geom_line() + theme_classic() + ylim(-6.5,6.5) +
  geom_line(data=resB_1[1:100,], aes(x=Time, y=Val), col="seagreen", size=1.1) +
  geom_line(data=resB_1[101:200,], aes(x=Time, y=Val), col="darkorange", size=1) +
  geom_line(data=resB_1[201:300,], aes(x=Time, y=Val), col="darkorange", size=1) +
  geom_line(data=quanTsFrame2, aes(x=Time, y=x, col=Interval), size=1, linetype="dashed") +
  scale_color_manual(values=c("royalblue","red3","royalblue", "red3"))+
  labs(x="Time", title="0.95 probability intervals for k=5")

# d) - State and model parameter inference



# Assignment 2
load("C:\\Users\\Gustav\\Documents\\AdvML\\Lab4\\CAPM.Rda")

# a) - Estimating the variance components by MLE

CAPMv <- CAPM_data[,c(10,17,18)]
Zt <- matrix(CAPM_data$MARKET - CAPM_data$RKFREE)

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
dlmLocalTrend <- buildLocalTrend(initVal, data = Zt) # Just to see that the build function works
MLEs <- dlmMLE(as.matrix(CAPMv), parm = initVal, build = buildLocalTrend, data = Zt)
dlmWithMLEs <- buildLocalTrend(MLEs$par, data = Zt)

Filter <- dlmFilter(y = as.matrix(CAPMv$IBM), dlmWithMLEs)
Filter_y <- data.frame(Filter$m)
plot(Filter_y[-1,1], type="l", ylim = c(0,1))
lines(Filter_y[-1,2])


Smooth <- dlmSmooth(y = as.matrix(CAPMv$IBM), dlmWithMLEs)
SmoothS <- Smooth$s
plot(Smooth$s[-1,1], type="l")
plot(Smooth$s[-1,2], type="l")


# c)
TimeB <- data.frame(Beta=(dlmSmooth(y = as.matrix(CAPMv[2:49,1]), dlmWithMLEs)$s)[,2])
TimeA <- data.frame(Beta=(dlmSmooth(y = as.matrix(CAPMv[51:121,1]), dlmWithMLEs)$s)[,2])

Sensitivity<- data.frame(Beta=c(rnorm(10000, mean(TimeB$Beta), sd(TimeB$Beta)),
rnorm(10000, mean(TimeA$Beta), sd(TimeA$Beta))), Period=rep(c("Before", "After"),
                                                            each=10000))

par(mfrow=c(2,1))
hist(Sensitivity[1:10000,1], xlim = c(-.1,.1), col="blue", main="Before 82")
hist(Sensitivity[10001:20000,1], xlim = c(-.1,.1), col="red", main="After 82")


ggplot(Sensitivity, aes(x=Beta, fill=Period)) + theme_classic() +
  geom_histogram(binwidth=0.01) + xlim(-0.1,0.1)
  scale_y_continuous(labels=percent)

ggplot(Sensitivity[1:10000,], aes(x=Beta)) + theme_classic() +
  geom_histogram(binwidth=0.01) 
ggplot(Sensitivity[10001:20000,], aes(x=Beta)) + theme_classic() +
  geom_histogram(binwidth=0.01) 





