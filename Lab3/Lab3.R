library(ggplot2)
# Lab 3
# Assignment 1
SqExpKernel <- function(x1, x2, hyperParam){
  K <- matrix(nrow=length(x1), ncol=length(x2))
  for (i in 1:length(x2)){
    K[, i] <- hyperParam[1]^2 * exp(-0.5 *( (x1-x2[i])/hyperParam[2]) ^2)
  }
  return(K)
}
PosteriorGP <- function(x, y, xStar, hyperParam, sigmaNoise){
  # Calculates f star bar
  K <- SqExpKernel(x, x, hyperParam) 
  L <- t(chol(K + sigmaNoise*diag(dim(K)[1]) )) 
  alpha <- solve(t(L),solve(L,y))
  # Posterior mean
  fStar <- (SqExpKernel(xStar, x, hyperParam)) %*% alpha 
  # Posterior variance
  v_f <- solve(L, t(SqExpKernel(xStar, x, hyperParam)))
  cov_fStar <- SqExpKernel(xStar, xStar, hyperParam) - (t(v_f) %*% v_f)  
  # Store all values in a list
  val_list <- list(fStar=fStar, xStar=xStar, cov_fStar=cov_fStar, xStar=xStar
  )
  return(val_list)
}
assign1B <- PosteriorGP(x=0.4, y=0.719, xStar=seq(-1,1, 0.01), hyperParam=c(1, 0.3),
            sigmaNoise=0.1^2)
Upper1B <- assign1B$fStar + 1.96 * sqrt(diag(assign1B$cov_fStar))
Lower1B <- assign1B$fStar - 1.96 * sqrt(diag(assign1B$cov_fStar))
plotD <- data.frame(fStar=assign1B$fStar, xStar=assign1B$xStar, Lower=Lower1B, Upper=Upper1B)
xy <- data.frame(x=0.4, y=0.719)
ggplot(plotD, aes(y=fStar, x=xStar)) + geom_point(col="darkorange") + ylim(-2,2.2) + 
  geom_line(data=plotD, aes(xStar, Lower), col="seagreen", size=1.1) +
  geom_line(data=plotD, aes(xStar, Upper), col="seagreen", size=1.1) +geom_point(data=xy, aes(x,y), size=3) + theme_classic()
assign1C <- PosteriorGP(x=c(0.4, -0.6), y=c(0.719, -0.044), xStar=seq(-1,1, 0.01), hyperParam=c(1, 0.3),sigmaNoise=0.1^2)
Upper1C <- assign1C$fStar + 1.96 * sqrt(diag(assign1C$cov_fStar))
Lower1C <- assign1C$fStar - 1.96 * sqrt(diag(assign1C$cov_fStar))

plotD <- data.frame(fStar=assign1C$fStar, xStar=assign1C$xStar, Lower=Lower1C, Upper=Upper1C)
xy <- data.frame(x=c(0.4, -0.6), y=c(0.719, -0.044))
ggplot(plotD, aes(y=fStar, x=xStar)) + geom_point(col="darkorange") + ylim(-2,2.2) + 
  geom_line(data=plotD, aes(xStar, Lower), col="seagreen", size=1.1) +
  geom_line(data=plotD, aes(xStar, Upper), col="seagreen", size=1.1) +geom_point(data=xy, aes(x,y), size=3) +theme_classic()
assign1D <- PosteriorGP(x=c(0.8, 0.4, -0.2, -0.6, -1), y=c(-0.664, 0.719, -0.94, -0.044, 0.768), xStar=seq(-1,1, 0.01),hyperParam=c(1, 0.3),sigmaNoise=0.1^2)
Upper1D <- assign1D$fStar + 1.96 * sqrt(diag(assign1D$cov_fStar))
Lower1D <- assign1D$fStar - 1.96 * sqrt(diag(assign1D$cov_fStar))

plotD <- data.frame(fStar=assign1D$fStar, xStar=assign1D$xStar, Lower=Lower1D, Upper=Upper1D)
xy <- data.frame(x=c(0.8, 0.4, -0.2, -0.6, -1), y=c(-0.664, 0.719, -0.94, -0.044, 0.768))
ggplot(plotD, aes(y=fStar, x=xStar)) + geom_point(col="darkorange") + ylim(-2,2.2) + 
  geom_line(data=plotD, aes(xStar, Lower), col="seagreen", size=1.1) +
  geom_line(data=plotD, aes(xStar, Upper), col="seagreen", size=1.1) +geom_point(data=xy, aes(x,y), size=3) +theme_classic()
assign1E <- PosteriorGP(x=c(0.8, 0.4, -0.2, -0.6, -1), y=c(-0.664, 0.719, -0.94, -0.044, 0.768), xStar=seq(-1,1, 0.01),hyperParam=c(1, 1),sigmaNoise=0.1^2)
Upper1E <- assign1E$fStar + 1.96 * sqrt(diag(assign1E$cov_fStar))
Lower1E <- assign1E$fStar - 1.96 * sqrt(diag(assign1E$cov_fStar))

plotD <- data.frame(fStar=assign1E$fStar, xStar=assign1E$xStar, Lower=Lower1E, Upper=Upper1E)
xy <- data.frame(x=c(0.8, 0.4, -0.2, -0.6, -1), y=c(-0.664, 0.719, -0.94, -0.044, 0.768))
ggplot(plotD, aes(y=fStar, x=xStar)) + geom_point(col="darkorange") + ylim(-2,2.2) + 
  geom_line(data=plotD, aes(xStar, Lower), col="seagreen", size=1.1) +
  geom_line(data=plotD, aes(xStar, Upper), col="seagreen", size=1.1) +geom_point(data=xy, aes(x,y), size=3) +theme_classic()
# Assignment 2
library(kernlab)
temps <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv",header=TRUE,   sep=";")
temps$time <- 1:2190
temps$day <- rep(1:365, 6)
temps$index <- rep(1:5, 438)
thinTemps <- subset(temps, temps$index == 1)[,1:4]
sqExpK <- function(sigmaf = 1, ell = 1) {
  rval <- function(x, xStar = NULL) {
    return(sigmaf * exp (-1/2 * ((x-xStar)/ell)^2 ))
  }
  class(rval) <- "kernel"
  return(rval)
} 
sqExpKFunc = sqExpK(sigmaf = 1, ell = 2) # sqExpKFunc is a kernel FUNCTION
sqExpKFunc(1,2) # Evaluating the kernel in x=1, x'=2
# Computing the whole covariance matrix K from the kernel.
x <- as.matrix(c(1,3,4))
xStar <- as.matrix(c(2,3,4))
K <- kernelMatrix(kernel = sqExpKFunc, x = x, y = xStar) # So this is K(X,Xstar)
K
quadLM <- lm(temp ~ time + time^2, data=temps)
sigma_n <- var(quadLM$residuals)

# Correct parametrizesed?
sqExpKFunc = sqExpK(sigmaf = 20^2, ell = 0.2)
reg <- gausspr(thinTemps$time, thinTemps$temp, kernel = sqExpKFunc, var=sigma_n)
postMean <- data.frame(pred=predict(reg), y=thinTemps$temp, x=thinTemps$time)
plotB <- ggplot(postMean, aes(x=x, y=y)) + geom_point(col="grey75") + geom_line(data=postMean, aes(x=x, y=pred), col="darkorange", size=1.1) + theme_classic()
plotB
covar_f <- PosteriorGP(x=thinTemps$time, y = thinTemps$temp, xStar = thinTemps$time, 
                       hyperParam=c(20^2, 0.2),sigmaNoise =sigma_n )
probBands <- data.frame(Upper=postMean[,1] + 1.96 * sqrt(diag(covar_f$cov_fStar)) , Lower = postMean[,1] - 1.96 * sqrt(diag(covar_f$cov_fStar)), x=thinTemps$time)
plotB + geom_line(data=probBands, aes(x=x, y=Upper), size=1.1, col="seagreen") + geom_line(data=probBands, aes(x=x, y=Lower), size=1.1, col="seagreen")
quadLM2 <- lm(temp ~ day + day^2, data=temps)
sigma_n2 <- var(quadLM2$residuals)
# Small difference, 67.3641 vs 64.10528
sqExpKFuncD = sqExpK(sigmaf = 20^2, ell = 1.2)
regD <- gausspr(thinTemps$day, thinTemps$temp, kernel = sqExpKFuncD, var=sigma_n2)
postMeanD <- data.frame(pred=predict(regD), y=thinTemps$temp, x=thinTemps$time)
postMT <- data.frame(rbind(postMean, postMeanD), model=c(rep("m1",438), rep("m2", 438)))
ggplot(postMT, aes(y=y, x=x)) + geom_point(col="grey75") + geom_line(aes(x=x, y=pred, col=model), size=1.1) +theme_classic() + scale_color_manual(values = c("darkorange", "seagreen"))
PeriodicK <- function(sigmaf = 1, ell_one = 1, ell_two=1, d=1){
  rval <- function(x, xStar = NULL) {
    return(sigmaf * exp(-(2*sin(pi*abs(x-xStar)/d)^2)/ ell_one^2) *
             exp(-0.5 * (abs(x-xStar)^2 / ell_two^2)))
  }
  class(rval) <- "kernel"
  return(rval)
}
PeriodicKFunc = PeriodicK(sigmaf = 20^2, ell_one = 1, ell_two=10, d=365/sd(thinTemps$time))
#PeriodicKFunc is a kernel FUNCTION
regE <- gausspr(thinTemps$time, thinTemps$temp, kernel = PeriodicKFunc, var=sigma_n)
postMeanE <- data.frame(pred=predict(regE), y=thinTemps$temp, x=thinTemps$time)
plotE <- ggplot(postMeanE, aes(x=x, y=y)) + geom_point(col="grey75") + geom_line(data=postMeanE, aes(x=x, y=pred), col="darkorange", size=1.1) + theme_classic()
postMT <- rbind(postMT, data.frame(postMeanE, model="m3"))
plotAll <- ggplot(postMT, aes(x=x, y=pred)) + geom_line(aes(col=model, linetype=model),size=1.1, alpha=0.8) + theme_classic() + scale_color_manual(values = c("darkorange", "seagreen", "skyblue")) + scale_linetype_manual(values = c("dashed", "longdash", "solid")) + theme(legend.justification = c(1,0), legend.position=c(1,0)) + scale_y_continuous(limits = c(-17,20))
library(gridExtra)
grid.arrange(plotE, plotAll, ncol=2)
# Assignment 3
library(AtmRay)
data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/banknoteFraud.csv",
                  header=FALSE, sep=",")
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])
set.seed(111); SelectTraining <- sample(1:dim(data)[1], size = 1000, replace = FALSE)
Train <- data[SelectTraining, ]
Test <- data[-SelectTraining, ]
# a)
cfA <- gausspr(fraud ~ varWave + skewWave, data=Train)
gridX <- seq(min(Train$varWave), max(Train$varWave), length=100)
gridY <- seq(min(Train$skewWave), max(Train$skewWave), length=100)
gridP <- meshgrid(gridX, gridY)
gridP <- cbind(c(gridP$x), c(gridP$y))
gridP <- data.frame(gridP)
names(gridP) <- names(Train)[1:2]
probPredsA <- predict(cfA, gridP, type="probabilities")
contour(x=gridX, y= gridY,z=matrix(probPredsA[,2],100),20)
points(Train[Train[,5]== 1,1],Train[Train[,5]== 1,2],col="blue")
points(Train[Train[,5]== 0,1],Train[Train[,5]== 0,2],col="red")
# predict on the training set
confMatA <- table(predict(cfA,Train[,1:2]), Train[,5]) # confusion matrix
confMatA
cat("Accuracy:", sum(diag(confMatA)) / sum(confMatA))
confMatB <- table(predict(cfA,Test[,1:2]), Test[,5]) # confusion matrix
confMatB
cat("Accuracy:", sum(diag(confMatB)) / sum(confMatB))
cfC <- gausspr(fraud ~ varWave + skewWave + kurtWave + entropyWave, data=Train)
confMatC <- table(predict(cfC,Test[,1:4]), Test[,5]) # confusion matrix
cat("Accuracy:", sum(diag(confMatC)) / sum(confMatC))
## NA
