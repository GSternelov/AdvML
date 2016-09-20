## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
library(HMM)
library(ggplot2)
library(scales)
library(gridExtra)

## ---- echo=FALSE---------------------------------------------------------
# 1
seq_eP <- rep(1:10,10)
transP <- matrix(data=0, nrow=10, ncol=10)
for(i in 1:10){
  transP[i, c(seq_eP[i], seq_eP[i+1])] <- 0.5
}
emissionP <- matrix(data=0, nrow=10, ncol=10)
  for(i in 1:10){
  j <- i+10
  emissionP[c(seq_eP[j-2],seq_eP[j-1],seq_eP[j], seq_eP[j+1],seq_eP[j+2]),i] <- 0.2
  }

HMM1 <- initHMM(States = c(1:10), Symbols = c(1:10), startProbs =  rep(0.1, 10),
        transProbs = transP, emissionProbs = emissionP)
# 2
set.seed(311015)
Sim1 <- simHMM(HMM1, 100)

## ---- echo=FALSE, fig.width=8, fig.height=3.5----------------------------
# 3
# Filtered prob dist
options(scipen=99)
alpha <- exp(forward(HMM1, Sim1$observation))
filtering1 <- matrix(0, ncol=100, nrow=10)
for(k in 1:100){
  filtering1[,k] <- alpha[,k] / colSums(alpha)[k]
}
filter <- data.frame(state=0, obs=1:100)
for(h in 1:100){
  filter[h,1] <- which.max(filtering1[,h])
}
filter_p <- ggplot(filter, aes(x=obs, y=state)) + geom_path() + labs(x="Path", y="State",title="Filtered distribution") + scale_y_continuous(breaks=c(1:10)) + theme_minimal()
# Smoothing
smoothing <- posterior(HMM1, Sim1[[2]])
smooth <- data.frame(state=0, obs=1:100)
for(h in 1:100){
  smooth[h,1] <- which.max(smoothing[,h])
}
smooth_p <- ggplot(smooth, aes(x=obs, y=state)) + geom_path() + labs(x="Path", y="", title= "Smoothed distribution") + scale_y_continuous(breaks=c(1:10)) + theme_minimal()
# Most probable path
path <- data.frame(state=as.numeric(viterbi(HMM1, Sim1[[2]])), obs=1:100)
path_p <- ggplot(path, aes(x=obs, y=state)) + geom_path() + labs(x="Path", y="", title="Most probable path") + scale_y_continuous(breaks=c(1:10)) + theme_minimal()
plots <- list(filter_p, smooth_p, path_p)
plots1 <- arrangeGrob(grobs = plots, nrow=1)
plot(plots1)

## ---- echo=FALSE, fig.width=7, fig.height=3.5----------------------------
filter$distribution <- "Filtered"
smooth$distribution <- "Smoothed"
path$distribution <- "Path"
allDist <- rbind(filter, smooth, path)
ggplot(allDist, aes(x=obs, y=state, col=distribution)) + geom_path() + labs(x="Path", y="State", title="Comparsion of the distributions") + scale_y_continuous(breaks=c(1:10)) + theme_minimal()

## ---- echo=FALSE---------------------------------------------------------
# 4
data.frame(AccuracyFiltered = sum(Sim1[[1]] == filter$state) / 100, AccuracySmoothed=sum(Sim1[[1]] == smooth$state) / 100, AccuracyPath= sum(Sim1[[1]] == path$state) / 100)

## ---- echo=FALSE---------------------------------------------------------

data.frame(Filtered=prop.table(table(filter$state))[1:10],Smoothed=prop.table(table(smooth$state))[1:10], Path=prop.table(table(path$state))[1:10])


## ---- echo=FALSE, fig.width=6, fig.height=3.5----------------------------
# 5
## Repeat the code above but with different sample sizes
accuracy <- data.frame(sampleS=seq(50,350,10), accF=0, accS=0, accP=0)
set.seed(0814)
for(t in 1:nrow(accuracy)){
  Sim2 <- simHMM(HMM1, accuracy[t,1])
  # Filtering
  alpha <- exp(forward(HMM1, Sim2[[2]]))
  filtering <- matrix(0, ncol=accuracy[t,1], nrow=10)
  for(k in 1:accuracy[t,1]){
    filtering[,k] <- alpha[,k] / colSums(alpha)[k]
  }
  filter <- data.frame(state=0, obs=1:accuracy[t,1])
  for(l in 1:accuracy[t,1]){
    filter[l,1] <- which.max(filtering[,l])
  }
  # Smoothing
  smoothing <- posterior(HMM1, Sim2[[2]])
  smooth <- data.frame(state=0, obs=1:accuracy[t,1])
  for(h in 1:accuracy[t,1]){
    smooth[h,1] <- which.max(smoothing[,h])
  }
  # Path
  path <- data.frame(state=as.numeric(viterbi(HMM1, Sim2[[2]])), obs=1:accuracy[t,1])
  # Accuracy
  accuracy[t,2] <- sum(Sim2[[1]] == filter$state) / accuracy[t,1]
  accuracy[t,3] <- sum(Sim2[[1]] == smooth$state) / accuracy[t,1]
  accuracy[t,4] <- sum(Sim2[[1]] == path$state) / accuracy[t,1]
}
accuracy2 <- tidyr::gather(accuracy, "method", "accuracy", 2:4)
ggplot(accuracy2, aes(x=sampleS, y=accuracy, col=method)) + geom_line() + theme_minimal() +
  scale_y_continuous(labels=percent) + scale_x_continuous(breaks=seq(50,350,50))

## ---- echo=FALSE, fig.width=6, fig.height=3.5, message=FALSE, warning=FALSE----
# 6
library(entropy)
entrF <- data.frame(sampleS=seq(20,350,10), entropy=0)
set.seed(311015)
for(t in 1:nrow(entrF)){
  Sim3 <- simHMM(HMM1, entrF[t,1])
  # Filtering
  alpha <- exp(forward(HMM1, Sim3[[2]]))
  filtering <- matrix(0, ncol=entrF[t,1], nrow=10)
  for(k in 1:entrF[t,1]){
    filtering[,k] <- alpha[,k] / colSums(alpha)[k]
  }  
  filter <- data.frame(state=0, obs=1:entrF[t,1])
  for(l in 1:entrF[t,1]){
    filter[l,1] <- which.max(filtering[,l])
  }
  entrF[t,2] <- entropy.empirical(table(filter$state))
}
ggplot(entrF, aes(x=sampleS, y=entropy)) + geom_line() + theme_minimal()+ scale_x_continuous(breaks=seq(0,350,50))

## ---- echo=FALSE---------------------------------------------------------
TP <- matrix(c(0,0,0,0,0,0,0.5,0.5,0,0),ncol=1)
Alph <- matrix(filtering1[,100])

TP * Alph / sum(TP * Alph )

## ----code=readLines(knitr::purl("C:\\Users\\Gustav\\Documents\\AdvML\\Lab2\\Lab2.Rmd",documentation = 1)), eval = FALSE, tidy=TRUE----
## 
## 

