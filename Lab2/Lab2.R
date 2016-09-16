
library(HMM)
library(ggplot2)
library(scales)
# 1
State <- c(1:10)
Symbol <- c(1:10)
startP <- rep(0.1, 10)

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

HMM1 <- initHMM(States = State, Symbols = Symbol, startProbs = startP,
        transProbs = transP, emissionProbs = emissionP)

# 2
n = 500
Sim1 <- simHMM(HMM1, n)

# 3
# Filtered prob dist
options(scipen=1)

alpha <- prop.table(exp(forward(HMM1, Sim1$observation)))
filtering <- matrix(0, ncol=n, nrow=10)
for(k in 1:n){
  filtering[,k] <- alpha[,k] / colSums(alpha)[k]
}
filter <- data.frame(state=0, obs=1:n)
for(h in 1:n){
  filter[h,1] <- which.max(filtering[,h])
}
ggplot(filter, aes(x=obs, y=state)) + geom_path() + labs(x="Path", y="State") +
  scale_y_continuous(breaks=c(1:10)) + theme_minimal()

# Smoothing
smoothing <- posterior(HMM1, Sim1[[2]])
smooth <- data.frame(state=0, obs=1:n)
for(h in 1:n){
  smooth[h,1] <- which.max(smoothing[,h])
}
ggplot(smooth, aes(x=obs, y=state)) + geom_path() + labs(x="Path", y="State") +
  scale_y_continuous(breaks=c(1:10)) + theme_minimal()


# Most probable path
path <- data.frame(state=as.numeric(viterbi(HMM1, Sim1[[2]])), obs=1:n)
ggplot(path, aes(x=obs, y=state)) + geom_path() + labs(x="Path", y="State") +
  scale_y_continuous(breaks=c(1:10)) + theme_minimal()

# 4
sum(Sim1[[1]] == filter$state) / n
sum(Sim1[[1]] == smooth$state) / n
sum(Sim1[[1]] == path$state) / n

prop.table(table(filter$state))
prop.table(table(smooth$state))
prop.table(table(path$state))


# 5
## Repeat the code above but with different sample sizes
accuracy <- data.frame(sampleS=seq(50,350,5), accF=0, accS=0, accP=0)

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
  scale_y_continuous(labels=percent)

colMeans(accuracy[,2:4])

# 6
library(entropy)

entrF <- data.frame(sampleS=seq(50,350,5), entropy=0)
for(t in 1:nrow(entrF)){
  Sim3 <- simHMM(HMM1, accuracy[t,1])
  # Filtering
  alpha <- exp(forward(HMM1, Sim3[[2]]))
  filtering <- matrix(0, ncol=accuracy[t,1], nrow=10)
  for(k in 1:accuracy[t,1]){
    filtering[,k] <- alpha[,k] / colSums(alpha)[k]
  }  
  filter <- data.frame(state=0, obs=1:accuracy[t,1])
  for(l in 1:accuracy[t,1]){
    filter[l,1] <- which.max(filtering[,l])
  }
  entrF[t,2] <- entropy.empirical(filter$state)
  
}
ggplot(entrF, aes(x=sampleS, y=entropy)) + geom_line() + theme_minimal()


# 7
TP <- matrix(c(0,0,0,0,0.5,0.5,0,0,0,0),ncol=1)
Alph <- matrix(filtering[,100])

TP * Alph / sum(TP * Alph )

(TP) %*% t(Alph)
Sim1$observation

