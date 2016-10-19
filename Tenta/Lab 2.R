library(HMM)
library(ggplot2)
library(scales)
library(gridExtra)
library(entropy)


# 1. Initialize the model
# Transition matrix is unchanged
# Emission matrix is modified, see below

seq_eP <- rep(1:10,10)
transP <- matrix(data=0, nrow=10, ncol=10)
for(i in 1:10){
  transP[i, c(seq_eP[i], seq_eP[i+1])] <- 0.5
}
# Emission matrix only give probabilities for states +-1
emissionP <- matrix(data=0, nrow=10, ncol=10)
for(i in 1:10){
  j <- i+10
  emissionP[c(seq_eP[j-1],seq_eP[j], seq_eP[j+1]),i] <- 0.2
}

HMM1 <- initHMM(States = c(1:10), Symbols = c(1:10), startProbs =  rep(0.1, 10),
                transProbs = transP, emissionProbs = emissionP)

# 100 observations are simulated from the model
set.seed(311015)
Sim1 <- simHMM(HMM1, 100)

# A plot over the simulated states and the observations
Sims1 <- data.frame(Time=rep(1:100, 2), Data=c(Sim1$states, Sim1$observation),
                    Var=rep(c("State", "Observation"), each=100))
ggplot(Sims1, aes(x=Time, y=Data, col=Var)) + geom_line(size=1) + theme_classic()
# State sequence never moves backwards but observation sometimes does


# The filtered distribution is computed
options(scipen=99)
alpha <- exp(forward(HMM1, Sim1$observation))
alpha <- prop.table(alpha)
filtFunc <- function(Data){
  Data <- Data / sum(Data)
}
filtering <- apply(alpha,2,filtFunc)
probState <- function(Data){
  Data <- which.max(Data)
}
Filter <- data.frame(state=apply(filtering, 2, probState), Time=1:100)

filter_p <- ggplot(Filter, aes(x=Time, y=state)) + geom_path() +
  labs(x="Path", y="State",title="Filtered distribution") +
  scale_y_continuous(breaks=c(1:10)) + theme_classic()

# Smoothed distribution
smoothing <- posterior(HMM1, Sim1[[2]])
smooth <- data.frame(state=apply(smoothing, 2, probState), Time=1:100) 
smooth_p <- ggplot(smooth, aes(x=Time, y=state)) + geom_path() +
  labs(x="Path", y="", title= "Smoothed distribution") +
  scale_y_continuous(breaks=c(1:10)) + theme_classic()

# The most probable path
path <- data.frame(state=as.numeric(viterbi(HMM1, Sim1[[2]])), Time=1:100)
path_p <- ggplot(path, aes(x=Time, y=state)) + geom_path() +
  labs(x="Path", y="", title="Most probable path") +
  scale_y_continuous(breaks=c(1:10)) + theme_classic()

# All distributions next to each other
plots <- list(filter_p, smooth_p, path_p)
plots1 <- arrangeGrob(grobs = plots, nrow=1)
plot(plots1)

# All distributions in the same plot
Filter$distribution <- "Filtered"
smooth$distribution <- "Smoothed"
path$distribution <- "Path"
allDist <- rbind(Filter, smooth, path)
ggplot(allDist, aes(x=Time, y=state, col=distribution)) + geom_path(size=1) +
  labs(x="Path", y="State", title="Comparsion of the distributions") +
  scale_y_continuous(breaks=c(1:10)) + theme_classic()

# The accuracy for the respective distribution
data.frame(AccuracyFiltered = sum(Sim1[[1]] == Filter$state) / 100,
           AccuracySmoothed=sum(Sim1[[1]] == smooth$state) / 100,
           AccuracyPath= sum(Sim1[[1]] == path$state) / 100)

# The frequency for each state for all distributions
# Not sure if this says something of importance
data.frame(Filtered=prop.table(table(Filter$state))[1:10],
           Smoothed=prop.table(table(smooth$state))[1:10],
           Path=prop.table(table(path$state))[1:10])


# The accuracy with different samples
# Samples 100 observations 1000 times and calculate the accuarcy for each sample
accuracy <- data.frame(sampleS=1:400, accF=0, accS=0, accP=0)
set.seed(0814)
for(t in 1:nrow(accuracy)){
  Sim2 <- simHMM(HMM1, 100)
  #Filtering
  alpha <- exp(forward(HMM1, Sim2$observation))
  alpha <- prop.table(alpha)
  filtering <- apply(alpha,2,filtFunc)
  Filter <- data.frame(state=apply(filtering, 2, probState), Time=1:100)
  # Smoothing
  smoothing <- posterior(HMM1, Sim2[[2]])
  smooth <- data.frame(state=apply(smoothing, 2, probState), Time=1:100) 
  # MPP
  path <- data.frame(state=as.numeric(viterbi(HMM1, Sim2[[2]])), Time=1:100)
  accuracy[t,2] <- sum(Sim2[[1]] == Filter$state) / 100
  accuracy[t,3] <- sum(Sim2[[1]] == smooth$state) /100
  accuracy[t,4] <- sum(Sim2[[1]] == path$state) / 100
}
# tidyr to transform data to ggplot format
accuracy2 <- tidyr::gather(accuracy, "method", "accuracy", 2:4)
ggplot(accuracy2, aes(x=sampleS, y=accuracy, col=method)) + geom_line( size=1)+
  theme_classic()+scale_y_continuous(labels=percent) +
  scale_x_continuous(breaks=seq(50,350,50))
# Smoothed the most accurate. Very even between filtered and MPP

# The entropy for each time step for the filtered distribution
set.seed(1897)
Sim3 <- simHMM(HMM1, 100)
alpha <- prop.table(exp(forward(HMM1, Sim3$observation)))
filtering <- apply(alpha,2,filtFunc)
Estimated_Entropy <- c()
for (i in 1:100) {
  Estimated_Entropy[i] <- entropy.empirical(filtering[, i])
}
gg_entropy <- data.frame(Entropy = Estimated_Entropy, Time = 1:100)
ggplot(gg_entropy, aes(x = Time, y = Entropy)) + geom_line()

# Higher value indicates higher uncertainty over where the robot is.
# More sure for lower values and knows exactly where it is if entropy = 0
# Can see that the certainty for the robot's position not increases with time.
# That we knew the position in the previous time step does not necessarily helps us
# in the next time step. 


# For the simulation in the previous assignment are the probabilities for state 101
# calculated. This is done by using the probabilities from the filtering for the 
# 100th simulation and the whole transition matrix

filtering[, 100] %*% transP/sum(filtering[, 100] %*% transP) 

# Highest probability for Z101 = 6. That is the individually most probable state
# for time step 101. 


