library(bnlearn)
# Assignment 1
data("learning.test")
data("alarm")
data("asia")

set.seed(1897)

init_graph <- random.graph(LETTERS[1:6])
res1 <- hc(learning.test, start = init_graph, restart = 100)
all.equal(cpdag(init_graph), cpdag(res1))
plot(init_graph, main="DAG 1")
plot(res1, main="DAG 2, with 100 restarts")
score(init_graph, learning.test)
score(res1, learning.test)
# Assignment 2
library(ggplot2)
res_asia <- hc(asia)
res1_asia <- hc(asia, start = res_asia, restart = 100)
p1 <- hc(asia, start=res_asia, restart = 100, score="bde", iss=10)
p2 <- hc(asia, start=res_asia, restart = 100, score="bde", iss=100)
p3 <- hc(asia, start=res_asia, restart = 100, score="bde", iss=1000)

par(mfrow=c(1,3))
plot(p1, main="ISS: 10")
plot(p2, main="ISS: 100")
plot(p3, main="ISS: 1000")
par(mfrow=c(1,1))
library(gRain)
# Assignment 3
#true data
Init2 <- hc(asia,restart = 100)
fit <- bn.fit(x = Init2,data =asia)
#exact value
grainOBJ<-as.grain(fit)
compGrainOBj<-compile(grainOBJ)
compGrainOBj2<-setFinding(compGrainOBj,nodes = c("A"),states = c("no"))
compGrainOBj4<-setFinding(compGrainOBj,nodes = c("A","B","E"),states = c("no","yes","no"))
compGrainOBj6<-setFinding(compGrainOBj,nodes = c("A","B","E","D","T"),states = c("yes","yes","no","no","no"))
prob1 <- querygrain(compGrainOBj2,type = "marginal") 
prob2 <- querygrain(compGrainOBj4,type = "marginal")
prob3 <- querygrain(compGrainOBj6,type = "marginal")
exaxt1 <- c(prob1$B[2],prob2$D[2],prob3$S[2])


#approximate data
#approximate answer will be bad when we have many nodes since it will give us very few observation
set.seed(11827)
compAprox2<-cpdist(fit,nodes =c("A","B"),evidence = (A=="no"))
set.seed(11827)
compAprox4<-cpdist(fit,nodes =c("A","B","E","D"),evidence = (A=="no"&B=="yes"&E=="no"))
set.seed(11827)
compAprox6<-cpdist(fit,nodes =c("A","B","E","D","T","S"),evidence = (A=="yes"&B=="yes"&E=="no"&D=="no"& T=="no"))

approx <- c(prop.table(table(compAprox2$B))[2], prop.table(table(compAprox4$D))[2], prop.table(table(compAprox6$S))[2])          

data.frame("Obs.Nodes" = c("A", "A,B,E", "A,B,E,D,T"),
           "Statement " = c("No", "No,Yes,No", "y,y,n,n,n"),
           "Next included" = c("B","D","S"), "Statement"=c(
             "Yes","Yes","Yes"), "True prob"=round(exaxt1,5), 
           "Approx prob"=round(approx,5), "Diff"=round(abs(exaxt1-approx), 5))

compAprox2<-cpdist(fit,nodes =c("A","B"),evidence = (A=="no"))
prop.table(table(compAprox2))
           
compAprox2<-cpdist(fit,nodes =c("A","B"),evidence = (A=="no"))
prop.table(table(compAprox2))
           
compAprox2<-cpdist(fit,nodes =c("A","B"),evidence = (A=="no"))
prop.table(table(compAprox2))
# Assignment 4

fraction <- data.frame(Every1=NA, Every5=NA, Every10=NA,
                       Every15=NA, Every20=NA)

Every <- c(1,5,10,15,20)
burn <- c(300, 500, 1000, 5000)
countCol <- 0

for(i in Every){
  countCol <- countCol + 1
  countRow <- 0
  for(j in burn) {
    countRow <- countRow + 1
    
    DAGSGraph<-random.graph(LETTERS[1:5],method = "melancon",num=30000,every=i,burn.in=j)
cpFraction<-list()
for(k in 1:length(DAGSGraph)){
  cpFraction[[k]]<-cpdag(DAGSGraph[[k]])  
}
  fraction[countRow, countCol] <- length(unique(cpFraction))  / 30000  
  }
}
row.names(fraction) <- c("Burn300", "Burn500", "Burn1000", "Burn5000")
fraction
## NA
set.seed(311015)
TestData <- data.frame(ground=rep(c("H","A"),500,each=10),
                       surface=rep(c("Artif", "Grass"),1000, each=5),
  goals=sample(c(0,2,3,4,5), 10000, replace=TRUE,prob=c(2/9,2/6,1/6,1/6,1/6)),
  bookings=sample(c(0,1,2,3), 10000, TRUE, rep(1/4,4)))

TestData$Res = 0
TestData$Red = 0           
set.seed(311015)
for (i in 1:nrow(TestData)){
  if(TestData$ground[i] == "H" & TestData$surface[i] == "Artif"){
    TestData$Res[i] = sample(c("1","X","2"),1,prob=c(2/3,1/6,1/6))    
  }
  if(TestData$ground[i] == "A" & TestData$surface[i] == "Artif"){
    TestData$Res[i] = sample(c("1","X","2"),1,prob=c(1/6,1/6,2/3)) 
  }
  if(TestData$ground[i] == "H" & TestData$surface[i] == "Grass"){
    TestData$Res[i] = sample(c("1","X","2"),1,prob=c(3/6,2/6,1/6)) 
  }
  if(TestData$ground[i] == "A" & TestData$surface[i] == "Grass"){
    TestData$Res[i] = sample(c("1","X","2"),1,prob=c(1/6,3/6,2/6)) 
  }
  if(TestData$goals[i] == 0){
    TestData$Res[i] = "X"
  }else{
    TestData$Res[i] = TestData$Res[i] 
  }
  if(TestData$bookings[i] > 2){
    TestData$Red[i] = 1
  }else{
    TestData$Red[i] = 0
  }
}

TestData$Res <- as.factor(TestData$Res)
TestData$ground <- seq_along(levels(TestData$ground))[TestData$ground]
TestData$surface <- seq_along(levels(TestData$surface))[TestData$surface]
TestData$Res <- seq_along(levels(TestData$Res))[TestData$Res]
TestData$ground <- as.factor(TestData$ground)
TestData$surface <- as.factor(TestData$surface)
TestData$Res <- as.factor(TestData$Res)
TestData$goals <- as.factor(TestData$goals)
TestData$bookings <- as.factor(TestData$bookings)
TestData$Red <- as.factor(TestData$Red)
# Lab 1

## Assignment 1
# Initialize a random graph with 5 nodes
set.seed(311015)
init_graph <- random.graph(colnames(TestData))
# Run the score based structure learning algorithm with 100 restarts
res1 <- hc(TestData, start = init_graph, restart = 100)
all.equal(cpdag(init_graph), cpdag(res1))
plot(init_graph, main = "DAG 1")
plot(res1, main = "DAG 2, with 100 restarts")
score(init_graph, TestData)
score(res1, TestData)

## Assignment 2
# Run the score based learning algorithm for 3 different imaginary sample sizes
# initial graph is the graph returned in assignment 1

p1 <- hc(TestData, start = res1, restart = 100, score = "bde", iss = 10)
p2 <- hc(TestData, start = res1, restart = 100, score = "bde", iss = 100)
p3 <- hc(TestData, start = res1, restart = 100, score = "bde", iss = 1000)
par(mfrow = c(1, 3))
plot(p1, main = "ISS: 10")
plot(p2, main = "ISS: 100")
plot(p3, main = "ISS: 1000")
par(mfrow = c(1, 1))

# In general, number of connections increases with iss. 
# Results in a model that finds non-existing connections. 
# Overfitted to this data, this sample. 


## Assignment 3
# creates object for approximate method from res1
res1Fit <- bn.fit(res1, TestData)

# approximate method, no evidence
evi_cb <- cpdist(res1Fit, nodes = c("surface"), evidence = TRUE)
prop.table(table(evi_cb))

# Exact method, no evidence
# object for exact method
res1Fit_grain <- as.grain(res1Fit)

res1Fit_grain <- compile(res1Fit_grain)
res1grain <- setFinding(res1Fit_grain,nodes=c("surface"),
                          states = TRUE)
# Which node that is shown below is determined by [i]
querygrain(res1grain)[2]

# Clearly see difference between methods when looking at ground or surface
# for which the proportions are known. 

# Approximate method, with evidence
evi_cb <- cpdist(res1Fit, nodes = c("goals"), evidence = (surface == "2"))
prop.table(table(evi_cb))

# Exact method, with evidence
ltFit_grain <- setEvidence(res1Fit_grain,nodes=c("surface"), states = c("2") )
querygrain(ltFit_grain)[3]

# Still similar results from the respective methods #

# Approxiamte method for several different seeds
set.seed(0814)
evi_cb <- cpdist(res1Fit, nodes = c("Red"), evidence = TRUE)
prop.table(table(evi_cb))

set.seed(311015)
evi_cb <- cpdist(res1Fit, nodes = c("Red"), evidence = TRUE)
prop.table(table(evi_cb))

set.seed(1991)
evi_cb <- cpdist(res1Fit, nodes = c("Red"), evidence = TRUE)
prop.table(table(evi_cb))

set.seed(1897)
evi_cb <- cpdist(res1Fit, nodes = c("Red"), evidence = TRUE)
prop.table(table(evi_cb))

# As expected, differet results are obtained for different seeds

# Approximate versus Exact when many observed nodes
evi_cb <- cpdist(res1Fit, nodes = c("Res"),evidence=c(ground=="2" & surface=="1" &
                                           goals=="3" & bookings=="2" & Red=="0"))
prop.table(table(evi_cb))

# Exact method, no evidence
ltFit_grain <- setFinding(res1Fit_grain,nodes=c("ground","surface","goals", "bookings", "red"),states = c("2","1", "3", "2", "0"))
querygrain(ltFit_grain)[4]

# Results differs more when more observed nodes
# Problem for approximate method: Fewer and fewer observations that matches the evidence
# when the number of observed nodes increases
