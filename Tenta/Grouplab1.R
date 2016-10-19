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
