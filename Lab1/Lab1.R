
library(bnlearn)
library(gRain)
library(ggplot2)

data(learning.test)
data(asia)
data(alarm)

# 1
# learning.test data
lt1 <- hc(learning.test, restart = 0)
plot(lt1)
vstructs(lt1)
cpdag(lt1)

score(lt1, data=learning.test ,type="bde")

lt2 <- hc(learning.test, start=lt1, restart=100)
plot(lt2)
all.equal(lt1, lt2)
vstructs(lt2)
cpdag(lt2)

# asia data set 
set.seed(3110)
as1 <- hc(asia, restart = 0)
plot(as1)
vstructs(as1)
cpdag(as1)

as2 <- hc(asia, start=as1, restart=100)
plot(as2)
all.equal(as1, as2)
vstructs(as2)
cpdag(as2)

# 2
# asia data set
as1_2 <- hc(asia, score = "bde", restart = 10)
vstructs(as1_2)

sSize <- seq(10, 1000, 10)
ScoreF <- data.frame(iss=seq(10, 1000, 10), score=NA)
for(i in 1:100){
  ScoreF[i,2] <- score(as1_2, data=asia ,type="bde", iss=ScoreF[i,1])
}
ggplot(ScoreF, aes(x=score))+geom_histogram(binwidth=400,alpha=0.5,fill="blue")+
  theme_light() + scale_x_continuous(breaks=c(seq(-11000,-14000,-400)))

as1_2 <- hc(asia, score = "bde", restart = 10, iss=10)
as2_2 <- hc(asia, score = "bde", restart = 10, iss=100)
as3_2 <- hc(asia, score = "bde", restart = 10, iss=500)
as4_2 <- hc(asia, score = "bde", restart = 10, iss=1000)
plot(as1_2, main="ISS=10")
plot(as2_2, main="ISS=100")
plot(as3_2, main="ISS=500")
plot(as4_2)

score(as1_2, data=asia ,type="bde", iss=1000)
par(mfrow=c(2,2))


# True DAG
res = empty.graph(names(learning.test))
modelstring(res) = "[A][C][F][B|A][D|A:C][E|B:F]"
plot(res)

# 3
plot(lt2)
ltFit <- bn.fit(lt2, learning.test)

cpdist(ltFit, nodes = c("A", "B", "C", "D", "E", "F"), evidence = (C == "c"))

evi_cb <- cpdist(ltFit, nodes = c("A"), evidence = TRUE)
prop.table(table(evi_cb))

# Jämföra med exakt och när jag har en eller två observationer (några)

ltFit_grain <- as.grain(ltFit)

ltFit_grain <- compile(ltFit_grain)
ltFit_grain <- setFinding(ltFit_grain,nodes=c("A"),
           states = TRUE)
querygrain(ltFit_grain)

# start over
# approximate method, no evidence
evi_cb <- cpdist(ltFit, nodes = c("A"), evidence = TRUE)
prop.table(table(evi_cb))

# Exact method, no evidence
ltFit_grain <- as.grain(ltFit)

ltFit_grain <- compile(ltFit_grain)
ltFit_grain <- setFinding(ltFit_grain,nodes=c("A"),
                          states = TRUE)
querygrain(ltFit_grain)[1]

# Approximate method, with evidence
evi_cb <- cpdist(ltFit, nodes = c("A"), evidence = (D == "c"))
prop.table(table(evi_cb))

# Exact method, with evidence
ltFit_grain <- as.grain(ltFit)
ltFit_grain <- compile(ltFit_grain)

ltFit_grain <- setEvidence(ltFit_grain,nodes=c("D"), states = c("c") )
querygrain(ltFit_grain)[1]

# Approxiamte method for several different seeds
set.seed(0814)
evi_cb <- cpdist(ltFit, nodes = c("A"), evidence = TRUE)
prop.table(table(evi_cb))

set.seed(311015)
evi_cb <- cpdist(ltFit, nodes = c("A"), evidence = TRUE)
prop.table(table(evi_cb))

set.seed(1991)
evi_cb <- cpdist(ltFit, nodes = c("A"), evidence = TRUE)
prop.table(table(evi_cb))

set.seed(1897)
evi_cb <- cpdist(ltFit, nodes = c("A"), evidence = TRUE)
prop.table(table(evi_cb))

# Approximate versus Exact when many observed nodes
# approximate method, no evidence
evi_cb <- cpdist(ltFit, nodes = c("B"),evidence=c(A=="b" & C=="a" & D=="c" 
                                                      & F=="a" & E=="a"))
prop.table(table(evi_cb))

# Exact method, no evidence
ltFit_grain <- as.grain(ltFit)

ltFit_grain <- compile(ltFit_grain)
ltFit_grain <- setFinding(ltFit_grain,nodes=c("A","C", "D", "F", "E"),
                          states = c("b", "a", "c", "a", "a"))
querygrain(ltFit_grain)[2]


# 4
gr_100 <- random.graph(nodes=c("A","B","C","D","E"), num = 30000, method="melancon",
                       burn.in=10000, every=0.1)

gr_list <- list()
for(i in 1:length(gr_100)){
  gr_list[[i]] <- cpdag(gr_100[[i]])  
}

length(unique(gr_list))

