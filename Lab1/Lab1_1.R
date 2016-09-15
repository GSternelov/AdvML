## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
library(bnlearn)
library(gRain)
library(ggplot2)

data(learning.test)
data(asia)
options("scipen"=10)

## ---- echo=FALSE, fig.height=3.5, fig.width=8, fig.align='center'--------
# 1
# learning.test data
set.seed(1897)
lt1 <- hc(learning.test, restart = 0)
par(mar=c(1.1,4.1,2.1,2.1))
plot(lt1, main="DAG1 (data: learning.test,no initial structure, 0 restarts)")
ltS1 <- round(score(lt1, data=learning.test ,type="bde"),2)
#vstructs(lt1)
#cpdag(lt1)

## ---- echo=FALSE, fig.height=3.5, fig.width=8, fig.align='center'--------
set.seed(1897)
lt2 <- hc(learning.test, start=lt1, restart=100)
par(mar=c(1.1,4.1,4.1,2.1))
plot(lt2, main="DAG2 (data: learning.test, DAG1 as initial structure, 100 restarts)")
ltS2 <- round(score(lt2, data=learning.test ,type="bde"),2)
#vstructs(lt2)
#cpdag(lt2)

## ---- echo=TRUE----------------------------------------------------------
all.equal(lt1, lt2)

## ---- echo=FALSE, fig.height=3.5, fig.width=8, fig.align='center'--------
# asia data set 
set.seed(3110)
as1 <- hc(asia, restart = 0)
par(mar=c(1.1,4.1,4.1,2.1))
plot(as1, main="DAG3 (data: Asia, no initial structure, 0 restarts)")
aS1 <- round(score(as1, data=asia ,type="bde"),2)
#vstructs(as1)
#cpdag(as1)

## ---- echo=FALSE, fig.height=3.5, fig.width=8, fig.align='center'--------
set.seed(3110)
as2 <- hc(asia, start=as1, restart=100)
par(mar=c(1.1,4.1,4.1,2.1))
plot(as2, main="DAG4 (data: Asia, DAG3 as initial structure, 100 restarts)")
aS2 <- round(score(as2, data=asia ,type="bde"),2)
#vstructs(as2)
#cpdag(as2)

## ----echo=TRUE-----------------------------------------------------------
all.equal(as1, as2)

## ---- echo=FALSE, fig.height=3.5, fig.width=6, fig.align='center'--------
# 2
# asia data set
as1_2 <- hc(asia, score = "bde", restart = 10)

sSize <- seq(10, 1000, 10)
ScoreF <- data.frame(iss=seq(10, 1000, 10), score=NA)
for(i in 1:100){
  ScoreF[i,2] <- score(as1_2, data=asia ,type="bde", iss=ScoreF[i,1])
}
ggplot(ScoreF, aes(x=score))+geom_histogram(binwidth=400,alpha=0.5,fill="blue")+
  theme_light() + scale_x_continuous(breaks=c(seq(-11000,-14000,-400))) + theme(axis.title.y = element_text(angle=0)) + ggtitle("DBeu score for different imaginary sample sizes")

## ---- echo=FALSE, fig.height=6, fig.width=8, fig.align='center'----------
par(mar=c(1.1,4.1,4.1,2.1))
par(mfrow=c(2,2))
as1_2 <- hc(asia, score = "bde", restart = 10, iss=10)
as2_2 <- hc(asia, score = "bde", restart = 10, iss=100)
as3_2 <- hc(asia, score = "bde", restart = 10, iss=500)
plot(as1_2, main="ISS=10")
plot(as2_2, main="ISS=100")
plot(as3_2, main="ISS=500")
# True DAG
res = empty.graph(names(learning.test))
modelstring(res) = "[A][C][F][B|A][D|A:C][E|B:F]"
plot(res, main="True DAG")


## ---- echo=FALSE---------------------------------------------------------
ltFit <- bn.fit(lt2, learning.test)
# approximate method, no evidence
evi_cb <- cpdist(ltFit, nodes = c("A"), evidence = TRUE)
prop.table(table(evi_cb))

# Exact method, no evidence
ltFit_grain <- as.grain(ltFit)
ltFit_grain <- compile(ltFit_grain)

ltFit_grain1 <- setFinding(ltFit_grain,nodes=c("A"),
                          states = TRUE)
querygrain(ltFit_grain1)[1]

## ---- echo=FALSE---------------------------------------------------------
# Approximate method, with evidence
evi_cb <- cpdist(ltFit, nodes = c("A"), evidence = (D == "c"))
prop.table(table(evi_cb))

# Exact method, with evidence
ltFit_grain2 <- setEvidence(ltFit_grain,nodes=c("D"), states = c("c") )
querygrain(ltFit_grain2)[1]

## ---- echo=FALSE---------------------------------------------------------
# Approximate method for several different seeds
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

## ---- echo=FALSE---------------------------------------------------------
# Approximate versus Exact when many observed nodes
# approximate method, no evidence
set.seed(0800014)
evi_cb <- cpdist(ltFit, nodes = c("B"),evidence=c(A=="b" & C=="a" & D=="c" 
                                                      & F=="a" & E=="a"))
prop.table(table(evi_cb))

# Exact method, no evidence
ltFit_grain3 <- setFinding(ltFit_grain,nodes=c("A","C", "D", "F", "E"),
                          states = c("b", "a", "c", "a", "a"))
querygrain(ltFit_grain3)[2]

## ---- echo=FALSE---------------------------------------------------------
gr_1 <- random.graph(nodes=c("A","B","C","D","E"), num = 30000, method="melancon")
gr_list <- list()
for(i in 1:length(gr_1)){
  gr_list[[i]] <- cpdag(gr_1[[i]])  
}
Lgr_1 <- length(unique(gr_list))

gr_2 <- random.graph(nodes=c("A","B","C","D","E"), num = 30000, method="melancon",
                       burn.in=1000)
gr_list <- list()
for(i in 1:length(gr_2)){
  gr_list[[i]] <- cpdag(gr_2[[i]])  
}
Lgr_2 <- length(unique(gr_list))

gr_3 <- random.graph(nodes=c("A","B","C","D","E"), num = 30000, method="melancon",
                       burn.in=1000, every=0.1)
gr_list <- list()
for(i in 1:length(gr_3)){
  gr_list[[i]] <- cpdag(gr_3[[i]])  
}
Lgr_3 <- length(unique(gr_list))

gr_4 <- random.graph(nodes=c("A","B","C","D","E"), num = 30000, method="melancon",
                       burn.in=1000, every=10)
gr_list <- list()
for(i in 1:length(gr_4)){
  gr_list[[i]] <- cpdag(gr_4[[i]])  
}
Lgr_4 <- length(unique(gr_list))

## ----code=readLines(knitr::purl("C:\\Users\\Gustav\\Documents\\AdvML\\Lab1\\Lab1_1.Rmd",documentation = 1)), eval = FALSE, tidy=TRUE----
## 
## 

