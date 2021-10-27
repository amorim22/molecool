library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

mammal2 <- read_csv("mammal.temp2.csv")
mammal2 <- mammal2 %>% mutate(temp2=T.high - T.low)

mammal2.log <- mammal2%>%
mutate_at(c("temp2", "T.high", "T.low", "mass.g"),log)

mammal2.T.high.allo <- nls(T.high~a*mass.g^b, start=list(b=0.1, a=0.1),data = mammal2.log)
mammal2.T.low.allo <- nls(T.low~a*mass.g^b, start=list(b=0.1, a=0.1),data = mammal2.log)
mammal2.temp2.allo <- nls(temp2~a*mass.g^b, start=list(b=0.1, a=0.1),data = mammal2.log) 

summary(mammal2.T.high.allo)
summary(mammal2.T.low.allo)
summary(mammal2.temp2.allo)

m.phy <- read.tree("mammal.tree.pruned.txt")
m.phy$tip.label <- gsub("(\\w+)_(\\w+)","\\1 \\2",m.phy$tip.label)

#PGLS under BM, T.high
pgls.BM1 <- gls(T.high ~mass.g, correlation = corBrownian(1,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

#PGLS under BM, T.low
pgls.BM2 <- gls(T.low ~mass.g, correlation = corBrownian(1,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

#PGLS under BM, temp2
pgls.BM3 <- gls(temp2 ~mass.g, correlation = corBrownian(1,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

#PGLS under OU, T.high
pgls.OU1 <- gls(T.high ~mass.g, correlation = corMartins(0,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

#PGLS under OU, T.low
pgls.OU2 <- gls(T.low ~mass.g, correlation = corMartins(0,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

#PGLS under OU, temp2
pgls.OU3 <- gls(temp2 ~mass.g, correlation = corMartins(0,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

mammal2.phylo.aic <- AICc(pgls.BM1,pgls.BM2,pgls.BM3,pgls.OU1,pgls.OU2,pgls.OU3)
aicw(mammal2.phylo.aic$AICc)

anova(pgls.BM1)
anova(pgls.BM2)
anova(pgls.BM3)
