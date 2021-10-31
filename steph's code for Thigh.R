library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

mammal2 <- read_csv("mammal.temp2.csv")
mammal2 <- mammal2 %>% 
  mutate(temp2=T.high - T.low)

mammal2.log <- mammal2%>%
  mutate_at(c("temp2", "T.high", "T.low", "mass.g"),log)

#pull(mammal2, var="temp2")

m.phy <- read.tree("mammal.tree.pruned.txt")
m.phy$tip.label <- gsub("(\\w+)_(\\w+)","\\1 \\2",m.phy$tip.label)

T.high.mod <- mammal2$T.high
names(T.high.mod) <- m.phy$tip.label
phylosig(m.phy,T.high.mod,method="lambda",test=T)

T.high.lm <- lm(log(T.high.mod)~log(mass.g),mammal2)
T.high.pgls <- gls(log(T.high.mod) ~ log(mass.g), mammal2, correlation=corBrownian(5, m.phy, form = ~species))
T.high.pgls2 <- gls(log(T.high.mod) ~ log(mass.g), mammal2, correlation=corMartins(5, m.phy, form = ~species))

mammal2$m.phy.pred <- predict(T.high.pgls)
mammal2$m.phy.pred2 <- predict(T.high.pgls2)

mammal2$pred <- predict(T.high.lm)

anova(T.high.lm)
anova(T.high.pgls)
anova(T.high.pgls2)

summary(T.high.lm)
summary(T.high.pgls)
summary(T.high.pgls2)
coef(T.high.lm)
coef(T.high.pgls)
coef(T.high.pgls2)

mammal2%>%
  ggplot(aes(log(mass.g),log(T.high.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=m.phy.pred))
mammal2%>%
  ggplot(aes(log(mass.g),log(T.high.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=pred))



