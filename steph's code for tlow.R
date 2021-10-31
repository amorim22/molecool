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
  mutate_at(c("temp2", "T.low", "T.low", "mass.g"),log)

#pull(mammal2, var="temp2")

m.phy <- read.tree("mammal.tree.pruned.txt")
m.phy$tip.label <- gsub("(\\w+)_(\\w+)","\\1 \\2",m.phy$tip.label)

T.low.mod <- mammal2$T.low
names(T.low.mod) <- m.phy$tip.label
phylosig(m.phy,T.low.mod,method="lambda",test=T)

T.low.lm <- lm(log(T.low.mod)~log(mass.g),mammal2)
T.low.pgls <- gls(log(T.low.mod) ~ log(mass.g), mammal2, correlation=corBrownian(5, m.phy, form = ~species))
T.low.pgls2 <- gls(log(T.low.mod) ~ log(mass.g), mammal2, correlation=corMartins(5.05, m.phy, form = ~species))

mammal2$m.phy.pred <- predict(T.low.pgls)
mammal2$m.phy.pred2 <- predict(T.low.pgls2)

mammal2$pred <- predict(T.low.lm)

anova(T.low.lm)
anova(T.low.pgls)
anova(T.low.pgls2)

summary(T.low.lm)
summary(T.low.pgls)
summary(T.low.pgls2)
coef(T.low.lm)
coef(T.low.pgls)
coef(T.low.pgls2)

mammal2%>%
  ggplot(aes(log(mass.g),log(T.low.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=m.phy.pred))
mammal2%>%
  ggplot(aes(log(mass.g),log(T.low.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=pred))



