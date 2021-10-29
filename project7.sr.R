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

temp2.mod <- mammal2$temp2
names(temp2.mod) <- m.phy$tip.label
phylosig(m.phy,temp2.mod,method="lambda",test=T)

temp2.lm <- lm(log(temp2.mod)~log(mass.g),mammal2)
temp2.pgls <- gls(log(temp2.mod) ~ log(mass.g), mammal2, correlation=corBrownian(5, m.phy, form = ~species))

mammal2$m.phy.pred <- predict(temp2.pgls)
mammal2$pred <- predict(temp2.lm)

anova(temp2.lm)
anova(temp2.pgls)

summary(temp2.lm)
summary(temp2.pgls)
coef(temp2.lm)
coef(temp2.pgls)

mammal2%>%
  ggplot(aes(log(mass.g),log(temp2.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=m.phy.pred))
mammal2%>%
  ggplot(aes(log(mass.g),log(temp2.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=pred))



