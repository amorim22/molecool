library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

#load data
mammal2 <- read_csv("mammal.temp2.csv")
mammal2 <- mammal2 %>% 
  mutate(temp2=abs(T.high - T.low))

mammal2.log <- mammal2%>%
  mutate_at(c("temp2", "T.high", "T.low", "mass.g"),log)

m.phy <- read.tree("mammal.tree.pruned.txt")
m.phy$tip.label <- gsub("(\\w+)_(\\w+)","\\1 \\2",m.phy$tip.label)

mammal2.T.high.lm <- lm(T.high~mass.g,mammal2)
mammal2.T.low.lm <- lm(T.low~mass.g,mammal2)
mammal2.temp2.lm <- lm(temp2~mass.g,mammal2)

mammal2.T.high.allo <- nls(T.high~a*mass.g^b, start=list(b=0.1, a=0.1),data = mammal2.log)
mammal2.T.low.allo <- nls(T.low~a*mass.g^b, start=list(b=0.1, a=0.1),data = mammal2.log)
mammal2.temp2.allo <- nls(temp2~a*mass.g^b, start=list(b=0.1, a=0.1),data = mammal2.log) 

mammal2.T.high.aic <- AICc(mammal2.T.high.lm,mammal2.T.high.allo)
mammal2.T.low.aic <- AICc(mammal2.T.low.lm,mammal2.T.low.allo)
mammal2.temp2.aic <- AICc(mammal2.temp2.lm,mammal2.temp2.allo)

mammal2.T.high.aicw <- aicw(mammal2.T.high.aic$AICc)
mammal2.T.low.aicw <- aicw(mammal2.T.low.aic$AICc)
mammal2.temp2.aicw <- aicw(mammal2.temp2.aic$AICc)

print(mammal2.T.high.aic)
print(mammal2.T.low.aic)
print(mammal2.temp2.aic)

print(mammal2.T.high.aicw)
print(mammal2.T.low.aicw)
print(mammal2.temp2.aicw)

#PGLS under BM, T.high
pgls.BM1 <- gls(T.high ~mass.g, correlation = corBrownian(5,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

#PGLS under BM, T.low
pgls.BM2 <- gls(T.low ~mass.g, correlation = corBrownian(5,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

#PGLS under BM, temp2
pgls.BM3 <- gls(temp2 ~mass.g, correlation = corBrownian(5,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

#PGLS under OU, T.high
pgls.OU1 <- gls(T.high ~mass.g, correlation = corMartins(5,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

#PGLS under OU, T.low
pgls.OU2 <- gls(T.low ~mass.g, correlation = corMartins(5.05,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

#PGLS under OU, temp2
pgls.OU3 <- gls(temp2 ~mass.g, correlation = corMartins(5,phy = m.phy,form=~species),data = mammal2.log, method = "ML")

mammal2.phylo.aic <- AICc(pgls.BM1,pgls.BM2,pgls.BM3,pgls.OU1,pgls.OU2,pgls.OU3)
mammal2.phylo.aicw <- aicw(mammal2.phylo.aic$AICc)
print(mammal2.phylo.aicw)

anova(pgls.OU1)
anova(pgls.OU2)
anova(pgls.OU3)

#T.high
T.high.mod <- mammal2$T.high
names(T.high.mod) <- m.phy$tip.label
phylosig(m.phy,T.high.mod,method="lambda",test=T)

T.high.lm <- lm(log(T.high.mod)~log(mass.g),mammal2)
T.high.pgls <- gls(log(T.high.mod) ~ log(mass.g), mammal2, correlation=corBrownian(5, m.phy, form = ~species))
T.high.pgls2 <- gls(log(T.high.mod) ~ log(mass.g), mammal2, correlation=corMartins(5, m.phy, form = ~species))

T.high.aic <- AICc(T.high.pgls, T.high.pgls2)
T.high.aicw <- aicw(T.high.aic$AICc)

print(T.high.aicw)

mammal2$m.phy.T.high.pred <- predict(T.high.pgls)
mammal2$m.phy.T.high.pred2 <- predict(T.high.pgls2)

mammal2$T.high.pred <- predict(T.high.lm)

anova(T.high.lm)
anova(T.high.pgls)

summary(T.high.lm)
summary(T.high.pgls)

coef(T.high.lm)
coef(T.high.pgls)

#T.low
T.low.mod <- mammal2$T.low
names(T.low.mod) <- m.phy$tip.label
phylosig(m.phy,T.low.mod,method="lambda",test=T)

T.low.lm <- lm(log(T.low.mod)~log(mass.g),mammal2)
T.low.pgls <- gls(log(T.low.mod) ~ log(mass.g), mammal2, correlation=corBrownian(5, m.phy, form = ~species))
T.low.pgls2 <- gls(log(T.low.mod) ~ log(mass.g), mammal2, correlation=corMartins(5.05, m.phy, form = ~species))

T.low.aic <- AICc(T.low.pgls, T.low.pgls2)
T.low.aicw <- aicw(T.low.aic$AICc)

print(T.low.aicw)

mammal2$m.phy.T.low.pred <- predict(T.low.pgls)
mammal2$m.phy.T.low.pred2 <- predict(T.low.pgls2)

mammal2$T.low.pred <- predict(T.low.lm)

anova(T.low.lm)
anova(T.low.pgls)

summary(T.low.lm)
summary(T.low.pgls)

coef(T.high.lm)
coef(T.low.pgls)

#temp2
m.phy <- read.tree("mammal.tree.pruned.txt")
m.phy$tip.label <- gsub("(\\w+)_(\\w+)","\\1 \\2",m.phy$tip.label)

temp2.mod <- mammal2$temp2
names(temp2.mod) <- m.phy$tip.label
phylosig(m.phy,temp2.mod,method="lambda",test=T)

temp2.lm <- lm(log(temp2.mod)~log(mass.g),mammal2)
temp2.pgls <- gls(log(temp2.mod) ~ log(mass.g), mammal2, correlation=corBrownian(5, m.phy, form = ~species))
temp2.pgls2 <- gls(log(temp2.mod) ~ log(mass.g), mammal2, correlation=corMartins(5.05, m.phy, form = ~species))

temp2.aic <- AICc(temp2.pgls, temp2.pgls2)
temp2.aicw <- aicw(temp2.aic$AICc)

print(temp2.aicw)

mammal2$m.phy.temp2.pred <- predict(temp2.pgls)
mammal2$m.phy.temp2.pred2 <- predict(temp2.pgls2)

mammal2$temp2.pred <- predict(temp2.lm)

anova(temp2.lm)
anova(temp2.pgls)

summary(temp2.lm)
summary(temp2.pgls)

coef(temp2.lm)
coef(temp2.pgls)

#graph zone
mammal2%>%
  ggplot(aes(log(mass.g),log(T.high.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=m.phy.T.high.pred))
mammal2%>%
  ggplot(aes(log(mass.g),log(T.high.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=m.phy.T.high.pred2))

mammal2%>%
  ggplot(aes(log(mass.g),log(T.low.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=m.phy.T.low.pred))
mammal2%>%
  ggplot(aes(log(mass.g),log(T.low.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=m.phy.T.low.pred2))

mammal2%>%
  ggplot(aes(log(mass.g),log(temp2.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=m.phy.temp2.pred))
mammal2%>%
  ggplot(aes(log(mass.g),log(temp2.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=m.phy.temp2.pred2))
mammal2%>%
  ggplot(aes(log(mass.g),log(temp2.mod)))+geom_point()+geom_smooth(method="lm",se = F)+geom_line(aes(y=temp2.pred))
