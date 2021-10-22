library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

mammal <- read_csv("mammal.temp.csv")
mammal <- mammal %>% mutate(temp=T.high - T.low)

print(mammal)

mammal.log <- mammal%>%
  mutate_at(c("temp", "mass.g"),log)

mammal.lm <- lm(temp~mass.g,mammal)

mammal.allo <- nls(temp~a*mass.g^b, start=list(b=0.1, a=0.1),data = mammal)

summary(mammal.allo)

mammal.aic <- AICc(mammal.lm,mammal.allo)

mammal.aicw <- aicw(mammal.aic$AICc)

print(mammal.aicw)

mammal.log%>%
  ggplot(aes(mass.g,temp,col=Order))+geom_point()+geom_smooth(method="lm")

mammal.log.Order.lm <- lm(temp~mass.g*Order,mammal.log)
summary(mammal.log.Order.lm)
anova(mammal.log.Order.lm)

mammal.log.lm  <- lm(temp~mass.g,mammal.log)
anova(mammal.log.lm)

mammal.log.aic <- AICc(mammal.log.lm,mammal.log.Order.lm)
aicw(mammal.log.aic$AICc)

mammal.log <- mammal.log %>%
  mutate(res=residuals(mammal.log.Order.lm))

mammal.log%>%
  ggplot(aes(Order,res))+geom_point()

p.Order <- mammal.log%>%
  ggplot(aes(x=Order,y=res))+geom_boxplot()
print(p.Order)

p.Order+geom_boxplot() +stat_summary(fun=mean, geom="point", size=2)
