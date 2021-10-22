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

mammal.lm <- lm(temp~mass.g,mammal)

mammal.allo <- nls(temp~a*mass.g^b, start=list(b=0.1, a=0.1),data = mammal)

summary(mammal.allo)

mammal.aic <- AICc(mammal.lm,mammal.allo)

mammal.aicw <- aicw(mammal.aic$AICc)

print(mammal.aicw)

mammal.Order.lm <- lm(temp~mass.g*Order,mammal)
summary(mammal.Order.lm)

anova(mammal.Order.lm)

anova(mammal.lm)

mammal.Order.aic <- AICc(mammal.lm,mammal.Order.lm)
aicw(mammal.Order.aic$AICc)

mammal <- mammal %>%
  mutate(res=residuals(mammal.Order.lm))

mammal%>%
  ggplot(aes(Order,res))+geom_point()

p.Order <- mammal%>%
  ggplot(aes(x=Order,y=res)) +geom_boxplot()
print(p.Order)

p.Order+ geom_boxplot() +stat_summary(fun=mean, geom="point", size=2)