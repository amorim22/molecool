library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

anole <- read_csv("anole.dat.csv") 
anole.eco <- read_csv("anole.eco.csv") 

anole2 <- anole%>% 
  left_join(anole.eco)%>%
  filter(!Ecomorph%in%c("U","CH"))%>%  
  na.omit()%>%  
  print()

anole.log <- anole2%>%  
  mutate_at(c("SVL", "HTotal","PH","ArbPD"),log)

anole2%>% 
  ggplot(aes(SVL,HTotal))+geom_point()+geom_smooth(method="lm")

anole.lm <- lm(HTotal~SVL,anole2)

coef(anole.lm)

anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_abline(slope=coef(anole.lm)[2],intercept=coef(anole.lm)[1],col="blue")

SVL2 <- seq(min(anole2$SVL),max(anole2$SVL),0.1)

pred.lm <-tibble(
  SVL=SVL2,
  H.pred=predict(anole.lm,newdata = data.frame(SVL=SVL2))
)

anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_point(data=pred.lm,aes(SVL,H.pred),col="blue")

summary(anole.lm)

anole.allo <- nls(HTotal~a*SVL^b, start=list(b=1, a=1),data = anole2)

summary(anole.allo)

anole.aic <- AICc(anole.lm,anole.allo)

anole.aicw <- aicw(anole.aic$AICc)

print(anole.aicw)

logLik(anole.lm)

logLik(anole.allo)

anole.log%>%
  ggplot(aes(HTotal,SVL,col=Ecomorph2))+geom_point()+geom_smooth(method="lm")

anole.log.eco.lm <- lm(HTotal~SVL*Ecomorph2,anole.log)
summary(anole.log.eco.lm)

anova(anole.log.eco.lm)

anole.log.lm  <- lm(HTotal~SVL,anole.log)
anova(anole.log.lm)

anole.log.aic <- AICc(anole.log.lm,anole.log.eco.lm)
aicw(anole.log.aic$AICc)

anole.log <- anole.log %>%
  mutate(res=residuals(anole.log.lm))

anole.log%>%
  ggplot(aes(Ecomorph2,res))+geom_point()

p.eco <- anole.log%>%
  ggplot(aes(x=Ecomorph2,y=res)) +geom_boxplot()
print(p.eco)

p.eco+ geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)

anole.tree <- read.tree("anole.tre")
plot(anole.tree,cex=0.4)

pgls.BM1 <- gls(HTotal ~SVL, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

pgls.BM2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

pgls.OU1 <- gls(HTotal ~SVL, correlation = corMartins(0,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

pgls.OU2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corMartins(0,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

anole.phylo.aic <- AICc(pgls.BM1,pgls.BM2,pgls.OU1,pgls.OU2)
aicw(anole.phylo.aic$AICc)

anova(pgls.BM2) 

anole.log <- anole.log%>%
  mutate(phylo.res=residuals(pgls.BM2))

p.eco.phylo <- anole.log%>%
  ggplot(aes(x=Ecomorph2,y=phylo.res)) +geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)

print(p.eco.phylo)

anole.log%>%
  dplyr::select(Ecomorph2,res,phylo.res)%>%
  pivot_longer(cols=c("res","phylo.res"))%>%
  print%>%
  ggplot(aes(x=Ecomorph2,y=value)) +geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)+facet_grid(name~.,scales = "free_y")+ylab("residual")

#CPK: Still lots of stuff you don't need. Most above is examples from the desription [-1]
anole.log

#CPK: 2/3 right? 
anole.log.PH.lm <- lm(HTotal~SVL+PH,anole.log)
summary(anole.log.PH.lm)
anova(anole.log.PH.lm)

anole.log.ArbPD.lm <- lm(HTotal~SVL+ArbPD,anole.log)
summary(anole.log.ArbPD.lm)
anova(anole.log.ArbPD.lm)

anole.log <- anole.log%>%
  mutate(PH.res=residuals(anole.log.PH.lm))

anole.log <- anole.log%>%
  mutate(ArbPD.res=residuals(anole.log.ArbPD.lm))

p.PH <- anole.log%>%
  ggplot(aes(x=Ecomorph2, y=PH.res))+geom_boxplot()+stat_summary(fun=mean, geom="point", size=3)
print(p.PH)

p.ArbPD <- anole.log%>%
  ggplot(aes(x=Ecomorph2, y=ArbPD.res))+geom_boxplot()+stat_summary(fun=mean, geom="point", size=3)
print(p.ArbPD)

#CPK didn't we need to plot the the HL against the residuals of the covariates (PH and PD). No mention of ecomorph. Prompt #3: Explore how both perch diameter and height effect the hindlimb-SVL relationship by plotting the residuals of your simple linear models against these covariates. This will require mutating a data tibble to include residuals from both models. Please produce two separate plots. Code to do so below. -1 pt for not including the right regression plots. It simple as this . . .

 anole.log%>%
  ggplot(aes(x=PH.res,y=HTotal))+geom_point()+geom_smooth(method="lm")


 anole.log%>%
   ggplot(aes(x=ArbPD.res,y=HTotal))+geom_point()+geom_smooth(method="lm")


 #CPK Well done with the gls modeling and AIC
pgls.BM3 <- gls(HTotal ~SVL * PH, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

pgls.BM4 <- gls(HTotal ~SVL * ArbPD, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

pgls.BM5 <- gls(HTotal ~SVL * PH * ArbPD, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

anole.PHArbPD.aic <- AICc(pgls.BM3,pgls.BM4,pgls.BM5)
aicw(anole.PHArbPD.aic$AICc)

anova(pgls.BM4)

anole.log <- anole.log%>%
  mutate(phylo2.res=residuals(pgls.BM4))

p.ArbPD.phylo <- anole.log%>%
  ggplot(aes(x=ArbPD,y=phylo2.res)) +geom_point() +geom_smooth(method="lm")
print(p.ArbPD.phylo)

anole.log%>%
  dplyr::select(ArbPD,res,phylo2.res)%>%
  pivot_longer(cols=c("res","phylo2.res"))%>%
  print%>%
  ggplot(aes(x=ArbPD,y=value)) +geom_point() +geom_smooth(method="lm")+stat_summary(fun=mean, geom="point", size=1)+facet_grid(name~.,scales = "free_y")+ylab("residual")

#CPK. But again, the last prompt didn't ask about ecomorph, just to produce a plot of your own design that concisely visualizes the effect of your covariate(s) on the hindlimb residuals of the best fitting PGLS model.As simple as this.

anole.log%>%
  ggplot(aes(x=phylo2.res,HTotal)) +geom_point() +geom_smooth(method="lm")

#CPK. Very good work. Just be careful not to include what we dont' need and address the prompts with the appropriate actions. [8/10]