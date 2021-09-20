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

#prepare data files

anole2 <- anole%>% #save filtered tibble to anole2
left_join(anole.eco)%>% #mutate operation 
  filter(!Ecomorph%in%c("U","CH"))%>% 
# filter to keep the data of ecomorphs that are not in c("U","CH") using !
  na.omit()%>%
#omit rows on tibble with missing values 
  print()

anole.log <- anole2%>%
  mutate_at(c("SVL","HTotal","PH","ArbPD"),log)
#mutate collective tibble with data - change column containing size and ecological data to log trnasformation of values 
#continuous data - log transformation convert data to proportional representations
#mutate_at - mutate columns in place, replace each with names wrapped in c()
#second parameter in mutate_at sepcify what we do to column - natural log 

anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_smooth(method="lm")
#does hind limb(HTotal) vary with size(SVL)
#smooth line on a linear model - method="lm"

anole.lm <- lm(HTotal~SVL,anole2)
coef(anole.lm)
#evaluate the fit using simple linear model using lm() function 
#HTotal~SVL: model formulae, variable to left is predicted by one or more variable to the right 
#after the ~, specify what data the variable are coming from - this case anole2

anole2%>%
ggplot(aes(SVL,HTotal))+geom_point()+geom_abline(slope=coef(anole.lm)[2],intercept=coef(anole.lm)[1],col="blue")

#now use range of snout vent length rather than have geom_smooth do it 
#make tibble with predictions

SVL2 <- seq(min(anole2$SVL),max(anole2$SVL),0.1)

pred.lm <- tibble(
  SVL=SVL2,
  H.pred=predict(anole.lm,newdata = data.frame(SVL=SVL2))
)

anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_point(data=pred.lm,aes(SVL,H.pred),col="blue")

summary(anole.lm)
#summary retrieve info abt the model - include estimated response of HTotal to SVL: ex. slope
#retrieve r^2 - pretty good model - 0.93 
#use an allometric model this time instead - use non-linear lesat squares: nls() - approximate relationship by minimizing sum of squares of the residuals 
#model is fit based on non-linear relationship 

anole.allo <- nls(HTotal~a*SVL^b, start = list(b=1, a=1),data = anole2)

summary(anole.allo)
#need to specify starting values for parameters in the model with a list - "start"
#"nls()" conducts a search of model parameter values to find those that reduce the sum of squares of the residuals 
#starting points estimated on basis that scaling coefficient b is probably around 1 - data fitted a linear model well, intercept a probably is close to but not zero 
#a and b values are significant, but no r^2 value - model not applicable to non-linear regression
#more than we want to know if the HTotal varies with SVL or not, we want to know if the model is linear or allometric 
#carry out likelihood ratio tests - compare likelihood of two nested models 
#statistical property of likelihood - defined as probability - model and set of parameter values, of obtaining particular set of data
#whichever model with higher likelihood describe our data better - better mathematical approximation
#to compare models of any type nested or not - Akaike Information Criterion (AIC) 
#description of model fit compares the likelihood of model against the number of parameters estimated 
#AIC - favor model that provides the best fit to the data with as few parameters as possible 
#use AIC for smaller sample size - AICc 
#also calculate relative support of each model using AIC weights (AICw)

anole.aic <- AICc(anole.lm,anole.allo)
anole.aicw <- aicw(anole.aic$AICc)
print(anole.aicw)
#anole.allo has lower AIC score - indicate better fit 
#relative fit is 3x better, but difference in AICc is somewhat smaller (change in AIC = 2.2)
#change of AIC of less than 4 indicate roughly equivalent models - little difference btw our linear and allometric 
#so we say that allometry and isometry are roughly equivalent models 

#now consider effect of ecomorph on the hindlimb-SVL relationship 
#use log-transformed data - visualize hindlimb-SVL relationship for each ecomorph in ggplot
anole.log%>%
  ggplot(aes(HTotal,SVL,col=Ecomorph2))+geom_point()+geom_smooth(method="lm")
#colors=Ecomorph2 - colors all the added geometries according to column values in Ecomorph2 
#group data to compute multiple regression lines for geom_smooth 

anole.log.eco.lm <- lm(HTotal~SVL*Ecomorph2,anole.log)
summary(anole.log.eco.lm)
anova(anole.log.eco.lm)
#ANOVA test of our models - more precisely defined, two-way analysis of covariance 
#assess effect of categorical variable (Ecomorph2) in the context of how HTotal covaries with SVL 
#results indicate that we should reject the null hypothesis that ecomorph groups do not have separate hindlimb-SVL relationships 
#want to know if adding Ecomorph2 parameter results in better fit compared to simple model 
#use log-transformed data and compare with AIC and AICw to more complicated model 

anole.log.lm <- lm(HTotal~SVL,anole.log)
anova(anole.log.lm)

anole.log.aic <- AICc(anole.log.lm,anole.log.eco.lm)
aicw(anole.log.aic$AICc)
#model with Ecomorph2 shows better fit 

#now see how much ecomorph varies among anoles in HTotal-SVL relationship
#residual - represent how much any one data point varies from prediction 
#have a global model predicting HTotal vs. SVL on all anole species, determine how much each species deviate from prediction
#residual of each species, look for pattern - plot residual vs. ecomorph 

#compute residual based on global anole.log.lm model
anole.log <- anole.log %>%
  mutate(res=residuals(anole.log.lm))
#used mutate to establish new res column that contains the residuals 
#redefined anole.log as a tibble containing the new column 

anole.log%>%
  ggplot(aes(Ecomorph2,res))+geom_point()
#plot residual against ecomorph2 
#establish x variable as ecomorph2, y as residual 
#model deviations are greatest in twig and trunk-ground ecomorphs 
#include median residual as well - more mutaiton of data or ggplot 

p.eco <- anole.log%>%
  ggplot(aes(x=Ecomorph2,y=res))+geom_boxplot()
  print(p.eco)
#visualize mean instead of median 
#add geometries to existing plot 

p.eco+geom_boxplot()+stat_summary(fun=mean,geom="point",size=3)
#added point with stat_summary(): apply summarizing function to already established group, plot values of summary with specified geometry
#applied function: fun=mean, specified point with size 3

#need to account for phylogeny - use phylogenetic comparative methods - PCM
#phylogenetic generalized least squares (PGLS) may be most flexible and widely used
#need to establish covriation matrix with tree, assume that characters have evolved in this tree under some evolutionary model

anole.tree <- read.tree("anole.tre")
plot(anole.tree,cex=0.4)

#evolutionary model - under what process or mode does a trait(S) evolve over a tree
#Brownian motion model of character evolution - random-walk 
#value of character walks get higher or lower form previous value - sum of changes describe change in value 
#how much the value can change over time - sigma squared - higher=greater change 
#or use Ornstein-Uhlenbeck(OU) model: trait evolving over tree are pulled toward an optimum, assumed to model process of stabilizing selection
#use similar ZO, sigma^2, and also optimum (theta) and strength of selection(alpha)
#alpha determines how strongly the trait is pulled back to optimum - when a=0, no selection - OU is mechanistically identical to BM process

#perform PGLS for simple regression models that doesn't include ecomorph and those that do 

#PGLS under BM, no ecomorph
pgls.BM1 <- gls(HTotal ~SVL, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#PGLS under BM, w ecomorph
pgls.BM2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")


#PGLS under OU, no ecomorph
pgls.OU1 <- gls(HTotal ~SVL, correlation = corMartins(0,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#PGLS under OU, w, ecomorph
pgls.OU2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corMartins(0,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#use the models to perform AIC operation - see which fit the best
anole.phylo.aic <- AICc(pgls.BM1,pgls.BM2,pgls.OU1,pgls.OU2)
aicw(anole.phylo.aic$AICc)
#including Ecomorph2, under BM is best fit 
#OU models that specify a pull to some global optimum rejected 
#traits evolved randomly, but randomly within each lineage 

#now consider whether ecomorph is significant factor in predicting HTotal-SVL relationship in a phylogenetic context - ANOVA on pgls.BM2
anova(pgls.BM2)
#variable remains significant - Ecomorph2 has significant effect on relationsihp-mega low p value 
#isn't as strong as a factor when we consider phylogeny 

#go back to residuals of HTotal-SVL, but now consider BM evolution of relationship over the tree 
#mutate and redifine anole.log data to include the column for phylogenetically corrected residuals, then plot against ecomorph 
anole.log <- anole.log%>%
  mutate(phylo.res=residuals(pgls.BM2))

p.eco.phylo <- anole.log%>%
  ggplot(aes(x=Ecomorph2, y=phylo.res))+geom_boxplot()+stat_summary(fun=mean, geom="point", size=3)

print(p.eco.phylo)

#now compare phylogenetically corrected residual with uncorrected oens 
#make anole.log tibble longer with respect to the two types of the residuals 
#column for each type of residual - make tibble longer, with a column for both residual values and another identifying to which type of residual the value belongs 
#ggplot faceting tool to include boxplot of both residual types in one plot 
#facets break plot up into grid according to values in data 

anole.log%>%
  dplyr::select(Ecomorph2,res,phylo.res)%>%
  pivot_longer(cols=c("res","phylo.res"))%>%
  print%>%
  ggplot(aes(x=Ecomorph2,y=value)) +geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)+facet_grid(name~.,scales = "free_y")+ylab("residual")
