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
