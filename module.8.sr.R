library(tidyverse)
library(dbplyr)
library(ggplot2)
library(MuMIn)
#compiling metadata
k<- list.files("~/Documents/Documents/molecool/Project 8 data",pattern=".csv", full.names = T)
print(k)
k.l <- list()
for(i in k){
  met.dat <- unlist(strsplit(i,"_"))
  group <- met.dat[1]
  subject <- met.dat[2]
  angle <- as.numeric(met.dat[3])
  experiment <- gsub(".csv","",met.dat[4])
  k.l[[i]] <- read_delim(i,delim = " ", col_names = c("Reading","Force","Unit"), id="Experiment") %>%
    mutate(group=group,subject=subject,angle=angle,experiment=experiment)
}

dat<-do.call(rbind,k.l)

if(TRUE){
  d.max1 <- dat%>%
    group_by(subject,experiment,angle)%>%
    mutate(max.force=max(abs(Force),na.rm = TRUE),n=n())
  
  ## computing the max among all of the "maxes"
  d.max2 <- d.max1%>%
    group_by(subject,experiment)%>%
    mutate(max.f2=max(abs(Force)))
  
  ## join the data
  d.max <- d.max2%>%
    left_join(d.max2)%>%
    mutate(normF=max.force/max.f2)
}

d.max%>%
  ggplot(aes(angle,normF,col=experiment))+geom_point()

good.data <- d.max%>%
  group_by(subject)%>%
  summarise(n=n)%>%
  mutate(good=n>=16)%>%
  select(subject,good)

d.max <- d.max%>%
  left_join(good.data)%>%
filter(good== TRUE)

#########################################


AICs <- d.max%>% 
  group_by(subject,experiment)%>% 
  summarize(m2=AICc(lm(normF~poly(angle,2))), 
            m3=AICc(lm(normF~poly(angle,3))), 
            m4=AICc(lm(normF~poly(angle,4))))%>% 
  pivot_longer(m2:m4,names_to="model",values_to="AICc")%>% 
  print()

x.pred <- seq(45,157.5,length.out = 1000)

fits <- d.max%>%
  group_by(subject,experiment)%>%
  summarize(
    m2=predict(lm(normF~poly(angle,2)),newdata=data.frame(angle=x.pred)), #second order
    m3=predict(lm(normF~poly(angle,3)),newdata=data.frame(angle=x.pred)), #third order
    m4=predict(lm(normF~poly(angle,4)),newdata=data.frame(angle=x.pred)) #fourth order
  )%>%
  pivot_longer(m2:m4,names_to="model")%>%
  group_by(subject,experiment,model)%>%
  summarize(theta_max=x.pred[which.max(value)])%>%
  print()

best.models <- fits%>%
  left_join(AICs)%>%
  group_by(subject,experiment)%>%
  mutate(best=AICc==min(AICc))%>%
  filter(best==TRUE)%>%
  select(-best)%>%
  print()

anova(lm(theta_max~experiment,best.models))

best.models%>%
  pivot_wider(id_cols=subject,names_from = experiment,values_from=theta_max)%>%
  mutate(shift=(fatigue-control))%>%
  ungroup()%>%
  summarise(mean.shift=mean(shift, na.rm=TRUE),se.shift=sd(shift, na.rm=TRUE)/sqrt(length(shift)))

