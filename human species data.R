library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

f <- list.files("./project6data",pattern=".csv",full.names = T)

f.l <- list()
for(i in f){
  f.i <- read_csv(i)
  if(ncol(f.i)==1) coln <- 1
  if(ncol(f.i)==2) coln <- 2
  m <- strsplit(i,"_")%>%unlist
  print(m)
  sub <- m[2]
  time <- m[3]%>%tolower()
  mass <- as.numeric(gsub(".csv","",m[4]))
  f.i <- f.i[,coln]
  colnames(f.i) <- "Temp"
  f.l[[i]] <- f.i%>%
    mutate(N=1:n(),Temp=as.numeric(Temp),subject=sub,Time=time,mass=mass)
}

dat <- do.call(rbind,f.l)

dat%>%
  group_by(subject,Time,mass)%>%
  na.omit()%>%
  filter(N>0.5*max(N))%>%
  summarize(m.temp=mean(Temp),n=n())%>%
  pivot_wider(names_from=Time,values_from = c(m.temp,n))%>%
  mutate(Tdifference=abs(m.temp_day-m.temp_night))%>%
  ggplot(aes(log(mass),log(Tdifference)))+geom_point()+geom_smooth(method='lm')

print(dat)

dat.log <- dat%>%
  mutate_at(c("Temp", "mass"),log)

dat.log.lm <- lm(Temp~mass,dat.log)

dat.log.allo <- nls(Temp~a*mass^b, start=list(b=0.1, a=0.1),data = dat.log)

summary(dat.log.lm)

summary(dat.log.allo)

dat.log.aic <- AICc(dat.log.lm,dat.log.allo)

dat.log.aicw <- aicw(dat.log.aic$AICc)

print(dat.log.aic)

print(dat.log.aicw)



