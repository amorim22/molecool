library(tidyverse)
pseed <- read_csv("pseed.fin.amps.csv")
spec(pseed)
pseed.bl <- read_csv("pseed.lengths.csv")
spec(pseed.bl)
speeds <- read_csv("pseed.calibration.csv")
spec(speeds)
pseed2 <- pseed%>%
  left_join(speeds,by=c("speed"="vol"))%>%
  print()
pseed.bl%>%
  print()
pseed2%>%
  select(fish)%>%
  unique()
pseed2 <- pseed2%>%
  left_join(pseed.bl,by="fish")%>%
  print()
pseed2 <- pseed2%>%
  mutate(bl.s=cm.s/bl)%>%
  print()
pseed2%>%
  ggplot(aes(x=bl.s,y=amp.bl))+geom_point()
pseed2%>%
  ggplot(aes(x=bl.s,y=amp.bl))+geom_point(alpha=0.01)
pseed2%>%
    filter(date=="2019-06-17-151149", fin=="L")%>%
    ggplot(aes(x=frame,y=amp.bl))+geom_point()
library(features)
exp1 <- pseed2%>%
  filter(date=="2019-06-17-151149", fin=="L")

f1 <-  features(x = exp1$frame,y=exp1$amp.bl)->f1
fget(f1)
pseed2%>%
  filter(date=="2019-06-17-151149", fin=="L")%>%
  ggplot(aes(x=frame,y=amp.bl))+geom_point()+geom_vline(xintercept = fget(f1)$crit.pts)
f2 <-  features(x = exp1$frame,y=exp1$amp.bl*100)
fget(f2)
f.tib <- fget(f2)[2:3]%>%
  as_tibble()%>%
  filter(curvature<0)%>%
  mutate(peaks=round(crit.pts,0))%>%
  print()
pseed2%>%
  filter(date=="2019-06-17-151149", fin=="L")%>%
  mutate(peak=frame %in% f.tib$peaks)%>%
  ggplot(aes(x=frame,y=amp.bl,col=peak))+geom_point()
pseed2%>%
    summarize(n=length(unique(date)))
find.peaks <- function(x,y,mult=100){
f <- fget(features(x = x, y=y*mult))[2:3]%>%
  as_tibble()%>%
  filter(curvature<0)%>%
  mutate(peaks=round(crit.pts,0))
  return(f$peaks)
}
pseed2%>%
  filter(date%in%unique(date)[1:3])%>%
  group_by(date,fin)%>%
  mutate(peak=frame %in% find.peaks(frame,amp.bl))%>%
  ggplot(aes(x=frame,y=amp.bl,alpha=peak,col=peak))+geom_point()+facet_grid(date~fin)
pseed.max <- pseed2%>%
  group_by(date,fin)%>%
  mutate(peak=frame %in% find.peaks(frame,amp.bl))%>%
  filter(peak==T) #new filter