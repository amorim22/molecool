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

#[CPK] This (and much more) is an example. Asked y'all not to include theses [- 1 point]
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
pseed.max%>%
  ggplot(aes(x=bl.s,y=amp.bl))+geom_point()+geom_smooth(method="lm")
amp.aov <-  aov(amp.bl~bl.s,pseed.max)
summary(amp.aov)
pseed.max %>%
  group_by(fish, bl.s) %>%
  summarize(mean.max=mean(amp.bl)) %>%
  ggplot(aes(x=bl.s,y=mean.max,col=fish))+geom_point()+geom_smooth(method="lm")
pseed2
pseed2 <- pseed2 %>%
  group_by(date,frame) %>%
  mutate(amp.sum=sum(amp.bl))
pseed2 %>%
  filter(fin=="R")
pseed.wide <- pseed2 %>%
  select(-amp)%>%
  pivot_wider(names_from = fin,values_from = amp.bl) %>%
  mutate(amp.sum=L+R)%>%
  print() 

#part 2 folks!
#[CPK] Keep library calls at top, the convention
library(tidyverse)
pseed.wide%>%
  filter(date=="2019-06-17-151149")%>%
  ggplot(aes(x=frame,y=amp.sum))+geom_point()

exp1 <- pseed.wide%>%
  filter(date=="2019-06-17-151149")

f1 <-  features(x = exp1$frame,y=exp1$amp.sum)->f1
fget(f1)

pseed.wide%>%
  filter(date=="2019-06-17-151149")%>%
  ggplot(aes(x=frame,y=amp.sum))+geom_point()+geom_vline(xintercept = fget(f1)$crit.pts)

f2 <-  features(x = exp1$frame,y=exp1$amp.sum*100)
fget(f2)

f.tib <- fget(f2)[2:3]%>%
  as_tibble()%>%
  filter(curvature<0)%>%
  mutate(peaks=round(crit.pts,0))%>%
  print()

pseed.wide%>%
  filter(date=="2019-06-17-151149")%>%
  mutate(peak=frame %in% f.tib$peaks)%>%
  ggplot(aes(x=frame,y=amp.sum,col=peak))+geom_point()

find.peaks <- function(x,y,mult=100){ 
  f <- fget(features(x = x,y=y*mult))[2:3]%>%
    as_tibble()%>% 
    filter(curvature<0)%>% 
    mutate(peaks=round(crit.pts,0))
  return(f$peaks) 
}

pseed.wide <- pseed.wide%>%
  group_by(date)%>%
  mutate(peak=frame %in% find.peaks(frame, amp.sum))

#[CPK] same operation as above
pseed.sum.max <- pseed.wide%>%
  group_by(date)%>%
  mutate(peak=frame %in% find.peaks(frame,amp.sum))%>%
  filter(peak==T) #all above was finding mean#[CPK] Nope, just the max

pseed.sum.max%>%
  ggplot(aes(x=bl.s,y=amp.sum))+geom_point()+geom_smooth(method="lm")

#[CPK] heres the mean
pseed.sum.max %>% 
  group_by(fish, bl.s) %>%
  summarize(amp.sum.mean=mean(amp.sum)) %>%
  ggplot(aes(x=bl.s,y=amp.sum.mean,col=fish))+geom_point()+geom_smooth(method="lm")

#[CPK] Repeating agina
pseed.sum.max <- pseed.sum.max %>% #this section adds the amp.sum.mean to pseed.sum.max
  group_by(fish, bl.s) %>%
  mutate(amp.sum.mean=mean(amp.sum))


amp.sum <- pseed.sum.max%>% pull(amp.sum)  #this is defining the amp.sum column in the tibble as a vector

#CPK: why have an agrument name amp.sum. Will this be the data in all cases. For me this causes a problem below. Commenting out and reworking.
#standard_error <- function(amp.sum) sd(amp.sum) / sqrt(length(amp.sum.mean)) #defines standard error function]

standard_error <- function(x) sd(x) / sqrt(length(x))

pseed.sum.max<- pseed.sum.max %>% #makes standard error column                      
  group_by(fish,bl.s)%>%  
  mutate(amp.sum.se = standard_error(amp.sum))
pseed.sum.max  


#CPK: Wouldn't a linear line through the points been better?
ggplot(pseed.sum.max, aes(x=bl.s, y=amp.sum.mean, colour=fish)) + #creates plot
  geom_errorbar(aes(ymin=amp.sum.mean-amp.sum.se, ymax=amp.sum.mean+amp.sum.se), width=.1) +
  geom_line() +
  geom_point()
#now the new data!
pseed.met.rate <- read_csv("pseed.met.rate.csv")

pseed.sum.max <- pseed.sum.max%>% #join tibbles
  left_join(pseed.met.rate,by=c("date"="date"))
print(pseed.sum.max)

#CPK: for me, bl.s.x.y doesn't exist. Probably left_join twice
pseed.sum.max <- pseed.sum.max%>% #neaten data
  select(-m.s.x,-cm.s.x,-L,-R,-fish.y,-m.s.x,-cm.s.y,-bl.s.x.y)



print(pseed.sum.max)

ggplot(pseed.sum.max, aes(x=amp.sum.mean, y=met.rate, colour=fish.x))+geom_point()

#CPK:well done! keep examples out of project code. stick to what you need. [9/10]
