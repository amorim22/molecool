## look at p50, when half of the population has arrived

## temperature, wind, precipitation as variables

## shift between species?

## carbon trapping in north-hemisphere, source of resources (food)

## look at variables from where they're are coming from as well as AL, NC and MA.
## see what they'll have to say about the time of these birds arriving in MA

## leave out weather data if including more variables.

## ideas: 10 total

Ruby-throated hummingbird
Great Crested Flycatcher
Yellow-bellied Sapsucker
Eastern Wood-Peewee
Scarlet tanager
Blue-winged warbler
Northern waterthrust
Swainsons Thrush
Broadwing Hawk
Black-billed Cuckoo


database: https://avibase.bsc-eoc.org/species.jsp?lang=EN&avibaseid=DEC3D8D0C74C598B

## using wind:
regress arrival time against this, direction

## sin(200 - wind vector)*(wind speed) in order to normalize data

ex: 200 degrees, subtract wind vector (aka direction), convert to radians

## wind 3,4,5 days before they arrived in Boston

## plotting residuals because some birds arrive later

## linear mixed effects

## is the arrival time changing? Deviation from arrival time.

## run AIC

## try different models + compare them

ACTUAL PROJECT

library(rgbif)
library(tidyverse)
library(MuMIn)
library(rnoaa)
library(data.table)
library(ggmap)
library(usmap)
library(magick)#for examples
library(cowplot)#for examples
library(lme4) #for linear mixed effect models
library(car) #for LME anova testing
library(data.table) #for frollmean function (and others)

species <- c("Sphyrapicus varius","Archilochus colubris","Vermivora cyanoptera","Parkesia noveboracensis","Buteo platypterus","Contopus virens","Coccyzus erythropthalmus","Myiarchus crinitus","Catharus ustulatus","Piranga olivacea")

y <- paste0("1990",",","2019")
m <- paste0("3",",","6")

dat.l <-list()

for(s in species){
  
  ## setting the limit=0 returns no records but can give you the number of observations if you access the meta data
  n.obs <-  occ_data(scientificName = s,year=y,month=m,limit=0,country="US",basisOfRecord = "HUMAN_OBSERVATION",stateProvince="Massachusetts")$meta$count 
  
  print(n.obs)
  
  
  dat.l[[paste0(s)]] <- occ_data(scientificName = s,year=y,month=m,
                                 limit=n.obs,country="US",
                                 basisOfRecord = "HUMAN_OBSERVATION",
                                 stateProvince="Massachusetts")[[2]]
  
  
}


dat <- rbindlist(dat.l,fill=T)

head(dat)

saveRDS(data,"massbird.data.RDS")

dat%>%
  group_by(year,species)%>%
  summarise(count=sum(individualCount,na.rm = T))%>%
  ggplot(aes(x=year,y=count,col=species))+geom_point()

options(noaakey = "xIQmhYFPwNFiyWfXvHvzdLNfxaGypXmE")

sts <- c(
  "GHCND:USW00013894", #Mobible, AL 2k away about 10 days away @200 km/day
  "GHCND:USW00013881", #Charlotte, NC 1000 km away about 6 days away @200 km/day
  "GHCND:USW00014739" #Boston
)

bos <- ncdc_stations(stationid = "GHCND:USW00014739")

print(bos)

library(rgdal)

sta.d <- bind_rows( #bind the rows
  lapply(sts,function(x) ncdc_stations(stationid = x)$data ) #use lapply to run through stations
)%>%
  left_join(usmap_transform(.[,c("longitude","latitude")]))%>% #join transformation of lat/long for projection with usmap
  mutate(name=str_sub(name, -5,-4))%>%#simplify the name column, grab just the state
  mutate(migr.day=c(10,5,0))%>% #so we can look at wind speed 0, 5 or 10 days before arrive in boston
  separate(id,into = c("station.type","id"))%>%#need to cut station type out from station id number
  print()

plot_usmap(
  include = c(.northeast_region,.south_region,.east_north_central)
)+geom_point(data=sta.d,aes(x=longitude.1,y=latitude.1,col=name),size=5)+geom_label(data=sta.d,aes(x=longitude.1,y=latitude.1,col=name,label=name),size=5,nudge_x = 1e6*0.25)+theme(legend.position = "none")

weather.d <- meteo_pull_monitors(sta.d$id,date_min = "2000-01-01")

head(weather.d)

weather.d <- weather.d%>%
  dplyr::mutate(year=as.integer(str_sub(date,1,4)), #add year
                date=as.Date(date))%>%
  group_by(year)%>% #group by year so we can compute julian day
  dplyr::mutate(j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01"))), #add julian day
                date2=date,
                wdir.rad=(200-abs(wdf2-200))*pi/180, #radians so we can use a trig function to compute wind vector, scale degrees first to 180 scale to 2x pi and subtract from 180 (wind comes out of a direction)
                wvec=cos(wdir.rad)*-1*awnd # we want a negative value for positive value for 2x pi
  )%>% #store day in new column
  select(id,year,date2,j.day,wvec)%>% #select the rows we need
  left_join(sta.d%>%select(id,name,migr.day))%>% #add the station id info (ie. name)
  mutate(j.day=j.day+migr.day)#make j.day ahead of BOS according to the migration days away so we can join weather along path

mc<- dat%>%
  group_by(year)%>%
  mutate(date=as.Date(paste0(year,"-",month,"-",day)),
         j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01")))
  )%>%
  group_by(species,year,j.day,date)%>%
  summarise(day.tot=sum(individualCount,na.rm=T))%>%
  group_by(species,year)%>%
  mutate(prop=cumsum(day.tot/sum(day.tot,na.rm = T)))%>%
  filter(year>1999)

mc.pred <- mc%>%
  group_by(year)%>%
  summarize(
    pred=predict(nls(prop~SSlogis(j.day,Asym, xmid, scal)),newdata=data.frame(j.day=min(j.day):max(j.day))),#predict the logistic curve for each species
    j.day=min(j.day):max(j.day),
  )%>%
  left_join(mc%>%select(j.day,date)) ## add date back to tibble

mc.arrive.date <-mc.pred%>%
  group_by(year)%>%
  filter(j.day==j.day[which.min(abs(pred-0.25))])


mc.arr.weath <- mc.arrive.date%>%
  left_join(weather.d)%>%
  left_join(mc%>%select(year,date,j.day))

head(mc.arr.weath)

weather.wk <-weather.d %>% 
  group_by(year,name) %>% 
  mutate(wk.wvec = frollmean(wvec, n=14,align="right")
  )%>%
  select(j.day,date2,name,wk.wvec)

mc.arr.weath2 <- mc.arrive.date%>%
  left_join(weather.wk)

head(mc.arr.weath2)

species.mass <- read_csv("species.mass.txt")

mc.arr.weath.mass <- mc.arr.weath%>%
  left_join(species.mass)%>%
  print()

mc.lmer <- lmer(j.day~wvec+(1|species),mc.arr.weath,na.action = "na.fail")

Anova(mc.lmer) #Anova from the car package

mc.lmer2 <- lmer(j.day~wk.wvec+(1|species),mc.arr.weath2,na.action = "na.fail")

Anova(mc.lmer2) 

library(kableExtra)

mc.arr.aic <- dredge(mc.lmer2,fixed = c("wk.wvec"),)

mc.kb <- kable(mc.arr.aic[1:2,],caption = "Fit values for nested models of the most complicated lme model")

kable_styling(mc.kb)

best.lmer <-  lmer(j.day~wk.wvec+(1|species),mc.arr.weath2,na.action = "na.fail")

Anova(best.lmer)