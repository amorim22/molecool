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

library(tidyverse)
dat <- readRDS("massbird.data.RDS")

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

mc <- dat%>%
  group_by(year)%>%
  mutate(date=as.Date(paste0(year,"-",month,"-",day)),
         j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01")))
  )%>%
  group_by(species,year,j.day,date)%>%
  summarise(day.tot=sum(individualCount,na.rm=T))%>%
  group_by(species,year)%>%
  mutate(prop=cumsum(day.tot/sum(day.tot,na.rm = T)))%>%
  filter(year>1999)

mc%>%
  ggplot(aes(j.day,prop))+geom_point()+facet_wrap(year~.)

mc.pred <- mc%>%
  group_by(year)%>%
  summarize(
    pred=predict(nls(prop~SSlogis(j.day,Asym, xmid, scal)),newdata=data.frame(j.day=min(j.day):max(j.day))),#predict the logistic curve for each species
    j.day=min(j.day):max(j.day),
  )%>%
  left_join(mc%>%select(j.day,date)) ## add date back to tibble

mc%>%
  ggplot(aes(j.day,prop))+geom_point(aes=0.3)+geom_line(data=mc.pred,aes(x=j.day,y=pred),col="blue",size=2)+facet_wrap(year~.)

mc.arrive.date <-mc.pred%>%
  group_by(year)%>%
  filter(j.day==j.day[which.min(abs(pred-0.25))])

mc.arrive.date%>%
  ggplot(aes(year,j.day))+geom_point()

weather.d <- weather.d%>%
  dplyr::mutate(year=as.integer(str_sub(date,1,4)), #add year
                date=as.Date(date))%>%
  group_by(year)%>% #group by year so we can compute julian day
  dplyr::mutate(j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01"))), #add julian day
                date2=date,
                wdir.rad=(200-abs(wdf2-180))*pi/180, #radians so we can use a trig function to compute wind vector, scale degrees first to 180 scale to 2x pi and subtract from 180 (wind comes out of a direction)
                wvec=cos(wdir.rad)*-1*awnd # we want a negative value for positive value for 2x pi
  )%>% #store day in new column
  select(id,year,date2,j.day,tmin,tmax,wvec)%>% #select the rows we need
  left_join(sta.d%>%select(id,name,migr.day))%>% #add the station id info (ie. name)
  mutate(j.day=j.day+migr.day)

mc.arr.weath <- mc.arrive.date%>%
  left_join(weather.d)%>%
  left_join(mc%>%select(year,date,j.day))

head(mc.arr.weath)

weather.wk <-weather.d %>% 
  group_by(year,name) %>% 
  mutate(wk.tmin = frollmean(tmin, n=14,align="right"),
         wk.tmax = frollmean(tmax, n=14,align="right"),
         wk.wvec = frollmean(wvec, n=14,align="right")
  )%>%
  select(j.day,date2,name,wk.tmin,wk.tmax,wk.wvec)

mc.arr.weath2 <- mc.arrive.date%>%
  left_join(weather.wk)

head(mc.arr.weath2)

mc.arr.weath2%>%
  group_by(year,species)%>%
  ggplot(aes(x=wk.wvec,y=date2, col=species))+geom_point()


