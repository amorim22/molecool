library(rgbif)
library(tidyverse)
library(MuMIn)
library(rnoaa)
library(data.table)
library(ggmap)
library(usmap)
library(magick)
library(cowplot)
library(lme4) 
library(car) 
library(data.table)

species <- c("Sphyrapicus varius","Archilochus colubris","Vermivora cyanoptera","Parkesia noveboracensis","Buteo platypterus","Contopus virens","Coccyzus erythropthalmus","Myiarchus crinitus","Catharus ustulatus","Piranga olivacea")
y <- paste0("1990",",","2019")
m <- paste0("3",",","6")

dat.l <-list()
for(s in species){
  n.obs <-  occ_data(scientificName=s,year=y,month=m,limit=0,country="US",basisOfRecord = "HUMAN_OBSERVATION",stateProvince="Massachusetts")$meta$count 
  print(n.obs)
  
  dat.l[[paste0(s)]] <- occ_data(scientificName = s,year=y,month=m,
                                 limit=n.obs,country="US",
                                 basisOfRecord = "HUMAN_OBSERVATION",
                                 stateProvince="Massachusetts")[[2]]
}
dat <- rbindlist(dat.l,fill=T)
head(dat)

saveRDS(data,"massbird.data.RDS")
library(tidyverse)


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
  select(id,year,date2,j.day,tmin,tmax,wvec)%>% #select the rows we need
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

species.mass <- read_csv("species.mass.txt")

mc.arr.weath.mass <- mc.arr.weath%>%
  left_join(species.mass)%>%
  print()



