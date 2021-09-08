library(ggplot2)
library(tidyverse)
setwd("~/Desktop/orgbiorepo/molecool/Project 1 JB")

## the <- command names the right as the left 
dat <- read.csv("scales.csv")

## prompt size of data set with:
dim(dat)

## reveal first lines of data with:
head(dat)

#with "$"
class(dat$N)
class(dat$quadrant)
class(dat$species)
class(dat$specimen)
#with index, produces the same result
class(dat[,1])
class(dat[,2])
class(dat[,3])
class(dat[,4])

## showing you cant do math based commands on non-numerical text:
mean(dat$N)
mean(dat$quadrant)

## hack for showing the class of all of the data at once:
sapply(dat,class)

##  ?  What species: Inspect levels, but make a factor first 
dat$species <- as.factor(dat$species)
species <- levels(dat$species)
species

# How many species
length(species)

## == asks does left = right
dat$species==species[1]

##   ?  nesting logical values within the “species” to return only “A. rupestris"
dat$species[dat$species==species[1]]

##shows the number of observations/punctures for each species
A.rup<-length(dat$species[dat$species==species[1]])
L.gib<-length(dat$species[dat$species==species[2]])
L.mac<-length(dat$species[dat$species==species[3]])
M.sal<-length(dat$species[dat$species==species[4]])
M.sax<-length(dat$species[dat$species==species[5]])
P.fla<-length(dat$species[dat$species==species[6]])
#combine the results with species
species.obs <- data.frame(sp=species,n=c(A.rup,L.gib,L.mac,M.sal,M.sax,P.fla))
species.obs



##USING TIDYVERSE : summarizing data made easy peasy

## the pipe %>% passes results from left to right
dat %>%
  group_by(species) %>%
  summarise(n = n())

## save summary to variable "species.n"
species.n<- dat %>%
  group_by(species) %>%
  summarise(n = n())
species.n

## How many specimens for each species
dat %>% 
  count(species,specimen) %>%
  print() %>%
  count(species,name = "n.specimens")



##LOOPS

## for loops - for “i” in 1 to 10, print each “i”
for(i in 1:10) print(i)

##make a box plot!
for(i in species){
  p <- dat %>%
    filter(species==i)%>%
    ggplot()+geom_boxplot(aes(x=quadrant,y=N))+ggtitle(i)
  print(p)
}

## save PDF file with all figures
pdf("species.quadrant.pdf")
for(i in species){
  p <- dat %>%
    filter(species==i)%>%
    ggplot()+geom_boxplot(aes(x=quadrant,y=N))+ggtitle(i)
  print(p)
}
dev.off()

list.files(pattern=".pdf")


