f <- list.files("./BIOL314020",pattern=".csv",full.names = T)

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
  filter(N>0.5*max(N))
  
dat <- dat%>%
  group_by(subject)%>%
  mutate(
      Max = max(Temp, na.rm = T),
      Min = min(Temp, na.rm = T))%>%
  mutate(Tdifference=(Max-Min))%>%
  filter(Tdifference<5)
dat%>% 
  ggplot(aes(x=mass,y=Tdifference))+geom_point()+geom_smooth(method='lm')

  

 