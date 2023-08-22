install.packages('maps')
library(maps)
library(ggplot2)
library(tidyverse)

world_map<-map_data('world')
country<-c('Australia','Canada','Denmark','Greece','Israel','Ireland','Malasiya', 'Switzerland','Poland','Norway','Russia','Spain','Sweden','India','Germany','UK','Japan','Italy','USA','China')
countryname<-c('Australia','Canada','Denmark','Greece','Israel','Ireland','Malasiya', 'Switzerland','Poland','Norway','Russia','Spain','Sweden','India','Germany','UK','Japan','Italy','USA','China')
number<-c(1,1,1,1,1,1,1,1,3,2, 2,4,4,5, 6,8,8,8,14,15)
freq<-c('1','1','1','1','1','1','1','1','3','2', '2','4','4','5', '6','8','8','8','14','15')
data<-data_frame(country,countryname,number)
data_world_map<-world_map %>% left_join(data,by=c('region'='country'))

ggplot(data_world_map,aes(x=long,y=lat,group=group))+
  geom_polygon(aes(fill=number),colour='white')+
  scale_x_continuous(breaks = seq(-180,210,45),labels=function(x){paste0(x,".")})+
  scale_y_continuous(breaks = seq(-60,100,30),labels=function(x){paste0(x,".")})+
  scale_fill_gradient(low='lightblue',high = 'steel blue')+
  labs( y='Latitude',x='Longitude')+
  #geom_text(data=data_world_map,aes(x=long-3,y=lat-5,label=number,group=group),col='orange')+
  theme_light()+theme(axis.text=element_text(size=28),axis.title.x=element_text(size=28),axis.title.y=element_text(size=28))


