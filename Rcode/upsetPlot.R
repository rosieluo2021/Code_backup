#multiple barchart
library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

##   test multiple group   ##

#test<-read.csv(file = '~/test.csv',header = TRUE)
#test$class[1:17]<-'lipid'
#test$class[11:13]<-'amino'
#test$class[13:17]<-'organic oxgen'

#reorganized data
#testnew<-test %>% gather(trends,realfrequency,incresed_frequency:decresed_frequency)
#testnew<-testnew[order(testnew$class),]

## real data ##

mappedhighfreqmets<-read.csv('~/drive/metaPD/results/metabolomics/data/v1/diagnosis/figure/mappedHighFreqMets_v2.csv')
mappedhighfreqmets<-mappedhighfreqmets[,c(1,4,19,20)]
#reorganized data
new<-mappedhighfreqmets %>% gather(trends,realfrequency,incresed_Realfrequency:decresed_Realfrequency)
new<-new[order(new$superclass),]

q1<-ggplot(data= filter(new, superclass =='Lipids and lipid-like molecules'), mapping=aes(x = all, y = realfrequency,fill=trends))+
  geom_bar(stat="identity",position=position_dodge(0.8))+scale_x_discrete(position = 'bottom')+labs(
    x = "Lipids",
    y = "Frequency" )+theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_blank(), 
      plot.title=element_text(size=15,face="bold"),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10),
      axis.title.x=element_text(vjust=1,size=15),
      axis.title.y=element_text(vjust=1,size=15),
      legend.title = element_text(size=10),
      legend.text = element_text(size=10))+scale_y_continuous(expand = c(0,0),breaks=seq(0,20,1))
q1



q2<-ggplot(data= filter(new, superclass =='Organic acids and derivatives'), mapping=aes(x = all, y = realfrequency,fill=trends))+
  geom_bar(stat="identity",position=position_dodge(0.8))+scale_x_discrete(position = 'bottom')+labs(
    x = "Organic acid",
    y = "Frequency")+theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_blank(), 
      plot.title=element_text(size=15,face="bold"),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10),
      axis.title.x=element_text(vjust=1,size=15),
      axis.title.y=element_text(vjust=1,size=15),
      legend.title = element_text(size=10),
      legend.text = element_text(size=10))+ scale_y_continuous(expand = c(0,0),breaks=seq(0,20,1))

q2


q3<-ggplot(data= filter(new, superclass =='Organic oxygen compounds'), mapping=aes(x = all, y = realfrequency,fill=trends))+
  geom_bar(stat="identity",position=position_dodge(0.8))+scale_x_discrete(position = 'bottom')+labs(
    x = "Organic oxygen compounds",
    y = "Frequency")+theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_blank(), 
      plot.title=element_text(size=15,face="bold"),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10),
      axis.title.x=element_text(vjust=1,size=15),
      axis.title.y=element_text(vjust=1,size=15),
      legend.title = element_text(size=10),
      legend.text = element_text(size=10))+ scale_y_continuous(expand = c(0,0),breaks=seq(0,20,1))


q4<-ggplot(data= filter(new, superclass =='Organoheterocyclic compounds'), mapping=aes(x = all, y = realfrequency,fill=trends))+
  geom_bar(stat="identity",position=position_dodge(0.8))+scale_x_discrete(position = 'bottom')+labs(
    x = "Organoheterocyclic compounds",
    y = "Frequency")+theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_blank(), 
      plot.title=element_text(size=15,face="bold"),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10),
      axis.title.x=element_text(vjust=1,size=15),
      axis.title.y=element_text(vjust=1,size=15),
      legend.title = element_text(size=10),
      legend.text = element_text(size=10))+ scale_y_continuous(expand = c(0,0),breaks=seq(0,20,1))


q5<-ggarrange(q1, q2, ncol = 1, nrow = 2) 
q6<-ggarrange(q3, q4, ncol = 2, nrow = 1) 
par(oma=c(1,1,1,1), mar=c(2,2,2,2))
ggarrange(q5, q6, ncol = 1, nrow = 2) 

################  upset plot #####################
#
#
#                 upset
#
#
##################################################

library(ComplexUpset)
specimens<-read.csv('~/drive/metaPD/papers/firstpaper/graphics/upset_specimen.csv')
specimens<-specimens[1:84,1:8] #only diagnosis
label<-c('Plasma','Serum','CSF','Urinary','Fecal','Brain','Other')
color<-c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2")
q1<-upset(specimens,c('Plasma','Serum','CSF','Urinary','Fecal','Brain','Other'),group_by='sets',
          queries=list(
            upset_query(group=label[1], color=color[1]),
            upset_query(group=label[2], color=color[2]),
            upset_query(group=label[3], color=color[3]),
            upset_query(group=label[4], color=color[4]),
            upset_query(group=label[5], color=color[5]),
            upset_query(group=label[6], color=color[6]),
            upset_query(group=label[7], color=color[7]),
            upset_query(set=label[1], fill=color[1]),
            upset_query(set=label[2], fill=color[2]),
            upset_query(set=label[3], fill=color[3]),
            upset_query(set=label[4], fill=color[4]),
            upset_query(set=label[5], fill=color[5]),
            upset_query(set=label[6], fill=color[6]),
            upset_query(set=label[7], fill=color[7])
          ),
          sort_intersections='descending',width_ratio=0.1,name='Number of study of each specimen',
          themes =upset_modify_themes(list(
            'Intersection size'=theme(
              axis.text = element_text(size=18),
              axis.title = element_text(size=18)
            )
          )),set_sizes=(upset_set_size()
            + theme(axis.text.x=element_text(size=10,angle=90))) 
)+theme(panel.grid.major = element_blank(),axis.text=element_text(size=18),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        legend.text=element_text(size=20),legend.title=element_text(colour="blue",size=20))
q1

allspecimens<-read.csv('~/drive/metaPD/papers/firstpaper/graphics/upset_allspecimen.csv')
allspecimens<-allspecimens[1:87,1:8] #included prognosis
label<-c('Plasma','Serum','CSF','Urinary','Fecal','Brain','Other')
color<-c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2")
q1<-upset(allspecimens,c('Plasma','Serum','CSF','Urinary','Fecal','Brain','Other'),group_by='sets',
          queries=list(
            upset_query(group=label[1], color=color[1]),
            upset_query(group=label[2], color=color[2]),
            upset_query(group=label[3], color=color[3]),
            upset_query(group=label[4], color=color[4]),
            upset_query(group=label[5], color=color[5]),
            upset_query(group=label[6], color=color[6]),
            upset_query(group=label[7], color=color[7]),
            upset_query(set=label[1], fill=color[1]),
            upset_query(set=label[2], fill=color[2]),
            upset_query(set=label[3], fill=color[3]),
            upset_query(set=label[4], fill=color[4]),
            upset_query(set=label[5], fill=color[5]),
            upset_query(set=label[6], fill=color[6]),
            upset_query(set=label[7], fill=color[7])
          ),
          sort_intersections='descending',width_ratio=0.1,name='Number of study of each specimen',
          themes =upset_modify_themes(list(
            'Intersection size'=theme(
              axis.text = element_text(size=18),
              axis.title = element_text(size=18),
            ))),set_sizes=(upset_set_size() + theme(axis.text.x=element_text(size=15,angle=90))) 
)+theme(panel.grid.major = element_blank(),axis.text=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))
q1

#themes=upset_default_themes(text=element_text(color='red'))


###################   platform


#########################################
platform<-read.csv('~/drive/metaPD/papers/firstpaper/graphics/upset_platform.csv')
platform<-platform[1:84,1:11]
label1<-c('LC-MS','LC-tandem MS','GC-MS','GC-tandem MS','NMR','MS','HPLC-ED','LCECA','MRS','other')
colnames(platform)<-c('study','LC-MS','LC-tandem MS','GC-MS','GC-tandem MS','NMR','MS','HPLC-ED','LCECA','MRS','other')
color1<-c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2","#7F7F7F","#BCBD22","#2297E6")
q2<-upset(platform, c('LC-MS','LC-tandem MS','GC-MS','GC-tandem MS','NMR','MS','HPLC-ED','LCECA','MRS','other'),group_by='sets',
          queries=list(
            upset_query(group=label1[1], color="#1F77B4"),
            upset_query(group=label1[2], color="#FF7F0E"),
            upset_query(group=label1[3], color="#2CA02C"),
            upset_query(group=label1[4], color="#D62728"),
            upset_query(group=label1[5], color="#9467BD"),
            upset_query(group=label1[6], color="#8C564B"),
            upset_query(group=label1[7], color="#E377C2"),
            upset_query(group=label1[8], color="#7F7F7F"),
            upset_query(group=label1[9], color="#BCBD22"),
            upset_query(group=label1[10], color="#2297E6"),
            upset_query(set=label1[1], fill=color1[1]),
            upset_query(set=label1[2], fill=color1[2]),
            upset_query(set=label1[3], fill=color1[3]),
            upset_query(set=label1[4], fill=color1[4]),
            upset_query(set=label1[5], fill=color1[5]),
            upset_query(set=label1[6], fill=color1[6]),
            upset_query(set=label1[7], fill=color1[7]),
            upset_query(set=label1[8], fill=color1[8]),
            upset_query(set=label1[9], fill=color1[9]),
            upset_query(set=label1[10], fill=color1[10])
          ), 
          sort_intersections='descending',width_ratio=0.1,name='Number of study of each specimen',
          themes =upset_modify_themes(list(
            'Intersection size'=theme(
              axis.text = element_text(size=18),
              axis.title = element_text(size=18),
            ))),set_sizes=(upset_set_size() + theme(axis.text.x=element_text(size=15,angle=90))) 
)+theme(panel.grid.major = element_blank(),axis.text=element_text(size=18),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))
q2

allplatform<-read.csv('~/drive/metaPD/papers/firstpaper/graphics/upset_allplatform.csv')
allplatform<-allplatform[1:87,1:11]
label1<-c('LC-MS','LC-tandem MS','GC-MS','GC-tandem MS','NMR','MS','HPLC-ED','LCECA','MRS','other')
colnames(allplatform)<-c('study','LC-MS','LC-tandem MS','GC-MS','GC-tandem MS','NMR','MS','HPLC-ED','LCECA','MRS','other')
color1<-c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2","#7F7F7F","#BCBD22","#2297E6")
q2<-upset(allplatform, c('LC-MS','LC-tandem MS','GC-MS','GC-tandem MS','NMR','MS','HPLC-ED','LCECA','MRS','other'),group_by='sets',
          queries=list(
            upset_query(group=label1[1], color="#1F77B4"),
            upset_query(group=label1[2], color="#FF7F0E"),
            upset_query(group=label1[3], color="#2CA02C"),
            upset_query(group=label1[4], color="#D62728"),
            upset_query(group=label1[5], color="#9467BD"),
            upset_query(group=label1[6], color="#8C564B"),
            upset_query(group=label1[7], color="#E377C2"),
            upset_query(group=label1[8], color="#7F7F7F"),
            upset_query(group=label1[9], color="#BCBD22"),
            upset_query(group=label1[10], color="#2297E6"),
            upset_query(set=label1[1], fill=color1[1]),
            upset_query(set=label1[2], fill=color1[2]),
            upset_query(set=label1[3], fill=color1[3]),
            upset_query(set=label1[4], fill=color1[4]),
            upset_query(set=label1[5], fill=color1[5]),
            upset_query(set=label1[6], fill=color1[6]),
            upset_query(set=label1[7], fill=color1[7]),
            upset_query(set=label1[8], fill=color1[8]),
            upset_query(set=label1[9], fill=color1[9]),
            upset_query(set=label1[10], fill=color1[10])
          ), 
          sort_intersections='descending',width_ratio=0.1,name='Number of study of each platform',
          themes =upset_modify_themes(list(
            'Intersection size'=theme(
              axis.text = element_text(size=18),
              axis.title = element_text(size=18),
            ))),set_sizes=(upset_set_size() + theme(axis.text.x=element_text(size=15,angle=90))) 
)+theme(panel.grid.major = element_blank(),axis.text=element_text(size=15),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))
q2



library(ggpubr)
ggarrange(q1, q2, ncol = 2, nrow = 1,widths = c(0.9, 1)) 
