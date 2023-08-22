library(tidyverse)
library(ggplot2)
#load data
melb<-read.csv('~/drive/metaPD/results/metabolomics/model/modelDecompSCT_new2/164/coremodel/subsystems.csv')
melb<-melb[,-3] #remove unnecessary column
melb_top<-melb[1:20,]

rxns<-melb_top[,1:3]
rxns$class<-rep('rxns',20)
rxns$order<-rep(20:1)
colnames(rxns)<-c('Subsystems','Num','Ratio','Class','order')
mets<-melb_top[,c(1,4,5)]
mets$class<-rep('mets',20)
mets$order<-rep(20:1)
colnames(mets)<-c('Subsystems','Num','Ratio','Class','order')
melb_top_new<-rbind(rxns,mets)


melb_pyramid <- melb_top_new %>% 
  mutate(
    Ratio = case_when(
      Class == "rxns" ~ -Ratio,
      Class == "mets" ~ Ratio
    ),
    order=factor(order)
  )


p=ggplot(melb_pyramid,
       aes(x = Ratio,
           y = order,
           fill = Class)) +
  geom_col()+scale_y_discrete(labels = rev(melb_pyramid$Subsystems[1:20]))+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(x='Reaction Ratio vs CoreMetabolites Ratio',y='Subsystems',title='Top 20 subsystems in the model')

datalabel<-paste0(melb_top_new$Ratio,'(',melb_top_new$Num,')')
p+geom_text(aes(label = datalabel,hjust=-0.01),size=4)+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=12),axis.text.x=element_blank())+
  scale_fill_discrete(name = "Class", labels = c('CoreMatabolite','Reaction'))+xlim(-0.9,0.2)
