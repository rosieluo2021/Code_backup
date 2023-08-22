library(tidyverse)
library(stringr)
library(ggplot2)
allsubs<- read.csv('~/Downloads/subs.csv',header = TRUE)
#
removedVar=c('Transport, nuclear','Transport, endoplasmic reticular','Transport, mitochondrial','Transport, extracellular','Exchange/demand reaction')
allsubs <- allsubs[!(allsubs$Subsystems %in% removedVar),]
allsubs=allsubs[1:20,]

## 在每个分组数据的后面插入几行缺失值
#empty_bar <- 3
#to_add <- data.frame(matrix(NA, empty_bar, ncol(allsubs)))
#colnames(to_add) <- colnames(allsubs)
#allsubs<-rbind(allsubs,to_add)
# 获取每个样本的名称在y轴的位置和倾斜角度
label_data <- allsubs
label_data$id <- seq(1, nrow(label_data))
number_of_bar <- nrow(label_data) # 计算条的数量
## 每个条上标签的轴坐标的倾斜角度
angle <- 90 - 360 * (label_data$id-0.2) /number_of_bar 
label_data$hjust <- ifelse( angle < -90, 1, 0) # 调整标签的对其方式
label_data$angle <- ifelse(angle < -90, angle+180, angle) ## 标签倾斜角度


## 可视化分组圆环条形图
p1 <- ggplot(label_data)+
  ## 添加条形图
  geom_bar(aes(x=as.factor(id), y=Num,fill=Subsystems),stat="identity",
           alpha=0.8)  +
  ylim(-200,200) + ## 设置y轴坐标表的取值范围,可流出更大的圆心空白
  ## 设置使用的主题并使用极坐标系可视化条形图
  theme_minimal() +
  theme(legend.position = "none", # 不要图例
        axis.text = element_blank(),# 不要x轴的标签
        axis.title = element_blank(), # 不要坐标系的名称
        panel.grid = element_blank(), # 不要网格线
        plot.margin = unit(rep(-1,4), "cm"))+ ## 整个图与周围的边距
  coord_polar() + ## 极坐标系
  ## 为条形图添加文本
  geom_text(data=label_data, 
            aes(x=id, y=Num+5, label=Subsystems,hjust=hjust),
            color="black",fontface="bold",alpha=0.8, size=4, 
            angle= label_data$angle, inherit.aes = FALSE)
p1
  

library(ggpubr)
ggarrange(p1, p2, ncol = 2, nrow = 1,widths = c(0.9, 1)) 



##

library(RColorBrewer)
library(ggplot2)
allsubs<- read.csv('~/Downloads/subs.csv',header = TRUE)
removedVar=c('Transport, nuclear','Transport, endoplasmic reticular','Transport, mitochondrial','Transport, extracellular','Exchange/demand reaction')
allsubs <- allsubs[!(allsubs$Subsystems %in% removedVar),]
topsubs<-allsubs[1:20,]

barplot(rev(topsubs$Proportion),horiz = T,xlim=c(-5,8),axes=F,col='Blue')
text(seq(from=0.7,length.out=135,by=1.2),x=-2,label=rev(topsubs$Subsystems))
axis(3,c(0,0.25,0.5,0.75,1))
