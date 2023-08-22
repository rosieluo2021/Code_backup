install.packages('tidyjson')
library(classyfireR)
library(magrittr)
remotes::install_github('aberHRML/classyfireR')
help(classyfireR)
inchi_key<-read.csv(file="C:/Users/xxerw/Desktop/inchi_key.csv",header=F)

class(inchi_key)
inchi<-as.character(inchi_key$V1)

class(inchi)
Classification_List <- purrr::map(inchi, get_classification)
Classification_List[[1]]@classification[["Classification"]]

result<-as.data.frame(matrix(NA,1,10))
colnames(result)<-c("kingdom","superclass","class","subclass","level 5","level 6","level 7","level 8","level 9","level 10")
library(dplyr)

for (i in 1:1100){
  if (is.null(Classification_List[[i]])==FALSE){
    result_1<-as.data.frame(t(Classification_List[[i]]@classification[["Classification"]]))
    colnames(result_1)<-Classification_List[[i]]@classification[["Level"]]
    result<-dplyr::bind_rows(result,result_1)
    }else{
      result_2<-as.data.frame(matrix(NA,1,10))
      result<-dplyr::bind_rows(result,result_2)
      }
  i=i+1
}
result<-result[,-c(11:20)]

des<-as.data.frame(matrix(NA,1,1))
for (i in 1:1100){
  if (is.null(Classification_List[[i]])==FALSE){
  des_1<-as.data.frame(t(Classification_List[[i]]@direct_parent[["description"]]))
  des<-dplyr::bind_rows(des,des_1)
  }else{
    des_2<-as.data.frame(matrix(NA,1,10))
    des<-dplyr::bind_rows(des,des_2)
  }
  i=i+1
}

classy_result<-dplyr::bind_cols(result,des)
classyfire_result<-classy_result[-1,-c(12:20)]
classyfire_result[,12]<-inchi_key$V1
write.csv(classyfire_result,file = "E:/PhD/1.metaPD/classyfire/classyfire_result.csv")

