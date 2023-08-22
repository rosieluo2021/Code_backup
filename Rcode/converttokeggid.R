#BridgeDbR
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BridgeDbR")
library(BridgeDbR)
#download database manually
location = normalizePath('~/Downloads/metabolites_20210109.bridge')
mapper <- loadDatabase(location)

getMatchingSources("C00803")
pubchem<-getSystemCode("PubChem-compound")
vmh<-getSystemCode("VMH metabolite")
kegg<-getSystemCode('KEGG Compound')

all_pubchemid<-unique(read.csv('~/Downloads/all_pubchem.csv',header = T))
source1<-data.frame(rep(pubchem,nrow(all_pubchemid)))
names(source1)<-'source'
names(all_pubchemid)<-'identifier'
input<-cbind(source1,all_pubchemid)
maptokegg<-maps(mapper,input,target = kegg)
write.csv(maptokegg,file='~/drive/metaPD/results/metabolomics/pubchemidtokegg.xlsx')

