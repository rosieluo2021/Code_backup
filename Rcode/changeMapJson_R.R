#install.packages('rjson')
library(rjson)
map<-fromJSON(file = '~/drive/metaPD/results/metabolomics/model/modelDecompSCT_new2/164/coremodel/formap/mapfile/new_map.json')
print(map)

for (i in 1:length(map[[2]]$reactions)) {
  map[[2]]$reactions[[i]]$name=map[[2]]$reactions[[i]]$bigg_id
}

# if need to change the mets ID
for (i in 1:length(map[[2]]$nodes)) {
  fields=map[[2]]$nodes[[i]]
  if (any(names(fields) == "bigg_id")){ 
  name=substring(map[[2]]$nodes[[i]]$bigg_id,1,nchar(map[[2]]$nodes[[i]]$bigg_id)-4)
  #name=paste0(name,'[c]')
  map[[2]]$nodes[[i]]$bigg_id=name
  }
}



datajson<- toJSON(map)
cat(datajson,file = "~/drive/metaPD/results/metabolomics/model/modelDecompSCT_new2/164/coremodel/formap/mapfile/new_map_change.json")
