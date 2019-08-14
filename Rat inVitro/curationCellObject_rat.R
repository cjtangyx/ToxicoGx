#extract from phenodata file
phenoData <- readRDS("rds/phenoData_rat.rds")
#subset necessary columns
curationCell <- subset(phenoData,select=c(cellid))#chose cellid instead of samplename
#add new column
curationCell$tggates.cellid <- curationCell$cellid
#change first column name
names(curationCell)[1] <- "unique.cellid"

#save object into data/
saveRDS(curationCell,file="rds/curationCell.rds")