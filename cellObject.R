# extract from phenodata file
phenodata <- read.csv("data/phenowLot.csv",sep=",",header=TRUE, stringsAsFactors=TRUE)

# subset necessary columns
cellObject <- subset(phenodata, select=c(X, ORGAN_ID, MATERIAL_ID, SPECIES, TEST_TYPE, CELL_LOT_TYPE))

# assign row names as barcode 
rownames(cellObject)<-X

# rename column names
names(cellObject)<-c("cellid","tissueid","materialid", "species","testType","batchid")

# save 
save(cellObject, file="cellObject.csv")
