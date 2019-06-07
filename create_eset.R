# To normalize gene expression data and export as expression set in .RDS file form.

# Install Affymetrix chip (cdf) package to manage data arrays of human hepatocytes. 
### Name of package: Human Genome U133 Plus 2.0 Array, found at http://www.affymetrix.com/support/technical/byproduct.affx?product=hg-u133-plus
install.packages("data/hgu133plus2hsensgcdf_22.0.0.tar.gz", repos = NULL, type = "source")

# Install Affymetrix package from BioConductor.
# Install PharmacoGx package from BioConductor.

# Load libraries installed.
setwd("data")
library("affy")
library("hgu133plus2hsensgcdf")
library("PharmacoGx")

# Group all .CEL files in one variable.
celfn <- list.celfiles("data/CELfiles", full.names = TRUE)

# Perfrom RMA Normalization.
eset <- just.rma(filenames = celfn, verbose=TRUE, cdfname = "hgu133plus2hsensgcdf")

# Save expression set (eset) as .rds file.
saveRDS(eset, file='eset.rds')

# Biobase
# pre-requisites
### info about each experiment
message("Read sample information")
library(Biobase)

path.sc=file.path("data")

# Read .txt file of feature info from metadata.
sampleinfo <- read.table(file.path(path.sc, "TGGATEsfeatureInfo.txt"), sep="\t")
sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA

# Read .csv file of phenoData from metadata.
annot <- read.csv("data/phenowLot.csv", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1)

pData(eset) <- as.data.frame(sampleinfo[match(gsub("[.]CEL[.]gz$", "", rownames(pData(eset))), rownames(sampleinfo)), , drop=FALSE])
colnames(exprs(eset)) <- rownames(pData(eset)) <- gsub("[.]CEL[.]gz$", "", colnames(exprs(eset)))
controls <- rownames(exprs(eset))[grep("AFFX", rownames(exprs(eset)))]
fData(eset) <- fData(eset)[which(!rownames(fData(eset)) %in% controls), , drop=FALSE]
exprs(eset) <- exprs(eset)[which(!rownames(exprs(eset)) %in% controls), , drop=FALSE]
ensemblIds <- sapply(strsplit(rownames(exprs(eset)), "_"), function (x) { return (x[[1]]) }) 
fData(eset) <- data.frame("Probe"=rownames(exprs(eset)), 
                          "EnsemblGeneId"=ensemblIds,
                          "EntrezGeneId"=annot[ensemblIds, "EntrezGene.ID"],
                          "Symbol"=annot[ensemblIds, "gene_name"],
                          "GeneBioType"=annot[ensemblIds, "gene_biotype"],
                          "BEST"=TRUE)
rownames(fData(eset)) <- rownames(exprs(eset))
pData(eset)[,"batchid"] <- NA
annotation(eset) <- "rna"
saveRDS(eset, file="eset_final.rds")
