## Created 2 different functions to import the counts from protein (CITE_count)or RNA libraries.

read_CITE_count <- function(dir = "./data/raw/CITEseq/CITEseq/CITE_counts/" , name) {
  m <- read.delim(paste0(dir, name, ".txt") )
  m <- m[1:78,]
  colnames(m) <- substring(colnames(m),6)
  return(m)
}

read_RNA_count <- function(dir = "./data/raw/COUNT_TABLES/COUNT_TABLES_" , name) {
  m <- read_delim(paste0(dir, name,"/COUNTS_",name,".tab"),
                  "\t", escape_double = FALSE, trim_ws = TRUE)
  # get umi count table
  m <- m[grepl("^ENSG|^ERCC",m$`#Sample:`),]
  m <- as.data.frame(m)
  rownames(m) <- m$`#Sample:`
  colnames(m) <- substring(colnames(m),1, 14)
  colnames(m) <- substring(colnames(m),6)
  m  <-  m[,-1]

  return(m)
}

## This specific experiment contains data from 10 plates.

## Read protein libraries (here I need CITEseq_platenumber names)
samplenames.citeseq_counts <- c(paste0("CITEseq_",  c(1578, 1579, 1580, 1582, 1584, 1585, 1586,1587, 1588, 1589)))

for (i in samplenames.citeseq_counts) {

  assign(i, read_CITE_count(name = i ))

}

## Read RNA libraries
samplenames.RNA <-  c(1578, 1579, 1580, 1582, 1584, 1585,1586, 1587, 1588, 1589)

for (i in samplenames.RNA) {

  assign(paste0("RNAseq_",i), read_RNA_count(name = i ))

}

## Combine all RNA & replace GeneIDs with gene names
RNA_all <- bind_cols(RNAseq_1578, RNAseq_1579, RNAseq_1580, RNAseq_1582, RNAseq_1584, RNAseq_1585, RNAseq_1586,  RNAseq_1587, RNAseq_1588, RNAseq_1589)

symbols <- make.unique(mapIds(org.Hs.eg.db, keys = sub("\\..*", "", rownames(RNA_all)), keytype = "ENSEMBL", column="SYMBOL"))
symbols[is.na(symbols)] <- "XXXX"
rownames(RNA_all)  <- symbols
RNA_all <- RNA_all[!grepl("^NA",rownames(RNA_all)),] ## remove lots of NA in genelist?!!

## Combine Protein from fixed plats and live plates in 2 objects

PROT_all_fix <- bind_cols(CITEseq_1586, CITEseq_1587, CITEseq_1588, CITEseq_1589)
PROT_all_live <- bind_cols(CITEseq_1578, CITEseq_1579, CITEseq_1580, CITEseq_1582, CITEseq_1584, CITEseq_1585)

PROT_all_fix <- PROT_all_fix[-c(77:78),]
PROT_all_live <- PROT_all_live[-c(54:78),]

row.names(PROT_all_fix) <- gsub('.{9}$', '', row.names(PROT_all_fix))
row.names(PROT_all_live) <- gsub('.{9}$', '', row.names(PROT_all_live))

## Create seurat objects & save

seu.RNA <- CreateSeuratObject(counts = RNA_all, project = "Bcells", min.cells = 5, min.features = 1)
seu.RNA[["percent.ERCC"]] <- PercentageFeatureSet(seu.RNA, pattern = "^ERCC")
seu.RNA[["percent.mt"]] <- PercentageFeatureSet(seu.RNA, pattern = "^MT")
seu.RNA[["percent.rb"]] <- PercentageFeatureSet(seu.RNA, pattern = "^RPL|^RPS")
saveRDS(seu.RNA, file = "output/seu.RNA.rds")

seu.PROT_fix <- CreateSeuratObject(counts = PROT_all_fix, project = "Bcells", min.cells = 0, min.features = 0, assay = "PROT")
seu.PROT_fix[["percent.Ig"]] <- PercentageFeatureSet(seu.PROT_fix, pattern = "^Ig", assay = "PROT")
seu.PROT_fix[["percent.HisH3"]] <- PercentageFeatureSet(seu.PROT_fix, pattern = "^His", assay = "PROT")
saveRDS(seu.PROT_fix, file = "output/seu.PROT_fix.rds")

seu.PROT_live <- CreateSeuratObject(counts = PROT_all_live, project = "Bcells", min.cells = 0, min.features = 0, assay = "PROT")
seu.PROT_live[["percent.Ig"]] <- PercentageFeatureSet(seu.PROT_live, pattern = "^Ig", assay = "PROT")
seu.PROT_live[["percent.HisH3"]] <- PercentageFeatureSet(seu.PROT_live, pattern = "^His", assay = "PROT")
saveRDS(seu.PROT_live, file = "output/seu.PROT_live.rds")
