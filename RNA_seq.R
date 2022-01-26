if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.14", force=TRUE)

require("DESeq2")

RNA_seq_data_normalization_and_clustering 

setwd("~/repos/me/R-scripts/data")

just.raw.counts = read.delim("Raw_counts_input.txt")
head(just.raw.counts)
dim(just.raw.counts)

just.raw.counts = read.delim("Raw_counts_input.txt" , row.names = 1)
head(just.raw.counts)

dim(just.raw.counts)

meta.data = read.delim(file = "meta_data.txt", row.names = 1)
head(meta.data)
count.data.set <- DESeqDataSetFromMatrix(countData=just.raw.counts, colData = meta.data,
                                         design = ~ condition)

count.data.set.object <- DESeq(count.data.set)

vsd <-vst(count.data.set.object)
#or, we can use 
#rld <-rlog(count.data.set.object)

norm.data = assay(vsd)

head(norm.data)

write.table(norm.data, sep="\t",file="Norm_data_all_genes_NO_counts_cut_off.txt"
            , row.names=TRUE,col.names=NA,quote=FALSE)

sampleDists <- dist(t(norm.data),  method = "euclidean")
reversed_rows_columns = (t(norm.data))
reversed_rows_columns[1:5,1:5]
sampleDists
clusters=hclust(sampleDists)
plot(clusters)
plotPCA(vsd, intgroup=c("condition")) 

require(ggplot2)
plotPCA(vsd, intgroup=c("condition")) +
  scale_colour_hue(breaks = c("E14.5","E14.5","Neonatal","Neonatal","Adult",
                              "Adult","TAC","TAC"))


