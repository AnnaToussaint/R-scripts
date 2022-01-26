if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", version = "3.14", force=TRUE)

BiocManager::install("Rsubread" , version = 3.14)
BiocManager::install("annotables")
BiocManager::install("parathyroidSE")
setwd("~/repos/me/R-scripts/data")
#load reference, build index, import reads, map reads to references 

library(Rsubread)
#single end sequence 
ref <- system.file("extdata", "reference.fa" , package="Rsubread")
buildindex(basename="reference_index", reference = ref)
reads <- system.file("extdata", "reads.txt.gz", package="Rsubread")

#align reads to reference
align.stat <- align(index = "reference_index",
                    readfile1 = reads, 
                    output_file = "alignResults.BAM",
                    phredOffset = 64)


#paired end sequence 
reads1 <- system.file("extdata","reads1.txt.gz" , package="Rsubread")
reads2 <- system.file("extdata", "reads2.txt.gz", package="Rsubread")

align.stat2 <- align(index = "reference_index", readfile1 = reads1, readfile2 = 
                       reads2, output_file = "alignresultsPE.BAM", phredOffset = 
                       64)

#annotation table
ann <- data.frame(
  GeneID=c("gene1","gene1","gene2","gene2"),
  Chr="chr_dummy",
  Start=c(100, 1000, 3000, 5000),
  End=c(500, 1800, 4000, 5500),
  Strand=c("+","+","-","-"),
  stringsAsFactors = FALSE) 
ann

#feature column 
fc_SE <- featureCounts("alignResults.BAM", annot.ext = ann)
fc_SE$counts
fc_PE <- featureCounts("alignresultsPE.BAM", annot.ext = ann, isPairedEnd = TRUE)
fc_PE

#using library data 

library("parathyroidSE")
data("parathyroidGenesSE")

library(DESeq2)

#create a count matrix for contrast
dummy_cm <- cbind(fc_SE$counts.fc_PE$counts)
colnames(dummy_cm) <-c("Sample1","Sample2")
countmatrix <- assay(parathyroidGenesSE)




