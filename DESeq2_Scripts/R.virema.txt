#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library('DESeq2')

Infile = paste(args[1], sep="")
counts <- read.delim(Infile, header=TRUE, row.names=1)
counts <- as.matrix(counts)
condition <- factor(c(rep(args[2], args[4]), rep(args[3], args[4])))
coldata <- data.frame(row.names=colnames(counts), condition)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- DESeq(dds, fitType="local")
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
Outfile = paste(args[1], "_ViReMa_DESeq2-results.csv", sep="")
write.csv(resdata, file=Outfile)

pdf(paste(args[2], "-vs-", args[3], "_PCA_plot_Genes.pdf", sep=""))
plotPCA(vsd, intgroup=c("condition"))
dev.off()




quit("yes")
