library("DESeq2")
library("ggplot2")
library("Glimma")

## Import count matrix generated with featureCounts
counts <- as.matrix(read.csv(
  file = "~/Data/ARCTIC/ARCTIC_CountMatrix.txt",
  sep = "\t",
  row.names = "Geneid"
))
## Regex to strip the colname to the sample basename
colnames(counts) <- sub("T[I]*_.+$", "", colnames(counts))

## Import clinical data
counts_clin <- read.csv(
  file = "~/Data/ARCTIC/ARCTIC_Clinical_data_reduced_MP.csv",
  row.names = 1
)
## Remove rows containing NA
counts_clin <- na.omit(counts_clin)

## Count rows in clinical csv and in countmatrix
nrow(counts_clin)
ncol(counts)

## Print clinical and countmatrix sample names
rownames(counts_clin)
colnames(counts)

## Intersect clin and count arrays
intersect <- Reduce(intersect, list(colnames(counts), rownames(counts_clin)))
counts <- counts[, intersect]
counts_clin <- counts_clin[intersect, ]

## Check if clin and count arrays contain the same data
all(ncol(counts[, intersect]) == nrow(counts_clin[intersect, ]))
all(nrow(counts_clin[intersect, ]) == ncol(counts[, intersect]))

factor(dds$Responder_type)

## Create DESeq dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = counts_clin,
  design = ~Responder_type
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

## Add responder types as factors
dds$Responder_type <- factor(dds$Responder_type, levels = c("0", "1", "2"))

dds$Responder_type <- relevel(dds$Responder_type, ref = "0")

dds$Responder_type <- droplevels(dds$Responder_type)

## Run differential expression
dds <- DESeq(dds)
res <- results(dds)
res

results(dds, name = "Responder_type_2_vs_0")
results(dds, name = "Responder_type_1_vs_0")

resultsNames(dds)

## Shrinking Log Fold Change
resultsNames(dds)
reslfc <- lfcShrink(dds,
  contrast = c("Responder_type", "1", "2"),
  type = "ashr"
)
reslfc

## Order results by smalles p value
res_ordered <- res[order(res$pvalue), ]

## Results summary
summary(res)

## Results summary with p value < 0.1
sum(res$padj < 0.1, na.rm = TRUE)

## Results summary with alpha of 0.05
res05 <- results(dds, alpha = 0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm = TRUE)

library(Glimma)

glimmaMDS(dds)

glimmaMA(dds)

## MA plot
plotMA(res, ylim = c(-2, 2))

## MA plot with shrunken LFC
plotMA(resLFC, ylim = c(-2, 2))

## Read counts for the gene with the lowest p value
plotCounts(dds, gene = which.min(res$padj), intgroup = "Responder_type")

d <- plotCounts(dds, gene = which.min(res$padj), intgroup = "Responder_type", returnData = TRUE)
library("ggplot2")
ggplot(d, aes(x = Responder_type, y = count)) +
  geom_point(position = position_jitter(w = 0.1, h = 0)) +
  scale_y_log10(breaks = c(25, 100, 400))

mcols(res)$description

sum(res$padj < 0.1, na.rm = TRUE)
sum(res05$padj < 0.05, na.rm = TRUE)

## MA-plot
plotMA(res, ylim = c(-2, 2))

## Plot counts
plotCounts(dds, gene = which.min(res$padj), intgroup = "Tumor_type")
