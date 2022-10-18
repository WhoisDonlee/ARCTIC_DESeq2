library("DESeq2")

counts <- as.matrix(read.csv(
  file = "~/Data/ARCTIC/ARCTIC_CountMatrix.txt",
  sep = "\t",
  row.names = "Geneid"
))
colnames(counts) <- sub("T[I]*_.+$", "", colnames(counts))

counts_clin <- read.csv(
  file = "~/Data/ARCTIC/ARCTIC_Clinical_data_reduced_MP.csv",
  row.names = 1
)
counts_clin <- na.omit(counts_clin)

nrow(counts_clin)

rownames(counts_clin)
colnames(counts)

intersect <- Reduce(intersect, list(colnames(counts), rownames(counts_clin)))
counts <- counts[, intersect]
counts_clin <- counts_clin[intersect, ]

all(ncol(counts[, intersect]) == nrow(counts_clin[intersect, ]))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = counts_clin,
  design = ~Responder_type
)

dds$Responder_type <- factor(dds$Responder_type, levels = c("0", "1", "2"))

dds <- DESeq(dds)

res <- results(dds)
res

resultsNames(dds)

summary(res)

res05 <- results(dds, alpha = 0.05)
summary(res05)

sum(res$padj < 0.1, na.rm = TRUE)
sum(res05$padj < 0.05, na.rm = TRUE)

## MA-plot
plotMA(res, ylim = c(-2, 2))

## Plot counts
plotCounts(dds, gene = which.min(res$padj), intgroup = "Tumor_type")
