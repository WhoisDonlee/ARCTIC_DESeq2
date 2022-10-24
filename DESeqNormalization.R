library("DESeq2")

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

## Create DESeq dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = counts_clin,
  design = ~Responder_type
)

## Add responder types as factors
dds$Responder_type <- factor(dds$Responder_type, levels = c("0", "1", "2"))

dds <- DESeq(dds)

## Run differential expression
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
