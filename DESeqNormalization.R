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

head(counts, 2)
