library("DESeq2")

counts <- as.matrix(read.csv(
  file = "~/Data/ARCTIC/ARCTIC_CountMatrix.txt",
  sep = "\t",
  row.names = "Geneid"
))

counts_clin <- read.csv(
  file = "~/Data/ARCTIC/ARCTIC_Clinical_data_reduced_MP.csv",
  row.names = 1
)

counts_clin

counts_clin[1]

colnames(counts) <- sub("T[I]*_.+$", "", colnames(counts))

all(counts_clin[1] %in% colnames(counts))


dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~condidtion
)

head(counts, 2)
