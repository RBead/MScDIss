setwd("data/scratch/bt23766/SRAdata/")

library(methylKit)

# Define the path to .cov files
input_dir4 <- "/data/home/bt23766/covfiles"

# List the .cov files
cov_files4 <- list.files(input_dir4, pattern = "*.cov.gz", full.names = TRUE)

# Extract the section starting with "SRR" and remove the part from "_trimmed" onwards
srr_list4 <- sub(".*/(SRR[0-9]+)_.*", "\\1", cov_files4)


# Define sample IDs
sample_ids4 <- srr_list4  # Example sample IDs

# Define treatment vector (all samples in the same treatment group, so 1)
treatment4 <- rep(1, length(cov_files4)) 


# read the files to a methylRawListDB object: myobjDB 
# and save in databases in folder methylDB
myobjDB4=methRead(location=as.list(cov_files4),
                  sample.id=as.list(sample_ids4),
                  assembly="hg38",
                  context="CpG",
                  treatment=treatment4,
                  dbtype = "tabix",
                  pipeline = "bismarkCoverage",
                  mincov = 2
)

head(myobjDB4)

# Filter
filtered.myobj4=filterByCoverage(myobjDB4,lo.count=10,lo.perc=NULL,
                                 hi.count=NULL,hi.perc=99.9)

# Normalize
normalizeCoverage()
newObj4 = normalizeCoverage(myobjDB4,method="median")

head(newObj4)

# Unite
meth4=unite(newObj4, destrand=FALSE, min.per.group=3L, chunk.size = 1e+06)
meth4



perc.meth4=percMethylation(meth4, rowids=TRUE)

head(perc.meth4)

perc_meth_df4 <- as.data.frame(perc.meth4)
perc_meth_df4
nrow(perc_meth_df4)

#below is just for density graph:

testmeth_df <- perc_meth_df4 %>% select(-position)
testmeth_df <- subset(perc_meth_df4, select = -c(position))
testmeth_df

meth_df_long <- melt(testmeth_df, variable.name = "Sample", value.name = "Methylation")

ggplot(meth_df_long, aes(x = Methylation, color = Sample)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Methylation for subset of Samples",
       x = "Methylation %",
       y = "Density") 



#genes extraction from list below


# Read the TSV file
data <- read.delim("gene-aging-mechanisms.tsv", header = FALSE, quote = "", sep = "\t")

# Extract gene names using regular expressions
gene_names <- gsub('^"([^"]+)",.*', '\\1', data$V1)

# Combine gene names with original data
filtered_data <- data.frame(GeneName = gene_names, AdditionalColumn = data$V2)

# Write filtered data to a new TSV file
write.table(filtered_data, "filtered_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(gene_names, "filtered_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# Read the file into R
newgene_data <- read.csv("filtered_genes.tsv", sep = "\t", header = FALSE)

# Remove the quotes
gene_data_clean <- gsub('"', '', newgene_data$V1)

# Write the cleaned data back to a TSV file
write.table(gene_data_clean, "cleangenes.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




#old impute below

BiocManager::install("impute")
library(impute)

knn_imputed <- impute.knn(perc.meth4)$data
head(knn_imputed)


options(expressions = 500000)

split_data <- split(perc.meth, ceiling(seq_along(rownames(perc.meth)) / 1000))

imputed_chunks <- lapply(split_data, function(chunk) {impute.knn(chunk)$data})

imputed_data <- do.call(rbind, imputed_chunks)





BiocManager::install("annotatr", force = TRUE)

BiocManager::install("GenomicRanges")

BiocManager::install("AnnotationHub", update = TRUE, ask = FALSE)


library(annotatr)
library(GenomicRanges)
library(AnnotationHub)

meth_df4 <- getData(meth4)


meth_df4
nrow(meth_df4)

#need to mix meth_df and percmeth to get %s in there
perc_meth_df4
row.names(perc_meth_df) <- NULL

combined_df4 <- bind_cols(meth_df4, perc_meth_df4)
combined_df4

srr_columns4 <- combined_df4[, grep("^SRR133", names(combined_df4))]
srr_columns4

mcols = srr_columns4

meth_gr4 <- GRanges(
  seqnames = meth_df4$chr,
  ranges = IRanges(start = meth_df4$start, end = meth_df4$end),
  strand = meth_df4$strand, srr_columns4)

meth_gr4 #has methylation and ranges

#now new bit below from here, with annotatr and filtering

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("AnnotationHub", "BiocFileCache", "dbplyr", "dplyr"), update = TRUE, ask = FALSE)

library(AnnotationHub)
AnnotationHub::clearCache()
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


annots = c('hg38_cpgs', 'hg38_basicgenes')
annotations <- build_annotations(genome = 'hg38', annotations = annots)


dm_annotated4 = annotate_regions(
  regions = meth_gr4,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

print(dm_annotated4)
dm_annotated4


df_dm_annotated4 = data.frame(dm_annotated4)

print(head(df_dm_annotated4))

# Filtering
subset_df4 <- subset(df_dm_annotated4, annot.symbol %in% genes_of_interest)
subset_df4


unique_genes4 <- unique(subset_df4$annot.symbol)
print(unique_genes4)
print(paste("Number of unique genes:", length(unique_genes)))

unique_start_df4 <- subset_df4 %>% distinct(start, .keep_all = TRUE)
unique_start_df4
nrow(unique_start_df4)


#initial knn imputation

imputationcols4 <- unique_start_df4[, grep("^SRR", names(unique_start_df4), value = TRUE)]

knn_imputednewtest4 <- impute.knn(as.matrix(imputationcols4))$data

unique_start_df4[, grep("^SRR", names(unique_start_df4), value = TRUE)] <- knn_imputednewtest

head(unique_start_df4)


#new mean impu test:

# Select the columns to impute
imputationcolsnew <- unique_start_df4[, grep("^SRR", names(unique_start_df4), value = TRUE)]

# Perform mean imputation

mean_imputednewtest <- apply(imputationcolsnew, 2, function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  return(col)
})

# Replace the original columns with the imputed data
unique_start_df4[, grep("^SRR", names(unique_start_df4), value = TRUE)] <- mean_imputednewtest

# Display the head of the dataframe to verify the changes
head(unique_start_df4)

# Write data
write.csv(unique_start_df4, "testfinaldata4.csv", row.names = FALSE)


#test imputation for RMSE:

# Introduce missing values (for demonstration)
set.seed(123) # For reproducibility
mask <- matrix(runif(length(imputationcolsnew)) < 0.1, nrow=nrow(imputationcolsnew), ncol=ncol(imputationcolsnew))
imputationcols_with_na <- imputationcolsnew
imputationcols_with_na[mask] <- NA

# Perform mean imputation
mean_imputed_test <- apply(imputationcols_with_na, 2, function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  return(col)
})

# Calculate RMSE
original_values <- imputationcolsnew[mask] # Original values before NA introduction
imputed_values <- mean_imputed_test[mask]  # Imputed values

rmse <- sqrt(mean((original_values - imputed_values)^2))
print(paste("RMSE:", rmse))
"RMSE: 34.0278671000422"
