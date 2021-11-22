library(tidyverse)
library(DESeq2)
library(DEGreport)
#library(ggpubr)


setwd("~/Desktop/Chapter 4 Rhizosphere/Rhizosphere Tnseq/Differential Abundance/DESeq2 Pooled Samples/")

#######################################################
### Rarefaction Analysis of Rhizosphere Tn-seq data ###
#######################################################

wigFiles <- list.files(pattern = "filt.3.wig")
wigData <- data_frame(filename = wigFiles) %>% mutate(file_contents = map(filename, ~read_delim(., delim = " ", skip = 1, col_names = c("pos", "count"))))
wigTable <- unnest(wigData, cols = c(file_contents))
wigTable <- wigTable %>% separate(col = filename, 
                                  sep = "[.]", remove = F, 
                                  into = c("library", NA, NA, NA))
wigTable <- wigTable %>% pivot_wider(values_from = count, names_from = library)
wigTable <- wigTable %>% select(-filename)
wigTable <- wigTable %>% rename(NSI1 = `4305P2-01`,
                                NSI3 = `4305P2-02`,
                                WTI1 = `4305P2-03`,
                                WTI3 = `4305P2-04`,
                                NS1NFR5 = `4305P2-05`,
                                NS3NFR5 = `4305P2-06`,
                                NS1GIFU = `4305P2-07`,
                                NS3GIFU = `4305P2-08`,
                                WT1NFR5 = `4305P2-09`,
                                WT3NFR5 = `4305P2-10`,
                                WT1GIFU = `4305P2-11`,
                                WT3GIFU = `4305P2-12`,
                                NS1COL0 = `4305P2-13`,
                                NS3COL0 = `4305P2-14`,
                                WT1COL0 = `4305P2-15`,
                                WT3COL0 = `4305P2-16`)
wigUniq <- as.data.frame(unique(wigTable$pos))
wigUniq <- wigUniq %>% transmute(pos = `unique(wigTable$pos)`)
wigUniq <- wigTable %>% select(pos, NSI1) %>% filter(NSI1 != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, NSI3) %>% filter(NSI3 != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, WTI1) %>% filter(WTI1 != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, WTI3) %>% filter(WTI3 != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, NS1NFR5) %>% filter(NS1NFR5 != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, NS3NFR5) %>% filter(NS3NFR5 != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, NS1GIFU) %>% filter(NS1GIFU != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, NS3GIFU) %>% filter(NS3GIFU != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, WT1NFR5) %>% filter(WT1NFR5 != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, WT3NFR5) %>% filter(WT3NFR5 != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, WT1GIFU) %>% filter(WT1GIFU != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, WT3GIFU) %>% filter(WT3GIFU != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, NS1COL0) %>% filter(NS1COL0 != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, NS3COL0) %>% filter(NS3COL0 != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, WT1COL0) %>% filter(WT1COL0 != "NA") %>% left_join(x = wigUniq)
wigUniq <- wigTable %>% select(pos, WT3COL0) %>% filter(WT3COL0 != "NA") %>% left_join(x = wigUniq)

rownames(wigUniq) <- paste("pos", wigUniq$pos, sep = "_")
wigUniq[is.na(wigUniq)] <- 0
wigUniq <- wigUniq %>% select(-pos)
write.table(x = wigUniq,
            file = "../rhizosphere.mutant.abundance.txt",
            row.names = T,
            col.names = T,
            sep = "\t",
            quote = F)


##################################################################
##              LoF Analysis Using DESeq2 LRT Analysis          ##
##################################################################

setwd("~/Desktop/Chapter 4 Rhizosphere/Rhizosphere Tnseq/Differential Abundance/DESeq2 Pooled Samples/")

# READ IN TRADIS OUTPUTS WITH ICE GENES FILTERED OUT FOR CHROMOSOMAL COMPARISON
tradisFiles <- list.files(pattern = "tradis_gene_insert_sites.ICE.Filt.csv.all.csv")
tradisData <- data_frame(filename = tradisFiles) %>%  mutate(file_contents = map(filename, ~read_csv(.)))
tradisTable <- unnest(tradisData)
tradisTable <- tradisTable %>% separate(col = filename,
                                        sep = "[.]", remove = F,
                                        into = c("treatment", NA, NA, NA, NA, NA, NA, NA, NA, NA))
tradisTable <- tradisTable %>% group_by(treatment)

#Building table of read counts
read_counts <- tradisTable %>%
        select(locus_tag, treatment, read_count) %>%
        spread(key = treatment, value = read_count)
read_counts <- read_counts %>% select(locus_tag,
                                      WTGR1, WTGR2, NSGR1, NSGR2,
                                      WTTY1, WTTY2, NSTY1, NSTY2,
                                      WT1INPUT = `4305P2-03`,
                                      WT3INPUT = `4305P2-04`,
                                      NS1INPUT = `4305P2-01`,
                                      NS3INPUT = `4305P2-02`,
                                      WT1GIFU = `4305P2-11`,
                                      WT3GIFU = `4305P2-12`,
                                      NS1GIFU = `4305P2-07`,
                                      NS3GIFU = `4305P2-08`,
                                      WT1NFR5 = `4305P2-09`,
                                      WT3NFR5 = `4305P2-10`,
                                      NS1NFR5 = `4305P2-05`,
                                      NS3NFR5 = `4305P2-06`,
                                      WT1COL0 = `4305P2-15`,
                                      WT3COL0 = `4305P2-16`,
                                      NS1COL0 = `4305P2-13`,
                                      NS3COL0 = `4305P2-14`)

#### Build a DESEQ2 Object to normalize the reads counts with ###
read_counts <- read_counts %>% select(locus_tag, 
                                      # WTGR1, WTGR2, NSGR1, NSGR2,
                                      # WTTY1, WTTY2, NSTY1, NSTY2,
                                      WT1INPUT, WT3INPUT, NS1INPUT, NS3INPUT,
                                      WT1GIFU, WT3GIFU, NS1GIFU, NS3GIFU,
                                      WT1NFR5, WT3NFR5, NS1NFR5, NS3NFR5,
                                      WT1COL0, WT3COL0, NS1COL0, NS3COL0)

tnseq_conditions <- c(
        # "WTGR", "WTGR", "NSGR", "NSGR",
        # "WTTY", "WTTY", "NSTY", "NSTY",
        "WTINPUT", "WTINPUT", "NSINPUT", "NSINPUT",
        "WTGIFU", "WTGIFU", "NSGIFU", "NSGIFU",
        "WTNFR5", "WTNFR5", "NSNFR5", "NSNFR5",
        "WTCOL0", "WTCOL0", "NSCOL0", "NSCOL0")

### Prepare DESeqDataSet object ###
# Tidy the count data
colnames(read_counts)
rhizo <- read.csv("~/Desktop/Chapter 4 Rhizosphere/Table 4.S1 Summary of Rhizosphere Tnseq Analysis.csv")
rhizo_NE_loci <- rhizo %>% filter(DESeq2.LRT == "TRUE") %>% select(locus_tag)

# filtering the counts table for gene in which both input libraries have greater than 0 counts
# filtered_read_counts <- read_counts %>%
#         filter(WT1INPUT > 0 & WT3INPUT > 0) %>%
#         filter(NS1INPUT > 0 & NS3INPUT > 0)
filtered_read_counts <- read_counts %>% filter(locus_tag %in% rhizo_NE_loci$locus_tag)

ID <- filtered_read_counts$locus_tag
filtered_read_counts <- filtered_read_counts %>% select(-locus_tag)
filtered_read_counts <- as.matrix(filtered_read_counts)
rownames(filtered_read_counts) <- ID

condition <- factor(tnseq_conditions)
condition <- DataFrame(condition)

# Make the dds Object
dds_filtered <- DESeqDataSetFromMatrix(countData = filtered_read_counts, colData = condition, ~ condition)

### DESeq2 analysis
design(dds_filtered) = ~ condition
dds_filtered <- estimateSizeFactors(dds_filtered)
sizeFactors(dds_filtered)

dds_filtered <- DESeq(dds_filtered, test = "LRT", reduced = ~ 1)
dds_results <- results(dds_filtered)
dds_results <- as.data.frame(dds_results)
dds_results$locus_tag <- row.names(dds_results)
dds_results <- dds_results %>%
        select(locus_tag, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
filt_norm_counts <- as.data.frame(counts(dds_filtered, normalized=T))
filt_norm_counts <- filt_norm_counts %>% rename_all(paste0, ".norm_counts")
filt_norm_counts$locus_tag <- row.names(filt_norm_counts)

### Clustering

sig_res_LRT <- dds_results %>%
        data.frame() %>%
        rownames_to_column() %>% 
        as_tibble(var='gene') %>% 
        filter(padj < 0.05)

sigLRT_genes <- sig_res_LRT %>% 
        pull(locus_tag)
length(sigLRT_genes)

clustering_sig_genes <- sig_res_LRT %>%
        arrange(padj) %>%
        head(n=1000)
for_clusters <- assay(dds_filtered)
for_clusters <- for_clusters[clustering_sig_genes$locus_tag, ]
library(DEGreport)
library(SummarizedExperiment)
meta <- colData(dds_filtered)
clusters <- degPatterns(as.matrix(for_clusters), metadata = meta, time = "condition", col = NULL, minc = 1)
sig_clusters <- as_tibble(clusters$df)
sig_clusters <- sig_clusters %>% rename(genes = "locus_tag")

### LoF Analysis output

rhizo <- left_join(rhizo, filt_norm_counts, "locus_tag")
rhizo <- left_join(rhizo, dds_results, "locus_tag")
rhizo <- left_join(rhizo, sig_clusters, "locus_tag")
#write_csv(x = rhizo, file = "~/Desktop/Chapter 4 Rhizosphere/Table 4.S2 Summary of Rhizosphere DESeq2 Tnseq Analysis.csv")

### Making the PCA PLOT with the NE genes ###

sig_loci <- rhizo %>% filter(padj < 0.05) %>% select(locus_tag)
setwd("~/Desktop/Chapter 4 Rhizosphere/Rhizosphere Tnseq/Differential Abundance/DESeq2 Pooled Samples/")
tradisFiles <- list.files(pattern = "tradis_gene_insert_sites.csv.all.csv")
tradisData <- data_frame(filename = tradisFiles) %>%  mutate(file_contents = map(filename, ~read_csv(.)))
tradisTable <- unnest(data = tradisData, cols = `file_contents`)
tradisTable <- tradisTable %>% separate(col = filename, 
                                        sep = "[.]", remove = F, 
                                        into = c("treatment", NA, NA, NA, NA, NA, NA, NA, NA, NA))
tradisTable <- tradisTable %>% group_by(treatment)

read_counts <- tradisTable %>% 
        select(locus_tag, treatment, read_count) %>% 
        spread(key = treatment, value = read_count)

read_counts <- read_counts %>% select(
        locus_tag,
        WTGR1, WTGR2, WTTY1, WTTY2,
        NSI1 = `4305P2-01`,
        NSI3 = `4305P2-02`,
        WTI1 = `4305P2-03`,
        WTI3 = `4305P2-04`,
        NS1NFR5 = `4305P2-05`,
        NS3NFR5 = `4305P2-06`,
        NS1GIFU = `4305P2-07`,
        NS3GIFU = `4305P2-08`,
        WT1NFR5 = `4305P2-09`,
        WT3NFR5 = `4305P2-10`,
        WT1GIFU = `4305P2-11`,
        WT3GIFU = `4305P2-12`,
        NS1COL0 = `4305P2-13`,
        NS3COL0 = `4305P2-14`,
        WT1COL0 = `4305P2-15`,
        WT3COL0 = `4305P2-16`
)


read_counts <- read_counts %>% select(locus_tag,
                                      #WTGR1, WTGR2, WTTY1, WTTY2,
                                      WTI1,
                                      WTI3,
                                      NSI1, NSI3,
                                      WT1NFR5,
                                      WT3NFR5,
                                      NS1NFR5,
                                      NS3NFR5,
                                      WT1GIFU,
                                      WT3GIFU,
                                      NS1GIFU,
                                      NS3GIFU,
                                      WT1COL0,
                                      WT3COL0,
                                      NS1COL0,
                                      NS3COL0)

condition <- c(
        #"WTGR", "WTGR", "WTTY", "WTTY",
        "WTINPUT", "WTINPUT",
        "NSINPUT", "NSINPUT",
        "WTNFR5", "WTNFR5",
        "NSNFR5", "NSNFR5",
        "WTGIFU", "WTGIFU",
        "NSGIFU", "NSGIFU",
        "WTCOL0", "WTCOL0",
        "NSCOL0", "NCOL0")

read_counts <- read_counts %>% filter(locus_tag %in% sig_loci$locus_tag)
ID <- read_counts$locus_tag
read_counts <- read_counts %>% select(-locus_tag)
read_counts <- as.matrix(read_counts)
rownames(read_counts) <- ID
condition <- factor(tnseq_conditions)
condition <- DataFrame(condition)

dds_filtered <- DESeqDataSetFromMatrix(countData = filtered_read_counts, colData = condition, ~ condition)

design(dds_filtered) = ~ condition
dds_filtered <- estimateSizeFactors(dds_filtered)
sizeFactors(dds_filtered)

dds_filtered <- DESeq(dds_filtered, test = "LRT", reduced = ~ 1)

vsdata <- varianceStabilizingTransformation(dds_filtered, blind = T)
str(vsdata$condition)
vsdata$condition <- factor(vsdata$condition, levels = c("WTINPUT", "NSINPUT",
                                                        "WTGIFU", "NSGIFU",
                                                        "WTNFR5", "NSNFR5",
                                                        "WTCOL0", "NSCOL0"))

PCA <- plotPCA(vsdata, intgroup="condition")
PCA <- PCA + theme_minimal(base_size = 10) 
PCA <- PCA + geom_point(size=6)
PCA
PCA + ggsave("~/Desktop/Chapter 4 Rhizosphere/PCA Plot of 114 sig diff VST Normalized Counts.svg", width = 12, height = 8)


# ### Prepare Master Table for Chapter 4 ###
# 
# setwd("~/Desktop/Chapter 4 Rhizosphere/Rhizosphere Tnseq/Differential Abundance/")
# # Read in the CE Table with functions.
# PLANTA <- read_csv(file = "../../Table S1 Rhizosphere TraDIS Essential Genes.csv",
#                   col_names = T,
#                   trim_ws = T)
# PLANTA <- PLANTA %>% select(locus_tag, state, FUN) %>% 
#         filter(state %in% c("PL-NS-UE", "PL-CORE-ES", "PL-WT-UE"))
# PLANTA <- PLANTA %>% transmute(locus_tag, state.planta = state, FUN.planta = FUN)
# 
# # Read in ANNOT files to append
# COGS <- read_tsv(file = "../Annotations/R7A2020.eggNOG.tsv", col_names = F, trim_ws = T, comment = "#")
# COGS <- COGS %>% select(locus_tag = X1, gene = X6, COG = X21, CA = X22)
# KEGG <- read_table(file = "../Annotations/R7A2020.KO.KofamKOALA.txt", col_names = F, comment = "#")
# KEGG <- KEGG %>% select(locus_tag = X2, KO = X3, KA = X7)
# #Filter duplicate KEGG entries for smallest E-value
# KEGG <- KEGG[order(KEGG$locus_tag, abs(KEGG$X6) ), ] ### sort first
# KEGG <- KEGG[ !duplicated(KEGG$locus_tag), ]  ### Keep highest
# KEGG <- KEGG %>%  select(locus_tag, KO, KA)
# 
# NCBI <- read_csv("DESeq2 Pooled Samples/4305P2-03.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv", col_names = T, trim_ws = )
# NCBI <- NCBI %>% select(locus_tag, NCBI = fcn)
# #Merge the annotations
# ANNOT <- left_join(NCBI, KEGG)
# ANNOT <- left_join(ANNOT, COGS)
# 
# # Read in TraDIS tables and filter and merge for ins, read counts, insertion index
# insIdxTable <- tradisTable %>% select(locus_tag, treatment, ins_index) %>% spread(key = treatment, value = ins_index)
# insIdxTable <- insIdxTable %>% rename_all(paste0, ".ins_index")
# insIdxTable$locus_tag <- insIdxTable$locus_tag.ins_index
# insIdxTable <- insIdxTable %>% select(-locus_tag.ins_index)
# 
# # Read in Table S1 from Chapter 3 for a backbone(?)
# chpt3 <- read_tsv("~/Desktop/Chapter 3 In Vitro/Table S1 Summary of In Vitro Tnseq Analysis.tsv",
#                   col_names = T,
#                   trim_ws = T)
# chpt3 <- chpt3 %>% select(locus_tag, 
#                           state.ES = state, 
#                           FUN.ES = FUN, 
#                           read_count.WTES, 
#                           ins_index.WTES, 
#                           read_count.NSES, 
#                           ins_index.NSES)
# # Build up the table from 6308 count backbone ANNOT, PLANTA, Abundance, Table S1
# 
# chpt4 <- as.data.frame(insIdxTable$locus_tag)
# chpt4 <- chpt4 %>% select(locus_tag = `insIdxTable$locus_tag`)
# chpt4 <- left_join(x = chpt4, y = ANNOT)
# chpt4 <- left_join(x = chpt4, y = chpt3)
# chpt4 <- left_join(x = chpt4, PLANTA)
# chpt4 <- left_join(x = chpt4, y = insIdxTable, "locus_tag")
# chpt4 <- left_join(x = chpt4, y = dds_results, "locus_tag")
# chpt4 <- left_join(x = chpt4, y = filt_norm_counts, "locus_tag")
# 
# write_delim(x = chpt4,
#             file = "../../Table S2 Summary of Rhizosphere DESeq2 Tnseq Analysis.tsv",
#             delim = "\t",
#             col_names = T,
#             append = F)
# # ^^^ Calculate the L2FC in excel for Table S1 an then make flags for the state of significance

################################
### ICESym Normalized Counts ###
################################
# READ IN TRADIS OUTPUTS WITH ICE GENES FOR CHROMOSOMAL COMPARISON
setwd("~/Desktop/Chapter 4 Rhizosphere/Rhizosphere Tnseq/Differential Abundance/DESeq2 Pooled Samples/")
tradisFiles <- list.files(pattern = "tradis_gene_insert_sites.csv.all.csv")
tradisData <- data_frame(filename = tradisFiles) %>%  mutate(file_contents = map(filename, ~read_csv(.)))
tradisTable <- unnest(data = tradisData, cols = `file_contents`)
tradisTable <- tradisTable %>% separate(col = filename, 
                                        sep = "[.]", remove = F, 
                                        into = c("treatment", NA, NA, NA, NA, NA, NA, NA, NA, NA))
tradisTable <- tradisTable %>% group_by(treatment)

read_counts <- tradisTable %>% 
        select(locus_tag, treatment, read_count) %>% 
        spread(key = treatment, value = read_count)

read_counts <- read_counts %>% select(
        locus_tag,
        WTGR1, WTGR2, WTTY1, WTTY2,
        NSI1 = `4305P2-01`,
        NSI3 = `4305P2-02`,
        WTI1 = `4305P2-03`,
        WTI3 = `4305P2-04`,
        NS1NFR5 = `4305P2-05`,
        NS3NFR5 = `4305P2-06`,
        NS1GIFU = `4305P2-07`,
        NS3GIFU = `4305P2-08`,
        WT1NFR5 = `4305P2-09`,
        WT3NFR5 = `4305P2-10`,
        WT1GIFU = `4305P2-11`,
        WT3GIFU = `4305P2-12`,
        NS1COL0 = `4305P2-13`,
        NS3COL0 = `4305P2-14`,
        WT1COL0 = `4305P2-15`,
        WT3COL0 = `4305P2-16`
)


read_counts <- read_counts %>% select(locus_tag,
                                      #WTGR1, WTGR2, WTTY1, WTTY2,
                                      WTI1,
                                      WTI3,
                                      NSI1, NSI3,
                                      WT1NFR5,
                                      WT3NFR5,
                                      NS1NFR5,
                                      NS3NFR5,
                                      WT1GIFU,
                                      WT3GIFU,
                                      NS1GIFU,
                                      NS3GIFU,
                                      WT1COL0,
                                      WT3COL0,
                                      NS1COL0,
                                      NS3COL0)

condition <- c(
        #"WTGR", "WTGR", "WTTY", "WTTY",
               "WTINPUT", "WTINPUT",
               "NSINPUT", "NSINPUT",
               "WTNFR5", "WTNFR5",
               "NSNFR5", "NSNFR5",
               "WTGIFU", "WTGIFU",
               "NSGIFU", "NSGIFU",
               "WTCOL0", "WTCOL0",
               "NSCOL0", "NCOL0")
condition <- as.data.frame(condition)

ID <- read_counts$locus_tag
read_counts <- read_counts %>% select(-locus_tag)
rownames(read_counts) <- ID

dds_all <- DESeqDataSetFromMatrix(countData = read_counts, colData = condition, ~ condition)
design(dds_all) = ~ condition
dds_all <- estimateSizeFactors(dds_all)
sizeFactors(dds_all)

dds_all <- DESeq(dds_all, test="LRT", reduced = ~ 1)
dds_all_results <- results(dds_all)
dds_all_results <- as.data.frame(dds_all_results)

all_norm_counts <- as.data.frame(counts(dds_all, normalized=T))


### Filter for ICESym genes
ICEloci <- read_table("~/ref/ICESym.locus.tags.txt", col_names = F)
ICEloci <- ICEloci %>% select(locus_tag = X1)
ICE_counts <- read_counts %>% filter(locus_tag %in% as.list(ICEloci$locus_tag))
ICErhizo <- read.csv("~/Desktop/Chapter 4 Rhizosphere/Table 4.S1 Rhizosphere TraDIS Essential Genes.csv")
ICE_NE_loci <- ICErhizo %>% filter(ICE.state == "ICE-NE") %>% select(locus_tag)
ICE_counts <- ICE_counts %>% filter(locus_tag %in% ICE_NE_loci$locus_tag)

ID <- ICE_counts$locus_tag
ICE_counts <- ICE_counts %>% select(-locus_tag)
rownames(ICE_counts) <- ID

# Make the dds Object
dds_ICE <- DESeqDataSetFromMatrix(countData = ICE_counts, colData = condition, ~ condition)

### DESeq2 analysis
design(dds_ICE) = ~ condition
dds_ICE <- estimateSizeFactors(dds_ICE)
sizeFactors(dds_ICE)

dds_ICE <- DESeq(dds_ICE, test = "LRT", reduced = ~ 1)
dds_ICE_results <- results(dds_ICE)
dds_ICE_results <- as.data.frame(dds_ICE_results)
dds_ICE_results$locus_tag <- row.names(dds_ICE_results)
dds_ICE_results <- dds_ICE_results %>% 
        select(locus_tag, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)


ICE_norm_counts <- as.data.frame(counts(dds_ICE, normalized=T))
ICE_norm_counts <- ICE_norm_counts %>% rename_all(paste0, ".norm_counts")
ICE_norm_counts$locus_tag <- row.names(ICE_norm_counts)

dds_ICE_results <- left_join(dds_ICE_results, ICE_norm_counts)
dds_ICE_results <- left_join(dds_ICE_results, ANNOT, "locus_tag")

### Clustering
library(DEGreport)
library(SummarizedExperiment)

clustering_sig_genes <- dds_ICE_results %>% filter(padj <0.05)

for_clusters <- assay(dds_ICE)
for_clusters <- for_clusters[clustering_sig_genes$locus_tag, ]

meta <- colData(dds_ICE)
clusters <- degPatterns(as.matrix(for_clusters), metadata = meta, time = "condition", col = NULL, minc = 1)
sig_clusters <- as_tibble(clusters$df)
sig_clusters <- sig_clusters %>% rename(genes = "locus_tag")

### LoF Analysis output

ICE_norm_counts <- left_join(ICE_norm_counts, dds_ICE_results, "locus_tag")
ICE_norm_counts <- left_join(ICE_norm_counts, sig_clusters, "locus_tag")

write_csv(
        x = ICE_norm_counts,
        file = "~/Desktop/Chapter 4 Rhizosphere/Table 4.S3  ICESym DESeq2 LoF Analysis.csv",
        col_names = T,
        append = F
)

vsdata <- varianceStabilizingTransformation(dds_ICE, blind=T)
PCA <- plotPCA(vsdata, intgroup="condition")
PCA <- PCA + theme_minimal(base_size = 10) 
PCA <- PCA + geom_point(size=1, alpha = .9)
PCA


###################################################
### Make LoF table for Figure 4.6 Metabolic Map ###
###################################################

setwd("~/Desktop/Chapter 4 Rhizosphere/Figure 6 Overview of Rhizosphere Central Carbon and Nitrogen Metabolism/")
mapTable <- read_tsv(file = "Table 3 Metabolic Map Table.txt", col_names = T, trim_ws = T)
chpt4 <- readxl::read_excel("../Table S2 Summary of Rhizosphere DESeq2 Tnseq Analysis.xlsx", sheet = "Chromosomal", col_names = T)
chpt4 <- chpt4 %>% select(locus_tag, baseMean, mean.INPUT, mean.GIFU, padj, L2FC, DESEQ2)

chpt4ICE <- readxl::read_excel("../Table S2 Summary of Rhizosphere DESeq2 Tnseq Analysis.xlsx", sheet = "ICESym", col_names = T)
chpt4ICE <- chpt4ICE %>% select(locus_tag, baseMean, mean.INPUT, mean.GIFU, padj, L2FC, DESEQ2)

chpt4 <- rbind(chpt4, chpt4ICE)

mapTable <- mapTable %>% select(locus_tag, Pathway, ID, state, FUN, gene, NCBI, KO, KA, COG, CA, start, end, gene_length, strand)

mapTable <- left_join(mapTable, chpt4, "locus_tag")

tradisPlanta <- readxl::read_excel("../Table S1 Rhizosphere TraDIS Essential Genes.xlsx", col_names = T)
tradisPlanta <- tradisPlanta %>% select(locus_tag, state.planta = state)
mapTable <- left_join(mapTable, tradisPlanta)

write_delim(x = mapTable,
            file = "../Table 4.3 Metabolic Map Summary Table.tsv",
            delim = "\t",
            col_names = T,
            append = F)

# Add on the amended table PLANTA states

# mapTable <- readxl::read_excel("Table 3 Metabolic Map Table.xlsx", col_names = T)
# tradisPlanta <- readxl::read_excel("../Table S1 Rhizosphere TraDIS Essential Genes.xlsx", col_names = T)
# tradisPlanta <- tradisPlanta %>% select(locus_tag, state.planta = state)
# mapTable <- left_join(mapTable, tradisPlanta)
# write_delim(x = mapTable,
#             file = "../Table 4.3 Metabolic Map Summary Table.tsv",
#             delim = "\t",
#             col_names = T,
#             append = F)

#############################################
### PCA Plot of all genes for rhizosphere ###
#############################################
setwd("~/Desktop/Chapter 4 Rhizosphere/Rhizosphere Tnseq/Differential Abundance/Origanl DESeq2 Analysis/deseq/")

# READ IN TRADIS OUTPUTS WITH ICE GENES FILTERED OUT FOR CHROMOSOMAL COMPARISON
tradisFiles <- list.files(pattern = "tradis_gene_insert_sites.csv.all.csv")
tradisData <- data_frame(filename = tradisFiles) %>%  mutate(file_contents = map(filename, ~read_csv(.)))
tradisTable <- unnest(tradisData)
tradisTable <- tradisTable %>% separate(col = filename, 
                                        sep = "[.]", remove = F, 
                                        into = c("treatment", NA, NA, NA, NA, NA, NA, NA, NA, NA))
tradisTable <- tradisTable %>% group_by(treatment)

#Building table of read counts
read_counts <- tradisTable %>% 
        select(locus_tag, treatment, read_count) %>% 
        spread(key = treatment, value = read_count)

read_counts <- read_counts %>% select(locus_tag,
                                      WT1INPUT, WT3INPUT, NS1INPUT, NS3INPUT,
                                      WT1GIFU, WT3GIFU, NS1GIFU, NS3GIFU,
                                      WT1NFR5, WT3NFR5, NS1NFR5, NS3NFR5,
                                      WT1COL0, WT3COL0, NS1COL0, NS3COL0)

tnseq_conditions <- c("WTINPUT", "WTINPUT", "NSINPUT", "NSINPUT",
                      "WTGIFU", "WTGIFU", "NSGIFU", "NSGIFU",
                      "WTNFR5", "WTNFR5", "NSNFR5", "NSNFR5",
                      "WTCOL0", "WTCOL0", "NSCOL0", "NSCOL0")
# Tidy the count data
colnames(read_counts)

# filtering the counts table for gene in which both input libraries have greater than 0 counts
filtered_read_counts <- read_counts %>% filter(NS1INPUT > 1 & NS3INPUT > 1) %>% filter(WT1INPUT > 1 & WT3INPUT > 1)
ID <- filtered_read_counts$locus_tag
filtered_read_counts <- filtered_read_counts %>% select(-locus_tag)
filtered_read_counts <- as.matrix(filtered_read_counts)
rownames(filtered_read_counts) <- ID

condition <- factor(tnseq_conditions)
condition <- DataFrame(condition)

# Make the dds Object
dds_filtered <- DESeqDataSetFromMatrix(countData = filtered_read_counts, colData = condition, ~ condition)
dds_filtered <- estimateSizeFactors(dds_filtered)

vsdata <- varianceStabilizingTransformation(dds_filtered, blind = T)
PCA <- plotPCA(vsdata, intgroup="condition", ntop = 250)
PCA <- PCA + theme_minimal(base_size = 10) 
PCA <- PCA + geom_point(size=8)
PCA <- PCA + ylim(-3.5, 3.5)
PCA <- PCA + xlim(-4.5, 2.5)
PCA

PCA + ggsave("../../../../Figure 5 Rhizosphere Differential Abundance Analysis with LoF Competitions/PCA Plot of top 250 variable VST Normalized Counts.svg", width = 12, height = 8)


#### Make genomic norm counts for merger with the rhizosphere Tn-seq summary table.
setwd("~/Desktop/Chapter 4 Rhizosphere/Rhizosphere Tnseq/Differential Abundance/DESeq2 Pooled Samples/")
tradisFiles <- list.files(pattern = "tradis_gene_insert_sites.csv.all.csv")
tradisData <- data_frame(filename = tradisFiles) %>%  mutate(file_contents = map(filename, ~read_csv(.)))
tradisTable <- unnest(data = tradisData, cols = `file_contents`)
tradisTable <- tradisTable %>% separate(col = filename, 
                                        sep = "[.]", remove = F, 
                                        into = c("treatment", NA, NA, NA, NA, NA, NA, NA, NA, NA))
tradisTable <- tradisTable %>% group_by(treatment)

read_counts <- tradisTable %>% 
        select(locus_tag, treatment, read_count) %>% 
        spread(key = treatment, value = read_count)

read_counts <- read_counts %>% select(
        locus_tag,
        WTGR1, WTGR2, WTTY1, WTTY2,
        NSI1 = `4305P2-01`,
        NSI3 = `4305P2-02`,
        WTI1 = `4305P2-03`,
        WTI3 = `4305P2-04`,
        NS1NFR5 = `4305P2-05`,
        NS3NFR5 = `4305P2-06`,
        NS1GIFU = `4305P2-07`,
        NS3GIFU = `4305P2-08`,
        WT1NFR5 = `4305P2-09`,
        WT3NFR5 = `4305P2-10`,
        WT1GIFU = `4305P2-11`,
        WT3GIFU = `4305P2-12`,
        NS1COL0 = `4305P2-13`,
        NS3COL0 = `4305P2-14`,
        WT1COL0 = `4305P2-15`,
        WT3COL0 = `4305P2-16`
)


read_counts <- read_counts %>% select(locus_tag,
                                      #WTGR1, WTGR2, WTTY1, WTTY2,
                                      WTI1,
                                      WTI3,
                                      NSI1, NSI3,
                                      WT1NFR5,
                                      WT3NFR5,
                                      NS1NFR5,
                                      NS3NFR5,
                                      WT1GIFU,
                                      WT3GIFU,
                                      NS1GIFU,
                                      NS3GIFU,
                                      WT1COL0,
                                      WT3COL0,
                                      NS1COL0,
                                      NS3COL0)

condition <- c(
        #"WTGR", "WTGR", "WTTY", "WTTY",
        "WTINPUT", "WTINPUT",
        "NSINPUT", "NSINPUT",
        "WTNFR5", "WTNFR5",
        "NSNFR5", "NSNFR5",
        "WTGIFU", "WTGIFU",
        "NSGIFU", "NSGIFU",
        "WTCOL0", "WTCOL0",
        "NSCOL0", "NCOL0")
condition <- as.data.frame(condition)

ID <- read_counts$locus_tag
read_counts <- read_counts %>% select(-locus_tag)
rownames(read_counts) <- ID

dds_all <- DESeqDataSetFromMatrix(countData = read_counts, colData = condition, ~ condition)
design(dds_all) = ~ condition
dds_all <- estimateSizeFactors(dds_all)
sizeFactors(dds_all)

dds_all <- DESeq(dds_all, test="LRT", reduced = ~ 1)
dds_all_results <- results(dds_all)
dds_all_results <- as.data.frame(dds_all_results)
dds_all_results <- dds_all_results %>% rename_all(paste0, ".all_deseq2")
dds_all_results$locus_tag <- rownames(dds_all_results)

all_norm_counts <- as.data.frame(counts(dds_all, normalized=T))
all_norm_counts <- all_norm_counts %>% rename_all(paste0, ".all_norm_counts")
all_norm_counts$locus_tag <- row.names(all_norm_counts)

rhizo <- read_csv("~/Desktop/Chapter 4 Rhizosphere/Table 4.S1 Summary of Rhizosphere Tnseq Analysis.csv", col_names = T)

rhizo <- left_join(rhizo, all_norm_counts, "locus_tag")
rhizo <- left_join(rhizo, dds_all_results, "locus_tag")
write_csv(rhizo, "~/Desktop/deseq_all.test.csv")