library(tidyverse)

setwd("/Users/bjp/Desktop/Chapter 3 Rhizosphere/Tnseq/Essential/")

list.files()

WT <- read_csv(file = "WT.tradis_insertion_plot.tradis_gene_insert_sites.csv.all.csv",col_names = T, trim_ws = T)
WT <- WT %>% rename_all(paste0, ".WT")
WT <- WT %>% rename(locus_tag = locus_tag.WT)
NS <- read_csv(file = "NS.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv", col_names = T, trim_ws = T)
NS <- NS %>% rename_all(paste0, ".NS")
NS <- NS %>% rename(locus_tag = locus_tag.NS)

#Prepare annotations from eggNOG and KEGG
COGS <- read_tsv(file = "../Annotations/R7A2020.eggNOG.tsv", col_names = F, trim_ws = T, comment = "#")
COGS <- COGS %>% rename(locus_tag = X1, gene = X6, COG = X21, CA = X22)
COGS <- COGS %>% select(locus_tag, gene, COG, CA)
KEGG <- read_table(file = "../Annotations/R7A2020.KO.KofamKOALA.txt", col_names = F, comment = "#")
KEGG <- KEGG %>% rename(locus_tag = X2, KO = X3, KA = X7)
#Filter duplicate KEGG entries for smallest E-value
KEGG <- KEGG[order(KEGG$locus_tag, abs(KEGG$X6) ), ] ### sort first
KEGG <- KEGG[ !duplicated(KEGG$locus_tag), ]  ### Keep highest
KEGG <- KEGG %>%  select(locus_tag, KO, KA)

ANNOT <- left_join(KEGG, COGS)

WT.ES <- read_table(file = "WT.TraDIS.ICEFilt.ES.Unique.txt", col_names = F)
WT.ES$state <- "WT-ES"
WT.ES <- WT.ES %>% rename(locus_tag = X1)

NS.ES <- read_table(file = "NS.TraDIS.ICEFilt.ES.Unique.txt", col_names = F)
NS.ES$state <- "NS-ES"
NS.ES <- NS.ES %>% rename(locus_tag = X1)

CORE.ES <- read_table(file = "CORE.TraDIS.ICEFilt.ES.Unique.txt", col_names = F)
CORE.ES$state <- "CORE-ES"
CORE.ES <- CORE.ES %>% rename(locus_tag = X1)

ICE.ES <- read_table(file = "ICE.TraDIS.ES.Unique.txt", col_names = F)
ICE.ES$state <- "ICE-ES"
ICE.ES <- ICE.ES %>% rename(locus_tag = X1)
ES <- rbind(CORE.ES, WT.ES, NS.ES, ICE.ES)

#Join the TraDIS data to the ES genes list

ES <- left_join(ES, WT)
ES <- left_join(ES, NS)
ES <- left_join(ES, ANNOT)

write_delim(ES, "Table S1 Essential Genes of M japonicum R7A.csv", delim = ",", col_names = T, append = F)


WTNSCorr <- ES %>% filter(state %in% c("NS-ES","WT-ES")) %>% select(locus_tag, state, ins_index.WT, ins_index.NS)

WTNSCorr %>% ggplot(aes(y=ins_index.WT, x=ins_index.NS)) + geom_point() + geom_smooth(method=lm)

#reading in the curated funcational categories table
COGTable <- read_tsv("ES genes functional summary.tsv", col_names = T, trim_ws = T)
COGTable <- COGTable %>% mutate(FUN = as.factor(FUN))
COGTable <- COGTable %>% group_by(state) %>% select(FUN) %>% table() %>% as.data.frame.array()
COGTable$state <- rownames(COGTable)
write_csv(COGTable, "COG Table.csv")
