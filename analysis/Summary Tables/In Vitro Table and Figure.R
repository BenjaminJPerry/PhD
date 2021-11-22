library(tidyverse)
setwd("/Users/bjp/Desktop/Chapter 3 In Vitro/In Vitro Tnseq/InVitro/")

################################################################################
############ NS vs WT Conditionally Essential Functional Categories ############
################################################################################
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
#Merge the annotations
ANNOT <- left_join(KEGG, COGS)


# Reading in TraDIS output files
WTINPUT <- read_csv(file = "WT-INPUT.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv",col_names = T, trim_ws = T)
WTINPUT <- WTINPUT %>% rename_all(paste0, ".WTINPUT")
WTINPUT <- WTINPUT %>% rename(locus_tag = locus_tag.WTINPUT)
NSINPUT <- read_csv(file = "NS-INPUT.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv",col_names = T, trim_ws = T)
NSINPUT <- NSINPUT %>% rename_all(paste0, ".NSINPUT")
NSINPUT <- NSINPUT %>% rename(locus_tag = locus_tag.NSINPUT)

WTTY <- read_csv(file = "WTTY.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv",col_names = T, trim_ws = T)
WTTY <- WTTY %>% rename_all(paste0, ".WTTY")
WTTY <- WTTY %>% rename(locus_tag = locus_tag.WTTY)
NSTY <- read_csv(file = "NSTY.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv",col_names = T, trim_ws = T)
NSTY <- NSTY %>% rename_all(paste0, ".NSTY")
NSTY <- NSTY %>% rename(locus_tag = locus_tag.NSTY)

WTGR <- read_csv(file = "WTGR.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv",col_names = T, trim_ws = T)
WTGR <- WTGR %>% rename_all(paste0, ".WTGR")
WTGR <- WTGR %>% rename(locus_tag = locus_tag.WTGR)
NSGR <- read_csv(file = "NSGR.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv",col_names = T, trim_ws = T)
NSGR <- NSGR %>% rename_all(paste0, ".NSGR")
NSGR <- NSGR %>% rename(locus_tag = locus_tag.NSGR)

#Prepare annotations from eggNOG and KEGG
 
FUNCS <- read_csv(file = "CE genes functional summary.csv", col_names = T)
FUNCS <- FUNCS %>% select(locus_tag, FUN)
FUNCS <- FUNCS %>% distinct(locus_tag, .keep_all = TRUE)

#Reading in the conditionally essential genes
WT.INPUT <- read_table(file = "TB.WT.CE.txt", col_names = F)
WT.INPUT$state <- "WT-TB-CE"
WT.INPUT <- WT.INPUT %>% rename(locus_tag = X1)
NS.INPUT <- read_table(file = "TB.NS.CE.txt", col_names = F)
NS.INPUT$state <- "NS-TB-CE"
NS.INPUT <- NS.INPUT %>% rename(locus_tag = X1)
CE.INPUT <- read_table(file = "TB.CE.txt", col_names = F)
CE.INPUT$state <- "TB-CE"
CE.INPUT <- CE.INPUT %>% rename(locus_tag = X1)
INPUT <- rbind(WT.INPUT, NS.INPUT, CE.INPUT)


WT.TY <- read_table(file = "TY.WT.CE.txt", col_names = F)
WT.TY$state <- "WT-TY-CE"
WT.TY <- WT.TY %>% rename(locus_tag = X1)
NS.TY <- read_table(file = "TY.NS.CE.txt", col_names = F)
NS.TY$state <- "NS-TY-CE"
NS.TY <- NS.TY %>% rename(locus_tag = X1)
CE.TY <- read_table(file = "TY.CE.txt", col_names = F)
CE.TY$state <- "TY-CE"
CE.TY <- CE.TY %>% rename(locus_tag = X1)
TY <- rbind(WT.TY, NS.TY, CE.TY)


WT.GR <- read_table(file = "GR.WT.CE.txt", col_names = F)
WT.GR$state <- "WT-GR-CE"
WT.GR <- WT.GR %>% rename(locus_tag = X1)
NS.GR <- read_table(file = "GR.NS.CE.txt", col_names = F)
NS.GR$state <- "NS-GR-CE"
NS.GR <- NS.GR %>% rename(locus_tag = X1)
CE.GR <- read_table(file = "GR.CE.txt", col_names = F)
CE.GR$state <- "GR-CE"
CE.GR <- CE.GR %>% rename(locus_tag = X1)
GR <- rbind(WT.GR, NS.GR, CE.GR)

CESPLIT <- rbind(INPUT, TY, GR)

CESPLIT <- left_join(CESPLIT, FUNCS)
CESPLIT <- left_join(CESPLIT, ANNOT)
CESPLIT <- left_join(CESPLIT, WTINPUT)
CESPLIT <- left_join(CESPLIT, WTGR)
CESPLIT <- left_join(CESPLIT, WTTY)
CESPLIT <- left_join(CESPLIT, NSINPUT)
CESPLIT <- left_join(CESPLIT, NSGR)
CESPLIT <- left_join(CESPLIT, NSTY)
write_delim(CESPLIT, "Table S2 In Vitro Conditionally Essential Genes of M japonicum R7A.csv", delim = ",", col_names = T, append = F)

#reading in the curated functional categories table
COGTable <- read_csv("Table S2 In Vitro Conditionally Essential Genes of M japonicum R7A.csv", col_names = T, trim_ws = T)
COGTable <- COGTable %>% mutate(FUN = as.factor(FUN))
COGTable <- COGTable %>% group_by(state) %>% select(FUN) %>% table() %>% as.data.frame.array()
COGTable$state <- rownames(COGTable)
write_csv(COGTable, "COG Table.csv", append = F)

################################################################################
############################ Make KEGG Mapping File ############################
################################################################################

# Reading in TraDIS output files
WTINPUT <- read_csv(file = "WT-INPUT.tradis_insertion_plot.tradis_gene_insert_sites.csv.all.csv",col_names = T, trim_ws = T)
WTINPUT <- WTINPUT %>% rename_all(paste0, ".WTTB")
WTINPUT <- WTINPUT %>% rename(locus_tag = locus_tag.WTTB)
NSINPUT <- read_csv(file = "NS-INPUT.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv",col_names = T, trim_ws = T)
NSINPUT <- NSINPUT %>% rename_all(paste0, ".NSTB")
NSINPUT <- NSINPUT %>% rename(locus_tag = locus_tag.NSTB)

WTTY <- read_csv(file = "WTTY.tradis_insertion_plot.tradis_gene_insert_sites.csv.all.csv",col_names = T, trim_ws = T)
WTTY <- WTTY %>% rename_all(paste0, ".WTTY")
WTTY <- WTTY %>% rename(locus_tag = locus_tag.WTTY)
NSTY <- read_csv(file = "NSTY.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv",col_names = T, trim_ws = T)
NSTY <- NSTY %>% rename_all(paste0, ".NSTY")
NSTY <- NSTY %>% rename(locus_tag = locus_tag.NSTY)

WTGR <- read_csv(file = "WTGR.tradis_insertion_plot.tradis_gene_insert_sites.csv.all.csv",col_names = T, trim_ws = T)
WTGR <- WTGR %>% rename_all(paste0, ".WTGR")
WTGR <- WTGR %>% rename(locus_tag = locus_tag.WTGR)
NSGR <- read_csv(file = "NSGR.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv",col_names = T, trim_ws = T)
NSGR <- NSGR %>% rename_all(paste0, ".NSGR")
NSGR <- NSGR %>% rename(locus_tag = locus_tag.NSGR)

WTES <- read_csv(file = "../Essential/WT.tradis_insertion_plot.tradis_gene_insert_sites.csv.all.csv",col_names = T, trim_ws = T)
WTES <- WTES %>% rename_all(paste0, ".WTES")
WTES <- WTES %>% rename(locus_tag = locus_tag.WTES)
NSES <- read_csv(file = "../Essential/NS.tradis_insertion_plot.tradis_gene_insert_sites.ICE.Filt.csv.all.csv",col_names = T, trim_ws = T)
NSES <- NSES %>% rename_all(paste0, ".NSES")
NSES <- NSES %>% rename(locus_tag = locus_tag.NSES)

#reading in the curated functional categories table
CETable <- read_csv("Table S2 In Vitro Conditionally Essential Genes of M japonicum R7A.csv", col_names = T, trim_ws = T)
CETable <- CETable %>% mutate(FUN = as.factor(FUN))
CETable <- CETable %>% select(locus_tag, state, FUN)
ESTable <- read_csv("../Essential/Table S1 Essential Genes of M japonicum R7A.csv", col_names = T, trim_ws = T)
ESTable <- ESTable %>% mutate(FUN = as.factor(FUN))
ESTable <- ESTable %>% select(locus_tag, state, FUN)
TnseqState <- rbind(ESTable, CETable)

tnseqFilter <- TnseqState$locus_tag
`%notin%` <- Negate(`%in%`)
TnseqTable <- WTINPUT %>% select(locus_tag) %>% filter(locus_tag %notin% tnseqFilter)
TnseqTable$state <- "NE"
TnseqTable$FUN <- NA

TnseqTable <- rbind(TnseqTable, TnseqState)

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
#Preparing NCBI Annotations Column
NCBI <- WTINPUT %>% select(locus_tag, fcn.WTTB) %>% transmute(locus_tag, NCBI = fcn.WTTB)
#Merge the annotations
ANNOT <- left_join(NCBI, KEGG)
ANNOT <- left_join(ANNOT, COGS)

#combine into megaTable
TnseqTable <- left_join(x = TnseqTable, y = ANNOT, by = "locus_tag")
TnseqTable <- left_join(TnseqTable, WTES, "locus_tag")
TnseqTable <- left_join(TnseqTable, NSES, "locus_tag")
TnseqTable <- left_join(TnseqTable, WTINPUT, "locus_tag")
TnseqTable <- left_join(TnseqTable, NSINPUT, "locus_tag")
TnseqTable <- left_join(TnseqTable, WTTY, "locus_tag")
TnseqTable <- left_join(TnseqTable, NSTY, "locus_tag")
TnseqTable <- left_join(TnseqTable, WTGR, "locus_tag")
TnseqTable <- left_join(TnseqTable, NSGR, "locus_tag")

write_delim(x = TnseqTable, file = "../../Table S1 Summary of In Vitro Tnseq Analysis.tsv", delim = "\t", col_names = T, append = F)
