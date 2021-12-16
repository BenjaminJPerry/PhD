# TODO: Essentiality Function ####
# tnseq <- function(locus_tags=NULL,NS_tnseq=NULL, WT_tnseq=NULL, NS_TraDIS_thresh="", WT_TraDIS_thresh="", length_perc=95, ES_diviser=2, GD_multiple=2){
#   
#   # Set default stares to NE
#   states$locus_tags <- locus_tags
#   states$states <- "NE"
#   # Layer on the states 
#   
#   # ii in NS and WT below 1/2 the TraDIS thresholds
#   ES$ES_state <- replace(x = ES$ES_state, ES$NS_ii < 0.05 & ES$WT_ii < 0.048, "CORE-ES")
#   # ii in NS less than TraDIS and 2x WT TraDIS thresholds
#   ES$ES_state <- replace(x = ES$ES_state, ES$NS_ii < 0.05 & ES$WT_ii > (0.048*2), "NS-ES")
#   # ii in WT less than TraDIS and 2x NS TraDIS thresholds
#   ES$ES_state <- replace(x = ES$ES_state, ES$NS_ii > (0.05*2) & ES$WT_ii < 0.048, "WT-ES")
#   # 
#   ES$ES_state <- replace(x = ES$ES_state, ES$NS_ii < 0.05 & ES$WT_ii <= (0.048*2) & ES$WT_ii >= 0.048, "GD")
#   
#   ES$ES_state <- replace(x = ES$ES_state, ES$NS_ii <= (0.05*2) & ES$NS_ii >= 0.05 & ES$WT_ii < 0.048 , "GD")
#   
#   ES$ES_state <- replace(x = ES$ES_state, ES$NS_ii <= (0.05*2) & ES$NS_ii >= 0.05 & ES$WT_ii <= 0.001, "WT-ES")
#   
#   ES$ES_state <- replace(x = ES$ES_state, ES$NS_ii <= 0.001 & ES$WT_ii <= (0.048*2) & ES$WT_ii >= 0.048, "NS-ES")
#   
#   return(states)
# }

### Code for chapter 3. In vitro Tn-seq

setwd("~/Desktop/PhD/")
set.seed(seed = 1953)

library(tidyverse)
library(cowplot)
library(ggpubr)

readTradis <- function(wrkingdir = "") {
  cwd<-getwd()
  setwd(dir = wrkingdir)
  tradisFiles <-
    list.files(pattern = ".tradis_gene_insert_sites.csv.all.csv")
  tradisData <-
    data_frame(filename = tradisFiles) %>% mutate(file_contents = map(filename, ~read_delim(., delim = ",")))
  tradisTable <- unnest(tradisData, cols = c(file_contents))
  tradisTable <- tradisTable %>% separate(
    col = filename,
    sep = "[.]",
    remove = F,
    into = c("library", NA, NA, NA, NA, NA)
  )
  tradisTable <-
    tradisTable %>% select(locus_tag, ins_index, library)
  tradisTable <-
    tradisTable %>% pivot_wider(values_from = ins_index, names_from = library, values_fill = NA)
  setwd(cwd)
  return(tradisTable)
}

### 1. Compile insertion index summary tables for analysis. ###

ES <- readTradis(wrkingdir = "~/Desktop/PhD/data/tradis/essential/")
INVITRO <- readTradis(wrkingdir = "~/Desktop/PhD/data/tradis/in vitro/")

### 2. Compile Annotations Table ###
# COGS <- read_tsv(file = "data/ref/R7A2020.eggNOG.tsv", col_names = F, trim_ws = T, comment = "#")
# COGS <- COGS %>% rename(locus_tag = X1, gene = X6, COG = X21, CA = X22)
# COGS <- COGS %>% select(locus_tag, gene, COG, CA)
# KEGG <- read_table(file = "data/ref/R7A2020.KO.KofamKOALA.txt", col_names = F, comment = "#")
# KEGG <- KEGG %>% rename(locus_tag = X2, KO = X3, KA = X7)
# #Filter duplicate KEGG entries for smallest E-value
# KEGG <- KEGG[order(KEGG$locus_tag, abs(KEGG$X6) ), ] ### sort first
# KEGG <- KEGG[ !duplicated(KEGG$locus_tag), ]  ### Keep highest
# KEGG <- KEGG %>%  select(locus_tag, KO, KA)
# ANNOT <- left_join(COGS,KEGG)

### 2.1 Take Mannual Annotations from Old Analysis ###
# NOTE: At this point the annotation table (ANNOT) was mannual curated. A Column "FUN" was created which
#       selected the single best COG category for genes of interest (ES/GD/etc.) but not all neutral genes.
#       Annotations were included for non protein features such as tRNAs as well. I've read in this table below.

ANNOT <- read_csv(file = "data/ref/MLR7A_2020.annotations.csv")
rm(COGS)
rm(KEGG)

### 3. Update the traDIS insertion index summary tables with essentiality states based on the thresholds calculated from the traDIS analysis. ###
###    Essentiality States are for combined data of the WT and the NS datasets
###    These values are printed out into a pdf file generate by the traDIS analysis pipeline. We then apply from rules to these thresholds which 
###    increase the robustnesss of the states. Specicially, we select a decreased threshold for ICE genes due to decreased ins density, and we apply
###    multipliers to the thresholds to identify cuts offs for GD genes.

# Update the ins_index table to have NS ICE genes registered as NA
ES <- left_join(ANNOT, ES)
ES <- ES %>% rename(NS_ii = NS, WT_ii = WT)
ES <- ES %>% mutate(NS_ii = replace(NS_ii, ICESym, NA))

# Create state variable to hold the gene state, default is NE
ES$ES_state <- "NE"

ES$ES_state <- replace(x = ES$ES_state,
                    ES$NS_ii < 0.05 & ES$WT_ii < 0.048,
                    "CORE-ES")
ES$ES_state <- replace(x = ES$ES_state,
                    ES$NS_ii < 0.05 & ES$WT_ii > (0.048*2),
                    "NS-ES")
ES$ES_state <- replace(x = ES$ES_state,
                    ES$NS_ii > (0.05*2) & ES$WT_ii < 0.048,
                    "WT-ES")
ES$ES_state <- replace(x = ES$ES_state,
                    ES$NS_ii < 0.05 & ES$WT_ii <= (0.048*2) & ES$WT_ii >= 0.048,
                    "GD")
ES$ES_state <- replace(x = ES$ES_state,
                    ES$NS_ii <= (0.05*2) & ES$NS_ii >= 0.05 & ES$WT_ii < 0.048 ,
                    "GD")
ES$ES_state <- replace(x = ES$ES_state,
                    ES$NS_ii <= (0.05*2) & ES$NS_ii >= 0.05 & ES$WT_ii <= 0.001,
                    "WT-ES")
ES$ES_state <- replace(x = ES$ES_state,
                    ES$NS_ii <= 0.001 & ES$WT_ii <= (0.048*2) & ES$WT_ii >= 0.048,
                    "NS-ES")
# Plot the correlated insertion indices
ES_plot <- ES %>% ggplot(aes(x=NS_ii, y=WT_ii, color=ES_state, alpha=0.1)) + 
  geom_point() + 
  theme_minimal(base_size = 20) + 
  guides(colour = guide_legend(override.aes = list(size=8))) +
  geom_hline(yintercept = 0.048) +
  geom_vline(xintercept = 0.05) + xlim(0,0.3) + ylim(0,0.3)
ggsave(plot = ES_plot, "results/chapter 3 in vitro/ES.ins_index.png", width = 12, height = 8)

ES_plot_zoom <- ES %>% ggplot(aes(x=NS_ii, y=WT_ii, color=ES_state, alpha=0.1)) + 
  geom_point() + 
  theme_minimal(base_size = 20) + 
  guides(colour = guide_legend(override.aes = list(size=8))) +
  geom_hline(yintercept = 0.048) +
  geom_vline(xintercept = 0.05) + xlim(0,0.1) + ylim(0,0.1)
ggsave(plot = ES_plot_zoom, "results/chapter 3 in vitro/ES.ins_index.zoom.png", width = 12, height = 8)

### 4. Update the traDIS insertion index summary tables for the conditionally essential states using similar rules to section 3. ###

# Update the insertion index table to have NS ICE genes registered as NA
INVITRO <- left_join(ANNOT, INVITRO)
INVITRO <- INVITRO %>% rename(NSTB_ii = `NS-INPUT`, NSGR_ii = NSGR, NSTY_ii = NSTY,
                              WTTB_ii = `WT-INPUT`, WTGR_ii = WTGR, WTTY_ii = WTTY)
INVITRO <- INVITRO %>% mutate(NSTB_ii = replace(NSTB_ii, ICESym, NA),NSGR_ii = replace(NSGR_ii, ICESym, NA),NSTY_ii = replace(NSTY_ii, ICESym, NA))

###    TY States - TY_states
INVITRO$TY_state <- "NE"
INVITRO$TY_state <- replace(x = INVITRO$TY_state,
                    INVITRO$NSTY_ii < 0.016 & INVITRO$WTTY_ii < 0.018,
                    "TY-ES")
INVITRO$TY_state <- replace(x = INVITRO$TY_state,
                    INVITRO$NSTY_ii < 0.016 & INVITRO$WTTY_ii > (0.018*2),
                    "NS-TY-ES")
INVITRO$TY_state <- replace(x = INVITRO$TY_state,
                    INVITRO$NSTY_ii > (0.016*2) & INVITRO$WTTY_ii < 0.018,
                    "WT-TY-ES")

INVITRO$TY_state <- replace(x = INVITRO$TY_state,
                    INVITRO$NSTY_ii < 0.016 & INVITRO$WTTY_ii <= (0.018*2) & INVITRO$WTTY_ii >= 0.018,
                    "TY-GD")
INVITRO$TY_state <- replace(x = INVITRO$TY_state,
                    INVITRO$NSTY_ii <= (0.016*2) & INVITRO$NSTY_ii >= 0.016 & INVITRO$WTTY_ii < 0.018,
                    "TY-GD")

INVITRO$TY_state <- replace(x = INVITRO$TY_state,
                    INVITRO$NSTY_ii <= 0.001 & INVITRO$WTTY_ii <= (0.018*2) & INVITRO$WTTY_ii >= 0.018,
                    "NS-TY-ES")
INVITRO$TY_state <- replace(x = INVITRO$TY_state,
                    INVITRO$NSTY_ii <= (0.016*2) & INVITRO$NSTY_ii >= 0.016 & INVITRO$WTTY_ii <= 0.001,
                    "WT-TY-ES")

INVITRO <- INVITRO %>% mutate(TY_state = as.factor(TY_state))
INVITRO$TY_state <- factor(INVITRO$TY_state, levels = c("TY-ES", "TY-GD", "NE", "NS-TY-ES", "WT-TY-ES"))

TY_plot <- INVITRO %>% ggplot(aes(x=NSTY_ii, y=WTTY_ii, color=TY_state, alpha=0.1)) + 
  geom_point() + 
  theme_minimal(base_size = 20) + 
  guides(colour = guide_legend(override.aes = list(size=8))) +
  geom_hline(yintercept = 0.018) +
  geom_vline(xintercept = 0.016) + xlim(0,0.1) + ylim(0,0.1)

ggsave(plot = TY_plot, "results/chapter 3 in vitro/TY.ins_index.png", width = 12, height = 8)

###    GR States - GR_states
INVITRO$GR_state <- "NE"
INVITRO$GR_state <- replace(x = INVITRO$GR_state,
                    INVITRO$NSGR_ii < 0.015 & INVITRO$WTGR_ii < 0.012,
                    "GR-ES")
INVITRO$GR_state <- replace(x = INVITRO$GR_state,
                    INVITRO$NSGR_ii < 0.015 & INVITRO$WTGR_ii > (0.012*2),
                    "NS-GR-ES")
INVITRO$GR_state <- replace(x = INVITRO$GR_state,
                    INVITRO$NSGR_ii > (0.015*2) & INVITRO$WTGR_ii < 0.012,
                    "WT-GR-ES")

INVITRO$GR_state <- replace(x = INVITRO$GR_state,
                    INVITRO$NSGR_ii < 0.015 & INVITRO$WTGR_ii <= (0.012*2) & INVITRO$WTGR_ii >= 0.012,
                    "GR-GD")
INVITRO$GR_state <- replace(x = INVITRO$GR_state,
                    INVITRO$NSGR_ii <= (0.015*2) & INVITRO$NSGR_ii >= 0.015 & INVITRO$WTGR_ii < 0.012,
                    "GR-GD")
INVITRO$GR_state <- replace(x = INVITRO$GR_state,
                    INVITRO$NSGR_ii <= 0.001 & INVITRO$WTGR_ii <= (0.012*2) & INVITRO$WTGR_ii >= 0.012,
                    "NS-GR-ES")
INVITRO$GR_state <- replace(x = INVITRO$GR_state,
                    INVITRO$NSGR_ii <= (0.015*2) & INVITRO$NSGR_ii >= 0.015 & INVITRO$WTGR_ii <= 0.001,
                    "WT-GR-ES")

INVITRO <- INVITRO %>% mutate(GR_state = as.factor(GR_state))
INVITRO$GR_state <- factor(INVITRO$GR_state, levels = c("GR-ES", "GR-GD", "NE", "NS-GR-ES", "WT-GR-ES"))

GR_plot <- INVITRO %>% ggplot(aes(x=NSGR_ii, y=WTGR_ii, color=GR_state, alpha=0.1)) + 
  geom_point() + 
  theme_minimal(base_size = 20) + 
  guides(colour = guide_legend(override.aes = list(size=8))) +
  geom_hline(yintercept = 0.012) +
  geom_vline(xintercept = 0.015) + xlim(0,0.1) + ylim(0,0.1)

ggsave(plot = GR_plot, "results/chapter 3 in vitro/GR.ins_index.png", width = 12, height = 8)

###    TB States - TB_states
INVITRO$TB_state <- "NE"
INVITRO$TB_state <- replace(x = INVITRO$TB_state,
                    INVITRO$NSTB_ii < 0.017 & INVITRO$WTTB_ii < 0.015,
                    "TB-ES")
INVITRO$TB_state <- replace(x = INVITRO$TB_state,
                    INVITRO$NSTB_ii < 0.017 & INVITRO$WTTB_ii > (0.015*2),
                    "NS-TB-ES")
INVITRO$TB_state <- replace(x = INVITRO$TB_state,
                    INVITRO$NSTB_ii > (0.017*2) & INVITRO$WTTB_ii < 0.015,
                    "WT-TB-ES")

INVITRO$TB_state <- replace(x = INVITRO$TB_state,
                    INVITRO$NSTB_ii < 0.017 & INVITRO$WTTB_ii <= (0.015*2) & INVITRO$WTTB_ii >= 0.015,
                    "TB-GD")
INVITRO$TB_state <- replace(x = INVITRO$TB_state,
                    INVITRO$NSTB_ii <= (0.017*2) & INVITRO$NSTB_ii >= 0.017 & INVITRO$WTTB_ii < 0.015,
                    "TB-GD")

INVITRO$TB_state <- replace(x = INVITRO$TB_state,
                    INVITRO$NSTB_ii <= 0.001 & INVITRO$WTTB_ii <= (0.015*2) & INVITRO$WTTB_ii >= 0.015,
                    "NS-TB-ES")
INVITRO$TB_state <- replace(x = INVITRO$TB_state,
                    INVITRO$NSTB_ii <= (0.017*2) & INVITRO$NSTB_ii >= 0.017 & INVITRO$WTTB_ii <= 0.001,
                    "WT-TB-ES")

INVITRO <- INVITRO %>% mutate(TB_state = as.factor(TB_state))
INVITRO$TB_state <- factor(INVITRO$TB_state, levels = c("TB-ES", "TB-GD", "NE", "NS-TB-ES", "WT-TB-ES"))

TB_plot <- INVITRO %>% ggplot(aes(x=NSTB_ii, y=WTTB_ii, color=TB_state, alpha=0.1)) + 
  geom_point() + 
  theme_minimal(base_size = 20) + 
  guides(colour = guide_legend(override.aes = list(size=8))) +
  geom_hline(yintercept = 0.015) +
  geom_vline(xintercept = 0.017) + xlim(0,0.1) + ylim(0,0.1)

ggsave(plot = TB_plot, "results/chapter 3 in vitro/TB.ins_index.png", width = 12, height = 8)


### Merge the ES and INVITRO dataframes into CHPT3 dataframe ###
CHPT3 <- ES %>% select(locus_tag, NS_ii, WT_ii, ES_state)
CHPT3 <- left_join(INVITRO, CHPT3, "locus_tag")

### 5. Update the essentiality states of ICE genes using reduced ii thresholds

# ICE Essentiality States

CHPT3$ICE_state <- "NE"

# TODO: Update the code ####

CHPT3$ICE_state <- replace(x = CHPT3$ICE_state, !CHPT3$ICESym, NA)

ICE$state <- replace(x = ICE$state, ICE$ES <= 0.03, "ICE-ES")
ICE$state <- replace(x = ICE$state, ICE$ES > 0.03 & ICE$GR == 0, "ICE-GR-ES")
ICE$state <- replace(x = ICE$state, ICE$ES > 0.03 & ICE$TY == 0, "ICE-TY-ES")
ICE$state <- replace(x = ICE$state, ICE$ES > 0.03 & ICE$GR == 0 & ICE$TY == 0, "ICE-AGAR-ES")
ICE$state <- replace(x = ICE$state, ICE$ES > 0.03 & ICE$TB == 0, "ICE-TB-ES")

write_tsv(ICE, "../ICE in vitro updated states.txt")

ICE <- ICE %>% mutate(state = as.factor(state))
ICE$state <- factor(ICE$state, levels = c("ICE-ES", "ICE-AGAR-ES", "NE", "ICE-TB-ES", "ICE-GR-ES"))

ICE %>% ggplot(aes(x=ES, y=TY, alpha=0.1, color=state)) + 
  geom_point() + 
  theme_minimal(base_size = 20) + 
  guides(colour = guide_legend(override.aes = list(size=8))) +
  geom_vline(xintercept = 0.03) +
  ylim(0,0.125) #+
ggsave("ICE.TY.ins_index.png", width = 12, height = 8) +
  ggsave("ICE.TY.ins_index.svg", width = 12, height = 8)

ICE %>% ggplot(aes(x=ES, y=TB, alpha=0.1, color=state)) + 
  geom_point() + 
  theme_minimal(base_size = 20) + 
  guides(colour = guide_legend(override.aes = list(size=8))) +
  geom_vline(xintercept = 0.03) +
  ylim(0,0.125) #+
ggsave("ICE.TB.ins_index.png", width = 12, height = 8) +
  ggsave("ICE.TB.ins_index.svg", width = 12, height = 8)

ICE %>% ggplot(aes(x=ES, y=GR, alpha=0.1, color=state)) + 
  geom_point() + 
  theme_minimal(base_size = 20) + 
  guides(colour = guide_legend(override.aes = list(size=8))) +
  geom_vline(xintercept = 0.03) +
  ylim(0,0.125) #+
ggsave("ICE.GR.ins_index.png", width = 12, height = 8) +
  ggsave("ICE.GR.ins_index.svg", width = 12, height = 8)

### Summary Plots ###
invitro_plot <- ggarrange(TB_plot, GR_plot, TY_plot, nrow = 1, align = "h", legend = "none")
invitro_plot
ggsave(plot = invitro_plot, "results/chapter 3 in vitro/invitro.ins_index.png", bg = "white")

essential_plot <- ggarrange(ES_plot, ES_plot_zoom, align = "h", legend = "none")
essential_plot
ggsave(plot = essential_plot, "results/chapter 3 in vitro/essential.ins_index.png", bg = "white")

