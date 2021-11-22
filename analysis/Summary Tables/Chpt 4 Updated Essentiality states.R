library(tidyverse)
setwd(dir = "~/Desktop/Chapter 3 In Vitro/In Vitro Tnseq/")

ES <- read_tsv(file = "ES.ins_indx.tsv", trim_ws = T)
ES$state <- "NE"
ES$state <- replace(x = ES$state,
                    ES$ins_index.NSES < 0.05 & ES$ins_index.WTES < 0.048,
                    "CORE-ES")
ES$state <- replace(x = ES$state,
                    ES$ins_index.NSES < 0.05 & ES$ins_index.WTES > (0.048*2),
                    "NS-ES")
ES$state <- replace(x = ES$state,
                    ES$ins_index.NSES > (0.05*2) & ES$ins_index.WTES < 0.048,
                    "WT-ES")
ES$state <- replace(x = ES$state,
                    ES$ins_index.NSES < 0.05 & ES$ins_index.WTES <= (0.048*2) & ES$ins_index.WTES >= 0.048,
                    "GD")
ES$state <- replace(x = ES$state,
                    ES$ins_index.NSES <= (0.05*2) & ES$ins_index.NSES >= 0.05 & ES$ins_index.WTES < 0.048 ,
                    "GD")
ES$state <- replace(x = ES$state,
                    ES$ins_index.NSES <= (0.05*2) & ES$ins_index.NSES >= 0.05 & ES$ins_index.WTES <= 0.001,
                    "WT-ES")
ES$state <- replace(x = ES$state,
                    ES$ins_index.NSES <= 0.001 & ES$ins_index.WTES <= (0.048*2) & ES$ins_index.WTES >= 0.048,
                    "NS-ES")

ES %>% ggplot(aes(x=ins_index.NSES, y=ins_index.WTES, color=state, alpha=0.1)) + 
        geom_point() + 
        theme_minimal(base_size = 20) + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_hline(yintercept = 0.048) +
        geom_vline(xintercept = 0.05) +
        ggsave("ES.ins_index.png", width = 12, height = 8)

ES %>% ggplot(aes(x=ins_index.NSES, y=ins_index.WTES, color=state, alpha=0.1)) + 
        geom_point() + 
        theme_minimal(base_size = 20) + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_hline(yintercept = 0.048) +
        geom_vline(xintercept = 0.05) + xlim(0,0.175) + ylim(0,0.175) +
        ggsave("ES.ins_index.zoom.png", width = 12, height = 8)


TY <- read_tsv(file = "TY.ins_indx.tsv", trim_ws = T)
TY$state <- "NE"
TY$state <- replace(x = TY$state,
                    TY$ins_index.NSTY < 0.016 & TY$ins_index.WTTY < 0.018,
                    "TY-ES")
TY$state <- replace(x = TY$state,
                    TY$ins_index.NSTY < 0.016 & TY$ins_index.WTTY > (0.018*2),
                    "NS-TY-ES")
TY$state <- replace(x = TY$state,
                    TY$ins_index.NSTY > (0.016*2) & TY$ins_index.WTTY < 0.018,
                    "WT-TY-ES")

TY$state <- replace(x = TY$state,
                    TY$ins_index.NSTY < 0.016 & TY$ins_index.WTTY <= (0.018*2) & TY$ins_index.WTTY >= 0.018,
                    "TY-GD")
TY$state <- replace(x = TY$state,
                    TY$ins_index.NSTY <= (0.016*2) & TY$ins_index.NSTY >= 0.016 & TY$ins_index.WTTY < 0.018,
                    "TY-GD")

TY$state <- replace(x = TY$state,
                    TY$ins_index.NSTY <= 0.001 & TY$ins_index.WTTY <= (0.018*2) & TY$ins_index.WTTY >= 0.018,
                    "NS-TY-ES")
TY$state <- replace(x = TY$state,
                    TY$ins_index.NSTY <= (0.016*2) & TY$ins_index.NSTY >= 0.016 & TY$ins_index.WTTY <= 0.001,
                    "WT-TY-ES")


TY <- TY %>% mutate(state = as.factor(state))
TY$state <- factor(TY$state, levels = c("TY-ES", "TY-GD", "NE", "NS-TY-ES", "WT-TY-ES"))

TY %>% ggplot(aes(x=ins_index.NSTY, y=ins_index.WTTY, color=state, alpha=0.1)) + 
        geom_point() + 
        theme_minimal(base_size = 20) + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_hline(yintercept = 0.018) +
        geom_vline(xintercept = 0.016) + xlim(0,0.175) + ylim(0,0.175) +
        ggsave("TY.ins_index.png", width = 12, height = 8)


TB <- read_tsv(file = "TB.ins_indx.tsv", trim_ws = T)
TB$state <- "NE"
TB$state <- replace(x = TB$state,
                    TB$ins_index.NSTB < 0.017 & TB$ins_index.WTTB < 0.015,
                    "TB-ES")
TB$state <- replace(x = TB$state,
                    TB$ins_index.NSTB < 0.017 & TB$ins_index.WTTB > (0.015*2),
                    "NS-TB-ES")
TB$state <- replace(x = TB$state,
                    TB$ins_index.NSTB > (0.017*2) & TB$ins_index.WTTB < 0.015,
                    "WT-TB-ES")

TB$state <- replace(x = TB$state,
                    TB$ins_index.NSTB < 0.017 & TB$ins_index.WTTB <= (0.015*2) & TB$ins_index.WTTB >= 0.015,
                    "TB-GD")
TB$state <- replace(x = TB$state,
                    TB$ins_index.NSTB <= (0.017*2) & TB$ins_index.NSTB >= 0.017 & TB$ins_index.WTTB < 0.015,
                    "TB-GD")

TB$state <- replace(x = TB$state,
                    TB$ins_index.NSTB <= 0.001 & TB$ins_index.WTTB <= (0.015*2) & TB$ins_index.WTTB >= 0.015,
                    "NS-TB-ES")
TB$state <- replace(x = TB$state,
                    TB$ins_index.NSTB <= (0.017*2) & TB$ins_index.NSTB >= 0.017 & TB$ins_index.WTTB <= 0.001,
                    "WT-TB-ES")



TB <- TB %>% mutate(state = as.factor(state))
TB$state <- factor(TB$state, levels = c("TB-ES", "TB-GD", "NE", "NS-TB-ES", "WT-TB-ES"))

TB %>% ggplot(aes(x=ins_index.NSTB, y=ins_index.WTTB, color=state, alpha=0.1)) + 
        geom_point() + 
        theme_minimal(base_size = 20) + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_hline(yintercept = 0.015) +
        geom_vline(xintercept = 0.017) + xlim(0,0.175) + ylim(0,0.175) +
        ggsave("TB.ins_index.png", width = 12, height = 8)


GR <- read_tsv(file = "GR.ins_indx.tsv", trim_ws = T)
GR$state <- "NE"
GR$state <- replace(x = GR$state,
                    GR$ins_index.NSGR < 0.015 & GR$ins_index.WTGR < 0.012,
                    "GR-ES")
GR$state <- replace(x = GR$state,
                    GR$ins_index.NSGR < 0.015 & GR$ins_index.WTGR > (0.012*2),
                    "NS-GR-ES")
GR$state <- replace(x = GR$state,
                    GR$ins_index.NSGR > (0.015*2) & GR$ins_index.WTGR < 0.012,
                    "WT-GR-ES")

GR$state <- replace(x = GR$state,
                    GR$ins_index.NSGR < 0.015 & GR$ins_index.WTGR <= (0.012*2) & GR$ins_index.WTGR >= 0.012,
                    "GR-GD")
GR$state <- replace(x = GR$state,
                    GR$ins_index.NSGR <= (0.015*2) & GR$ins_index.NSGR >= 0.015 & GR$ins_index.WTGR < 0.012,
                    "GR-GD")
GR$state <- replace(x = GR$state,
                    GR$ins_index.NSGR <= 0.001 & GR$ins_index.WTGR <= (0.012*2) & GR$ins_index.WTGR >= 0.012,
                    "NS-GR-ES")
GR$state <- replace(x = GR$state,
                    GR$ins_index.NSGR <= (0.015*2) & GR$ins_index.NSGR >= 0.015 & GR$ins_index.WTGR <= 0.001,
                    "WT-GR-ES")

GR %>% ggplot(aes(x=ins_index.NSGR, y=ins_index.WTGR, color=state, alpha=0.1)) + 
        geom_point() + 
        theme_minimal(base_size = 20) + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_hline(yintercept = 0.012) +
        geom_vline(xintercept = 0.015) + xlim(0,0.175) + ylim(0,0.175) +
        ggsave("GR.ins_index.png", width = 12, height = 8)

ES <- ES %>% rename(ES.state = state)
TB <- TB %>% rename(TB.state = state)
TY <- TY %>% rename(TY.state = state)
GR <- GR %>% rename(GR.state = state)

ES <- tibble(ES)
TB <- tibble(TB)
TY <- tibble(TY)
GR <- tibble(GR)

write_tsv(ES, "ES.ins_index.states.tsv")
write_tsv(TB, "TB.ins_index.states.tsv")
write_tsv(TY, "TY.ins_index.states.tsv")
write_tsv(GR, "GR.ins_index.states.tsv")

### ICE Buisness, reducing the threshold to compensate for the decreased insertion density
### ICE Essentiality States
hist(ICE$ES, breaks = 25, plot = T)
hist(ICE$TY, breaks = 50, plot = F)
hist(ICE$TB, breaks = 50, plot = T)
hist(ICE$GR, breaks = 50, plot = F)

ICE <- read_tsv(file = "ICE.ins_indx.tsv")
ICE$state <- "NE"
ICE$state <- replace(x = ICE$state,
                    ICE$ES <= 0.03,
                    "ICE-ES")
ICE$state <- replace(x = ICE$state,
                     ICE$ES > 0.03 & ICE$GR == 0,
                     "ICE-GR-ES")
ICE$state <- replace(x = ICE$state,
                     ICE$ES > 0.03 & ICE$TY == 0,
                     "ICE-TY-ES")
ICE$state <- replace(x = ICE$state,
                     ICE$ES > 0.03 & ICE$GR == 0 & ICE$TY ==0,
                     "ICE-AGAR-ES")
ICE$state <- replace(x = ICE$state,
                     ICE$ES > 0.03 & ICE$TB == 0,
                     "ICE-TB-ES")

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


# Chapter 4 - Rhizosphere ES genes
rhizo <- read_csv(
        file = "~/Desktop/Chapter 4 Rhizosphere/Rhizosphere Tnseq/RHIZO.ins_indices.csv")

rhizo$state <- "NE"
rhizo$state <- replace(x = rhizo$state, rhizo$NS.RHIZO < 0.033 & rhizo$WT.RHIZO < 0.033, "RHIZO-ES")
rhizo$state <- replace(x = rhizo$state, rhizo$NS.INPUT < 0.017 & rhizo$WT.INPUT < 0.015, "INPUT-ES")
rhizo$state <- replace(x = rhizo$state, rhizo$NS.INPUT < 0.017 & rhizo$WT.INPUT >= 0.015, "INPUT-GD")
rhizo$state <- replace(x = rhizo$state, rhizo$NS.INPUT >= 0.017 & rhizo$WT.INPUT < 0.015, "INPUT-GD")
rhizo$state <- replace(x = rhizo$state, rhizo$NS.INPUT <= 0.001 & rhizo$WT.INPUT >= 0.015,"INPUT-ES")
rhizo$state <- replace(x = rhizo$state, rhizo$NS.INPUT >= 0.017 & rhizo$WT.INPUT <= 0.001, "INPUT-ES")
rhizo$state <- replace(x = rhizo$state, rhizo$state == "NE" & rhizo$NS.RHIZO >= 0.033 & rhizo$WT.RHIZO <= 0.033, "RHIZO-GD")
rhizo$state <- replace(x = rhizo$state, rhizo$state == "NE" & rhizo$NS.RHIZO <= 0.033 & rhizo$WT.RHIZO >= 0.033, "RHIZO-GD")
rhizo$state <- replace(x = rhizo$state, rhizo$state == "INPUT-GD" & rhizo$NS.RHIZO >= 0.033 & rhizo$WT.RHIZO >= 0.033, "INPUT-GD-RR")
rhizo$state <- replace(x = rhizo$state, rhizo$state == "INPUT-ES" & rhizo$NS.RHIZO >= 0.033 & rhizo$WT.RHIZO >= 0.033, "INPUT-ES-RR")
rhizo$state <- replace(x = rhizo$state, rhizo$state == "INPUT-GD" & rhizo$NS.RHIZO < 0.033 & rhizo$WT.RHIZO < 0.033, "INPUT-GD-RE")

### RHIZO and Background States
rhizo <- rhizo %>% mutate(state = as.factor(state))
rhizo$state <- factor(rhizo$state, levels = c("INPUT-ES",
                                              "INPUT-GD",
                                              "NE",
                                              "RHIZO-ES",
                                              "RHIZO-GD",
                                              "INPUT-GD-RR",
                                              "INPUT-ES-RR",
                                              "INPUT-GD-RE"))

rhizo %>% ggplot(aes(y = WT.INPUT, x = NS.INPUT, alpha = 0.05, color = state)) +
        geom_point() +
        theme_minimal() + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_vline(xintercept = 0.017) + 
        geom_hline(yintercept = 0.015) + xlim(0,0.175) + ylim(0,0.175) #+
        ggsave("Desktop/Chapter 4 Rhizosphere/RHIZO.INPUT.ins_index.png", width = 12, height = 8) +
        ggsave("Desktop/Chapter 4 Rhizosphere/RHIZO.INPUT.ins_index.svg", width = 12, height = 8)


rhizo %>% ggplot(aes(y = WT.RHIZO, x = NS.RHIZO, alpha = 0.05, color = state)) +
        geom_point() +
        theme_minimal() + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_vline(xintercept = 0.033) + 
        geom_hline(yintercept = 0.033) + xlim(0,0.25) + ylim(0,0.25) #+
        ggsave("../../Chapter 4 Rhizosphere/RHIZO.RHIZO.ins_index.png", width = 12, height = 8) +
        ggsave("../../Chapter 4 Rhizosphere/RHIZO.RHIZO.ins_index.svg", width = 12, height = 8)

### Plant specific rhizosphere states
rhizo$GIFU.state <- "NE"
rhizo$GIFU.state <- replace(x = rhizo$GIFU.state, rhizo$NS.RHIZO < 0.033 & rhizo$WT.RHIZO < 0.033, "RHIZO-ES")
rhizo$GIFU.state <- replace(x = rhizo$GIFU.state, rhizo$NS.INPUT < 0.017 & rhizo$WT.INPUT < 0.015, "INPUT-ES")
rhizo$GIFU.state <- replace(x = rhizo$GIFU.state, rhizo$NS.INPUT < 0.017 & rhizo$WT.INPUT >= 0.015, "INPUT-GD")
rhizo$GIFU.state <- replace(x = rhizo$GIFU.state, rhizo$NS.INPUT >= 0.017 & rhizo$WT.INPUT < 0.015, "INPUT-GD")
rhizo$GIFU.state <- replace(x = rhizo$GIFU.state, rhizo$NS.INPUT <= 0.001 & rhizo$WT.INPUT >= 0.015,"INPUT-ES")
rhizo$GIFU.state <- replace(x = rhizo$GIFU.state, rhizo$NS.INPUT >= 0.017 & rhizo$WT.INPUT <= 0.001, "INPUT-ES")
rhizo$GIFU.state <- replace(x = rhizo$GIFU.state, rhizo$GIFU.state == "NE" & rhizo$NS.RHIZO >= 0.033 & rhizo$WT.RHIZO <= 0.033, "RHIZO-GD")
rhizo$GIFU.state <- replace(x = rhizo$GIFU.state, rhizo$GIFU.state == "NE" & rhizo$NS.RHIZO <= 0.033 & rhizo$WT.RHIZO >= 0.033, "RHIZO-GD")


rhizo$GIFU.state <- replace(x = rhizo$GIFU.state, rhizo$state == "NE" & rhizo$NS.GIFU >= 0.013 & rhizo$WT.GIFU <= 0.013, "GIFU-GD")
rhizo$GIFU.state <- replace(x = rhizo$GIFU.state, rhizo$state == "NE" & rhizo$NS.GIFU <= 0.013 & rhizo$WT.GIFU >= 0.013, "GIFU-GD")
rhizo$GIFU.state <- replace(x = rhizo$GIFU.state, rhizo$state == "NE" & rhizo$NS.GIFU < 0.013 & rhizo$WT.GIFU < 0.013, "GIFU-ES")

rhizo <- rhizo %>% mutate(GIFU.state = as.factor(GIFU.state))
rhizo$GIFU.state <- factor(rhizo$GIFU.state, levels = c("RHIZO-ES",
                                                    "RHIZO-GD",
                                                    "NE",
                                                    "INPUT-ES",
                                                    "INPUT-GD",
                                                    "GIFU-ES",
                                                    "GIFU-GD"))

rhizo %>% ggplot(aes(y = WT.GIFU, x = NS.GIFU, alpha = 0.05, color = GIFU.state)) +
        geom_point() +
        theme_minimal() + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_vline(xintercept = 0.013) + 
        geom_hline(yintercept = 0.013) + xlim(0,0.15) + ylim(0,0.16) +
        ggsave("~/Desktop/Chapter 4 Rhizosphere/RHIZO.GIFU.ins_index.png", width = 12, height = 8) +
        ggsave("~/Desktop/Chapter 4 Rhizosphere/RHIZO.GIFU.ins_index.svg", width = 12, height = 8)

rhizo$NFR5.state <- "NE"
rhizo$NFR5.state <- replace(x = rhizo$NFR5.state, rhizo$NS.RHIZO < 0.033 & rhizo$WT.RHIZO < 0.033, "RHIZO-ES")
rhizo$NFR5.state <- replace(x = rhizo$NFR5.state, rhizo$NS.INPUT < 0.017 & rhizo$WT.INPUT < 0.015, "INPUT-ES")
rhizo$NFR5.state <- replace(x = rhizo$NFR5.state, rhizo$NS.INPUT < 0.017 & rhizo$WT.INPUT >= 0.015, "INPUT-GD")
rhizo$NFR5.state <- replace(x = rhizo$NFR5.state, rhizo$NS.INPUT >= 0.017 & rhizo$WT.INPUT < 0.015, "INPUT-GD")
rhizo$NFR5.state <- replace(x = rhizo$NFR5.state, rhizo$NS.INPUT <= 0.001 & rhizo$WT.INPUT >= 0.015,"INPUT-ES")
rhizo$NFR5.state <- replace(x = rhizo$NFR5.state, rhizo$NS.INPUT >= 0.017 & rhizo$WT.INPUT <= 0.001, "INPUT-ES")
rhizo$NFR5.state <- replace(x = rhizo$NFR5.state, rhizo$NFR5.state == "NE" & rhizo$NS.RHIZO >= 0.033 & rhizo$WT.RHIZO <= 0.033, "RHIZO-GD")
rhizo$NFR5.state <- replace(x = rhizo$NFR5.state, rhizo$NFR5.state == "NE" & rhizo$NS.RHIZO <= 0.033 & rhizo$WT.RHIZO >= 0.033, "RHIZO-GD")


rhizo$NFR5.state <- replace(x = rhizo$NFR5.state, rhizo$state == "NE" & rhizo$NS.NFR5 >= 0.013 & rhizo$WT.NFR5 <= 0.010, "NFR5-GD")
rhizo$NFR5.state <- replace(x = rhizo$NFR5.state, rhizo$state == "NE" & rhizo$NS.NFR5 <= 0.013 & rhizo$WT.NFR5 >= 0.010, "NFR5-GD")
rhizo$NFR5.state <- replace(x = rhizo$NFR5.state, rhizo$state == "NE" & rhizo$NS.NFR5 < 0.013 & rhizo$WT.NFR5 < 0.010, "NFR5-ES")

rhizo <- rhizo %>% mutate(NFR5.state = as.factor(NFR5.state))
rhizo$NFR5.state <- factor(rhizo$NFR5.state, levels = c("RHIZO-ES",
                                                        "RHIZO-GD",
                                                        "NE",
                                                        "INPUT-ES",
                                                        "INPUT-GD",
                                                        "NFR5-ES",
                                                        "NFR5-GD"))

rhizo %>% ggplot(aes(y = WT.NFR5, x = NS.NFR5, alpha = 0.05, color = NFR5.state)) +
        geom_point() +
        theme_minimal() + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_vline(xintercept = 0.013) + 
        geom_hline(yintercept = 0.010) + xlim(0,0.15) + ylim(0,0.16) +
        ggsave("~/Desktop/Chapter 4 Rhizosphere/RHIZO.NFR5.ins_index.png", width = 12, height = 8) +
        ggsave("~/Desktop/Chapter 4 Rhizosphere/RHIZO.NFR5.ins_index.svg", width = 12, height = 8)


rhizo$COL0.state <- "NE"
rhizo$COL0.state <- replace(x = rhizo$COL0.state, rhizo$NS.RHIZO < 0.033 & rhizo$WT.RHIZO < 0.033, "RHIZO-ES")
rhizo$COL0.state <- replace(x = rhizo$COL0.state, rhizo$NS.INPUT < 0.017 & rhizo$WT.INPUT < 0.015, "INPUT-ES")
rhizo$COL0.state <- replace(x = rhizo$COL0.state, rhizo$NS.INPUT < 0.017 & rhizo$WT.INPUT >= 0.015, "INPUT-GD")
rhizo$COL0.state <- replace(x = rhizo$COL0.state, rhizo$NS.INPUT >= 0.017 & rhizo$WT.INPUT < 0.015, "INPUT-GD")
rhizo$COL0.state <- replace(x = rhizo$COL0.state, rhizo$NS.INPUT <= 0.001 & rhizo$WT.INPUT >= 0.015,"INPUT-ES")
rhizo$COL0.state <- replace(x = rhizo$COL0.state, rhizo$NS.INPUT >= 0.017 & rhizo$WT.INPUT <= 0.001, "INPUT-ES")
rhizo$COL0.state <- replace(x = rhizo$COL0.state, rhizo$COL0.state == "NE" & rhizo$NS.RHIZO >= 0.033 & rhizo$WT.RHIZO <= 0.033, "RHIZO-GD")
rhizo$COL0.state <- replace(x = rhizo$COL0.state, rhizo$COL0.state == "NE" & rhizo$NS.RHIZO <= 0.033 & rhizo$WT.RHIZO >= 0.033, "RHIZO-GD")

rhizo$COL0.state <- replace(x = rhizo$COL0.state, rhizo$state == "NE" & rhizo$NS.COL0 >= 0.009 & rhizo$WT.COL0 <= 0.010, "COL0-GD")
rhizo$COL0.state <- replace(x = rhizo$COL0.state, rhizo$state == "NE" & rhizo$NS.COL0 <= 0.009 & rhizo$WT.COL0 >= 0.010, "COL0-GD")
rhizo$COL0.state <- replace(x = rhizo$COL0.state, rhizo$state == "NE" & rhizo$NS.COL0 < 0.009 & rhizo$WT.COL0 < 0.010, "COL0-ES")

rhizo <- rhizo %>% mutate(COL0.state = as.factor(COL0.state))
rhizo$COL0.state <- factor(rhizo$COL0.state, levels = c("RHIZO-ES",
                                                        "RHIZO-GD",
                                                        "NE",
                                                        "INPUT-ES",
                                                        "INPUT-GD",
                                                        "COL0-ES",
                                                        "COL0-GD"))

rhizo %>% ggplot(aes(y = WT.COL0, x = NS.COL0, alpha = 0.05, color = COL0.state)) +
        geom_point() +
        theme_minimal() + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_vline(xintercept = 0.009) + 
        geom_hline(yintercept = 0.010) + xlim(0,0.15) + ylim(0,0.16) +
        ggsave("~/Desktop/Chapter 4 Rhizosphere/RHIZO.COL0.ins_index.png", width = 12, height = 8) +
        ggsave("~/Desktop/Chapter 4 Rhizosphere/RHIZO.COL0.ins_index.svg", width = 12, height = 8)

write_tsv(rhizo, "~/Desktop/Chapter 4 Rhizosphere/Rhizosphere updated states.txt")


### ICESym States
ICErhizo <- read_csv(file = "../../Chapter 4 Rhizosphere/Rhizosphere Tnseq/ICE.RHIZO.ins_indices.csv")

ICErhizo$state <- replace(x = ICErhizo$state,
                          ICErhizo$ES <= 0.03,
                          "ICE-INPUT-ES")
ICErhizo$state <- replace(x = ICErhizo$state,
                          ICErhizo$INPUT <= 0.008 & ICErhizo$ES > 0.03,
                          "ICE-INPUT-GD")
ICErhizo$state <- replace(x = ICErhizo$state,
                          ICErhizo$state == "NE" & ICErhizo$ES > 0.03 & ICErhizo$RHIZO <= 0.01,
                          "ICE-RHIZO-GD")
ICErhizo$state <- replace(x = ICErhizo$state,
                          ICErhizo$state == "NE" & ICErhizo$GIFU <= 0.001 & ICErhizo$NFR5 <= 0.001,
                          "ICE-RHIZO-LJAP")
ICErhizo$state <- replace(x = ICErhizo$state,
                          ICErhizo$state == "NE" & ICErhizo$ES > 0.03 & ICErhizo$GIFU <= 0.001,
                          "ICE-RHIZO-GIFU")
ICErhizo$state <- replace(x = ICErhizo$state,
                          ICErhizo$state == "NE" & ICErhizo$ES > 0.03 & ICErhizo$NFR5 <= 0.001,
                          "ICE-RHIZO-NFR5")
ICErhizo$state <- replace(x = ICErhizo$state,
                          ICErhizo$state == "NE" & ICErhizo$ES > 0.03 & ICErhizo$COL0 <= 0.001,
                          "ICE-RHIZO-COL0")

ICErhizo <- ICErhizo %>% mutate(state = as.factor(state))
ICErhizo$state <- factor(ICErhizo$state, levels = c("ICE-INPUT-ES",
                                                    "ICE-INPUT-GD",
                                                    "ICE-RHIZO-GD",
                                                    "NE",
                                                    "ICE-AGAR-ES",
                                                    "ICE-GR-ES",
                                                    "ICE-RHIZO-LJAP",
                                                    "ICE-RHIZO-COL0",
                                                    "ICE-RHIZO-GIFU",
                                                    "ICE-RHIZO-NFR5"))

ICErhizo %>% ggplot(aes(x = ES, y = INPUT, alpha = 0.05, color=state)) +
        geom_point() +
        theme_minimal() + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_vline(xintercept = 0.03) +
        geom_hline(yintercept = 0.008) +
        ggsave("../../Chapter 4 Rhizosphere/ICESym.INPUT.ins_index.png", width = 12, height = 8) +
        ggsave("../../Chapter 4 Rhizosphere/ICESym.INPUT.ins_index.svg", width = 12, height = 8)


ICErhizo %>% ggplot(aes(x = ES, y = RHIZO, alpha = 0.05, color=state)) +
        geom_point() +
        theme_minimal() + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_vline(xintercept = 0.03) +
        geom_hline(yintercept = 0.01) +
        ggsave("../../Chapter 4 Rhizosphere/ICESym.RHIZO.ins_index.png", width = 12, height = 8) +
        ggsave("../../Chapter 4 Rhizosphere/ICESym.RHIZO.ins_index.svg", width = 12, height = 8)

ICErhizo %>% ggplot(aes(x = ES, y = GIFU, alpha = 0.05, color=state)) +
        geom_jitter() +
        theme_minimal() + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_vline(xintercept = 0.03) +
        geom_hline(yintercept = 0.001) +
        ggsave("../../Chapter 4 Rhizosphere/ICESym.RHIZO.GIFU.ins_index.png", width = 12, height = 8) +
        ggsave("../../Chapter 4 Rhizosphere/ICESym.RHIZO.GIFU.ins_index.svg", width = 12, height = 8)

ICErhizo %>% ggplot(aes(x = ES, y = NFR5, alpha = 0.05, color=state)) +
        geom_jitter() +
        theme_minimal() + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_vline(xintercept = 0.03) +
        geom_hline(yintercept = 0.001) +
        ggsave("../../Chapter 4 Rhizosphere/ICESym.RHIZO.NFR5.ins_index.png", width = 12, height = 8) +
        ggsave("../../Chapter 4 Rhizosphere/ICESym.RHIZO.NFR5.ins_index.svg", width = 12, height = 8)

ICErhizo %>% ggplot(aes(x = ES, y = COL0, alpha = 0.05, color=state)) +
        geom_jitter() +
        theme_minimal() + 
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_vline(xintercept = 0.03) +
        geom_hline(yintercept = 0.001) +
        ggsave("../../Chapter 4 Rhizosphere/ICESym.RHIZO.COL0.ins_index.png", width = 12, height = 8) +
        ggsave("../../Chapter 4 Rhizosphere/ICESym.RHIZO.COL0.ins_index.svg", width = 12, height = 8)

write_tsv(ICErhizo, "../../Chapter 4 Rhizosphere/Rhizosphere ICESym updated states.txt")

### Read in tables and parse our data for venn diagrams and functional summaries
#invitro <- invitro[!duplicated(invitro$locus_tag), ]
#write_csv(invitro, file = "Desktop/Chapter 3 In Vitro/Table 3.S1 in vitro TraDIS essentiality results.csv")
library(tidyverse)

invitro <- read_csv("Desktop/Chapter 3 In Vitro/Table 3.S1 in vitro TraDIS essentiality results.csv")
rhizo <- read.csv("~/Desktop/Chapter 4 Rhizosphere/Table 4.S1 Summary of Rhizosphere Tnseq Analysis.csv")

COGTable <- invitro %>% mutate(FUN = as.factor(FUN))
COGTable <- COGTable %>% group_by(state) %>% select(FUN) %>% table() %>% as.data.frame.array()
COGTable$state <- rownames(COGTable)
write_csv(COGTable, "Desktop/Chapter 3 In Vitro/Figure 3 In Vitro Conditionally Essential Genes/Conditional Updated COG Table.csv")

COGTable <- rhizo %>% mutate(FUN = as.factor(FUN))
COGTable <- COGTable %>% group_by(state) %>% select(FUN) %>% table() %>% as.data.frame.array()
COGTable$state <- rownames(COGTable)
write_csv(COGTable, "Desktop/Chapter 4 Rhizosphere/Figure 2 Rhizospheres Conditionally Essential Genes/Rhizosphere Updated COG Table.csv")

TBES <- invitro %>% select(locus_tag, state, TB.state) %>%
        filter(!state %in% c("CORE-ES", "NS-ES", "WT-ES")) %>%
        filter(TB.state %in% c("TB-ES", "WT-TB-ES", "NS-TB-ES")) %>%
        select(locus_tag)
write_tsv(TBES, "Desktop/TBES.tsv", col_names = F)
TYES <- invitro %>% select(locus_tag, state, TY.state) %>%
        filter(!state %in% c("CORE-ES", "NS-ES", "WT-ES")) %>%
        filter(TY.state %in% c("TY-ES", "WT-TY-ES", "NS-TY-ES")) %>% 
        select(locus_tag)
write_tsv(TYES, "Desktop/TYES.tsv", col_names = F)
GRES <- invitro %>% select(locus_tag, state, GR.state) %>%
        filter(!state %in% c("CORE-ES", "NS-ES", "WT-ES")) %>%
        filter(GR.state %in% c("GR-ES", "WT-GR-ES", "NS-GR-ES")) %>% 
        select(locus_tag)
write_tsv(GRES, "Desktop/GRES.tsv", col_names = F)

TBGD <- invitro %>% select(locus_tag, state, TB.state) %>%
        filter(!state %in% c("CORE-ES", "NS-ES", "WT-ES")) %>%
        filter(TB.state %in% c("TB-GD")) %>%
        select(locus_tag)
write_tsv(TBGD, "Desktop/TBGD.tsv", col_names = F)
TYGD <- invitro %>% select(locus_tag, state, TY.state) %>%
        filter(!state %in% c("CORE-ES", "NS-ES", "WT-ES")) %>%
        filter(TY.state %in% c("TY-GD")) %>% 
        select(locus_tag)
write_tsv(TYGD, "Desktop/TYGD.tsv", col_names = F)
GRGD <- invitro %>% select(locus_tag, state, GR.state) %>%
        filter(!state %in% c("CORE-ES", "NS-ES", "WT-ES")) %>%
        filter(GR.state %in% c("GR-GD")) %>% 
        select(locus_tag)
write_tsv(GRGD, "Desktop/GRGD.tsv", col_names = F)

COGTable <- invitro %>% mutate(FUN = as.factor(FUN))
TBCOGTable <- COGTable %>% group_by(TB.state) %>% select(FUN) %>% table() %>% as.data.frame.array()
TBCOGTable$state <- rownames(TBCOGTable)
write_csv(TBCOGTable, "Desktop/Chapter 3 In Vitro/Figure 3 In Vitro Conditionally Essential Genes/TB Updated COG Table.csv")
TYCOGTable <- COGTable %>% group_by(TY.state) %>% select(FUN) %>% table() %>% as.data.frame.array()
TYCOGTable$state <- rownames(TYCOGTable)
write_csv(TYCOGTable, "Desktop/Chapter 3 In Vitro/Figure 3 In Vitro Conditionally Essential Genes/TY Updated COG Table.csv")
GRCOGTable <- COGTable %>% group_by(GR.state) %>% select(FUN) %>% table() %>% as.data.frame.array()
GRCOGTable$state <- rownames(GRCOGTable)
write_csv(GRCOGTable, "Desktop/Chapter 3 In Vitro/Figure 3 In Vitro Conditionally Essential Genes/GR Updated COG Table.csv")


COGTable <- rhizo %>% mutate(FUN = as.factor(FUN))

RHIZOCOGTable <- COGTable %>% group_by(rhizo.state) %>% select(FUN) %>% table() %>% as.data.frame.array()
RHIZOCOGTable$state <- rownames(RHIZOCOGTable)
write_csv(RHIZOCOGTable, "~/Desktop/RHIZO COG Table.csv")

GIFUCOGTable <- COGTable %>% group_by(GIFU.state) %>% select(FUN) %>% table() %>% as.data.frame.array()
GIFUCOGTable$state <- rownames(GIFUCOGTable)
write_csv(GIFUCOGTable, "~/Desktop/GIFU COG Table.csv")

NFR5COGTable <- COGTable %>% group_by(NFR5.state) %>% select(FUN) %>% table() %>% as.data.frame.array()
NFR5COGTable$state <- rownames(NFR5COGTable)
write_csv(NFR5COGTable, "~/Desktop/NFR5 COG Table.csv")

COL0COGTable <- COGTable %>% group_by(COL0.state) %>% select(FUN) %>% table() %>% as.data.frame.array()
COL0COGTable$state <- rownames(COL0COGTable)
write_csv(COL0COGTable, "~/Desktop/COL0 COG Table.csv")

GIFU_ES <- rhizo %>% filter(GIFU.state == "GIFU-ES") %>% select(locus_tag)
GIFU_GD <- rhizo %>% filter(GIFU.state == "GIFU-GD") %>% select(locus_tag)

NFR5_ES <- rhizo %>% filter(NFR5.state == "NFR5-ES") %>% select(locus_tag)
NFR5_GD <- rhizo %>% filter(NFR5.state == "NFR5-GD") %>% select(locus_tag)

COL0_ES <- rhizo %>% filter(COL0.state == "COL0-ES") %>% select(locus_tag)
COL0_GD <- rhizo %>% filter(COL0.state == "COL0-GD") %>% select(locus_tag)

write_csv(GIFU_ES, "~/Desktop/GIFU ES.txt", col_names = F)
write_csv(GIFU_GD, "~/Desktop/GIFU GD.txt", col_names = F)
write_csv(NFR5_ES, "~/Desktop/NFR5 ES.txt", col_names = F)
write_csv(NFR5_GD, "~/Desktop/NFR5 GD.txt", col_names = F)
write_csv(COL0_ES, "~/Desktop/COL0 ES.txt", col_names = F)
write_csv(COL0_GD, "~/Desktop/COL0 GD.txt", col_names = F)
