### Code for chapter 3. In Vitro Tn-seq

setwd("~/Desktop/PhD/")
set.seed(seed = 1953)

library(tidyverse)



readTradis <- function(wrkingdir = "") {
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
  return(tradisTable)
}

# 1. Compile insertion index summary tables for analysis.

ES <- readTradis(wrkingdir = "~/Desktop/PhD/data/tradis/essential/")
INVITRO <- readTradis(wrkingdir = "~/Desktop/PhD/data/tradis/in vitro/")



# 2. 



