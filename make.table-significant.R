library(feather)
library(dplyr)
library(readr)
library(xlsx)
library(pander)

options(stringsAsFactors = FALSE)

mwas.tab <- read_feather('table_mwas.ft')

mwas.boot.tab <- mwas.tab %>% na.omit() %>%
    mutate(p.val = pnorm((lodds -boot.lodds)/boot.lodds.se, lower.tail = FALSE))

n.tot <- mwas.tab %>% filter(pheno == 'NP') %>% nrow()

cutoff <- 0.05 / n.tot / 3
significant <- mwas.boot.tab %>% filter(p.val < cutoff, lodds > 0) %>% select(cg) %>% unique()
mwas.sig <- mwas.boot.tab %>% filter(cg %in% significant$cg)

## output to feather for visualization
.write.tab <- function(...) {
    ## write.table(..., row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
    write_tsv(..., col_names = TRUE)
}

pheno.names <- c('NP', 'NFT', 'Cog')
.write.tab(mwas.tab %>% filter(pheno %in% pheno.names), path = gzfile('table_mwas.txt.gz'))
.write.tab(mwas.sig, path = gzfile('table_mwas_significant.txt.gz'))
