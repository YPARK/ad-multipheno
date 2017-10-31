library(feather)
library(dplyr)
library(readr)
library(xlsx)
library(pander)

options(stringsAsFactors = FALSE)

mwas.tab <- read_feather('table_mwas.ft')

mwas.boot.tab <- mwas.tab %>%
    mutate(p.val = pnorm((lodds -boot.lodds)/boot.lodds.se, lower.tail = FALSE)) %>%
    mutate(p.val.2 = 2*pnorm(abs(theta -boot.theta)/boot.theta.se, lower.tail = FALSE)) %>%
    mutate(p.val = pmax(p.val, p.val.2)) %>%
    select(-p.val.2)

n.tot <- mwas.tab %>% filter(pheno == 'NP') %>% nrow()

cutoff <- 0.05 / n.tot / 3
significant <- mwas.boot.tab %>% filter(p.val < cutoff) %>% select(cg) %>% unique()
mwas.sig <- mwas.boot.tab %>% filter(cg %in% significant$cg)

## output to feather for visualization
.write.tab <- function(...) {
    ## write.table(..., row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
    write_tsv(..., col_names = TRUE)
}

pheno.names <- c('NP', 'NFT', 'Cog')
.write.tab(mwas.boot.tab %>% filter(pheno %in% pheno.names), path = gzfile('table_mwas.txt.gz'))
.write.tab(mwas.sig, path = gzfile('table_mwas_significant.txt.gz'))
