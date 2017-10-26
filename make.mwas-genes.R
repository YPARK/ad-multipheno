library(feather)
library(dplyr)
library(readr)
library(xlsx)
library(pander)

options(stringsAsFactors = FALSE)

mwas.tab <- read_feather('table_mwas.ft')
m2t.tab <- read_tsv('m2t.txt.gz', col_names = TRUE)

mwas.boot.tab <- mwas.tab %>% na.omit() %>%
    mutate(p.val = pnorm((lodds -boot.lodds)/boot.lodds.se, lower.tail = FALSE))

n.tot <- mwas.tab %>% filter(pheno == 'NP') %>% nrow()

cutoff <- 0.05 / n.tot / 3

mwas.gene.sig <- mwas.boot.tab %>% filter(p.val < cutoff, lodds > 0) %>%
    left_join(m2t.tab, by = c('chr', 'cg', 'cg.loc')) %>% na.omit() %>%
    dplyr::select(chr, ld.lb, ld.ub, cg, cg.loc, pheno, theta, theta.se, lodds,
                  boot.lodds, boot.lodds.se, p.val,
                  ensg, hgnc, tss, tes, m2t.theta, m2t.theta.var) %>%
    dplyr::mutate(m2t.theta.var = sqrt(m2t.theta.var)) %>%
    dplyr::rename(m2t.theta.se = m2t.theta.var)

mwas.sig <- mwas.boot.tab %>% filter(p.val < cutoff, lodds > 0)

## output to feather for visualization
.write.tab <- function(...) {
    ## write.table(..., row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
    write_tsv(..., col_names = TRUE)
}

pheno.names <- c('NP', 'NFT', 'Cog')
.write.tab(mwas.tab %>% filter(pheno %in% pheno.names), path = gzfile('table_mwas.txt.gz'))
.write.tab(mwas.gene.sig, path = gzfile('table_mwas_significant_gene.txt.gz'))
.write.tab(mwas.sig, path = gzfile('table_mwas_significant.txt.gz'))

## output to Excel
mwas.sig <- mwas.boot.tab %>% filter(p.val < cutoff)
xlsx.file <- 'table_mwas_significant.xlsx'

write.xlsx(mwas.gene.sig %>% dplyr::select(-ld.lb, -ld.ub),
           file = xlsx.file, sheetName = 'linked.genes', append = FALSE)
