library(dplyr)
library(readr)
library(feather)

mwas.tab <- read_feather('../table_mwas.ft')
mwas.boot.tab <- mwas.tab %>% na.omit() %>%
    mutate(p.val = pnorm((lodds -boot.lodds)/boot.lodds.se, lower.tail = FALSE))

n.tot <- mwas.tab %>% filter(pheno == 'NP') %>% nrow()
cutoff <- 0.05 / n.tot / 3

mwas.sig <- mwas.boot.tab %>% filter(p.val < cutoff, lodds > 0) %>%
    dplyr::mutate(lb = cg.loc - 1, ub = cg.loc) %>%
    dplyr::select(chr, lb, ub, cg, pheno, theta, theta.se, lodds, p.val)

write.table(mwas.sig, file = gzfile('mwas.bed.gz'),
            row.names = FALSE, col.names = FALSE,
            sep = '\t', quote = FALSE)
