library(ggplot2)
library(dplyr)
library(readr)
library(feather)

options(stringsAsFactors = FALSE)

pheno.names <- c('NP', 'NFT', 'Cog')
pheno.colors <- c('#FFAA00','#AA55FF','#00DD99')
pheno.shapes <- c(21,22,25)
pheno.names <- read_tsv('phenotype.names.txt', col_names = FALSE)[,1]

mwas.tab <- read_feather('table_mwas.ft')

mwas.boot.tab <- mwas.tab %>% na.omit() %>%
    mutate(p.val = pnorm((lodds -boot.lodds)/boot.lodds.se, lower.tail = FALSE))

n.tot <- mwas.tab %>% filter(pheno == 'NP') %>% nrow()
cutoff <- 0.05 / n.tot / 3

lodds.cutoff <- log(.9) - log(.1)

mwas.sig <- mwas.boot.tab %>%
    filter(p.val < cutoff)

multi.pheno.tab <- mwas.sig %>% group_by(cg) %>%
    summarize(n = length(unique(pheno))) %>%
    filter(n > 1)

multi.mwas.sig <- mwas.tab %>% filter(cg %in% multi.pheno.tab$cg) %>%
    mutate(trunc.z = pmax(pmin(theta/theta.se, 4), -4))

ggplot(multi.mwas.sig, aes(x = pheno, y = cg, fill = trunc.z)) +
    geom_tile() +
    scale_fill_gradient2(low = 'blue', high = 'red')


## look at



