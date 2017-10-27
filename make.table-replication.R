## show repliated probes
source('util.R')
library(ggplot2)
library(pander)

load('replication.summary.rdata') ## full replication analysis
mwas.tot <- read_tsv('table_mwas.txt.gz')
mwas.sig <- read_tsv('table_mwas_significant.txt.gz')
dejager <- read_tsv('dejager.cpgs.txt', col_names = 'cg')

discovery.mwas <- melt.mwas.tab(gauss.rep.data %>% filter(type == 'large'))
replication.mwas <- melt.mwas.tab(gauss.rep.data %>% filter(type == 'small'))

gene.cols <- c('cg', 'chr', 'cg.loc', 'ensg', 'tss', 'tes', 'hgnc', 'strand', 'dist')
genes.tab <- read_tsv('nearest.genes.gz', col_names = gene.cols)

cg.key <- c('cg', 'chr', 'cg.loc', 'pheno')

matched <- discovery.mwas %>% select(cg, chr, cg.loc, pheno, theta) %>%
    left_join(y = replication.mwas %>% select(cg, chr, cg.loc, pheno, theta),
              by = cg.key,
              suffix = c('.discovery', '.replication'))

concordant <- matched %>%
    filter(sign(theta.discovery) == sign(theta.replication)) %>%
    select(cg, chr, cg.loc, pheno)

lodds.cutoff <- log(.95) - log(.05)

out <- concordant %>% left_join(mwas.sig, by = cg.key) %>%
    na.omit() %>% filter(pheno %in% c('NP', 'NFT', 'Cog')) %>%
    left_join(discovery.mwas, by = cg.key, suffix = c('', '.discovery')) %>%
    left_join(replication.mwas, by = cg.key, suffix = c('', '.replication'))

np.out <- out %>%
    filter(pheno == 'NP', lodds > lodds.cutoff, lodds.discovery > 0, lodds.replication > 0) %>%
    arrange(desc(lodds), desc(lodds.discovery), desc(lodds.replication))

nft.out <- out %>%
    filter(pheno == 'NFT', lodds > lodds.cutoff, lodds.discovery > 0, lodds.replication > 0) %>%
    arrange(desc(lodds), desc(lodds.discovery), desc(lodds.replication))

cog.out <- out %>%
    filter(pheno == 'Cog', lodds > lodds.cutoff, lodds.discovery > 0, lodds.replication > 0) %>%
    arrange(desc(lodds), desc(lodds.discovery), desc(lodds.replication))



np.out.10 <- np.out %>% head(10) %>% arrange(chr, cg.loc)
nft.out.10 <- out %>% head(10) %>% arrange(chr, cg.loc)
cog.out.10 <- out %>% head(10) %>% arrange(chr, cg.loc)






################################################################
## just show replication of De Jager CpGs

dejager.tab <- rbind(dejager %>% left_join(discovery.mwas, by = 'cg') %>%
                     filter(pheno %in% c('NP', 'NFT', 'Cog')) %>%
                     mutate(batch = 'discovery'),
                     dejager %>% left_join(replication.mwas, by = 'cg') %>%
                     filter(pheno %in% c('NP', 'NFT', 'Cog')) %>%
                     mutate(batch = 'replication'),
                     dejager %>% left_join(mwas.tot, by = 'cg') %>%
                     select(cg, chr, cg.loc, pheno, theta, theta.se, lodds) %>%
                     mutate(batch = 'joint') %>% na.omit()) %>%
    mutate(trunc.z = pmin(pmax(theta/theta.se, -3), 3))

dejager.tab$batch <- factor(dejager.tab$batch, c('replication', 'discovery', 'joint'))
dejager.tab$pheno <- factor(dejager.tab$pheno, c('NP', 'NFT', 'Cog'))

gg.plot <- function(...) {
    ggplot(...) + theme_bw() + theme(panel.background = element_blank())
}

.aes <- aes(y = batch, x = cg, fill = trunc.z, size = 1/(1+exp(-lodds)))

plt <- gg.plot() +
    geom_point(data = dejager.tab, .aes, pch = 22) +
    facet_grid(pheno ~ chr, scales = 'free', space = 'free') +
    scale_fill_gradientn('effect z-score', colors = c('blue', 'gray', 'yellow')) +
    scale_size_continuous('PIP', range = c(0, 4), breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.9)) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
          axis.title = element_blank())

