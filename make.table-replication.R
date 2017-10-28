## show repliated probes
source('util.R')
library(ggplot2)
library(pander)
library(xlsx)

load('replication.summary.rdata') ## full replication analysis
mwas.tot <- read_tsv('table_mwas.txt.gz')
mwas.sig <- read_tsv('table_mwas_significant.txt.gz')
dejager <- read_tsv('dejager.cpgs.txt', col_names = 'cg')

discovery.mwas <- melt.mwas.tab(gauss.rep.data %>% filter(type == 'large'))
replication.mwas <- melt.mwas.tab(gauss.rep.data %>% filter(type == 'small'))

gene.cols <- c('cg', 'chr', 'cg.loc', 'ensg', 'hgnc', 'tss', 'tes', 'strand', 'dist')
genes.tab <- read_tsv('nearest.genes.gz', col_names = gene.cols)

cg.key <- c('cg', 'chr', 'cg.loc', 'pheno')

matched <- discovery.mwas %>% select(cg, chr, cg.loc, pheno, theta) %>%
    left_join(y = replication.mwas %>% select(cg, chr, cg.loc, pheno, theta),
              by = cg.key,
              suffix = c('.discovery', '.replication'))

concordant <- matched %>%
    filter(sign(theta.discovery) == sign(theta.replication)) %>%
    select(cg, chr, cg.loc, pheno)

out <- concordant %>% left_join(mwas.sig, by = cg.key) %>%
    na.omit() %>% filter(pheno %in% c('NP', 'NFT', 'Cog')) %>%
    left_join(discovery.mwas, by = cg.key, suffix = c('', '.disc')) %>%
    left_join(replication.mwas, by = cg.key, suffix = c('', '.repl'))

format.tab <- function(tab) {
    ret <- tab %>%
        mutate(location = format(cg.loc, big.mark = ','), p.val = signif(p.val, 2)) %>%
            mutate(joint = signif(theta, 2) %&&% ' (' %&&% signif(theta.se,1) %&&% ')') %>%
                mutate(discovery = signif(theta.disc, 2) %&&% ' (' %&&% signif(theta.se.disc,1) %&&% ')') %>%
                    mutate(replication = signif(theta.repl, 2) %&&% ' (' %&&% signif(theta.se.repl,1) %&&% ')')

    ret <- ret %>% select(cg, chr, location, pheno, p.val,
                          joint, discovery, replication, hgnc)
    
    return(ret)
}

take.repl.qc <- function(tab) {

    lodds.cutoff <- log(.95) - log(.05)
    dist.cutoff <- 1

    ret <- tab %>%
        filter(lodds > lodds.cutoff, lodds.disc > 0, lodds.repl > 0) %>%
        arrange(desc(lodds), desc(lodds.disc), desc(lodds.repl))

    ret <- ret %>%
        left_join(genes.tab, by = c('cg', 'chr', 'cg.loc')) %>%
            filter(dist < dist.cutoff) %>%
                select(-dist, -tss, -tes, -strand) %>%
                    as.data.frame()

    return(ret)
}

np.out <- take.repl.qc(out %>% filter(pheno == 'NP'))
nft.out <- take.repl.qc(out %>% filter(pheno == 'NFT'))
cog.out <- take.repl.qc(out %>% filter(pheno == 'Cog'))

np.out.10 <- np.out %>% head(10) %>% arrange(chr, cg.loc) %>% format.tab()
nft.out.10 <- nft.out %>% head(10) %>% arrange(chr, cg.loc) %>% format.tab()
cog.out.10 <- cog.out %>% head(10) %>% arrange(chr, cg.loc) %>% format.tab()

small.out.tab <- rbind(np.out.10, nft.out.10, cog.out.10)

tab.ret <- pandoc.table.return(small.out.tab, caption = 'Top CpGs',
                               style = 'simple', digits = 2,
                               split.tables = 200)

cat(tab.ret, file = 'table_top_cpgs.md')

write.xlsx(np.out, file = 'table_replication.xlsx', sheetName = 'NP', append = FALSE, row.names = FALSE)
write.xlsx(nft.out, file = 'table_replication.xlsx', sheetName = 'NFT', append = TRUE, row.names = FALSE)
write.xlsx(cog.out, file = 'table_replication.xlsx', sheetName = 'Cog', append = TRUE, row.names = FALSE)

write.xlsx(np.out %>% format.tab(), file = 'table_replication_formatted.xlsx', sheetName = 'NP', append = FALSE, row.names = FALSE)
write.xlsx(nft.out %>% format.tab(), file = 'table_replication_formatted.xlsx', sheetName = 'NFT', append = TRUE, row.names = FALSE)
write.xlsx(cog.out %>% format.tab(), file = 'table_replication_formatted.xlsx', sheetName = 'Cog', append = TRUE, row.names = FALSE)



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
    scale_size_continuous('PIP', range = c(0, 2), breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.9)) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 6),
          axis.title = element_blank(),
          strip.text = element_text(size = 8),
          strip.background = element_blank())

pdf(file = 'figure-replication/fig_dejager.pdf', width = 8, height = 4)
print(plt)
dev.off()
