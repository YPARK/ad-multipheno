library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(piano)
library(pander)
library(xlsx)

options(stringsAsFactors = FALSE)
source('util.R')

parse.gsc <- function(gsc.in) {

    rm.header <- function(s) paste(strsplit(s, '[_]')[[1]][-1], collapse = '_')
    nn <- length(gsc.in$gsc)
    gsc.names <- names(gsc.in$gsc)
    ret <- lapply(1:nn, function(g) data.frame(hgnc = gsc.in$gsc[[g]], pathway = gsc.names[g]))
    ret <- do.call(rbind, ret)
    ret <- ret %>% mutate(pathway = sapply(as.character(pathway), rm.header)) %>%
        as.data.frame()

    return(ret)
}

pheno.names <- c('NP', 'NFT', 'Cog')
pheno.colors <- c('#FFAA00','#AA55FF','#00DD99')
pheno.shapes <- c(21,22,25)
pheno.names <- read_tsv('phenotype.names.txt', col_names = 'var')
pheno.names <- pheno.names$var

mwas.tab <- read_tsv('table_mwas.txt.gz')

cpg2gene <- read_tsv('nearest.genes.gz', col_names = FALSE)
colnames(cpg2gene) <- c('cg', 'chr', 'cg.loc', 'ensg', 'hgnc', 'tss', 'tes', 'strand', 'dist')

mwas.gene.tab <- mwas.tab %>%
    left_join(cpg2gene, by = c('cg', 'chr', 'cg.loc')) %>%
    filter(dist < 1) %>%
    select(-dist, -boot.theta, -boot.theta.se, -boot.lodds, -boot.lodds.se) %>%
    na.omit()


## pathway x cpg (gene)
kegg.tab <- loadGSC('genesets/c2.cp.kegg.v6.0.symbols.gmt') %>% parse.gsc() %>% mutate(gs.type = 'KEGG')
reactome.tab <- loadGSC('genesets/c2.cp.reactome.v6.0.symbols.gmt') %>% parse.gsc() %>% mutate(gs.type = 'REACTOME')

fdr.cutoff <- .2

path.sig <- rbind(read_tsv('pathway/full-KEGG_Cog.txt.gz', col_names = TRUE) %>%
    rename(pathway = gs) %>%
    filter(q.val.up < fdr.cutoff | q.val.dn < fdr.cutoff) %>% mutate(gs.type = 'KEGG'),
    read_tsv('pathway/full-KEGG_NP.txt.gz', col_names = TRUE) %>%
    rename(pathway = gs) %>%
    filter(q.val.up < fdr.cutoff | q.val.dn < fdr.cutoff) %>% mutate(gs.type = 'KEGG'),
    read_tsv('pathway/full-KEGG_NFT.txt.gz', col_names = TRUE) %>%
    rename(pathway = gs) %>%
    filter(q.val.up < fdr.cutoff | q.val.dn < fdr.cutoff) %>% mutate(gs.type = 'KEGG'),
    read_tsv('pathway/full-REACTOME_Cog.txt.gz', col_names = TRUE) %>%
    rename(pathway = gs) %>%
    filter(q.val.up < fdr.cutoff | q.val.dn < fdr.cutoff) %>% mutate(gs.type = 'REACTOME'),
    read_tsv('pathway/full-REACTOME_NP.txt.gz', col_names = TRUE) %>%
    rename(pathway = gs) %>%
    filter(q.val.up < fdr.cutoff | q.val.dn < fdr.cutoff) %>% mutate(gs.type = 'REACTOME'),
    read_tsv('pathway/full-REACTOME_NFT.txt.gz', col_names = TRUE) %>%
    rename(pathway = gs) %>%
    filter(q.val.up < fdr.cutoff | q.val.dn < fdr.cutoff) %>% mutate(gs.type = 'REACTOME'))

## pathway x phenotype
tab.ret <- pandoc.table.return(path.sig %>% as.data.frame(),
                               caption = 'Enriched pathways',
                               style = 'simple', digits = 2,
                               split.tables = 200)

cat(tab.ret, file = 'table_pathways.md')

write.xlsx(data.frame(path.sig), file = 'table_pathways.xlsx')


## pathway x gene
gene.summary.tab <- mwas.gene.tab %>% 
    filter(lodds > 0) %>%
    mutate(z = theta/theta.se) %>%
    group_by(hgnc, pheno) %>%
    summarize(z = max(z)) %>%
    mutate(trunc.z = pmin(pmax(z, -3), 3))

full.tab <- path.sig %>% full_join(rbind(kegg.tab, reactome.tab)) %>%
    full_join(gene.summary.tab) %>% na.omit()

temp <- full.tab %>% select(pathway, trunc.z, hgnc) %>% spread(key = pathway, value = trunc.z)
hgnc <- temp$hgnc
temp.mat <- temp %>% select(-hgnc) %>% as.matrix()
temp.mat[!is.finite(temp.mat)] <- 0
hgnc.order <- row.order(temp.mat)

temp <- full.tab %>% select(pathway, trunc.z, hgnc) %>% spread(key = hgnc, value = trunc.z)
path <- temp$pathway
temp.mat <- temp %>% select(-pathway) %>% as.matrix()
temp.mat[!is.finite(temp.mat)] <- 0
path.order <- row.order(temp.mat)

full.tab$pathway <- factor(full.tab$pathway, path[path.order])
full.tab$hgnc <- factor(full.tab$hgnc, hgnc[hgnc.order])

plt <- gg.plot(full.tab, aes(x=hgnc, y=pathway, fill=trunc.z)) +
    geom_tile() +
    facet_grid(gs.type ~ pheno, scales = 'free', space = 'free') +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.background = element_rect(fill = '#DDDDDD'),
          strip.background = element_blank(),
          legend.position = 'bottom') +
    scale_fill_gradient2('gene-level z', low='blue', high='red') +
    xlab('genes bearing significant CpGs')

ggsave(filename = 'fig_pathway_genes.pdf', plot = plt, width = 8, height = 8)

