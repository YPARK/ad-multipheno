source('util.R')
library(ggplot2)
library(ggrepel)
library(pander)
library(xlsx)
library(dplyr)
library(tidyr)

################################################################
## read all data
mwas.tab <- read_tsv('table_mwas.txt.gz') %>% rename(variable = pheno)
pve.tab <- read_tsv('table_pve_tot.txt.gz')

var.names <- read_tsv('phenotype.names.txt', col_names = 'var')
var.names <- var.names$var
pve.melt <- pve.tab %>% gather(key = 'variable', value = 'pve', c(var.names, 'genetic'))

full.tab <- mwas.tab %>% left_join(pve.melt)

n.tot <- 413223
cutoff <- 0.05 / n.tot / 3
log.cutoff <- -log10(cutoff)

gene.cols <- c('cg', 'chr', 'cg.loc', 'ensg', 'hgnc', 'tss', 'tes', 'strand', 'dist')
genes.tab <- read_tsv('nearest.genes.gz', col_names = gene.cols)

################################################################
## discrepancy between log-odds and -log10 P

plt <- gg.plot(full.tab)

p1 <- plt +
    geom_hex(aes(x = 1/(1+exp(-lodds)), y = -log10(p.val)), bins = 50) +
    scale_fill_gradientn(colors = c('blue', 'gray', 'yellow'), trans = 'log10') +
    geom_hline(yintercept = log.cutoff, lty = 2, color = 'red') +
    xlim(c(0, 1)) +
    xlab('Posterior inclusion probability')+
    ylab('-log10 p-value of parametric bootstrap')

p2 <- plt +
    geom_hex(aes(x = 1/(1+exp(-lodds)), y = pve), bins = 50) +
    scale_fill_gradientn(colors = c('blue', 'gray', 'yellow'), trans = 'log10') +
    xlim(c(0, 1)) +
    xlab('Posterior inclusion probability') +
    ylab('Proportion of DNAme variance explained')

p3 <- plt +
    geom_hex(aes(-log10(p.val), y = pve), bins = 50) +
    scale_fill_gradientn(colors = c('blue', 'gray', 'yellow'), trans = 'log10') +
    scale_x_sqrt() +
    scale_y_sqrt(breaks=c(1e-3, 1e-2, 0.1, 0.2, 0.3, 0.4, 0.5)) +
    geom_vline(xintercept = log.cutoff, lty = 2, color = 'red') +
    xlab('-log10 p-value of parametric bootstrap') +
    ylab('Proportion of DNAme variance explained')

ggsave(filename = 'fig_pip_pvalue.pdf', plot=p1, width = 4, height = 4)
ggsave(filename = 'fig_pip_pve.pdf', plot=p2, width = 4, height = 4)
ggsave(filename = 'fig_pve_pvalue.pdf', plot=p3, width = 4, height = 4)

################################################################
## show local views on significant CpGs
pheno.names <- c('NP', 'NFT', 'Cog')
pheno.colors <- c('#FFAA00','#AA55FF','#00DD99')
pheno.shapes <- c(21,22,25)

tab.to.show <- read_tsv('table_pve_mult_top.txt.gz')
cg.to.show <- tab.to.show$cg

.cg <- cg.to.show[2]
.loc <- mwas.tab %>% filter(cg == .cg) %>% select(chr, cg.loc) %>% unique()
.chr <- .loc$chr
lb <- .loc$cg.loc - 5e5
ub <- .loc$cg.loc + 5e5

.neigh.tab <- full.tab %>%
    filter(chr == .chr, cg.loc > lb, cg.loc < ub) %>%
    filter(variable %in% pheno.names) %>%
    mutate(p.val.trunc = pmin(pmax(-log10(p.val), 0), 10)) %>%
    mutate(z.trunc = pmin(pmax(theta/theta.se, -3), 3)) %>%
    mutate(pip = 1/(1+exp(-lodds)))

.neigh.tab$variable <- factor(.neigh.tab$variable, pheno.names)
.neigh.tab$p.val.trunc[is.na(.neigh.tab$p.val.trunc)] <- 1

.neigh.genes <- genes.tab %>% filter(cg %in% unique(.neigh.tab$cg)) %>%
    select(hgnc, tss, tes, strand) %>% unique() %>%
    arrange(tss)

#y=(theta - boot.theta)/boot.theta.se,

.aes <- aes(x = cg.loc, y = pve, shape = variable, fill = z.trunc, size = p.val.trunc)
.aes.lab <- aes(x = cg.loc, y = pve, label = cg)

##            size = 1/(1+exp(-(lodds - boot.lodds))))

## .neig %>% filter(p.val.trunc > log.cutoff)

.neigh.tab.lab <- .neigh.tab %>% na.omit() %>%
    filter(p.val.trunc > log.cutoff, pve > 1e-2)

gg.plot() +
    geom_point(data = .neigh.tab, .aes) +
    geom_text_repel(data = .neigh.tab.lab, .aes.lab) +
    scale_color_manual(values = pheno.colors) +
    scale_shape_manual(values = pheno.shapes) +
    scale_fill_gradientn(colors = c('blue', 'gray', 'yellow')) +
    scale_size_continuous(breaks=c(0, 1, 3, signif(log.cutoff, 2), 10, 20), range = c(0.5, 2)) +
    facet_grid(variable ~ ., scales='free') +
    xlim(lb, ub)


