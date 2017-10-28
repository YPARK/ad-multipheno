source('util.R')
library(ggplot2)
library(pander)
library(xlsx)
library(dplyr)
library(tidyr)

repl.tab <- read_tsv('table_replication.txt.gz')
load('summary-G.rdata')

var.names <- read_tsv('phenotype.names.txt', col_names = 'var')
var.names <- var.names$var

colnames(gauss.var.decomp) <- c('cg', 'chr', 'cg.loc',
                                'var.obs', 'total',
                                var.names,
                                'genetic', 'MF')

pve.stat.tot <- gauss.var.decomp %>%
    gather(variable, pve, c(var.names, 'genetic')) %>%
    group_by(variable) %>%
    summarize(pve.mean = mean(pve),
              pve.sd = sd(pve),
              pve.se = sd(pve) / sqrt(n()))

cg.pheno <- repl.tab %>%
    select(cg, pheno) %>%
    unique() %>%
    group_by(cg) %>%
    summarize(pheno = paste(pheno, collapse = '|'))

var.tab <- cg.pheno %>%
    left_join(gauss.var.decomp, by = 'cg') %>%
    arrange(chr, cg.loc, pheno)

pve.tab <- var.tab %>%
    select(-chr, -cg.loc, -var.obs) %>%
    gather(variable, pve, c('genetic', var.names))

pve.stat <- pve.tab %>% group_by(variable, pheno) %>%
    summarize(pve.mean = mean(pve),
              pve.sd = sd(pve),
              pve.se = sd(pve) / sqrt(n())) %>%
    arrange(pheno)

pve.stat$variable <- factor(pve.stat$variable,
                            c(var.names, 'genetic'))

pve.stat.tot$variable <- factor(pve.stat.tot$variable,
                                c(var.names, 'genetic'))

pve.stat$pheno <- factor(pve.stat$pheno, c('NP', 'NFT', 'Cog'))

gg.plot <- function(...) {
    ggplot(...) + theme_bw() + theme(panel.background = element_blank())
}

plt.pve.stat <- gg.plot() +
    geom_bar(data = pve.stat.tot, aes(x = variable, y = pve.mean),
             fill = '#009900', stat = 'identity')

plt.pve.stat <- plt.pve.stat +
    geom_errorbar(data = pve.stat.tot, aes(x = variable,
                      ymin = pmax(pve.mean - 2 * pve.se, 0),
                      ymax = pmin(pve.mean + 2 * pve.se, 1)),
                  color = '#009900', width = .1)

plt.pve.stat <-
    plt.pve.stat +
    geom_linerange(data = pve.stat, aes(x = variable,
                      ymin = pve.mean - 2 * pve.se,
                      ymax = pve.mean + 2 * pve.se)) +
    geom_point(data = pve.stat, aes(x = variable, y = pve.mean)) +
    facet_grid(pheno ~ .) +
    ylab('proportion of variance explained (DNAme)') +
    theme(axis.text.x = element_text(size=8, angle = 45, hjust = 1, vjust = 1))

pdf('figure-replication/fig_pve_stat_tot.pdf', width = 3, height = 7)
print(plt.pve.stat)
dev.off()

################################################################
pve.tab$variable <- factor(pve.tab$variable, c(var.names, 'genetic'))
pve.tab$pheno <- factor(pve.tab$pheno, c('NP', 'NFT', 'Cog'))

temp <- pve.tab %>% spread(key = variable, value = pve) %>%
    select(-pheno, -total, -MF)

cg.list <- temp$cg
cg.order <- row.order(temp %>% select(-cg) %>% as.matrix())
cg.ordered <- cg.list[cg.order]

pve.tab$cg <- factor(pve.tab$cg, cg.ordered)

plt.pve.heatmap <- gg.plot(pve.tab, aes(y = cg, x = variable, fill = pve)) +
    geom_tile() +
    facet_grid(pheno ~ ., scales = 'free', space = 'free') +
    scale_fill_gradientn('PVE',
                         colors = c('white', 'gray40', 'gray20', 'black'),
                         breaks = c(0, 0.1, 0.25, 0.5), trans = 'sqrt') +
    theme(legend.position = 'bottom',
          legend.text = element_text(size=6),
          axis.text.y = element_text(size=4), axis.title = element_blank(),
          axis.text.x = element_text(size=8, angle = 45, hjust = 1, vjust = 1))

pdf('figure-replication/fig_pve_stat_heatmap.pdf', width = 3, height = 7)
print(plt.pve.heatmap)
dev.off()

################################################################
write.xlsx(pve.tab %>% spread(key = variable, value = pve) %>% as.data.frame(),
           file = 'table_replication_pve.xlsx', row.names = FALSE)
write_tsv(pve.tab, path = gzfile('table_replication_pve.txt.gz'))
