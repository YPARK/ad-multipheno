source('util.R')
library(ggplot2)
library(pander)
library(xlsx)
library(dplyr)
library(tidyr)

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

mwas.sig <- read_tsv('table_mwas_significant.txt.gz')

cg.pheno <- mwas.sig %>%
    select(cg, pheno) %>%
    unique() %>%
    group_by(cg) %>%
    summarize(pheno = paste(sort(pheno, decreasing=TRUE), collapse = '+'))

var.tab <- cg.pheno %>%
    left_join(gauss.var.decomp, by = 'cg')

pve.tab <- var.tab %>%
    select(-chr, -cg.loc, -var.obs) %>%
    gather(variable, pve, c('genetic', var.names)) %>%
    na.omit()

pve.stat <- pve.tab %>% group_by(variable, pheno) %>%
    summarize(pve.mean = mean(pve),
              pve.sd = sd(pve),
              pve.se = sd(pve) / sqrt(n()),
              n = n()) %>%
    arrange(pheno) %>% filter(n > 1)

pve.stat$pheno <- factor(pve.stat$pheno, c('NP', 'NFT', 'Cog', 'NP+NFT', 'NP+Cog', 'NFT+Cog'))

pve.stat$variable <- factor(pve.stat$variable, c(var.names, 'genetic'))

pve.stat.tot$variable <- factor(pve.stat.tot$variable, c(var.names, 'genetic'))

gg.plot <- function(...) {
    ggplot(...) + theme_bw() + theme(panel.background = element_blank())
}

plt.pve.stat <- gg.plot() +
    geom_bar(data = pve.stat.tot, aes(x = variable, y = pve.mean),
             fill = '#009900', stat = 'identity')

plt.pve.stat <- plt.pve.stat +
    geom_errorbar(data = pve.stat.tot, aes(x = variable,
                      ymin = pmax(pve.mean - 2 * pve.se, 1e-4),
                      ymax = pmin(pve.mean + 2 * pve.se, 1)),
                  color = '#009900', width = .1)

plt.pve.stat <-
    plt.pve.stat +
    geom_linerange(data = pve.stat, aes(x = variable,
                      ymin = pmax(pve.mean - 2 * pve.se, 1e-4),
                      ymax = pmin(pve.mean + 2 * pve.se, 1))) +
    geom_point(data = pve.stat, aes(x = variable, y = pve.mean)) +
    facet_wrap(~pheno, scales = 'free', nrow = 3, dir = 'v') +
    ylab('proportion of variance explained (DNAme)') +
    theme(axis.text.x = element_text(size=8, angle = 45, hjust = 1, vjust = 1))

plt.pve.stat <- 
    plt.pve.stat +
    scale_y_continuous(trans = 'sqrt', breaks = c(1e-3, 1e-2, 0.025, 0.05, 0.1, 0.2))

pdf('fig_pve_stat_tot.pdf', width = 5, height = 7)
print(plt.pve.stat)
dev.off()

################################################################
## show double phenotype stuff
pve.tab$variable <- factor(pve.tab$variable, c(var.names, 'genetic'))
pve.tab$pheno <- factor(pve.tab$pheno, c('NP', 'NFT', 'Cog', 'NP+NFT', 'NP+Cog', 'NFT+Cog'))

pve.tab.mult <- pve.tab %>% filter(pheno %in% c('NP+NFT', 'NP+Cog', 'NFT+Cog'))

temp <- pve.tab.mult %>% spread(key = variable, value = pve) %>%
    select(-pheno, -total, -MF)

cg.list <- temp$cg
cg.order <- row.order(temp %>% select(-cg) %>% as.matrix())
cg.ordered <- cg.list[cg.order]

pve.tab.mult$cg <- factor(pve.tab.mult$cg, cg.ordered)

plt.pve.heatmap <-
    gg.plot(pve.tab.mult, aes(y = cg, x = variable, fill = pve)) +
    geom_tile() +
    facet_grid(pheno ~ ., scales = 'free', space = 'free') +
    scale_fill_gradientn('PVE',
                         colors = c('white', 'gray40', 'gray20', 'black'),
                         breaks = c(0, 0.1, 0.25, 0.5), trans = 'sqrt') +
    theme(legend.position = 'bottom',
          legend.text = element_text(size=6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8, angle = 45, hjust = 1, vjust = 1))

pdf('fig_pve_stat_heatmap.pdf', width = 2.5, height = 7)
print(plt.pve.heatmap)
dev.off()

################################################################
write_tsv(pve.tab, path = gzfile('table_pve.txt.gz'))
