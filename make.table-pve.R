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

n.tot <- 413223
cutoff <- 0.05 / n.tot / 3

cg.pheno <- mwas.sig %>%
    filter(p.val < cutoff) %>%
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
    arrange(pheno) %>% filter(n > 1) %>% select(-n) %>%
    mutate(cpg.data = 'MWAS')

pve.df <- rbind(pve.stat,
                pve.stat %>% select(variable, pheno) %>%
                left_join(pve.stat.tot, by = 'variable') %>%
                mutate(cpg.data = 'Total'))

pve.df$pheno <- factor(pve.df$pheno,
                       c('NP', 'NFT', 'Cog', 'NP+NFT', 'NP+Cog', 'NFT+Cog',
                         'NP+NFT+Cog'))

pve.df$variable <- factor(pve.df$variable, c(var.names, 'genetic'))

gg.plot <- function(...) {
    ggplot(...) + theme_bw() + theme(panel.background = element_blank())
}

plt.pve.stat <-
    gg.plot(pve.df, aes(x = variable, fill = cpg.data,
                        ymin = pmax(pve.mean - 2 * pve.se, 1e-4), 
                        ymax = pmin(pve.mean + 2 * pve.se, 1))) +
    geom_bar(aes(y = pve.mean), stat = 'identity', position = 'dodge') +
    geom_errorbar(position = position_dodge(.9), color = 'black', width = .3, size =2) +
    facet_wrap(~pheno, scales = 'free', nrow = 2, dir = 'h') +
    ylab('proportion of variance explained (DNAme)') +
    theme(axis.text.x = element_text(size=8, angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual('CpGs', values = c('#770000', 'gray')) +
    scale_y_continuous(trans = 'sqrt', breaks = c(1e-3, 1e-2, 0.025, 0.05, 0.1, 0.2))

ggsave(filename = 'fig_pve_stat_tot.pdf', plot = plt.pve.stat, width = 8, height = 5)

################################################################
## show double phenotype stuff
pve.tab$variable <- factor(pve.tab$variable, rev(c(var.names, 'genetic')))

pve.tab$pheno <- factor(pve.tab$pheno,
                        c('NP', 'NFT', 'Cog', 'NP+NFT', 'NP+Cog', 'NFT+Cog', 'NP+NFT+Cog'))

pve.tab.mult <- pve.tab %>% filter(pheno %in% c('NP+NFT', 'NP+Cog', 'NFT+Cog', 'NP+NFT+Cog'))

temp <- pve.tab.mult %>% spread(key = variable, value = pve) %>%
    select(-pheno, -total, -MF)

cg.list <- temp$cg
cg.order <- row.order(temp %>% select(-cg) %>% as.matrix())
cg.ordered <- cg.list[cg.order]

pve.tab.mult$cg <- factor(pve.tab.mult$cg, cg.ordered)

plt.pve.heatmap <-
    gg.plot(pve.tab.mult, aes(x = cg, y = variable, fill = pve)) +
    geom_tile() + xlab('CpGs') +
    facet_grid(. ~ pheno, scales = 'free', space = 'free') +
    scale_fill_gradientn('PVE',
                         colors = c('white', 'gray40', 'gray20', 'black'),
                         breaks = c(0, 0.1, 0.25, 0.5), trans = 'sqrt') +
    theme(legend.position = 'bottom',
          legend.text = element_text(size=6),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=8, hjust = 1, vjust = 1))

ggsave(filename = 'fig_pve_stat_heatmap.pdf', plot = plt.pve.heatmap, width = 8, height = 3)

################################################################
## make tables for top doubly-associated CpGs

gene.cols <- c('cg', 'chr', 'cg.loc', 'ensg', 'hgnc', 'tss', 'tes', 'strand', 'dist')
genes.tab <- read_tsv('nearest.genes.gz', col_names = gene.cols)
pve.tab.mult$cg <- as.character(pve.tab.mult$cg)
pve.tab.mult$variable <- factor(pve.tab.mult$variable, c(var.names, 'genetic'))

mwas.compact <- mwas.sig %>% rename(variable=pheno) %>%
    select(cg, variable, theta, theta.se, lodds) %>%
    mutate(theta = signif(theta, 2), theta.se = signif(theta.se, 2),
           pip = signif(1/(1+exp(-lodds)), 2))

mwas.compact.formatted <- mwas.compact %>%
    mutate(stat = theta %&&% ' (' %&&% signif(theta.se, 1) %&&% ', ' %&&% signif(pip, 1) %&&% ')') %>%
    select(-theta, -theta.se, -pip, -lodds) %>%
    mutate(variable = variable %&&% '_effect')

mwas.compact.formatted$variable <- factor(mwas.compact.formatted$variable,
                                          c('NP', 'NFT', 'Cog') %&&% '_effect')

mwas.compact.tab <- mwas.compact.formatted %>%
    spread(key = variable, value = stat, fill = '.')

pve.mult.genes <- 
    pve.tab.mult %>%
    left_join(genes.tab, by = 'cg') %>%
    na.omit() %>%
    rename(gene = hgnc) %>%
    filter(dist < 1) %>%
    select(cg, chr, cg.loc, gene, pheno, variable, pve) %>%
    group_by(cg, chr, cg.loc, pheno, variable) %>%
    summarize(pve = signif(mean(pve), 2), gene = paste(gene, collapse = '|')) %>%
    mutate(pve = signif(pve, 2)) %>%
    spread(key = variable, value = pve)


pve.z.tab <- pve.tab.mult %>% mutate(variable = as.character(variable)) %>%
    left_join(mwas.sig %>% rename(variable = pheno)) %>%
    na.omit()

.aes <- aes(x = theta/theta.se, y = pve, size = 1/(1+exp(-lodds)))

plt.effect.dir <-
    gg.plot() +
    geom_point(data = pve.z.tab %>% filter(p.val < cutoff), .aes, pch = 21) +
    scale_y_sqrt() +
    scale_size_continuous('PIP', breaks = c(0.6, 0.7, 0.8, 0.9, 0.95), range = c(0, 3)) +
    facet_wrap(~variable) + xlab('effect size z-score') + ylab('DNAme variance explained')

ggsave(filename = 'fig_mult_pheno_direction.pdf', plot = plt.effect.dir, width = 8, height = 3)

################################################################
mult.tab.out.top <- 
    pve.mult.genes %>% mutate(AD = NP + NFT + Cog) %>% arrange(desc(AD)) %>% head(30) %>%
    arrange(pheno, chr, cg.loc) %>% as.data.frame() %>%
    rename(location = cg.loc) %>%
    mutate(location = format(location, big.mark = ',')) %>%
    select(-NeuN, -PMI, -genetic, -AD) %>%
    left_join(mwas.compact.tab) %>% 
    as.data.frame()


mult.tab.out <- 
    pve.mult.genes %>% mutate(AD = NP + NFT + Cog) %>% arrange(desc(AD)) %>%
    arrange(pheno, chr, cg.loc) %>% as.data.frame() %>%
    rename(location = cg.loc) %>%
    mutate(location = format(location, big.mark = ',')) %>%
    select(-NeuN, -PMI, -genetic, -AD) %>% 
    left_join(mwas.compact.tab) %>% 
    as.data.frame()


## tables with 30 highest PVE genes
tab.ret <- pandoc.table.return(mult.tab.out.top,
                               caption = 'Top CpGs associated with multiple phenotypes',
                               style = 'simple', digits = 2,
                               split.tables = 200)

cat(tab.ret, file = 'table_pve_mult_cpgs.md')



################################################################
write_tsv(mult.tab.out.top, path = gzfile('table_pve_mult_top.txt.gz'))
write_tsv(mult.tab.out, path = gzfile('table_pve_mult.txt.gz'))
write.xlsx(mult.tab.out.top, file = 'table_pve_mult.xlsx', sheetName = 'Top 30', row.names = FALSE)
write.xlsx(mult.tab.out, file = 'table_pve_mult.xlsx', sheetName = 'All significant', row.names = FALSE, append = TRUE)

write_tsv(pve.tab, path = gzfile('table_pve.txt.gz'))
