#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)

library(piano)
library(cba)
library(dplyr)
library(tidyr)
library(readr)
library(feather)

`%&&%` <- function(a, b) paste(a, b, sep = '')

rm.header <- function(s) paste(strsplit(s, '[_]')[[1]][-1], collapse = '_')

.write.tab <- function(...) {
    write_tsv(..., col_names = TRUE)
}

pheno.names <- c('NP', 'NFT', 'Cog')
pheno.colors <- c('#FFAA00','#AA55FF','#00DD99')
pheno.shapes <- c(21,22,25)

var.names <- read_tsv('phenotype.names.txt', col_names = FALSE)[,1]

cpg2gene <- read_tsv('nearest.genes.gz', col_names = FALSE)
colnames(cpg2gene) <- c('cg', 'chr', 'cg.loc', 'ensg', 'hgnc', 'tss', 'tes', 'strand', 'dist')

mwas.tab <- read_feather('table_mwas.ft')

mwas.boot.tab <- mwas.tab %>% na.omit() %>%
    mutate(p.val = pnorm((lodds -boot.lodds)/boot.lodds.se, lower.tail = FALSE))

n.tot <- mwas.tab %>% filter(pheno == 'NP') %>% nrow()
cutoff <- 0.05 / n.tot / 3

mwas.sig <- mwas.boot.tab %>%
    filter(p.val < cutoff)

################################################################
take.gsa <- function(gsc.file, .pheno, n.perm = 5e4, dist.cutoff = 1e4, full.list = FALSE) {

    gsc <- loadGSC(gsc.file)

    if(full.list) {
        gene.tab <- mwas.tab %>%
            filter(pheno == .pheno)
    } else {
        sig.cg <- mwas.sig %>%
            dplyr::select(cg) %>% unique()
        gene.tab <- mwas.tab %>%
            filter(pheno == .pheno, cg %in% sig.cg$cg)
    }

    gene.tab <- gene.tab %>%
        left_join(cpg2gene %>% dplyr::filter(dist < dist.cutoff), by = 'cg') %>%
            na.omit() %>%
                group_by(hgnc) %>%
                    slice(which.max(lodds)) %>%
                        mutate(z = theta / (theta.se + 1e-4)) %>%
                            filter(is.finite(z))

    gsc.stat <- gene.tab$z
    names(gsc.stat) <- gene.tab$hgnc

    gsa <- runGSA(geneLevelStats = gsc.stat, gsc = gsc, geneSetStat = 'median',
                  nPerm = n.perm, gsSizeLim = c(10, 500), verbose = TRUE)

    .tab <- GSAsummaryTable(gsa, save = FALSE)
    gsa.tab <- data.frame(gs = .tab[, 'Name'],
                          p.val.up = .tab[, 'p (dist.dir.up)'],
                          q.val.up = .tab[, 'p adj (dist.dir.up)'],
                          p.val.dn = .tab[, 'p (dist.dir.dn)'],
                          q.val.dn = .tab[, 'p adj (dist.dir.dn)'],
                          pheno = .pheno) %>%
                              dplyr::mutate(gs = sapply(as.character(gs), rm.header)) %>%
                                  arrange(p.val.up)

    return(gsa.tab)
}

take.pair.gsa <- function(gsc.file, pheno.pair, n.perm = 5e4, dist.cutoff = 1e4) {

    gsc <- loadGSC(gsc.file)

    sig.cg <- mwas.sig %>%
        dplyr::filter(pheno %in% pheno.pair) %>%
            dplyr::select(cg) %>% unique()

    gene.pair.tab <-
        mwas.tab %>% filter(cg %in% sig.cg$cg) %>%
            dplyr::filter(pheno %in% pheno.pair) %>%
                dplyr::select(cg, pheno, lodds) %>%
                    left_join(cpg2gene %>% dplyr::filter(dist < dist.cutoff), by = 'cg') %>%
                        dplyr::group_by(hgnc, pheno) %>%
                            dplyr::summarize(lodds = max(lodds)) %>%
                                tidyr::spread(key = 'pheno', value = 'lodds') %>%
                                    na.omit()

    gsc.stat <- gene.pair.tab[, pheno.pair[1]] - gene.pair.tab[, pheno.pair[2]]
    gsc.stat <- gsc.stat[, 1]
    names(gsc.stat) <- gene.pair.tab$hgnc

    gsa <- runGSA(geneLevelStats = gsc.stat, gsc = gsc, geneSetStat = 'median',
                  nPerm = n.perm, gsSizeLim = c(10, 500), verbose = TRUE)

    .tab <- GSAsummaryTable(gsa, save = FALSE)

    gsa.tab <- data.frame(.tab[, 'Name'],
                          .tab[, 'p (dist.dir.up)'],
                          .tab[, 'p adj (dist.dir.up)'],
                          .tab[, 'p (dist.dir.dn)'],
                          .tab[, 'p adj (dist.dir.dn)'],
                          pheno.pair[1],
                          pheno.pair[2])

    colnames(gsa.tab) <- c('gs', 'p.val.' %&&% pheno.pair[1], 'q.val.' %&&% pheno.pair[1],
                           'p.val.' %&&% pheno.pair[2], 'q.val.' %&&% pheno.pair[2],
                           'pheno.1', 'pheno.2')

    gsa.tab <- gsa.tab %>%
        dplyr::mutate(gs = sapply(as.character(gs), rm.header))

    return(gsa.tab)
}



################################################################
dir.create('pathway', showWarnings = FALSE)

for(.pheno in pheno.names) {
    gsc.file <- 'genesets/c2.cp.kegg.v6.0.symbols.gmt'
    out.file <- 'pathway/' %&&% 'KEGG_' %&&% .pheno %&&% '.txt.gz'
    if(!file.exists(out.file)) {
        gsa.tab <- take.gsa(gsc.file, .pheno, n.perm = 5e4, dist.cutoff = 1e4)
        .write.tab(gsa.tab, path = gzfile(out.file))
    }

    gsc.file <- 'genesets/c2.cp.reactome.v6.0.symbols.gmt'
    out.file <- 'pathway/' %&&% 'REACTOME_' %&&% .pheno %&&% '.txt.gz'
    if(!file.exists(out.file)) {
        gsa.tab <- take.gsa(gsc.file, .pheno, n.perm = 5e4, dist.cutoff = 1e4)
        .write.tab(gsa.tab, path = gzfile(out.file))
    }

    gsc.file <- 'genesets/c2.cp.kegg.v6.0.symbols.gmt'
    out.file <- 'pathway/' %&&% 'full-KEGG_' %&&% .pheno %&&% '.txt.gz'
    if(!file.exists(out.file)) {
        gsa.tab <- take.gsa(gsc.file, .pheno, n.perm = 5e4, dist.cutoff = 1e4, full.list = TRUE)
        .write.tab(gsa.tab, path = gzfile(out.file))
    }

    gsc.file <- 'genesets/c2.cp.reactome.v6.0.symbols.gmt'
    out.file <- 'pathway/' %&&% 'full-REACTOME_' %&&% .pheno %&&% '.txt.gz'
    if(!file.exists(out.file)) {
        gsa.tab <- take.gsa(gsc.file, .pheno, n.perm = 5e4, dist.cutoff = 1e4, full.list = TRUE)
        .write.tab(gsa.tab, path = gzfile(out.file))
    }
}

for(.p1 in pheno.names) {
    for(.p2 in pheno.names) {
        if(.p1 != .p2) {
            .pp <- c(.p1, .p2)
            .pp.str <- paste(.pp, collapse = '-')

            gsc.file <- 'genesets/c2.cp.kegg.v6.0.symbols.gmt'
            out.file <- 'pathway/' %&&% 'KEGG_' %&&% .pp.str %&&% '.txt.gz'

            if(!file.exists(out.file)) {
                gsa.tab <- take.pair.gsa(gsc.file, pheno.pair = .pp, n.perm = 5e4, dist.cutoff = 1e4)
                .write.tab(gsa.tab, path = gzfile(out.file))
            }

            gsc.file <- 'genesets/c2.cp.reactome.v6.0.symbols.gmt'
            out.file <- 'pathway/' %&&% 'REACTOME_' %&&% .pp.str %&&% '.txt.gz'
            if(!file.exists(out.file)) {
                gsa.tab <- take.pair.gsa(gsc.file, pheno.pair = .pp, n.perm = 5e4, dist.cutoff = 1e4)
                .write.tab(gsa.tab, path = gzfile(out.file))
            }
        }
    }
}
