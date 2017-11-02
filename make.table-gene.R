source('util.R')
library(ggplot2)
library(pander)
library(xlsx)
library(dplyr)
library(tidyr)
library(qvalue)
library(goseq)

.files <- 'genetest/genetest-chr' %&&% 1:22 %&&% '.txt.gz'

gene.test.cols <- c('ensg', 'chr', 'tss', 'tes', 'hgnc', 'pheno', 'T.obs',
               'T.null.max', 'T.null.95', 'T.null.50',
               'n.false', 'n.tot', 'n.cpgs')

gene.test.tab <- do.call(rbind, lapply(.files, read.table, col.names = gene.test.cols)) %>%
    mutate(p.val = (n.false + 1)/(n.tot + 1))

qobj <- qvalue(gene.test.tab$p.val, pi0.method = 'smooth', pfdr = FALSE)
gene.test.tab <- data.frame(gene.test.tab, q.val = qobj$qvalues)

gene.cols <- c('cg', 'chr', 'cg.loc', 'ensg', 'hgnc', 'tss', 'tes', 'strand', 'dist')
gene.info.tab <- read.table('nearest.genes.gz', col.names = gene.cols)

gene.length.tab <- gene.info.tab %>% select(ensg, tss, tes) %>% unique() %>%
    mutate(len = tes - tss)


n.tot <- gene.test.tab %>% select(ensg, pheno) %>% unique() %>% nrow()
cutoff <- 0.05 / n.tot

genes <- gene.test.tab$ensg %>% unique()
n.genes <- length(genes)

.pheno <- 'NFT'

de.genes <- gene.test.tab %>%
    dplyr::filter(pheno == .pheno, q.val <= 5e-2) %>%
    dplyr::select(ensg)

gene.vector <- as.integer(genes %in% de.genes$ensg)
names(gene.vector) <- genes
head(genes)

pwf <- nullp(gene.vector, 'hg19', 'ensGene')
ret.0 <- goseq(pwf, 'hg19', 'ensGene')
qobj <- qvalue(pmin(ret.0$over_represented_pvalue, 1),
               pi0.method = 'smooth', pfdr = FALSE)

ret <- data.frame(ret.0, q.val = qobj$qvalues) %>%
    dplyr::rename(p.val = over_represented_pvalue) %>%
    dplyr::select(category, p.val, q.val, ontology, term, numDEInCat, numInCat) %>%
    arrange(q.val, ontology)



temp <- getgo(de.genes$ensg, 'hg19', 'ensGene')

names(temp)

## TODO: show gene x GO 

