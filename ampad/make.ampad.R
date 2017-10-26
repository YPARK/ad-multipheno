library(xlsx)
library(dplyr)
library(readr)
library(feather)
library(biomaRt)
library(valr)

genes <- read.xlsx('AMPAD.xlsx', sheetIndex = 1)
genes <- genes[, 1:2] %>% na.omit()
genes <- genes[-1, ]
colnames(genes) <- c('hgnc', 'ensg')

grch37 <- useEnsembl(biomart = 'ensembl',
                     dataset = 'hsapiens_gene_ensembl',
                     GRCh = 37)

bm.attr <- c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name',
             'transcript_start', 'transcript_end')
bm.tab <- getBM(attributes = bm.attr,
                filters = 'ensembl_gene_id',
                values = genes$ensg,
                mart = grch37)

genes.annot <- bm.tab %>%
    group_by(chromosome_name, hgnc_symbol, ensembl_gene_id) %>%
    summarize(tss = min(transcript_start),
              tes = min(transcript_end)) %>%
    mutate(start = tss - 1e5, end = tes + 1e5) %>%
    rename(chrom = chromosome_name,
           ensg = ensembl_gene_id,
           hgnc = hgnc_symbol)


mwas.tab <- read_feather('table_mwas.ft')
mwas.boot.tab <- mwas.tab %>% na.omit() %>%
    mutate(p.val = pnorm((lodds -boot.lodds)/boot.lodds.se, lower.tail = FALSE))

n.tot <- mwas.tab %>% filter(pheno == 'NP') %>% nrow()
cutoff <- 0.05 / n.tot / 3

mwas.sig <- mwas.boot.tab %>% filter(p.val < cutoff, lodds > 0) %>%
    dplyr::mutate(chrom = chr, start = cg.loc - 1, end = cg.loc)


mwas.bed <- bed_sort(mwas.sig)
gene.bed <- bed_sort(genes.annot)




## bedr(
