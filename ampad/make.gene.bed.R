library(xlsx)
library(biomaRt)
library(dplyr)

genes <- read.xlsx('AMPAD.xlsx', sheetIndex = 1)
genes <- genes[, 1:2] %>% na.omit()
genes <- genes[-1, ]
colnames(genes) <- c('hgnc', 'ensg')

grch37 <- useEnsembl(biomart = 'ensembl',
                     dataset = 'hsapiens_gene_ensembl',
                     GRCh = 37)

bm.attr <- c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name',
             'transcript_start', 'transcript_end', 'strand')

bm.tab <- getBM(attributes = bm.attr,
                filters = 'ensembl_gene_id',
                values = genes$ensg,
                mart = grch37) %>%
    rename(chr = chromosome_name, hgnc = hgnc_symbol, ensg = ensembl_gene_id)

window <- 1e5

genes.annot <- bm.tab %>%
    group_by(chr, hgnc, ensg, strand) %>%
    summarize(left = min(c(transcript_start, transcript_end)),
              right = max(c(transcript_start, transcript_end))) %>%
    mutate(lb = left - window, ub = right + window) %>%
    select(chr, lb, ub, ensg, hgnc, left, right, strand) %>%
    unique()

write.table(genes.annot, file = gzfile('AMPAD.genes.bed.gz'),
            row.names = FALSE, col.names = FALSE,
            sep = '\t', quote = FALSE)
