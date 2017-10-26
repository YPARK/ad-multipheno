library(xlsx)
library(readr)
library(dplyr)

tab <- read_tsv('intersect.bed.gz', col_names = FALSE)

colnames(tab) <- c('chr', 'cg.loc.1', 'cg.loc', 'cg', 'pheno',
                   'effect', 'effect.se', 'effect.pip', 'p.val',
                   'chr.1', 'lb', 'ub', 'ensg', 'hgnc', 'left', 'right',
                   'strand')

tab <- tab %>% mutate(effect.pip = 1/(1 + exp(-effect.pip))) %>%
    mutate(effect.pip = signif(effect.pip, 2))

.left <- which(tab$cg.loc < tab$left)
.right <- which(tab$cg.loc > tab$right)
.middle <- which(tab$cg.loc < tab$right & tab$cg.loc > tab$left)
tab$dist <- 0
tab$dist[.left] <- tab$left[.left] - tab$cg.loc[.left]
tab$dist[.right] <- tab$cg.loc[.right] - tab$right[.right]

ret <- tab %>%
    dplyr::select(hgnc, ensg, cg, pheno, p.val, effect, effect.se, effect.pip,
                  chr, left, right, cg.loc, strand,
                  dist)

ret <- ret %>% mutate(strand = ifelse(strand > 0, '+', '-'))


ret1 <- as.data.frame(ret %>% filter(dist < 1e4))
ret2 <- as.data.frame(ret)


.col <- function(...) paste(..., collapse = '|')

ret3 <- ret %>%
    arrange(chr, cg.loc) %>%
    group_by(ensg, hgnc, chr) %>%
    summarize(cg = .col(cg),
              pheno = .col(pheno),
              cg.loc = .col(cg.loc),
              effect = .col(signif(effect, 2)),
              effect.se = .col(signif(effect.se, 2)),
              effect.pip = .col(signif(effect.pip, 2))) %>%
    as.data.frame()

ret4 <- ret %>%
    group_by(cg, cg.loc, chr) %>%
    arrange(chr, cg.loc) %>%
    summarize(hgnc = .col(hgnc),
              ensg = .col(ensg),
              pheno = .col(pheno),
              effect = .col(signif(effect, 2)),
              effect.se = .col(signif(effect.se, 2)),
              effect.pip = .col(signif(effect.pip, 2))) %>%
    as.data.frame()

write.xlsx(ret1, file = 'AMPAD-mwas.xlsx', sheetName = 'dist < 10k', row.names = FALSE)
write.xlsx(ret2, file = 'AMPAD-mwas.xlsx', sheetName = 'dist < 100k', row.names = FALSE, append = TRUE)

write.xlsx(ret3, file = 'AMPAD-mwas.xlsx', sheetName = 'gene-centric', row.names = FALSE, append = TRUE)

write.xlsx(ret4, file = 'AMPAD-mwas.xlsx', sheetName = 'cpg-centric', row.names = FALSE, append = TRUE)


