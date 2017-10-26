## store all the summary data as feather format
library(readr)
library(feather)
library(reshape2)
library(dplyr)
options(stringsAsFactors = FALSE)

`%c%` <- function(mat, pos) mat[, pos, drop = FALSE]


load('summary.rdata')
phenotypes <- read.table('phenotype.names.txt')[, 1]
n.pheno <- length(phenotypes)

colnames(gauss.boot) <- c('cg', 'chr', 'cg.loc', 'theta', 'theta.var', 'lodds', 'lodds.var', 'pheno')
gauss.boot$cg <- as.character(gauss.boot$cg)
gauss.boot$pheno <- as.character(gauss.boot$pheno)

boot.tab <- gauss.boot %>% 
    dplyr::rename(theta.se = theta.var) %>% dplyr::mutate(theta.se = sqrt(theta.se)) %>%
    dplyr::rename(lodds.se = lodds.var) %>% dplyr::mutate(lodds.se = sqrt(lodds.se)) %>%
    dplyr::rename(boot.theta = theta, boot.theta.se = theta.se,
                  boot.lodds = lodds, boot.lodds.se = lodds.se)

write_feather(gauss.boot, path = 'table_bootstrap.ft')

theta.cols <- 3 + 1:n.pheno
theta.var.cols <- 3 + n.pheno + 1:n.pheno
lodds.cols <- 3 + 2*n.pheno + 1:n.pheno

cg.key <- c('cg', 'chr', 'cg.loc')

take.pheno.cols <- function(.dat, val.name, .cols) {
    .temp <- .dat %c% c(1:3, .cols)
    colnames(.temp) <- .var.names <- c(cg.key, phenotypes)
    ret <- melt(.temp, id.vars = cg.key, variable.name = 'pheno',
                  value.name = val.name, as.is = TRUE, factorsAsStrings = TRUE) %>%
                      as.data.frame()

    ret$cg <- as.character(ret$cg)
    ret$pheno <- as.character(ret$pheno)
    ret$chr <- as.integer(ret$chr)
    ret$cg.loc <- as.integer(ret$cg.loc)

    rm(.temp)
    return(ret)
}

mwas.tab <- take.pheno.cols(gauss.data, 'theta', theta.cols) %>%
    left_join(take.pheno.cols(gauss.data, 'theta.var', theta.var.cols), by = c(cg.key, 'pheno')) %>%
    left_join(take.pheno.cols(gauss.data, 'lodds', lodds.cols), by = c(cg.key, 'pheno')) %>%
    dplyr::rename(theta.se = theta.var) %>% dplyr::mutate(theta.se = sqrt(theta.se)) %>%
    left_join(boot.tab, by = c(cg.key, 'pheno'))

write_feather(mwas.tab, path = 'table_mwas.ft')
        
        

