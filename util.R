options(stringsAsFactors = FALSE)
library(readr)
library(feather)
library(reshape2)
library(dplyr)

`%c%` <- function(mat, pos) mat[, pos, drop = FALSE]

`%&&%` <- function(a, b) paste(a, b, sep = '')

melt.mwas.tab <- function(data) {

    cg.key <- c('cg', 'chr', 'cg.loc')
    phenotypes <- read.table('phenotype.names.txt')[, 1]
    n.pheno <- length(phenotypes)

    theta.cols <- 3 + 1:n.pheno
    theta.var.cols <- 3 + n.pheno + 1:n.pheno
    lodds.cols <- 3 + 2*n.pheno + 1:n.pheno

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
    
    ret <- take.pheno.cols(data, 'theta', theta.cols) %>%
        left_join(take.pheno.cols(data, 'theta.var', theta.var.cols), by = c(cg.key, 'pheno')) %>%
            left_join(take.pheno.cols(data, 'lodds', lodds.cols), by = c(cg.key, 'pheno')) %>%
                dplyr::rename(theta.se = theta.var) %>% dplyr::mutate(theta.se = sqrt(theta.se))
    
    return(ret)
}

