## store all the summary data as feather format
source('util.R')

load('summary-G.rdata')

colnames(gauss.boot) <- c('cg', 'chr', 'cg.loc', 'theta', 'theta.var', 'lodds', 'lodds.var', 'pheno')
gauss.boot$cg <- as.character(gauss.boot$cg)
gauss.boot$pheno <- as.character(gauss.boot$pheno)

boot.tab <- gauss.boot %>% 
    dplyr::rename(theta.se = theta.var) %>% dplyr::mutate(theta.se = sqrt(theta.se)) %>%
    dplyr::rename(lodds.se = lodds.var) %>% dplyr::mutate(lodds.se = sqrt(lodds.se)) %>%
    dplyr::rename(boot.theta = theta, boot.theta.se = theta.se,
                  boot.lodds = lodds, boot.lodds.se = lodds.se)

cg.key <- c('cg', 'chr', 'cg.loc')

mwas.tab <- melt.mwas.tab(gauss.data) %>%
    left_join(boot.tab, by = c(cg.key, 'pheno'))

write_feather(mwas.tab, path = 'table_mwas.ft')
        
        

