all: table_mwas_significant.xlsx \
  nearest.genes.gz


## 1. assign genes by mediation analysis
table_mwas_significant.xlsx: make.mwas-genes.R m2t.txt.gz table_bootstrap.ft table_mwas.ft
	Rscript --vanilla $<

table_mwas.txt.gz: table_mwas_significant.xlsx
table_mwas_significant_gene.txt.gz: table_mwas_significant.xlsx
table_mwas_significant.txt.gz: table_mwas_significant.xlsx


## 2. more extensively assign nearest CpGs to codig genes
nearest.genes.gz: temp/coding.genes.bed.gz
	cat table_mwas.txt.gz | gzip -d | tail -n +2 | awk -F'\t' '{ bed = $$2 FS ($$3 - 1) FS $$3 FS $$1; BED[bed] ++  } END { for(b in BED) print b }' | sort -k1,1 -k2,2n | bedtools closest -a stdin -b $< -t all -d | awk -F'\t' '{ print $$4 FS $$1 FS $$3 FS $$8 FS $$10 FS $$6 FS $$7 FS $$9 FS $$11 }' | gzip > $@

temp/coding.genes.bed.gz: coding.genes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $< | gzip -d | awk -F'\t' '{ print $$1 FS $$2 FS $$3 FS $$5 FS $$4 FS $$6 }' | sort -k1,1 -k2,2n | gzip > $@


## convert large rdata to easily loadable feather format
table_mwas.ft: make.feather.R
	Rscript --vanilla $<

table_bootstrap.ft: table_mwas.ft
