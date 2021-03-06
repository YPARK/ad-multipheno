all: table_mwas.txt.gz \
  table_replication_formatted.xlsx \
  nearest.genes.gz \
  table_pve.txt.gz \
  table_pathways.md

## 1. Take significant CpGs while converting large rdata to easily
## loadable feather format
table_mwas.txt.gz: make.table-significant.R table_mwas.ft
	Rscript --vanilla $<

table_mwas_significant.txt.gz: table_mwas.txt.gz

table_mwas.ft: make.feather.R
	Rscript --vanilla $<

table_bootstrap.ft: table_mwas.ft

## 2. more extensively assign nearest CpGs to codig genes
nearest.genes.gz: temp/coding.genes.bed.gz
	cat table_mwas.txt.gz | gzip -d | tail -n +2 | awk -F'\t' '{ bed = $$2 FS ($$3 - 1) FS $$3 FS $$1; BED[bed] ++  } END { for(b in BED) print b }' | sort -k1,1 -k2,2n | bedtools closest -a stdin -b $< -t all -d | awk -F'\t' '{ print $$4 FS $$1 FS $$3 FS $$8 FS $$10 FS $$6 FS $$7 FS $$9 FS $$11 }' | gzip > $@

temp/coding.genes.bed.gz: coding.genes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $< | gzip -d | awk -F'\t' '{ print $$1 FS $$2 FS $$3 FS $$5 FS $$4 FS $$6 }' | sort -k1,1 -k2,2n | gzip > $@

## 3. take nicely replicated CpGs
table_replication.txt.gz: table_replication_formatted.xlsx
table_top_replication_cpgs.md: table_replication_formatted.xlsx
table_replication.xlsx: table_replication_formatted.xlsx
figure-replication/fig_dejager.pdf: table_replication_formatted.xlsx
table_replication_formatted.xlsx: make.table-replication.R table_mwas.txt.gz nearest.genes.gz dejager.cpgs.txt
	Rscript --vanilla $<

## 4. Proportion of variance explained
fig_pve_stat_tot.pdf: table_pve.txt.gz
fig_pve_stat_heatmap.pdf: table_pve.txt.gz
table_pve_tot.txt.gz: table_pve.txt.gz
table_pve.txt.gz: make.table-pve.R table_mwas_significant.txt.gz
	Rscript --vanilla $<


