all: AMPAD-mwas.xlsx

# 1. generate gene list and intersect with MWAS

intersect.bed.gz: AMPAD.genes.bed.gz mwas.bed.gz
	bedtools intersect -a mwas.bed.gz -b AMPAD.genes.bed.gz  -wa -wb | gzip > $@

mwas.bed.gz: make.mwas.bed.R
	Rscript --vanilla make.mwas.bed.R

AMPAD.genes.bed.gz: make.gene.bed.R
	Rscript --vanilla make.gene.bed.R

# 2. take Excel table
AMPAD-mwas.xlsx: intersect.bed.gz make.report.R
	Rscript --vanilla make.report.R

