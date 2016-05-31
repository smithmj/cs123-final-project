#!/bin/bash

invcf=$1
prefix=$2

annovarpath=/cs123-final-project/software/annovar
humandb=/cs123-final-project/software/annovar/humandb

# Convert vcf to annovar format 
$annovarpath/convert2annovar.pl -format vcf4 $invcf -outfile ${prefix}.annovar

# Annotate the variants

# Annotate with refgene info
$annovarpath/annotate_variation.pl --geneanno --buildver hg19 -dbtype gene ${prefix}.avinput $humandb

# Annotate with GWAS info
$annovarpath/annotate_variation.pl -regionanno --buildver hg19 -dbtype gwascatalog ${prefix}.avinput $humandb

# Annotate with dbSNP info
$annovarpath/annotate_variation.pl -regionanno --buildver hg19 -dbtype gwascatalog ${prefix}.avinput $humandb


$annovarpath/table_annovar.pl ${prefix}.annovar $humandb --buildver hg19 -protocol refGene,snp135,gwascatalog,ljb2_all -operation g,f,r,f
