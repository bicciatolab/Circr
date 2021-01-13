File explanation:
Genome assembly (XXX) specific files
 - XXX.ensGene.gtf	esnembl gene file; downloaded from UCSC
 - XXX.fa.gz		genome fasta file; downloaded from UCSC
 - XXX.rRNA.bed		rRNA coordinates; downloaded from UCSC
 - XXX.AGO.bed		collection of AGO peaks; multiple sources
 - XXX.INT.bed		collection of validated miRNA-RNA interactions; multiple sources

Common files
 - mature.fa		miRBase v22 miRNA mature sequences for all organisms
 - miR_family_info.txt	input file required for TargetScan
 - org.3lc.tax.txt	table connecting 3 letter code (e.g. mmu) and taxonomy ID (needed for TargetScan)

Folder organization

support_files
   |_ fruitfly
   |	|_ dm3
   |	|   |_ dm3.ensGene.gtf
   |	|   |_ dm3.fa.gz
   |	|   |_ dm3.rRNA.bed
   |	|
   |	|_ dm6
   |	    |_ dm6.ensGene.gtf
   |	    |_ dm6.fa.gz
   |	    |_ dm6.rRNA.bed
   |
   |_ human
   |	|_ hg19
   |	|   |_ hg19.ensGene.gtf
   |	|   |_ hg19.fa.gz
   |	|   |_ hg19.rRNA.bed
   |	|_ hg38
   |	    |_ hg38.ensGene.gtf
   |	    |_ hg38.fa.gz
   |	    |_ hg38.rRNA.bed
   |
   |_ miRNA
   |	|_ mature_v22.fa
   |	|_ miR_family_info.txt
   |	|_ org.3lc.tax.txt
   |
   |_ mouse
   |	|_ mm9
   |	|   |_ mm9.ensGene.gtf
   |	|   |_ mm9.fa.gz
   |	|   |_ mm9.rRNA.bed
   |	|
   |	|_ mm10
   |	    |_ mm10.AGO.bed
   |	    |_ mm10.ensGene.gtf
   |	    |_ mm10.fa.gz
   |	    |_ mm10.INT.bed
   |	    |_ mm10.rRNA.bed
   |
   |_ worm
	|_ ce10
	|   |_ ce10.ensGene.gtf
	|   |_ ce10.fa.gz
	|   |_ ce10.rRNA.bed
	|_ ce11
	    |_ ce11.ensGene.gtf
	    |_ ce11.fa.gz
	    |_ ce11.rRNA.bed

