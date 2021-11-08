# Analyze RNAseq data from grapevine leaves with the intention of relating gene expression to methylation patterns.

## fastqc and trim reads

Leaf

	cd /data/vondras/methylome.project/rnaseq/fq/original
	ls *_leaf.fq.gz | sed 's/_leaf.fq.gz//g' | sed 's/_R.*//g'| sort | uniq > leaf.list
	bash launcher.sh leaf.list leaf.sh 3
	cat leaf.list.*out leaf.list.*err

	#!/bin/bash

	WD=/data/vondras/methylome.project/rnaseq
	FQ=$WD/fq
	FQC=$WD/fastqc
	TRIMMO=/data/Assembly_tools/Tools/Trimmomatic-0.36/trimmomatic-0.36.jar
	
	for line in $( cat $1 ); do
        	echo ${line}
	
		# FastQC, PreTrim
        	cd $FQ/original
        	fastqc $FQ/original/${line}_R1_leaf.fq.gz $FQ/original/${line}_R2_leaf.fq.gz

		# Trim and FastQC
		trim_galore --length 80 --fastqc --gzip --paired $FQ/original/${line}_R1_leaf.fq.gz $FQ/original/${line}_R2_leaf.fq.gz --output_dir $FQ/trimmed
		rm $FQ/original/${line}_R1_leaf.fq.gz $FQ/original/${line}_R2_leaf.fq.gz
	done
	mv $FQ/trimmed/*fastqc* $FQC/trimmed/
	mv $FQ/original/*fastqc* $FQC/original/
	
Berry
	
	cd /data/vondras/methylome.project/rnaseq/fq/original
	ls *_berry.fq.gz > berry.list
	bash launcher.sh berry.list berry.sh 8
	cat berry.list.*out berry.list.*err

	#!/bin/bash

	WD=/data/vondras/methylome.project/rnaseq
	FQ=$WD/fq
	FQC=$WD/fastqc
	TRIMMO=/data/Assembly_tools/Tools/Trimmomatic-0.36/trimmomatic-0.36.jar
	
    	for line in $( cat $1 ); do
        	echo ${line}
	
		# FastQC, PreTrim
        	cd $FQ/original
        	fastqc $FQ/original/${line}

		# Trim and FASTQC
		trim_galore --length 80 --fastqc --gzip $FQ/original/${line} --output_dir $FQ/trimmed	 
		rm $FQ/pretrim/${line}	
	done
	mv $FQ/trimmed/*fastqc* $FQC/trimmed/
	mv $FQ/original/*fastqc* $FQC/original/


**Check trimming before and after here**


## HISAT2 & Salmon

[Salmon FAQ](https://combine-lab.github.io/salmon/)

[Salmon Documentation](https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-mapping-based-mode)





## Make indexes for hisat2 and salmon. Make genome bed annotation.

	cd /data/vondras/methylome.project/rnaseq/index

	#!/bin/bash

	INDEX=/data/vondras/methylome.project/rnaseq/index
	GENOME=/DATA7/All_vinifera_genomes/Genomes/VITVvi_vCabSauv08_v1.1.pseudomolecules.all.fasta
	ANNOTATION=/DATA7/All_vinifera_genomes/Annotations/VITVvi_vCabSauv08_v1.1.pseudomolecules.all.gff3

	GFFREAD=/DATA/users/aminio/Assembly_tools/Tools/cufflinks-2.2.1.Linux_x86_64/gffread
	GTFBED=/DATA7/Resources/Tools/augustus/scripts/gtf2bed.pl
	SALMON=/DATA7/Resources/Tools/salmon-1.5.1_linux_x86_64/bin/salmon	

	# hisat2 indexes
	hisat2-build $GENOME VITVvi_vCabSauv08_v1.1.pseudomolecules.all
	
	# index genome
	ln -s $GENOME . ; samtools faidx VITVvi_vCabSauv08_v1.1.pseudomolecules.all.fasta # Index genome
		
	# Create mRNA fasta using genome and genome annotation
	$GFFREAD -w VITVvi_vCabSauv08_v1.1.pseudomolecules.all.mRNA.salmon.fasta -g VITVvi_vCabSauv08_v1.1.pseudomolecules.all.fasta $ANNOTATION 	
	
	# make bed annotation
	$GFFREAD -o - -T $ANNOTATION > VITVvi_vCabSauv08_v1.1.pseudomolecules.all.salmon.gtf
	$GTFBED < VITVvi_vCabSauv08_v1.1.pseudomolecules.all.salmon.gtf > VITVvi_vCabSauv08_v1.1.pseudomolecules.all.salmon.bed ; rm VITVvi_vCabSauv08_v1.1.pseudomolecules.all.salmon.gtf
		
	# salmon indexes and decoys only used if salmon both maps and quatifies
	grep "^>" $GENOME | sed 's/>//g' > decoys.txt
 	TRANSCRIPTOME=/data/vondras/methylome.project/rnaseq/index/VITVvi_vCabSauv08_v1.1.pseudomolecules.all.mRNA.salmon.fasta
	cat $TRANSCRIPTOME $GENOME | pigz -9vc > VITVvi_vCabSauv08_v1.1.pseudomolecules.mrna.genome.salmon.fasta.gz
	$SALMON index -t VITVvi_vCabSauv08_v1.1.pseudomolecules.mrna.genome.salmon.fasta.gz -d decoys.txt -p 12 -i VITVvi_vCabSauv08_v1.1.pseudomolecules.mrna.genome.salmon.index --gencode --keepDuplicates > salmon_index.log 2> salmon_index.err
	



## Use salmon to map and quantify

	cd /data/vondras/methylome.project/rnaseq/align/
	nohup bash ap4.sh &> ap4.out &

	#!/bin/bash
	
	SALMON=/DATA7/Resources/Tools/salmon-1.5.1_linux_x86_64/bin/salmon
	SINDEX=/data/vondras/methylome.project/rnaseq/index/VITVvi_vCabSauv08_v1.1.pseudomolecules.mrna.genome.salmon.index
	FQ=/data/vondras/methylome.project/rnaseq/fq/trimmed
	ALIGN=/data/vondras/methylome.project/rnaseq/align
	
	# Leaf (paired-end)
	for line in $( ls $FQ/*_leaf* | sed -e 's/.*.\///g' -e 's/_R.*//g' | sort | uniq ); do
		$SALMON quant -i $SINDEX -l A -1  $FQ/${line}_R1_leaf_val_1.fq.gz -2 $FQ/${line}_R2_leaf_val_2.fq.gz --validateMappings -o ${line}_leaf_salmon_direct.salmon_quant -p 12 --seqBias
	done

	# Berry (single-end)
	for line in $( ls $FQ/*_berry_trimmed.fq.gz | sed 's/_berry_trimmed.fq.gz//g' | sed 's/.*.\///g' | sort | uniq ); do
		$SALMON quant -i $SINDEX -l A -r $FQ/${line}_berry_trimmed.fq.gz --validateMappings -o ${line}_berry_salmon_direct.salmon_quant -p 12 --seqBias
	done

check alignment rate

	grep 'Logs will be written to \|Mapping rate' salmon.out | paste - - | sed 's/ /\t/g' | cut -f 6,14 | sed 's/_salmon_direct.salmon_quant\/logs//g' 


## Combine counts tables from different libraries into single table

	grep '.*' *salmon_direct.salmon_quant/*.sf | sed -e 's/\:/\t/g' -e 's/_salmon_direct.salmon_quant\/quant.sf//g' | grep -v "NumReads" | sed '1ilibrary\ttranscript\tlength\teffective_length\ttpm\tnumber_reads' > combined_transcripts.quant.sf





## Keep haplotypes separate, but merge transcripts belonging to same gene in one haplotype using tximport and prep for DESeq2

[Salmon documentation recommends using tximport](https://combine-lab.github.io/salmon/getting_started/#after-quantification)

[tximport documentation](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)


	cd /data/vondras/methylome.project/rnaseq/align; 	
	
	cut -f 1 combined_transcripts.quant.sf | sed '1d' | awk '{print $1"_salmon_direct.salmon_quant/quant.sf" "\t" $0}' | sort | uniq | sed '1ifile_path\tsample' > tximport_filenames.txt
	
	zcat /data/vondras/methylome.project/rnaseq/index/VITVvi_vCabSauv08_v1.1.pseudomolecules.mrna.genome.salmon.fasta.gz | grep ">" | grep "gene="| sed -e 's/>//g' -e 's/gene=//g' | cut -f 1,2 -d ' ' | sort | uniq | sed '1itranscript\tgene' > tx2gene
	
	R
	
	library(tximport); library(readr)
	
	# make sample key pointing at salmon quantitation files
	files_samples <- read.table( "tximport_filenames.txt", header = TRUE)
	files <- files_samples$file_path
	names(files) <- files_samples$sample
	all(file.exists(files))

	# make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID
	tx2gene <- read.table( "tx2gene", header = TRUE)
	
	# import and convert to gene-level info
	txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
 
 	dim(txi$counts)
	names(txi)
 	txi$counts[1:3,1:3]
 	txi$counts[57000:57003,1:3]

	summary(txi$counts[grepl("Primary", row.names(txi$counts)),])
	summary(txi$counts[grepl("Haplotig", row.names(txi$counts)),])
	
	tmp <-summary(txi$counts[grepl("Primary", row.names(txi$counts)),])
	tmp <- gsub("^.*.\\:", "", tmp[grepl("Median", tmp)])
	summary(as.numeric(gsub(" ", "", tmp)))
	#   Min.   1st Qu.  Median    Mean    3rd Qu.    Max. 
	#  93,455  359,114  468,380  495,665  625,157    1,335,393  # Max
	# 3.000   6.000   7.240   9.687  10.000  33.000 
		
	tmp <-summary(txi$counts[grepl("Haplotig", row.names(txi$counts)),])
	tmp <- gsub("^.*.\\:", "", tmp[grepl("Median", tmp)])
	summary(as.numeric(gsub(" ", "", tmp)))
	#   Min.   1st Qu.  Median   Mean     3rd Qu.    Max. 
	# 214,998   440,539  540,612  596,780  667,084   1,678,706  # Max
	# 4.00    7.50    9.40   12.89   14.60   44.40 
	
	rm(tmp)
	
	library(DESeq2)
	
	sampleTable <- data.frame(
		tissue= gsub(".*._", "", names(data.frame(txi$counts))), 
		cultivar= gsub("\\..*.", "", names(data.frame(txi$counts)))
		)
	rownames(sampleTable) <- names(data.frame(txi$counts))
	row.names(sampleTable) <- gsub("\\.", "-", row.names(sampleTable))
	table( row.names(sampleTable)==colnames(txi$counts))
	sampleTable$condition <- as.factor(paste(sampleTable$cultivar, sampleTable$tissue, sep= "_"))	
	dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
	
Save tables

	save(txi, files, tx2gene, sampleTable, dds, file="tximport_gene_level_salmon_quant.RData")
	

