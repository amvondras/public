# Get data

	WD=/data/vondras/methylome.project
	FQ=$WD/methylation/fq
	cd $FQ/pretrim
	
	# Cabernet franc

	CF=$WD/00-RawData/00_01-Reads/00_01_01-Methylation/00_01_01_02-CabFranc
	ln -s $CF/CF_1_S1_R1_001.fastq.gz CF01_R1.fq.gz
	ln -s $CF/CF_1_S1_R2_001.fastq.gz CF01_R2.fq.gz
	ln -s $CF/CF_3_S2_R1_001.fastq.gz CF03_R1.fq.gz
	ln -s $CF/CF_3_S2_R2_001.fastq.gz CF03_R2.fq.gz
	ln -s $CF/R320-L1-P2-TGACCA-READ1-Sequences.txt.gz CF04_R1.fq.gz
	ln -s $CF/R320-L1-P2-TGACCA-READ2-Sequences.txt.gz  CF04_R2.fq.gz

	# Sauvignon blanc

	SB=$WD/00-RawData/00_01-Reads/00_01_01-Methylation/00_01_01_03-SauvBlanc
	ln -s $SB/R320-L1-P1-ATCACG-READ1-Sequences.txt.gz SB01_R1.fq.gz
	ln -s $SB/R320-L1-P1-ATCACG-READ2-Sequences.txt.gz SB01_R2.fq.gz
	ln -s $SB/SB_14_S4_R1_001.fastq.gz SB14_R1.fq.gz
	ln -s $SB/SB_14_S4_R2_001.fastq.gz SB14_R2.fq.gz
	ln -s $SB/SB_6_S3_R1_001.fastq.gz SB06_R1.fq.gz
	ln -s $SB/SB_6_S3_R2_001.fastq.gz SB06_R2.fq.gz

	# Cabernet sauvignon

	CS=$WD/00-RawData/00_01-Reads/00_01_01-Methylation/00_01_01_01-CabSauv
	ln -s $CS/CS_6_S1_L001_R1_001.fastq.gz CS06_R1.fq.gz
	ln -s $CS/CS_6_S1_L001_R2_001.fastq.gz CS06_R2.fq.gz
	ln -s $CS/R320-L1-P3-ACTTGA-READ1-Sequences.txt.gz CS08_R1.fq.gz
	ln -s $CS/R320-L1-P3-ACTTGA-READ2-Sequences.txt.gz CS08_R2.fq.gz
	ln -s $CS/CS_47_S2_L001_R1_001.fastq.gz CS47_R1.fq.gz
	ln -s $CS/CS_47_S2_L001_R2_001.fastq.gz CS47_R2.fq.gz
	
	CS=$WD/00-RawData/00_01-Reads/00_01_01-Methylation/00_01_01_04-CabSauv_BSseq_protocol
	ln -s $CS/P17-TCATTC-READ1-Sequences.txt.gz CS08cv_R1.fq.gz
	ln -s $CS/P17-TCATTC-READ2-Sequences.txt.gz CS08cv_R2.fq.gz
	
	
# Check read length

	#!/bin/bash
	WD=/data/vondras/methylome.project
	FQ=$WD/methylation/fq/pretrim
	cd $FQ	
	for i in $(ls *.fq.gz); do echo $i; zcat $i | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' ; done


# FASTQC and Trim

	cd /data/vondras/methylome.project/methylation/fq/pretrim
	ls *.gz | sed -e 's/.fq.gz//g' -e 's/_R.//g' | sort | uniq > ../../files.list

	#!/bin/bash
	WD=/data/vondras/methylome.project/methylation
	FQ=$WD/fq	
	FQC=$WD/fastqc
	for i in $(cat $1); do 
		fastqc $FQ/pretrim/${i}_R1.fq.gz $FQ/pretrim/${i}_R2.fq.gz -o $FQC/pretrim
		trim_galore --quality 20 --length 80 --fastqc --fastqc_args "-o $FQC/trimmed" --dont_gzip --paired $FQ/pretrim/${i}_R1.fq.gz $FQ/pretrim/${i}_R2.fq.gz --output_dir $FQ/trimmed
	done


# Downsample

Number of reads per library

	CS08cv_R1_val_1.fq	 41,118,142 	Do not downsample
	SB01_R1_val_1.fq	 82,462,299 	DS #
	CF03_R1_val_1.fq	 88,045,633 	
	SB14_R1_val_1.fq	 96,856,402 	
	CF04_R1_val_1.fq	 110,756,543 	
	CS08_R1_val_1.fq	 112,933,656 	
	SB06_R1_val_1.fq	 113,545,565 	
	CF01_R1_val_1.fq	 117,168,372 	
	CS06_R1_val_1.fq	 154,458,619 	
	CS47_R1_val_1.fq	 154,940,833 	
		
		
Doesn't get downsampled (CS08cv), remove from file list and zip it up

	#!/bin/bash
	WD=/data/vondras/methylome.project
	FQ=$WD/methylation/fq
        mv $FQ/trimmed/CS08cv*fq $FQ/downsampled/
        cd $FQ/downsampled
	gzip CS08cv_R1_val_1.fq 
        gzip CS08cv_R2_val_2.fq

Downsample

	#!/bin/bash
	FQ=/data/vondras/methylome.project/methylation/fq
	for i in $(cat $1); do
		seqtk sample -s100 $FQ/trimmed/${i}_R1_val_1.fq 82462299  | gzip > $FQ/downsampled/${i}_R1.fq.gz
		seqtk sample -s100 $FQ/trimmed/${i}_R2_val_2.fq 82462299 | gzip > $FQ/downsampled/${i}_R2.fq.gz
		rm $FQ/trimmed/${i}_R1_val_1.fq $FQ/trimmed/${i}_R2_val_2.fq
	done	


# Bismark

[Bismark_User_Guide.pdf](https://github.com/amvondras/methylome/files/6216666/Bismark_User_Guide.pdf)
 
*Files produced by bismark:*

1. BAM alignment file: sample__bismark_bt2_pe.bam
2. general report: sample_bismark_bt2_PE_report.txt


*Files produced by bismark_methylation_extractor:*
	
1. bedgraph: sample.bedGraph.gz
2. context-specific chg, chh, and cpg files: CHG/CHH/CpG_context_sample.txt.gz
3. cx: sample.CX_report.txt.gz
4. mbias: sample.M-bias.txt
5. general report: sample_splitting_report.txt
6. coverage: sample.bismark.cov.gz


*(I) Prepare reference (CS and lambda), making sure it will be the same one used for RNAseq and small RNA analysis*

	cd data/vondras/methylome.project/reference

	HAPONE=/data/vondras/methylome.project/00-RawData/00_02-Genomes/VITVvi_vCabSauv08_v1.1.pseudomolecules.hap1.fasta
	HAPTWO=/data/vondras/methylome.project/00-RawData/00_02-Genomes/VITVvi_vCabSauv08_v1.1.pseudomolecules.hap2.fasta
	GENOME=/DATA7/All_vinifera_genomes/Genomes/VITVvi_vCabSauv08_v1.1.pseudomolecules.all.fasta
	LAMBDA=/DATA/users/darcantu/FD_methylation/lambda/lambda.fasta
	programs=/DATA/users/darcantu/FD_methylation/programs/
		
	ln -s $HAPONE cs.h1/
	ln -s $HAPTWO cs.h2/
	ln -s $GENOME cs.h1.h2/ # For later use by Salmon, RNAseq analysis
	cp $LAMBDA lambda/

	$programs/Bismark/bismark_genome_preparation cs.h1
	$programs/Bismark/bismark_genome_preparation cs.h2
	$programs/Bismark/bismark_genome_preparation lambda




*(II) Run bismark and deduplicate the bams produced*

	#!/bin/bash

	WD=/data/vondras/methylome.project
	FQ=$WD/methylation/fq
	programs=/DATA/users/darcantu/FD_methylation/programs
	

	for i in $( ls $FQ/downsampled/*.gz | sed 's/.*.\///g' | sed 's/_.*//g' | sort | uniq  ); do
	echo ${i}

	# on lambda
	cd $WD/methylation/bismark/lambda
	ln -s $FQ/downsampled/${i}_R1.fq.gz ; ln -s $FQ/downsampled/${i}_R2.fq.gz
	REFLAM=/data/vondras/methylome.project/reference/lambda
	$programs/Bismark/bismark --bowtie2 -o $WD/methylation/bismark/lambda $REFLAM -1 ${i}_R1.fq.gz -2 ${i}_R2.fq.gz --parallel 40 -p 4 -N 1
	rm ${i}_R1.fq.gz ${i}_R2.fq.gz

	# on CS
	cd $WD/methylation/bismark/cs
	ln -s $FQ/downsampled/${i}_R1.fq.gz ; ln -s $FQ/downsampled/${i}_R2.fq.gz
		
        # on CS Hap1 
        REFCS=/data/vondras/methylome.project/reference/cs.h1
        $programs/Bismark/bismark --bowtie2 -o $WD/methylation/bismark/cs $REFCS -1 ${i}_R1.fq.gz -2 ${i}_R2.fq.gz --parallel 40 -p 4 -N 1
        mv ${i}_R1_bismark_bt2_pe.bam  ${i}.h1.bam
        mv ${i}_R1_bismark_bt2_PE_report.txt ${i}_bismark_report.h1.txt
	$programs/Bismark/deduplicate_bismark --paired --bam --output_dir bam ${i}.h1.bam
                
        # on CS Hap2
        REFCS=/data/vondras/methylome.project/reference/cs.h2
        $programs/Bismark/bismark --bowtie2 -o $WD/methylation/bismark/on_cs $REFCS -1 ${i}_R1.fq.gz -2 ${i}_R2.fq.gz --parallel 40 -p 4 -N 1
        mv ${i}_R1_bismark_bt2_pe.bam  ${i}.h2.bam
	mv ${i}_R1_bismark_bt2_PE_report.txt  ${i}_bismark_report.h2.txt
	$programs/Bismark/deduplicate_bismark --paired --bam --output_dir bam ${i}.h2.bam
               
	rm ${i}_R1.fq.gz ${i}_R2.fq.gz
	done



*Mapping efficiency:*
	
	cd /data/vondras/methylome.project/methylation/bismark/cs/bismark_report
	cat *.txt | grep 'Bismark report for\|efficiency' | paste - -

*Conversion efficiency:*

`Conversion efficiency = 1 - {no. of methylated Cs in CpG context/(no of methylated Cs in CpG context + no of unmethylated Cs in CpG context)}`
	
	#!/bin/bash

	cd /data/vondras/methylome.project/methylation/bismark/lambda
	for line in $( ls *report.txt | sed 's/_.*.//g' ); do echo ${line} ; grep "C's in CpG context" ${line}_R1_bismark_bt2_PE_report.txt | cut -f 2 | paste - - ; done	


*(III) Extract methylation*
	
	#!/bin/bash
	
	ALIGN=/data/vondras/methylome.project/methylation/bismark/cs/bam
	programs=/DATA/users/darcantu/FD_methylation/programs

	cd $ALIGN
        
        REFCS=/data/vondras/methylome.project/reference/cs.h1
	for i in $(ls *.h1.deduplicated.bam | sed 's/.deduplicated.bam//g' ); do
		$programs/Bismark/bismark_methylation_extractor --paired-end --gzip --multicore 40 --buffer_size 12G --scaffolds --CX --cytosine_report --comprehensive --bedGraph --ignore 15 --ignore_r2 15 --ignore_3prime 2 --ignore_3prime_r2 2 --genome_folder $REFCS  -o /data/vondras/methylome.project/methylation/bismark/cs ${i}.deduplicated.bam 	
	done

        REFCS=/data/vondras/methylome.project/reference/cs.h2
        for i in $(ls *.h2.deduplicated.bam | sed 's/.deduplicated.bam//g' ); do
	        $programs/Bismark/bismark_methylation_extractor --paired-end --gzip --multicore 40 --buffer_size 12G --scaffolds --CX --cytosine_report --comprehensive --bedGraph --ignore 15 --ignore_r2 15 --ignore_3prime 2 --ignore_3prime_r2 2 --genome_folder $REFCS -o /data/vondras/methylome.project/methylation/bismark/cs ${i}.deduplicated.bam 	
        done


