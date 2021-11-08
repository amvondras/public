# Analysis of small RNA sequencing data from Cabernet sauvignon, Sauvignon blanc, and Cabernet franc

[sRNA analysis training](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/srna/tutorial.html#read-quality-checking)


## Index references for mapping

*Index rRNA*

- I retrieved rRNA sequences from https://www.ncbi.nlm.nih.gov/nuccore using these search criteria
- rRNA[All Fields] AND (plants[filter] AND biomol_rrna[PROP] AND genbank[filter] AND is_nuccore[filter]) 
- (("Vitis vinifera"[Organism] OR Vitis vinifera[All Fields]) AND ribosomal rna[All Fields]) AND "Vitis vinifera"[porgn] AND (plants[filter] AND genbank[filter] AND is_nuccore[filter])
- I downloaded these fasta files:[rrna sequences.zip](https://github.com/amvondras/methylome/files/6269038/rrna.sequences.zip)
- I edited manually to remove `patterns`

	grep ">" rrna.fasta | grep -i -v 'rrna\|ribo' | grep -v "5S" > patterns # remove the two sequences not named rRNAs manually


        #!/bin/bash
        INDEX=/data/vondras/methylome.project/reference
        cd $INDEX/rRNA
        bowtie2-build rrna.fasta rrna

*Index genome references*

        #!/bin/bash
        REFERENCE=/data/vondras/methylome.project/reference
	cd $REFERENCE
        bowtie2-build $REFERENCE/cs.h1/VITVvi_vCabSauv08_v1.1.pseudomolecules.hap1.fasta VITVvi_vCabSauv08_v1.1.pseudomolecules.hap1
	bowtie2-build $REFERENCE/cs.h2/VITVvi_vCabSauv08_v1.1.pseudomolecules.hap2.fasta VITVvi_vCabSauv08_v1.1.pseudomolecules.hap2

## Get data

	cd /data/vondras/methylome.project/smallrna/fq/original/berry
        for i in $( ls /data/vondras/methylome.project/00-RawData/00_01-Reads/00_01_02-SmallRNA/*.gz ); do ln -s "$i" "`paste <(echo ${i} | sed -e 's/.*\///g' -e 's/_.*.fastq.gz//g' -e 's/-/_/g') <(echo "_berry.fq.gz") | sed 's/\t//g' `"; done # Berry
	
	cd /data/vondras/methylome.project/smallrna/fq/original/leaf
        for i in $( ls /DATA12/seq/Trio_RNAseq_smallRNA/small_RNA/*R1*.fastq.gz ); do ln -s "$i" "`paste <(echo ${i} | sed -e 's/.*\///g' -e 's/s_S.._L.*._R.//g' -e 's/_001.fastq.gz//g' -e 's/-//g') <(echo "_leaf.fq.gz") | sed 's/\t//g'`"; done # Leaf, R1 only

	cd /data/vondras/methylome.project/smallrna/fq/original/new_leaf
	for i in $( ls /DATA/seq/CS_CF_CB_smallRNA/clones/*.gz | grep -v "Undet" ); do ln -s "$i" "`paste <(echo ${i} | sed -e 's/^.*.\///g' -e 's/S_S.*//g') <(echo "_nleaf.fq.gz") | sed 's/\t//g' `"; done # New Leaf libraries

	cd /data/vondras/methylome.project/smallrna/fq/original/
	mv */* .
	rm -rf berry leaf new_leaf

## Retrieve FASTQC already run on original data

	FQC=/data/vondras/methylome.project/smallrna/fastqc/original
	
	cd $FQC/berry
	for i in $( ls /data/vondras/methylome.project/00-RawData/00_01-Reads/00_01_02-SmallRNA/fastqc/*.html ); do ln -s "$i" "`paste <(echo ${i} | sed -e 's/.*\///g' -e 's/_.*//g' | sed -e 's/-/_/g') <(echo "_berry.fastqc.html") | sed 's/\t//g' `"; done # Berry

	cd $FQC/leaf 
	for i in $( ls /DATA12/seq/Trio_RNAseq_smallRNA/small_RNA/*R1*_fastqc.html ); do ln -s "$i" "`paste <(echo ${i} | sed -e 's/.*\///g' -e 's/s_S.._L.*._R.//g' -e 's/_001_fastqc.html//g' -e 's/-//g') <(echo "_leaf.fastqc.html") | sed 's/\t//g'`"; done # Leaf, R1 only

	cd $FQC/new_leaf
	for i in $( ls /DATA/seq/CS_CF_CB_smallRNA/clones/*.html | grep -v "Undet" | grep -v "lane" ); do ln -s "$i" "`paste <(echo ${i} | sed -e 's/^.*.\///g' -e 's/S_S.*//g') <(echo "_nleaf.fastqc.html") | sed 's/\t//g' `"; done # New Leaf libraries



## Trim and FASTQC; Only use R1 leaf sequencing data.
        
        FQ=/data/vondras/methylome.project/smallrna/fq
        cd $FQ/original
	ls *.fq.gz | sed 's/.fq.gz//g' | sort | uniq > file.list

	#!/bin/bash
        FQ=/data/vondras/methylome.project/smallrna/fq
        FQC=/data/vondras/methylome.project/smallrna/fastqc

	for i in $(cat $1); do       	
	        trim_galore --phred33 --small_rna --gzip --quality 0 --length 18 --max_length 27 --fastqc --fastqc_args "--outdir $FQC/trimmed" --output_dir $FQ/trimmed $FQ/original/${i}.fq.gz 
        done

Run and save outputs

	bash launcher.sh file.list trim.sh 9
	rm *.out ; cat file*err > trim.out ; rm file.list.* file.list.*	

## Remove rRNA and align to CS

	#!/bin/bash
	REFERENCE=/data/vondras/methylome.project/reference
	ALIGN=/data/vondras/methylome.project/smallrna/align
	FQ=/data/vondras/methylome.project/smallrna/fq
        FQC=/data/vondras/methylome.project/smallrna/fastqc

 	for i in $( ls $FQ/trimmed/*.gz | sed 's/^.*.\///g' | sed 's/_trimmed.fq.gz//g' | sort ); do 
	
	echo -e '$i\t#AMV'
	echo "#AMV"
	
	echo "rRNA alignment"
	# Remove rRNA
	bowtie2 -x $REFERENCE/rRNA/rrna -U $FQ/trimmed/${i}_trimmed.fq.gz --un-gz $FQ/rrna_filter/${i}.fq.gz -S $FQ/rrna_filter/${i}.sam
	#samtools view -bS $FQ/rrna_filter/${i}.sam | samtools sort -o $FQ/rrna_filter/${i}.bam
	rm $FQ/rrna_filter/${i}.sam
	fastqc $FQ/rrna_filter/${i}.fq.gz --outdir $FQC/rrna_filter

	echo "hap1 alignment"
	# Align to CS
	bowtie2 -x $REFERENCE/cs.h1/VITVvi_vCabSauv08_v1.1.pseudomolecules.hap1 -U $FQ/rrna_filter/${i}.fq.gz -S $ALIGN/${i}.h1.sam --al-gz $FQ/genome_aligned/${i}.h1.fq.gz 
	
	echo "hap2 alignment"
	bowtie2 -x $REFERENCE/cs.h2/VITVvi_vCabSauv08_v1.1.pseudomolecules.hap2 -U $FQ/rrna_filter/${i}.fq.gz -S $ALIGN/${i}.h2.sam --al-gz $FQ/genome_aligned/${i}.h2.fq.gz

	# SAM to sorted BAM
	cd $ALIGN
	samtools view -bS $ALIGN/${i}.h1.sam | samtools sort -o $ALIGN/${i}.h1.bam
	samtools view -bS $ALIGN/${i}.h2.sam | samtools sort -o $ALIGN/${i}.h2.bam
	rm $ALIGN/${i}.h1.sam $ALIGN/${i}.h2.sam 

	done

## Check how many reads are left after the rRNA filter

	#!/bin/bash
	cd /data/vondras/methylome.project/smallrna/fq/rrna_filter
	for i in $(ls *.gz | sed 's/.fq.gz//g' ); do paste <(echo $i) <(echo $(zcat ${i}.fq.gz |wc -l)/4|bc); done


## Check alignment rate

	cd /data/vondras/methylome.project/smallrna

	cat nohup.out | grep '_\|overall'| grep -v 'Approx\|Started\|complete\|bam_sort_core' | sed 's/overall alignment rate//g' | tr '\n' '\t' | sed -e 's/CS/\nCS/g' -e 's/CF/\nCF/g' -e 's/SB/\nSB/g' | sort -k2,2n


## Build counts tables of unique small RNA sequences

	cd /data/vondras/methylome.project/smallrna/counts
	
	#!/bin/bash
	for i in $( ls ../fq/genome_aligned/*.fq.gz | sed -e 's/..\/fq\/genome_aligned\///g' -e 's/.fq.gz//g' | sort | uniq ); do
		echo ${i}
		seqtk seq -a ../fq/genome_aligned/${i}.fq.gz | grep -v ">" | sort | uniq -c | awk -v var=$i '{print var "\t" $1 "\t" $2}' > ${i}.count
	done
	cat *.count > all_long.count
	rm *h1.count *h2.count


	cd /data/vondras/methylome.project/smallrna/counts; R
	library(reshape2)
	strSort <- function(x)
        sapply(lapply(strsplit(x, NULL), sort), paste, collapse="")
	
	dat <- read.delim("all_long.count", header=FALSE);dat <- dat[,c(3,1,2)]
	
	names(dat) <- c("sequence", "sample", "count")
	dat$hap <- gsub("^.*\\.", "", dat$sample)
	dat$sam <- gsub("\\.h.", "", dat$sample); dat$sam <- gsub("_.*._", "_", dat$sam)
	dat$complexity <- nchar(gsub('([[:alpha:]])\\1+', '\\1', strSort(dat$sequence)))
	dat <- dat[dat$complexity > 1 & dat$complexity < 5,]

	seqs <- list()
	for (set in unique(dat$sam)){
		seqs[[set]] <- unique(dat[dat$sam==set,"sequence"])
	}; rm(set)
	save(dat,seqs, file="seq_count.RData")
	
## How many unique sequences are present

	cd /data/vondras/methylome.project/smallrna/counts; R
	load(file="seq_count.RData")
	uniqueseqs <- data.frame();for(i in names(seqs)){uniqueseqs <- rbind(uniqueseqs, data.frame(sample=i, count=length(seqs[[i]])) )}
	uniqueseqs$group <- gsub("^.*._", "", uniqueseqs$sample)
	uniqueseqs$cultivar <- gsub("_.*.", "", uniqueseqs$sample);
		uniqueseqs$cultivar <- gsub("CS.*.", "CS", uniqueseqs$cultivar);
		uniqueseqs$cultivar <- gsub("CF.*.", "CF", uniqueseqs$cultivar);
		uniqueseqs$cultivar <- gsub("SB.*.", "SB", uniqueseqs$cultivar);
	uniqueseqs <- uniqueseqs[order(uniqueseqs$group, uniqueseqs$cultivar),]
	uniqueseqs


	# Update excel sheet with unique seq info per library
	uniqueseqs <- data.frame();
	for(i in unique(dat$sample)){ uniqueseqs <- rbind(uniqueseqs, data.frame(sample=i,count=length(unique(dat[dat$sample==i, "sequence"] )) ))}
	uniqueseqs <- uniqueseqs[order(uniqueseqs$count),]


## Make simple figure(s)summarizing pipeline steps for leaf libraries

	cd /data/vondras/methylome.project/smallrna; R
	library(readxl); library(ggplot2); library(reshape2)
	
	dat <- read_excel("leaf_RNAseq_SmallRNAseq.libraries.info.xlsx", sheet = "Read bases counts")
	dat$`Original or rerun` <- gsub("Original_reanalyzed_w_reruns", "Original", dat$`Original or rerun`)
	names(dat)
	
	tmp <- dat[,c(1,3,8:10,17:18,21:22)]
	names(tmp) <- c("Run", "Sample", "Raw", "Trimmed and sized", "rRNA-filtered", "Hap1-aln", "Hap2-aln", "Unique H1 Seqs", "Unique H2 Seqs")
	tmp <- melt(tmp, id.vars = names(tmp)[1:2] )
	
	pdf("pdf/srna_pipeline_steps.pdf", width=5, height=9)	
	ggplot(tmp, aes(Sample, value)) + geom_bar(stat="identity", position = position_dodge(.9)) + 
			theme_classic() +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			facet_grid(variable~Run, scales="free", labeller = labeller(groupwrap = label_wrap_gen(10))) + ylab("Number of reads (row 1:5) or unique sequences (row 6:7)")
	dev.off()
	
	pdf("pdf/srna_pipeline_steps_side_by_side.pdf", width=9, height=5)	
	ggplot(tmp, aes(variable, value, fill=Sample)) + geom_bar(stat="identity", position = position_dodge(.9)) + 
			theme_classic() +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			facet_grid(.~Run, scales="free", labeller = labeller(groupwrap = label_wrap_gen(10))) + ylab("Number of reads (row 1:5) or unique sequences (row 6:7)")
	dev.off()
	
	tmp <- dat[,c(1,3,14,15)]
	names(tmp) <- c("Run", "Sample", "Hap1", "Hap2")
	tmp <- melt(tmp, id.vars = names(tmp)[1:2] )
	
	pdf("pdf/srna_pipeline_alignment.pdf", width=5, height=5)
	ggplot(tmp, aes(Sample, value, fill=Run)) + geom_bar(stat="identity", position = position_dodge(.9)) + 
			theme_classic() +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			facet_grid(variable~.) + ylab("Percent reads aligned")
	dev.off()

[srna_pipeline_alignment.pdf](https://github.com/amvondras/trio/files/7341510/srna_pipeline_alignment.pdf)

[srna_pipeline_steps_side_by_side.pdf](https://github.com/amvondras/trio/files/7341511/srna_pipeline_steps_side_by_side.pdf)

[srna_pipeline_steps.pdf](https://github.com/amvondras/trio/files/7341512/srna_pipeline_steps.pdf)


![small rna update 001](https://user-images.githubusercontent.com/33852065/137215863-b1fd4f85-e5dd-42c1-83f3-d5114829431b.jpeg)
![small rna update 002](https://user-images.githubusercontent.com/33852065/137215871-b00cebb8-cebb-4412-b498-f35548a2bf2f.jpeg)
![small rna update 003](https://user-images.githubusercontent.com/33852065/137215888-13a7af8c-b81a-4229-a8d8-65d4c281eeab.jpeg)

## Get intersections, leaf

	cd /data/vondras/methylome.project/smallrna/counts; R
	load(file="seq_count.RData")
	library(UpSetR)	
	names(seqs)
 	
	# Get core smallRNAs to each cultivar, based on OG data, rerun, and both runs together
	core_leaf <- list(
	CF = list(
		original = 		Reduce(intersect, list(seqs[["CF1_leaf"]], seqs[["CF3_leaf"]], seqs[["CF4_leaf"]])),
		rerun = 		Reduce(intersect, list(seqs[["CF1_nleaf"]], seqs[["CF3_nleaf"]],seqs[["CF4_nleaf"]])),	
		runs_combined = 	Reduce(intersect,list( unique(c(seqs[["CF1_leaf"]],seqs[["CF1_nleaf"]])), 
								unique(c(seqs[["CF3_leaf"]],seqs[["CF3_nleaf"]])), 
								unique(c(seqs[["CF4_leaf"]],seqs[["CF4_nleaf"]]))  )), 
		runs_reproducible= 	Reduce(intersect, list(seqs[["CF1_leaf"]], seqs[["CF3_leaf"]],seqs[["CF4_leaf"]], 
								seqs[["CF1_nleaf"]], seqs[["CF3_nleaf"]],seqs[["CF4_nleaf"]]))),
	CS = list(
		original = 		Reduce(intersect, list(seqs[["CS47_leaf"]], seqs[["CS8_leaf"]], seqs[["CS6_leaf"]])),	
		rerun = 		Reduce(intersect, list(seqs[["CS47_nleaf"]], seqs[["CS8_nleaf"]], seqs[["CS6_nleaf"]])),		
		runs_combined = 	Reduce(intersect,list( unique(c(seqs[["CS47_leaf"]],seqs[["CS47_nleaf"]])), 
								unique(c(seqs[["CS8_leaf"]],seqs[["CS8_nleaf"]])), 
								unique(c(seqs[["CS6_leaf"]],seqs[["CS6_nleaf"]]))  )), 
		runs_reproducible= 	Reduce(intersect, list(seqs[["CS47_leaf"]], seqs[["CS8_leaf"]], seqs[["CS6_leaf"]], 
								seqs[["CS47_nleaf"]], seqs[["CS8_nleaf"]], seqs[["CS6_nleaf"]]))),
	SB = list(
		original = 		Reduce(intersect, list(seqs[["SB1_leaf"]], seqs[["SB6_leaf"]], seqs[["SB14_leaf"]])),	
		rerun = 		Reduce(intersect, list(seqs[["SB1_nleaf"]], seqs[["SB6_nleaf"]],seqs[["SB14_nleaf"]])),
		runs_combined = 	Reduce(intersect,list( unique(c(seqs[["SB1_leaf"]],seqs[["SB1_nleaf"]])), 
								unique(c(seqs[["SB6_leaf"]],seqs[["SB6_nleaf"]])), 
								unique(c(seqs[["SB14_leaf"]],seqs[["SB14_nleaf"]]))  )), 
		runs_reproducible= 	Reduce(intersect, list(seqs[["SB1_leaf"]], seqs[["SB6_leaf"]],seqs[["SB14_leaf"]], 
								seqs[["SB1_nleaf"]], seqs[["SB6_nleaf"]],seqs[["SB14_nleaf"]]))))
	
	for(set in names(core)){
		for(i in names(core[[set]]) ){
			print(paste( set, i, length(core[[set]][[i]]) ))
		}	
	}; rm(set, i)

## Plot small RNA sequences in leaves and berries per cultivar	

	pdf("unique_small_rna_sequences_berry.pdf")
	berry <- seqs[grepl("berry", names(seqs))]
	mymax <- max(length(berry[["CF_berry"]]), length(core_leaf[["CS_berry"]]), length(core_leaf[["SB_berry"]])  )
	print(upset( fromList(berry), sets=c("CF_berry","CS_berry","SB_berry"), 
		#order.by = c("degree","freq"), 
		keep.order=TRUE, set_size.show=TRUE, 
		#text.scale= c(0.9, 2, 0, 0, 2, 2), 
		set_size.scale_max= mymax + 0.60 * mymax, 
		mainbar.y.label="Number shared berry small RNA sequences", 
		intersections = list(list("CS_berry", "CF_berry", "SB_berry"), list("CS_berry", "CF_berry"),list("CS_berry", "SB_berry"),list("CF_berry", "SB_berry"),list("CS_berry"), list("CF_berry"), list("SB_berry") )
	))
	dev.off(); rm(berry, mymax)

	pdf("unique_small_rna_sequences_leaf.pdf", width=8, height=5)
	for (set in names(core_leaf[[1]])){
		tmp <- list(CF=core_leaf[["CF"]][[set]], CS=core_leaf[["CS"]][[set]], SB=core_leaf[["SB"]][[set]])
		mymax <- max(length(core_leaf[["CF"]][[set]]), length(core_leaf[["CS"]][[set]]), length(core_leaf[["SB"]][[set]])  )
		print(upset(
			fromList(tmp), #order.by = c("degree","freq"), 
			keep.order=TRUE, set_size.show=TRUE, 
			#text.scale= c(0.9, 2, 0, 0, 2, 2), 
			set_size.scale_max= mymax + 0.50 * mymax, 
			sets.x.label=set, 
			mainbar.y.label="N shared leaf small RNA sequences", 
			intersections = list(list("CS", "CF", "SB"), list("CS", "CF"),list("CS", "SB"),list("CF", "SB"),list("CS"), list("CF"), list("SB") )
		)); rm(tmp, mymax)
	}; 
	dev.off(); rm(set)




[unique_small_rna_sequences_berry.pdf](https://github.com/amvondras/trio/files/7356104/unique_small_rna_sequences_berry.pdf)

[unique_small_rna_sequences_leaf.pdf](https://github.com/amvondras/trio/files/7356105/unique_small_rna_sequences_leaf.pdf)

