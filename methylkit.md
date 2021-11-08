# Do differential methylation analysis

Bismark analyzed BSseq reads on each haplotype separately. MethylKit can compare samples along a single haplotype in one context at a time, but not across haplotypes. Where are cultivars differentially methylated along a single haplotype? Then, use synteny information to determine whether cultivars are differentially methylated over the same region on the opposite haplotype.


With pipeline='bismarkCytosineReport', the function expects cytosine report files from Bismark, which have chr,start, strand, number of cytosines (methylated bases) , number of thymines (unmethylated bases), context and trinucletide context format. These are the CX files.


	VITVvi_vCabSauv08_v1.1.hap1.chr01	2	+	0	0	CHG	CAG

## methylKit doesn't like all contexts combined, so split CX files

	cd /data/vondras/methylome.project/methylation/methylkit

        #!/bin/bash
        CX=/data/vondras/methylome.project/methylation/bismark/cs/cx
        for i in $( ls $CX/*.gz | sed -e 's/^.*.\///g' -e 's/.deduplicated.CX_report.txt.gz//g' ); do
		echo ${i}
	        zcat $CX/${i}.deduplicated.CX_report.txt.gz | awk '$6=="CHG"'| gzip > chg/chg_${i}.gz 
	        zcat $CX/${i}.deduplicated.CX_report.txt.gz | awk '$6=="CG"' | gzip > cg/cg_${i}.gz 
	        zcat $CX/${i}.deduplicated.CX_report.txt.gz | awk '$6=="CHH"'| gzip > chh/chh_${i}.gz 
        done

## MethylKit

_Import, format, and save data_
  
  	cd /data/vondras/methylome.project/methylation/methylkit; R
        setwd("/data/vondras/methylome.project/methylation/methylkit")
        library(methylKit)

  	# First haplotype
	file.list.cg=list(
               "cg/cg_CF01.h1.gz", "cg/cg_CF03.h1.gz", "cg/cg_CF04.h1.gz",
               "cg/cg_CS06.h1.gz", "cg/cg_CS08.h1.gz", "cg/cg_CS47.h1.gz",
               "cg/cg_SB01.h1.gz", "cg/cg_SB06.h1.gz", "cg/cg_SB14.h1.gz")   
	       
	file.list.chg=list(
               "chg/chg_CF01.h1.gz", "chg/chg_CF03.h1.gz", "chg/chg_CF04.h1.gz",
               "chg/chg_CS06.h1.gz", "chg/chg_CS08.h1.gz", "chg/chg_CS47.h1.gz",
               "chg/chg_SB01.h1.gz", "chg/chg_SB06.h1.gz", "chg/chg_SB14.h1.gz") 	   
	
	file.list.chh=list(
		"chh/chh_CF01.h1.gz", "chh/chh_CF03.h1.gz", "chh/chh_CF04.h1.gz",
 		"chh/chh_CS06.h1.gz", "chh/chh_CS08.h1.gz", "chh/chh_CS47.h1.gz",
		"chh/chh_SB01.h1.gz", "chh/chh_SB06.h1.gz", "chh/chh_SB14.h1.gz") 	   

	# Second haplotype	
	file.list.cg.h2=list(
               "cg/cg_CF01.h2.gz", "cg/cg_CF03.h2.gz", "cg/cg_CF04.h2.gz",
               "cg/cg_CS06.h2.gz", "cg/cg_CS08.h2.gz", "cg/cg_CS47.h2.gz",
               "cg/cg_SB01.h2.gz", "cg/cg_SB06.h2.gz", "cg/cg_SB14.h2.gz")   
	       
	file.list.chg.h2=list(
               "chg/chg_CF01.h2.gz", "chg/chg_CF03.h2.gz", "chg/chg_CF04.h2.gz",
               "chg/chg_CS06.h2.gz", "chg/chg_CS08.h2.gz", "chg/chg_CS47.h2.gz",
               "chg/chg_SB01.h2.gz", "chg/chg_SB06.h2.gz", "chg/chg_SB14.h2.gz") 	   
	
	file.list.chh.h2=list(
		"chh/chh_CF01.h2.gz", "chh/chh_CF03.h2.gz", "chh/chh_CF04.h2.gz",
 		"chh/chh_CS06.h2.gz", "chh/chh_CS08.h2.gz", "chh/chh_CS47.h2.gz",
		"chh/chh_SB01.h2.gz", "chh/chh_SB06.h2.gz", "chh/chh_SB14.h2.gz")
		
	dat <- list(
 		cg_h1 = methRead(file.list.cg, 
                        sample.id=list( 
				"CF01","CF03","CF04",
                        	"CS06","CS08","CS47",
                        	"SB01","SB06","SB14"), 
                        treatment=c(0,0,0,1,1,1,2,2,2),mincov = 0,
                        assembly="CS", 
                        context=c("CpG"), 
                        pipeline="bismarkCytosineReport"),

        	chg_h1 = methRead(file.list.chg, 
                          sample.id=list( 
				"CF01","CF03","CF04",
                        	"CS06","CS08","CS47",
                        	"SB01","SB06","SB14"),   
                        treatment=c(0,0,0,1,1,1,2,2,2),mincov = 0,
                        assembly="CS", 
                        context=c("CHG"), 
                        pipeline="bismarkCytosineReport"),

        	chh_h1 = methRead(file.list.chh, 
                         sample.id=list( 
				"CF01","CF03","CF04",
                        	"CS06","CS08","CS47",
                        	"SB01","SB06","SB14"),    
                        treatment=c(0,0,0,1,1,1,2,2,2),mincov = 0,
                        assembly="CS", 
                        context=c("CHH"), 
                        pipeline="bismarkCytosineReport"),
			
 		cg_h2 = methRead(file.list.cg.h2, 
                        sample.id=list( 
				"CF01","CF03","CF04",
                        	"CS06","CS08","CS47",
                        	"SB01","SB06","SB14"), 
                        treatment=c(0,0,0,1,1,1,2,2,2),mincov = 0,
                        assembly="CS", 
                        context=c("CpG"), 
                        pipeline="bismarkCytosineReport"),

        	chg_h2 = methRead(file.list.chg.h2, 
                          sample.id=list( 
				"CF01","CF03","CF04",
                        	"CS06","CS08","CS47",
                        	"SB01","SB06","SB14"),  
                        treatment=c(0,0,0,1,1,1,2,2,2),mincov = 0,
                        assembly="CS", 
                        context=c("CHG"), 
                        pipeline="bismarkCytosineReport"),

        	chh_h2 = methRead(file.list.chh.h2, 
                         sample.id=list( 
				"CF01","CF03","CF04",
                        	"CS06","CS08","CS47",
                        	"SB01","SB06","SB14"),    
                        treatment=c(0,0,0,1,1,1,2,2,2),mincov = 0,
                        assembly="CS", 
                        context=c("CHH"), 
                        pipeline="bismarkCytosineReport"))

	save(dat, file="methylkit.RData")
	
	
_Check coverage and % methylation distribution_

	cd /data/vondras/methylome.project/methylation/methylkit; R
	setwd("/data/vondras/methylome.project/methylation/methylkit")
	library(methylKit)
	load("methylkit.RData")
	
	for(set in names(dat)){
		print(set)
		for(mysam in 1:9){
			getCoverageStats(dat[[set]][[mysam]],plot=FALSE,both.strands=FALSE)
		}	
	}; rm(set)

	
	
	
	pdf("getMethylationStats.pdf")
	for(set in names(dat)){
		print(set)
		print(getMethylationStats(dat[[set]][[2]],plot=TRUE,both.strands=FALSE))
		print(getCoverageStats(dat[[set]][[2]],plot=TRUE,both.strands=FALSE))
	}; rm(set)
	dev.off()
	
_Filter, unite, and save_

	dat_f <- list(); dat_u <- list()
	for(set in names(dat)){
		print(set)
		dat_f[[set]] <-filterByCoverage(dat[[set]], lo.count=10)
		dat_u[[set]] <- unite(dat_f[[set]], destrand=TRUE)
       	}; rm(set) 
	
	dat_u_split <- list()
	for(set in names(dat_u)){
		print(set)
		dat_u_split[[set]] = list(
			all = dat_u[[set]], 
			CS_CF = reorganize(dat_u[[set]],sample.ids=c("CS06","CS08","CS47","CF01","CF03","CF04"),treatment=c(0,0,0,1,1,1) ),
			CS_SB = reorganize(dat_u[[set]],sample.ids=c("CS06","CS08","CS47","SB01","SB06","SB14"),treatment=c(0,0,0,1,1,1) ),
			CF_SB = reorganize(dat_u[[set]],sample.ids=c("CF01","CF03","CF04","SB01","SB06","SB14"),treatment=c(0,0,0,1,1,1) )
		)
	}; rm(set, dat_f, dat_u)
	dat <- dat_u_split; rm(dat_u_split)
	save(dat, file="methylkit_filtered_united.RData")


_PCAs and clustering_

	cd /data/vondras/methylome.project/methylation/methylkit; R
	setwd("/data/vondras/methylome.project/methylation/methylkit")
	library(methylKit)
	load("methylkit_filtered_united.RData")

	pdf("pdf/pcas_and_clustering_haps_separate.pdf")
	for(set in names(dat)){
		print(set)
		getCorrelation(dat[[set]][[1]],plot=TRUE)
        	PCASamples(dat[[set]][[1]])
        	clusterSamples(dat[[set]][[1]], dist="correlation", method="ward", plot=TRUE)
	}	; rm(set)
	dev.off()


_Differential methylation along each haplotype at base resolution_

	# https://groups.google.com/g/methylkit_discussion/c/mTpy8I9dJUk?pli=1

	cd /data/vondras/methylome.project/methylation/methylkit; R
	setwd("/data/vondras/methylome.project/methylation/methylkit")
	library(methylKit)
	load("methylkit_filtered_united.RData")

	meth_diff <- list()
	for(set in names(dat)){
		print(set)
		for(p in names(dat[[1]])){
			print(p)
			meth_diff[[set]][[p]] <- calculateDiffMeth(dat[[set]][[p]], overdispersion="MN", mc.cores=24 )
			a <- nrow(meth_diff[[set]][[p]] ) ; b <- nrow(meth_diff[[set]][[p]][ meth_diff[[set]][[p]]$pvalue < 0.01 , ] )
			print(a) ; print(b) ; print(round(b / a * 100, 2)); rm(a, b)
		}
	}; rm(set, p)

	save(meth_diff, file="methylkit_diff.RData")



_Calculate methylation in windows / "tiles" and do differential analysis_

	cd /data/vondras/methylome.project/methylation/methylkit; R
	setwd("/data/vondras/methylome.project/methylation/methylkit")
	library(methylKit)
	load("methylkit.RData")
	
	tiles <- list(); tiles_u <- list(); tiles_diff <- list(); tiles_diff_od <- list()
	for(set in names(dat)){
		print(set)
		tiles[[set]] <- tileMethylCounts(dat[[set]], win.size=1000, step.size=1000, cov.bases = 10)
		tiles_u[[set]] <- unite(tiles[[set]], destrand=TRUE)
		
		tiles_split[[set]] = list(
			all = tiles_u[[set]], 
			CS_CF = reorganize(tiles_u[[set]],sample.ids=c("CS06","CS08","CS47","CF01","CF03","CF04"),treatment=c(0,0,0,1,1,1) ),
			CS_SB = reorganize(tiles_u[[set]],sample.ids=c("CS06","CS08","CS47","SB01","SB06","SB14"),treatment=c(0,0,0,1,1,1) ),
			CF_SB = reorganize(tiles_u[[set]],sample.ids=c("CF01","CF03","CF04","SB01","SB06","SB14"),treatment=c(0,0,0,1,1,1) ))
		
		for(p in names(tiles_split[[1]])){
			tiles_diff[[set]][[p]] <- calculateDiffMeth(tiles_split[[set]][[p]], overdispersion="MN", mc.cores=24 )
			a <- nrow(tiles_diff[[set]][[p]] ) ; b <- nrow(tiles_diff[[set]][[p]][ tiles_diff[[set]][[p]]$pvalue < 0.01 , ] )
			print(a) ; print(b) ; print(round(b / a * 100, 2)); rm(a, b)
		}		
	}; rm(set, tiles_u)
	save(tiles_split, tiles_diff, file="methylkit_tiles.RData")

_Annotate differentially methylated bases / regions with [synteny](synteny/synteny_between_haplotypes_genes_mcscanx.md) and genomic feature information_

https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#4_Annotating_differentially_methylated_bases_or_regions
	
https://bioconductor.org/packages/release/bioc/vignettes/genomation/inst/doc/GenomationManual.html#3_Data_input	
	
	cd /data/vondras/methylome.project/methylation/methylkit
	ln -s ../../synteny/mcscanx/formatted_mcscanx.txt .
	
	cd /data/vondras/methylome.project/methylation/methylkit; R; 
	library(methylKit); 
	library(genomation);
	library(rtracklayer)
	
	load("methylkit_tiles.RData") # meth_diff
	load("methylkit_diff.RData") ; rm(tiles_split) # tiles_diff, tiles_split

	setwd("/data/vondras/methylome.project/methylation/methylkit/grangebeds")
	
	syn <- readGeneric("../formatted_mcscanx.txt", header = TRUE, chr = 2, start = 3, end = 4, keep.all.metadata = TRUE)
	write.table(syn, file="h1_formatted_mcscanx.bed", sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

	syn <- readGeneric("../formatted_mcscanx.txt", header = TRUE, chr = 6, start = 7, end = 8, keep.all.metadata = TRUE)
	write.table(syn, file="h2_formatted_mcscanx.bed", sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

	methylkit_tiles_GRanges <- list()	
	for(set in names(tiles_diff)){
		print(set)
		for(p in names(tiles_diff[[set]])){
			print(p)
			methylkit_tiles_GRanges[[set]][[p]] <- as(tiles_diff[[set]][[p]], "GRanges")
			write.table(methylkit_tiles_GRanges[[set]][[p]], file=paste("tiles_",set, "_", p, ".bed", sep=""), sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)	
		}
	}; rm(set, p)	
		
On command line:

	#!/bin/bash	
	for i in $(ls tiles*h1*.bed | sed 's/.bed//g' ); do
		echo $i
		bedtools intersect -a ${i}.bed -b h1_formatted_mcscanx.bed -wa -wb > ${i}_ann.bed
	done
	for i in $(ls tiles*h2*.bed | sed 's/.bed//g' ); do
        	echo $i
        	bedtools intersect -a ${i}.bed -b h2_formatted_mcscanx.bed -wa -wb > ${i}_ann.bed
	done	

Relate hap1 and hap2 in R using the ann.bed files and identify syntenic genes differentially methylated on both haplotypes.

	cd /data/vondras/methylome.project/methylation/methylkit/grangebeds; R
	list.files()
	
	dat = list(
                cg_all 		= list(read.delim("tiles_cg_h1_all_ann.bed", header=FALSE),  read.delim("tiles_cg_h2_all_ann.bed", header=FALSE)),
		cg_cf_sb 	= list(read.delim("tiles_cg_h1_CF_SB_ann.bed", header=FALSE),read.delim("tiles_cg_h2_CF_SB_ann.bed", header=FALSE)),
		cg_cs_cf 	= list(read.delim("tiles_cg_h1_CS_CF_ann.bed", header=FALSE),read.delim("tiles_cg_h2_CS_CF_ann.bed", header=FALSE)),
		cg_cs_sb 	= list(read.delim("tiles_cg_h1_CS_SB_ann.bed", header=FALSE),read.delim("tiles_cg_h2_CS_SB_ann.bed", header=FALSE)),
		chg_all 	= list(read.delim("tiles_chg_h1_all_ann.bed", header=FALSE),  read.delim("tiles_chg_h2_all_ann.bed", header=FALSE)),
		chg_cf_sb 	= list(read.delim("tiles_chg_h1_CF_SB_ann.bed", header=FALSE),read.delim("tiles_chg_h2_CF_SB_ann.bed", header=FALSE)),
		chg_cs_cf 	= list(read.delim("tiles_chg_h1_CS_CF_ann.bed", header=FALSE),read.delim("tiles_chg_h2_CS_CF_ann.bed", header=FALSE)),
		chg_cs_sb 	= list(read.delim("tiles_chg_h1_CS_SB_ann.bed", header=FALSE),read.delim("tiles_chg_h2_CS_SB_ann.bed", header=FALSE)),
		chh_all 	= list(read.delim("tiles_chh_h1_all_ann.bed", header=FALSE),  read.delim("tiles_chh_h2_all_ann.bed", header=FALSE)),
		chh_cf_sb 	= list(read.delim("tiles_chh_h1_CF_SB_ann.bed", header=FALSE),read.delim("tiles_chh_h2_CF_SB_ann.bed", header=FALSE)),
		chh_cs_cf 	= list(read.delim("tiles_chh_h1_CS_CF_ann.bed", header=FALSE),read.delim("tiles_chh_h2_CS_CF_ann.bed", header=FALSE)),
		chh_cs_sb 	= list(read.delim("tiles_chh_h1_CS_SB_ann.bed", header=FALSE),read.delim("tiles_chh_h2_CS_SB_ann.bed", header=FALSE)))
	for(set in names(dat)){ names(dat[[set]]) <- c("hap1", "hap2") }; rm(set)
	
	for(set in names(dat)){	
		names(dat[[set]][[1]])<- c(
		"h1_data_chr", "h1_data_start", "h1_data_end", "h1_data_width", 
		"h1_data_strand", "h1_data_pvalue", "h1_data_qvalue", "h1_data_meth.diff", 	
		"h1_chr", "h1_start", "h1_stop", "h1_width", "h1_strand", "h1_gene", 
		"h2_gene","h2_chr", "h2_start", "h2_stop", 
		"Block_ID", "Position_in_Block", "Block_score", "Block_e.value",
		"Block_loci_count", "Block_Direction", "table" )
		
		names(dat[[set]][[2]])<-c(
		"h2_data_chr", "h2_data_start", "h2_data_end", "h2_data_width", 
		"h2_data_strand", "h2_data_pvalue", "h2_data_qvalue", "h2_data_meth.diff",		
		"h2_chr", "h2_start", "h2_stop", "h2_width", "h2_strand", 
		"h1_gene","h1_chr", "h1_start", "h1_stop",
		"h2_gene", 
		"Block_ID", "Position_in_Block", "Block_score", "Block_e.value",
		"Block_loci_count", "Block_Direction", "table" )
	}; rm(set)
	
	for(set in names(dat)){	
		dat[[set]][[1]] <- dat[[set]][[1]][ ,c( "h1_data_chr", "h1_data_start", "h1_data_end", "h1_data_pvalue", "h1_data_meth.diff", 
	"h1_chr", "h1_start", "h1_stop", "h1_gene",           
	"h2_chr", "h2_start", "h2_stop", "h2_gene" )]
	
		dat[[set]][[2]] <- dat[[set]][[2]][ ,c( "h2_data_chr", "h2_data_start", "h2_data_end", "h2_data_pvalue", "h2_data_meth.diff", 
	 "h1_chr", "h1_start", "h1_stop", "h1_gene",
	 "h2_chr", "h2_start", "h2_stop", "h2_gene" )]
	
	}; rm(set)
	
	for(set in names(dat)){	
		dat[[set]][["merged"]] <- unique(merge(dat[[set]][[1]], dat[[set]][[2]], by=c("h1_chr", "h1_start", "h1_stop","h1_gene","h2_chr", "h2_start", "h2_stop", "h2_gene"), all=TRUE))		
		dat[[set]][["merged_ic"]] <- unique(dat[[set]][["merged"]][!complete.cases(dat[[set]][["merged"]]), ]) 
		dat[[set]][["merged_cc"]] <- unique(dat[[set]][["merged"]][complete.cases(dat[[set]][["merged"]]), ])
		dat[[set]][["merged_cc"]]$dm_count <- rowSums(dat[[set]][["merged_cc"]][,c("h1_data_pvalue", "h2_data_pvalue") ] < 0.05) 
	}; rm(set)
	
	res <- dat; rm(dat)
	save(res, file="/data/vondras/methylome.project/methylation/methylkit/methylkit_diff_synteny.RData")


Append gene functional annotation to results

	cd /data/vondras/methylome.project/reference/functional
	ln -s /DATA7/All_vinifera_genomes/Functional_annotations/VITVvi_vCabSauv08_v1.1_functional_annotation.txt .

	# format gene annotation for general use
	cat VITVvi_vCabSauv08_v1.1_functional_annotation.txt | sed 's/ /__/g' | cut -f 1,3 | sed 's/ /\t/g' | sed 's/__/ /g' | sed 's/\.p..\t/\t/g' | sort | uniq > CS_genefunctionalannotation.txt
	
	# format transcript annotation for general use
	cat VITVvi_vCabSauv08_v1.1_functional_annotation.txt | sed 's/ /__/g' | cut -f 1,3 | sed 's/ /\t/g' | sed 's/__/ /g' | sort | uniq > CS_transcriptfunctionalannotation.txt
	
	# format GO annotation for general use
	# gene : go 
	cat VITVvi_vCabSauv08_v1.1_functional_annotation.txt | sed 's/ /__/g' | cut -f 1,4 | sed 's/__/ /g'| sed '1d' | awk 'BEGIN{OFS=" "}{sign=""; {for (i=2;i<=NF; i++){print $1 "\t" $i}}}' | sed -e 's/,//g' -e 's/"//g' -e 's/\.p..\t/\t/g' | awk '{print $2 "\t" $1}' | sort | uniq | sed '1igo\tgene'  > gene.go.key

	# go : go-id
	 paste <( cat VITVvi_vCabSauv08_v1.1_functional_annotation.txt  | cut -f 4,5 | awk '$1!="-" && $2 !="-"' |  cut -f 1 | sed 's/ /__/g' | sed -e 's/,__/\n/g' -e 's/;__/\n/g'  -e 's/"//g' ) <( cat VITVvi_vCabSauv08_v1.1_functional_annotation.txt  | cut -f 4,5 | awk '$1!="-" && $2 !="-"'  | cut -f 2 |  sed -e 's/"//g' -e 's/, F\:/\nF\:/g' -e 's/, P\:/\nP\:/g' -e 's/, C\:/\nC\:/g'  -e 's/; F\:/\nF\:/g' -e 's/; P\:/\nP\:/g' -e 's/; C\:/\nC\:/g'   | less -S) | sed '1d' | sort | uniq |sed '1igo\tname' > go.name.key

	R
	functional <- list(gene_annotation = read.delim("/data/vondras/methylome.project/reference/functional/CS_genefunctionalannotation.txt"), 
		transcript_annotation = read.delim("/data/vondras/methylome.project/reference/functional/CS_transcriptfunctionalannotation.txt"),
		go_gene = read.delim("/data/vondras/methylome.project/reference/functional/gene.go.key"),
		go_name = read.delim("/data/vondras/methylome.project/reference/functional/go.name.key"))
	save(functional, file="/data/vondras/methylome.project/reference/functional/annotation.RData")


Summarize differentially methylated pairs of syntenic genes

	cd /data/vondras/methylome.project/methylation/methylkit/; R
	load("methylkit_diff_synteny.RData")
	load("/data/vondras/methylome.project/reference/functional/annotation.RData")

	library(clusterProfiler)

	for(set in names(res)){
		print(set)
		tmp <- unique(merge(res[[set]][["merged_cc"]] , 
			data.frame( "h1_gene" = functional$gene_annotation$SeqName, "h1_fun" = functional$gene_annotation$Description ),
			all.x=TRUE, all.y=FALSE))	
		tmp <- unique(merge(tmp, 
			data.frame( "h2_gene" = functional$gene_annotation$SeqName, "h2_fun" = functional$gene_annotation$Description ),
			all.x=TRUE, all.y=FALSE))	
		res[[set]][["merged_cc_annotated"]] <- tmp ; rm(tmp)
		tmp <- enricher(unique(c( 
			res[[set]][["merged_cc"]][res[[set]][["merged_cc"]]$dm_count > 0, "h1_gene"] , 
			res[[set]][["merged_cc"]][res[[set]][["merged_cc"]]$dm_count > 0, "h2_gene"] )), 
			pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = functional$go_gene, TERM2NAME = functional$go_name)	 
		res[[set]][["enrichres"]]  <- tmp@result ; rm(tmp)
	}; rm(set)

	for(set in names(res)){print(set); print(res[[set]][["enrichres"]][res[[set]][["enrichres"]]$p.adjust < 0.1 ,"Description"]); print("#########")};rm(set)
	
	save(res, file="/data/vondras/methylome.project/methylation/methylkit/methylkit_diff_synteny.RData")


Visualize differentially methylated pairs of syntenic genes

	cd /data/vondras/methylome.project/methylation/methylkit/; R
	load("methylkit_diff_synteny.RData")
	library(GenomicFeatures); 
	
	names(res[[1]])
	#[1] "hap1"                "hap2"                "merged"             
	#[4] "merged_ic"           "merged_cc"           "merged_cc_annotated"
	#[7] "enrichres"
 
 	names(res[[1]][["merged_cc"]])
 	#[1] "h1_chr"            "h1_start"          "h1_stop"          
 	#[4] "h1_gene"           "h2_chr"            "h2_start"         
 	#[7] "h2_stop"           "h2_gene"           "h1_data_chr"      
	#[10] "h1_data_start"     "h1_data_end"       "h1_data_pvalue"   
	#[13] "h1_data_meth.diff" "h2_data_chr"       "h2_data_start"    
	#[16] "h2_data_end"       "h2_data_pvalue"    "h2_data_meth.diff"
	#[19] "dm_count"  
