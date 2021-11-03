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



	
	
