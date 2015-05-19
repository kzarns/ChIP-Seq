#!/usr/bin/Rscript
##################################################################
## This script uses a cpp_result file to collect the reannotation
## results and reformat them to be used with make_heatmap
##################################################################

require(reshapes2) || stop("The reshape2 library is not available!")

##################################################################
## Functions

strand_cleaner <- function(i){
        i_clean_noname <- i[, -c(4,5, 7:9)]
        i_clean <- rename(i_clean_noname, c('db.desc1' = 'X.hg19.knownGene.name', 
                'db.desc2' = 'hg19.kgXref.geneSymbol', 'db.chr' = 
                'hg19.knownGene.chrom'))
}

#################################################################

## read in files and original gene table
file_list <- list.files(pattern='*_deduplicated.cpp_result')
first_frame_list <- lapply(file_list, read.delim2)
hg_original <- read.delim2('hg19.gene_list')

## Trim frames, seperate Description vector, and correct format
first_frame_list <- lapply(first_frame_list, 
				function(x) {x[, -c(6:8, 10, 11)]})
first_frame_list <- lapply(first_frame_list, function(x) 
		{cbind(x, colsplit(x$db.desc2, ':', c('Description', 'start')))})
first_frame_list <- lapply(first_frame_list, function(x) {x$start <-
				 as.integer(x$start); return(x)})

## Add a column containing the difference between start and db start
first_frame_list <- lapply(first_frame_list, function(x) 
			{x$diff <- x$start - x$db.start; return(x)})
## Seperate plus and minus strand genes
plus_strand_genes <- first_frame_list[[1]]
plus_strand_genes <- plus_strand_genes[!plus_strand_genes$diff == 200, ]
minus_strand_genes <- first_frame_list[[2]]
minus_strand_genes <- minus_strand_genes[!minus_strand_genes$diff == 100, ]

## Match appropriate columns to filter plus and minus genes
hg_plus <- hg_original[hg_original$hg19.kgXref.geneSymbol 
			%in% plus_strand_genes$db.desc2, ]
hg_minus <- hg_original[hg_original$hg19.kgXref.geneSymbol 
			%in% minus_strand_genes$db.desc2, ]

plus_strand_genes <- strand_cleaner(plus_strand_genes)
minus_strand_genes <- strand_cleaner(minus_strand_genes)

## Merge old and new annotation tables, then remove old annotations
## and add a properly formatted strand designation column
plus_merge <- merge(plus_strand_genes, hg_plus)
minus_merge <- merge(minus_strand_genes, hg_minus)
plus_merge <- plus_merge[, -5]
minus_merge <- minus_merge[, -5]
plus_merge <- rename(plus_merge, c('q.start' = 'hg19.knownGenetxStart'))
minus_merge <- rename(minus_merge, c('q.start' = 'hg19.knownGenetxStart'))
plus_merge$X <- rep('plus', nrow(plus_merge))
minus_merge$X <- rep('minus', nrow(minus_merge))
new_tss <- rbind(plus_merge, minus_merge)

write.table(new_tss, file = 'new_hg19.gene_list', sep = "\t",
                 row.names = FALSE, quote = FALSE)
