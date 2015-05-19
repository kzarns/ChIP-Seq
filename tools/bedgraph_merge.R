#!/usr/bin/Rscript
##################################################################
## This script will normalize colections of  bedgraph files in a 
## folder based on the number of mapped reads present in those 
## bedgraph files.  files are normalized down to the lowest number 
## of reads to minimilize the noise in downstream analysis.
##################################################################

require(plyr) || stop("The plyr library is not available!")

##################################################################
## Functions

write.table.with.header <- function(x, file, header, ...){
    cat(header, '\n',  file = paste(gsub("^(.*)[.].*", "\\1", file), 
	'norm_merge.bedgraph', sep=""))
    write.table(x, file = paste(gsub("^(.*)[.].*", "\\1", file), 
	'norm_merge.bedgraph', sep=""), append = T, ...)
}

frame.merge <- function(x, y){
        ddply(merge(x, y, all=TRUE), .(chr, start, end), 
                summarise, count=sum(count))
}

##################################################################

file_list <- list.files(pattern = '*.bedgraph')
frame_list <- lapply(file_list, read.delim, skip=1, header=FALSE, col.names =
			c('chr', 'start', 'end', 'count'))

## Sum the counts column of all frames and save a vector of those counts.
mapped_reads <- lapply(frame_list, function(x) sum(x$count))

## Find the smallest and create a list of normalization factors
norm_factors <- lapply(mapped_reads, function(x) x/min(unlist(mapped_reads)))

## Divide the lists by their nf
norm_frame_list <- Map(function(x, i) {x[, 4]=x[, 4]/i; return(x)}, 
			frame_list, norm_factors)

## merge all data frames and sum counts in overlapping rows
merge_frame <- Reduce(frame.merge, norm_frame_list)

## Create a generic header
header <- "track type=bedGraph visibility=full color=179,27,27 altColor=179,27,27 priority=20"

## Write back to beadgraph file
Map(write.table.with.header, merge_frame, file_list, header=header, sep="\t", 
	row.names=FALSE, col.names=FALSE, quote=FALSE)

