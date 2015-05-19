#!/usr/bin/Rscript
##################################################################
## This script will generate a heatmap and metagene plot from a 
## .heatmap matrix (make_heatmap) 
##################################################################

require(gplots) || stop("The gplots library is not available!")
require(RColorBrewer) || stop("The RcolorBrewer library is not available!")
require(reshape2) || stop("The reshape2 library is not available!")

##################################################################
## Functions

heatmap_gen_list <- function(x, i) {
    pdf(paste(gsub("^(.*)[_].*", "\\1", i), '_orig_no_orph_heatmap.pdf', sep = ""),
        width = 4, height = 50)
    heatmap.2(data.matrix(x[, 4:(ncol(x)-2)]), margins = c(4, 1.75),
              dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none', col = palette,
              breaks = c_break, labRow = '', labCol = '', main = paste(gsub("^(.*)[_].*", "\\1", i)),
              lmat = lmat, lwid = lwid, lhei = lhei, key.title = NA,
              key.ylab = NA, key.xlab = 'Aligned Reads')
    dev.off()
    plot.new()
}

metagene_gen_list <- function(x, i) {
    pdf(paste(gsub("^(.*)[.].*", "\\1", i), '_metagene.pdf', sep = ""))
    plot(colSums(Filter(is.numeric, x[, -(ncol(x)-1)])), 
         type = 'l', col = 'darkgreen', ann = FALSE)
    title(main = paste(gsub("^(.*)[.].*", "\\1", i)), xlab = 'Distance from TSS', 
          ylab = 'Total Reads Across all Genes')
    dev.off()
    plot.new()
}

##################################################################

file_list <- list.files(pattern='*.heatmap')
frame_list <- lapply(file_list, read.delim2, skip = 8, dec = '.')

## optional list of genes to order the heatmap. This file should a text  
## file listing one gene per line
order_file <- 'orig_gene_list_no_orph.txt'

## Drop unused columns from the data frames
drops <- c('Gene.Start', 'Gene.End', 'Strand')
frame_list <- lapply(frame_list, function(x) {x[, !(names(x) %in% drops)]})

## Sum all bins for each gene and attach as a column 
vector_list <- lapply(frame_list, function(x) {apply(x[, -(1:3)], 1, sum)})
frame_list <- Map(function(x, i) {cbind(x, i)}, frame_list, vector_list)

## Split description column and add gene names back to the frame 
frame_list <- lapply(frame_list, function(x) {cbind(x, colsplit(x$Description, 
			':', c('Description', 'Start'))$Description)})
frame_list <- lapply(frame_list, function(x) {colnames(x)[ncol(x)] <- 'Names'; return(x)})

## This code will slice out and order the data frames using a predifined 
## set of genes.  Will also save the genes removed from the analysis
order_genes <- read.table(order_file)
order_genes <- as.character(order_genes$x)
order_orphans_list <- lapply(frame_list, function(x) {x[!x$Names %in% order_genes, ]})
frame_list <- lapply(frame_list, function(x) {x[x$Names %in% order_genes, ]})
final_list <- lapply(frame_list, function(x) {x[match(order_genes, x$Names), ]})

## Check for complete cases to avoid errors in plotting
final_list <- lapply(frame_list, function(x) {x[complete.cases(x), ];return(x)})

#Write final data frame to csv
Map(function(x, i) {write.csv(x, file = paste(gsub("^(.*)[_].*", "\\1", i), 
			'_cleaned.csv', sep = ""))}, final_list, file_list)
 
## Write orphans from indexing as well
Map(function(x, i) {write.csv(x, file = paste(gsub("^(.*)[_].*", "\\1", i), 
			'orphans.csv', sep = ""))}, unique_orphans_list, file_list)

## This code sets some of the color options for heatmap creation
c_break <- c(seq(0, 7, length = 5), seq(7, 200, length = 10), seq(200, 250, length = 5))
palette <- colorRampPalette(c("white", "darkgreen", "black"))(n = 19)
lmat = rbind(c(0,3),c(2,1),c(0,4))
lwid = c(0.25, 4)
lhei = c(0.1, 4, 0.1)

## Heatmap generation.  Best used with 10nt bins
Map(heatmap_gen_list, final_list, file_list)

## Metagene generation.  Best used with 1nt bins 
Map(metagene_gen_list, final_list, file_list)


