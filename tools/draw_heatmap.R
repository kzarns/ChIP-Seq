#! /export/apps/r-3.0.2/bin/Rscript
 

usage <- "USAGE: Rscript draw_heatmap.R <in_file_name> <out_file_name>"
args <- commandArgs(TRUE)
if("--help" %in% args){
	print(usage)
	q()
}
infile <- args[1]
outfile <- args[2]
print(sprintf("infile: %s, outfile: %s", infile, outfile))

require("RColorBrewer")
wr <- colorRampPalette(c("white", "pink", "red", "darkred", "black", bias=1))(n = 299)
png(outfile)
data <-read.delim(infile, header=TRUE, sep="\t", skip=8)
data_matrix <- data.matrix(data[,-1:-6])
data_heatmap <- heatmap(data_matrix, Rowv=NA, Colv=NA, col=wr, margins=c(5,6))

q(status=1)
