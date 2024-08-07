# Load required libraries
library(rtracklayer)
library(Signac)
library(RColorBrewer)
library(ggplot2)

# Create a BigwigTrack object
# It represents a track in a genomic browser displaying coverage information from a bigWig file
# Parameters:
#   region: A GRanges object specifying the genomic region to display
#   plot.bw: The path to the bigWig file containing coverage data
#   smooth: The smoothing window size for the coverage data (in base pairs)
#   extend.upstream: The distance to extend the region upstream
#   extend.downstream: The distance to extend the region downstream
#   type: The type of track, in this case, "coverage"
#   y_label: The label for the y-axis
#   ymax: The maximum y-value for the track, set to NULL for automatic scaling
#   max.downsample: The maximum number of data points to display on the track
#   downsample.rate: The downsample rate for data points
#   bigwig.scale: The scaling method for the bigWig data
plot.bigwigs = function(chr, pos1 , pos2 , pos , plot.bw, bwcolor, bigwig.scale){

  bigwig.plot = BigwigTrack(
    region = GRanges(seqnames = chr, ranges = IRanges(start = pos1, end = pos2)),
    bigwig = plot.bw,
    smooth = 100,
    type = "coverage",
    y_label = "",
    ymax = NULL,
    max.downsample = 1000,
    downsample.rate = 0.1,
    bigwig.scale = bigwig.scale
  ) + xlab("") + theme_bw() + scale_fill_manual(values = rep(bwcolor, length(plot.bw) )) +  theme_classic() +
    theme(legend.position = "none")

  if(!is.na(pos)){
    bigwig.plot = bigwig.plot + geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
  }

  return(bigwig.plot)
}
