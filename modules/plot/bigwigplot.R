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

# make y labels smaller if too big
wrap_string <- function(str, x) {
  # Split the string into chunks of x characters
  split_str <- substring(str, seq(1, nchar(str), by = x), seq(x, nchar(str) + x - 1, by = x))
  # Combine the chunks with a newline separator
  paste(split_str, collapse = "-\n")
}


plot.bigwigs = function(chr, pos1 , pos2 , pos , plot.bw, bwcolor, bigwig.scale){

  names(plot.bw) = lapply(names(plot.bw), function(x) wrap_string(x,30))

  if(!is.na(pos) && SuREbigwig == TRUE){
    GTs = query_snps[7:12]
    GTs = gsub(0, query_snps$REF, GTs)
    GTs = gsub(1, query_snps$ALT, GTs)
    names(plot.bw) = paste0(names(plot.bw),"\n","(", GTs, ")")
  }

  empty_bw <- c()
  non_empty_bw <- list()
  for(bw in names(plot.bw)){
    if( nrow(data.frame(import.bw(plot.bw[[bw]], which = GRanges(seqnames = chr, ranges = IRanges(start = pos1, end = pos2))))) == 0 ){
      empty_bw <- c(empty_bw,bw)
    } else {
      non_empty_bw[[bw]] <- plot.bw[[bw]]
    }
  }
  
  if(length(non_empty_bw)!=0){
    non_empty_bw.plots <- BigwigTrack(
        region = GRanges(seqnames = chr, ranges = IRanges(start = pos1, end = pos2)),
        bigwig = non_empty_bw,
        smooth = 100,
        type = "coverage",
        y_label = "",
        ymax = NULL,
        max.downsample = 1000,
        downsample.rate = 0.1,
        bigwig.scale = bigwig.scale) + 
        xlab("") +
        scale_fill_manual(values = rep(bwcolor, length(non_empty_bw) )) +  
        theme_bw() +
        theme(legend.position = "none")

    if (!is.na(pos)) {
      non_empty_bw.plots <- non_empty_bw.plots + geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
    }
      } 

  empty_bw.plots <- list()
  for(bw in empty_bw){
    pl <- ggplot() +
            geom_line(aes(x = c(pos1, pos2), y = c(0, 1)), linetype = "dashed", color = "transparent") +
            labs(fill = "") +
            xlab("") +
            ylab(bw) +
            guides(color = FALSE) +
            coord_cartesian(xlim = c(pos1, pos2)) +
            theme_bw() +
            theme(
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),  # Remove y-axis ticks
              axis.line.y = element_blank(),   # Remove y-axis line
              panel.border = element_rect(color = "black")  # Keep the box around the plot
            )
    
    if (!is.na(pos)) {
      pl <- pl + geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
    }

    empty_bw.plots[[bw]] <- pl

  }

  if(length(empty_bw)>0){
    if(length(non_empty_bw)>0){
      bigwig.plot <- wrap_plots(
            c(list(non_empty_bw.plots), empty_bw.plots),
            ncol = 1, heights = c(0.2 * length(non_empty_bw), rep(0.1,length(empty_bw))))
    } else {
      bigwig.plot <- wrap_plots(
            empty_bw.plots,
            ncol = 1, heights = rep(0.1,length(empty_bw)) )
    }
  } else {
      bigwig.plot <- non_empty_bw.plots
  }

  return(bigwig.plot)
}
