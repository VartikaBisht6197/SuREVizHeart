# This script provides functions to visualize genomic regions from BED files
# using ggplot2 and patchwork packages in R.

# Load necessary libraries
library(patchwork)
library(GenomicRanges)
library(data.table)

# Define a function to plot a single BED file
plot.single.bed <- function(chr, pos1, pos2, pos, bed_file, file_name) {
    # Read and process the BED file
    bed_data <- fread(bed_file)[, 1:3]
    colnames(bed_data) <- c("chr", "start", "end")
    bed_data <- makeGRangesFromDataFrame(bed_data)

    # Define the region of interest
    region_of_interest <- GRanges(seqnames = chr, ranges = IRanges(start = pos1, end = pos2))

    # Subset the bed_data to include only the specified region
    bed_data <- as.data.frame(bed_data[queryHits(findOverlaps(bed_data, region_of_interest)), ])

    # Plot the BED file
    if(nrow(bed_data)>0){
    bed.plot <- ggplot(bed_data, aes(x = start, xend = end, y = 1, yend = 1)) +
        geom_segment(size = 2) +
        theme_bw() +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_text(vjust = 1, hjust = 0.5, size = 14)
        ) +
        labs(y = file_name, x = "") +
        coord_cartesian(xlim = c(pos1, pos2)) } else {
           bed.plot = ggplot() +
            geom_line(aes(x = c(pos1, pos2), y = c(0, 1)), linetype = "dashed") +
            labs(fill = "") +
            xlab("") +
            ylab("") +
            guides(color = FALSE) +
            coord_cartesian(xlim = c(pos1, pos2)) +
            theme_bw() +
            theme(
              axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
              axis.text.y = element_blank(), axis.ticks.y = element_blank())
        }

    # Add vertical line if pos specified
    if (!is.na(pos)) {
        bed.plot <- bed.plot + geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
    }

    return(bed.plot)
}

# Define a function to plot multiple BED files
plot.beds <- function(chr, pos1, pos2, pos, bed_files) {
    beds <- list()
    for (i in 1:length(bed_files)) {
        beds <- c(beds, list(plot.single.bed(chr, pos1, pos2, pos, as.character(bed_files)[i], names(bed_files)[i])))
    }

    return(wrap_plots(beds, ncol = 1))
}
