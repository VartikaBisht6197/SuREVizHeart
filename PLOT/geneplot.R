################################################################################
# Gene Plotting Script
# --------------------
# This R script generates a gene plot for visualization of genomic features around
# a specified genomic position. It uses the 'ggplot2' and 'plotly' libraries for
# visualization and assumes the existence of a 'genedata' dataset with relevant
# genomic information.
#
# Usage:
# source("geneplot.R"); plot.genes(chr, pos, flank, genedata)
#
# Parameters:
#   chr      : Chromosome
#   pos      : Genomic position
#   flank    : Flanking region size
#   genedata : Gene dataset (assumed to be loaded externally)
#
# Example:
# source("geneplot.R"); plot.genes("chr19", 56252947, 100000, genedata)
################################################################################

# Required Libraries: ggplot2, GenomicRanges
require(ggplot2)
require(GenomicRanges)

# Function definition
plot.genes = function(chr, pos1, pos2, pos, genedata) {
  # Convert the 'strand' column to a factor with meaningful levels for proper legends
  genedata$strand = factor(genedata$strand)
  levels(genedata$strand) = list("Plus Strand Gene" = "+", "Minus Strand Gene" = "-")
  
  # Create a ggplot for gene visualization
  gene.plot = ggplot(genedata, aes( text = paste0("Gene Name: ", GENE , "\n",
                                                  "Strand: ", strand, "\n",
                                                  "Type: ", type , "\n") )) +
    geom_segment(data = genedata, mapping = aes(x = start, xend = end, y = GENE, yend = GENE, color = strand)) +
  
    # Customize color scales for strand
    # For the gene segment, you have color and fill but for rectangles you have fill
    scale_color_manual(values = c("blue", "red"),labels = c("Plus Strand Gene","Minus Strand Gene")) +
    scale_fill_manual(values = c( "blue", "red"),labels = c("Plus Strand Gene","Minus Strand Gene")) +
    
    # Customize labels and appearance
    # Gene plot is on the top, we do not want x label.
    # No fill or color legend but show the guides for only fill
    labs(fill = "") +
    xlab(chr) +
    ylab("Genes") +
    guides(color = FALSE) +
    theme_bw() +
    # Set the x-axis limits to focus on the specified region
    # coord_cartesian does this without altering the data which is plotted
    coord_cartesian(xlim = c(pos1, pos2))
    
  # Add extra rectangles to show exons
  # This can be changed based on what characteristic you want to highlight
  if (sum(genedata$type %in% c("exon")) != 0){
    gene.plot = gene.plot + geom_rect(data = genedata[genedata$type %in% c("exon"), ],
              aes(xmin = start, xmax = end, ymin = GENE, ymax = GENE, color = strand, fill = strand),
              size = 5)}  
  
  return(gene.plot)
}
