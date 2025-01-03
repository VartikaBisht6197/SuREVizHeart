#-----------------------------------------------------------------------------------------------
# Author: Vartika Bisht
# Date: 19.08.2024
#
# Description:
# This script, named geneplot.R, is designed to generate a gene plot for visualizing genomic
# features around a specified genomic position. It utilizes the 'ggplot2' and 'plotly' libraries
# for static and interactive visualizations, respectively. The script assumes the existence of
# a dataset named 'genedata' that contains genomic information and annotations.
#-----------------------------------------------------------------------------------------------

# Function definition
plot.genes <- function(chr, pos1, pos2, pos, genedata) {
  # Convert the 'strand' column to a factor with meaningful levels for proper legends
  genedata$strand <- factor(genedata$strand)
  levels(genedata$strand) <- list("Plus Strand Gene" = "+", "Minus Strand Gene" = "-")

  # Determine the strands present in the data
  strand_levels <- levels(genedata$strand)

  # Dynamically set color mapping based on the strands present
  if (all(strand_levels %in% unique(genedata$strand))) {
    # Both strands are present
    color_values <- c("Plus Strand Gene" = "blue", "Minus Strand Gene" = "red")
  } else if ("Plus Strand Gene" %in% unique(genedata$strand)) {
    # Only Plus Strand Gene is present
    color_values <- c("Plus Strand Gene" = "blue")
  } else {
    # Only Minus Strand Gene is present
    color_values <- c("Minus Strand Gene" = "red")
  }

  # Create a ggplot object for gene visualization
  gene.plot <- ggplot(genedata, aes(text = paste0(
    "Gene Name: ", GENE, "\n",
    "Strand: ", strand, "\n",
    "Type: ", type, "\n"
  ))) +
    geom_segment(aes(x = start, xend = end, y = GENE, yend = GENE, color = strand)) +

    # Customize color scales for strand
    scale_color_manual(values = color_values) +
    scale_fill_manual(values = color_values) +

    # Customize labels and appearance
    labs(fill = "") +
    xlab(chr) +
    ylab("Genes") +
    guides(color = FALSE) +
    theme_bw() +
    coord_cartesian(xlim = c(pos1, pos2))

  # Add extra rectangles to show exons if available
  if (sum(genedata$type %in% c("exon")) != 0) {
    gene.plot <- gene.plot + geom_rect(
      data = genedata[genedata$type %in% c("exon"), ],
      aes(
        xmin = start, xmax = end, ymin = GENE, ymax = GENE,
        color = strand, fill = strand
      ),
      size = 5
    )
  }

  return(gene.plot)
}


# Source the script to query the gene database
source(file.path(DBQueryScriptsDIR, "DBquery.gene.r"))

# Define the path to the gene SQLite database
dbSQL <- file.path(dataDIR, "gencode.v46.annotation.hg38.genes.db")

# Query gene data based on the specified chromosome and position range
genedata <- suREVizGeneQuery(chr, pos1, pos2, dbSQL)

# Generate the gene plot
gene.plot <- plot.genes(chr, pos1, pos2, pos, genedata)

# Highlight the specific genomic position if provided
if (!is.na(pos)) {
  gene.plot <- gene.plot + geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
}

# Create an empty Plotly plot if no data is found
if (nrow(gene.plot$data) == 0) {
  gene.plot <- ggplot() +
    geom_line(
      aes(x = c(pos1, pos2), y = c(0, 1)),
      linetype = "dashed", color = "transparent"
    ) +
    labs(fill = "") +
    xlab("") +
    ylab("Genes") +
    guides(color = FALSE) +
    coord_cartesian(xlim = c(pos1, pos2)) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) 

  if (!is.na(pos)) {
    gene.plot <- gene.plot + geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
  }
  gene.plotly <- ggplotly(gene.plot)
} else {
  # Convert ggplot object to interactive Plotly plot
  gene.plotly <- ggplotly(gene.plot, tooltip = "text", showlegend = FALSE, dynamicTicks = FALSE)
}
