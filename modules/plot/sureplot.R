#-----------------------------------------------------------------------------------------------
# Author: Vartika Bisht
# Date: 19.08.2024
#
# Description:
# This script, named sureplot.R, is designed to generate a SuRE (Survey of Regulatory Elements)
# signal plot for a specified genomic variant. It utilizes the 'ggplot2' library for creating
# static plots and 'plotly' for interactive visualizations. The script assumes the presence of
# a SQLite database named 'SuREViz.Main.db' and a function `suREVizMainQuery` for querying the
# dataset. The `snps_plot` function visualizes SuRE signals and highlights specific SNPs based
# on user inputs and significance thresholds.
#
#-----------------------------------------------------------------------------------------------

# Function to create SuRE signal plot
snps_plot <- function(chr, pos1, pos2, pos, suredata, pvalcutoff, SureSignalcutoff, query_snps) {
  # Modify the suredata dataframe within the function based on significance thresholds
  suredata$raQTL <- ifelse(max(suredata$ref.mean, suredata$alt.mean) > SureSignalcutoff & suredata$wilcox.p.value <= pvalcutoff, TRUE, FALSE)
  suredata$max.allele <- ifelse(suredata$ref.mean > suredata$alt.mean,
    "Loss of function\n( REF > ALT )",
    "Gain of function\n( ALT > REF )"
  )
  suredata$gnomad312_AF <- round(as.numeric(suredata$gnomad312_AF) * 100, 2)
  suredata$sig <- ifelse(suredata$wilcox.p.value < pvalcutoff,
    -log10(suredata$wilcox.p.value) * 100,
    -log10(suredata$wilcox.p.value) / 100
  )
  suredata$min.allele.exp <- signif(mapply(min, suredata$ref.mean, suredata$alt.mean), 3)
  suredata$max.allele.exp <- signif(mapply(max, suredata$ref.mean, suredata$alt.mean), 3)

  # Create ggplot object for SuRE signal
  p_signal <- ggplot(suredata, aes(
    x = POS,
    y = ref.mean,
    xend = POS,
    yend = alt.mean,
    alpha = sig,
    color = max.allele,
    size = -log10(wilcox.p.value),
    text = paste0(
      "Click to centre on this variant", "\n", "\n",
      "Position: ", paste0(CHROM, "-", POS, "-", REF, "-", ALT), "\n",
      "p-value: ", format(wilcox.p.value, scientific = TRUE, digits = 2), "\n",
      "Alternate allele frequency: ", gnomad312_AF, "%\n",
      "Reference allele: ", REF, "\n",
      "Alternate allele: ", ALT, "\n",
      max.allele, "\n",
      "Min allele expression: ", min.allele.exp, "\n",
      "Max allele expression: ", max.allele.exp
    )
  )) +
    geom_segment() +
    guides(size = "none", alpha = "none") +
    scale_color_manual(
      values = c("#3c5379", "#64BEA5"),
      breaks = c("Gain of function\n( ALT > REF )", "Loss of function\n( REF > ALT )"),
      guide = guide_legend(title = "")
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.title.x = element_blank()
    ) +
    ylab("SuRE Signal") +
    xlim(c(pos1, pos2)) +
    coord_cartesian(xlim = c(pos1, pos2))

  # Highlight specific SNPs if available
  if (!is.na(pos) && nrow(query_snps) != 0) {
    p_signal <- p_signal +
      geom_segment(aes(x = pos, xend = pos, y = query_snps$ref.mean, yend = query_snps$alt.mean),
        color = "yellow", size = 0.5, alpha = 1
      ) +
      geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
  }

  return(p_signal)
}

# Source the script to query the main SQL database
source(file.path(DBQueryScriptsDIR, "DBquery.main.r"))

# Define the path to the SQLite database
dbSQL <- file.path(dataDIR, "SuREViz.Main.db")

# Query SuRE data based on the specified chromosome and position range
suredata <- suREVizMainQuery(chr, pos1, pos2, dbSQL)

if (nrow(suredata) == 0) {
  # Create an empty Plotly plot if no data is found
  snp.plot <- ggplot() +
    geom_line(aes(x = c(pos1, pos1), y = c(0, 1)), linetype = "dashed") +
    labs(fill = "") +
    xlab("") +
    ylab("SuRE Signal") +
    guides(color = FALSE) +
    coord_cartesian(xlim = c(pos1, pos1)) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  # Add the SNP position if specified
  if (!is.na(pos)) {
    snp.plot <- snp.plot + geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
  }
  snp.plot <- ggplotly(snp.plot)
} else {
  # Generate and convert ggplot object to Plotly interactive plot
  static.snp.plot <- snps_plot(chr, pos1, pos2, pos, suredata, pvalcutoff, SureSignalcutoff, query_snps)
  snp.plot <- ggplotly(static.snp.plot, tooltip = "text", showlegend = TRUE, dynamicTicks = TRUE)
}
