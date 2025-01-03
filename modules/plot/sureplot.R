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
  suredata$AF <- round(as.numeric(suredata$AF) * 100, 2)
  suredata$sig <- ifelse(suredata$wilcox.p.value < pvalcutoff,
    -log10(suredata$wilcox.p.value) * 100,
    -log10(suredata$wilcox.p.value) / 100
  )
  suredata$min.allele.exp <- signif(mapply(min, suredata$ref.mean, suredata$alt.mean), 3)
  suredata$max.allele.exp <- signif(mapply(max, suredata$ref.mean, suredata$alt.mean), 3)

  # Modify the query_snps dataframe within the function based on significance thresholds if variant view 
  if (!is.na(pos) && nrow(query_snps) != 0) {
  query_snps$raQTL <- ifelse(max(query_snps$ref.mean, query_snps$alt.mean) > SureSignalcutoff & query_snps$wilcox.p.value <= pvalcutoff, TRUE, FALSE)
  query_snps$max.allele <- ifelse(query_snps$ref.mean > query_snps$alt.mean,
    "Loss of function\n( REF > ALT )",
    "Gain of function\n( ALT > REF )"
  )
  query_snps$AF <- round(as.numeric(query_snps$AF) * 100, 2)
  query_snps$sig <- ifelse(query_snps$wilcox.p.value < pvalcutoff,
    -log10(query_snps$wilcox.p.value) * 100,
    -log10(query_snps$wilcox.p.value) / 100
  )
  query_snps$min.allele.exp <- signif(mapply(min, query_snps$ref.mean, query_snps$alt.mean), 3)
  query_snps$max.allele.exp <- signif(mapply(max, query_snps$ref.mean, query_snps$alt.mean), 3) 
  }

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
      "Alternate allele frequency: ", AF , "%\n",
      "Reference allele: ", REF, "\n",
      "Alternate allele: ", ALT, "\n",
      max.allele, "\n",
      "Max allele expression: ", max.allele.exp, "\n",
      "Min allele expression: ", min.allele.exp
    )
  )) +
  geom_segment() +
  guides(size = "none", alpha = "none") +
  scale_color_manual(
    values = c("#3c5379", "#64BEA5"),
    breaks = c(
      "Gain of function\n( ALT > REF )",
      "Loss of function\n( REF > ALT )"
    ),
    guide = guide_legend(title = "")
  ) +
  scale_alpha_continuous(
    limits = c(alpha_min, alpha_max)
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
      geom_segment(aes(x = query_snps$POS, xend = query_snps$POS, y = query_snps$ref.mean, yend = query_snps$alt.mean,
      text = paste0(
      "Click to centre on this variant", "\n", "\n",
      "Position: ", paste0(query_snps$CHROM, "-", query_snps$POS, "-", query_snps$REF, "-", query_snps$ALT), "\n",
      "p-value: ", format(query_snps$wilcox.p.value, scientific = TRUE, digits = 2), "\n",
      "Alternate allele frequency: ", query_snps$AF, "%\n",
      "Reference allele: ", query_snps$REF, "\n",
      "Alternate allele: ", query_snps$ALT, "\n",
      query_snps$max.allele, "\n",
      "Max allele expression: ", query_snps$max.allele.exp, "\n",
      "Min allele expression: ", query_snps$min.allele.exp
    )),
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
  static.snp.plot <- ggplot() +
    geom_line(
      aes(x = c(pos1, pos2), y = c(0, 1)),
      linetype = "dashed", color = "transparent"
    ) +
    labs(fill = "") +
    xlab("") +
    ylab("SuRE Signal") +
    guides(color = FALSE) +
    coord_cartesian(xlim = c(pos1, pos2)) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    annotate(
      "text",
      x = mean(c(pos1, pos2)), y = 0.5,
      label = "No data available to plot",
      size = 5, hjust = 0.5, vjust = 0.5
    )
  # Add the SNP position if specified
  if (!is.na(pos)) {
    static.snp.plot <- static.snp.plot + geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
  }
  snp.plot <- ggplotly(static.snp.plot)
} else {
  # Generate and convert ggplot object to Plotly interactive plot
  static.snp.plot <- snps_plot(chr, pos1, pos2, pos, suredata, pvalcutoff, SureSignalcutoff, query_snps)
  snp.plot <- ggplotly(static.snp.plot, tooltip = "text",dynamicTicks = TRUE)
}
