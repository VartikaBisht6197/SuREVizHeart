# Get SNP.PLOT
source(file.path(plotDIR, "sureplot.R"), local = TRUE)
source(file.path(DBQueryScriptsDIR,"DBquery.main.r"))
suredata <- suREVizMainQuery(chr, pos1, pos2, file.path(DBDIR,"SuREViz.Main.db"))

if(nrow(suredata) == 0){
  # It is an empty plot which would give errors so plot an empty plotly
  # Create an empty Plotly plot
  snp.plot = ggplot() + geom_line(aes(x = c(pos1, pos1), y = c(0, 1)), linetype = "dashed") +
    labs(fill = "") + xlab("") + ylab("SuRE Signal") + guides(color = FALSE) +
    coord_cartesian(xlim = c(pos1, pos1)) + theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  # Add the SNP seleted
  if(!is.na(pos)){snp.plot =  snp.plot + geom_vline(xintercept = pos, colour = "black", linetype = "dashed")}
  snp.plot = ggplotly(snp.plot)
  
} else {
  # Convert ggplot object to plotly
  static.snp.plot <- snps_plot(chr, pos1, pos2, pos, suredata, pvalcutoff, SureSignalcutoff, query_snps)
  snp.plot <- ggplotly(static.snp.plot,tooltip = "text", showlegend = TRUE, dynamicTicks = TRUE)
}
