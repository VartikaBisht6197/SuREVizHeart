source(file.path(DBQueryScriptsDIR,"DBquery.gene.r"))
genedata <- suREVizGeneQuery(chr, pos1, pos2, file.path(DBDIR, "gencode.v46.annotation.hg38.genes.db"))
source(file.path(plotDIR,"geneplot.R"))
gene.plot <- plot.genes(chr, pos1, pos2, pos, genedata) 
if (!is.na(pos)) {
  gene.plot <- gene.plot + geom_vline(xintercept = pos, colour = "black", linetype = "dashed")
}

if(nrow(gene.plot$data) == 0){
  # It is an empty plot which would give errors so plot an empty plotly
  # Create an empty Plotly plot
  gene.plotly = ggplot() + 
    labs(fill = "") + xlab("") + ylab("Gene") + guides(color = FALSE) +
    coord_cartesian(xlim = c(pos1, pos2)) + theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  # Add the SNP seleted
  if(!is.na(pos)){gene.plotly =  gene.plotly + geom_vline(xintercept = pos, colour = "black", linetype = "dashed")}
  gene.plotly = ggplotly(gene.plotly)
  
} else {
  # Convert ggplot object to plotly
  gene.plotly <- ggplotly(gene.plot, tooltip = "text", showlegend = F, dynamicTicks = FALSE)
}
