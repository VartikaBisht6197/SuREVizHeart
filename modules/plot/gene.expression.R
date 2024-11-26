
extra_genes = c("TBX18", "GATA4", "HNF4A")
protien.coding.genes <- fread(file.path(dataDIR, "Homo.Sapiens.Gene.List.txt"), header = FALSE)$V1
get.expression.plot = function(genes){
  genes = genes[genes %in% protien.coding.genes]
  if(length(genes)>0){
    genes = c(genes,extra_genes)
    normcount = read.table(file.path(dataDIR, "Tissue_DevelopmentStages_CellLine_normcount.txt"))
    # Subset to genes of interest
    cnts = normcount[genes,]
    g = pheatmap(t(cnts), cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = FALSE, silent=TRUE)
  } else {
    g = ggplot() + 
      theme_minimal() +
      annotate("text", x = 0.5, y = 0.5, label = "No protien coding gene within the defined window.", size = 6, hjust = 0.5, vjust = 0.5) +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
      )
  }
  return(g)
}


