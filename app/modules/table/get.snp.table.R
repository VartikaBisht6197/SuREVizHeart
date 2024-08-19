
highlight = NA
if( nrow(suredata) != 0 ){
  SNP.info = suredata
  if(!is.na(pos)){
    # highlight the pos
    highlight = which(SNP.info$POS == pos)
  }
  # colnames(SNP.info) <- c("Chromosome", "Position in human hg38",
  # "Reference allele", "Alternate allele",
  # "rsID", "Population allele frequency gnomAD 3.1.2",
  #  "Genotype in SuRE57", "Genotype in SuRE59",
  # "Genotype in SuRE67", "Genotype in SuRE68", "Genotype in SuRE86",
  # "Alternate allele coverage", "Reference allele coverage",
  # "Reference allele mean expression", "Alternate allele mean expression",
  # "p-value")
  
  colnames(SNP.info) <- c("Chromosome", "Position in human hg38",
                          "Reference allele", "Alternate allele",
                          "rsID", "Population allele frequency gnomAD 3.1.2", 
                          "Genotype in Patient 1", "Genotype in Patient 2",
                          "Genotype in Patient 3", "Genotype in Patient 4", 
                          "Genotype in Patient 5", "Alternate allele coverage", 
                          "Reference allele coverage", "Reference allele mean expression", 
                          "Alternate allele mean expression", "p-value", "Description")
  
  
}

