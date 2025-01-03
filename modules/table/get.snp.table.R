
highlight = NA
SNP.info <- suredata
col.order <- c(
  "Chromosome", "Position",
  "REF", "ALT",
  "rsID", "Population allele frequency gnomAD 3.1.2",
  "SuREX38", "SuREX57", "SuREX59", "SuREX67", "SuREX68", "SuREX86",
  "REF allele coverage", "ALT allele coverage",
  "REF allele mean expression", "ALT allele mean expression",
  "p-value", "Description"
)

if( nrow(suredata) != 0 ){
  if(!is.na(pos)){
    colnames(SNP.info) <- c(
      "Chromosome", "Position",
      "REF", "ALT",
      "rsID", "Population allele frequency gnomAD 3.1.2",
      "SuREX38", "SuREX57", "SuREX59",
      "SuREX67", "SuREX68",
      "SuREX86", "ALT allele coverage",
      "REF allele coverage", "REF allele mean expression",
      "ALT allele mean expression", "p-value", "Description"
    )

    SNP.info <- SNP.info[, col.order]

    # highlight the pos
    highlight = which(SNP.info$Position == pos)
  }  
  message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", "SuRE Table : Functionally assessed variants in the region of interest printed in tabular format.  ✅"))
} else {
  SNP.info <- data.frame(matrix(data = "", nrow = 1, ncol = length(col.order)))
  colnames(SNP.info) <- col.order

  message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", "SuRE Table : No functionally assessed variants in the region of interest."))
}
