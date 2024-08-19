#-----------------------------------------------------------------------------------------------
# Author: Vartika Bisht
# Date: 19.08.2024
#
# Description:
# This script is designed to handle and reformat input queries for genomic data analysis. It processes
# user input to determine whether the query is in the format of chromosome:position or gene name, and
# then adjusts the search parameters accordingly. The script performs the following tasks:
# 1. **Reformat Input**: Parses the search query to determine if it is in chromosome:position or gene name format.
# 2. **Variant View**:
#    - **Chromosome:Position**: Checks if the specified position is present in the SuRE database. If not found,
#      it displays an error alert. If found, it sets the region around the position for further analysis.
# 3. **Gene View**:
#    - **Gene Name**: Converts the gene name to uppercase and searches for its coordinates. If the gene is found,
#      it adjusts the region around the gene based on its strand. If not found, it displays an error alert.
# 4. **Invalid Input Handling**: Provides a warning alert for any input format that does not match the expected
#    formats.
#
# Script Workflow:
# 1. **Parse the Query**: Splits the search query based on delimiters and determines the format (chr:pos or gene).
# 2. **Process Based on Format**:
#    - **Chromosome:Position**: Sets up the search region and queries the SuRE database.
#    - **Gene Name**: Converts the gene name, queries for its coordinates, and adjusts the search region based on
#      the gene's strand.
# 3. **Error Handling**:
#    - Displays appropriate alerts for invalid inputs or missing data.
#-----------------------------------------------------------------------------------------------

##################
# Reformat Input #
##################

# Parse the search query and extract relevant information
query_parts <- strsplit(search_query, "[:,]")
query_snps <- NULL

# Reformat input
# query_parts should be chr:pos or gene name
# Check if chr:pos is a data point in our data or not
# If not, show pos - flank, pos + flank region

if (length(query_parts[[1]]) == 2) {
  # chr:pos format
  variant_view <- TRUE
  chr <- query_parts[[1]][1]
  pos <- as.numeric(query_parts[[1]][2])
  pos1 <- pos - flank
  pos2 <- pos + flank

  # Query SuRE database for SNPs at the specified position
  source(file.path(DBQueryScriptsDIR, "DBquery.find.in.SuRE.r"))
  query_snps <- inSuRE(chr, pos, file.path(DBDIR, "SuREViz.Main.db"))

  if (nrow(query_snps) == 0) {
    # Alert user if variant is not found in SuRE database
    shinyalert(
      title = "Variant not in SuRE Database",
      type = "error",
      text = ""
    )
    alert_confirmed(FALSE)
    return(NULL)
  } else {
    alert_confirmed(TRUE)
  }
} else if (length(query_parts[[1]]) == 1) {
  # Gene name format
  variant_view <- FALSE
  gene_name <- toupper(query_parts[[1]][1]) # Convert gene name to uppercase

  # Query gene coordinates from database
  source(file.path(DBQueryScriptsDIR, "DBquery.find.gene.cord.r"))
  gene <- GeneQuery(gene_name, file.path(DBDIR, "gencode.v46.annotation.hg38.genes.db"))

  if (nrow(gene) == 0) {
    # Alert user if gene is not found
    shinyalert(
      type = "error",
      title = "Gene not found",
      text = "Please refer to names as specified in https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz"
    )
    alert_confirmed(FALSE)
    return(NULL)
  } else {
    # Adjust region based on gene strand
    chr <- unique(gene$CHROM)
    pos <- NA
    if (unique(gene$strand) == "+" && length(unique(gene$strand)) == 1) {
      pos1 <- min(gene$start) - flank
      pos2 <- max(gene$end) + flank
    } else {
      pos1 <- min(gene$end) - flank
      pos2 <- max(gene$start) + flank
    }
    alert_confirmed(TRUE)
  }
} else {
  # Alert user for invalid input format
  shinyalert(
    title = "Invalid input format",
    text = "Please provide input in one of the following formats: chr:pos, or gene."
  )
  alert_confirmed(FALSE)
  return(NULL)
}
