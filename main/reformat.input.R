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
hg38.chrom.sizes = fread(file.path(dataDIR, "hg38.chrom.sizes"))
query_parts <- strsplit(search_query, "[:,]")
query_snps <- NULL

# Function to validate input
validate_number <- function(chr,pos) {
  # Check if the input is numeric and matches the desired pattern (digits only)
  if (grepl("^[0-9]+$", pos) && (as.numeric(pos)>=1) && (as.numeric(pos)<=hg38.chrom.sizes[hg38.chrom.sizes$V1 == chr,]$V2)) {
    return(TRUE) # Input is valid
  } else {
    return(FALSE) # Input is invalid
  }
}

# Function to validate chromosome name
validate_chromosome <- function(chromosome) {
  # Check if the chromosome matches the pattern "chr_" followed by 1-22 or X
  if (grepl("^chr(1[0-9]|2[0-2]|[1-9]|X)$", chromosome)) {
    return(TRUE)  # Valid chromosome
  } else {
    return(FALSE) # Invalid chromosome
  }
}

# Reformat input
# query_parts should be chr:pos or gene name
# Check if chr:pos is a data point in our data or not
# If not, show pos - flank, pos + flank region

if (length(query_parts[[1]]) == 2) {
  # chr:pos format
  variant_view <- TRUE
  chr <- query_parts[[1]][1]
  pos <- query_parts[[1]][2]

  if(!(validate_number(chr,pos) && validate_chromosome(chr)) ){

    # Alert user for invalid input format
    message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", "Invalid input format. Please provide input in one of the following formats: chr:pos, or gene. ❌"))
    shinyalert(
      html = TRUE,
      title = "Invalid Input Format",
      text = paste0(
        "<div style='text-align: left; font-family: Arial, sans-serif; font-size: 14px;'>",
        "<strong>Oops! Your input appears to be in an incorrect format.</strong><br><br>",
        "To search for a variant, please follow the required format:<br>",
        "<ul style='margin-left: 20px;'>",
        "  <li><strong>chr:pos</strong></li>",
        "    <ul style='margin-left: 20px;'>",
        "      <li><code>chr</code>: Chromosome, formatted as <code>chrN</code> (e.g., <code>chr1</code>, <code>chrX</code>, etc.).</li>",
        "      <li><code>pos</code>: A valid numeric position on the chromosome.</li>",
        "    </ul>",
        "</ul><br>",
        "<strong>Important tips:</strong><br>",
        "<ul style='margin-left: 20px;'>",
        "  <li>Ensure there are no unnecessary spaces before or after your input.</li>",
        "  <li>The position must lie within the valid range for the chromosome size. You can refer to the size definitions at ",
        "<a href='http://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes' target='_blank'>",
        "UCSC Genome Browser</a>.</li>",
        "  <li>Double-check for typos or incorrect characters.</li>",
        "</ul><br>",
        "<strong>Your Input:</strong> <code>", chr, ":", pos, "</code><br><br>",
        "Please correct your input and try again.<br>",
        "</div>"
      ),
      type = "error"
    )
    input_pass_test(FALSE)
    return(NULL)

  }else{

    pos <- as.numeric(pos)
    pos1 <- pos - flank
    pos2 <- pos + flank

    if (pos1 < 1) {
      pos1 = 1
    }

    if (pos2 > hg38.chrom.sizes[hg38.chrom.sizes$V1 == chr,]$V2){
      pos2 = as.numeric(hg38.chrom.sizes[hg38.chrom.sizes$V1 == chr, ]$V2)
    }

    # Query SuRE database for SNPs at the specified position
    source(file.path(DBQueryScriptsDIR, "DBquery.find.in.SuRE.r"))
    query_snps <- inSuRE(chr, pos, file.path(dataDIR, "SuREViz.Main.db"))

    if (nrow(query_snps) == 0) {
      # Alert user if variant is not found in SuRE database
      # Inform the user through a shinyalert
      shinyalert(
        html = TRUE,
        title = "Variant Not Found",
        type = "info",
        text = paste0(
          "<div style='text-align: left; font-family: Arial, sans-serif; font-size: 14px;'>",
          "Unfortunately, we couldn't locate the variant at the specified position: <strong>", chr, ":", pos, "</strong> in the SuRE database.<br><br>",
          "But no worries! We're automatically expanding the search to include the surrounding region to help you find relevant results.<br>",
          "</div>"
        )
      )
      # Adjust the region to flank regions
      pos1 <- pos - flank
      pos2 <- pos + flank
      pos <- NA # Set pos to NA as requested

      # Continue with the adjusted region (render or other processing)
      # You can now render everything or proceed with the region you just calculated
      input_pass_test(TRUE)

      return(NULL)
    } else {
      input_pass_test(TRUE)
    }
  }
  
} else if (length(query_parts[[1]]) == 1) {
  # Gene name format
  variant_view <- FALSE
  gene_name <- toupper(query_parts[[1]][1]) # Convert gene name to uppercase

  # Query gene coordinates from database
  source(file.path(DBQueryScriptsDIR, "DBquery.find.gene.cord.r"), local = TRUE)
  gene <- GeneQuery(gene_name, file.path(dataDIR, "gencode.v46.annotation.hg38.genes.db"))

  if (nrow(gene) == 0) {
    # Alert user if gene is not found
    message("Either the search query is a single string but not a gene name or the gene specified is not defined in GENECODE GTF GRCh38.p14 Human Release 46 ❌")
    shinyalert(
      html = TRUE,
      type = "error",
      title = "Invalid Input Format",
      text = paste0(
        "<div style='text-align: left; font-family: Arial, sans-serif; font-size: 14px;'>",
        "<strong>Oops! Your input appears to be in an incorrect format.</strong><br><br>",
        "<strong>To search for a gene:</strong><br>",
        "<ul style='margin-left: 20px;'>",
        "  <li>Use the gene names as specified in <strong>GENCODE GTF GRCh38.p14 Human Release 46</strong>.</li>",
        "  <li>Ensure there are no extra spaces at the beginning or end of your input.</li>",
        "</ul><br>",
        "Please double-check your input and try again.<br>",
        "</div>"
      )
    )
    input_pass_test(FALSE)
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
    input_pass_test(TRUE)
  }
} else {
  # Alert user for invalid input format
  message( "Invalid input format.Please provide input in one of the following formats: chr:pos, or gene. ❌")
  shinyalert(
    html = TRUE,
    title = "Invalid Input Format",
    text = paste0(
      "<div style='text-align: left; font-family: Arial, sans-serif; font-size: 14px;'>",
      "<strong>Invalid Input Format</strong><br><br>",
      "Please provide your input in one of the following formats:<br>",
      "<ul style='margin-left: 20px;'>",
      "  <li><strong>chr:pos</strong> - Specify the chromosome and position (e.g., <code>chr12:128797635</code>).</li>",
      "  <li><strong>gene</strong> - Provide a valid gene name.</li>",
      "</div>"
    ),
    type = "error"
  )

  input_pass_test(FALSE)
  return(NULL)
}
