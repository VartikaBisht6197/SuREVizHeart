library(formattable) 

##################
# Reformat Input #
##################

# Parse the search query and extract relevant information
query_parts <- strsplit(search_query, "[:,]")
query_snps <- NULL
# Reformat input
# I have query_parts and flank
# query_parts should be chr:pos or gene name
# check if chr:pos is a data point in our data or not
# If not then just show pos - flank , pos + flank region

if (length(query_parts[[1]]) == 2) {
  # chr:pos format
  variant_view <- TRUE
  chr <- query_parts[[1]][1]
  pos <- as.numeric(query_parts[[1]][2])
  pos1 <- pos - flank
  pos2 <- pos + flank
  
  source(file.path(DBQueryScriptsDIR,"DBquery.find.in.SuRE.r"))
  query_snps = inSuRE(chr,pos,file.path(DBDIR, "SuREViz.Main.db"))
  
  if(nrow(query_snps)==0){
    shinyalert(title = "Variant not in SuRE Database", 
               type = "error",
               text = "")
    alert_confirmed(FALSE)
    return(NULL)
  }else{
    alert_confirmed(TRUE)
  }
  
} else if (length(query_parts[[1]]) == 1) {
  # Gene name format
  variant_view <- FALSE
  gene_name <- query_parts[[1]][1]
  # Convert gene_name to uppercase
  gene_name <- toupper(query_parts[[1]][1])
  # Process the gene name as needed
  source(file.path(DBQueryScriptsDIR,"DBquery.find.gene.cord.r"))
  gene = GeneQuery(gene_name, file.path(DBDIR, "gencode.v46.annotation.hg38.genes.db"))
  #if gene is not found, show shiny alert warning
  if (nrow(gene) == 0) {
    shinyalert(type = "error",
               title = "Gene not found", 
               text = "Please refer to names as specified in https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz")
    alert_confirmed(FALSE)
    return(NULL)
  } else {
    # if gene is found, if strand is positive then start-end otherwise end-start
    chr = unique(gene$CHROM)
    pos = NA
    if(unique(gene$strand)=="+" && length(unique(gene$strand))==1){
      pos1 = min(gene$start) - flank
      pos2 = max(gene$end) + flank
    } else {
      pos1 = min(gene$end) - flank
      pos2 = max(gene$start) + flank
    }
    alert_confirmed(TRUE)
  }
} else {
  # Show shiny alert warning for invalid input
  shinyalert(title = "Invalid input format", 
              text = "Please provide input in one of the following formats: chr:pos, or gene.")
  alert_confirmed(FALSE)
  return(NULL)
  
}