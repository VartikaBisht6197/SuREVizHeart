# gnomAD info
gnomad_complete_url <- NA
if(!is.na(pos)){
  source(file.path(DBQueryScriptsDIR,"DBquery.clin.r"))
  if( !is.na(query_snps$rsID) ){
    
    log.this.info("gnomAD URL View : Variant found in gnomAD database. ✅")

    # Define the base URL for CLIVVAR variation data
    base_url <- "https://gnomad.broadinstitute.org/variant/"
    
    # Define the clinker ID you want to query (e.g., 1221899)
    clinker_id <- paste(substr(query_snps$CHROM,4,nchar(query_snps$CHROM)),query_snps$POS,query_snps$REF,query_snps$ALT,sep = "-")
    
    # Construct the complete URL
    gnomad_complete_url <- paste0(base_url, clinker_id)
  } else {
    log.this.info("gnomAD URL View : Variant not found in gnomAD database. ⚠️")
    # Point to the local HTML file
    gnomad_complete_url <- file.path(dataDIR, "not_in_gnomAD.html")
  }
} else {
    log.this.info("gnomAD URL View : No variant selected. ⚠️")
    # Point to the local HTML file
    gnomad_complete_url <- file.path(dataDIR, "variant_view_false.html")
    
  }