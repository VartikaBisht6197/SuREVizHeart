
# CLIN info
clin_complete_url <- NA
  if(!is.na(pos)){
    source(file.path(DBQueryScriptsDIR,"DBquery.clin.r"))
    clindata <- suREVizClinVarQuery(query_snps$CHROM, query_snps$POS, query_snps$REF , query_snps$ALT , file.path(DBDIR,"SuREViz.ClinVar.db"))
    if(nrow(clindata)>0){
      
      # Define the base URL for CLIVVAR variation data
      base_url <- "https://www.ncbi.nlm.nih.gov/clinvar/variation/"
      
      # Define the clinker ID you want to query (e.g., 1221899)
      clinker_id <- clindata$ID
      
      # Construct the complete URL
      clin_complete_url <- paste0(base_url, clinker_id)
    } else {
      # Point to the local HTML file
      clin_complete_url <- file.path(DBDIR, "not_in_clinvar.html")
    }
  } else {
    # Point to the local HTML file
    clin_complete_url <- file.path(DBDIR, "variant_view_false.html")
    
  }