
# CLIN info
clin_complete_url <- NA
  if(!is.na(pos)){
    source(file.path(DBQueryScriptsDIR,"DBquery.clin.r"))
    clindata <- suREVizClinVarQuery(query_snps$CHROM, query_snps$POS, query_snps$REF , query_snps$ALT , file.path(dataDIR,"SuREViz.ClinVar.db"))
    if(nrow(clindata)>0){

      message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", "ClinVar URL View : Variant found in clinvar database. ✅"))
      
      # Define the base URL for CLIVVAR variation data
      base_url <- "https://www.ncbi.nlm.nih.gov/clinvar/variation/"
      
      # Define the clinker ID you want to query (e.g., 1221899)
      clinker_id <- clindata$ID
      
      # Construct the complete URL
      clin_complete_url <- paste0(base_url, clinker_id)
    } else {
      message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", "ClinVar URL View : Variant not found in clinvar database. ⚠️"))
      # Point to the local HTML file
      clin_complete_url <- file.path(appDIR, "www", "not_in_clinvar.html")
    }
  } else {
    message(paste(format(Sys.time(), "%d/%m/%Y %H:%M:%S"), ":", "ClinVar URL View : No variant selected. ⚠️"))
    # Point to the local HTML file
    clin_complete_url <- file.path(appDIR, "www", "variant_view_false.html")
    
  }