jaspar_output = list("plot" = NULL,
                    "table" = data.frame())
jaspar.data = NULL

if( !is.na(pos) && nrow(suredata)!=0 ){
  
  # Run grep command to check for rows
  grep_output <- system(paste0("grep -w ", pos, " ", file.path(DBDIR, "JASPAR2022.HLHS.raQTL.TF.filter.txt")), intern = TRUE)
  
  # Check if grep command returns any rows
  if (length(grep_output) > 0) {
    # If there are rows, read data into jaspar.data
    
    jaspar.data <- fread(cmd = paste0("grep -w ", pos, " ", file.path(DBDIR, "JASPAR2022.HLHS.raQTL.TF.filter.txt")), sep = "\t")
    colnames(jaspar.data) <- colnames(fread(cmd = paste("head -n1",file.path(DBDIR, "JASPAR2022.HLHS.raQTL.TF.filter.txt")), header = TRUE))
    
    if(nrow(jaspar.data)==1 & (jaspar.data$CHROM == chr & jaspar.data$POS == pos & jaspar.data$REF == query_snps$REF & jaspar.data$ALT == query_snps$ALT)) {
      JASPAR.PPM = read_meme(file.path(DBDIR,"data","JASPAR2022_CORE_vertebrates_non-redundant_v2.meme"))
      names(JASPAR.PPM) = unlist(lapply(JASPAR.PPM,function(x) x@name))
      
      source(file.path(plotDIR,"jasparplot.R"), local = TRUE)
      jaspar_output = jaspar.analysis(jaspar.data, JASPAR.PPM) 
      motifs = (nrow(jaspar_output$table) + 2)/plots_in_one_set  }
  } else {
    if(query_snps$raQTL == "raQTL"){
      jaspar_output$plot = ggplot() + 
        theme_minimal() +
        annotate("text", x = 0.5, y = 0.5, label = "The variant does not affect any TFBS\n(as defined by JASPAR 2022 CORE) significantly.", size = 6, hjust = 0.5, vjust = 0.5) +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()
        )
    } else {
      jaspar_output$plot = ggplot() + 
        theme_minimal() +
        annotate("text", x = 0.5, y = 0.5, label = "The variant selected is not an raQTL\nThis view is only possible for raQTLs.", size = 6, hjust = 0.5, vjust = 0.5) +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()
        )
    }
    
  }
} else {
  jaspar_output$plot = ggplot() + 
    theme_minimal() +
    annotate("text", x = 0.5, y = 0.5, label = "Display is only possible for the variant view.", size = 6, hjust = 0.5, vjust = 0.5) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
      ) 
}

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
   
# gnomAD info
gnomad_complete_url <- NA
if(!is.na(pos)){
  source(file.path(DBQueryScriptsDIR,"DBquery.clin.r"))
  if( !is.na(query_snps$rsID) ){
    
    # Define the base URL for CLIVVAR variation data
    base_url <- "https://gnomad.broadinstitute.org/variant/"
    
    # Define the clinker ID you want to query (e.g., 1221899)
    clinker_id <- paste(substr(query_snps$CHROM,4,nchar(query_snps$CHROM)),query_snps$POS,query_snps$REF,query_snps$ALT,sep = "-")
    
    # Construct the complete URL
    gnomad_complete_url <- paste0(base_url, clinker_id)
  } else {
    # Point to the local HTML file
    gnomad_complete_url <- file.path(DBDIR, "not_in_gnomAD.html")
  }
} else {
    # Point to the local HTML file
    gnomad_complete_url <- file.path(DBDIR, "variant_view_false.html")
    
  }

