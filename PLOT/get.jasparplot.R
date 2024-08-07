if( !is.na(pos) && nrow(suredata)!=0 ){
  # Run grep command to check for rows
  grep_output <- system(paste0("grep -w ", pos, " ", file.path(DBDIR, "JASPAR.TF.raQTL.txt")), intern = TRUE)
  
  # Check if grep command returns any rows
  if (length(grep_output) > 0) {
    # If there are rows, read data into jaspar.data
    jaspar.data <- fread(cmd = paste0("grep -w ", pos, " ", file.path(DBDIR, "JASPAR.TF.raQTL.txt")), sep = "\t")[, 1:14]
    colnames(jaspar.data) = c("CHROM", "POS", "REF", "ALT", "start", "end", "strand", "pos", "motif_id", "motif_alt_id", "refs.score", "alts.score", "refs.pval", "alts.pval")
    
    if(nrow(jaspar.data)==1 & (jaspar.data$CHROM == chr & jaspar.data$POS == pos & jaspar.data$REF == query_snps$REF & jaspar.data$ALT == query_snps$ALT)) {
      JASPAR.PPM = read_meme(file.path(DBDIR,"data","JASPAR2022_CORE_vertebrates_non-redundant_v2.meme"))
      names(JASPAR.PPM) = unlist(lapply(JASPAR.PPM,function(x) x@name))

      source(file.path(plotDIR,"jasparplot.R"))
      source(file.path(plotDIR,"TF.RefAlt.binding.plot.R"))
      result_list = jaspar.analysis(jaspar.data, JASPAR.PPM)
      
      # Display extra plots
      if (!is.null(result_list) && length(result_list) > 0) {
        rv <- reactiveValues()  # Create a reactiveValues object to store plots
        
        output$tableTab4 <- renderUI({
          plot_list <- lapply( names(result_list), function(plot_name){plotOutput(plot_name)} )
        })
        
        observe({
          # Need local so that each item gets its own number. Without it, the value
          # of i in the renderPlot() will be the same across all instances, because
          # of when the expression is evaluated.
          for (plot_name in names(result_list)) {
            local({
              plot_name_local <- plot_name
              output[[plot_name_local]] <- renderPlot({
                result_list[[plot_name_local]]
              })
            })
          }
        })
        
        
      }}
  }else{
    output$tableTab4 <- renderUI({
      ggplotly(ggplot() + theme_minimal())
    })
    
    
  }}else{
    output$tableTab4 <- renderUI({
      ggplotly(ggplot() + theme_minimal())
    })
    
  }

