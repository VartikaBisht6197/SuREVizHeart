library(kableExtra)

# Output plots and tables
output$plotTab1 <- renderPlotly({
      # Plotly plot
      p <- subplot(snp.plot, gene.plotly, nrows = 2, heights = c(0.8, 0.2), shareX = TRUE) %>%
       layout(xaxis = list(autorange = FALSE, range = c(pos1, pos2)))
      # Register the plotly_click event
      event_register(p, "plotly_click")
})

# Output bigwigs
# Bigwigs Quite big
# output$plotTab2 <- renderPlot({
#   bigwig.plot
# },  res = 100)


output$plotTab3 <- renderPlot({
  gene.expression.plots
},  res = 100)

# Output userdefined bigwigs
output$plotTab5 <- renderPlot({
  user.defined.plot
}, res = 100)
print(motifs)
output$jasparplot <- renderPlot(
  {
    # Ensure the plot dimensions are consistent
    jaspar_output$plot
  },
  res = 100,
  height = motifvh*motifs,
  width = 1000
) # Adjust height and width as needed


output$fimo_table <- renderText({
  
  if(nrow(jaspar_output$table)!=0){
    jaspar_output$table %>%
      knitr::kable("html", escape = FALSE) %>%
      kable_styling(
        full_width = FALSE,
        bootstrap_options = c("striped", "hover", "condensed", "responsive")
      ) %>%
      # Center align all columns
      column_spec(1, bold = TRUE, color = "black", 
                  background = "#64BEA5",
                  border_right = TRUE, 
                  extra_css = "text-align: center;") %>%
      column_spec(2:ncol(jaspar_output$table), 
                  border_left = TRUE, border_right = TRUE, 
                  extra_css = "text-align: center;") %>%
      row_spec(0, bold = TRUE, color = "black", background = "#64BEA5") %>%
      row_spec(1:nrow(jaspar_output$table), background = "white") %>%
      kable_minimal()
      } else {
        jaspar_output$table %>%
          knitr::kable("html", escape = FALSE) 
    }
  
})


output$clintable <- renderUI({
  if (is.na(clin_complete_url)) {
    # Return NULL if clin_complete_url is NA (this should never happen with the current logic)
    return(NULL)
  } else if (file.exists(clin_complete_url)) {
    # Render local HTML file if clin_complete_url points to a file
    includeHTML(clin_complete_url)
  } else {
    # Create an iframe element that embeds the URL
    tags$iframe(
      src = clin_complete_url,
      width = "100%", # Adjust width as needed
      height = "800px", # Adjust height as needed
      frameborder = "0",
      # Add the 'sandbox' attribute to allow popups
      sandbox = "allow-same-origin allow-scripts allow-popups"
    )
  }
})


output$gnomad <- renderUI({
  if (is.na(gnomad_complete_url)) {
    # Return NULL if gnomad_complete_url is NA
    return(NULL)
  } else if (file.exists(gnomad_complete_url)) {
    # Render local HTML file if clin_complete_url points to a file
    includeHTML(gnomad_complete_url)
  } else {
    # Create an iframe element that embeds the URL
    tags$iframe(
      src = gnomad_complete_url,
      width = "100%", # Adjust width as needed
      height = "800px", # Adjust height as needed
      frameborder = "0",
      # Add the 'sandbox' attribute to allow popups
      sandbox = "allow-same-origin allow-scripts allow-popups"
    )
  }
})


output$tableTab3 <- renderText({
  if (nrow(suredata) == 0) {
    showNotification("Data not found")
    return("Data not found")  # Return a message or an empty string
  }
  
  if(!is.na(highlight)){
    (SNP.info) %>%
      knitr::kable("html") %>% 
      kable_styling(full_width = TRUE,
                    bootstrap_options = c("striped", "hover", "condensed")
      ) %>% kable_minimal() %>%
    row_spec(highlight, background = "yellow")
    
  }else{
    (SNP.info) %>%
      knitr::kable("html") %>% 
      kable_styling(full_width = TRUE,
                    bootstrap_options = c("striped", "hover", "condensed")
      ) %>% kable_minimal()
  }
  
})



# Print the combined plots along with the table
# Create table using plot_ly
# Create the table using plot_ly
table_plot <- plot_ly(
  type = "table",
  cells = list(
    values = list(
      colnames(SNP.info),
      as.character(t(SNP.info))
    )
  )
)
