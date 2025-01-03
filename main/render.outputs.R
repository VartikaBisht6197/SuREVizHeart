#-----------------------------------------------------------------------------------------------
# Author: Vartika Bisht
# Date: 19.08.2024
#
# Description:
# This script is responsible for rendering and displaying various output plots and tables in the
# Shiny application. It includes functionalities for:
# 1. Displaying Plotly plots for SNP and gene visualizations.
# 2. Rendering static plots for BigWig and gene expression data.
# 3. Showing user-defined BigWig plots.
# 4. Displaying JASPAR motif plots with adjustable dimensions.
# 5. Rendering FIMO table data with styled HTML tables.
# 6. Embedding ClinVar and gnomAD pages either from local files or external URLs.
# 7. Displaying a table of SNP information with optional row highlighting.
# 8. Generating and displaying a table plot using Plotly for SNP information.
#
# Script Workflow:
# 1. **Output Plots and Tables**:
#    - **Plotly Plot**: Renders a combined Plotly plot with SNP and gene visualizations, allowing
#      interactivity with the `plotly_click` event.
#    - **BigWig Plot**: Displays the static BigWig plot.
#    - **Gene Expression Plot**: Renders the gene expression plots.
#    - **User Defined BigWig Plot**: Displays plots for user-defined BigWig files.
#    - **JASPAR Plot**: Renders JASPAR motif plots with specific dimensions.
#    - **FIMO Table**: Renders the FIMO table using styled HTML.
#    - **ClinVar Page**: Embeds or includes ClinVar data from local files or URLs.
#    - **gnomAD Page**: Embeds or includes gnomAD data from local files or URLs.
#    - **SNP Information Table**: Renders and optionally highlights rows in a table of SNP information.
#    - **Table Plot**: Creates and displays a Plotly table for SNP information.
#
#-----------------------------------------------------------------------------------------------

# Output plots and tables

# Render combined Plotly plot for SNP and gene visualizations
output$plotTab1 <- renderPlotly({
   p <- subplot(
          style(snp.plot %>%
                config(
                  displayModeBar = TRUE, # Show the modebar
                  modeBarButtons = list(list("resetViews"), list("toImage"))
                ), showlegend = FALSE),
          style(gene.plotly %>%
                config(
                  displayModeBar = TRUE, # Show the modebar
                 modeBarButtons = list( list("resetViews"), list("toImage"))
                ), showlegend = FALSE),
          nrows = 2, heights = c(0.8, 0.2), shareX = TRUE
        ) %>%
       layout(xaxis = list(autorange = FALSE, range = c(pos1, pos2)))
  # Register the plotly_click event for interactivity
  event_register(p, "plotly_click")
})

# Output static BigWig plot
output$plotTab2 <- renderPlot(
  {
    bigwig.plot
  },
  res = 100
)

# Output static gene expression plot
output$plotTab3 <- renderPlot(
  {
    gene.expression.plots
  },
  res = 100,
  height = ngenes * 75
)

# Output user-defined BigWig plot
output$plotTab5 <- renderPlot(
  {
    user.defined.plot
  },
  res = 100,
  height = 500 * user_define_plots
)

# Output JASPAR plot with adjustable dimensions
output$jasparplot <- renderPlot(
  {
    # Ensure the plot dimensions are consistent with the motif height and count
    jaspar_output$plot
  },
  res = 100,
  height = 600 * motifs
) # Adjust height and width as needed

# Output FIMO table with styled HTML
output$fimo_table <- renderText({
  if (nrow(jaspar_output$table) != 0) {
    jaspar_output$table %>%
      knitr::kable("html", escape = FALSE) %>%
      kable_styling(
        full_width = FALSE,
        bootstrap_options = c("striped", "hover", "condensed", "responsive")
      ) %>%
      column_spec(1,
        bold = TRUE, color = "black",
        background = "#64BEA5",
        border_right = TRUE,
        extra_css = "text-align: center;"
      ) %>%
      column_spec(2:ncol(jaspar_output$table),
        border_left = TRUE, border_right = TRUE,
        extra_css = "text-align: center;"
      ) %>%
      row_spec(0, bold = TRUE, color = "black", background = "#64BEA5") %>%
      row_spec(1:nrow(jaspar_output$table), background = "white") %>%
      kable_minimal()
  } else {
    jaspar_output$table %>%
      knitr::kable("html", escape = FALSE)
  }
})

# Output ClinVar page
output$clintable <- renderUI({
  if (is.na(clin_complete_url)) {
    return(NULL) # Return NULL if URL is NA
  } else if (file.exists(clin_complete_url)) {
    includeHTML(clin_complete_url) # Render local HTML file
  } else {
    tags$iframe(
      src = clin_complete_url,
      width = "100%", # Adjust width as needed
      height = "800px", # Adjust height as needed
      frameborder = "0",
      sandbox = "allow-same-origin allow-scripts allow-popups" # Add sandbox attributes
    )
  }
})

# Output gnomAD page
output$gnomad <- renderUI({
  if (is.na(gnomad_complete_url)) {
    return(NULL) # Return NULL if URL is NA
  } else if (file.exists(gnomad_complete_url)) {
    includeHTML(gnomad_complete_url) # Render local HTML file
  } else {
    tags$iframe(
      src = gnomad_complete_url,
      width = "100%", # Adjust width as needed
      height = "800px", # Adjust height as needed
      frameborder = "0",
      sandbox = "allow-same-origin allow-scripts allow-popups" # Add sandbox attributes
    )
  }
})

# Output SNP information table with optional row highlighting
output$tableTab3 <- renderText({

  if (!is.na(highlight)) {
    (SNP.info) %>%
      knitr::kable("html") %>%
      kable_styling(
        full_width = TRUE,
        bootstrap_options = c("striped", "hover", "condensed")
      ) %>%
      kable_minimal() %>%
      row_spec(highlight, background = "yellow")
  } else {
    (SNP.info) %>%
      knitr::kable("html") %>%
      kable_styling(
        full_width = TRUE,
        bootstrap_options = c("striped", "hover", "condensed")
      ) %>%
      kable_minimal()
  }
})

# Output table plot using Plotly for SNP information
table_plot <- plot_ly(
  type = "table",
  cells = list(
    values = list(
      colnames(SNP.info),
      as.character(t(SNP.info))
    )
  )
)
