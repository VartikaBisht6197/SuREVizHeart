#-----------------------------------------------------------------------------------------------
# Define Tab Content Function
# Author: Vartika Bisht
# Date: 19.08.2024
#-----------------------------------------------------------------------------------------------
# This function defines the content for each tab in the Shiny application's dashboard. It uses
# `tabItems` and `tabItem` functions from the `shinydashboard` package to create multiple tabs,
# each containing a `box` with various output components. Each tab represents a different
# section of the application, displaying plots, HTML content, or tables.
#-----------------------------------------------------------------------------------------------

tabContent <- function() {
    tabItems(
        # Tab: Functional Impact Assessment
        tabItem(
            tabName = "impact_assessment",
            box(
                width = 12, status = "primary", solidHeader = TRUE, title = "Functional Impact Assessment",
                includeHTML(file.path(appDIR, "www", "SuREViz.html")), # Include HTML content for impact assessment
                div(
                    style = "overflow: auto;",
                    withSpinner(plotlyOutput("plotTab1", height = "1200px", width = "100%")) # Plotly output with spinner for loading
                )
            )
        ),

        # Tab: Variant Data Overview
        tabItem(
            tabName = "variant_summary",
            box(
                width = 12, status = "primary", solidHeader = TRUE, title = "Variant Data Overview",
                includeHTML(file.path(appDIR, "www", "variant.summary.html")), # Include HTML content for variant summary
                div(
                    style = "overflow: auto;",
                    withSpinner(htmlOutput("tableTab3")) # HTML output with spinner for loading
                )
            )
        ),

        # Tab: Gene Expression Overview
        tabItem(
            tabName = "gene_expression_summary",
            box(
                width = 12, status = "primary", solidHeader = TRUE, title = "Gene Expression Overview",
                includeHTML(file.path(appDIR, "www", "gene.expression.summary.html")), # Include HTML content for gene expression
                div(
                    style = "overflow: auto;",
                    withSpinner(plotOutput("plotTab3", width = "100%", height = "800px")) # Plot output with spinner for loading
                )
            )
        ),

        # Tab: Integrating Complementary Datasets with SuRE Data
        tabItem(
            tabName = "sure_profiles",
            box(
                width = 12, status = "primary", solidHeader = TRUE, title = "SuRE Profiles",
                includeHTML(file.path(appDIR, "www", "SuREprofileView.html")), # Include HTML content for SuRE profile view
                div(
                    style = "overflow: auto;",
                    withSpinner(plotOutput("plotTab2", width = "100%", height = "1800px")) # Plot output with spinner for loading
                )
            )
        ),

        # Tab: Transcription Factor Disruption Analysis
        tabItem(
            tabName = "jaspar",
            box(
                width = 12, status = "primary", solidHeader = TRUE, title = "Transcription Factor Binding Site Impact (TFBSi)",
                includeHTML(file.path(appDIR, "www", "TF.impact.html")), # Include HTML content for TF disruption analysis
                div(
                    style = "overflow: auto; width: 100%; text-align: center; margin-bottom: 20px;",
                    withSpinner(plotOutput("jasparplot", width = "100%", height = paste0(motifvh, "px"))) # Plot output with spinner for loading
                ),
                div(
                    style = "overflow: auto; width: 100%; margin-top: 10px; text-align: center; margin-top: 20px;",
                    withSpinner(tableOutput("fimo_table")) # Table output with spinner for loading
                )
            )
        ),

        # Tab: Genome Aggregation Database (gnomAD) Viewer
        tabItem(
            tabName = "gnomad",
            box(
                width = 12, status = "primary", solidHeader = TRUE, title = "Genome Aggregation Database (gnomAD) Viewer",
                includeHTML(file.path(appDIR, "www", "population.genetics.html")), # Include HTML content for gnomAD viewer
                div(
                    style = "overflow: auto;",
                    withSpinner(htmlOutput("gnomad")) # HTML output with spinner for loading
                )
            )
        ),

        # Tab: ClinVar Viewer
        tabItem(
            tabName = "clinvar",
            box(
                width = 12, status = "primary", solidHeader = TRUE, title = "ClinVar Viewer",
                includeHTML(file.path(appDIR, "www", "ClinVar.html")), # Include HTML content for ClinVar viewer
                div(
                    style = "overflow: auto;",
                    withSpinner(htmlOutput("clintable")) # HTML output with spinner for loading
                )
            )
        ),

        # Tab: Uploaded Tracks Viewer
        tabItem(
            tabName = "comparative_study",
            box(
                width = 12, status = "primary", solidHeader = TRUE, title = "Uploaded Tracks Viewer",
                includeHTML(file.path(appDIR, "www", "UserUpload.html")), # Include HTML content for variant summary
                div(
                    style = "overflow: auto;",
                    withSpinner(plotOutput("plotTab5", width = "100%", height = "3000px")) # Plot output with spinner for loading
                )
            )
        )
    )
}
