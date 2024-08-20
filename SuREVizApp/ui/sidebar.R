#-----------------------------------------------------------------------------------------------
# Define Sidebar Function
# Author: Vartika Bisht
# Date: 19.08.2024
#-----------------------------------------------------------------------------------------------
# This function defines the sidebar layout for the Shiny application using the `dashboardSidebar`
# function from the `shinydashboard` package. It sets up the sidebar menu with various navigation
# items and sub-items, each linked to a specific tab in the dashboard. Icons are used to represent
# different sections visually.
#-----------------------------------------------------------------------------------------------

sidebar <- function() {
    dashboardSidebar(
        width = 300, # Set the width of the sidebar
        sidebarMenu(
            # Menu item for Functional Impact Assessment
            menuItem(
                "Functional Impact Assessment", # Title displayed in the sidebar
                tabName = "impact_assessment", # Tab name associated with this menu item
                icon = icon("chart-bar") # Icon displayed next to the menu item
            ),

            # Menu item for Variant Data Overview
            menuItem(
                "Variant Data Overview",
                tabName = "variant_summary",
                icon = icon("table")
            ),

            # Menu item for Gene Expression Overview
            menuItem(
                "Gene Expression Overview",
                tabName = "gene_expression_summary",
                icon = icon("chart-simple")
            ),

            # Menu item for Integrating Complementary Datasets
            menuItem(
                "Integrating Complementary Datasets",
                tabName = "sure_profiles",
                icon = icon("dna")
            ),

            # Menu item for Highlighted Variant with sub-items
            menuItem(
                "Highlighted Variant",
                tabName = "highlighted_variant", # This main menu item links to a tab name
                icon = icon("vial"), # Icon for the main menu item
                style = "background-color: #e9e9e9;", # Background color for the highlighted menu item
                menuSubItem(
                    "Transcription Factor Disruption Analysis", # Sub-menu item title
                    tabName = "jaspar" # Tab name associated with this sub-item
                ),
                menuSubItem(
                    "Genome Aggregation Database (gnomAD) Viewer",
                    tabName = "gnomad"
                ),
                menuSubItem(
                    "ClinVar Viewer",
                    tabName = "clinvar"
                )
            ),

            # Menu item for Uploaded Tracks Viewer
            menuItem(
                "Uploaded Tracks Viewer",
                tabName = "comparative_study",
                icon = icon("project-diagram")
            )
        )
    )
}
