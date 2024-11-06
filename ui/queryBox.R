#-----------------------------------------------------------------------------------------------
# Define Query Box Function
# Author: Vartika Bisht
# Date: 19.08.2024
#-----------------------------------------------------------------------------------------------
# This function creates a UI component, `queryBox`, using `fluidRow` and `box` elements from the
# `shiny` and `shinydashboard` packages. The component includes inputs for querying the genome,
# a file upload button, and options for rendering and downloading data. The layout is organized
# into rows and columns to ensure a clean and user-friendly interface.
#-----------------------------------------------------------------------------------------------

queryBox <- function() {
    fluidRow(
        box(
            width = 12, # Box occupies full width of the container
            status = "primary", # Box style: primary color
            solidHeader = TRUE, # Solid header for the box
            title = "Query the human genome", # Title of the box

            # First row with query inputs
            fluidRow(
                column(
                    6,
                    textInput(
                        "search_query", # Input ID for search query
                        "Search", # Label for the input
                        value = "chr12:128797635", # Default value
                        placeholder = "Enter chr:pos or gene name" # Placeholder text
                    )
                ),
                column(
                    6,
                    shinyWidgets::sliderTextInput(
                        "flank", # Input ID for flank selection
                        "Flank", # Label for the input
                        selected = "10kb", # Default selected value
                        choices = lapply(seq(1000, 100000, 1000), to_kb) # Dynamic choices
                    )
                )
            ),

            # Second row with action buttons and file upload
            fluidRow(
                column(
                    12,
                    style = "display: flex; align-items: center;", # Style for aligning items horizontally
                    actionButton(
                        "runButton", # Button ID for rendering
                        "Render", # Button label
                        icon = icon("play"), # Button icon
                        class = "btn-primary", # Button style class
                        style = "margin-left: 10px;" # Inline style for spacing
                    ),
                    shinyFilesButton(
                        "file", # Button ID for file upload
                        "Upload Files", # Button label
                        title = "Select files", # Button title
                        multiple = TRUE, # Allow multiple file selection
                        class = "btn-primary", # Button style class
                        icon = icon("upload"), # Button icon
                        style = "margin-left: 10px;" # Inline style for spacing
                    ),
                    downloadButton(
                        "downloadData", # Button ID for downloading data
                        "Download data/error report", # Button label
                        class = "btn-primary", # Button style class
                        style = "margin-left: 10px;", # Inline style for spacing
                        icon = icon("download") # Button icon
                    )
                ),
                column(
                    12,
                    uiOutput("selected_files"), # Dynamic UI output for displaying selected files
                    style = "margin-top: 10px;" # Inline style for spacing
                )
            )
        )
    )
}
