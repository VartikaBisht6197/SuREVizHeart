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

            # Second row with action buttons
            fluidRow(
                # Second row with action buttons (Render and Download)
                column(
                    12,
                    actionButton(
                        "runButton", # Button ID for rendering
                        "Render", # Button label
                        icon = icon("play"), # Button icon
                        class = "btn-primary", # Button style class
                        style = "margin-top: 10px;" # Inline style for spacing
                    ),
                    downloadButton(
                        "downloadData", # Button ID for downloading data
                        "Download Files", # Button label
                        class = "btn-primary", # Button style class
                        style = "margin-top: 10px;", # Inline style for spacing
                        icon = icon("download") # Button icon
                    )
                )
            ),

            # Third row for uploads and UI for uploaded files
            fluidRow(
                style = "margin-top: 30px;", # Adds 30px space above this row
                column(
                    12,
                    # Add the line and message above the file upload box
                    div(
                        style = "border-top: 2px solid #64BEA5; padding-top: 10px; margin-bottom: 15px;  text-align: center;",
                        span(
                            "Please refer to the 'Uploaded Tracks Viewer' tab for expected format of uploaded files",
                            style = "font-size: 14px; color: #64BEA5; font-weight: bold;"
                        )
                    ),
                    # File upload container with dotted box
                    div(
                        class = "file-upload-container",
                        fileInput(
                            "file",
                            label = NULL,
                            multiple = TRUE,
                            placeholder = "Upload files (max 10 files can be uploaded)",
                            accept = c(".bed", ".bw", ".bigwig", ".bigWig", ".tsv")
                        ),
                        # Display uploaded files
                        uiOutput("selected_files")
                    )
                )
            )

        )
    )
}
